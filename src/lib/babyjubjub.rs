use std::marker::PhantomData;
use plonk::composer::Variable;
use plonk::{Composer};
use crate::lib::mux::Mux;
use crate::utils::Field;

pub struct Babyjubjub<F: Field>{
    _marker: PhantomData<F>,
}

impl<F: Field> Babyjubjub<F> {
    /// input: P1(x1, y1), P2(x2, y2)
    /// output: (x, y) = P1 + P2
    /// xout == (x1y2 + y1x2) / (1 + 168696x1y2y1x2)
    /// yout == (168700 * x1 * y2 - y1 * x2 + (y1 - 168700 * x1)*(x2 + y2)) / (1 - 168696 * x1 * y2 * y1 * x2)
    /// = y1y2 - 168700x1x2 / ()
    pub fn add(
        cs: &mut Composer<F>,
        x1: Variable,
        y1: Variable,
        x2: Variable,
        y2: Variable,
    ) -> (Variable, Variable) {
        let x1y2 = cs.mul(x1, y2);
        let x2y1 = cs.mul(x2, y1);
        let x1x2 = cs.mul(x1, x2);
        let y1y2 = cs.mul(y1, y2);

        let x1y2_val = cs.get_assignment(x1y2);
        let x2y1_val = cs.get_assignment(x2y1);
        let x1x2_val = cs.get_assignment(x1x2);
        let y1y2_val = cs.get_assignment(y1y2);

        //x out
        let tmp1_val:F = F::one() + F::from(168696 as u128) * x1y2_val * x2y1_val;
        let tmp1 = cs.alloc(tmp1_val);
        cs.poly_gate(
            vec![(x1y2, F::zero()), (x2y1, F::zero()), (tmp1, -F::one())],
            F::from(168696 as u128),
            F::one(),
        );

        let xout_val:F = (x1y2_val + x2y1_val) / tmp1_val;
        let xout = cs.alloc(xout_val);
        cs.poly_gate(
            vec![(tmp1, F::zero()), (xout, F::zero()), (x1y2, -F::one()), (x2y1, -F::one())],
            F::one(),
            F::zero(),
        );

        //y out
        let tmp2_val:F = F::one() - F::from(168696 as u128) * x1y2_val * x2y1_val;
        let tmp2 = cs.alloc(tmp2_val);
        cs.poly_gate(
            vec![(x1y2, F::zero()), (x2y1, F::zero()), (tmp2, -F::one())],
            -F::from(168696 as u128),
            F::one(),
        );

        let yout_val:F = (y1y2_val - F::from(168700 as u128) * x1x2_val) / tmp2_val;
        let yout = cs.alloc(yout_val);
        cs.poly_gate(
            vec![(tmp2, F::zero()), (yout, F::zero()), (y1y2, -F::one()), (x1x2, F::from(168700 as u128))],
            F::one(),
            F::zero(),
        );

        (xout, yout)
    }

    /// input: P1(x1, y1)
    /// output: (x, y) = 2 * P1
    /// xout == (2 * x1y1) / (1 + 168696x1x1y1y1)
    /// yout == y1y1 - 168700x1x1 / (1 - 168696x1x1y1y1)
    pub fn double(
        cs: &mut Composer<F>,
        x1: Variable,
        y1: Variable,
    ) -> (Variable, Variable) {
        //单独实现也只省一个门，算了
        Babyjubjub::add(cs, x1, y1, x1, y1)
    }

    /// input: P1(x1, y1)
    /// output: (x, y) = n * P1
    /// n is bit format
    pub fn mul_scalar(
        cs: &mut Composer<F>,
        x1: Variable,
        y1: Variable,
        n: &[Variable],
    ) -> (Variable, Variable) {
        assert_eq!(n.len(), 254);
        let y_zero = cs.alloc(F::one());

        let mut sum_x = Composer::<F>::null();
        let mut sum_y = y_zero;

        let mut xi = x1;
        let mut yi = y1;

        for &i in n {
            let x = Mux::mux1(cs, Composer::<F>::null(), xi, i);
            let y = Mux::mux1(cs, y_zero, yi, i);

            let (tmpl, tmpr) = Babyjubjub::add(cs, x, y, sum_x, sum_y);
            sum_x = tmpl;
            sum_y = tmpr;

            let (tmpl, tmpr) = Babyjubjub::double(cs, xi, yi);
            xi = tmpl;
            yi = tmpr;

        }

        (sum_x, sum_y)
    }

    /// input: P1(x1, y1)
    /// check if P1 is on babyjub
    #[allow(dead_code)]
    pub fn check_point(
        cs: &mut Composer<F>,
        x1: Variable,
        y1: Variable,
    ) -> bool {
        let x1x1 = cs.mul(x1, x1);
        let y1y1 = cs.mul(y1, y1);

        //168700x2 + y2 === 1 + 168696x2*y2;
        cs.poly_gate(
            vec![(x1x1, -F::one()), (y1y1, -F::from(168700 as u128))],
            F::from(168696 as u128),
            F::one(),
        );

        true
    }
}


#[cfg(test)]
mod tests {
    use std::str::FromStr;
    use plonk::MFr;
    use ark_std::{test_rng};
    use plonk::*;

    use ark_poly::{EvaluationDomain};
    use num_bigint::{BigUint};
    use crate::lib::babyjubjub::Babyjubjub;

    use crate::lib::mimc::MiMC;
    use crate::lib::num_to_bits::NumToBits;

    #[test]
    fn test_add() -> Result<(), Error>{
        let mut cs = {
            // x^3 + x + pi = 35
            let mut cs = Composer::new(4);
            let pi = cs.alloc_input(MFr::from(5 as u64));
            let x = cs.alloc(MFr::from(3));
            let y = cs.mul(x, x);
            let z = cs.mul(x, y);
            let u = cs.add(x, z);
            let v = cs.add(pi, u);
            cs.enforce_constant(v, MFr::from(BigUint::parse_bytes(b"23", 16).unwrap()));


            let p1x = cs.alloc(MFr::from(0));
            let p1y = cs.alloc(MFr::from(1));
            let p2x = cs.alloc(MFr::from(0));
            let p2y = cs.alloc(MFr::from(1));

            let (right_x, right_y) = Babyjubjub::add(&mut cs, p1x, p1y, p2x, p2y);
            println!("{}", cs.get_assignment(right_x));
            println!("{}", cs.get_assignment(right_y));

            let Bx = cs.alloc(
                MFr::from_str(
                    "995203441582195749578291179787384436505546430278305826713579947235728471134"
                ).unwrap_or_default()
            );
            let By = cs.alloc(
                MFr::from_str(
                    "5472060717959818805561601436314318772137091100104008585924551046643952123905"
                ).unwrap_or_default()
            );
            println!("B8X{}", MFr::from_str(
                "5299619240641551281634865583518297030282874472190772894086521144482721001553"
            ).unwrap_or_default());
            println!("B8Y{}", MFr::from_str(
                "16950150798460657717958625567821834550301663161624707787222815936182638968203"
            ).unwrap_or_default());

            let (Ax2, Ay2) = Babyjubjub::double(&mut cs, Bx, By);
            let (Ax4, Ay4) = Babyjubjub::double(&mut cs, Ax2, Ay2);
            let (Ax8, Ay8) = Babyjubjub::double(&mut cs, Ax4, Ay4);
            println!("{}", cs.get_assignment(Ax8));
            println!("{}", cs.get_assignment(Ay8));

            cs
        };

        let rng = &mut test_rng();

        let pk = cs.compute_prover_key::<GeneralEvaluationDomain<MFr>>()?;
        let mut prover = plonk::prover::Prover::new(&pk);
        prover.prove(&mut cs, rng)?;

        Ok(())
    }

    #[test]
    fn test_mul() -> Result<(), Error>{
        let mut cs = {
            // x^3 + x + pi = 35
            let mut cs = Composer::new(4);
            let pi = cs.alloc_input(MFr::from(5 as u64));
            let x = cs.alloc(MFr::from(3 as u128));
            let y = cs.mul(x, x);
            let z = cs.mul(x, y);
            let u = cs.add(x, z);
            let v = cs.add(pi, u);
            cs.enforce_constant(v, MFr::from(BigUint::parse_bytes(b"23", 16).unwrap()));


            let Px = cs.alloc(
                MFr::from_str(
                    "17777552123799933955779906779655732241715742912184938656739573121738514868268"
                ).unwrap_or_default()
            );
            let Py = cs.alloc(
                MFr::from_str(
                    "2626589144620713026669568689430873010625803728049924121243784502389097019475"
                ).unwrap_or_default()
            );

            //mul 3
            println!("3{}", MFr::from_str(
                "19372461775513343691590086534037741906533799473648040012278229434133483800898"
            ).unwrap_or_default());
            println!("3{}", MFr::from_str(
                "9458658722007214007257525444427903161243386465067105737478306991484593958249"
            ).unwrap_or_default());

            let h2bits = NumToBits::num_to_bits(&mut cs, x, 254);
            let (Ax2, Ay2) = Babyjubjub::mul_scalar(&mut cs, Px, Py, &h2bits);
            println!("{}", cs.get_assignment(Ax2));
            println!("{}", cs.get_assignment(Ay2));


            //mul rr
            let rr = cs.alloc(
                MFr::from_str(
                    "14035240266687799601661095864649209771790948434046947201833777492504781204499"
                ).unwrap_or_default()
            );
            println!("*rr{}", MFr::from_str(
                "17070357974431721403481313912716834497662307308519659060910483826664480189605"
            ).unwrap_or_default());
            println!("*rr{}", MFr::from_str(
                "4014745322800118607127020275658861516666525056516280575712425373174125159339"
            ).unwrap_or_default());

            let h2bits = NumToBits::num_to_bits(&mut cs, rr, 254);
            let (Ax2, Ay2) = Babyjubjub::mul_scalar(&mut cs, Px, Py, &h2bits);
            println!("{}", cs.get_assignment(Ax2));
            println!("{}", cs.get_assignment(Ay2));


            cs
        };

        let rng = &mut test_rng();

        let pk = cs.compute_prover_key::<GeneralEvaluationDomain<MFr>>()?;
        let mut prover = plonk::prover::Prover::new(&pk);
        prover.prove(&mut cs, rng)?;

        Ok(())
    }
}