#![allow(non_snake_case)]
use std::marker::PhantomData;

use plonk::{Composer};
use plonk::composer::Variable;
use crate::lib::babyjubjub::Babyjubjub;
use crate::lib::comparators::enforce_equal_if_enabled;
use crate::lib::mimc::MiMC;
use crate::lib::num_to_bits::NumToBits;
use crate::utils::Field;

pub struct EdDSA<F: Field>{
    _marker: PhantomData<F>,
}

impl<F: Field> EdDSA<F> {
    pub fn mimc_sponge_verifier(
        cs: &mut Composer<F>,
        Ax: Variable,
        Ay: Variable,
        S: Variable,
        R8x: Variable,
        R8y: Variable,
        M: Variable,
        selector: Variable,
    ){
        // Ensure S<Subgroup Order


        // Calculate the h = H(R,A, msg)
        let hash = MiMC::MiMC_sponge(
            cs,
            &[R8x, R8y, Ax, Ay, M],
            Composer::<F>::null(),
            1,
        );

        // Calculate second part of the right side:  right2 = h*8*A
        let (Ax2, Ay2) = Babyjubjub::double(cs, Ax, Ay);
        let (Ax4, Ay4) = Babyjubjub::double(cs, Ax2, Ay2);
        let (Ax8, Ay8) = Babyjubjub::double(cs, Ax4, Ay4);

        // check that A is not zero (Ax != 0)
        assert!(!cs.get_assignment(Ax8).is_zero());
        let Ax8_inv_val = cs.get_assignment(Ax8).inverse().unwrap();
        let Ax8_inv = cs.alloc(Ax8_inv_val);
        cs.poly_gate(
            vec![(Ax8_inv, F::zero()), (Ax8, F::zero())],
            F::one(),
            -F::one()
        );

        let h2bits = NumToBits::num_to_bits(cs, hash[0], 254);

        let (right2x, right2y) = Babyjubjub::mul_scalar(cs, Ax8, Ay8, &h2bits);

        // Compute the right side: right =  R8 + right2
        let (right_x, right_y) = Babyjubjub::add(cs, right2x, right2y, R8x, R8y);

        // Calculate left side of equation left = S*B8

        let B8x = cs.alloc(
            F::from_str(
                "5299619240641551281634865583518297030282874472190772894086521144482721001553"
            ).unwrap_or_default()
        );
        let B8y = cs.alloc(
            F::from_str(
                "16950150798460657717958625567821834550301663161624707787222815936182638968203"
            ).unwrap_or_default()
        );

        let S2bits = NumToBits::num_to_bits(cs, S, 254);

        let (left_x, left_y) = Babyjubjub::mul_scalar(cs, B8x, B8y, &S2bits);

        // Do the comparation left == right if enabled;
        enforce_equal_if_enabled(cs, left_x, right_x, selector);
        enforce_equal_if_enabled(cs, left_y, right_y, selector);
    }
}

#[cfg(test)]
mod tests {
    use plonk::MFr;
    use ark_std::{test_rng};
    use plonk::*;

    use ark_ff::{One};
    use num_bigint::{BigUint};

    use crate::lib::eddsa::EdDSA;

    #[test]
    fn verify_signature() -> Result<(), Error>{
        let mut cs = {
            let mut cs = Composer::new(4);

            let m = MFr::from(BigUint::parse_bytes(b"42649378395939397566720", 10).unwrap());
            let r8x = MFr::from(BigUint::parse_bytes(b"11384336176656855268977457483345535180380036354188103142384839473266348197733", 10).unwrap());
            let r8y = MFr::from(BigUint::parse_bytes(b"15383486972088797283337779941324724402501462225528836549661220478783371668959", 10).unwrap());
            let s = MFr::from(BigUint::parse_bytes(b"1868336918738674306327358602987493427631678603535639134028485964115448322340", 10).unwrap());
            let ax = MFr::from(BigUint::parse_bytes(b"13277427435165878497778222415993513565335242147425444199013288855685581939618", 10).unwrap());
            let ay = MFr::from(BigUint::parse_bytes(b"13622229784656158136036771217484571176836296686641868549125388198837476602820", 10).unwrap());

            let m = cs.alloc(m);
            let r8x = cs.alloc(r8x);
            let r8y = cs.alloc(r8y);
            let s = cs.alloc(s);
            let ax = cs.alloc(ax);
            let ay = cs.alloc(ay);

            let selector = cs.alloc(MFr::one());

            EdDSA::mimc_sponge_verifier(
                &mut cs,
                ax,
                ay,
                s,
                r8x,
                r8y,
                m,
                selector,
            );

            cs
        };

        let rng = &mut test_rng();

        let pk = cs.compute_prover_key::<GeneralEvaluationDomain<MFr>>()?;
        let mut prover = plonk::prover::Prover::new(&pk);
        prover.prove(&mut cs, rng)?;

        Ok(())
    }
}