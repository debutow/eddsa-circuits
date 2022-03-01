use ark_ff::{BigInteger, FftField, FpParameters, PrimeField};
use crate::utils::Field;
use std::cmp;
use std::marker::PhantomData;
use ark_std::Zero;

use num_bigint::BigUint;
use crate::{Composer, Error};
use crate::composer::Variable;


pub trait FpField: FftField + PrimeField {
}

pub(crate) struct NonnativeField<F: Field, FP: FpField> {
    pub A: F,
    ///A^2
    pub A2: F,
    ///A^3
    pub A3: F,

    pub P: BigNum<F>,
    ///P - R
    pub P_mod_R: F,
    ///T - P
    pub minus_P_mod_T: BigNum<F>,
    ///A^4 mod P
    pub A4_mod_P: BigNum<F>,
    ///A^5 mod P
    pub A5_mod_P: BigNum<F>,
    ///A^6 mod P
    pub A6_mod_P: BigNum<F>,

    // pub cs: &Composer<F>,

    _pfp: PhantomData<FP>,
}

#[derive(Clone)]
pub struct BigNum<F: Field>{
    ///a = a.0 + a.1 * A + a.2 * A^2 + a.3 * A^3
    value: [F; 4],

    //store the bigUint form for easy calculating
    value_biguint: BigUint,

    ///bits of value
    value_bits: [u32; 4], //必须用统一的算法得出

    value_alloc: [Variable; 4],

    num_bits: u32, //必须用统一的算法得出

}

// const A: dyn Field = Field::zero();

impl<F: Field> BigNum<F> {
    pub fn new(value: &[F; 4]) -> BigNum<F> {


        let Abigu: BigUint = BigUint::from(1u128 << 68);

        let mut vbignum: BigUint = Abigu.clone() * value[3].into_repr().into();
        vbignum += value[2].into_repr().into();
        vbignum *= Abigu.clone();
        vbignum += value[1].into_repr().into();
        vbignum *= Abigu.clone();
        vbignum += value[0].into_repr().into();

        BigNum{
            value: value.clone(),
            value_biguint: vbignum,
            value_bits: [68,68,68,68],
            value_alloc: [Composer::<F>::null(),Composer::<F>::null(),Composer::<F>::null(),Composer::<F>::null()],
            num_bits: 272,
        }
    }

    pub fn new_from_biguint(value: &BigUint) -> BigNum<F> {
        let Abigu: BigUint = BigUint::from(1u128 << 68);

        let a0 = value.clone() % Abigu.clone();
        let tmp = value.clone() / Abigu.clone();
        let a1 = tmp.clone() % Abigu.clone();
        let tmp = tmp.clone() / Abigu.clone();
        let a2 = tmp.clone() % Abigu.clone();
        let tmp = tmp.clone() / Abigu.clone();
        let a3 = tmp.clone() % Abigu.clone();

        let values = [F::from(a0), F::from(a1), F::from(a2), F::from(a3)];

        // let num_bits = BigNum::cal_num_bits(&values);

        BigNum{
            value: values,
            value_biguint: value.clone(),
            value_bits: [68,68,68,68],
            value_alloc: [Composer::<F>::null(),Composer::<F>::null(),Composer::<F>::null(),Composer::<F>::null()],
            num_bits: 272,
        }
    }

    ///1. sort num_bits, 2. calculate from small
    pub fn cal_add_bits(value: &[u32]) -> u32 {

        let mut sorted = value.iter().map(|x| {*x}).collect::<Vec<_>>();

        // let mut sorted = value.clone();
        sorted.sort();
        // println!("{:?}", sorted);

        let mut num = sorted[1] + 1;
        for i in 2..sorted.len() {
            num = cmp::max(num, sorted[i]) + 1;
        }

        // println!("{}", num);
        num
    }

    pub fn cal_num_bits(value: &[u32; 4]) -> u32 {
        let mut tmp = value.clone();
        tmp[1] = value[1] + 68;
        tmp[2] = value[2] + 68*2;
        tmp[3] = value[3] + 68*3;

        BigNum::<F>::cal_add_bits(&tmp)
    }

    ///only once
    pub fn alloc_vars(&mut self, cs: &mut Composer<F>) {
        if self.value_alloc[0] != Composer::<F>::null() {return; }

        for i in 0..4 {
            self.value_alloc[i] = cs.alloc(self.value[i]);
        }

    }


    pub fn get_value(&self) -> [F; 4] {
        self.value
    }

    ///will not change anything else
    pub fn set_value(&mut self, value: &[F; 4]) {
        self.value = value.clone();

        //should not
        // let res = BigNum::<F>::cal_num_bits(value);
        // self.set_num_bits(res);
    }

    pub fn get_num_bits(&self) -> u32 {
        self.num_bits
    }

    pub fn set_num_bits(&mut self, num_bits: u32) {
        self.num_bits = num_bits;
    }

    pub fn get_value_bits(&self) -> [u32; 4] {
        self.value_bits
    }

    pub fn set_value_bits(&mut self, value_bits: [u32; 4]) {
        self.value_bits = value_bits;
    }


}

impl<F: Field, FP: FpField> NonnativeField<F, FP> {
    pub fn new(cs: &mut Composer<F>) -> NonnativeField<F, FP> {
        let Abigu: BigUint = BigUint::from(1u128 << 68);
        println!("A {}", Abigu);
        let A2bigu: BigUint = Abigu.clone() * Abigu.clone();
        let A3bigu: BigUint = A2bigu.clone() * Abigu.clone();
        let Tbigu = A3bigu.clone() * Abigu.clone();
        println!("T {}", Tbigu);


        let fp = FP::Params::MODULUS;
        let fr = F::Params::MODULUS;

        let p_mod_r: BigUint = fp.into() - fr.into();

        let t_sub_p: BigUint = Tbigu.clone() - fp.into();
        let mut t_sub_p_mybigint = BigNum::<F>::new_from_biguint(&t_sub_p.clone());
        t_sub_p_mybigint.alloc_vars(cs);
        t_sub_p_mybigint.set_value_bits(
            [t_sub_p_mybigint.value[0].into_repr().num_bits(),
                t_sub_p_mybigint.value[1].into_repr().num_bits(),
                t_sub_p_mybigint.value[2].into_repr().num_bits(),
                t_sub_p_mybigint.value[3].into_repr().num_bits()]
        );
        t_sub_p_mybigint.set_num_bits(BigNum::<F>::cal_num_bits(&t_sub_p_mybigint.value_bits));

        let mut p_mybigint = BigNum::<F>::new_from_biguint(&fp.into());
        p_mybigint.alloc_vars(cs);
        p_mybigint.set_value_bits(
            [p_mybigint.value[0].into_repr().num_bits(),
            p_mybigint.value[1].into_repr().num_bits(),
            p_mybigint.value[2].into_repr().num_bits(),
            p_mybigint.value[3].into_repr().num_bits()]
        );
        p_mybigint.set_num_bits(BigNum::<F>::cal_num_bits(&p_mybigint.value_bits));

        let mut A4_mod_P = BigNum::<F>::new_from_biguint(&(Tbigu.clone() % fp.into()));
        let mut A5_mod_P = BigNum::<F>::new_from_biguint(&(Tbigu.clone() * Abigu.clone() % fp.into()));
        let mut A6_mod_P = BigNum::<F>::new_from_biguint(&(Tbigu.clone() * Abigu.clone() * Abigu.clone() % fp.into()));
        A4_mod_P.alloc_vars(cs);
        A5_mod_P.alloc_vars(cs);
        A6_mod_P.alloc_vars(cs);
        A4_mod_P.set_value_bits(
            [A4_mod_P.value[0].into_repr().num_bits(),
                A4_mod_P.value[1].into_repr().num_bits(),
                A4_mod_P.value[2].into_repr().num_bits(),
                A4_mod_P.value[3].into_repr().num_bits()]
        );
        A5_mod_P.set_value_bits(
            [A5_mod_P.value[0].into_repr().num_bits(),
                A5_mod_P.value[1].into_repr().num_bits(),
                A5_mod_P.value[2].into_repr().num_bits(),
                A5_mod_P.value[3].into_repr().num_bits()]
        );
        A6_mod_P.set_value_bits(
            [A6_mod_P.value[0].into_repr().num_bits(),
                A6_mod_P.value[1].into_repr().num_bits(),
                A6_mod_P.value[2].into_repr().num_bits(),
                A6_mod_P.value[3].into_repr().num_bits()]
        );
        A4_mod_P.set_num_bits(BigNum::<F>::cal_num_bits(&A4_mod_P.value_bits));
        A5_mod_P.set_num_bits(BigNum::<F>::cal_num_bits(&A5_mod_P.value_bits));
        A6_mod_P.set_num_bits(BigNum::<F>::cal_num_bits(&A6_mod_P.value_bits));


        NonnativeField{
            A: F::from(1u128 << 68),
            A2: F::from(A2bigu),
            A3: F::from(A3bigu),
            P: p_mybigint,
            P_mod_R: F::from(p_mod_r),
            minus_P_mod_T: t_sub_p_mybigint,
            A4_mod_P,
            A5_mod_P,
            A6_mod_P,
            _pfp: Default::default()
        }
    }

    pub fn is_divide_exactly_by_T(&self, bignum: &BigNum<F>, cs: &mut Composer<F>) -> Result<bool, Error> {

        let b_u0 = 1 + cmp::max(bignum.value_bits[0], bignum.value_bits[1] + 68);
        let b_u1 = 1 + cmp::max(bignum.value_bits[2], bignum.value_bits[3] + 68);

        if b_u0 < 254 && b_u1 < 254 {
            let u0: F = bignum.value[0] + bignum.value[1] * self.A;
            let u1: F = bignum.value[2] + bignum.value[3] * self.A;
            //在F上没法判断是否整除，因为一定会有域上的结果。所以转换到Bigint做除法
            let u0bigu = u0.into_repr().into();
            let u1bigu = u1.into_repr().into();

            if !(u0bigu.clone() % self.A2.into_repr().into()).is_zero() {
                return Ok(false);
            }
            let v0bigu = u0bigu.clone() / self.A2.into_repr().into();
            let B_v0 = cmp::max(b_u0 - 68*2, 0);

            if !((u1bigu.clone() + v0bigu.clone()) % self.A2.into_repr().into()).is_zero() {
                return Ok(false);
            }
            let v1bigu = (u1bigu.clone() + v0bigu.clone()) / self.A2.into_repr().into();
            let B_u1v0 = cmp::max(b_u1, B_v0) + 1;
            let B_v1 = cmp::max(B_u1v0 - 68*2, 0);

            //约束u0 = t0 + t1*A
            //约束u1 = t2 + t3*A
            //约束u0 = v0*A2 和v0范围
            //约束u1 + v0 = v1*A2 和v1范围
            let u0_alloc = cs.alloc(u0);
            let u1_alloc = cs.alloc(u1);
            let v0 = F::from(v0bigu);
            let v1 = F::from(v1bigu);
            let v0_alloc = cs.alloc(v0);
            let v1_alloc = cs.alloc(v1);

            cs.poly_gate(
                vec![(bignum.value_alloc[0], F::one()),
                     (bignum.value_alloc[1], self.A),
                     (u0_alloc, -F::one()),
                     (Composer::<F>::null(), F::zero())],
                F::zero(), F::zero());
            cs.poly_gate(
                vec![(bignum.value_alloc[2], F::one()),
                     (bignum.value_alloc[3], self.A),
                     (u1_alloc, -F::one()),
                     (Composer::<F>::null(), F::zero())],
                F::zero(), F::zero());
            cs.poly_gate(
                vec![(u0_alloc, -F::one()),
                     (v0_alloc, self.A2),
                     (Composer::<F>::null(), F::zero()),
                     (Composer::<F>::null(), F::zero())],
                F::zero(), F::zero());
            cs.poly_gate(
                vec![(v0_alloc, -F::one()),
                     (v1_alloc, self.A2),
                     (u1_alloc, -F::one()),
                     (Composer::<F>::null(), F::zero())],
                F::zero(), F::zero());

            cs.enforce_range(v0_alloc, (B_v0 as usize + 1) / 2 * 2)?;
            cs.enforce_range(v1_alloc, (B_v1 as usize + 1) / 2 * 2)?;


        } else {

            let t0bigu = bignum.value[0].into_repr().into();
            let t1bigu = bignum.value[1].into_repr().into();
            let t2bigu = bignum.value[2].into_repr().into();
            let t3bigu = bignum.value[3].into_repr().into();

            if !(t0bigu.clone() % self.A.into_repr().into()).is_zero() { return Ok(false); }
            let w0bigu = t0bigu.clone() / self.A.into_repr().into();
            let B_w0 = cmp::max(bignum.value_bits[0] - 68, 0);

            if !((t1bigu.clone() + w0bigu.clone()) % self.A.into_repr().into()).is_zero() { return Ok(false); }
            let v0bigu = (t1bigu.clone() + w0bigu.clone()) / self.A.into_repr().into();
            let B_t1w0 = cmp::max(bignum.value_bits[1], B_w0) + 1;
            let B_v0 = cmp::max(B_t1w0 - 68, 0);

            if !((v0bigu.clone() + t2bigu.clone()) % self.A.into_repr().into()).is_zero() { return Ok(false); }
            let w1bigu = (v0bigu.clone() + t2bigu.clone()) / self.A.into_repr().into();
            let B_t2v0 = cmp::max(bignum.value_bits[2], B_v0) + 1;
            let B_w1 = cmp::max(B_t2v0 - 68, 0);

            if !((w1bigu.clone() + t3bigu.clone()) % self.A.into_repr().into()).is_zero() { return Ok(false); }
            let v1bigu = (w1bigu.clone() + t3bigu.clone()) / self.A.into_repr().into();
            let B_t3w1 = cmp::max(bignum.value_bits[3], B_w1) + 1;
            let B_v1 = cmp::max(B_t3w1 - 68, 0);

            //约束t0 = w0*A
            //约束w0 + t1 = v0*A 和w0,v0范围
            //约束v0 + t2 = w1*A
            //约束w1 + t3 = v1*A 和w1,v1范围
            let w0 = F::from(w0bigu);
            let w1 = F::from(w1bigu);
            let w0_alloc = cs.alloc(w0);
            let w1_alloc = cs.alloc(w1);
            let v0 = F::from(v0bigu);
            let v1 = F::from(v1bigu);
            let v0_alloc = cs.alloc(v0);
            let v1_alloc = cs.alloc(v1);

            cs.poly_gate(
                vec![(bignum.value_alloc[0], -F::one()),
                     (w0_alloc, self.A),
                     (Composer::<F>::null(), F::zero()),
                     (Composer::<F>::null(), F::zero())],
                F::zero(), F::zero());
            cs.poly_gate(
                vec![(bignum.value_alloc[1], -F::one()),
                     (v0_alloc, self.A),
                     (w0_alloc, -F::one()),
                     (Composer::<F>::null(), F::zero())],
                F::zero(), F::zero());
            cs.poly_gate(
                vec![(bignum.value_alloc[2], -F::one()),
                     (w1_alloc, self.A),
                     (v0_alloc, -F::one()),
                     (Composer::<F>::null(), F::zero())],
                F::zero(), F::zero());
            cs.poly_gate(
                vec![(bignum.value_alloc[3], -F::one()),
                     (v1_alloc, self.A),
                     (w1_alloc, -F::one()),
                     (Composer::<F>::null(), F::zero())],
                F::zero(), F::zero());

            cs.enforce_range(v0_alloc, (B_v0 as usize +1) / 2 * 2)?;
            cs.enforce_range(w0_alloc, (B_w0 as usize +1) / 2 * 2)?;
            cs.enforce_range(v1_alloc, (B_v1 as usize +1) / 2 * 2)?;
            cs.enforce_range(w1_alloc, (B_w1 as usize +1) / 2 * 2)?;

        }


        Ok(true)
    }

    ///c is a const
    pub fn is_eq_c_mod_P(&self, bignum: &BigNum<F>, c: FP, cs: &mut Composer<F>) -> Result<bool, Error> {

        let fp = FP::Params::MODULUS;

        let p_sub_c: BigUint = fp.into() - c.into_repr().into();
        let p_sub_c_mybigint = BigNum::<F>::new_from_biguint(&p_sub_c);

        //calculate p (where a-c-(p-1)P=0)
        if !((bignum.value_biguint.clone() + p_sub_c.clone()) % fp.into()).is_zero() { return Ok(false); }
        let p: BigUint = (bignum.value_biguint.clone() + p_sub_c.clone()) / fp.into();
        println!("p {}", p);
        let mut p_mybigint = BigNum::<F>::new_from_biguint(&p);
        println!("p_mybigint");
        for item in p_mybigint.value {
            println!("{}", item);
        }
        p_mybigint.alloc_vars(cs);

        let mut B_p = 0;
        if cmp::max(bignum.num_bits, p_sub_c.bits() as u32) + 1 >= 254 {
            B_p = cmp::max(bignum.num_bits, p_sub_c.bits() as u32) + 1 - 254
        }
        let mut B_pi = [0u32, 0u32, 0u32, 0u32];
        if B_p < 68 {
            B_pi[0] = B_p;
        } else {
            B_pi[0] = 68;
            let tmp = B_p - 68;
            if tmp < 68 {
                B_pi[1] = tmp;
            } else {
                B_pi[1] = 68;
                let tmp = tmp - 68;
                if tmp < 68 {
                    B_pi[2] = tmp;
                } else {
                    B_pi[2] = 68;
                    B_pi[3] = tmp - 68;
                }
            }
        }
        p_mybigint.set_value_bits(B_pi);
        p_mybigint.set_num_bits(BigNum::<F>::cal_num_bits(&p_mybigint.value_bits));


        let ar: BigUint = bignum.value_biguint.clone() % F::Params::MODULUS.into();
        let cr: BigUint = c.into_repr().into() % F::Params::MODULUS.into();
        let pr: BigUint = p.clone() % F::Params::MODULUS.into();

        let ar = F::from(ar);
        let cr = F::from(cr);
        let pr = F::from(pr);
        assert!((ar - cr - (pr - F::one()) * self.P_mod_R).is_zero());

        let psubcr = p_sub_c.clone() % F::Params::MODULUS.into();

        //约束ar + cr’ - pr * self.P_mod_R = 0
        // let tmp0: F = bignum.value[0] + p_sub_c_mybigint.value[0] - p_mybigint.value[0] * self.P_mod_R;
        // let tmp1: F = bignum.value[1] + p_sub_c_mybigint.value[1] - p_mybigint.value[1] * self.P_mod_R;
        // let tmp2: F = bignum.value[2] + p_sub_c_mybigint.value[2] - p_mybigint.value[2] * self.P_mod_R;
        // let tmp3: F = bignum.value[3] + p_sub_c_mybigint.value[3] - p_mybigint.value[3] * self.P_mod_R;
        // let tmp4: F = tmp2 * self.A2 + tmp1 * self.A + tmp0;
        // let tmp5: F = tmp3 * self.A3 + tmp4;
        // println!("tmp5 {}", tmp5);
        let tmp0: F = bignum.value[0] + bignum.value[1] * self.A + bignum.value[2] * self.A2;
        let tmp1: F = tmp0 + bignum.value[3] * self.A3 - p_mybigint.value[0] * self.P_mod_R;
        let tmp2: F = tmp1 - p_mybigint.value[1]  * self.A*self.P_mod_R - p_mybigint.value[2]  * self.A2*self.P_mod_R;
        // let tmp3: F = tmp2 - p_mybigint.value[3]  * self.A3*self.P_mod_R + F::from(p_sub_c.clone()) ;

        let tmp0_alloc = cs.alloc(tmp0);
        let tmp1_alloc = cs.alloc(tmp1);
        let tmp2_alloc = cs.alloc(tmp2);
        // let tmp3_alloc = cs.alloc(tmp3);

        cs.poly_gate(
            vec![(tmp0_alloc, -F::one()),
                 (bignum.value_alloc[1], self.A),
                 (bignum.value_alloc[2], self.A2),
                 (bignum.value_alloc[0], F::one())],
            F::zero(), F::zero());
        cs.poly_gate(
            vec![(tmp1_alloc, -F::one()),
                 (tmp0_alloc, F::one()),
                 (bignum.value_alloc[3], self.A3),
                 (p_mybigint.value_alloc[0], -self.P_mod_R)],
            F::zero(), F::zero());
        cs.poly_gate(
            vec![(tmp2_alloc, -F::one()),
                 (tmp1_alloc, F::one()),
                 (p_mybigint.value_alloc[1], -self.A*self.P_mod_R),
                 (p_mybigint.value_alloc[2], -self.A2*self.P_mod_R)],
            F::zero(), F::zero());
        cs.poly_gate(
            vec![(Composer::<F>::null(), -F::one()),
                 (tmp2_alloc, F::one()),
                 (Composer::<F>::null(), F::zero()),
                 (p_mybigint.value_alloc[3], -self.A3*self.P_mod_R)],
            F::zero(), F::from(psubcr));

        cs.enforce_range(p_mybigint.value_alloc[0], (p_mybigint.value_bits[0] as usize +1) / 2 * 2)?;
        cs.enforce_range(p_mybigint.value_alloc[1], (p_mybigint.value_bits[1] as usize +1) / 2 * 2)?;
        cs.enforce_range(p_mybigint.value_alloc[2], (p_mybigint.value_bits[2] as usize +1) / 2 * 2)?;
        cs.enforce_range(p_mybigint.value_alloc[3], (p_mybigint.value_bits[3] as usize +1) / 2 * 2)?;

        //a + c’ + p * self.minus_P_mod_T = 0 mod T
        //calculate t
        //t0 = a0 + psc0 + p0*Pt0
        //t1 = a1 + psc1 + p0*Pt1 + p1*Pt0
        //t2 = a2 + psc2 + p0*Pt2 + p1*Pt1 + p2*Pt0
        //t3 = a3 + psc3 + p0*Pt3 + p1*Pt2 + p2*Pt1 + p3*Pt0
        let t0 = bignum.value[0] + p_sub_c_mybigint.value[0]
            + p_mybigint.value[0] * self.minus_P_mod_T.value[0];

        let t1 = bignum.value[1] + p_sub_c_mybigint.value[1]
            + p_mybigint.value[0] * self.minus_P_mod_T.value[1]
            + p_mybigint.value[1] * self.minus_P_mod_T.value[0];

        let t2 = bignum.value[2] + p_sub_c_mybigint.value[2]
            + p_mybigint.value[0] * self.minus_P_mod_T.value[2]
            + p_mybigint.value[1] * self.minus_P_mod_T.value[1]
            + p_mybigint.value[2] * self.minus_P_mod_T.value[0];

        let t3 = bignum.value[3] + p_sub_c_mybigint.value[3]
            + p_mybigint.value[0] * self.minus_P_mod_T.value[3]
            + p_mybigint.value[1] * self.minus_P_mod_T.value[2]
            + p_mybigint.value[2] * self.minus_P_mod_T.value[1]
            + p_mybigint.value[3] * self.minus_P_mod_T.value[0];
        // println!("t0123 {} {} {} {}", t0, t1, t2, t3);

        let mut t_mybignum = BigNum::<F>::new(&[t0,t1,t2,t3]);
        t_mybignum.alloc_vars(cs);
        t_mybignum.set_value_bits(
            [
                BigNum::<F>::cal_add_bits(&[bignum.value_bits[0],
                    p_sub_c_mybigint.value_bits[0],
                    p_mybigint.value_bits[0]+self.minus_P_mod_T.value_bits[0],
                ]),
                BigNum::<F>::cal_add_bits(&[bignum.value_bits[1],
                    p_sub_c_mybigint.value_bits[1],
                    p_mybigint.value_bits[0]+self.minus_P_mod_T.value_bits[1],
                    p_mybigint.value_bits[1]+self.minus_P_mod_T.value_bits[0],
                ]),
                BigNum::<F>::cal_add_bits(&[bignum.value_bits[2],
                    p_sub_c_mybigint.value_bits[2],
                    p_mybigint.value_bits[0]+self.minus_P_mod_T.value_bits[2],
                    p_mybigint.value_bits[1]+self.minus_P_mod_T.value_bits[1],
                    p_mybigint.value_bits[2]+self.minus_P_mod_T.value_bits[0],
                ]),
                BigNum::<F>::cal_add_bits(&[bignum.value_bits[3],
                    p_sub_c_mybigint.value_bits[3],
                    p_mybigint.value_bits[0]+self.minus_P_mod_T.value_bits[3],
                    p_mybigint.value_bits[1]+self.minus_P_mod_T.value_bits[2],
                    p_mybigint.value_bits[2]+self.minus_P_mod_T.value_bits[1],
                    p_mybigint.value_bits[3]+self.minus_P_mod_T.value_bits[0],
                ]),
            ]
        );
        t_mybignum.set_num_bits(BigNum::<F>::cal_num_bits(&t_mybignum.value_bits));

        //约束ti
        cs.poly_gate(
            vec![(t_mybignum.value_alloc[0], -F::one()),
                 (bignum.value_alloc[0], F::one()),
                 (Composer::<F>::null(), F::zero()),
                 (p_mybigint.value_alloc[0], self.minus_P_mod_T.value[0])],
            F::zero(), p_sub_c_mybigint.value[0]);
        cs.poly_gate(
            vec![(t_mybignum.value_alloc[1], -F::one()),
                 (bignum.value_alloc[1], F::one()),
                 (p_mybigint.value_alloc[1], self.minus_P_mod_T.value[0]),
                 (p_mybigint.value_alloc[0], self.minus_P_mod_T.value[1])],
            F::zero(), p_sub_c_mybigint.value[1]);

        let tmp0: F = p_mybigint.value[1] * self.minus_P_mod_T.value[1]
            + p_mybigint.value[2] * self.minus_P_mod_T.value[0];
        let tmp0_alloc = cs.alloc(tmp0);
        cs.poly_gate(
            vec![(tmp0_alloc, -F::one()),
                 (p_mybigint.value_alloc[1], self.minus_P_mod_T.value[1]),
                 (p_mybigint.value_alloc[2], self.minus_P_mod_T.value[0]),
                 (Composer::<F>::null(), F::zero())],
            F::zero(), F::zero());
        cs.poly_gate(
            vec![(t_mybignum.value_alloc[2], -F::one()),
                 (bignum.value_alloc[2], F::one()),
                 (p_mybigint.value_alloc[0], self.minus_P_mod_T.value[2]),
                 (tmp0_alloc, F::one())],
            F::zero(), p_sub_c_mybigint.value[2]);

        let tmp1: F = p_mybigint.value[1] * self.minus_P_mod_T.value[2]
            + p_mybigint.value[2] * self.minus_P_mod_T.value[1]
            + p_mybigint.value[3] * self.minus_P_mod_T.value[0];
        let tmp1_alloc = cs.alloc(tmp1);
        cs.poly_gate(
            vec![(tmp1_alloc, -F::one()),
                 (p_mybigint.value_alloc[1], self.minus_P_mod_T.value[2]),
                 (p_mybigint.value_alloc[2], self.minus_P_mod_T.value[1]),
                 (p_mybigint.value_alloc[3], self.minus_P_mod_T.value[0])],
            F::zero(), F::zero());
        cs.poly_gate(
            vec![(t_mybignum.value_alloc[3], -F::one()),
                 (bignum.value_alloc[3], F::one()),
                 (p_mybigint.value_alloc[0], self.minus_P_mod_T.value[3]),
                 (tmp1_alloc, F::one())],
            F::zero(), p_sub_c_mybigint.value[3]);


        self.is_divide_exactly_by_T(&t_mybignum, cs)

    }

    ///ci = ai + bi
    pub fn a_add_b_mod_P(&self, a_bignum: &BigNum<F>, b_bignum: &BigNum<F>, cs: &mut Composer<F>) -> Result<BigNum<F>, Error> {
        let mut a_new = a_bignum.clone();
        let mut b_new = b_bignum.clone();
        for i in 0..4 {
            if a_bignum.value_bits[i] > 253 {
                println!("add:a need mod P");
                a_new = self.a_mod_P(a_bignum, cs)?;
            }
            break;
        }
        for i in 0..4 {
            if b_bignum.value_bits[i] > 253 {
                println!("add:b need mod P");
                b_new = self.a_mod_P(b_bignum, cs)?;
            }
            break;
        }

        let mut cvalue = [F::zero(); 4];
        // let c = a_bignum.value.iter().zip(b_bignum.value.iter()).map(|(&a, &b)| {a + b}).collect();
        for i in 0..4 {
            cvalue[i] = a_new.value[i] + b_new.value[i];
        }

        let mut res = BigNum::<F>::new(&cvalue);
        res.alloc_vars(cs);
        res.set_value_bits(
            [cmp::max(a_new.value_bits[0], b_new.value_bits[0])+1,
                cmp::max(a_new.value_bits[1], b_new.value_bits[1])+1,
                cmp::max(a_new.value_bits[2], b_new.value_bits[2])+1,
                cmp::max(a_new.value_bits[3], b_new.value_bits[3])+1,]
        );
        res.set_num_bits(BigNum::<F>::cal_num_bits(&res.value_bits));

        cs.add_gate(a_new.value_alloc[0], b_new.value_alloc[0], res.value_alloc[0]);
        cs.add_gate(a_new.value_alloc[1], b_new.value_alloc[1], res.value_alloc[1]);
        cs.add_gate(a_new.value_alloc[2], b_new.value_alloc[2], res.value_alloc[2]);
        cs.add_gate(a_new.value_alloc[3], b_new.value_alloc[3], res.value_alloc[3]);

        Ok(res)
    }

    ///Bn = max(
    ///     max(B_a3 - B_P^3, ... , B_a0 - B_P^3) + 1,
    ///     0
    /// )
    pub fn minus_a_mod_P(&self, a_bignum: &BigNum<F>, cs: &mut Composer<F>) -> Result<BigNum<F>, Error> {
        let mut a_new = a_bignum.clone();
        for i in 0..4 {
            // (bn+68 < 254 && Bai-Bpi < bn-1) => Bai < 235
            if a_bignum.value_bits[i] > 235 {
                a_new = self.a_mod_P(a_bignum, cs)?;
            }
            break;
        }

        let tmp1 = cmp::max(a_new.value[0].into_repr().num_bits() as i64 - self.P.value[0].into_repr().num_bits() as i64,
                            a_new.value[1].into_repr().num_bits() as i64 - self.P.value[1].into_repr().num_bits() as i64);
        let tmp2 = cmp::max(a_new.value[2].into_repr().num_bits() as i64 - self.P.value[2].into_repr().num_bits() as i64,
                            a_new.value[3].into_repr().num_bits() as i64 - self.P.value[3].into_repr().num_bits() as i64);
        let tmp3 = cmp::max(tmp1, tmp2) + 1;
        let bn = cmp::max(tmp3, 0) as u32;

        //可能溢出(bn至少是68-50=18)
        // let n = F::from(1u128 << (bn as u128));
        println!("bn {}", bn);
        assert!(bn < 186);

        let mut nslice = Vec::new();
        for _ in 0..(bn / 32) {
            nslice.push(0u32);
        }
        nslice.push(1u32 << (bn%32));
        let bign = BigUint::from_slice(&nslice);
        let n = F::from(bign);

        for i in 0..4 {
            assert!(n*self.P.value[i] > a_new.value[i]);
        }

        let mut res = BigNum::<F>::new(&[
            n*self.P.value[0] - a_new.value[0],
            n*self.P.value[1] - a_new.value[1],
            n*self.P.value[2] - a_new.value[2],
            n*self.P.value[3] - a_new.value[3]
        ]);
        res.alloc_vars(cs);
        res.set_value_bits(
            [ bn + self.P.value_bits[0],
                bn + self.P.value_bits[1],
                bn + self.P.value_bits[2],
                bn + self.P.value_bits[3],]
        );
        res.set_num_bits(BigNum::<F>::cal_num_bits(&res.value_bits));

        // let n_alloc = cs.alloc(n);
        for i in 0..4 {
            cs.poly_gate(
                vec![(res.value_alloc[i], -F::one()),
                     (a_new.value_alloc[i], -F::one()),
                     (Composer::<F>::null(), F::zero()),
                     (Composer::<F>::null(), F::zero())],
                F::zero(), n*self.P.value[i]
            );
        }

        Ok(res)
    }

    ///c = a + (-b)
    pub fn a_sub_b_mod_P(&self, a_bignum: &BigNum<F>, b_bignum: &BigNum<F>, cs: &mut Composer<F>) -> Result<BigNum<F>, Error> {
        let minus_b = self.minus_a_mod_P(b_bignum, cs)?;

        self.a_add_b_mod_P(a_bignum, &minus_b, cs)
    }

    ///both are var BigNum
    pub fn is_eq_mod_P(&self, a_bignum: &BigNum<F>, b_bignum: &BigNum<F>, cs: &mut Composer<F>) -> Result<bool, Error> {

        //不能这么比，这只是最大范围（但理论上应该不影响
        return if a_bignum.num_bits > b_bignum.num_bits {
            let c = self.a_sub_b_mod_P(a_bignum, b_bignum, cs)?;
            self.is_eq_c_mod_P(&c, FP::zero(), cs)
        } else {
            let c = self.a_sub_b_mod_P(b_bignum, a_bignum, cs)?;
            self.is_eq_c_mod_P(&c, FP::zero(), cs)
        }

    }

    ///when B_ai > 253, need this
    pub fn a_mod_P(&self, a_bignum: &BigNum<F>, cs: &mut Composer<F>) -> Result<BigNum<F>, Error> {
        println!("MOD MOD MOD");
        let a_mod_p = a_bignum.value_biguint.clone() % self.P.value_biguint.clone();

        let mut res = BigNum::<F>::new_from_biguint(&a_mod_p);
        res.alloc_vars(cs);

        cs.enforce_range(res.value_alloc[0], (res.value_bits[0] as usize +1) / 2 * 2)?;
        cs.enforce_range(res.value_alloc[1], (res.value_bits[1] as usize +1) / 2 * 2)?;
        cs.enforce_range(res.value_alloc[2], (res.value_bits[2] as usize +1) / 2 * 2)?;
        cs.enforce_range(res.value_alloc[3], (res.value_bits[3] as usize +1) / 2 * 2)?;

        assert!(
            self.is_eq_mod_P(a_bignum, &res, cs)?
        );

        Ok(res)
    }

    //
    pub fn a_mul_b_mod_P(&self, a_bignum: &BigNum<F>, b_bignum: &BigNum<F>, cs: &mut Composer<F>) -> Result<BigNum<F>, Error> {
        let mut a_new = a_bignum.clone();
        let mut b_new = b_bignum.clone();
        let mut Bai_max = 0;
        let mut Bbi_max = 0;
        let mut a_is_modded = false;
        let mut b_is_modded = false;
        for i in 0..4 {
            Bai_max = cmp::max(Bai_max, a_bignum.value_bits[i]);
            Bbi_max = cmp::max(Bbi_max, b_bignum.value_bits[i]);
        }
        println!("a b max {} {}", Bai_max, Bbi_max);

        if Bai_max + 68 + 68 >= 254 {
            println!("1mul:a need to mod P");
            a_new = self.a_mod_P(a_bignum, cs)?;
            a_is_modded = true;
        }
        if Bbi_max + 68 + 68 >= 254 {
            println!("1mul:b need to mod P");
            b_new = self.a_mod_P(b_bignum, cs)?;
            b_is_modded = true;
        }
        if Bai_max + Bbi_max + 68 >= 254 {
            if Bai_max > Bbi_max{
                if !a_is_modded {
                    println!("2mul:a need to mod P");
                    a_new = self.a_mod_P(a_bignum, cs)?;
                }
            } else {
                if !b_is_modded {
                    println!("2mul:b need to mod P");
                    b_new = self.a_mod_P(b_bignum, cs)?;
                }
            }
        }


        //c=A4, d=A5, e=A6
        //t0 = a3*b3*e0 + (a3*b2 + a2*b3)*d0 + (a3*b1 + a2*b2 +a1*b3)*c0 + a0*b0
        //t1 = a3*b3*e1 + (a3*b2 + a2*b3)*d1 + (a3*b1 + a2*b2 +a1*b3)*c1 + a1*b0 + a0*b1
        //t2 = a3*b3*e2 + (a3*b2 + a2*b3)*d2 + (a3*b1 + a2*b2 +a1*b3)*c2 + a2*b0 + a1*b1 + a0*b2
        //t3 = a3*b3*e3 + (a3*b2 + a2*b3)*d3 + (a3*b1 + a2*b2 +a1*b3)*c3 + a3*b0 + a2*b1 + a1*b2 + a0*b3
        let t0 = a_new.value[3] * b_new.value[3] * self.A6_mod_P.value[0]
            + (a_new.value[3] * b_new.value[2] + a_new.value[2] * b_new.value[3]) * self.A5_mod_P.value[0]
            + (a_new.value[3] * b_new.value[1] + a_new.value[2] * b_new.value[2] + a_new.value[1] * b_new.value[3]) * self.A4_mod_P.value[0]
            + a_new.value[0] * b_new.value[0];

        let t1 = a_new.value[3] * b_new.value[3] * self.A6_mod_P.value[1]
            + (a_new.value[3] * b_new.value[2] + a_new.value[2] * b_new.value[3]) * self.A5_mod_P.value[1]
            + (a_new.value[3] * b_new.value[1] + a_new.value[2] * b_new.value[2] + a_new.value[1] * b_new.value[3]) * self.A4_mod_P.value[1]
            + a_new.value[1] * b_new.value[0]
            + a_new.value[0] * b_new.value[1];

        let t2 = a_new.value[3] * b_new.value[3] * self.A6_mod_P.value[2]
            + (a_new.value[3] * b_new.value[2] + a_new.value[2] * b_new.value[3]) * self.A5_mod_P.value[2]
            + (a_new.value[3] * b_new.value[1] + a_new.value[2] * b_new.value[2] + a_new.value[1] * b_new.value[3]) * self.A4_mod_P.value[2]
            + a_new.value[2] * b_new.value[0]
            + a_new.value[1] * b_new.value[1]
            + a_new.value[0] * b_new.value[2];

        let t3 = a_new.value[3] * b_new.value[3] * self.A6_mod_P.value[3]
            + (a_new.value[3] * b_new.value[2] + a_new.value[2] * b_new.value[3]) * self.A5_mod_P.value[3]
            + (a_new.value[3] * b_new.value[1] + a_new.value[2] * b_new.value[2] + a_new.value[1] * b_new.value[3]) * self.A4_mod_P.value[3]
            + a_new.value[3] * b_new.value[0]
            + a_new.value[2] * b_new.value[1]
            + a_new.value[1] * b_new.value[2]
            + a_new.value[0] * b_new.value[3];
        // println!("t0123 {} {} {} {}", t0, t1, t2, t3);

        let mut t_mybignum = BigNum::<F>::new(&[t0,t1,t2,t3]);
        t_mybignum.alloc_vars(cs);
        t_mybignum.set_value_bits(
            [
                BigNum::<F>::cal_add_bits(&[a_new.value_bits[3]+b_new.value_bits[3]+self.A6_mod_P.value_bits[0],
                    a_new.value_bits[3]+b_new.value_bits[2]+self.A5_mod_P.value_bits[0],
                    a_new.value_bits[2]+b_new.value_bits[3]+self.A5_mod_P.value_bits[0],
                    a_new.value_bits[3]+b_new.value_bits[1]+self.A4_mod_P.value_bits[0],
                    a_new.value_bits[2]+b_new.value_bits[2]+self.A4_mod_P.value_bits[0],
                    a_new.value_bits[1]+b_new.value_bits[3]+self.A4_mod_P.value_bits[0],
                    a_new.value_bits[0]+b_new.value_bits[0],
                ]),
                BigNum::<F>::cal_add_bits(&[a_new.value_bits[3]+b_new.value_bits[3]+self.A6_mod_P.value_bits[1],
                    a_new.value_bits[3]+b_new.value_bits[2]+self.A5_mod_P.value_bits[1],
                    a_new.value_bits[2]+b_new.value_bits[3]+self.A5_mod_P.value_bits[1],
                    a_new.value_bits[3]+b_new.value_bits[1]+self.A4_mod_P.value_bits[1],
                    a_new.value_bits[2]+b_new.value_bits[2]+self.A4_mod_P.value_bits[1],
                    a_new.value_bits[1]+b_new.value_bits[3]+self.A4_mod_P.value_bits[1],
                    a_new.value_bits[1]+b_new.value_bits[0],
                    a_new.value_bits[0]+b_new.value_bits[1],
                ]),
                BigNum::<F>::cal_add_bits(&[a_new.value_bits[3]+b_new.value_bits[3]+self.A6_mod_P.value_bits[2],
                    a_new.value_bits[3]+b_new.value_bits[2]+self.A5_mod_P.value_bits[2],
                    a_new.value_bits[2]+b_new.value_bits[3]+self.A5_mod_P.value_bits[2],
                    a_new.value_bits[3]+b_new.value_bits[1]+self.A4_mod_P.value_bits[2],
                    a_new.value_bits[2]+b_new.value_bits[2]+self.A4_mod_P.value_bits[2],
                    a_new.value_bits[1]+b_new.value_bits[3]+self.A4_mod_P.value_bits[2],
                    a_new.value_bits[2]+b_new.value_bits[0],
                    a_new.value_bits[1]+b_new.value_bits[1],
                    a_new.value_bits[0]+b_new.value_bits[2],
                ]),
                BigNum::<F>::cal_add_bits(&[a_new.value_bits[3]+b_new.value_bits[3]+self.A6_mod_P.value_bits[3],
                    a_new.value_bits[3]+b_new.value_bits[2]+self.A5_mod_P.value_bits[3],
                    a_new.value_bits[2]+b_new.value_bits[3]+self.A5_mod_P.value_bits[3],
                    a_new.value_bits[3]+b_new.value_bits[1]+self.A4_mod_P.value_bits[3],
                    a_new.value_bits[2]+b_new.value_bits[2]+self.A4_mod_P.value_bits[3],
                    a_new.value_bits[1]+b_new.value_bits[3]+self.A4_mod_P.value_bits[3],
                    a_new.value_bits[3]+b_new.value_bits[0],
                    a_new.value_bits[2]+b_new.value_bits[1],
                    a_new.value_bits[1]+b_new.value_bits[2],
                    a_new.value_bits[0]+b_new.value_bits[3],
                ]),
            ]
        );
        t_mybignum.set_num_bits(BigNum::<F>::cal_num_bits(&t_mybignum.value_bits));

        //约束ti
        let t0tmp0: F = a_new.value[3] * b_new.value[3] * self.A6_mod_P.value[0];
        let t0tmp1: F = t0tmp0 + a_new.value[3] * b_new.value[2] * self.A5_mod_P.value[0];
        let t0tmp2: F = t0tmp1 + a_new.value[2] * b_new.value[3] * self.A5_mod_P.value[0];
        let t0tmp3: F = t0tmp2 + a_new.value[3] * b_new.value[1] * self.A4_mod_P.value[0];
        let t0tmp4: F = t0tmp3 + a_new.value[2] * b_new.value[2] * self.A4_mod_P.value[0];
        let t0tmp5: F = t0tmp4 + a_new.value[1] * b_new.value[3] * self.A4_mod_P.value[0];
        // let t0tmp6: F = t0tmp5 + a_new.value[0] * b_new.value[0];
        // println!("t0tmp6 {}", t0tmp6);
        let t0tmp0_alloc = cs.alloc(t0tmp0);
        let t0tmp1_alloc = cs.alloc(t0tmp1);
        let t0tmp2_alloc = cs.alloc(t0tmp2);
        let t0tmp3_alloc = cs.alloc(t0tmp3);
        let t0tmp4_alloc = cs.alloc(t0tmp4);
        let t0tmp5_alloc = cs.alloc(t0tmp5);
        // let t0tmp6_alloc = cs.alloc(t0tmp6);

        let t1tmp0: F = a_new.value[3] * b_new.value[3] * self.A6_mod_P.value[1];
        let t1tmp1: F = t1tmp0 + a_new.value[3] * b_new.value[2] * self.A5_mod_P.value[1];
        let t1tmp2: F = t1tmp1 + a_new.value[2] * b_new.value[3] * self.A5_mod_P.value[1];
        let t1tmp3: F = t1tmp2 + a_new.value[3] * b_new.value[1] * self.A4_mod_P.value[1];
        let t1tmp4: F = t1tmp3 + a_new.value[2] * b_new.value[2] * self.A4_mod_P.value[1];
        let t1tmp5: F = t1tmp4 + a_new.value[1] * b_new.value[3] * self.A4_mod_P.value[1];
        let t1tmp6: F = t1tmp5 + a_new.value[1] * b_new.value[0];
        // let t1tmp7: F = t1tmp6 + a_new.value[0] * b_new.value[1];
        // println!("t1tmp7 {}", t1tmp7);
        let t1tmp0_alloc = cs.alloc(t1tmp0);
        let t1tmp1_alloc = cs.alloc(t1tmp1);
        let t1tmp2_alloc = cs.alloc(t1tmp2);
        let t1tmp3_alloc = cs.alloc(t1tmp3);
        let t1tmp4_alloc = cs.alloc(t1tmp4);
        let t1tmp5_alloc = cs.alloc(t1tmp5);
        let t1tmp6_alloc = cs.alloc(t1tmp6);
        // let t1tmp7_alloc = cs.alloc(t1tmp7);

        let t2tmp0: F = a_new.value[3] * b_new.value[3] * self.A6_mod_P.value[2];
        let t2tmp1: F = t2tmp0 + a_new.value[3] * b_new.value[2] * self.A5_mod_P.value[2];
        let t2tmp2: F = t2tmp1 + a_new.value[2] * b_new.value[3] * self.A5_mod_P.value[2];
        let t2tmp3: F = t2tmp2 + a_new.value[3] * b_new.value[1] * self.A4_mod_P.value[2];
        let t2tmp4: F = t2tmp3 + a_new.value[2] * b_new.value[2] * self.A4_mod_P.value[2];
        let t2tmp5: F = t2tmp4 + a_new.value[1] * b_new.value[3] * self.A4_mod_P.value[2];
        let t2tmp6: F = t2tmp5 + a_new.value[2] * b_new.value[0];
        let t2tmp7: F = t2tmp6 + a_new.value[1] * b_new.value[1];
        // let t2tmp8: F = t2tmp7 + a_new.value[0] * b_new.value[2];
        // println!("t2tmp8 {}", t2tmp8);
        let t2tmp0_alloc = cs.alloc(t2tmp0);
        let t2tmp1_alloc = cs.alloc(t2tmp1);
        let t2tmp2_alloc = cs.alloc(t2tmp2);
        let t2tmp3_alloc = cs.alloc(t2tmp3);
        let t2tmp4_alloc = cs.alloc(t2tmp4);
        let t2tmp5_alloc = cs.alloc(t2tmp5);
        let t2tmp6_alloc = cs.alloc(t2tmp6);
        let t2tmp7_alloc = cs.alloc(t2tmp7);
        // let t2tmp8_alloc = cs.alloc(t2tmp8);

        let t3tmp0: F = a_new.value[3] * b_new.value[3] * self.A6_mod_P.value[3];
        let t3tmp1: F = t3tmp0 + a_new.value[3] * b_new.value[2] * self.A5_mod_P.value[3];
        let t3tmp2: F = t3tmp1 + a_new.value[2] * b_new.value[3] * self.A5_mod_P.value[3];
        let t3tmp3: F = t3tmp2 + a_new.value[3] * b_new.value[1] * self.A4_mod_P.value[3];
        let t3tmp4: F = t3tmp3 + a_new.value[2] * b_new.value[2] * self.A4_mod_P.value[3];
        let t3tmp5: F = t3tmp4 + a_new.value[1] * b_new.value[3] * self.A4_mod_P.value[3];
        let t3tmp6: F = t3tmp5 + a_new.value[3] * b_new.value[0];
        let t3tmp7: F = t3tmp6 + a_new.value[2] * b_new.value[1];
        let t3tmp8: F = t3tmp7 + a_new.value[1] * b_new.value[2];
        // let t3tmp9: F = t3tmp8 + a_new.value[0] * b_new.value[3];
        // println!("t3tmp9 {}", t3tmp9);
        let t3tmp0_alloc = cs.alloc(t3tmp0);
        let t3tmp1_alloc = cs.alloc(t3tmp1);
        let t3tmp2_alloc = cs.alloc(t3tmp2);
        let t3tmp3_alloc = cs.alloc(t3tmp3);
        let t3tmp4_alloc = cs.alloc(t3tmp4);
        let t3tmp5_alloc = cs.alloc(t3tmp5);
        let t3tmp6_alloc = cs.alloc(t3tmp6);
        let t3tmp7_alloc = cs.alloc(t3tmp7);
        let t3tmp8_alloc = cs.alloc(t3tmp8);
        // let t3tmp9_alloc = cs.alloc(t3tmp9);

        cs.poly_gate(
            vec![(a_new.value_alloc[3], F::zero()),
                 (b_new.value_alloc[3], F::zero()),
                 (Composer::<F>::null(), F::zero()),
                 (t0tmp0_alloc, -F::one())],
            self.A6_mod_P.value[0], F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[3], F::zero()),
                 (b_new.value_alloc[2], F::zero()),
                 (t0tmp0_alloc, F::one()),
                 (t0tmp1_alloc, -F::one())],
            self.A5_mod_P.value[0], F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[2], F::zero()),
                 (b_new.value_alloc[3], F::zero()),
                 (t0tmp1_alloc, F::one()),
                 (t0tmp2_alloc, -F::one())],
            self.A5_mod_P.value[0], F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[3], F::zero()),
                 (b_new.value_alloc[1], F::zero()),
                 (t0tmp2_alloc, F::one()),
                 (t0tmp3_alloc, -F::one())],
            self.A4_mod_P.value[0], F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[2], F::zero()),
                 (b_new.value_alloc[2], F::zero()),
                 (t0tmp3_alloc, F::one()),
                 (t0tmp4_alloc, -F::one())],
            self.A4_mod_P.value[0], F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[1], F::zero()),
                 (b_new.value_alloc[3], F::zero()),
                 (t0tmp4_alloc, F::one()),
                 (t0tmp5_alloc, -F::one())],
            self.A4_mod_P.value[0], F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[0], F::zero()),
                 (b_new.value_alloc[0], F::zero()),
                 (t0tmp5_alloc, F::one()),
                 (t_mybignum.value_alloc[0], -F::one())],
            F::one(), F::zero());

        cs.poly_gate(
            vec![(a_new.value_alloc[3], F::zero()),
                 (b_new.value_alloc[3], F::zero()),
                 (Composer::<F>::null(), F::zero()),
                 (t1tmp0_alloc, -F::one())],
            self.A6_mod_P.value[1], F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[3], F::zero()),
                 (b_new.value_alloc[2], F::zero()),
                 (t1tmp0_alloc, F::one()),
                 (t1tmp1_alloc, -F::one())],
            self.A5_mod_P.value[1], F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[2], F::zero()),
                 (b_new.value_alloc[3], F::zero()),
                 (t1tmp1_alloc, F::one()),
                 (t1tmp2_alloc, -F::one())],
            self.A5_mod_P.value[1], F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[3], F::zero()),
                 (b_new.value_alloc[1], F::zero()),
                 (t1tmp2_alloc, F::one()),
                 (t1tmp3_alloc, -F::one())],
            self.A4_mod_P.value[1], F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[2], F::zero()),
                 (b_new.value_alloc[2], F::zero()),
                 (t1tmp3_alloc, F::one()),
                 (t1tmp4_alloc, -F::one())],
            self.A4_mod_P.value[1], F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[1], F::zero()),
                 (b_new.value_alloc[3], F::zero()),
                 (t1tmp4_alloc, F::one()),
                 (t1tmp5_alloc, -F::one())],
            self.A4_mod_P.value[1], F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[1], F::zero()),
                 (b_new.value_alloc[0], F::zero()),
                 (t1tmp5_alloc, F::one()),
                 (t1tmp6_alloc, -F::one())],
            F::one(), F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[0], F::zero()),
                 (b_new.value_alloc[1], F::zero()),
                 (t1tmp6_alloc, F::one()),
                 (t_mybignum.value_alloc[1], -F::one())],
            F::one(), F::zero());

        cs.poly_gate(
            vec![(a_new.value_alloc[3], F::zero()),
                 (b_new.value_alloc[3], F::zero()),
                 (Composer::<F>::null(), F::zero()),
                 (t2tmp0_alloc, -F::one())],
            self.A6_mod_P.value[2], F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[3], F::zero()),
                 (b_new.value_alloc[2], F::zero()),
                 (t2tmp0_alloc, F::one()),
                 (t2tmp1_alloc, -F::one())],
            self.A5_mod_P.value[2], F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[2], F::zero()),
                 (b_new.value_alloc[3], F::zero()),
                 (t2tmp1_alloc, F::one()),
                 (t2tmp2_alloc, -F::one())],
            self.A5_mod_P.value[2], F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[3], F::zero()),
                 (b_new.value_alloc[1], F::zero()),
                 (t2tmp2_alloc, F::one()),
                 (t2tmp3_alloc, -F::one())],
            self.A4_mod_P.value[2], F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[2], F::zero()),
                 (b_new.value_alloc[2], F::zero()),
                 (t2tmp3_alloc, F::one()),
                 (t2tmp4_alloc, -F::one())],
            self.A4_mod_P.value[2], F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[1], F::zero()),
                 (b_new.value_alloc[3], F::zero()),
                 (t2tmp4_alloc, F::one()),
                 (t2tmp5_alloc, -F::one())],
            self.A4_mod_P.value[2], F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[2], F::zero()),
                 (b_new.value_alloc[0], F::zero()),
                 (t2tmp5_alloc, F::one()),
                 (t2tmp6_alloc, -F::one())],
            F::one(), F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[1], F::zero()),
                 (b_new.value_alloc[1], F::zero()),
                 (t2tmp6_alloc, F::one()),
                 (t2tmp7_alloc, -F::one())],
            F::one(), F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[0], F::zero()),
                 (b_new.value_alloc[2], F::zero()),
                 (t2tmp7_alloc, F::one()),
                 (t_mybignum.value_alloc[2], -F::one())],
            F::one(), F::zero());

        cs.poly_gate(
            vec![(a_new.value_alloc[3], F::zero()),
                 (b_new.value_alloc[3], F::zero()),
                 (Composer::<F>::null(), F::zero()),
                 (t3tmp0_alloc, -F::one())],
            self.A6_mod_P.value[3], F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[3], F::zero()),
                 (b_new.value_alloc[2], F::zero()),
                 (t3tmp0_alloc, F::one()),
                 (t3tmp1_alloc, -F::one())],
            self.A5_mod_P.value[3], F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[2], F::zero()),
                 (b_new.value_alloc[3], F::zero()),
                 (t3tmp1_alloc, F::one()),
                 (t3tmp2_alloc, -F::one())],
            self.A5_mod_P.value[3], F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[3], F::zero()),
                 (b_new.value_alloc[1], F::zero()),
                 (t3tmp2_alloc, F::one()),
                 (t3tmp3_alloc, -F::one())],
            self.A4_mod_P.value[3], F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[2], F::zero()),
                 (b_new.value_alloc[2], F::zero()),
                 (t3tmp3_alloc, F::one()),
                 (t3tmp4_alloc, -F::one())],
            self.A4_mod_P.value[3], F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[1], F::zero()),
                 (b_new.value_alloc[3], F::zero()),
                 (t3tmp4_alloc, F::one()),
                 (t3tmp5_alloc, -F::one())],
            self.A4_mod_P.value[3], F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[3], F::zero()),
                 (b_new.value_alloc[0], F::zero()),
                 (t3tmp5_alloc, F::one()),
                 (t3tmp6_alloc, -F::one())],
            F::one(), F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[2], F::zero()),
                 (b_new.value_alloc[1], F::zero()),
                 (t3tmp6_alloc, F::one()),
                 (t3tmp7_alloc, -F::one())],
            F::one(), F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[1], F::zero()),
                 (b_new.value_alloc[2], F::zero()),
                 (t3tmp7_alloc, F::one()),
                 (t3tmp8_alloc, -F::one())],
            F::one(), F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[0], F::zero()),
                 (b_new.value_alloc[3], F::zero()),
                 (t3tmp8_alloc, F::one()),
                 (t_mybignum.value_alloc[3], -F::one())],
            F::one(), F::zero());

        Ok(t_mybignum)
    }

    ///c = a*b mod P (abc are var BigNum
    pub fn is_a_mul_b_eq_c_mod_P(&self, a_bignum: &BigNum<F>, b_bignum: &BigNum<F>, c_bignum: &BigNum<F>, cs: &mut Composer<F>) -> Result<bool, Error> {
        let mut a_new = a_bignum.clone();
        let mut b_new = b_bignum.clone();

        let mut Bai_max = 0;
        let mut Bbi_max = 0;
        let mut a_is_modded = false;
        let mut b_is_modded = false;
        for i in 0..4 {
            Bai_max = cmp::max(Bai_max, a_bignum.value_bits[i]);
            Bbi_max = cmp::max(Bbi_max, b_bignum.value_bits[i]);
        }

        if Bai_max + 68 >= 254 {
            a_new = self.a_mod_P(a_bignum, cs)?;
            a_is_modded = true;
        }
        if Bbi_max + 68 >= 254 {
            b_new = self.a_mod_P(b_bignum, cs)?;
            b_is_modded = true;
        }
        if Bai_max + Bbi_max >= 254 {
            if Bai_max > Bbi_max{
                if !a_is_modded {
                    a_new = self.a_mod_P(a_bignum, cs)?;
                }
            } else {
                if !b_is_modded {
                    b_new = self.a_mod_P(b_bignum, cs)?;
                }
            }
        }

        let fp = FP::Params::MODULUS;


        let minus_c = self.minus_a_mod_P(c_bignum, cs)?;

        //calculate p (where a + c' -p*P=0)
        if !((a_new.value_biguint.clone()*b_new.value_biguint.clone() + minus_c.value_biguint.clone()) % fp.into()).is_zero() { return Ok(false); }
        let p: BigUint = (a_new.value_biguint.clone()*b_new.value_biguint.clone() + minus_c.value_biguint.clone()) / fp.into();
        println!("pp {}", p);
        println!("pbits {}", p.bits());
        let mut p_mybigint = BigNum::<F>::new_from_biguint(&p);
        p_mybigint.alloc_vars(cs);

        let mut B_p = 0;
        if cmp::max(a_new.num_bits + b_new.num_bits, minus_c.num_bits as u32) + 1 >= 254 {
            B_p = cmp::max(a_new.num_bits + b_new.num_bits, minus_c.num_bits as u32) + 1 - 254
        }
        let mut B_pi = [0u32, 0u32, 0u32, 0u32];
        //todo 当某些Bpi为0时可以优化，不增加pi（上面某节也是）
        if B_p < 68 {
            B_pi[0] = B_p;
        } else {
            B_pi[0] = 68;
            let tmp = B_p - 68;
            if tmp < 68 {
                B_pi[1] = tmp;
            } else {
                B_pi[1] = 68;
                let tmp = tmp - 68;
                if tmp < 68 {
                    B_pi[2] = tmp;
                } else {
                    B_pi[2] = 68;
                    B_pi[3] = tmp - 68;
                }
            }
        }

        p_mybigint.set_value_bits(B_pi);
        p_mybigint.set_num_bits(BigNum::<F>::cal_num_bits(&p_mybigint.value_bits));

        cs.enforce_range(p_mybigint.value_alloc[0], (p_mybigint.value_bits[0] as usize +1) / 2 * 2)?;
        cs.enforce_range(p_mybigint.value_alloc[1], (p_mybigint.value_bits[1] as usize +1) / 2 * 2)?;
        cs.enforce_range(p_mybigint.value_alloc[2], (p_mybigint.value_bits[2] as usize +1) / 2 * 2)?;
        cs.enforce_range(p_mybigint.value_alloc[3], (p_mybigint.value_bits[3] as usize +1) / 2 * 2)?;

        let ar: BigUint = a_new.value_biguint.clone() % F::Params::MODULUS.into();
        let br: BigUint = b_new.value_biguint.clone() % F::Params::MODULUS.into();
        let mcr: BigUint = minus_c.value_biguint.clone() % F::Params::MODULUS.into();
        let pr: BigUint = p.clone() % F::Params::MODULUS.into();

        let ar = F::from(ar);
        let br = F::from(br);
        let mcr = F::from(mcr);
        let pr = F::from(pr);

        assert!((ar*br + mcr - pr * self.P_mod_R).is_zero());

        let tmp0: F = a_new.value[0] + a_new.value[1]*self.A + a_new.value[2]*self.A2;
        let tmp1: F = tmp0 + a_new.value[3]*self.A3;

        let tmp2: F = b_new.value[0] + b_new.value[1]*self.A + b_new.value[2]*self.A2;
        let tmp3: F = tmp2 + b_new.value[3]*self.A3;

        let tmp4: F = tmp1 * tmp3 + minus_c.value[0];
        let tmp5: F = tmp4 + minus_c.value[1]*self.A + minus_c.value[2]*self.A2;
        let tmp6: F = tmp5 + minus_c.value[3]*self.A3 - p_mybigint.value[0]*self.P_mod_R;
        let tmp7: F = tmp6 - p_mybigint.value[1]*self.A*self.P_mod_R - p_mybigint.value[2]*self.A2*self.P_mod_R;
        // let tmp8: F = tmp7 - p_mybigint.value[3]*self.A3*self.P_mod_R;
        // println!("tmp8 {}", tmp8);

        let tmp0_alloc = cs.alloc(tmp0);
        let tmp1_alloc = cs.alloc(tmp1);
        let tmp2_alloc = cs.alloc(tmp2);
        let tmp3_alloc = cs.alloc(tmp3);
        let tmp4_alloc = cs.alloc(tmp4);
        let tmp5_alloc = cs.alloc(tmp5);
        let tmp6_alloc = cs.alloc(tmp6);
        let tmp7_alloc = cs.alloc(tmp7);
        // let tmp8_alloc = cs.alloc(tmp8);

        cs.poly_gate(
            vec![(a_new.value_alloc[0], F::one()),
                 (a_new.value_alloc[1], self.A),
                 (a_new.value_alloc[2], self.A2),
                 (tmp0_alloc, -F::one())],
            F::zero(), F::zero());
        cs.poly_gate(
            vec![(tmp0_alloc, F::one()),
                 (a_new.value_alloc[3], self.A3),
                 (tmp1_alloc, -F::one())],
            F::zero(), F::zero());
        cs.poly_gate(
            vec![(b_new.value_alloc[0], F::one()),
                 (b_new.value_alloc[1], self.A),
                 (b_new.value_alloc[2], self.A2),
                 (tmp2_alloc, -F::one())],
            F::zero(), F::zero());
        cs.poly_gate(
            vec![(tmp2_alloc, F::one()),
                 (b_new.value_alloc[3], self.A3),
                 (tmp3_alloc, -F::one())],
            F::zero(), F::zero());
        cs.poly_gate(
            vec![(tmp1_alloc, F::zero()),
                 (tmp3_alloc, F::zero()),
                 (minus_c.value_alloc[0], F::one()),
                 (tmp4_alloc, -F::one())],
            F::one(), F::zero());
        cs.poly_gate(
            vec![(tmp4_alloc, F::one()),
                 (minus_c.value_alloc[1], self.A),
                 (minus_c.value_alloc[2], self.A2),
                 (tmp5_alloc, -F::one())],
            F::zero(), F::zero());
        cs.poly_gate(
            vec![(tmp5_alloc, F::one()),
                 (minus_c.value_alloc[3], self.A3),
                 (p_mybigint.value_alloc[0], -self.P_mod_R),
                 (tmp6_alloc, -F::one())],
            F::zero(), F::zero());
        cs.poly_gate(
            vec![(tmp6_alloc, F::one()),
                 (p_mybigint.value_alloc[1], -self.A*self.P_mod_R),
                 (p_mybigint.value_alloc[2], -self.A2*self.P_mod_R),
                 (tmp7_alloc, -F::one())],
            F::zero(), F::zero());
        cs.poly_gate(
            vec![(tmp7_alloc, F::one()),
                 (p_mybigint.value_alloc[3], -self.A3*self.P_mod_R),
                 (Composer::<F>::null(), -F::one())],
            F::zero(), F::zero());


        //calculate t
        //t0 = a0*b0 + mc0 + p0*Pt0
        //t1 = a0*b1 +a1*b0 + mc1 + p0*Pt1 + p1*Pt0
        //t2 = a0*b2 +a1*b1 + a2*b0 + mc2 + p0*Pt2 + p1*Pt1 + p2*Pt0
        //t3 = a0*b3 +a1*b2 + a2*b1 + a3*b0 + mc3 + p0*Pt3 + p1*Pt2 + p2*Pt1 + p3*Pt0
        let t0 = a_new.value[0] * b_new.value[0]
            + minus_c.value[0]
            + p_mybigint.value[0] * self.minus_P_mod_T.value[0];

        let t1 = a_new.value[0] * b_new.value[1]
            + a_new.value[1] * b_new.value[0]
            + minus_c.value[1]
            + p_mybigint.value[0] * self.minus_P_mod_T.value[1]
            + p_mybigint.value[1] * self.minus_P_mod_T.value[0];

        let t2 = a_new.value[0] * b_new.value[2]
            + a_new.value[1] * b_new.value[1]
            + a_new.value[2] * b_new.value[0]
            + minus_c.value[2]
            + p_mybigint.value[0] * self.minus_P_mod_T.value[2]
            + p_mybigint.value[1] * self.minus_P_mod_T.value[1]
            + p_mybigint.value[2] * self.minus_P_mod_T.value[0];

        let t3 = a_new.value[0] * b_new.value[3]
            + a_new.value[1] * b_new.value[2]
            + a_new.value[2] * b_new.value[1]
            + a_new.value[3] * b_new.value[0]
            + minus_c.value[3]
            + p_mybigint.value[0] * self.minus_P_mod_T.value[3]
            + p_mybigint.value[1] * self.minus_P_mod_T.value[2]
            + p_mybigint.value[2] * self.minus_P_mod_T.value[1]
            + p_mybigint.value[3] * self.minus_P_mod_T.value[0];
        // println!("t0123 {} {} {} {}", t0, t1, t2, t3);

        let mut t_mybignum = BigNum::<F>::new(&[t0,t1,t2,t3]);
        t_mybignum.alloc_vars(cs);
        t_mybignum.set_value_bits(
            [
                BigNum::<F>::cal_add_bits(&[a_new.value_bits[0]+b_new.value_bits[0],
                    minus_c.value_bits[0],
                    p_mybigint.value_bits[0] * self.minus_P_mod_T.value_bits[0],
                ]),
                BigNum::<F>::cal_add_bits(&[a_new.value_bits[0]+b_new.value_bits[1],
                    a_new.value_bits[1]+b_new.value_bits[0],
                    minus_c.value_bits[1],
                    p_mybigint.value_bits[0] * self.minus_P_mod_T.value_bits[1],
                    p_mybigint.value_bits[1] * self.minus_P_mod_T.value_bits[0],
                ]),
                BigNum::<F>::cal_add_bits(&[a_new.value_bits[0]+b_new.value_bits[2],
                    a_new.value_bits[1]+b_new.value_bits[1],
                    a_new.value_bits[0]+b_new.value_bits[2],
                    minus_c.value_bits[2],
                    p_mybigint.value_bits[0] * self.minus_P_mod_T.value_bits[2],
                    p_mybigint.value_bits[1] * self.minus_P_mod_T.value_bits[1],
                    p_mybigint.value_bits[2] * self.minus_P_mod_T.value_bits[0],
                ]),
                BigNum::<F>::cal_add_bits(&[a_new.value_bits[0]+b_new.value_bits[3],
                    a_new.value_bits[1]+b_new.value_bits[2],
                    a_new.value_bits[2]+b_new.value_bits[1],
                    a_new.value_bits[3]+b_new.value_bits[0],
                    minus_c.value_bits[3],
                    p_mybigint.value_bits[0] * self.minus_P_mod_T.value_bits[3],
                    p_mybigint.value_bits[1] * self.minus_P_mod_T.value_bits[2],
                    p_mybigint.value_bits[2] * self.minus_P_mod_T.value_bits[1],
                    p_mybigint.value_bits[3] * self.minus_P_mod_T.value_bits[0],
                ]),
            ]
        );
        t_mybignum.set_num_bits(BigNum::<F>::cal_num_bits(&t_mybignum.value_bits));


        //约束ti
        let t0tmp0: F = a_new.value[0] * b_new.value[0] + minus_c.value[0];
        // let t0tmp1: F = t0tmp0 + p_mybigint.value[0] * self.minus_P_mod_T.value[0];
        // println!("t0tmp1 {}", t0tmp1);
        let t1tmp0: F = a_new.value[0] * b_new.value[1] + minus_c.value[1];
        let t1tmp1: F = a_new.value[1] * b_new.value[0] + t1tmp0;
        // let t1tmp2: F = t1tmp1 + p_mybigint.value[0] * self.minus_P_mod_T.value[1]
        //     + p_mybigint.value[1] * self.minus_P_mod_T.value[0];
        // println!("t1tmp2 {}", t1tmp2);
        let t2tmp0: F = a_new.value[0] * b_new.value[2] + minus_c.value[2];
        let t2tmp1: F = a_new.value[1] * b_new.value[1] + t2tmp0;
        let t2tmp2: F = a_new.value[2] * b_new.value[0] + t2tmp1;
        let t2tmp3: F = t2tmp2 + p_mybigint.value[0] * self.minus_P_mod_T.value[2]
            + p_mybigint.value[1] * self.minus_P_mod_T.value[1];
        // let t2tmp4: F = t2tmp3 + p_mybigint.value[2] * self.minus_P_mod_T.value[0];
        // println!("t2tmp4 {}", t2tmp4);
        let t3tmp0: F = a_new.value[0] * b_new.value[3] + minus_c.value[3];
        let t3tmp1: F = a_new.value[1] * b_new.value[2] + t3tmp0;
        let t3tmp2: F = a_new.value[2] * b_new.value[1] + t3tmp1;
        let t3tmp3: F = a_new.value[3] * b_new.value[0] + t3tmp2;
        let t3tmp4: F = t3tmp3 + p_mybigint.value[0] * self.minus_P_mod_T.value[3]
            + p_mybigint.value[1] * self.minus_P_mod_T.value[2];
        // let t3tmp5: F = t3tmp4 + p_mybigint.value[2] * self.minus_P_mod_T.value[1]
        //     + p_mybigint.value[3] * self.minus_P_mod_T.value[0];
        // println!("t3tmp5 {}", t3tmp5);

        let t0tmp0_alloc = cs.alloc(t0tmp0);
        // let t0tmp1_alloc = cs.alloc(t0tmp1);

        let t1tmp0_alloc = cs.alloc(t1tmp0);
        let t1tmp1_alloc = cs.alloc(t1tmp1);
        // let t1tmp2_alloc = cs.alloc(t1tmp2);

        let t2tmp0_alloc = cs.alloc(t2tmp0);
        let t2tmp1_alloc = cs.alloc(t2tmp1);
        let t2tmp2_alloc = cs.alloc(t2tmp2);
        let t2tmp3_alloc = cs.alloc(t2tmp3);
        // let t2tmp4_alloc = cs.alloc(t2tmp4);

        let t3tmp0_alloc = cs.alloc(t3tmp0);
        let t3tmp1_alloc = cs.alloc(t3tmp1);
        let t3tmp2_alloc = cs.alloc(t3tmp2);
        let t3tmp3_alloc = cs.alloc(t3tmp3);
        let t3tmp4_alloc = cs.alloc(t3tmp4);
        // let t3tmp5_alloc = cs.alloc(t3tmp5);


        cs.poly_gate(
            vec![(a_new.value_alloc[0], F::zero()),
                 (b_new.value_alloc[0], F::zero()),
                 (t0tmp0_alloc, -F::one()),
                 (minus_c.value_alloc[0], F::one())],
            F::one(), F::zero());
        cs.poly_gate(
            vec![(t0tmp0_alloc, F::one()),
                 (p_mybigint.value_alloc[0], self.minus_P_mod_T.value[0]),
                 (t_mybignum.value_alloc[0], -F::one())],
            F::zero(), F::zero());

        cs.poly_gate(
            vec![(a_new.value_alloc[0], F::zero()),
                 (b_new.value_alloc[1], F::zero()),
                 (t1tmp0_alloc, -F::one()),
                 (minus_c.value_alloc[1], F::one())],
            F::one(), F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[1], F::zero()),
                 (b_new.value_alloc[0], F::zero()),
                 (t1tmp1_alloc, -F::one()),
                 (t1tmp0_alloc, F::one())],
            F::one(), F::zero());
        cs.poly_gate(
            vec![(t1tmp1_alloc, F::one()),
                 (p_mybigint.value_alloc[0], self.minus_P_mod_T.value[1]),
                 (p_mybigint.value_alloc[1], self.minus_P_mod_T.value[0]),
                 (t_mybignum.value_alloc[1], -F::one())],
            F::zero(), F::zero());

        cs.poly_gate(
            vec![(a_new.value_alloc[0], F::zero()),
                 (b_new.value_alloc[2], F::zero()),
                 (t2tmp0_alloc, -F::one()),
                 (minus_c.value_alloc[2], F::one())],
            F::one(), F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[1], F::zero()),
                 (b_new.value_alloc[1], F::zero()),
                 (t2tmp1_alloc, -F::one()),
                 (t2tmp0_alloc, F::one())],
            F::one(), F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[2], F::zero()),
                 (b_new.value_alloc[0], F::zero()),
                 (t2tmp2_alloc, -F::one()),
                 (t2tmp1_alloc, F::one())],
            F::one(), F::zero());
        cs.poly_gate(
            vec![(t2tmp2_alloc, F::one()),
                 (p_mybigint.value_alloc[0], self.minus_P_mod_T.value[2]),
                 (p_mybigint.value_alloc[1], self.minus_P_mod_T.value[1]),
                 (t2tmp3_alloc, -F::one())],
            F::zero(), F::zero());
        cs.poly_gate(
            vec![(t2tmp3_alloc, F::one()),
                 (p_mybigint.value_alloc[2], self.minus_P_mod_T.value[0]),
                 (t_mybignum.value_alloc[2], -F::one())],
            F::zero(), F::zero());

        cs.poly_gate(
            vec![(a_new.value_alloc[0], F::zero()),
                 (b_new.value_alloc[3], F::zero()),
                 (t3tmp0_alloc, -F::one()),
                 (minus_c.value_alloc[3], F::one())],
            F::one(), F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[1], F::zero()),
                 (b_new.value_alloc[2], F::zero()),
                 (t3tmp1_alloc, -F::one()),
                 (t3tmp0_alloc, F::one())],
            F::one(), F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[2], F::zero()),
                 (b_new.value_alloc[1], F::zero()),
                 (t3tmp2_alloc, -F::one()),
                 (t3tmp1_alloc, F::one())],
            F::one(), F::zero());
        cs.poly_gate(
            vec![(a_new.value_alloc[3], F::zero()),
                 (b_new.value_alloc[0], F::zero()),
                 (t3tmp3_alloc, -F::one()),
                 (t3tmp2_alloc, F::one())],
            F::one(), F::zero());
        cs.poly_gate(
            vec![(t3tmp3_alloc, F::one()),
                 (p_mybigint.value_alloc[0], self.minus_P_mod_T.value[3]),
                 (p_mybigint.value_alloc[1], self.minus_P_mod_T.value[2]),
                 (t3tmp4_alloc, -F::one())],
            F::zero(), F::zero());
        cs.poly_gate(
            vec![(t3tmp4_alloc, F::one()),
                 (p_mybigint.value_alloc[2], self.minus_P_mod_T.value[1]),
                 (p_mybigint.value_alloc[3], self.minus_P_mod_T.value[0]),
                 (t_mybignum.value_alloc[3], -F::one())],
            F::zero(), F::zero());


        assert!(self.is_divide_exactly_by_T(&t_mybignum, cs)?);

        Ok(true)
    }

    ///c = a/b mod P (abc are var BigNum
    pub fn is_a_div_b_eq_c_mod_P(&self, a_bignum: &BigNum<F>, b_bignum: &BigNum<F>, c_bignum: &BigNum<F>, cs: &mut Composer<F>) -> Result<bool, Error> {
        self.is_a_mul_b_eq_c_mod_P(c_bignum, b_bignum, a_bignum, cs)
    }

    ///first get (a - b)
    /// then give c, witch c(a-b) == 1
    pub fn is_a_not_eq_b_mod_P(&self, a_bignum: &BigNum<F>, b_bignum: &BigNum<F>, cs: &mut Composer<F>) -> Result<bool, Error> {
        let fp = FP::Params::MODULUS;

        if a_bignum.value_biguint.clone() % fp.into() == b_bignum.value_biguint.clone() % fp.into() { return Ok(false); }

        let asubb = self.a_sub_b_mod_P(a_bignum, b_bignum, cs)?;

        //todo
        // let test1 = FP::from(1u128);
        // let test3 = FP::from(3u128);
        // let res: FP = test1.div(test3);
        // println!("{}", res);
        // let asubb_modp = self.a_mod_P(&asubb, cs);

        let asubb_modp: BigUint = asubb.value_biguint.clone() % fp.into();
        let asubb_p = FP::from(asubb_modp.clone());

        let one = FP::from(1u128);
        let c: FP = one.div(asubb_p);
        println!("c {}", c);

        let one_mybignum = BigNum::<F>::new(&[F::one(), F::zero(), F::zero(), F::zero()]);
        let mut c_mybignum = BigNum::<F>::new_from_biguint(&c.into_repr().into());
        c_mybignum.alloc_vars(cs);
        // c_mybignum.set_value_bits(); //不需要，但要约束range
        cs.enforce_range(c_mybignum.value_alloc[0], (c_mybignum.value_bits[0] as usize +1) / 2 * 2)?;
        cs.enforce_range(c_mybignum.value_alloc[1], (c_mybignum.value_bits[1] as usize +1) / 2 * 2)?;
        cs.enforce_range(c_mybignum.value_alloc[2], (c_mybignum.value_bits[2] as usize +1) / 2 * 2)?;
        cs.enforce_range(c_mybignum.value_alloc[3], (c_mybignum.value_bits[3] as usize +1) / 2 * 2)?;


        //延迟mod P
        // let c_mul_asubb = self.a_mul_b_mod_P(&c_mybignum, &asubb, cs);

        self.is_a_mul_b_eq_c_mod_P(&c_mybignum, &asubb, &one_mybignum, cs)

    }
}


#[cfg(test)]
mod tests {
    use ark_bn254::Fq;
    use ark_bn254::FqParameters;
    use ark_bn254::Fr;
    use ark_ff::{FpParameters, One, Zero};
    use ark_std::test_rng;
    use num_bigint::BigUint;
    use crate::{Composer, Error};
    use crate::nonnative_field::{BigNum, FpField, NonnativeField};
    use crate::prover::Prover;
    use crate::utils::Field;
    use crate::GeneralEvaluationDomain;

    // use super::*;

    // impl Field for Fr {
    //     fn coset_generator(index: usize) -> Self {
    //         match index % 4 {
    //             1 => Fr::from(7_u64),
    //             2 => Fr::from(13_u64),
    //             3 => Fr::from(17_u64),
    //             _ => Self::one(),
    //         }
    //     }
    // }

    impl FpField for Fq {

    }

    #[test]
    fn test_nnf() -> Result<(), Error>  {
        let mut cs = {
            // x^3 + x + pi = 35
            let mut cs = Composer::new(4);
            let pi = cs.alloc_input(Fr::from(5));
            let x = cs.alloc(Fr::from(3));
            let y = cs.mul(x, x);
            let z = cs.mul(x, y);
            let u = cs.add(x, z);
            let v = cs.add(pi, u);
            cs.enforce_constant(v, Fr::from(35));

            cs
        };

        let nft = NonnativeField::<Fr, Fq>::new(&mut cs);

        let mut testp = BigNum::<Fr>::new_from_biguint(&FqParameters::MODULUS.into());
        testp.alloc_vars(&mut cs);
        //68 67 62 50
        println!("P {:?}", nft.P.value_bits);
        //62 68 68 48
        println!("A4 {:?}", nft.A4_mod_P.value_bits);
        //67 68 67 49
        println!("A5 {:?}", nft.A5_mod_P.value_bits);
        //68 66 68 49
        println!("A6 {:?}", nft.A6_mod_P.value_bits);

        let mut testb = BigNum::<Fr>::new(&[Fr::zero(),Fr::zero(),Fr::zero(),nft.A]);
        testb.alloc_vars(&mut cs);

        let res = nft.is_eq_c_mod_P(&testp.clone(), Fq::zero(), &mut cs)?;
        println!("{}", res);

        let testadd =  nft.a_add_b_mod_P(&testb.clone(), &testp.clone(), &mut cs)?;

        let testminus = nft.minus_a_mod_P(&testadd, &mut cs)?;
        println!("testminus");
        for item in testminus.value {
            println!("{}", item);
        }

        let testsub = nft.a_sub_b_mod_P(&testb, &testminus, &mut cs)?;
        println!("testsub");
        for item in testsub.value {
            println!("{}", item);
        }

        let test2b = nft.a_add_b_mod_P(&testb, &testb, &mut cs)?;
        println!("test2b");
        for item in test2b.value {
            println!("{}", item);
        }
        let testeq = nft.is_eq_mod_P(&testsub, &test2b, &mut cs)?;
        println!("{}", testeq);

        let testmul = nft.a_mul_b_mod_P(&testb, &test2b, &mut cs)?;
        println!("testmul");
        for item in testmul.value {
            println!("{}", item);
        }

        let testmul2 = nft.a_mul_b_mod_P(&testb, &testmul, &mut cs)?;
        println!("testmul2");
        for item in testmul2.value {
            println!("{}", item);
        }

        let mul2modp = nft.a_mod_P(&testmul2, &mut cs)?;
        println!("mul2modp");
        for item in mul2modp.value {
            println!("{}", item);
        }

        let res = nft.is_a_mul_b_eq_c_mod_P(&testb, &testmul, &mul2modp, &mut cs)?;
        println!("{}", res);


        let rng = &mut test_rng();

        let pk = cs.compute_prover_key::<GeneralEvaluationDomain<Fr>>()?;
        let mut prover = Prover::new(&pk);
        prover.prove(&mut cs, rng)?;

        Ok(())
    }

    #[test]
    fn test_unrelated_matters() {
    }

}
