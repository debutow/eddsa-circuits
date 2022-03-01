use std::marker::PhantomData;
use ark_ff::BigInteger;
use crate::utils::Field;
use plonk::Composer;
use plonk::composer::Variable;

pub struct NumToBits<F: Field>{
    _marker: PhantomData<F>,
}

impl<F: Field> NumToBits<F> {
    ///little endian (no need range check
    pub fn num_to_bits(
        cs: &mut Composer<F>,
        num: Variable,
        n: usize,
    ) -> Vec<Variable> {
        //range need an even 'n'
        // cs.enforce_range(num, (n+1) / 2 * 2 );

        let mut out = Vec::new();

        let value_num = cs.get_assignment(num);
        let out_bits = value_num.into_repr().to_bits_le();
        for i in n..out_bits.len() {
            assert_eq!(out_bits[i], false);
        }

        //must write constraint in circuit
        let mut tmp_sum = F::zero();
        let mut var_tmp_old = Composer::<F>::null();
        let mut tmp = F::one();
        for i in 0..n {
            if out_bits[i] {
                out.push(cs.alloc(F::one()));
                tmp_sum += tmp;
            } else {
                out.push(cs.alloc(F::zero()));
            }

            let var_tmp = cs.alloc(tmp_sum);
            //new = old + 'Const' * out[i]
            cs.poly_gate(
                vec![(out[i], tmp), (var_tmp_old, F::one()), (var_tmp, -F::one())],
                F::zero(),
                F::zero()
            );
            var_tmp_old = var_tmp;

            tmp = tmp.double();
        }
        cs.enforce_eq(var_tmp_old, num);

        out
    }

    ///output n selectors, only out[num_value] = 1, others = 0
    #[allow(dead_code)]
    pub fn num_to_selectors(
        cs: &mut Composer<F>,
        num: Variable,
        n: usize,
    ) -> Vec<Variable> {
        //convert to u64
        let num_uint: Vec<u64> = cs.get_assignment(num).into_repr().as_ref().into();
        assert_eq!(num_uint[1], 0);
        assert_eq!(num_uint[2], 0);
        assert_eq!(num_uint[3], 0);
        assert!(num_uint[0] < n as u64);

        let mut out = Vec::new();

        for i in 0..n {
            // let out_i = cs.alloc(F::zero());
            let out_i = Composer::<F>::null();
            cs.enforce_bool(out_i);
            out.push(out_i);

            //(num - CONSTi) * out == 0
            cs.poly_gate(
                vec![(num, F::zero()), (out_i, -F::from(i as u128))],
                F::one(),
                F::zero(),
            );
        }

        //the selector eq to 1
        let out_i = cs.alloc(F::one());
        cs.enforce_bool(out_i);
        out[num_uint[0] as usize] = out_i;

        //(num - CONSTi) * out == 0
        cs.poly_gate(
            vec![(num, F::zero()), (out_i, -F::from(num_uint[0] as u128))],
            F::one(),
            F::zero(),
        );

        out
    }
}