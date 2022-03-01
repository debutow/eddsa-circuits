use std::marker::PhantomData;
use crate::utils::Field;
use plonk::Composer;
use plonk::composer::Variable;

pub struct Mux<F: Field>{
    _marker: PhantomData<F>,
}

impl<F: Field> Mux<F> {
    ///out == var_0 + selector * (var_1 - var_0)  (selector must be 0 or 1)
    pub fn mux1(
        cs: &mut Composer<F>,
        var_0: Variable,
        var_1: Variable,
        selector: Variable
    ) -> Variable {
        let value_s = cs.get_assignment(selector);
        let value_0 = cs.get_assignment(var_0);
        let value_1 = cs.get_assignment(var_1);

        let tmp: F = value_0 - value_0 * value_s;
        let var_tmp = cs.alloc(tmp);
        cs.poly_gate(
            vec![(var_0, F::one()), (selector, F::zero()), (var_tmp, -F::one())],
            -F::one(),
            F::zero(),
        );

        let out: F = tmp + value_1 * value_s;
        let var_out = cs.alloc(out);
        cs.poly_gate(
            vec![(var_1, F::zero()), (selector, F::zero()), (var_tmp, F::one()), (var_out, -F::one())],
            F::one(),
            F::zero(),
        );

        var_out
    }

    ///determine n vars with 1 selector
    #[allow(dead_code)]
    pub fn mux1_n(
        cs: &mut Composer<F>,
        var_0_n: &[Variable],
        var_1_n: &[Variable],
        selector: Variable
    ) -> Vec<Variable> {
        assert_eq!(var_0_n.len(), var_1_n.len());

        let mut var_out = Vec::new();

        for i in 0..var_0_n.len() {
            var_out.push(Mux::mux1(cs, var_0_n[i], var_1_n[i], selector));
        }

        var_out
    }

    ///choose 1 from n, with n selectors(only one is 1)
    #[allow(dead_code)]
    pub fn muxn_n(
        cs: &mut Composer<F>,
        var_n: &[Variable],
        selectors: &[Variable]
    ) -> Variable {
        assert_eq!(var_n.len(), selectors.len());

        let mut tmp = Composer::<F>::null();
        let mut tmp_out = Composer::<F>::null();

        for i in 0..selectors.len() {
            let tmp_out_value: F = cs.get_assignment(selectors[i])
                    * cs.get_assignment(var_n[i])
                    + cs.get_assignment(tmp);
            tmp_out = cs.alloc(tmp_out_value);

            // var[i] * selectors[i] + tmp = tmp_out
            cs.poly_gate(
                vec![(selectors[i], F::zero()), (var_n[i], F::zero()),
                     (tmp, F::one()), (tmp_out, -F::one())],
                F::one(),
                F::zero(),
            );

            tmp = tmp_out;

        }

        tmp_out
    }
}