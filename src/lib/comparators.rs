use crate::utils::Field;
use plonk::Composer;
use plonk::composer::Variable;


/// if selector = 1, var_0 must eq var_1
pub fn enforce_equal_if_enabled<F: Field>(
    cs: &mut Composer<F>,
    var_0: Variable,
    var_1: Variable,
    selector: Variable
) {
    //value_0 - value_1
    let tmp: F = cs.get_assignment(var_0) - cs.get_assignment(var_1);
    let var_tmp = cs.alloc(tmp);
    cs.poly_gate(
        vec![(var_0, F::one()), (var_1, -F::one()), (var_tmp, -F::one())],
        F::zero(),
        F::zero(),
    );

    //s * (value_0 - value_1) === 0
    cs.poly_gate(
        vec![(selector, F::zero()), (var_tmp, F::zero())],
        F::one(),
        F::zero(),
    );
}
