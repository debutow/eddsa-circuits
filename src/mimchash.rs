#![allow(non_snake_case)]
use std::marker::PhantomData;
use crate::lib::mimc::{ROUND_KEYS, ROUNDS};
use crate::utils::Field;

pub struct MiMCHash<F: Field>{
    _marker: PhantomData<F>,
}

impl<F: Field> MiMCHash<F> {
    /// MiMC (follow circomlib's MiMCsponge
    /// rounds = 220
    /// x_l[i+1] = (k + x_l[i] + c[i])**5 + x_r[i]
    /// x_r[i+1] = x_l[i]
    /// if last round:
    /// x_l[i+1] = x_l[i]
    /// x_r[i+1] = x_r[i] + (k + x_l[i] + c[i])**5
    fn MiMC_feistel(
        l_data: F,
        r_data: F,
        k: F,
    ) -> (F, F){

        let mut x_l_i = l_data;
        let mut x_r_i = r_data;

        let mut x_l_ipp = F::zero();
        let mut x_r_ipp = F::zero();

        for i in 0..ROUND_KEYS.len()-1 {
            let t: F = k + x_l_i + F::from_str(ROUND_KEYS[i]).unwrap_or_default();

            let x_lp_value: F = t.square().square() * t + x_r_i;
            x_l_ipp = x_lp_value;

            x_r_ipp = x_l_i;

            //update tmp vars
            x_l_i = x_l_ipp;
            x_r_i = x_r_ipp;
        }

        //last round
        let t: F = k + x_l_i + F::from_str(ROUND_KEYS[ROUNDS-1]).unwrap_or_default();

        let x_ro_value: F = t.square().square() * t + x_r_i;
        let x_r_out = x_ro_value;

        let x_l_out = x_l_i;

        (x_l_out, x_r_out)
    }

    /// MiMC (follow circomlib's MiMCsponge
    pub fn MiMC_sponge(
        inputs: &[F],
        // outputs: &mut [Variable],
        k: F,
        // n_inputs: usize,
        n_outputs: usize,
    ) -> Vec<F>{
        assert!(n_outputs >= 1);
        let mut outputs = vec![];

        let (mut s_lout,mut s_rout) = MiMCHash::MiMC_feistel(inputs[0], F::zero(), k);

        for i in 1..inputs.len() {
            let tmp: F = s_lout + inputs[i];
            let (tmpl, tmpr) = MiMCHash::MiMC_feistel(tmp, s_rout, k);
            s_lout = tmpl;
            s_rout = tmpr;
        }

        outputs.push(s_lout);

        for i in 0..n_outputs-1 {
            let (tmpl, tmpr) = MiMCHash::MiMC_feistel(s_lout, s_rout, k);
            s_lout = tmpl;
            s_rout = tmpr;
            outputs.push(s_lout);
        }

        outputs
    }
}

#[cfg(test)]
mod tests {
    use plonk::MFr;
    use num_traits::Zero;
    use crate::mimchash::MiMCHash;

    #[test]
    fn test_hash() {
            let x1 = MFr::from(1 as u64);
            let x2 = MFr::from(2);
            let x3 = MFr::from(3 as u64);
            let x4 = MFr::from(4);

            let hash = MiMCHash::MiMC_sponge(
                &[x1, x2],
                MFr::zero(),
                1,
            );
            let hash = hash[0];
            println!("{}", hash);

            let hash = MiMCHash::MiMC_sponge(
                &[x1, x2, x3, x4],
                MFr::zero(),
                1,
            );
            let hash = hash[0];
            println!("{}", hash);

        // "0x2bcea035a1251603f1ceaf73cd4ae89427c47075bb8e3a944039ff1e3d6d2a6f"
        // "0x3e86bdc4eac70bd601473c53d8233b145fe8fd8bf6ef25f0b217a1da305665c"
    }
}