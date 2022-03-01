use crate::Composer;
use ark_ff::{PrimeField, Field};
use crate::composer::Variable;

#[derive(Debug)]
pub struct MimcC<F: Field>{
    n_rounds: usize,
    mimc_c: Vec<F>,
}

impl<F: Field> MimcC<F>{
    pub fn new() -> Self{
        MimcC {
            n_rounds: 0,
            mimc_c: Vec::new(),
        }
    }

    pub fn init_mimc_c(&mut self, n_rounds: usize, mimc_c: Vec<F>){
        assert!(mimc_c.len() >= 1);
        assert!(n_rounds >= 1);
        self.mimc_c = mimc_c;
        self.n_rounds = n_rounds;
    }
}

impl<F: PrimeField + Field> Composer<F> {

    fn create_mimc_gates(
        &mut self,
        l: (Variable, F), // w_l, q_l
        r: (Variable, F), // w_r, q_r
        o: (Variable, F), // w_o, q_o
        aux: Option<(Variable, F)>, //w0
        q_m: F,
        q_c: F,
        pi: F,
        q_arith: F,
        q_mimc: F,
    ) {
        let aux = match aux {
            Some(aux) => aux,
            None => (self.null_var, F::zero()),
        };
        let index = self.n;
        self.permutation.insert_gate(aux.0, l.0, r.0, o.0, index);

        self.w_0.push(aux.0);
        self.w_1.push(l.0);
        self.w_2.push(r.0);
        self.w_3.push(o.0);
        self.pi.push(pi);

        self.q_0.push(aux.1);
        self.q_1.push(l.1);
        self.q_2.push(r.1);
        self.q_3.push(o.1);
        self.q_m.push(q_m);
        self.q_c.push(q_c);
        self.q_arith.push(q_arith);
        self.q_range.push(F::zero());
        self.q_mimc.push(q_mimc);
        //self.q_mimc_c.push(q_mimc_c);
        self.q_lookup.push(F::zero());

        self.n += 1;
    }

    //Return the hash value
    //input a Vec<F> as data
    //put q_mimc_c to 'w3'
    pub fn create_mimc_hash(
        &mut self,
        var: Variable,
        original_data: Vec<F>,
    ) -> F {
        let n_data =  original_data.len();
        let n_mimc_c = self.mimc_c_container.mimc_c.len();
        assert!(n_mimc_c >= 1);

        //let q_mimc_ci = self.q_mimc_c.get(0).unwrap().clone();
        let q_mimc_ci = self.mimc_c_container.mimc_c.get(0).unwrap().clone();
        let w1 = original_data.get(0).unwrap().clone();

        let tmp :F = F::zero() + w1.clone() + q_mimc_ci.clone();
        let w2 = tmp.square() * tmp.clone();

        let v1 = self.alloc_and_assign(w1);
        let v2 = self.alloc_and_assign(w2);
        let v3 = self.alloc_and_assign(q_mimc_ci.clone());
        let v0 = self.alloc_and_assign(F::zero());
        self.create_mimc_gates(
            (v1, F::zero()),
            (v2, F::zero()),
            (v3, F::zero()),
            Option::from((v0, F::zero())),
            F::zero(),
            F::zero(),
            F::zero(),
            F::zero(),
            F::one(),
        );
        let mut w0_omega:F  = w2.square() * tmp;
        for n_datum in 1..n_data {
            let q_mimc_ci = self.mimc_c_container.mimc_c.get(n_datum % n_mimc_c).unwrap().clone();
            let w1 = original_data.get(n_datum).unwrap().clone();

            let tmp :F = w0_omega.clone() + w1.clone() + q_mimc_ci.clone();
            let w2 = tmp.square() * tmp.clone();

            let v1 = self.alloc_and_assign(w1);
            let v2 = self.alloc_and_assign(w2);
            let v3 = self.alloc_and_assign(q_mimc_ci.clone());
            let v0 = self.alloc_and_assign(w0_omega.clone());
            self.create_mimc_gates(
                (v1, F::zero()),
                (v2, F::zero()),
                (v3, F::zero()),
                Option::from((v0, F::zero())),
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
                F::one(),
            );
            w0_omega = w2.square() * tmp;
        }
        //store hash. then need an arith gate to make it equal to another var
        let v0 = self.alloc_and_assign(w0_omega.clone());
        //this q_mimc must equal 0
        self.create_mimc_gates(
            (var, -F::one()),
            (var, F::zero()),
            (var, F::zero()),
            Option::from((v0, F::one())),
            F::zero(),
            -F::zero(),
            F::zero(),
            F::one(),
            F::zero(),
        );
        w0_omega
    }

    /// MiMC (follow circomlib's MiMCsponge
    /// rounds = 220
    /// i=0: c[i]=0
    ///
    /// x_l[i+1] = (k + x_l[i] + c[i])**5 + x_r[i]
    /// x_r[i+1] = x_l[i]
    ///
    ///
    ///(no use below)
    /// rounds = 220
    /// x_l[i+1] = (x_l[i] + c[i])**3 + x_r[i]
    /// x_r[i+1] = x_l[i]
    /// hash = x_l[rounds]
    pub fn create_mimc_hash_feistel(
        &mut self,
        var: Variable,
        l_data: F,
        r_data: F,
        k: F,
    ) -> F {
        let n_rounds =  self.mimc_c_container.n_rounds;
        let n_mimc_c = self.mimc_c_container.mimc_c.len();
        assert!(n_mimc_c >= 1);
        assert!(n_rounds >= 1);

        let q_mimc_ci = self.mimc_c_container.mimc_c.get(0).unwrap().clone();
        let mut w0 = l_data.clone();
        let w1 = r_data.clone();
        let tmp :F = w0.clone() + q_mimc_ci.clone();
        let w3 :F = tmp.square() * tmp.clone();

        let v0 = self.alloc_and_assign(w0.clone());
        let v1 = self.alloc_and_assign(w1);
        let v2 = self.alloc_and_assign(q_mimc_ci.clone());
        let v3 = self.alloc_and_assign(w3);
        self.create_mimc_gates(
            (v1, F::zero()),
            (v2, F::zero()),
            (v3, F::zero()),
            Option::from((v0.clone(), F::zero())),
            F::zero(),
            F::zero(),
            F::zero(),
            F::zero(),
            F::one(),
        );
        let mut w0_omega:F  = w3.clone() + w1.clone();
        let mut last_v0 = v0;
        for n_round in 1..n_rounds {
            let q_mimc_ci = self.mimc_c_container.mimc_c.get(n_round % n_mimc_c).unwrap().clone();
            let w1 = w0.clone();
            let tmp :F = w0_omega.clone() + q_mimc_ci.clone();
            let w3 :F = tmp.square() * tmp.clone();

            let v0 = self.alloc_and_assign(w0_omega.clone());
            let v2 = self.alloc_and_assign(q_mimc_ci.clone());
            let v3 = self.alloc_and_assign(w3);
            self.create_mimc_gates(
                (last_v0, F::zero()),
                (v2, F::zero()),
                (v3, F::zero()),
                Option::from((v0, F::zero())),
                F::zero(),
                F::zero(),
                F::zero(),
                F::zero(),
                F::one(),
            );
            w0_omega = w3.clone() + w1.clone();
            last_v0 = v0;
        }

        let v0 = self.alloc_and_assign(w0_omega.clone());
        //this q_mimc must equal 0
        self.create_mimc_gates(
            (var, -F::one()),
            (var, F::zero()),
            (var, F::zero()),
            Option::from((v0, F::one())),
            F::zero(),
            -F::zero(),
            F::zero(),
            F::one(),
            F::zero(),
        );
        w0_omega
    }
}
