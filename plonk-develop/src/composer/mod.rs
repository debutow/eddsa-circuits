use ark_std::{cfg_iter, cfg_iter_mut, cmp::max, format, vec, vec::Vec};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::prover::ProverKey;
use crate::{Domain, Error, Field, Map};

mod arithmetic;

mod permutation;

mod lookup;
pub use lookup::Table;

mod range;

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
struct Wire {
    pub col: usize, //在哪一列
    pub row: usize,
}

impl Wire {
    fn new(col: usize, row: usize) -> Self {
        Self { col, row }
    }
}

#[derive(Debug, Eq, PartialEq, Clone, Copy, Hash, Ord, PartialOrd, Default)]
pub struct Variable(usize);

#[derive(Debug, Default)]
pub struct Composer<F: Field> {
    pub program_width: usize,

    size: usize,
    is_finalized: bool,

    wires: Map<String, Vec<Variable>>,
    selectors: Map<String, Vec<F>>,
    public_input: Vec<Variable>,

    assignments: Vec<F>,
    epicycles: Vec<Vec<Wire>>,
    tables: Vec<Table<F>>,
}

/// basics
impl<F: Field> Composer<F> {
    const SELECTOR_LABELS: [&'static str; 6] =
        ["q_m", "q_c", "q_arith", "q_range", "q_lookup", "q_table"];

    pub fn new(program_width: usize) -> Composer<F> {
        let mut cs = Composer::default();

        for j in 0..program_width {
            cs.wires.insert(format!("w_{}", j), Vec::new());
            cs.selectors.insert(format!("q_{}", j), Vec::new());
        }
        for label in Self::SELECTOR_LABELS {
            cs.selectors.insert(label.to_string(), Vec::new());
        }
        cs.program_width = program_width;
        cs.alloc(F::zero());

        cs
    }

    #[inline]
    pub fn size(&self) -> usize {
        self.size
    }

    #[inline]
    pub fn input_size(&self) -> usize {
        self.public_input.len()
    }

    pub fn alloc(&mut self, value: F) -> Variable {
        let var = self.alloc_variable(value);

        var
    }

    ///alloc public input variable
    pub fn alloc_input(&mut self, value: F) -> Variable {
        let var = self.alloc_variable(value);
        self.public_input.push(var);

        var
    }

    //change to pub
    pub fn null() -> Variable {
        Variable(0)
    }

    #[inline]
    fn alloc_variable(&mut self, value: F) -> Variable {
        let var = Variable(self.epicycles.len());
        //每个变量都会新建一个epicycles（copy constraint）
        self.epicycles.push(Vec::new());
        self.assignments.push(value);

        var
    }

    #[inline]
    //只负责加一行，其它值都用0填充，由别的（调用它的函数）改成所需
    fn insert_gate(&mut self, mut gate_wires: Vec<Variable>) -> usize {
        assert!(gate_wires.len() <= self.program_width);
        while gate_wires.len() < self.program_width {
            gate_wires.push(Self::null());
        }

        let i = self.size;
        for j in 0..self.program_width {
            let wire_label = format!("w_{}", j);
            self.wires.get_mut(&wire_label).unwrap().push(gate_wires[j]);
            self.epicycles[gate_wires[j].0].push(Wire::new(j, i));

            let scaling_label = format!("q_{}", j);
            self.selectors
                .get_mut(&scaling_label)
                .unwrap()
                .push(F::zero());
        }

        for label in Self::SELECTOR_LABELS {
            self.selectors.get_mut(label).unwrap().push(F::zero());
        }
        self.size += 1;

        i
    }

    fn get_assignments(&self, vars: &[Variable]) -> Vec<F> {
        cfg_iter!(vars).map(|&v| self.assignments[v.0]).collect()
    }

    pub fn get_assignment(&self, var: Variable) -> F {
        self.assignments[var.0]
    }
}

//都要先.finalize()把pi放进wires里
/// synthesis
impl<F: Field> Composer<F> {
    pub fn compute_prover_key<D: Domain<F>>(&mut self) -> Result<ProverKey<F, D>, Error> {
        self.finalize();

        let size = max(self.size(), self.sorted_size());

        let mut prover_key = ProverKey::new(size + 1, self.input_size(), self.program_width)?;

        for (k, q) in self.selectors.iter() {
            prover_key.insert(&k, q.clone());
        }

        let sigmas = self.compute_sigmas(prover_key.domain);
        for (col, sigma) in sigmas.iter().enumerate() {
            prover_key.insert(&format!("sigma_{}", col), sigma.clone());
        }

        // the last column for table indices
        let table_values = self.compute_table_values();
        for (col, table_value) in table_values.iter().enumerate() {
            prover_key.insert(&format!("table_{}", col), table_value.clone());
        }

        Ok(prover_key)
    }

    pub(crate) fn compute_public_input(&mut self) -> Vec<F> {
        self.finalize();

        cfg_iter!(self.public_input)
            .map(|v| self.assignments[v.0])
            .collect()
    }

    pub(crate) fn compute_wire_values(&mut self) -> Result<Map<String, Vec<F>>, Error> {
        self.finalize();

        let mut wires = Map::new();

        let assign = |v: &Variable| self.assignments[v.0];
        for (l, w) in self.wires.iter() {
            wires.insert(l.to_string(), cfg_iter!(w).map(assign).collect());
        }
        let sorted_values = self.compute_sorted_values();
        for (col, sorted_value) in sorted_values.iter().enumerate() {
            wires.insert(format!("s_{}", col), sorted_value.clone());
        }

        Ok(wires)
    }

    fn finalize(&mut self) {
        if self.is_finalized {
            return;
        };

        let input_size = self.input_size();
        self.size += input_size;

        for epicycle in self.epicycles.iter_mut() {
            *epicycle = cfg_iter_mut!(epicycle)
                //将所有wire整体下移input size行
                .map(|w| Wire::new(w.col, w.row + input_size))
                .collect()
        }

        for (i, var) in self.public_input.iter().enumerate() {
            self.epicycles[var.0].push(Wire::new(0, i));
            for col in 1..self.program_width {
                //往var0所在的epicycle装点东西。也就是说pi所有var对应的行的其它地方要连到0上。
                self.epicycles[Self::null().0].push(Wire::new(col, i));
            }
        }

        let mut wires = Map::new();
        //更新，将pi放到wires w_0的最前面，其它填0
        for (label, wire) in self.wires.iter_mut() {
            let mut vars = if label == "w_0" {
                self.public_input.clone()
            } else {
                vec![Self::null(); input_size]
            };
            vars.append(wire);
            wires.insert(label.to_string(), vars);
        }
        self.wires = wires;

        let mut selectors = Map::new();
        for (label, selector) in self.selectors.iter_mut() {
            let mut values = if label == "q_0" {
                vec![F::one(); input_size]
            } else {
                vec![F::zero(); input_size]
            };
            values.append(selector);
            selectors.insert(label.to_string(), values);
        }
        self.selectors = selectors;

        self.is_finalized = true;
    }
}
