use ark_poly::GeneralEvaluationDomain;
use ark_std::test_rng;
use ff::to_hex;
use num_bigint::{BigInt, BigUint};
use num_traits::One;
use plonk::{Composer, MFr};
use crate::eddsa_mimc::{new_key, verify};
use crate::lib::eddsa::EdDSA;

mod lib;
mod utils;
mod eddsa_mimc;
mod mimchash;

fn main() {
    //generate eddsa signature
    println!("generate eddsa signature:");
    let privatekey = new_key();
    let pubkey = privatekey.public();
    let msg = BigInt::parse_bytes(b"123456789012345678901234567890", 10).unwrap();

    //eddsa-mimcsponge sign
    let sig = privatekey.sign(msg.clone()).unwrap();
    println!("done");

    let v = verify(pubkey.clone(), sig.clone(), msg.clone());
    assert_eq!(v, true);

    //verify in circuits

    //create circuits
    println!("create circuits:");
    let mut cs = Composer::new(4);

    //witness in circuits
    //message
    let m = MFr::from(msg.to_biguint().unwrap());
    //pubkey
    let ax = MFr::from(BigUint::parse_bytes(to_hex(&pubkey.x).as_bytes(), 16).unwrap());
    let ay = MFr::from(BigUint::parse_bytes(to_hex(&pubkey.y).as_bytes(), 16).unwrap());
    //signature
    let r8x = MFr::from(BigUint::parse_bytes(to_hex(&sig.r_b8.x).as_bytes(), 16).unwrap());
    let r8y = MFr::from(BigUint::parse_bytes(to_hex(&sig.r_b8.y).as_bytes(), 16).unwrap());
    let s = MFr::from(sig.s.to_biguint().unwrap());

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
    println!("done");

    //circuits verify
    println!("circuits verify:");
    let rng = &mut test_rng();

    let prover_key = cs.compute_prover_key::<GeneralEvaluationDomain<MFr>>().unwrap();
    let mut prover = plonk::prover::Prover::new(&prover_key);
    prover.prove(&mut cs, rng).unwrap();
    println!("done.");
}
