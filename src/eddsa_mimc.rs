// copied from babyjubjub-rs, change hash function to mimc-sponge
// BabyJubJub elliptic curve implementation in Rust.
// For LICENSE check https://github.com/arnaucube/babyjubjub-rs
#![allow(dead_code)]
use plonk::MFr;

extern crate num;
extern crate num_bigint;
extern crate num_traits;

extern crate arrayref;
extern crate generic_array;

pub extern crate ff;
extern crate rand;

extern crate lazy_static;
use ff::*;
use num_traits::{One, Zero};
use std::cmp::min;
use std::convert::TryInto;
use ark_ff::PrimeField;
// use ark_ff::PrimeField;
// use ark_ff::{BigInteger, Field, PrimeField};
use arrayref::array_ref;

use ff::Field;
use ff::PrimeField as OtherPrimeField;

#[cfg(feature = "default")]
extern crate blake_hash; // compatible version with Blake used at circomlib
#[cfg(feature = "default")]
use blake_hash::Digest;

use generic_array::GenericArray;
use lazy_static::lazy_static;
use num_bigint::{BigInt, Sign, ToBigInt, RandBigInt, BigUint};
use crate::mimchash::MiMCHash;

// pub type Fr = plonk::MFr;

pub type Fr = poseidon_rs::Fr; // alias

#[cfg(feature = "default")]
fn blh(b: &[u8]) -> Vec<u8> {
    let hash = blake_hash::Blake512::digest(&b);
    hash.to_vec()
}


pub fn modulus(a: &BigInt, m: &BigInt) -> BigInt {
    ((a % m) + m) % m
}

pub fn modinv(a: &BigInt, q: &BigInt) -> Result<BigInt, String> {
    let big_zero: BigInt = Zero::zero();
    if a == &big_zero {
        return Err("no mod inv of Zero".to_string());
    }

    let mut mn = (q.clone(), a.clone());
    let mut xy: (BigInt, BigInt) = (Zero::zero(), One::one());

    while mn.1 != big_zero {
        xy = (xy.1.clone(), xy.0 - (mn.0.clone() / mn.1.clone()) * xy.1);
        mn = (mn.1.clone(), modulus(&mn.0, &mn.1));
    }

    while xy.0 < Zero::zero() {
        xy.0 = modulus(&xy.0, q);
    }
    Ok(xy.0)
}

pub fn modsqrt(a: &BigInt, q: &BigInt) -> Result<BigInt, String> {
    // Tonelli-Shanks Algorithm (https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm)
    //
    // This implementation is following the Go lang core implementation https://golang.org/src/math/big/int.go?s=23173:23210#L859
    // Also described in https://www.maa.org/sites/default/files/pdf/upload_library/22/Polya/07468342.di020786.02p0470a.pdf
    // -> section 6

    let zero: BigInt = Zero::zero();
    let one: BigInt = One::one();
    if legendre_symbol(&a, q) != 1 || a == &zero || q == &2.to_bigint().unwrap() {
        return Err("not a mod p square".to_string());
    } else if q % 4.to_bigint().unwrap() == 3.to_bigint().unwrap() {
        let r = a.modpow(&((q + one) / 4), &q);
        return Ok(r);
    }

    let mut s = q - &one;
    let mut e: BigInt = Zero::zero();
    while &s % 2 == zero {
        s >>= 1;
        e += &one;
    }

    let mut n: BigInt = 2.to_bigint().unwrap();
    while legendre_symbol(&n, q) != -1 {
        n = &n + &one;
    }

    let mut y = a.modpow(&((&s + &one) >> 1), q);
    let mut b = a.modpow(&s, q);
    let mut g = n.modpow(&s, q);
    let mut r = e;

    loop {
        let mut t = b.clone();
        let mut m: BigInt = Zero::zero();
        while t != one {
            t = modulus(&(&t * &t), q);
            m += &one;
        }

        if m == zero {
            return Ok(y);
        }

        t = g.modpow(&(2.to_bigint().unwrap().modpow(&(&r - &m - 1), q)), q);
        g = g.modpow(&(2.to_bigint().unwrap().modpow(&(r - &m), q)), q);
        y = modulus(&(y * t), q);
        b = modulus(&(b * &g), q);
        r = m.clone();
    }
}

pub fn concatenate_arrays<T: Clone>(x: &[T], y: &[T]) -> Vec<T> {
    x.iter().chain(y).cloned().collect()
}

pub fn legendre_symbol(a: &BigInt, q: &BigInt) -> i32 {
    // returns 1 if has a square root modulo q
    let one: BigInt = One::one();
    let ls: BigInt = a.modpow(&((q - &one) >> 1), &q);
    if ls == q - one {
        return -1;
    }
    1
}

lazy_static! {
    static ref D: Fr = Fr::from_str("168696").unwrap();
    static ref D_BIG: BigInt = BigInt::parse_bytes(b"168696", 10).unwrap();
    static ref A: Fr = Fr::from_str("168700").unwrap();
    static ref A_BIG: BigInt = BigInt::parse_bytes(b"168700", 10).unwrap();
    pub static ref Q: BigInt = BigInt::parse_bytes(
        b"21888242871839275222246405745257275088548364400416034343698204186575808495617",10
    )
        .unwrap();
    static ref B8: Point = Point {
        x: Fr::from_str(
               "5299619240641551281634865583518297030282874472190772894086521144482721001553",
           )
            .unwrap(),
            y: Fr::from_str(
                "16950150798460657717958625567821834550301663161624707787222815936182638968203",
            )
                .unwrap(),
    };
    static ref ORDER: Fr = Fr::from_str(
        "21888242871839275222246405745257275088614511777268538073601725287587578984328",
    )
        .unwrap();

    // SUBORDER = ORDER >> 3
    static ref SUBORDER: BigInt = &BigInt::parse_bytes(
        b"21888242871839275222246405745257275088614511777268538073601725287587578984328",
        10,
    )
        .unwrap()
        >> 3;
    // static ref POSEIDON: poseidon_rs::Poseidon = Poseidon::new();
}

#[derive(Clone, Debug)]
pub struct PointProjective {
    pub x: Fr,
    pub y: Fr,
    pub z: Fr,
}

impl PointProjective {
    pub fn affine(&self) -> Point {
        if self.z.is_zero() {
            return Point {
                x: Fr::zero(),
                y: Fr::zero(),
            };
        }

        let zinv = self.z.inverse().unwrap();
        let mut x = self.x;
        x.mul_assign(&zinv);
        let mut y = self.y;
        y.mul_assign(&zinv);

        Point { x, y }
    }

    #[allow(clippy::many_single_char_names)]
    pub fn add(&self, q: &PointProjective) -> PointProjective {
        // add-2008-bbjlp https://hyperelliptic.org/EFD/g1p/auto-twisted-projective.html#doubling-dbl-2008-bbjlp
        let mut a = self.z;
        a.mul_assign(&q.z);
        let mut b = a;
        b.square();
        let mut c = self.x;
        c.mul_assign(&q.x);
        let mut d = self.y;
        d.mul_assign(&q.y);
        let mut e = *D;
        e.mul_assign(&c);
        e.mul_assign(&d);
        let mut f = b;
        f.sub_assign(&e);
        let mut g = b;
        g.add_assign(&e);
        let mut x1y1 = self.x;
        x1y1.add_assign(&self.y);
        let mut x2y2 = q.x;
        x2y2.add_assign(&q.y);
        let mut aux = x1y1;
        aux.mul_assign(&x2y2);
        aux.sub_assign(&c);
        aux.sub_assign(&d);
        let mut x3 = a;
        x3.mul_assign(&f);
        x3.mul_assign(&aux);
        let mut ac = *A;
        ac.mul_assign(&c);
        let mut dac = d;
        dac.sub_assign(&ac);
        let mut y3 = a;
        y3.mul_assign(&g);
        y3.mul_assign(&dac);
        let mut z3 = f;
        z3.mul_assign(&g);

        PointProjective {
            x: x3,
            y: y3,
            z: z3,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Point {
    pub x: Fr,
    pub y: Fr,
}

impl Point {
    pub fn projective(&self) -> PointProjective {
        PointProjective {
            x: self.x,
            y: self.y,
            z: Fr::one(),
        }
    }

    pub fn mul_scalar(&self, n: &BigInt) -> Point {
        let mut r: PointProjective = PointProjective {
            x: Fr::zero(),
            y: Fr::one(),
            z: Fr::one(),
        };
        let mut exp: PointProjective = self.projective();
        let (_, b) = n.to_bytes_le();

        for i in 0..n.bits() {
            if test_bit(&b, i.try_into().unwrap()) {
                r = r.add(&exp);
            }
            exp = exp.add(&exp);
        }

        r.affine()
    }

    pub fn compress(&self) -> [u8; 32] {
        let p = &self;
        let mut r: [u8; 32] = [0; 32];
        let x_big = BigInt::parse_bytes(to_hex(&p.x).as_bytes(), 16).unwrap();
        let y_big = BigInt::parse_bytes(to_hex(&p.y).as_bytes(), 16).unwrap();
        // let x_big: BigUint = p.x.into_repr().into();
        // let x_big = x_big.to_bigint().unwrap();
        // let y_big: BigUint = p.y.into_repr().into();
        // let y_big = y_big.to_bigint().unwrap();
        let (_, y_bytes) = y_big.to_bytes_le();
        let len = min(y_bytes.len(), r.len());
        r[..len].copy_from_slice(&y_bytes[..len]);
        if x_big > (&Q.clone() >> 1 as i32) {
            r[31] |= 0x80;
        }
        r
    }

    pub fn equals(&self, p: Point) -> bool {
        if self.x == p.x && self.y == p.y {
            return true;
        }
        false
    }
}

///第i位=1时，输出1
pub fn test_bit(b: &[u8], i: usize) -> bool {
    b[i / 8] & (1 << (i % 8)) != 0
}

pub fn decompress_point(bb: [u8; 32]) -> Result<Point, String> {
    // https://tools.ietf.org/html/rfc8032#section-5.2.3
    let mut sign: bool = false;
    let mut b = bb;
    if b[31] & 0x80 != 0x00 {
        sign = true;
        b[31] &= 0x7F;
    }
    let y: BigInt = BigInt::from_bytes_le(Sign::Plus, &b[..]);
    if y >= Q.clone() {
        return Err("y outside the Finite Field over R".to_string());
    }
    let one: BigInt = One::one();

    // x^2 = (1 - y^2) / (a - d * y^2) (mod p)
    let den = modinv(
        &modulus(
            &(&A_BIG.clone() - modulus(&(&D_BIG.clone() * (&y * &y)), &Q)),
            &Q,
        ),
        &Q,
    )?;
    let mut x: BigInt = modulus(&((one - modulus(&(&y * &y), &Q)) * den), &Q);
    x = modsqrt(&x, &Q)?;

    if sign && (x <= (&Q.clone() >> 1)) || (!sign && (x > (&Q.clone() >> 1))) {
        x *= -(1.to_bigint().unwrap());
    }
    x = modulus(&x, &Q);
    let x_fr: Fr = Fr::from_str(&x.to_string()).unwrap();
    let y_fr: Fr = Fr::from_str(&y.to_string()).unwrap();
    Ok(Point { x: x_fr, y: y_fr })
}

#[derive(Debug, Clone)]
pub struct Signature {
    pub r_b8: Point,
    pub s: BigInt,
}

impl Signature {
    pub fn compress(&self) -> [u8; 64] {
        let mut b: Vec<u8> = Vec::new();
        b.append(&mut self.r_b8.compress().to_vec());
        let (_, s_bytes) = self.s.to_bytes_le();
        let mut s_32bytes: [u8; 32] = [0; 32];
        let len = min(s_bytes.len(), s_32bytes.len());
        s_32bytes[..len].copy_from_slice(&s_bytes[..len]);
        b.append(&mut s_32bytes.to_vec());
        let mut r: [u8; 64] = [0; 64];
        r[..].copy_from_slice(&b[..]);
        r
    }
}

pub fn decompress_signature(b: &[u8; 64]) -> Result<Signature, String> {
    let r_b8_bytes: [u8; 32] = *array_ref!(b[..32], 0, 32);
    let s: BigInt = BigInt::from_bytes_le(Sign::Plus, &b[32..]);
    let r_b8 = decompress_point(r_b8_bytes);
    match r_b8 {
        Result::Err(err) => Err(err),
        Result::Ok(res) => Ok(Signature { r_b8: res, s }),
    }
}

pub struct PrivateKey {
    key: [u8; 32],
}

impl PrivateKey {
    pub fn import(b: Vec<u8>) -> Result<PrivateKey, String> {
        if b.len() != 32 {
            return Err(String::from("imported key can not be bigger than 32 bytes"));
        }
        let mut sk: [u8; 32] = [0; 32];
        sk.copy_from_slice(&b[..32]);
        Ok(PrivateKey { key: sk })
    }

    pub fn scalar_key(&self) -> BigInt {
        // not-compatible with circomlib implementation, but using Blake2b
        // let mut hasher = Blake2b::new();
        // hasher.update(sk_raw_bytes);
        // let mut h = hasher.finalize();

        // compatible with circomlib implementation
        let hash: Vec<u8> = blh(&self.key.to_vec());
        let mut h: Vec<u8> = hash[..32].to_vec();

        h[0] &= 0xF8;
        h[31] &= 0x7F;
        h[31] |= 0x40;

        let sk = BigInt::from_bytes_le(Sign::Plus, &h[..]);
        sk >> 3
    }

    pub fn public(&self) -> Point {
        B8.mul_scalar(&self.scalar_key())
    }

    ///use mimc-sponge
    pub fn sign(&self, msg: BigInt) -> Result<Signature, String> {
        if msg > Q.clone() {
            return Err("msg outside the Finite Field".to_string());
        }
        // let (_, sk_bytes) = self.key.to_bytes_le();
        // let mut hasher = Blake2b::new();
        // hasher.update(sk_bytes);
        // let mut h = hasher.finalize(); // h: hash(sk), s: h[32:64]
        let mut h: Vec<u8> = blh(&self.key.to_vec());

        let (_, msg_bytes) = msg.to_bytes_le();
        let mut msg32: [u8; 32] = [0; 32];
        msg32[..msg_bytes.len()].copy_from_slice(&msg_bytes[..]);
        let msg_fr: Fr = Fr::from_str(&msg.to_string()).unwrap();

        // println!("{}", msg.to_string());

        // https://tools.ietf.org/html/rfc8032#section-5.1.6
        let s = GenericArray::<u8, generic_array::typenum::U32>::from_mut_slice(&mut h[32..64]);
        let r_bytes = concatenate_arrays(s, &msg32);
        let r_hashed: Vec<u8> = blh(&r_bytes);
        let mut r = BigInt::from_bytes_le(Sign::Plus, &r_hashed[..]);
        r = modulus(&r, &SUBORDER);
        let r_b8: Point = B8.mul_scalar(&r);
        let a = &self.public();

        // let hm_input = vec![r_b8.x, r_b8.y, a.x, a.y, msg_fr];
        // let hm = POSEIDON.hash(hm_input)?;

        let hm_input = vec![
            MFr::from(BigUint::parse_bytes(to_hex(&r_b8.x).as_bytes(), 16).unwrap()),
            MFr::from(BigUint::parse_bytes(to_hex(&r_b8.y).as_bytes(), 16).unwrap()),
            MFr::from(BigUint::parse_bytes(to_hex(&a.x).as_bytes(), 16).unwrap()),
            MFr::from(BigUint::parse_bytes(to_hex(&a.y).as_bytes(), 16).unwrap()),
            MFr::from(BigUint::parse_bytes(to_hex(&msg_fr).as_bytes(), 16).unwrap()),
        ];

        let hm = MiMCHash::MiMC_sponge(&hm_input, MFr::zero(), 1);
        let hm = hm[0];

        let mut s = &self.scalar_key() << 3 as i32;
        // let hm_b = BigInt::parse_bytes(to_hex(&hm).as_bytes(), 16).unwrap();
        let hm_big: BigUint = hm.into_repr().into();
        let hm_b = hm_big.to_bigint().unwrap();

        s = hm_b * s;
        s = r + s;
        s %= &SUBORDER.clone();

        Ok(Signature { r_b8, s })
    }


}

pub fn new_key() -> PrivateKey {
    // https://tools.ietf.org/html/rfc8032#section-5.1.5
    let mut rng = rand::thread_rng();
    let sk_raw = rng.gen_biguint(1024).to_bigint().unwrap();
    let (_, sk_raw_bytes) = sk_raw.to_bytes_be();
    PrivateKey::import(sk_raw_bytes[..32].to_vec()).unwrap()
}

pub fn verify(pk: Point, sig: Signature, msg: BigInt) -> bool {
    if msg > Q.clone() {
        return false;
    }
    let msg_fr: Fr = Fr::from_str(&msg.to_string()).unwrap();
    // let hm_input = vec![sig.r_b8.x, sig.r_b8.y, pk.x, pk.y, msg_fr];
    // let hm = match POSEIDON.hash(hm_input) {
    //     Result::Err(_) => return false,
    //     Result::Ok(hm) => hm,
    // };

    let hm_input = vec![
        MFr::from(BigUint::parse_bytes(to_hex(&sig.r_b8.x).as_bytes(), 16).unwrap()),
        MFr::from(BigUint::parse_bytes(to_hex(&sig.r_b8.y).as_bytes(), 16).unwrap()),
        MFr::from(BigUint::parse_bytes(to_hex(&pk.x).as_bytes(), 16).unwrap()),
        MFr::from(BigUint::parse_bytes(to_hex(&pk.y).as_bytes(), 16).unwrap()),
        MFr::from(BigUint::parse_bytes(to_hex(&msg_fr).as_bytes(), 16).unwrap()),
    ];

    let hm = MiMCHash::MiMC_sponge(&hm_input, MFr::zero(), 1);
    let hm = hm[0];

    let l = B8.mul_scalar(&sig.s);
    // let hm_b = BigInt::parse_bytes(to_hex(&hm).as_bytes(), 16).unwrap();
    let hm_big: BigUint = hm.into_repr().into();
    let hm_b = hm_big.to_bigint().unwrap();


    let r = sig
        .r_b8
        .projective()
        .add(&pk.mul_scalar(&((8 as i32).to_bigint().unwrap() * hm_b)).projective());
    l.equals(r.affine())
}

#[cfg(test)]
mod tests {
    use ark_ff::PrimeField;
    use rustc_hex::ToHex;
    use super::*;

    #[test]
    fn test_new_key_sign_verify_0() {
        let sk = new_key();
        let pk = sk.public();
        let msg = 5.to_bigint().unwrap();
        let sig = sk.sign(msg.clone()).unwrap();
        let v = verify(pk, sig, msg);
        assert_eq!(v, true);
    }

    #[test]
    fn test_new_key_sign_verify_1() {
        let sk = new_key();
        let pk = sk.public();
        let msg = BigInt::parse_bytes(b"123456789012345678901234567890", 10).unwrap();
        let sig = sk.sign(msg.clone()).unwrap();
        let v = verify(pk, sig, msg);
        assert_eq!(v, true);
    }

    #[test]
    fn test_circomlib_testvector() {

        let sk_raw_bytes =
            hex::decode("0001020304050607080900010203040506070809000102030405060708090001")
                .unwrap();

        // test blake compatible with circomlib implementation
        let h: Vec<u8> = blh(&sk_raw_bytes);
        assert_eq!(h.to_hex(), "c992db23d6290c70ffcc02f7abeb00b9d00fa8b43e55d7949c28ba6be7545d3253882a61bd004a236ef1cdba01b27ba0aedfb08eefdbfb7c19657c880b43ddf1");

        // test private key
        let sk = PrivateKey::import(
            hex::decode("0001020304050607080900010203040506070809000102030405060708090001")
                .unwrap(),
        )
            .unwrap();
        assert_eq!(
            sk.scalar_key().to_string(),
            "6466070937662820620902051049739362987537906109895538826186780010858059362905"
        );
        println!("{:?}", sk.key);

        // test public key
        let pk = sk.public();

        assert_eq!(
            pk.x.to_string(),
            "Fr(0x1d5ac1f31407018b7d413a4f52c8f74463b30e6ac2238220ad8b254de4eaa3a2)"
        );
        assert_eq!(
            pk.y.to_string(),
            "Fr(0x1e1de8a908826c3f9ac2e0ceee929ecd0caf3b99b3ef24523aaab796a6f733c4)"
        );

        // test signature & verification
        let msg = BigInt::from_bytes_le(Sign::Plus, &hex::decode("00010203040506070809").unwrap());
        println!("msg {:?}", msg.to_string());
        let sig = sk.sign(msg.clone()).unwrap();
        assert_eq!(
            sig.r_b8.x.to_string(),
            "Fr(0x192b4e51adf302c8139d356d0e08e2404b5ace440ef41fc78f5c4f2428df0765)"
        );
        assert_eq!(
            sig.r_b8.y.to_string(),
            "Fr(0x2202bebcf57b820863e0acc88970b6ca7d987a0d513c2ddeb42e3f5d31b4eddf)"
        );
        assert_eq!(
            sig.s.to_string(),
            "1868336918738674306327358602987493427631678603535639134028485964115448322340"
        );
        let v = verify(pk, sig, msg);
        assert_eq!(v, true);
    }
}