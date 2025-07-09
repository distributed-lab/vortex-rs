use p3_dft::{Radix2DFTSmallBatch, TwoAdicSubgroupDft};
use p3_field::integers::QuotientMap;
use p3_field::{PrimeCharacteristicRing, PrimeField32};
use p3_koala_bear::KoalaBear;
use p3_maybe_rayon::prelude::*;
use rayon::current_num_threads;
use std::slice::Iter;

pub struct RSis {
    a: Vec<Vec<KoalaBear>>,
    ag: Vec<Vec<KoalaBear>>,
    max_nb_elements_to_hash: usize,
}

const KOALA_BEAR_BITS: usize = 31;
const KOALA_BEAR_BYTES: usize = 4;

const LOG_TWO_DEGREE: usize = 9;
const LOG_TWO_BOUND: usize = 16;
const LIMB_SIZE: usize = LOG_TWO_BOUND / 8;

pub const DEGREE: usize = 1 << LOG_TWO_DEGREE;

impl RSis {
    pub fn new(
        seed: u64,
        //log_two_degree: usize,
        //log_two_bound: usize,
        max_nb_elements_to_hash: usize,
    ) -> Self {
        // assert!(
        //     log_two_bound <= 64 && log_two_bound <= KOALA_BEAR_BITS,
        //     "logTwoBound too large"
        // );
        // assert_eq!(log_two_bound % 8, 0, "logTwoBound must be a multiple of 8");

        // assert_eq!(
        //     KOALA_BEAR_BYTES % nb_bytes_per_limb,
        //     0,
        //     "nbBytesPerLimb must divide field size"
        // );
        let mut n = (KOALA_BEAR_BYTES / LIMB_SIZE) * max_nb_elements_to_hash;

        if n % DEGREE == 0 {
            n /= DEGREE;
        } else {
            n /= DEGREE;
            n += 1;
        }

        let mut r = Self {
            a: vec![vec![]; n],
            ag: vec![vec![]; n],
            max_nb_elements_to_hash,
        };

        let n = n;

        r.a = (0..n)
            .into_par_iter()
            .chunks(current_num_threads())
            .map(|indexes| {
                let mut res = Vec::with_capacity(indexes.len());
                for i in indexes {
                    let mut a = Vec::with_capacity(DEGREE);
                    for j in 0..DEGREE {
                        a.push(derive_random_element_from_seed(seed, i as u64, j as u64));
                    }
                    res.push(a);
                }
                res
            })
            .flatten()
            .collect();

        r.ag = (0..n)
            .into_par_iter()
            .chunks(current_num_threads())
            .map(|indexes| {
                let mut res = Vec::with_capacity(indexes.len());
                let dft = Radix2DFTSmallBatch::new(DEGREE);

                for i in indexes {
                    let ag = dft.dft(r.a[i].clone());
                    res.push(ag);
                }
                res
            })
            .flatten()
            .collect();

        r
    }

    pub fn hash(&self, v: &Vec<KoalaBear>, dft: &Radix2DFTSmallBatch<KoalaBear>) -> Vec<KoalaBear> {
        assert!(
            v.len() <= self.max_nb_elements_to_hash,
            "can't hash more than configured elements with params provided in constructor"
        );
        let mut res = vec![KoalaBear::ZERO; DEGREE];

        for i in 0..self.ag.len() {
            res = self.inner_hash(res, &mut v.iter(), i, &dft);
        }

        dft.idft(res)
    }

    fn inner_hash(
        &self,
        mut state: Vec<KoalaBear>,
        vit: &mut Iter<KoalaBear>,
        pol_id: usize,
        dft: &Radix2DFTSmallBatch<KoalaBear>,
    ) -> Vec<KoalaBear> {
        assert_eq!(state.len(), DEGREE, "invalid state size");

        let mut k = vec![KoalaBear::ZERO; DEGREE];
        let mut zero = 0u32;

        struct LimbIterator {
            buf: [u8; KOALA_BEAR_BYTES],
            j: usize,
        }

        impl LimbIterator {
            fn next(&mut self, it: &mut Iter<KoalaBear>) -> Option<u32> {
                if self.j == KOALA_BEAR_BYTES {
                    let l = it.next();
                    if l.is_none() {
                        return None;
                    }

                    self.j = 0;
                    self.buf = l.unwrap().as_canonical_u32().to_le_bytes();
                }

                let bytes: [u8; 2] = [self.buf[self.j], self.buf[self.j + 1]];
                let r = u16::from_le_bytes(bytes) as u32;
                self.j += 2;
                Some(r)
            }
        }

        let mut it = LimbIterator {
            buf: [0u8; KOALA_BEAR_BYTES],
            j: KOALA_BEAR_BYTES,
        };

        for j in 0..DEGREE {
            let l = it.next(vit);
            if l.is_none() {
                break;
            }

            let l = l.unwrap();
            zero = zero | l;
            k[j] = KoalaBear::from_int(l);
        }

        if zero == 0 {
            // means m[i*r.Degree : (i+1)*r.Degree] == [0...0]
            // we can skip this, FFT(0) = 0
            return state;
        }

        let k = dft.dft(k);
        for i in 0..DEGREE {
            state[i] += k[i] * self.ag[pol_id][i];
        }

        state
    }
}

fn derive_random_element_from_seed(seed: u64, i: u64, j: u64) -> KoalaBear {
    let mut buf = [0u8; 3 + 3 * 8];
    buf[..3].copy_from_slice(b"SIS");
    buf[3..11].copy_from_slice(&seed.to_be_bytes()); // bytes 3–10
    buf[11..19].copy_from_slice(&i.to_be_bytes()); // bytes 11–18
    buf[19..27].copy_from_slice(&j.to_be_bytes()); // bytes 19–26

    let res = blake3::hash(&buf);

    let mut first16 = [0u8; 16];
    first16.copy_from_slice(&res.as_bytes()[..16]);
    KoalaBear::from_int(u128::from_be_bytes(first16))
}
