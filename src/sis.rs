use p3_dft::{Radix2DFTSmallBatch, TwoAdicSubgroupDft};
use p3_field::integers::QuotientMap;
use p3_field::{PrimeCharacteristicRing, PrimeField32, TwoAdicField};
use p3_koala_bear::KoalaBear;
use p3_maybe_rayon::prelude::*;
use p3_util::reverse_slice_index_bits;
use rayon::current_num_threads;
use std::cmp::min;
use std::iter;
use std::slice::Iter;

pub struct RSis {
    pub(crate) a: Vec<Vec<KoalaBear>>,
    pub(crate) ag: Vec<Vec<KoalaBear>>,
    pub(crate) max_nb_elements_to_hash: usize,
    pub(crate) coset: Vec<KoalaBear>,
    pub(crate) ag_shuffled: Vec<Vec<KoalaBear>>,
    pub(crate) twiddles: Vec<Vec<KoalaBear>>,
}

const KOALA_BEAR_BITS: usize = 31;
const KOALA_BEAR_BYTES: usize = 4;

const LOG_TWO_DEGREE: usize = 9;
const LOG_TWO_BOUND: usize = 16;
const LIMB_SIZE: usize = LOG_TWO_BOUND / 8;

pub const DEGREE: usize = 1 << LOG_TWO_DEGREE;

pub const MUL_GENERATOR: u32 = 3;

pub const HASH_STEP: usize = 256;

#[cfg(all(
    feature = "nightly-features",
    target_arch = "x86_64",
    target_feature = "avx512f"
))]
mod ffi {
    include!(concat!(env!("CARGO_MANIFEST_DIR"), "/src/bindings.rs"));
}
#[cfg(all(
    feature = "nightly-features",
    target_arch = "x86_64",
    target_feature = "avx512f"
))]
use ffi::*;

impl RSis {
    pub fn new(seed: u64, max_nb_elements_to_hash: usize) -> Self {
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
            coset: vec![],
            ag_shuffled: vec![],
            twiddles: vec![],
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
            .map_init(
                || Radix2DFTSmallBatch::new(DEGREE),
                |dft, indexes| {
                    let mut res = Vec::with_capacity(indexes.len());
                    for i in indexes {
                        let ag = dft.dft(r.a[i].clone());
                        res.push(ag);
                    }
                    res
                },
            )
            .flatten()
            .collect();

        #[cfg(all(
            feature = "nightly-features",
            target_arch = "x86_64",
            target_feature = "avx512f"
        ))]
        {
            r.ag_shuffled = (0..n)
                .into_par_iter()
                .chunks(current_num_threads())
                .map(|indexes| {
                    let mut res = Vec::with_capacity(indexes.len());

                    for i in indexes {
                        let mut ag_i = r.ag[i].clone();
                        let ag_i_slice = GoSlice {
                            data: ag_i.as_mut_ptr().cast(),
                            len: ag_i.len() as _,
                            cap: ag_i.capacity() as _,
                        };

                        unsafe {
                            SisShuffle_avx512(ag_i_slice);
                        }

                        res.push(ag_i);
                    }

                    res
                })
                .flatten()
                .collect();

            r.coset = vec![KoalaBear::ONE; DEGREE];
            let mul_gen = KoalaBear::from_int(MUL_GENERATOR);
            for i in 1..DEGREE {
                r.coset[i] = r.coset[i - 1] * mul_gen;
            }

            r.twiddles = compute_twiddles(DEGREE);
        }

        r
    }

    pub fn hash(&self, v: &Vec<KoalaBear>, dft: &Radix2DFTSmallBatch<KoalaBear>) -> Vec<KoalaBear> {
        assert!(
            v.len() <= self.max_nb_elements_to_hash,
            "can't hash more than configured elements with params provided in constructor"
        );
        let mut res = vec![KoalaBear::ZERO; DEGREE];

        #[cfg(all(
            feature = "nightly-features",
            target_arch = "x86_64",
            target_feature = "avx512f"
        ))]
        {
            let mut pol_id = 0;

            let mut twiddles = convert_2d_arr_to_go(&self.twiddles);

            for j in (0..v.len()).step_by(HASH_STEP) {
                let start = j;
                let end = min(j + HASH_STEP, v.len());
                let mut v_ = v[start..end].to_vec();

                let v__slice = GoSlice {
                    data: v_.as_mut_ptr().cast(),
                    len: v_.len() as _,
                    cap: v_.capacity() as _,
                };

                let coset_slice = GoSlice {
                    data: self.coset.as_ptr() as *mut _,
                    len: self.coset.len() as _,
                    cap: self.coset.capacity() as _,
                };

                let twiddles_slice = GoSlice {
                    data: twiddles.as_mut_ptr().cast(),
                    len: twiddles.len() as _,
                    cap: twiddles.capacity() as _,
                };

                let ag_shuffled_slice = GoSlice {
                    data: self.ag_shuffled[pol_id].as_ptr() as *mut _,
                    len: self.ag_shuffled[pol_id].len() as _,
                    cap: self.ag_shuffled[pol_id].capacity() as _,
                };

                let res_slice = GoSlice {
                    data: res.as_mut_ptr().cast(),
                    len: res.len() as _,
                    cap: res.capacity() as _,
                };

                unsafe {
                    Sis512_16_avx512(
                        v__slice,
                        coset_slice,
                        twiddles_slice,
                        ag_shuffled_slice,
                        res_slice,
                    );
                }
                pol_id += 1;
            }

            let res_slice = GoSlice {
                data: res.as_mut_ptr().cast(),
                len: res.len() as _,
                cap: res.capacity() as _,
            };

            unsafe {
                SisUnshuffle_avx512(res_slice);
            }
        }

        #[cfg(all(not(target_feature = "avx512f")))]
        {
            for i in 0..self.ag.len() {
                res = self.inner_hash(res, &mut v.iter(), i, &dft);
            }
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

fn roots_of_unity_table(n: usize) -> Vec<Vec<KoalaBear>> {
    let lg_n = n.ilog2() as usize;
    let generator = KoalaBear::two_adic_generator(lg_n);

    let mut t = vec![vec![]; lg_n];
    t[0] = vec![KoalaBear::ONE; 1 + (1 << (lg_n - 1))];
    for i in 1..t[0].len() {
        t[0][i] = t[0][i - 1] * generator;
    }

    for i in 1..lg_n {
        t[i] = vec![KoalaBear::ZERO; 1 + (1 << (lg_n - i - 1))];
        let mut k = 0;
        for j in 0..t[i].len() {
            t[i][j] = t[0][k];
            k += 1 << i;
        }
    }

    t
}

fn compute_twiddles(fft_len: usize) -> Vec<Vec<KoalaBear>> {
    // roots_of_unity_table(fft_len) returns a vector of twiddles of length log_2(fft_len).
    let mut new_twiddles = roots_of_unity_table(fft_len);

    new_twiddles.iter_mut().for_each(|ts| {
        let ln = ts.len();
        reverse_slice_index_bits(ts[0..ln - 1].as_mut());
    });

    new_twiddles
}
