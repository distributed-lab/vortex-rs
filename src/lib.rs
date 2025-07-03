use crate::hash::{PoseidonHash, hash_poseidon2, Digest};
use crate::merkle_tree::MerkleTree;
use crate::rs::encode_reed_solomon;
use p3_dft::{Radix2Dit, Radix2DitParallel};
use p3_field::PrimeCharacteristicRing;
use p3_koala_bear::{KoalaBear, Poseidon2KoalaBear};
use p3_maybe_rayon::prelude::*;
use p3_symmetric::Permutation;
use rand::SeedableRng;
use rand::prelude::SmallRng;
use rayon::current_num_threads;

pub mod hash;
pub mod merkle_tree;
pub mod rs;

pub fn commit(
    perm: &PoseidonHash,
    nb_row: usize,
    nb_col: usize,
    w: Vec<Vec<KoalaBear>>,
) -> MerkleTree {
    let w_: Vec<Vec<KoalaBear>> = w
        .into_par_iter().chunks(current_num_threads())
        .map(|chunk| {
            let mut res = Vec::with_capacity(chunk.len());
            for wi in chunk {
                res.push(encode_reed_solomon(wi, 2, &Radix2Dit::default()))
            }
            res
        })
        .flatten()
        .collect();

    let hash: Vec<Digest> = (0..nb_col)
        .into_par_iter()
        .chunks(current_num_threads())
        .map(|indexes| {
            let mut res = Vec::with_capacity(indexes.len());

            for i in indexes {
                let mut left = [KoalaBear::ZERO; 8];
                let mut right = [KoalaBear::ZERO; 8];
                for j in (0..nb_row).step_by(8) {
                    right[0] = w_[j][i];
                    right[1] = w_[j + 1][i];
                    right[2] = w_[j + 2][i];
                    right[3] = w_[j + 3][i];
                    right[4] = w_[j + 4][i];
                    right[5] = w_[j + 5][i];
                    right[6] = w_[j + 6][i];
                    right[7] = w_[j + 7][i];
                    left = hash_poseidon2(&perm, left, right);
                }
                res.push(left)
            }

            res
        })
        .flatten()
        .collect();

    MerkleTree::build(perm, hash)
}

#[cfg(test)]
mod tests {
    use super::*;
    use p3_dft::{Radix2DitParallel, TwoAdicSubgroupDft};
    use p3_koala_bear::KoalaBear;
    use rand::prelude::SmallRng;
    use rand::{Rng, SeedableRng};
    use std::cmp::Ordering;

    #[test]
    fn fft_works() {
        let dft: Radix2DitParallel<KoalaBear> = Radix2DitParallel::default();
        let mut rng = SmallRng::seed_from_u64(1);

        let sz = 1 << 10;

        let mut v = vec![];
        for i in 0..sz {
            v.push(KoalaBear::new(rng.random()));
        }

        let res = dft.dft(v.clone());
        let vres = dft.idft(res);

        for i in 0..vres.len() {
            assert_eq!(v[i].cmp(&vres[i]), Ordering::Equal)
        }
    }

    #[test]
    fn test_rs_encode_works() {
        let dft: Radix2DitParallel<KoalaBear> = Radix2DitParallel::default();
        let mut rng = SmallRng::seed_from_u64(1);

        const RHO: usize = 2;
        let sz = 1 << 10;
        //let shift = KoalaBear::two_adic_generator(11); // sz in bits + rho in bits

        let mut input = vec![];
        for i in 0..sz {
            input.push(KoalaBear::new(rng.random()));
        }

        let mut result = input.clone();
        result = dft.idft(input.clone());
        while result.len() < input.len() * RHO {
            result.push(KoalaBear::new(0));
        }
        result = dft.dft(result);
    }

    fn new_row(sz: usize, rng: &mut SmallRng) -> Vec<KoalaBear> {
        let mut v = vec![];
        for i in 0..sz {
            v.push(KoalaBear::new(rng.random()))
        }
        v
    }

    #[test]
    fn test_rs_encode_matrix() {
        let nb_col = 1 << 11;
        let nb_row = 1 << 19;

        let mut rng = SmallRng::seed_from_u64(1);
        let perm: PoseidonHash = Poseidon2KoalaBear::new_from_rng_128(&mut rng);

        let mut w = vec![];
        for i in 0..nb_row {
            w.push(new_row(nb_col, &mut rng))
        }

        commit(&perm, nb_row, nb_col, w);
    }
}
