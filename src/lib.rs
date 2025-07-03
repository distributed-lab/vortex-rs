use crate::hash::{Digest, PoseidonHash, hash_poseidon2};
use crate::merkle_tree::{MerkleTree, verify_merkle_proof};
use crate::rs::{encode_reed_solomon, encode_reed_solomon_ext};
use p3_dft::{Radix2DFTSmallBatch, Radix2Dit, Radix2DitParallel};
use p3_field::PrimeCharacteristicRing;
use p3_field::extension::BinomialExtensionField;
use p3_koala_bear::{KoalaBear, Poseidon2KoalaBear};
use p3_maybe_rayon::prelude::*;
use p3_symmetric::Permutation;
use rand::SeedableRng;
use rand::prelude::SmallRng;
use rayon::current_num_threads;

pub mod hash;
pub mod merkle_tree;
pub mod rs;

pub type KoalaBearExt = BinomialExtensionField<KoalaBear, 4>;

pub const RS_RATE: usize = 2;

pub fn commit(
    perm: &PoseidonHash,
    nb_row: usize,
    nb_col: usize,
    w: Vec<Vec<KoalaBear>>,
) -> (MerkleTree, Vec<Vec<KoalaBear>>) {
    let w_: Vec<Vec<KoalaBear>> = w
        .into_par_iter()
        .chunks(current_num_threads())
        .map(|chunk| {
            let mut res = Vec::with_capacity(chunk.len());
            let dft = Radix2DFTSmallBatch::default();
            for wi in chunk {
                res.push(encode_reed_solomon(wi, RS_RATE, &dft))
            }
            res
        })
        .flatten()
        .collect();

    let hash: Vec<Digest> = (0..nb_col * RS_RATE)
        .into_par_iter()
        .chunks(current_num_threads())
        .map(|indexes| {
            let mut res = vec![Digest::default(); indexes.len()];

            for (idx, i) in indexes.into_iter().enumerate() {
                let mut buf = Digest::default();
                for j in (0..nb_row).step_by(8) {
                    buf[0] = w_[j][i];
                    buf[1] = w_[j + 1][i];
                    buf[2] = w_[j + 2][i];
                    buf[3] = w_[j + 3][i];
                    buf[4] = w_[j + 4][i];
                    buf[5] = w_[j + 5][i];
                    buf[6] = w_[j + 6][i];
                    buf[7] = w_[j + 7][i];
                    res[idx] = hash_poseidon2(&perm, res[idx], buf);
                }
            }

            res
        })
        .flatten()
        .collect();

    (MerkleTree::build(perm, hash), w_)
}

pub struct OpenProof {
    pub columns: Vec<Vec<KoalaBear>>,
    pub merkle_proofs: Vec<Vec<Digest>>,
    pub lin_comb: Vec<KoalaBearExt>,
    pub column_ids: Vec<usize>,
    pub beta: KoalaBearExt,
}

fn eval_lin_comb(
    w: &Vec<Vec<KoalaBear>>,
    nb_row: usize,
    nb_col: usize,
    beta: KoalaBearExt,
) -> Vec<KoalaBearExt> {
    let collected: Vec<Vec<KoalaBearExt>> = (0..nb_row)
        .into_par_iter()
        .chunks(current_num_threads())
        .map(|indexes| {
            let mut res = vec![KoalaBearExt::ZERO; nb_col];

            let mut multiplier = beta.exp_u64(indexes[0] as u64);
            for i in indexes {
                for j in 0..nb_col {
                    res[j] += KoalaBearExt::from(w[i][j]) * multiplier;
                }
                multiplier *= beta;
            }
            res
        })
        .collect();

    let result = (0..nb_col)
        .into_par_iter()
        .chunks(current_num_threads())
        .map(|indexes| {
            let mut result = vec![KoalaBearExt::ZERO; indexes.len()];
            for (idx, i) in indexes.into_iter().enumerate() {
                for j in 0..collected.len() {
                    result[idx] += collected[j][i];
                }
            }
            result
        })
        .flatten()
        .collect();

    result
}
pub fn eval(
    w: &Vec<Vec<KoalaBear>>,
    nb_row: usize,
    nb_col: usize,
    coin: KoalaBearExt,
) -> Vec<KoalaBearExt> {
    (0..nb_row)
        .into_par_iter()
        .chunks(current_num_threads())
        .map(|indexes| {
            let mut res = vec![KoalaBearExt::ZERO; indexes.len()];
            let mut x = vec![KoalaBearExt::ONE; nb_col];
            for i in 1..nb_col {
                x[i] = x[i - 1] * coin;
            }

            for (idx, i) in indexes.into_iter().enumerate() {
                for j in 0..nb_col {
                    res[idx] += KoalaBearExt::from(w[i][j]) * x[j];
                }
            }

            res
        })
        .flatten()
        .collect()
}

pub fn open(
    w: &Vec<Vec<KoalaBear>>,
    w_: &Vec<Vec<KoalaBear>>,
    nb_row: usize,
    nb_col: usize,
    mt: &MerkleTree,
    beta: KoalaBearExt,
    column_ids: Vec<usize>,
) -> OpenProof {
    let open_columns = column_ids
        .par_iter()
        .chunks(current_num_threads())
        .map(|indexes| {
            let mut res = vec![vec![KoalaBear::ZERO; nb_row]; indexes.len()];

            for (idx, col) in indexes.into_iter().enumerate() {
                for i in 0..nb_row {
                    res[idx][i] = w_[i][*col]
                }
            }
            res
        })
        .flatten()
        .collect();

    OpenProof {
        columns: open_columns,
        merkle_proofs: column_ids.iter().map(|col| mt.open(*col)).collect(),
        lin_comb: eval_lin_comb(w, nb_row, nb_col, beta),
        column_ids,
        beta,
    }
}
pub fn verify(
    perm: &PoseidonHash,
    proof: OpenProof,
    nb_row: usize,
    nb_col: usize,
    root: Digest,
    y: Vec<KoalaBearExt>,
    coin: KoalaBearExt,
) {
    let mut betas = vec![KoalaBearExt::ONE; nb_row];
    for i in 1..nb_row {
        betas[i] = betas[i - 1] * proof.beta;
    }

    let mut x = vec![KoalaBearExt::ONE; nb_col];
    for i in 1..nb_col {
        x[i] = x[i - 1] * coin;
    }

    let ux: KoalaBearExt = (0..nb_col)
        .into_par_iter()
        .map(|i| proof.lin_comb[i] * x[i])
        .sum();

    let beta_y: KoalaBearExt = (0..nb_row).into_par_iter().map(|i| y[i] * betas[i]).sum();
    assert_eq!(beta_y, ux, "failed to verify evaluation");

    let dft = Radix2DFTSmallBatch::default();
    let u_ = encode_reed_solomon_ext(proof.lin_comb, RS_RATE, &dft);

    proof
        .columns
        .into_par_iter()
        .enumerate()
        .chunks(current_num_threads())
        .for_each(|chunks| {
            for (idx, column) in chunks {
                let mut column_hash = Digest::default();
                let mut buf = Digest::default();

                for i in (0..column.len()).step_by(8) {
                    buf[0] = column[i];
                    buf[1] = column[i + 1];
                    buf[2] = column[i + 2];
                    buf[3] = column[i + 3];
                    buf[4] = column[i + 4];
                    buf[5] = column[i + 5];
                    buf[6] = column[i + 6];
                    buf[7] = column[i + 7];
                    column_hash = hash_poseidon2(&perm, column_hash, buf);
                }

                assert!(
                    verify_merkle_proof(
                        proof.column_ids[idx],
                        column_hash,
                        root,
                        proof.merkle_proofs[idx].clone(),
                        perm,
                    ),
                    "Failed to verify merkle proof"
                );

                let mut beta_column = KoalaBearExt::ZERO;
                for i in 0..nb_row {
                    beta_column += KoalaBearExt::from(column[i]) * betas[i];
                }

                assert_eq!(
                    beta_column, u_[proof.column_ids[idx]],
                    "failed to verify RS linearity"
                )
            }
        });
}

#[cfg(test)]
mod tests {
    use super::*;
    use p3_dft::{Radix2DitParallel, TwoAdicSubgroupDft};
    use p3_koala_bear::KoalaBear;
    use rand::distr::{StandardUniform, Uniform};
    use rand::prelude::SmallRng;
    use rand::seq::IteratorRandom;
    use rand::{Rng, SeedableRng};
    use std::cmp::Ordering;
    use std::time::{Instant, SystemTime, UNIX_EPOCH};

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

        let start = Instant::now();

        commit(&perm, nb_row, nb_col, w);

        let elapsed = start.elapsed();
        println!("Commit: {:?}", elapsed);
    }

    #[test]
    fn test_rs_open_matrix() {
        let nb_col = 1 << 11;
        let nb_row = 1 << 19;
        let nb_columns_to_open = 256;

        let mut rng = SmallRng::seed_from_u64(1);
        let perm: PoseidonHash = Poseidon2KoalaBear::new_from_rng_128(&mut rng);

        let mut w = vec![];
        for i in 0..nb_row {
            w.push(new_row(nb_col, &mut rng))
        }

        let start = Instant::now();

        let (mt, w_) = commit(&perm, nb_row, nb_col, w.clone());

        let elapsed = start.elapsed();
        println!("Commit: {:?}", elapsed);

        let beta: KoalaBearExt = rng.sample(StandardUniform {});
        println!("{:?}", beta);
        let columns: Vec<usize> =
            (0..nb_col * RS_RATE).choose_multiple(&mut rng, nb_columns_to_open);

        let start = Instant::now();

        open(&w, &w_, nb_row, nb_col, &mt, beta, columns);

        let elapsed = start.elapsed();
        println!("Open: {:?}", elapsed);
    }

    #[test]
    fn test_rs_verify_matrix() {
        let nb_col = 1 << 11;
        let nb_row = 1 << 19;
        let nb_columns_to_open = 256;

        let mut rng = SmallRng::seed_from_u64(
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_secs(),
        );
        let perm: PoseidonHash = Poseidon2KoalaBear::new_from_rng_128(&mut rng);

        let mut w = vec![];
        for i in 0..nb_row {
            w.push(new_row(nb_col, &mut rng))
        }

        // ------------------------
        let start = Instant::now();

        let (mt, w_) = commit(&perm, nb_row, nb_col, w.clone());

        let elapsed = start.elapsed();
        println!("Commit: {:?}", elapsed);

        // ------------------------

        let beta: KoalaBearExt = rng.sample(StandardUniform {});
        let columns: Vec<usize> =
            (0..nb_col * RS_RATE).choose_multiple(&mut rng, nb_columns_to_open);

        // ------------------------
        let start = Instant::now();

        let coin: KoalaBearExt = rng.sample(StandardUniform {});
        let y = eval(&w, nb_row, nb_col, coin);

        let elapsed = start.elapsed();
        println!("Eval: {:?}", elapsed);

        // ------------------------
        let start = Instant::now();

        let proof = open(&w, &w_, nb_row, nb_col, &mt, beta, columns);

        let elapsed = start.elapsed();
        println!("Open: {:?}", elapsed);

        // ------------------------
        let start = Instant::now();
        verify(&perm, proof, nb_row, nb_col, mt.root(), y, coin);

        let elapsed = start.elapsed();
        println!("Verify: {:?}", elapsed);
    }
}
