use crate::PoseidonHash;
use crate::hash::{Digest, hash_poseidon2};
use p3_koala_bear::KoalaBear;
use p3_maybe_rayon::prelude::*;
use rayon::current_num_threads;

pub struct MerkleTree {
    pub(crate) digest_layers: Vec<Vec<Digest>>,
}

impl MerkleTree {
    pub fn build(perm: &PoseidonHash, leaves: Vec<Digest>) -> Self {
        let depth = leaves.len().ilog2() as usize;
        let mut layers = Vec::with_capacity(depth + 1);
        layers.push(leaves);

        for i in 0..depth {
            let prev = &layers[i];

            if prev.len() / 2 >= 64 {
                let next: Vec<Digest> = (0..prev.len() / 2)
                    .into_par_iter()
                    .chunks(current_num_threads())
                    .map(|indexes| {
                        let mut res = vec![Digest::default(); indexes.len()];
                        for (idx, j) in indexes.into_iter().enumerate() {
                            res[idx] = hash_poseidon2(perm, prev[2 * j], prev[2 * j + 1]);
                        }
                        res
                    })
                    .flatten()
                    .collect();

                layers.push(next);
            } else {
                let mut next =  vec![Digest::default(); prev.len() / 2];
                for j in 0..(prev.len() / 2) {
                    next[j] = hash_poseidon2(perm, prev[2 * j], prev[2 * j + 1]);
                }

                layers.push(next);
            }
        }

        assert_eq!(layers.last().unwrap().len(), 1);

        MerkleTree {
            digest_layers: layers,
        }
    }

    pub fn open(&self, mut index: usize) -> Vec<Digest> {
        let mut proof = Vec::with_capacity(self.digest_layers.len());

        for level in (0..self.digest_layers.len() - 1) {
            let sibling = index ^ 1;
            proof.push(self.digest_layers[level][sibling]);
            index >>= 1;
        }

        proof
    }

    pub fn root(&self) -> Digest {
        self.digest_layers.last().unwrap()[0]
    }
}

pub fn verify_merkle_proof(
    mut index: usize,
    leaf: Digest,
    root: Digest,
    proof: Vec<Digest>,
    perm: &PoseidonHash,
) -> bool {
    let mut current = leaf;

    for sibling in proof {
        let (left, right) = if index & 1 == 0 {
            (current, sibling)
        } else {
            (sibling, current)
        };
        current = hash_poseidon2(perm, left, right);
        index >>= 1;
    }

    current == root
}

#[cfg(test)]
mod tests {
    use crate::hash::{Digest, PoseidonHash};
    use crate::merkle_tree::{MerkleTree, verify_merkle_proof};
    use p3_koala_bear::{KoalaBear, Poseidon2KoalaBear};
    use rand::rngs::SmallRng;
    use rand::{Rng, SeedableRng};

    fn random_leaf(rng: &mut SmallRng) -> Digest {
        let mut leaf = [KoalaBear::default(); 8];
        for i in 0..8 {
            leaf[i] = KoalaBear::new(rng.random());
        }
        leaf
    }

    #[test]
    fn test_merkle_tree() {
        let mut rng = SmallRng::seed_from_u64(1);
        let perm: PoseidonHash = Poseidon2KoalaBear::new_from_rng_128(&mut rng);

        // Generate a power-of-two number of leaves (e.g., 8)
        let leaves: Vec<Digest> = (0..256).map(|_| random_leaf(&mut rng)).collect();

        // Build Merkle Tree
        let tree = MerkleTree::build(&perm, leaves.clone());
        let root = tree.digest_layers.last().unwrap()[0];

        // Check proof verification for every leaf
        for (i, leaf) in leaves.iter().enumerate() {
            let proof = tree.open(i);
            let valid = verify_merkle_proof(i, *leaf, root, proof, &perm);
            assert!(valid, "Merkle proof failed for index {}", i);
        }
    }
}
