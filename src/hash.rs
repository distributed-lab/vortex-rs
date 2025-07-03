use p3_koala_bear::{KoalaBear, Poseidon2ExternalLayerKoalaBear, Poseidon2InternalLayerKoalaBear};
use p3_poseidon2::Poseidon2;
use p3_symmetric::Permutation;

pub type Digest = [KoalaBear; 8];

pub type PoseidonHash = Poseidon2<
    KoalaBear,
    Poseidon2ExternalLayerKoalaBear<16>,
    Poseidon2InternalLayerKoalaBear<16>,
    16,
    3,
>;

pub fn hash_poseidon2(
    perm: &PoseidonHash,
    left: Digest,
    right: Digest,
) -> Digest {
    let mut input = [
        left[0], left[1], left[2], left[3], left[4], left[5], left[6], left[7], right[0], right[1],
        right[2], right[3], right[4], right[5], right[6], right[7],
    ];

    perm.permute_mut(&mut input);
    [
        input[0], input[1], input[2], input[3], input[4], input[5], input[6], input[7],
    ]
}

#[cfg(test)]
mod tests {
    use p3_koala_bear::{KoalaBear, Poseidon2KoalaBear};
    use rand::prelude::SmallRng;
    use rand::{Rng, SeedableRng};
    use crate::hash::{hash_poseidon2, PoseidonHash};

    #[test]
    fn test_hash_poseidon_works() {
        let mut rng = SmallRng::seed_from_u64(1);
        let perm: PoseidonHash = Poseidon2KoalaBear::new_from_rng_128(&mut rng);

        let mut left = [KoalaBear::new(0); 8];
        for i in 0..8 {
            left[i] = (KoalaBear::new(rng.random()));
        }

        let mut right = [KoalaBear::new(0); 8];
        for i in 0..8 {
            right[i] = KoalaBear::new(rng.random());
        }

        let result = hash_poseidon2(&perm, left, right);

        for i in 0..8 {
            print!("{} ", left[i]);
        }
    }
}
