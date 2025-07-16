use crate::KoalaBearExt;
use p3_dft::{Radix2DFTSmallBatch, Radix2DitParallel, TwoAdicSubgroupDft};
use p3_field::PrimeCharacteristicRing;
use p3_koala_bear::KoalaBear;

pub fn encode_reed_solomon(
    input: Vec<KoalaBear>,
    rho: usize,
    dft: &Radix2DFTSmallBatch<KoalaBear>,
) -> Vec<KoalaBear> {
    let word_sz = input.len();
    let codeword_sz = word_sz * rho;
    let mut result = input;
    result = dft.idft(result);
    result.resize(codeword_sz, KoalaBear::ZERO);
    dft.dft(result)
}

pub fn encode_reed_solomon_ext(
    input: Vec<KoalaBearExt>,
    rho: usize,
    dft: &Radix2DitParallel<KoalaBearExt>,
) -> Vec<KoalaBearExt> {
    let word_sz = input.len();
    let codeword_sz = word_sz * rho;
    let mut result = input;
    result = dft.idft(result);
    result.resize(codeword_sz, KoalaBearExt::ZERO);
    dft.dft(result)
}
