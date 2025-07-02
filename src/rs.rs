use p3_dft::{Radix2Dit, Radix2DitParallel, TwoAdicSubgroupDft};
use p3_koala_bear::KoalaBear;

type Dft = Radix2Dit<KoalaBear>;
pub fn encode_reed_solomon(input: Vec<KoalaBear>, rho: usize, dft: &Dft) -> Vec<KoalaBear> {
    let word_sz = input.len();
    let codeword_sz = word_sz * rho;
    let mut result = input;
    result = dft.idft(result);
    result.resize(codeword_sz, KoalaBear::new(0));
    dft.dft(result)
}
