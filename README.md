# Vortex PCS Rust implementation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Pull Requests welcome](https://img.shields.io/badge/PRs-welcome-ff69b4.svg?style=flat-square)](https://github.com/distributed-lab/vortex-rs/issues)
<a href="https://github.com/distributed-lab/vortex-rs">
<img src="https://img.shields.io/github/stars/distributed-lab/vortex-rs?style=social"/>
</a>

⚠️ __Please note - this crypto library has not been audited, so use it at your own risk.__

## Abstract

This implementation leverages Poseidon2 hash function for both hashing columns and Merkle tree. It also leverages
KoalaBear prime field. The field, hash and DFT implementations are taken
from [Plonky3](https://github.com/Plonky3/Plonky3) repository.

The original paper can be found [here](https://eprint.iacr.org/2024/185).

## Usage

Use the following functions to commit/evaluate/open and verify:

1. Commit
   ```rust
      pub fn commit(perm: &PoseidonHash, nb_row: usize, nb_col: usize, w: Vec<Vec<KoalaBear>>) -> (MerkleTree, Vec<Vec<KoalaBear>>)
    ```
2. Eval
   ```rust
      pub fn eval(w: &Vec<Vec<KoalaBear>>, nb_row: usize, nb_col: usize, coin: KoalaBearExt) -> Vec<KoalaBearExt>
    ```
3. Open
   ```rust
      pub fn open(w: &Vec<Vec<KoalaBear>>, w_: &Vec<Vec<KoalaBear>>, nb_row: usize, nb_col: usize, mt: &MerkleTree, beta: KoalaBearExt, column_ids: Vec<usize>) -> OpenProof
    ```
4. Verify proof
   ```rust
      pub fn verify(perm: &PoseidonHash, proof: OpenProof, nb_row: usize, nb_col: usize, root: Digest, y: Vec<KoalaBearExt>, coin: KoalaBearExt)
    ```

The verification function asserts if proof is invalid. The implementation is fixed over KoalaBear field and its
4-degree extension. Check tests in [lib.rs](./src/lib.rs) for more details.

## Benchmarks

Run command:

```shell
cargo test test_rs_encode_matrix --features nightly-features --release -- --show-output
cargo test test_rs_open_matrix --features nightly-features --release -- --show-output
cargo test test_rs_verify_matrix --features nightly-features --release -- --show-output
```

The following benches taken on the M3 Pro 36GB MakBook comparing to the Golang implementation (if changed
to Poseidon2)
from [gnark-crypto](https://github.com/Consensys/gnark-crypto/blob/master/field/koalabear/vortex/prover_test.go#L232)

All test are performed for $2^{19}$ polynomials of $2^{11}$ degree according
to [official benchmarks](https://hackmd.io/@YaoGalteland/SJ1WmzgTJg).

|            | Gnark (sec) | Rust (sec) |
|------------|-------------|------------|
| Commit     | 70-75       | 31-35      |
| Open Proof | 2           | 1.2-1.5    |
| Verify     | 28          | 1.6-1.7    |

## Definition

Imagine we have a list of polynomials $f_0,\dots,f_{k-1}$. We want to commit them and evaluate at the same point at the
same time. We describe each polynomial as vector $a_i$ where

$$
f_i(x) = \sum_{j=0}^{n-1} x^j\cdot a_{i,j}
$$

For simplicity, let's define the interpolation function $Int$ that works as follows (it takes the vector elements as
polynomial coefficients and evaluates it at point $x$):

$$
Int_{a_i}(x) = f_i(x)
$$

Then, we organize these vectors into the matrix $W \in \mathbb{F}^{k\times n}$.

![](./assets/vortex1.png)

We extend our matrix with additional columns by replacing each word $a_i$ with its codeword $a_i'$,
resulting in a matrix $W' \in \mathbb{F}^{k\times m}$, where $m > n$.

Then we hash each column, receiving $m$ values of $h_i$ — we will use these values as our commitment to the polynomials
$f_i$.

![](./assets/vortex1_5.png)

### Open

Given input $r$ from the verifier, the prover responds with values $y_0,\dots,y_{k-1}$ where

$$
y_i = f_i(x)
$$

### Prove & Verification

* The verifier samples a challenge $\beta$
* The prover responds with $u = B\cdot W$, where $B = (1, \beta, \beta^1,\dots,\beta^{k-1})$. Note that naturally, each
  element in $u$ equals to the sum of corresponding polynomials' coefficients over corresponding weight -- polynomial
  $i$ will be multiplied by $\beta^i$.
  ![](./assets/vortex2.png)
* Then, the verifier samples $t$ indexes $q_1,\dots,q_t$ where $q_i \in [m]$
* The prover opens the corresponding columns $s_1,\dots,s_t$ from the matrix $W'$
  ![](./assets/vortex3.png)
* The verifier computes the Reed-Solomon encoding of $u$ named $u'$.
* The verifier checks:
    * $hash(s_i) == h_{q_i}$ for each $i\in [t]$
    * $B\cdot s_i == u'_{q_i}$ → this follows from the linearity of Reed-Solomon
* The verifier checks that $Int_u(x) = B\cdot y$. This check follows from the following observation: $a\cdot Int_c(x) =
  Int_{a\cdot c}(x)$, so:
    * $1 \cdot Int_{w_0}(x) = 1 \cdot y_0$ → $Int_{u_0}(x) = 1 \cdot y_1$
    * $\beta \cdot Int_{w_1}(x) = \beta \cdot y_1$ → $Int_{u_1}(x) = \beta \cdot y_1$
    * $\beta^2 \cdot Int_{w_2}(x) = \beta^2 \cdot y_2$ → $Int_{u_2}(x) = \beta^2 \cdot y_2$
    * etc.

The parameter $t$ is selected according to the security parameters.