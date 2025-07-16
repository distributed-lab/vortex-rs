package main

import "C"

import "sis-binding/asm"

// Re‑export the Element type so callers see the same layout.
type Element = asm.Element

// -----------------------------------------------------------------------------
// Wrappers that the Go toolchain will turn into C symbols
// -----------------------------------------------------------------------------

//export SisShuffle_avx512
func SisShuffle_avx512(a []Element) {
	asm.SisShuffle_avx512(a)
}

//export SisUnshuffle_avx512
func SisUnshuffle_avx512(a []Element) {
	asm.SisUnshuffle_avx512(a)
}

//export Sis512_16_avx512
func Sis512_16_avx512(k256, cosets []Element,
	twiddles [][]Element, rag, res []Element) {
	asm.Sis512_16_avx512(k256, cosets, twiddles, rag, res)
}

// main() is mandatory in every package built with -buildmode=c-shared
func main() {}
