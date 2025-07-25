//go:build linux && amd64

// Code generated by gnark-crypto/generator. DO NOT EDIT.
// Refer to the generator for more documentation.
// Some sub-functions are derived from Plonky3:
// https://github.com/Plonky3/Plonky3/blob/36e619f3c6526ee86e2e5639a24b3224e1c1700f/monty-31/src/x86_64_avx512/packing.rs#L319

#include "textflag.h"
#include "funcdata.h"
#include "go_asm.h"

#define BUTTERFLYD1Q(in0, in1, in2, in3, in4) \
	VPADDD  in0, in1, in3 \
	VPSUBD  in1, in0, in1 \
	VPSUBD  in2, in3, in0 \
	VPMINUD in3, in0, in0 \
	VPADDD  in2, in1, in4 \
	VPMINUD in4, in1, in1 \

#define BUTTERFLYD2Q(in0, in1, in2, in3, in4) \
	VPSUBD  in1, in0, in4 \
	VPADDD  in0, in1, in3 \
	VPADDD  in2, in4, in1 \
	VPSUBD  in2, in3, in0 \
	VPMINUD in3, in0, in0 \

#define BUTTERFLYD2Q2Q(in0, in1, in2, in3) \
	VPSUBD in1, in0, in3 \
	VPADDD in0, in1, in0 \
	VPADDD in2, in3, in1 \

#define MULD(in0, in1, in2, in3, in4, in5, in6, in7, in8, in9) \
	VPSRLQ    $32, in0, in2 \
	VPSRLQ    $32, in1, in3 \
	VPMULUDQ  in0, in1, in4 \
	VPMULUDQ  in2, in3, in5 \
	VPMULUDQ  in4, in9, in6 \
	VPMULUDQ  in5, in9, in7 \
	VPMULUDQ  in6, in8, in6 \
	VPADDQ    in4, in6, in4 \
	VPMULUDQ  in7, in8, in7 \
	VPADDQ    in5, in7, in5 \
	VMOVSHDUP in4, K3, in5  \
	VPSUBD    in8, in5, in7 \
	VPMINUD   in5, in7, in0 \

#define PERMUTE8X8(in0, in1, in2) \
	VSHUFI64X2 $0x000000000000004e, in1, in0, in2 \
	VPBLENDMQ  in0, in2, K1, in0                  \
	VPBLENDMQ  in2, in1, K1, in1                  \

#define PERMUTE4X4(in0, in1, in2, in3) \
	VMOVDQA64 in2, in3          \
	VPERMI2Q  in1, in0, in3     \
	VPBLENDMQ in0, in3, K2, in0 \
	VPBLENDMQ in3, in1, K2, in1 \

#define PERMUTE2X2(in0, in1, in2) \
	VSHUFPD   $0x0000000000000055, in1, in0, in2 \
	VPBLENDMQ in0, in2, K3, in0                  \
	VPBLENDMQ in2, in1, K3, in1                  \

#define PERMUTE1X1(in0, in1, in2) \
	VPSHRDQ   $32, in1, in0, in2 \
	VPBLENDMD in0, in2, K3, in0  \
	VPBLENDMD in2, in1, K3, in1  \

#define LOAD_Q(in0, in1) \
	MOVD         $const_q, AX       \
	VPBROADCASTD AX, in0            \
	MOVD         $const_qInvNeg, AX \
	VPBROADCASTD AX, in1            \

#define LOAD_MASKS() \
	MOVQ  $0x0000000000000f0f, AX \
	KMOVQ AX, K1                  \
	MOVQ  $0x0000000000000033, AX \
	KMOVQ AX, K2                  \
	MOVQ  $0x0000000000005555, AX \
	KMOVD AX, K3                  \

#define BUTTERFLY_MULD(in0, in1, in2, in3, in4, in5, in6, in7, in8, in9, in10, in11, in12, in13, in14) \
BUTTERFLYD2Q(in0, in1, in2, in3, in4)                       \
MULD(in5, in6, in7, in8, in9, in10, in11, in12, in13, in14) \

TEXT ·Sis512_16_avx512(SB), $1024-120
	// refer to the code generator for comments and documentation.
	LOAD_Q(Z0, Z1)
	LOAD_MASKS()
	MOVQ k256+0(FP), R15
	MOVQ cosets+24(FP), CX
	MOVQ twiddles+48(FP), SI
	MOVQ rag+72(FP), R8
	MOVQ res+96(FP), R9
	MOVQ 0(SI), DI               // twiddles[0]
	MOVQ R15, DX
	MOVQ CX, BX
	ADDQ $0x0000000000000200, DX
	ADDQ $0x0000000000000400, BX

#define FROMMONTGOMERY(in0, in1, in2, in3, in4, in5, in6, in7, in8, in9, in10, in11) \
	VPSRLQ    $32, in0, in3     \
	VPSRLQ    $32, in1, in7     \
	VPMULUDQ  in0, in11, in4    \
	VPMULUDQ  in3, in11, in5    \
	VPMULUDQ  in1, in11, in8    \
	VPMULUDQ  in7, in11, in9    \
	VPMULUDQ  in4, in10, in4    \
	VPMULUDQ  in8, in10, in8    \
	VPMULUDQ  in5, in10, in5    \
	VPMULUDQ  in9, in10, in9    \
	VPANDD.Z  in0, in0, K3, in2 \
	VPANDD.Z  in1, in1, K3, in6 \
	VPADDQ    in2, in4, in2     \
	VPADDQ    in6, in8, in6     \
	VPADDQ    in3, in5, in3     \
	VPADDQ    in7, in9, in7     \
	VMOVSHDUP in6, K3, in7      \
	VMOVSHDUP in2, K3, in3      \
	VPSUBD    in10, in3, in5    \
	VPSUBD    in10, in7, in9    \
	VPMINUD   in3, in5, in0     \
	VPMINUD   in7, in9, in1     \

	VMOVDQU32     0(R15), Z16
	VMOVDQU32     0(DX), Z14
	FROMMONTGOMERY(Z16, Z14, Z4, Z5, Z6, Z7, Z8, Z9, Z10, Z11, Z0, Z1)
	VEXTRACTI64X4 $1, Z16, Y17
	VPMOVZXWD     Y16, Z16
	VPMOVZXWD     Y17, Z17
	VMOVDQU32     0(CX), Z12
	MULD(Z16, Z12, Z2, Z3, Z4, Z5, Z6, Z7, Z0, Z1)
	VMOVDQU32     64(CX), Z13
	MULD(Z17, Z13, Z8, Z9, Z10, Z11, Z6, Z7, Z0, Z1)
	VEXTRACTI64X4 $1, Z14, Y15
	VPMOVZXWD     Y14, Z14
	VPMOVZXWD     Y15, Z15
	VMOVDQU32     0(BX), Z12
	MULD(Z14, Z12, Z2, Z3, Z4, Z5, Z6, Z7, Z0, Z1)
	VMOVDQU32     64(BX), Z13
	MULD(Z15, Z13, Z8, Z9, Z10, Z11, Z6, Z7, Z0, Z1)
	BUTTERFLYD2Q(Z16, Z14, Z0, Z4, Z6)
	BUTTERFLYD2Q(Z17, Z15, Z0, Z5, Z7)
	VMOVDQU32     0(DI), Z12
	MULD(Z14, Z12, Z2, Z3, Z4, Z5, Z6, Z7, Z0, Z1)
	VMOVDQU32     64(DI), Z13
	MULD(Z15, Z13, Z8, Z9, Z10, Z11, Z6, Z7, Z0, Z1)
	VMOVDQU32     Z14, 0(SP)
	VMOVDQU32     Z15, 64(SP)
	VMOVDQU32     64(R15), Z18
	VMOVDQU32     64(DX), Z14
	FROMMONTGOMERY(Z18, Z14, Z4, Z5, Z6, Z7, Z8, Z9, Z10, Z11, Z0, Z1)
	VEXTRACTI64X4 $1, Z18, Y19
	VPMOVZXWD     Y18, Z18
	VPMOVZXWD     Y19, Z19
	VMOVDQU32     128(CX), Z12
	MULD(Z18, Z12, Z2, Z3, Z4, Z5, Z6, Z7, Z0, Z1)
	VMOVDQU32     192(CX), Z13
	MULD(Z19, Z13, Z8, Z9, Z10, Z11, Z6, Z7, Z0, Z1)
	VEXTRACTI64X4 $1, Z14, Y15
	VPMOVZXWD     Y14, Z14
	VPMOVZXWD     Y15, Z15
	VMOVDQU32     128(BX), Z12
	MULD(Z14, Z12, Z2, Z3, Z4, Z5, Z6, Z7, Z0, Z1)
	VMOVDQU32     192(BX), Z13
	MULD(Z15, Z13, Z8, Z9, Z10, Z11, Z6, Z7, Z0, Z1)
	BUTTERFLYD2Q(Z18, Z14, Z0, Z4, Z6)
	BUTTERFLYD2Q(Z19, Z15, Z0, Z5, Z7)
	VMOVDQU32     128(DI), Z12
	MULD(Z14, Z12, Z2, Z3, Z4, Z5, Z6, Z7, Z0, Z1)
	VMOVDQU32     192(DI), Z13
	MULD(Z15, Z13, Z8, Z9, Z10, Z11, Z6, Z7, Z0, Z1)
	VMOVDQU32     Z14, 128(SP)
	VMOVDQU32     Z15, 192(SP)
	VMOVDQU32     128(R15), Z20
	VMOVDQU32     128(DX), Z14
	FROMMONTGOMERY(Z20, Z14, Z4, Z5, Z6, Z7, Z8, Z9, Z10, Z11, Z0, Z1)
	VEXTRACTI64X4 $1, Z20, Y21
	VPMOVZXWD     Y20, Z20
	VPMOVZXWD     Y21, Z21
	VMOVDQU32     256(CX), Z12
	MULD(Z20, Z12, Z2, Z3, Z4, Z5, Z6, Z7, Z0, Z1)
	VMOVDQU32     320(CX), Z13
	MULD(Z21, Z13, Z8, Z9, Z10, Z11, Z6, Z7, Z0, Z1)
	VEXTRACTI64X4 $1, Z14, Y15
	VPMOVZXWD     Y14, Z14
	VPMOVZXWD     Y15, Z15
	VMOVDQU32     256(BX), Z12
	MULD(Z14, Z12, Z2, Z3, Z4, Z5, Z6, Z7, Z0, Z1)
	VMOVDQU32     320(BX), Z13
	MULD(Z15, Z13, Z8, Z9, Z10, Z11, Z6, Z7, Z0, Z1)
	BUTTERFLYD2Q(Z20, Z14, Z0, Z4, Z6)
	BUTTERFLYD2Q(Z21, Z15, Z0, Z5, Z7)
	VMOVDQU32     256(DI), Z12
	MULD(Z14, Z12, Z2, Z3, Z4, Z5, Z6, Z7, Z0, Z1)
	VMOVDQU32     320(DI), Z13
	MULD(Z15, Z13, Z8, Z9, Z10, Z11, Z6, Z7, Z0, Z1)
	VMOVDQU32     Z14, 256(SP)
	VMOVDQU32     Z15, 320(SP)
	VMOVDQU32     192(R15), Z22
	VMOVDQU32     192(DX), Z14
	FROMMONTGOMERY(Z22, Z14, Z4, Z5, Z6, Z7, Z8, Z9, Z10, Z11, Z0, Z1)
	VEXTRACTI64X4 $1, Z22, Y23
	VPMOVZXWD     Y22, Z22
	VPMOVZXWD     Y23, Z23
	VMOVDQU32     384(CX), Z12
	MULD(Z22, Z12, Z2, Z3, Z4, Z5, Z6, Z7, Z0, Z1)
	VMOVDQU32     448(CX), Z13
	MULD(Z23, Z13, Z8, Z9, Z10, Z11, Z6, Z7, Z0, Z1)
	VEXTRACTI64X4 $1, Z14, Y15
	VPMOVZXWD     Y14, Z14
	VPMOVZXWD     Y15, Z15
	VMOVDQU32     384(BX), Z12
	MULD(Z14, Z12, Z2, Z3, Z4, Z5, Z6, Z7, Z0, Z1)
	VMOVDQU32     448(BX), Z13
	MULD(Z15, Z13, Z8, Z9, Z10, Z11, Z6, Z7, Z0, Z1)
	BUTTERFLYD2Q(Z22, Z14, Z0, Z4, Z6)
	BUTTERFLYD2Q(Z23, Z15, Z0, Z5, Z7)
	VMOVDQU32     384(DI), Z12
	MULD(Z14, Z12, Z2, Z3, Z4, Z5, Z6, Z7, Z0, Z1)
	VMOVDQU32     448(DI), Z13
	MULD(Z15, Z13, Z8, Z9, Z10, Z11, Z6, Z7, Z0, Z1)
	VMOVDQU32     Z14, 384(SP)
	VMOVDQU32     Z15, 448(SP)
	VMOVDQU32     256(R15), Z24
	VMOVDQU32     256(DX), Z14
	FROMMONTGOMERY(Z24, Z14, Z4, Z5, Z6, Z7, Z8, Z9, Z10, Z11, Z0, Z1)
	VEXTRACTI64X4 $1, Z24, Y25
	VPMOVZXWD     Y24, Z24
	VPMOVZXWD     Y25, Z25
	VMOVDQU32     512(CX), Z12
	MULD(Z24, Z12, Z2, Z3, Z4, Z5, Z6, Z7, Z0, Z1)
	VMOVDQU32     576(CX), Z13
	MULD(Z25, Z13, Z8, Z9, Z10, Z11, Z6, Z7, Z0, Z1)
	VEXTRACTI64X4 $1, Z14, Y15
	VPMOVZXWD     Y14, Z14
	VPMOVZXWD     Y15, Z15
	VMOVDQU32     512(BX), Z12
	MULD(Z14, Z12, Z2, Z3, Z4, Z5, Z6, Z7, Z0, Z1)
	VMOVDQU32     576(BX), Z13
	MULD(Z15, Z13, Z8, Z9, Z10, Z11, Z6, Z7, Z0, Z1)
	BUTTERFLYD2Q(Z24, Z14, Z0, Z4, Z6)
	BUTTERFLYD2Q(Z25, Z15, Z0, Z5, Z7)
	VMOVDQU32     512(DI), Z12
	MULD(Z14, Z12, Z2, Z3, Z4, Z5, Z6, Z7, Z0, Z1)
	VMOVDQU32     576(DI), Z13
	MULD(Z15, Z13, Z8, Z9, Z10, Z11, Z6, Z7, Z0, Z1)
	VMOVDQU32     Z14, 512(SP)
	VMOVDQU32     Z15, 576(SP)
	VMOVDQU32     320(R15), Z26
	VMOVDQU32     320(DX), Z14
	FROMMONTGOMERY(Z26, Z14, Z4, Z5, Z6, Z7, Z8, Z9, Z10, Z11, Z0, Z1)
	VEXTRACTI64X4 $1, Z26, Y27
	VPMOVZXWD     Y26, Z26
	VPMOVZXWD     Y27, Z27
	VMOVDQU32     640(CX), Z12
	MULD(Z26, Z12, Z2, Z3, Z4, Z5, Z6, Z7, Z0, Z1)
	VMOVDQU32     704(CX), Z13
	MULD(Z27, Z13, Z8, Z9, Z10, Z11, Z6, Z7, Z0, Z1)
	VEXTRACTI64X4 $1, Z14, Y15
	VPMOVZXWD     Y14, Z14
	VPMOVZXWD     Y15, Z15
	VMOVDQU32     640(BX), Z12
	MULD(Z14, Z12, Z2, Z3, Z4, Z5, Z6, Z7, Z0, Z1)
	VMOVDQU32     704(BX), Z13
	MULD(Z15, Z13, Z8, Z9, Z10, Z11, Z6, Z7, Z0, Z1)
	BUTTERFLYD2Q(Z26, Z14, Z0, Z4, Z6)
	BUTTERFLYD2Q(Z27, Z15, Z0, Z5, Z7)
	VMOVDQU32     640(DI), Z12
	MULD(Z14, Z12, Z2, Z3, Z4, Z5, Z6, Z7, Z0, Z1)
	VMOVDQU32     704(DI), Z13
	MULD(Z15, Z13, Z8, Z9, Z10, Z11, Z6, Z7, Z0, Z1)
	VMOVDQU32     Z14, 640(SP)
	VMOVDQU32     Z15, 704(SP)
	VMOVDQU32     384(R15), Z28
	VMOVDQU32     384(DX), Z14
	FROMMONTGOMERY(Z28, Z14, Z4, Z5, Z6, Z7, Z8, Z9, Z10, Z11, Z0, Z1)
	VEXTRACTI64X4 $1, Z28, Y29
	VPMOVZXWD     Y28, Z28
	VPMOVZXWD     Y29, Z29
	VMOVDQU32     768(CX), Z12
	MULD(Z28, Z12, Z2, Z3, Z4, Z5, Z6, Z7, Z0, Z1)
	VMOVDQU32     832(CX), Z13
	MULD(Z29, Z13, Z8, Z9, Z10, Z11, Z6, Z7, Z0, Z1)
	VEXTRACTI64X4 $1, Z14, Y15
	VPMOVZXWD     Y14, Z14
	VPMOVZXWD     Y15, Z15
	VMOVDQU32     768(BX), Z12
	MULD(Z14, Z12, Z2, Z3, Z4, Z5, Z6, Z7, Z0, Z1)
	VMOVDQU32     832(BX), Z13
	MULD(Z15, Z13, Z8, Z9, Z10, Z11, Z6, Z7, Z0, Z1)
	BUTTERFLYD2Q(Z28, Z14, Z0, Z4, Z6)
	BUTTERFLYD2Q(Z29, Z15, Z0, Z5, Z7)
	VMOVDQU32     768(DI), Z12
	MULD(Z14, Z12, Z2, Z3, Z4, Z5, Z6, Z7, Z0, Z1)
	VMOVDQU32     832(DI), Z13
	MULD(Z15, Z13, Z8, Z9, Z10, Z11, Z6, Z7, Z0, Z1)
	VMOVDQU32     Z14, 768(SP)
	VMOVDQU32     Z15, 832(SP)
	VMOVDQU32     448(R15), Z30
	VMOVDQU32     448(DX), Z14
	FROMMONTGOMERY(Z30, Z14, Z4, Z5, Z6, Z7, Z8, Z9, Z10, Z11, Z0, Z1)
	VEXTRACTI64X4 $1, Z30, Y31
	VPMOVZXWD     Y30, Z30
	VPMOVZXWD     Y31, Z31
	VMOVDQU32     896(CX), Z12
	MULD(Z30, Z12, Z2, Z3, Z4, Z5, Z6, Z7, Z0, Z1)
	VMOVDQU32     960(CX), Z13
	MULD(Z31, Z13, Z8, Z9, Z10, Z11, Z6, Z7, Z0, Z1)
	VEXTRACTI64X4 $1, Z14, Y15
	VPMOVZXWD     Y14, Z14
	VPMOVZXWD     Y15, Z15
	VMOVDQU32     896(BX), Z12
	MULD(Z14, Z12, Z2, Z3, Z4, Z5, Z6, Z7, Z0, Z1)
	VMOVDQU32     960(BX), Z13
	MULD(Z15, Z13, Z8, Z9, Z10, Z11, Z6, Z7, Z0, Z1)
	BUTTERFLYD2Q(Z30, Z14, Z0, Z4, Z6)
	BUTTERFLYD2Q(Z31, Z15, Z0, Z5, Z7)
	VMOVDQU32     896(DI), Z12
	MULD(Z14, Z12, Z2, Z3, Z4, Z5, Z6, Z7, Z0, Z1)
	VMOVDQU32     960(DI), Z13
	MULD(Z15, Z13, Z8, Z9, Z10, Z11, Z6, Z7, Z0, Z1)
	VMOVDQU32     Z14, 896(SP)
	VMOVDQU32     Z15, 960(SP)
	ADDQ          $24, SI
	MOVQ          SI, DI
	MOVQ          $2, R10

fft256_2:
	MOVQ         0(SI), R11
	VMOVDQU32    0(R11), Z2
	VMOVDQU32    64(R11), Z3
	VMOVDQU32    128(R11), Z4
	VMOVDQU32    192(R11), Z5
	VMOVDQU32    256(R11), Z6
	VMOVDQU32    320(R11), Z7
	VMOVDQU32    384(R11), Z8
	VMOVDQU32    448(R11), Z9
	BUTTERFLYD2Q(Z16, Z24, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z17, Z25, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z18, Z26, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z19, Z27, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z20, Z28, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z21, Z29, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z22, Z30, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z23, Z31, Z0, Z14, Z11)
	MULD(Z24, Z2, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z25, Z3, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z26, Z4, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z27, Z5, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z28, Z6, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z29, Z7, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z30, Z8, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z31, Z9, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	ADDQ         $24, SI
	MOVQ         0(SI), R11
	VMOVDQU32    0(R11), Z2
	VMOVDQU32    64(R11), Z3
	VMOVDQU32    128(R11), Z4
	VMOVDQU32    192(R11), Z5
	BUTTERFLYD2Q(Z16, Z20, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z17, Z21, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z18, Z22, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z19, Z23, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z24, Z28, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z25, Z29, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z26, Z30, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z27, Z31, Z0, Z14, Z11)
	MULD(Z20, Z2, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z21, Z3, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z22, Z4, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z23, Z5, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z28, Z2, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z29, Z3, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z30, Z4, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z31, Z5, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	ADDQ         $24, SI
	MOVQ         0(SI), R11
	VMOVDQU32    0(R11), Z2
	VMOVDQU32    64(R11), Z3
	BUTTERFLYD2Q(Z16, Z18, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z17, Z19, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z20, Z22, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z21, Z23, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z24, Z26, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z25, Z27, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z28, Z30, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z29, Z31, Z0, Z14, Z11)
	MULD(Z18, Z2, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z19, Z3, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z22, Z2, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z23, Z3, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z26, Z2, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z27, Z3, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z30, Z2, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z31, Z3, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	ADDQ         $24, SI
	MOVQ         0(SI), R11
	VMOVDQU32    0(R11), Z2
	BUTTERFLYD2Q(Z16, Z17, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z18, Z19, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z20, Z21, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z22, Z23, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z24, Z25, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z26, Z27, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z28, Z29, Z0, Z14, Z11)
	BUTTERFLYD2Q(Z30, Z31, Z0, Z14, Z11)
	MULD(Z17, Z2, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z19, Z2, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z21, Z2, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z23, Z2, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z25, Z2, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z27, Z2, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z29, Z2, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z31, Z2, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	ADDQ         $24, SI
	MOVQ         ·vInterleaveIndices+0(SB), R12
	VMOVDQU64    0(R12), Z6
	MOVQ         0(SI), R11
	VMOVDQU32    0(R11), Y2
	VINSERTI64X4 $1, Y2, Z2, Z2
	MOVQ         24(SI), R11
	VMOVDQU32    0(R11), X3
	VINSERTI64X2 $1, X3, Z3, Z3
	VINSERTI64X2 $0x0000000000000002, X3, Z3, Z3
	VINSERTI64X2 $0x0000000000000003, X3, Z3, Z3
	MOVQ         48(SI), R11
	VPBROADCASTD 0(R11), Z4
	VPBROADCASTD 4(R11), Z5
	VPBLENDMD    Z4, Z5, K3, Z4
	PERMUTE8X8(Z16, Z17, Z10)
	BUTTERFLYD2Q(Z16, Z17, Z0, Z14, Z11)
	PERMUTE8X8(Z18, Z19, Z10)
	BUTTERFLYD2Q(Z18, Z19, Z0, Z14, Z11)
	PERMUTE8X8(Z20, Z21, Z10)
	BUTTERFLYD2Q(Z20, Z21, Z0, Z14, Z11)
	PERMUTE8X8(Z22, Z23, Z10)
	BUTTERFLYD2Q(Z22, Z23, Z0, Z14, Z11)
	PERMUTE8X8(Z24, Z25, Z10)
	BUTTERFLYD2Q(Z24, Z25, Z0, Z14, Z11)
	PERMUTE8X8(Z26, Z27, Z10)
	BUTTERFLYD2Q(Z26, Z27, Z0, Z14, Z11)
	PERMUTE8X8(Z28, Z29, Z10)
	BUTTERFLYD2Q(Z28, Z29, Z0, Z14, Z11)
	PERMUTE8X8(Z30, Z31, Z10)
	BUTTERFLYD2Q(Z30, Z31, Z0, Z14, Z11)
	MULD(Z17, Z2, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z19, Z2, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z21, Z2, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z23, Z2, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z25, Z2, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z27, Z2, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z29, Z2, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z31, Z2, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	PERMUTE4X4(Z16, Z17, Z6, Z10)
	BUTTERFLYD2Q(Z16, Z17, Z0, Z14, Z11)
	PERMUTE4X4(Z18, Z19, Z6, Z10)
	BUTTERFLYD2Q(Z18, Z19, Z0, Z14, Z11)
	PERMUTE4X4(Z20, Z21, Z6, Z10)
	BUTTERFLYD2Q(Z20, Z21, Z0, Z14, Z11)
	PERMUTE4X4(Z22, Z23, Z6, Z10)
	BUTTERFLYD2Q(Z22, Z23, Z0, Z14, Z11)
	PERMUTE4X4(Z24, Z25, Z6, Z10)
	BUTTERFLYD2Q(Z24, Z25, Z0, Z14, Z11)
	PERMUTE4X4(Z26, Z27, Z6, Z10)
	BUTTERFLYD2Q(Z26, Z27, Z0, Z14, Z11)
	PERMUTE4X4(Z28, Z29, Z6, Z10)
	BUTTERFLYD2Q(Z28, Z29, Z0, Z14, Z11)
	PERMUTE4X4(Z30, Z31, Z6, Z10)
	BUTTERFLYD2Q(Z30, Z31, Z0, Z14, Z11)
	MULD(Z17, Z3, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z19, Z3, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z21, Z3, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z23, Z3, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z25, Z3, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z27, Z3, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z29, Z3, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	MULD(Z31, Z3, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	PERMUTE2X2(Z16, Z17, Z10)
	BUTTERFLYD2Q(Z16, Z17, Z0, Z14, Z11)
	PERMUTE2X2(Z18, Z19, Z10)
	BUTTERFLYD2Q(Z18, Z19, Z0, Z14, Z11)
	PERMUTE2X2(Z20, Z21, Z10)
	BUTTERFLYD2Q(Z20, Z21, Z0, Z14, Z11)
	PERMUTE2X2(Z22, Z23, Z10)
	BUTTERFLYD2Q(Z22, Z23, Z0, Z14, Z11)
	PERMUTE2X2(Z24, Z25, Z10)
	BUTTERFLYD2Q(Z24, Z25, Z0, Z14, Z11)
	PERMUTE2X2(Z26, Z27, Z10)
	BUTTERFLYD2Q(Z26, Z27, Z0, Z14, Z11)
	PERMUTE2X2(Z28, Z29, Z10)
	BUTTERFLYD2Q(Z28, Z29, Z0, Z14, Z11)
	PERMUTE2X2(Z30, Z31, Z10)
	BUTTERFLYD2Q(Z30, Z31, Z0, Z14, Z11)
	MULD(Z17, Z4, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	PERMUTE1X1(Z16, Z17, Z10)
	MULD(Z19, Z4, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	PERMUTE1X1(Z18, Z19, Z10)
	MULD(Z21, Z4, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	PERMUTE1X1(Z20, Z21, Z10)
	MULD(Z23, Z4, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	PERMUTE1X1(Z22, Z23, Z10)
	MULD(Z25, Z4, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	PERMUTE1X1(Z24, Z25, Z10)
	MULD(Z27, Z4, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	PERMUTE1X1(Z26, Z27, Z10)
	MULD(Z29, Z4, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	PERMUTE1X1(Z28, Z29, Z10)
	MULD(Z31, Z4, Z13, Z15, Z10, Z11, Z12, Z14, Z0, Z1)
	PERMUTE1X1(Z30, Z31, Z10)
	BUTTERFLYD2Q2Q(Z16, Z17, Z0, Z10)
	BUTTERFLYD2Q2Q(Z18, Z19, Z0, Z10)
	BUTTERFLYD2Q2Q(Z20, Z21, Z0, Z10)
	BUTTERFLYD2Q2Q(Z22, Z23, Z0, Z10)
	BUTTERFLYD2Q2Q(Z24, Z25, Z0, Z10)
	BUTTERFLYD2Q2Q(Z26, Z27, Z0, Z10)
	BUTTERFLYD2Q2Q(Z28, Z29, Z0, Z10)
	BUTTERFLYD2Q2Q(Z30, Z31, Z0, Z10)
	VMOVDQU32    0(R8), Z7
	MULD(Z16, Z7, Z9, Z2, Z3, Z4, Z5, Z10, Z0, Z1)
	VMOVDQU32    64(R8), Z8
	MULD(Z17, Z8, Z11, Z13, Z15, Z12, Z5, Z10, Z0, Z1)
	VMOVDQU32    128(R8), Z7
	MULD(Z18, Z7, Z9, Z2, Z3, Z4, Z5, Z10, Z0, Z1)
	VMOVDQU32    192(R8), Z8
	MULD(Z19, Z8, Z11, Z13, Z15, Z12, Z5, Z10, Z0, Z1)
	VMOVDQU32    256(R8), Z7
	MULD(Z20, Z7, Z9, Z2, Z3, Z4, Z5, Z10, Z0, Z1)
	VMOVDQU32    320(R8), Z8
	MULD(Z21, Z8, Z11, Z13, Z15, Z12, Z5, Z10, Z0, Z1)
	VMOVDQU32    384(R8), Z7
	MULD(Z22, Z7, Z9, Z2, Z3, Z4, Z5, Z10, Z0, Z1)
	VMOVDQU32    448(R8), Z8
	MULD(Z23, Z8, Z11, Z13, Z15, Z12, Z5, Z10, Z0, Z1)
	VMOVDQU32    512(R8), Z7
	MULD(Z24, Z7, Z9, Z2, Z3, Z4, Z5, Z10, Z0, Z1)
	VMOVDQU32    576(R8), Z8
	MULD(Z25, Z8, Z11, Z13, Z15, Z12, Z5, Z10, Z0, Z1)
	VMOVDQU32    640(R8), Z7
	MULD(Z26, Z7, Z9, Z2, Z3, Z4, Z5, Z10, Z0, Z1)
	VMOVDQU32    704(R8), Z8
	MULD(Z27, Z8, Z11, Z13, Z15, Z12, Z5, Z10, Z0, Z1)
	VMOVDQU32    768(R8), Z7
	MULD(Z28, Z7, Z9, Z2, Z3, Z4, Z5, Z10, Z0, Z1)
	VMOVDQU32    832(R8), Z8
	MULD(Z29, Z8, Z11, Z13, Z15, Z12, Z5, Z10, Z0, Z1)
	VMOVDQU32    896(R8), Z7
	MULD(Z30, Z7, Z9, Z2, Z3, Z4, Z5, Z10, Z0, Z1)
	VMOVDQU32    960(R8), Z8
	MULD(Z31, Z8, Z11, Z13, Z15, Z12, Z5, Z10, Z0, Z1)
	VPADDD       0(R9), Z16, Z16
	VPADDD       64(R9), Z17, Z17
	VPSUBD       Z0, Z16, Z5
	VPSUBD       Z0, Z17, Z10
	VPMINUD      Z5, Z16, Z16
	VPMINUD      Z10, Z17, Z17
	VMOVDQU32    Z16, 0(R9)
	VMOVDQU32    Z17, 64(R9)
	VPADDD       128(R9), Z18, Z18
	VPADDD       192(R9), Z19, Z19
	VPSUBD       Z0, Z18, Z5
	VPSUBD       Z0, Z19, Z10
	VPMINUD      Z5, Z18, Z18
	VPMINUD      Z10, Z19, Z19
	VMOVDQU32    Z18, 128(R9)
	VMOVDQU32    Z19, 192(R9)
	VPADDD       256(R9), Z20, Z20
	VPADDD       320(R9), Z21, Z21
	VPSUBD       Z0, Z20, Z5
	VPSUBD       Z0, Z21, Z10
	VPMINUD      Z5, Z20, Z20
	VPMINUD      Z10, Z21, Z21
	VMOVDQU32    Z20, 256(R9)
	VMOVDQU32    Z21, 320(R9)
	VPADDD       384(R9), Z22, Z22
	VPADDD       448(R9), Z23, Z23
	VPSUBD       Z0, Z22, Z5
	VPSUBD       Z0, Z23, Z10
	VPMINUD      Z5, Z22, Z22
	VPMINUD      Z10, Z23, Z23
	VMOVDQU32    Z22, 384(R9)
	VMOVDQU32    Z23, 448(R9)
	VPADDD       512(R9), Z24, Z24
	VPADDD       576(R9), Z25, Z25
	VPSUBD       Z0, Z24, Z5
	VPSUBD       Z0, Z25, Z10
	VPMINUD      Z5, Z24, Z24
	VPMINUD      Z10, Z25, Z25
	VMOVDQU32    Z24, 512(R9)
	VMOVDQU32    Z25, 576(R9)
	VPADDD       640(R9), Z26, Z26
	VPADDD       704(R9), Z27, Z27
	VPSUBD       Z0, Z26, Z5
	VPSUBD       Z0, Z27, Z10
	VPMINUD      Z5, Z26, Z26
	VPMINUD      Z10, Z27, Z27
	VMOVDQU32    Z26, 640(R9)
	VMOVDQU32    Z27, 704(R9)
	VPADDD       768(R9), Z28, Z28
	VPADDD       832(R9), Z29, Z29
	VPSUBD       Z0, Z28, Z5
	VPSUBD       Z0, Z29, Z10
	VPMINUD      Z5, Z28, Z28
	VPMINUD      Z10, Z29, Z29
	VMOVDQU32    Z28, 768(R9)
	VMOVDQU32    Z29, 832(R9)
	VPADDD       896(R9), Z30, Z30
	VPADDD       960(R9), Z31, Z31
	VPSUBD       Z0, Z30, Z5
	VPSUBD       Z0, Z31, Z10
	VPMINUD      Z5, Z30, Z30
	VPMINUD      Z10, Z31, Z31
	VMOVDQU32    Z30, 896(R9)
	VMOVDQU32    Z31, 960(R9)
	DECQ         R10
	TESTQ        R10, R10
	JEQ          done_1
	MOVQ         DI, SI
	VMOVDQU32    0(SP), Z16
	VMOVDQU32    64(SP), Z17
	VMOVDQU32    128(SP), Z18
	VMOVDQU32    192(SP), Z19
	VMOVDQU32    256(SP), Z20
	VMOVDQU32    320(SP), Z21
	VMOVDQU32    384(SP), Z22
	VMOVDQU32    448(SP), Z23
	VMOVDQU32    512(SP), Z24
	VMOVDQU32    576(SP), Z25
	VMOVDQU32    640(SP), Z26
	VMOVDQU32    704(SP), Z27
	VMOVDQU32    768(SP), Z28
	VMOVDQU32    832(SP), Z29
	VMOVDQU32    896(SP), Z30
	VMOVDQU32    960(SP), Z31
	ADDQ         $0x0000000000000400, R8
	ADDQ         $0x0000000000000400, R9
	JMP          fft256_2

done_1:
	RET

TEXT ·SisShuffle_avx512(SB), NOSPLIT, $0-24
	MOVQ      a+0(FP), R15
	MOVQ      a_len+8(FP), DX
	SHRQ      $5, DX
	LOAD_MASKS()
	MOVQ      ·vInterleaveIndices+0(SB), CX
	VMOVDQU64 0(CX), Z3

loop_3:
	TESTQ     DX, DX
	JEQ       done_4
	DECQ      DX
	VMOVDQU32 0(R15), Z1  // load a[i]
	VMOVDQU32 64(R15), Z2 // load a[i+16]
	PERMUTE8X8(Z1, Z2, Z0)
	PERMUTE4X4(Z1, Z2, Z3, Z0)
	PERMUTE2X2(Z1, Z2, Z0)
	PERMUTE1X1(Z1, Z2, Z0)
	VMOVDQU32 Z1, 0(R15)  // store a[i]
	VMOVDQU32 Z2, 64(R15) // store a[i+16]
	ADDQ      $128, R15
	JMP       loop_3

done_4:
	RET

TEXT ·SisUnshuffle_avx512(SB), NOSPLIT, $0-24
	MOVQ      a+0(FP), R15
	MOVQ      a_len+8(FP), DX
	SHRQ      $5, DX
	LOAD_MASKS()
	MOVQ      ·vInterleaveIndices+0(SB), CX
	VMOVDQU64 0(CX), Z3

loop_5:
	TESTQ      DX, DX
	JEQ        done_6
	DECQ       DX
	VMOVDQU32  0(R15), Z1  // load a[i]
	VMOVDQU32  64(R15), Z2 // load a[i+16]
	VPUNPCKLDQ Z2, Z1, Z0
	VPUNPCKHDQ Z2, Z1, Z2
	VMOVDQA32  Z0, Z1
	PERMUTE4X4(Z1, Z2, Z3, Z0)
	PERMUTE8X8(Z1, Z2, Z0)
	VMOVDQU32  Z1, 0(R15)  // store a[i]
	VMOVDQU32  Z2, 64(R15) // store a[i+16]
	ADDQ       $128, R15
	JMP        loop_5

done_6:
	RET
