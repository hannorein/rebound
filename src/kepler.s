.text
.p2align 4
.globl	reb_whfast512_kepler_step
reb_whfast512_kepler_step:
	vmovapd	384(%rdi), %zmm14
	vmovapd	448(%rdi), %zmm12
	vmulpd	%zmm14, %zmm14, %zmm5
	vmovapd	576(%rdi), %zmm15
	vmovapd	640(%rdi), %zmm13
	vmulpd	%zmm15, %zmm14, %zmm4
	vmulpd	%zmm15, %zmm15, %zmm0
	vfmadd231pd	%zmm12, %zmm12, %zmm5
	vmovapd	512(%rdi), %zmm6
	vmovapd	704(%rdi), %zmm11
	vfmadd231pd	%zmm13, %zmm12, %zmm4
	vfmadd231pd	%zmm13, %zmm13, %zmm0
	vfmadd231pd	%zmm6, %zmm6, %zmm5
	vbroadcastsd	.LC1(%rip), %zmm6
	vmovapd	(%rdi), %zmm3
	vfmadd231pd	512(%rdi), %zmm11, %zmm4
	vfmadd231pd	%zmm11, %zmm11, %zmm0
	vsqrtpd	%zmm5, %zmm5
	vdivpd	%zmm5, %zmm6, %zmm16
	vmulpd	64(%rdi), %zmm16, %zmm2
	vaddpd	%zmm3, %zmm3, %zmm1
	vbroadcastsd	.LC3(%rip), %zmm17
	vbroadcastsd	.LC5(%rip), %zmm26
	vfmsub132pd	%zmm16, %zmm0, %zmm1
	vmulpd	%zmm4, %zmm2, %zmm0
	vbroadcastsd	.LC7(%rip), %zmm24
	vbroadcastsd	.LC13(%rip), %zmm22
	vbroadcastsd	.LC9(%rip), %zmm25
	vbroadcastsd	.LC11(%rip), %zmm23
	vmulpd	%zmm17, %zmm0, %zmm0
	vbroadcastsd	.LC17(%rip), %zmm20
	vbroadcastsd	.LC15(%rip), %zmm21
	vbroadcastsd	.LC21(%rip), %zmm18
	vbroadcastsd	.LC19(%rip), %zmm19
	vfnmadd132pd	%zmm16, %zmm6, %zmm0
	vmovapd	%zmm1, %zmm7
	vmovapd	%zmm5, %zmm28
	vfnmadd132pd	%zmm5, %zmm3, %zmm7
	vbroadcastsd	.LC25(%rip), %zmm31
	vmulpd	%zmm2, %zmm0, %zmm0
	vbroadcastsd	.LC23(%rip), %zmm30
	vbroadcastsd	.LC27(%rip), %zmm29
	vmulpd	%zmm0, %zmm0, %zmm9
	vfmsub213pd	64(%rdi), %zmm0, %zmm28
	vmulpd	%zmm9, %zmm1, %zmm2
	vmovapd	%zmm2, %zmm27
	vfnmadd132pd	%zmm26, %zmm24, %zmm27
	vmovapd	%zmm2, %zmm8
	vfnmadd132pd	%zmm25, %zmm23, %zmm8
	vmovapd	%zmm2, %zmm10
	vfnmadd132pd	%zmm2, %zmm22, %zmm27
	vfnmadd132pd	%zmm2, %zmm21, %zmm8
	vfnmadd132pd	%zmm2, %zmm20, %zmm27
	vfnmadd132pd	%zmm2, %zmm19, %zmm8
	vfnmadd132pd	%zmm2, %zmm18, %zmm27
	vfnmadd132pd	%zmm2, %zmm17, %zmm8
	vmulpd	%zmm27, %zmm0, %zmm27
	vfnmadd132pd	%zmm8, %zmm6, %zmm10
	vmulpd	%zmm8, %zmm9, %zmm8
	vfnmadd132pd	%zmm27, %zmm0, %zmm2
	vmulpd	%zmm27, %zmm9, %zmm9
	vmovapd	%zmm4, %zmm27
	vfmadd231pd	%zmm8, %zmm4, %zmm28
	vmulpd	%zmm10, %zmm4, %zmm10
	vfmadd132pd	%zmm2, %zmm5, %zmm27
	vfmadd132pd	%zmm7, %zmm28, %zmm9
	vfmadd132pd	%zmm7, %zmm10, %zmm2
	vfmadd132pd	%zmm7, %zmm27, %zmm8
	vmovapd	%zmm5, %zmm28
	vmulpd	%zmm9, %zmm2, %zmm2
	vmulpd	%zmm8, %zmm8, %zmm10
	vmulpd	%zmm29, %zmm9, %zmm9
	vmulpd	%zmm31, %zmm10, %zmm10
	vfnmadd132pd	%zmm30, %zmm10, %zmm2
	vsqrtpd	%zmm2, %zmm2
	vaddpd	%zmm2, %zmm8, %zmm8
	vfmsub132pd	%zmm8, %zmm9, %zmm0
	vdivpd	%zmm8, %zmm0, %zmm0
	vmulpd	%zmm0, %zmm0, %zmm9
	vfmsub213pd	64(%rdi), %zmm0, %zmm28
	vmulpd	%zmm9, %zmm1, %zmm2
	vmovapd	%zmm2, %zmm10
	vfnmadd132pd	%zmm26, %zmm24, %zmm10
	vmovapd	%zmm2, %zmm8
	vfnmadd132pd	%zmm25, %zmm23, %zmm8
	vmovapd	%zmm2, %zmm27
	vfnmadd132pd	%zmm2, %zmm22, %zmm10
	vfnmadd132pd	%zmm2, %zmm21, %zmm8
	vfnmadd132pd	%zmm2, %zmm20, %zmm10
	vfnmadd132pd	%zmm2, %zmm19, %zmm8
	vfnmadd132pd	%zmm2, %zmm18, %zmm10
	vfnmadd132pd	%zmm2, %zmm17, %zmm8
	vmulpd	%zmm10, %zmm0, %zmm10
	vfnmadd132pd	%zmm8, %zmm6, %zmm27
	vmulpd	%zmm8, %zmm9, %zmm8
	vfnmadd132pd	%zmm10, %zmm0, %zmm2
	vmulpd	%zmm10, %zmm9, %zmm9
	vmovapd	%zmm4, %zmm10
	vfmadd231pd	%zmm8, %zmm4, %zmm28
	vmulpd	%zmm27, %zmm4, %zmm27
	vfmadd132pd	%zmm2, %zmm5, %zmm10
	vfmadd132pd	%zmm7, %zmm28, %zmm9
	vfmadd132pd	%zmm7, %zmm27, %zmm2
	vfmadd231pd	%zmm8, %zmm7, %zmm10
	vmulpd	%zmm2, %zmm9, %zmm2
	vmulpd	%zmm10, %zmm10, %zmm8
	vmulpd	%zmm29, %zmm9, %zmm9
	vbroadcastsd	.LC41(%rip), %zmm29
	vmulpd	%zmm31, %zmm8, %zmm8
	vbroadcastsd	.LC31(%rip), %zmm31
	vfnmadd132pd	%zmm30, %zmm8, %zmm2
	vbroadcastsd	.LC29(%rip), %zmm8
	vbroadcastsd	.LC39(%rip), %zmm30
	vsqrtpd	%zmm2, %zmm2
	vaddpd	%zmm2, %zmm10, %zmm10
	vfmsub132pd	%zmm10, %zmm9, %zmm0
	vbroadcastsd	.LC33(%rip), %zmm9
	vdivpd	%zmm10, %zmm0, %zmm0
	vmulpd	%zmm0, %zmm0, %zmm27
	vmulpd	%zmm27, %zmm1, %zmm2
	vmovapd	%zmm2, %zmm10
	vfnmadd132pd	%zmm8, %zmm31, %zmm10
	vfnmadd213pd	.LC34(%rip), %zmm2, %zmm9
	vfnmadd213pd	.LC36(%rip), %zmm2, %zmm10
	vfnmadd132pd	%zmm2, %zmm30, %zmm9
	vbroadcastsd	.LC43(%rip), %zmm28
	kmovb	2368(%rdi), %k1
	vfnmadd132pd	%zmm2, %zmm29, %zmm10
	vfnmadd132pd	%zmm2, %zmm28, %zmm9
	vfnmadd132pd	%zmm2, %zmm26, %zmm10
	vfnmadd132pd	%zmm2, %zmm25, %zmm9
	vfnmadd132pd	%zmm2, %zmm24, %zmm10
	vfnmadd132pd	%zmm2, %zmm23, %zmm9
	vfnmadd132pd	%zmm2, %zmm22, %zmm10
	vfnmadd132pd	%zmm2, %zmm21, %zmm9
	vfnmadd132pd	%zmm2, %zmm20, %zmm10
	vfnmadd132pd	%zmm2, %zmm19, %zmm9
	vfnmadd132pd	%zmm2, %zmm18, %zmm10
	vfnmadd132pd	%zmm2, %zmm17, %zmm9
	vmulpd	%zmm10, %zmm0, %zmm10
	vmulpd	%zmm9, %zmm27, %zmm9
	vfnmadd132pd	%zmm10, %zmm0, %zmm2
	vmulpd	%zmm10, %zmm27, %zmm27
	vmulpd	%zmm2, %zmm4, %zmm2
	vfmadd231pd	%zmm9, %zmm7, %zmm2
	vmulpd	%zmm2, %zmm0, %zmm0
	vaddpd	%zmm5, %zmm2, %zmm2
	vdivpd	%zmm2, %zmm6, %zmm2
	vfnmadd132pd	%zmm4, %zmm0, %zmm9
	vfnmadd132pd	%zmm7, %zmm9, %zmm27
	vbroadcastsd	.LC33(%rip), %zmm9
	vaddpd	64(%rdi), %zmm27, %zmm27
	vmulpd	%zmm27, %zmm2, %zmm2
	vmulpd	%zmm2, %zmm2, %zmm0
	vmulpd	%zmm0, %zmm1, %zmm1
	vfnmadd132pd	%zmm1, %zmm31, %zmm8
	vfnmadd213pd	.LC34(%rip), %zmm1, %zmm9
	vfnmadd213pd	.LC36(%rip), %zmm1, %zmm8
	vmovapd	%zmm8, %zmm10
	vfnmadd132pd	%zmm1, %zmm29, %zmm10
	vmovapd	%zmm9, %zmm8
	vfnmadd132pd	%zmm1, %zmm30, %zmm8
	vmovapd	%zmm10, %zmm9
	vfnmadd132pd	%zmm1, %zmm26, %zmm9
	vfnmadd132pd	%zmm1, %zmm28, %zmm8
	vfnmadd132pd	%zmm1, %zmm24, %zmm9
	vfnmadd132pd	%zmm1, %zmm25, %zmm8
	vfnmadd132pd	%zmm1, %zmm22, %zmm9
	vfnmadd132pd	%zmm1, %zmm23, %zmm8
	vfnmadd132pd	%zmm1, %zmm20, %zmm9
	vfnmadd132pd	%zmm1, %zmm21, %zmm8
	vfnmadd132pd	%zmm1, %zmm18, %zmm9
	vfnmadd132pd	%zmm1, %zmm19, %zmm8
	vmulpd	%zmm9, %zmm2, %zmm9
	vfnmadd132pd	%zmm1, %zmm17, %zmm8
	vfnmadd132pd	%zmm9, %zmm2, %zmm1
	vmulpd	%zmm8, %zmm0, %zmm2
	vmulpd	%zmm9, %zmm0, %zmm0
	vmulpd	%zmm1, %zmm4, %zmm4
	vfnmadd213pd	64(%rdi), %zmm3, %zmm0
	vfmadd132pd	%zmm2, %zmm4, %zmm7
	vmulpd	%zmm2, %zmm3, %zmm2
	vmulpd	%zmm1, %zmm3, %zmm3
	vaddpd	%zmm5, %zmm7, %zmm7
	vmulpd	%zmm16, %zmm2, %zmm4
	vmulpd	%zmm16, %zmm3, %zmm3
	vmovapd	512(%rdi), %zmm5
	vdivpd	%zmm7, %zmm6, %zmm6
	vmovapd	%zmm4, %zmm1
	vmulpd	%zmm6, %zmm3, %zmm3
	vmulpd	%zmm2, %zmm6, %zmm6
	vmovapd	%zmm4, %zmm2
	vfnmadd132pd	%zmm14, %zmm14, %zmm2{%k1}{z}
	vfnmadd132pd	%zmm12, %zmm12, %zmm1{%k1}{z}
	vfnmadd132pd	%zmm5, %zmm5, %zmm4{%k1}{z}
	vfmadd231pd	%zmm15, %zmm0, %zmm2{%k1}{z}
	vfmadd231pd	%zmm13, %zmm0, %zmm1{%k1}{z}
	vfnmadd132pd	%zmm6, %zmm15, %zmm15{%k1}{z}
	vfmadd132pd	%zmm11, %zmm4, %zmm0{%k1}{z}
	vfnmadd132pd	%zmm6, %zmm13, %zmm13{%k1}{z}
	vfnmadd132pd	%zmm6, %zmm11, %zmm11{%k1}{z}
	vmovapd	%zmm2, 384(%rdi)
	vfnmadd132pd	%zmm3, %zmm15, %zmm14{%k1}{z}
	vmovapd	%zmm1, 448(%rdi)
	vfnmadd132pd	%zmm3, %zmm13, %zmm12{%k1}{z}
	vfnmadd132pd	%zmm5, %zmm11, %zmm3{%k1}{z}
	vmovapd	%zmm0, 512(%rdi)
	vmovapd	%zmm14, 576(%rdi)
	vmovapd	%zmm12, 640(%rdi)
	vmovapd	%zmm3, 704(%rdi)
	vzeroupper
	ret

.section	.rodata.cst8,"aM",@progbits,8
.align 8
.LC1:
	.long	0
	.long	1072693248
	.align 8
.LC3:
	.long	0
	.long	1071644672
	.align 8
.LC5:
	.long	1744127204
	.long	1046144581
	.align 8
.LC7:
	.long	-1521039564
	.long	1053236707
	.align 8
.LC9:
	.long	-1216831652
	.long	1049787983
	.align 8
.LC11:
	.long	436314138
	.long	1056571808
	.align 8
.LC13:
	.long	436314138
	.long	1059717536
	.align 8
.LC15:
	.long	381774871
	.long	1062650220
	.align 8
.LC17:
	.long	286331153
	.long	1065423121
	.align 8
.LC19:
	.long	1431655765
	.long	1067799893
	.align 8
.LC21:
	.long	1431655765
	.long	1069897045
	.align 8
.LC23:
	.long	0
	.long	1077149696
	.align 8
.LC25:
	.long	0
	.long	1076887552
	.align 8
.LC27:
	.long	0
	.long	1075052544
	.align 8
.LC29:
	.long	1182875991
	.long	1013118107
	.align 8
.LC31:
	.long	1882238282
	.long	1021924039
	.align 8
.LC33:
	.long	1673100695
	.long	1017545336

.section	.rodata
.align 64
.LC34:
	.long	-416040929
	.long	1026222067
	.long	-416040929
	.long	1026222067
	.long	-416040929
	.long	1026222067
	.long	-416040929
	.long	1026222067
	.long	-416040929
	.long	1026222067
	.long	-416040929
	.long	1026222067
	.long	-416040929
	.long	1026222067
	.long	-416040929
	.long	1026222067
	.align 64
.LC36:
	.long	-416040929
	.long	1030416371
	.long	-416040929
	.long	1030416371
	.long	-416040929
	.long	1030416371
	.long	-416040929
	.long	1030416371
	.long	-416040929
	.long	1030416371
	.long	-416040929
	.long	1030416371
	.long	-416040929
	.long	1030416371
	.long	-416040929
	.long	1030416371

.section	.rodata.cst8
.align 8
.LC39:
	.long	-1463780195
	.long	1034500468
	.align 8
.LC41:
	.long	329805065
	.long	1038488134
	.align 8
.LC43:
	.long	-268904296
	.long	1042411224
	
    .section	.note.GNU-stack,"",@progbits
