	.file	"tmp.c"
	.text
	.p2align 4
	.globl	mm_stiefel_Gs13_avx512
	.type	mm_stiefel_Gs13_avx512, @function
mm_stiefel_Gs13_avx512:
.LFB6404:
	.cfi_startproc
	vmulpd	%zmm1, %zmm1, %zmm2
	vbroadcastsd	.LC1(%rip), %zmm3
	vmovapd	%zmm3, (%rdx)
	vbroadcastsd	.LC3(%rip), %zmm3
	vmulpd	%zmm2, %zmm0, %zmm0
	vmovapd	%zmm3, (%rsi)
	vmovapd	(%rdx), %zmm3
	vfnmadd213pd	.LC5(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rdx)
	vmovapd	(%rsi), %zmm3
	vfnmadd213pd	.LC7(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rsi)
	vmovapd	(%rdx), %zmm3
	vfnmadd213pd	.LC9(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rdx)
	vmovapd	(%rsi), %zmm3
	vfnmadd213pd	.LC11(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rsi)
	vmovapd	(%rdx), %zmm3
	vfnmadd213pd	.LC13(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rdx)
	vmovapd	(%rsi), %zmm3
	vfnmadd213pd	.LC15(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rsi)
	vmovapd	(%rdx), %zmm3
	vfnmadd213pd	.LC17(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rdx)
	vmovapd	(%rsi), %zmm3
	vfnmadd213pd	.LC19(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rsi)
	vmovapd	(%rdx), %zmm3
	vfnmadd213pd	.LC21(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rdx)
	vmovapd	(%rsi), %zmm3
	vfnmadd213pd	.LC23(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rsi)
	vmovapd	(%rdx), %zmm3
	vfnmadd213pd	.LC25(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rdx)
	vmovapd	(%rsi), %zmm3
	vfnmadd213pd	.LC27(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rsi)
	vmovapd	(%rdx), %zmm3
	vfnmadd213pd	.LC29(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rdx)
	vmovapd	(%rsi), %zmm3
	vfnmadd213pd	.LC31(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rsi)
	vmovapd	(%rdx), %zmm3
	vfnmadd213pd	.LC33(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rdx)
	vmovapd	(%rsi), %zmm3
	vfnmadd213pd	.LC35(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rsi)
	vmulpd	(%rdx), %zmm1, %zmm3
	vfnmadd132pd	%zmm3, %zmm1, %zmm0
	vmovapd	%zmm3, (%rdx)
	vmovapd	%zmm0, (%rdi)
	vmulpd	(%rdx), %zmm2, %zmm0
	vmovapd	%zmm0, (%rdx)
	vmulpd	(%rsi), %zmm2, %zmm2
	vmovapd	%zmm2, (%rsi)
	ret
	.cfi_endproc
.LFE6404:
	.size	mm_stiefel_Gs13_avx512, .-mm_stiefel_Gs13_avx512
	.p2align 4
	.globl	mm_stiefel_Gs03_avx512
	.type	mm_stiefel_Gs03_avx512, @function
mm_stiefel_Gs03_avx512:
.LFB6405:
	.cfi_startproc
	vmulpd	%zmm1, %zmm1, %zmm2
	vbroadcastsd	.LC17(%rip), %zmm3
	vmovapd	%zmm3, (%rcx)
	vbroadcastsd	.LC19(%rip), %zmm3
	vmulpd	%zmm2, %zmm0, %zmm0
	vmovapd	%zmm3, (%rdx)
	vmovapd	(%rcx), %zmm3
	vfnmadd213pd	.LC21(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rcx)
	vmovapd	(%rdx), %zmm3
	vfnmadd213pd	.LC23(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rdx)
	vmovapd	(%rcx), %zmm3
	vfnmadd213pd	.LC25(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rcx)
	vmovapd	(%rdx), %zmm3
	vfnmadd213pd	.LC27(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rdx)
	vmovapd	(%rcx), %zmm3
	vfnmadd213pd	.LC29(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rcx)
	vmovapd	(%rdx), %zmm3
	vfnmadd213pd	.LC31(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rdx)
	vmovapd	(%rcx), %zmm3
	vfnmadd213pd	.LC33(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rcx)
	vmovapd	(%rdx), %zmm3
	vfnmadd213pd	.LC35(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rdx)
	vfnmadd213pd	.LC37(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rdi)
	vmulpd	(%rcx), %zmm1, %zmm3
	vfnmadd132pd	%zmm3, %zmm1, %zmm0
	vmovapd	%zmm3, (%rcx)
	vmovapd	%zmm0, (%rsi)
	vmulpd	(%rcx), %zmm2, %zmm0
	vmovapd	%zmm0, (%rcx)
	vmulpd	(%rdx), %zmm2, %zmm2
	vmovapd	%zmm2, (%rdx)
	ret
	.cfi_endproc
.LFE6405:
	.size	mm_stiefel_Gs03_avx512, .-mm_stiefel_Gs03_avx512
	.p2align 4
	.globl	reb_whfast512_kepler_step
	.type	reb_whfast512_kepler_step, @function
reb_whfast512_kepler_step:
.LFB6406:
	.cfi_startproc
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
	vbroadcastsd	.LC37(%rip), %zmm6
	vmovapd	(%rdi), %zmm3
	vfmadd231pd	512(%rdi), %zmm11, %zmm4
	vfmadd231pd	%zmm11, %zmm11, %zmm0
	vsqrtpd	%zmm5, %zmm5
	vdivpd	%zmm5, %zmm6, %zmm16
	vmulpd	64(%rdi), %zmm16, %zmm2
	vaddpd	%zmm3, %zmm3, %zmm1
	vbroadcastsd	.LC35(%rip), %zmm17
	vbroadcastsd	.LC17(%rip), %zmm26
	vfmsub132pd	%zmm16, %zmm0, %zmm1
	vmulpd	%zmm4, %zmm2, %zmm0
	vbroadcastsd	.LC21(%rip), %zmm24
	vbroadcastsd	.LC25(%rip), %zmm22
	vbroadcastsd	.LC19(%rip), %zmm25
	vbroadcastsd	.LC23(%rip), %zmm23
	vmulpd	%zmm17, %zmm0, %zmm0
	vbroadcastsd	.LC29(%rip), %zmm20
	vbroadcastsd	.LC27(%rip), %zmm21
	vbroadcastsd	.LC33(%rip), %zmm18
	vbroadcastsd	.LC31(%rip), %zmm19
	vfnmadd132pd	%zmm16, %zmm6, %zmm0
	vmovapd	%zmm1, %zmm7
	vmovapd	%zmm5, %zmm28
	vfnmadd132pd	%zmm5, %zmm3, %zmm7
	vbroadcastsd	.LC41(%rip), %zmm31
	vmulpd	%zmm2, %zmm0, %zmm0
	vbroadcastsd	.LC39(%rip), %zmm30
	vbroadcastsd	.LC43(%rip), %zmm29
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
	vbroadcastsd	.LC13(%rip), %zmm29
	vmulpd	%zmm31, %zmm8, %zmm8
	vbroadcastsd	.LC5(%rip), %zmm31
	vfnmadd132pd	%zmm30, %zmm8, %zmm2
	vbroadcastsd	.LC1(%rip), %zmm8
	vbroadcastsd	.LC11(%rip), %zmm30
	vsqrtpd	%zmm2, %zmm2
	vaddpd	%zmm2, %zmm10, %zmm10
	vfmsub132pd	%zmm10, %zmm9, %zmm0
	vbroadcastsd	.LC3(%rip), %zmm9
	vdivpd	%zmm10, %zmm0, %zmm0
	vmulpd	%zmm0, %zmm0, %zmm27
	vmulpd	%zmm27, %zmm1, %zmm2
	vmovapd	%zmm2, %zmm10
	vfnmadd132pd	%zmm8, %zmm31, %zmm10
	vfnmadd213pd	.LC6(%rip), %zmm2, %zmm9
	vfnmadd213pd	.LC8(%rip), %zmm2, %zmm10
	vfnmadd132pd	%zmm2, %zmm30, %zmm9
	vbroadcastsd	.LC15(%rip), %zmm28
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
	vbroadcastsd	.LC3(%rip), %zmm9
	vaddpd	64(%rdi), %zmm27, %zmm27
	vmulpd	%zmm27, %zmm2, %zmm2
	vmulpd	%zmm2, %zmm2, %zmm0
	vmulpd	%zmm0, %zmm1, %zmm1
	vfnmadd132pd	%zmm1, %zmm31, %zmm8
	vfnmadd213pd	.LC6(%rip), %zmm1, %zmm9
	vfnmadd213pd	.LC8(%rip), %zmm1, %zmm8
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
	.cfi_endproc
.LFE6406:
	.size	reb_whfast512_kepler_step, .-reb_whfast512_kepler_step
	.globl	invfactorial
	.section	.rodata
	.align 32
	.type	invfactorial, @object
	.size	invfactorial, 280
invfactorial:
	.long	0
	.long	1072693248
	.long	0
	.long	1072693248
	.long	0
	.long	1071644672
	.long	1431655765
	.long	1069897045
	.long	1431655765
	.long	1067799893
	.long	286331153
	.long	1065423121
	.long	381774871
	.long	1062650220
	.long	436314138
	.long	1059717536
	.long	436314138
	.long	1056571808
	.long	-1521039564
	.long	1053236707
	.long	-1216831652
	.long	1049787983
	.long	1744127204
	.long	1046144581
	.long	-268904296
	.long	1042411224
	.long	329805065
	.long	1038488134
	.long	-1463780195
	.long	1034500468
	.long	-416040929
	.long	1030416371
	.long	-416040929
	.long	1026222067
	.long	1882238282
	.long	1021924039
	.long	1673100695
	.long	1017545336
	.long	1182875991
	.long	1013118107
	.long	-1543372251
	.long	1008620587
	.long	-153291406
	.long	1003953038
	.long	1840773203
	.long	999345721
	.long	320223257
	.long	994533812
	.long	426964344
	.long	989801712
	.long	-585736279
	.long	984871884
	.long	-60141991
	.long	979930757
	.long	-1025716573
	.long	974985905
	.long	641009757
	.long	969974154
	.long	-1958520658
	.long	964844025
	.long	1346885134
	.long	959681324
	.long	-410782275
	.long	954479826
	.long	-410782275
	.long	949236946
	.long	1423773010
	.long	943953938
	.long	834731386
	.long	938635522
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC1:
	.long	1182875991
	.long	1013118107
	.align 8
.LC3:
	.long	1673100695
	.long	1017545336
	.align 8
.LC5:
	.long	1882238282
	.long	1021924039
	.section	.rodata
	.align 64
.LC6:
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
	.set	.LC7,.LC6
	.align 64
.LC8:
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
	.set	.LC9,.LC8
	.section	.rodata.cst8
	.align 8
.LC11:
	.long	-1463780195
	.long	1034500468
	.align 8
.LC13:
	.long	329805065
	.long	1038488134
	.align 8
.LC15:
	.long	-268904296
	.long	1042411224
	.align 8
.LC17:
	.long	1744127204
	.long	1046144581
	.align 8
.LC19:
	.long	-1216831652
	.long	1049787983
	.align 8
.LC21:
	.long	-1521039564
	.long	1053236707
	.align 8
.LC23:
	.long	436314138
	.long	1056571808
	.align 8
.LC25:
	.long	436314138
	.long	1059717536
	.align 8
.LC27:
	.long	381774871
	.long	1062650220
	.align 8
.LC29:
	.long	286331153
	.long	1065423121
	.align 8
.LC31:
	.long	1431655765
	.long	1067799893
	.align 8
.LC33:
	.long	1431655765
	.long	1069897045
	.align 8
.LC35:
	.long	0
	.long	1071644672
	.align 8
.LC37:
	.long	0
	.long	1072693248
	.align 8
.LC39:
	.long	0
	.long	1077149696
	.align 8
.LC41:
	.long	0
	.long	1076887552
	.align 8
.LC43:
	.long	0
	.long	1075052544
	.ident	"GCC: (Debian 12.2.0-14+deb12u1) 12.2.0"
	.section	.note.GNU-stack,"",@progbits
