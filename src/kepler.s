.include "header.s"
.text
.p2align 4
# High accuracy: (Gs1, Gs2, Gs3)
mm_stiefel_Gs13_avx512:
	vmulpd	%zmm1, %zmm1, %zmm2     # X^2
	vbroadcastsd	.LC1(%rip), %zmm3
	vbroadcastsd	.LC3(%rip), %zmm5
	vmulpd	%zmm2, %zmm0, %zmm0
	vfnmadd213pd	.LC5(%rip){1to8}, %zmm0, %zmm3
	vfnmadd213pd	.LC7(%rip){1to8}, %zmm0, %zmm5
	vfnmadd213pd	.LC9(%rip){1to8}, %zmm0, %zmm3
	vfnmadd213pd	.LC11(%rip){1to8}, %zmm0, %zmm5
	vfnmadd213pd	.LC13(%rip){1to8}, %zmm0, %zmm3
	vfnmadd213pd	.LC15(%rip){1to8}, %zmm0, %zmm5
	vfnmadd213pd	.LC17(%rip){1to8}, %zmm0, %zmm3
	vfnmadd213pd	.LC19(%rip){1to8}, %zmm0, %zmm5
	vfnmadd213pd	.LC21(%rip){1to8}, %zmm0, %zmm3
	vfnmadd213pd	.LC23(%rip){1to8}, %zmm0, %zmm5
	vfnmadd213pd	.LC25(%rip){1to8}, %zmm0, %zmm3
	vfnmadd213pd	.LC27(%rip){1to8}, %zmm0, %zmm5
	vfnmadd213pd	.LC29(%rip){1to8}, %zmm0, %zmm3
	vfnmadd213pd	.LC31(%rip){1to8}, %zmm0, %zmm5
	vfnmadd213pd	.LC33(%rip){1to8}, %zmm0, %zmm3
	vfnmadd213pd	.LC35(%rip){1to8}, %zmm0, %zmm5
	vmulpd	%zmm3, %zmm1, %zmm3
	vfnmadd132pd	%zmm3, %zmm1, %zmm0
	vmovapd	%zmm0, (%rdi)
	vmulpd	%zmm3, %zmm2, %zmm0
	vmovapd	%zmm0, (%rdx)
	vmulpd	%zmm5, %zmm2, %zmm2
	vmovapd	%zmm2, (%rsi)
	ret

# Low accuracy: (Gs0, Gs1, Gs2, Gs3)
.p2align 4
mm_stiefel_Gs03_avx512:
	vmulpd	%zmm1, %zmm1, %zmm2
	vbroadcastsd	.LC17(%rip), %zmm3
	vbroadcastsd	.LC19(%rip), %zmm4
	vmulpd	%zmm2, %zmm0, %zmm0
	vfnmadd213pd	.LC21(%rip){1to8}, %zmm0, %zmm5
	vfnmadd213pd	.LC23(%rip){1to8}, %zmm0, %zmm4
	vfnmadd213pd	.LC25(%rip){1to8}, %zmm0, %zmm5
	vfnmadd213pd	.LC27(%rip){1to8}, %zmm0, %zmm4
	vfnmadd213pd	.LC29(%rip){1to8}, %zmm0, %zmm5
	vfnmadd213pd	.LC31(%rip){1to8}, %zmm0, %zmm4
	vfnmadd213pd	.LC33(%rip){1to8}, %zmm0, %zmm5
	vmovapd	%zmm4, %zmm3
	vfnmadd213pd	.LC35(%rip){1to8}, %zmm0, %zmm4
	vfnmadd213pd	.LC37(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm3, (%rdi)
	vmulpd	%zmm5, %zmm1, %zmm3
	vfnmadd132pd	%zmm3, %zmm1, %zmm0
	vmovapd	%zmm0, (%rsi)
	vmulpd	%zmm3, %zmm2, %zmm0
	vmovapd	%zmm0, (%rcx)
	vmulpd	%zmm4, %zmm2, %zmm2
	vmovapd	%zmm2, (%rdx)
	ret

.p2align 4
.globl reb_whfast512_kepler_step
reb_whfast512_kepler_step:
	leaq	8(%rsp), %r10
	andq	$-64, %rsp
	pushq	-8(%r10)
	movq	%rdi, %rax
	pushq	%rbp
	movq	%rsp, %rbp
	pushq	%r10
	subq	$296, %rsp
.set X, %zmm13
.set Y, %zmm11
.set Z, %zmm15
.set VX, %zmm14
.set VY, %zmm12
.set VZ, %zmm10
	vmovapd	P512_X(%rdi), X
	vmovapd	P512_Y(%rdi), Y
	vmovapd	P512_Z(%rdi), Z
	vmulpd	X, X, %zmm6
	vfmadd231pd	Y, Y, %zmm6
	vfmadd231pd	Z, Z, %zmm6                 # r^2
	vsqrtpd	%zmm6, %zmm6                    # r
	vbroadcastsd	.LC37(%rip), %zmm9      # 1.0  #TODO: Keep in memory or try loading 8 doubles in one go
	vdivpd	%zmm6, %zmm9, %zmm16            # 1/r
	vmovapd	P512_VX(%rdi), VX
	vmovapd	P512_VY(%rdi), VY
	vmovapd	P512_VZ(%rdi), VZ
	vmulpd	VX, VX, %zmm0
	vfmadd231pd	VY, VY, %zmm0
	vfmadd231pd	VZ, VZ, %zmm0               # v^2
	vmovapd	P512_M(%rdi), %zmm4             # M
	vaddpd	%zmm4, %zmm4, %zmm17            # beta
	vmovapd	P512_DT(%rdi), %zmm8            # dt
	vfmsub132pd	%zmm16, %zmm0, %zmm17
	vmulpd	VX, X, %zmm0
	leaq	-240(%rbp), %rdx #             Gs2
	leaq	-304(%rbp), %rsi #             Gs1
	leaq	-176(%rbp), %rcx #             Gs3
	leaq	-112(%rbp), %rdi #             Gs0
	vmovapd	%zmm17, %zmm7
	vfmadd231pd	VY, Y, %zmm0
	vfnmadd132pd	%zmm6, %zmm4, %zmm7
	vmovapd	%zmm0, %zmm5
	vfmadd231pd	VZ, Z, %zmm5
	vmulpd	%zmm16, %zmm8, %zmm0
	vmovapd	%zmm5, %zmm18
	vmulpd	%zmm5, %zmm0, %zmm1
	vmovapd	%zmm5, %zmm21
	vmulpd	.LC35(%rip){1to8}, %zmm1, %zmm1
	vfnmadd132pd	%zmm16, %zmm9, %zmm1
	vmulpd	%zmm0, %zmm1, %zmm1
	vmovapd	%zmm17, %zmm0
            subq     $64, %rsp
            vmovdqu64 %zmm4, (%rsp)
            subq     $64, %rsp
            vmovdqu64 %zmm5, (%rsp)
	call	mm_stiefel_Gs03_avx512
            vmovdqu64 (%rsp), %zmm5
            addq     $64, %rsp
            vmovdqu64 (%rsp), %zmm4
            addq     $64, %rsp
	vmovapd	-304(%rbp), %zmm0
	vmovapd	%zmm6, %zmm3
	vfmadd132pd	%zmm0, %zmm6, %zmm18
	vfmsub132pd	%zmm1, %zmm8, %zmm3
	vmovapd	-240(%rbp), %zmm2
	vbroadcastsd	.LC41(%rip), %zmm19
	vfmadd231pd	%zmm2, %zmm5, %zmm3
	vfmadd132pd	%zmm7, %zmm18, %zmm2
	vmulpd	-112(%rbp), %zmm5, %zmm18
	vfmadd231pd	-176(%rbp), %zmm7, %zmm3
	vmulpd	%zmm2, %zmm2, %zmm20
	vfmadd132pd	%zmm7, %zmm18, %zmm0
	vbroadcastsd	.LC39(%rip), %zmm18
	vmulpd	%zmm19, %zmm20, %zmm20
	vmulpd	%zmm3, %zmm0, %zmm0
	vfnmadd132pd	%zmm18, %zmm20, %zmm0
	vbroadcastsd	.LC43(%rip), %zmm20
	vmulpd	%zmm20, %zmm3, %zmm3
	vsqrtpd	%zmm0, %zmm0
	vaddpd	%zmm0, %zmm2, %zmm2
	vmovapd	%zmm17, %zmm0
	vfmsub132pd	%zmm2, %zmm3, %zmm1
	vdivpd	%zmm2, %zmm1, %zmm1
            subq     $64, %rsp
            vmovdqu64 %zmm4, (%rsp)
            subq     $64, %rsp
            vmovdqu64 %zmm5, (%rsp)
	call	mm_stiefel_Gs03_avx512
            vmovdqu64 (%rsp), %zmm5
            addq     $64, %rsp
            vmovdqu64 (%rsp), %zmm4
            addq     $64, %rsp
	vmovapd	-304(%rbp), %zmm0
	vmovapd	%zmm6, %zmm3
	vfmadd132pd	%zmm0, %zmm6, %zmm21
	vfmsub132pd	%zmm1, %zmm8, %zmm3
	vmovapd	-240(%rbp), %zmm2
	movq	%rsi, %r8
	movq	%rdx, %rsi
	movq	%rcx, %rdx
	movq	%r8, %rdi
	vfmadd231pd	%zmm2, %zmm5, %zmm3
	vfmadd132pd	%zmm7, %zmm21, %zmm2
	vmulpd	-112(%rbp), %zmm5, %zmm21
	vfmadd231pd	-176(%rbp), %zmm7, %zmm3
	vfmadd132pd	%zmm7, %zmm21, %zmm0
	vmulpd	%zmm2, %zmm2, %zmm21
	vmulpd	%zmm0, %zmm3, %zmm0
	vmulpd	%zmm19, %zmm21, %zmm19
	vmulpd	%zmm20, %zmm3, %zmm3
	vfnmadd132pd	%zmm18, %zmm19, %zmm0
	vsqrtpd	%zmm0, %zmm0
	vaddpd	%zmm0, %zmm2, %zmm2
	vmovapd	%zmm17, %zmm0
	vfmsub132pd	%zmm2, %zmm3, %zmm1
	vdivpd	%zmm2, %zmm1, %zmm1
            subq     $64, %rsp
            vmovdqu64 %zmm4, (%rsp)
            subq     $64, %rsp
            vmovdqu64 %zmm5, (%rsp)
	call	mm_stiefel_Gs13_avx512
            vmovdqu64 (%rsp), %zmm5
            addq     $64, %rsp
            vmovdqu64 (%rsp), %zmm4
            addq     $64, %rsp
	vmulpd	-304(%rbp), %zmm5, %zmm2
	vmovapd	-240(%rbp), %zmm0
	vfmadd231pd	%zmm0, %zmm7, %zmm2
	vmulpd	%zmm2, %zmm1, %zmm1
	vfnmadd132pd	%zmm5, %zmm1, %zmm0
	vaddpd	%zmm6, %zmm2, %zmm1
	vdivpd	%zmm1, %zmm9, %zmm1
	vfnmadd231pd	-176(%rbp), %zmm7, %zmm0
	vaddpd	%zmm0, %zmm8, %zmm0
	vmulpd	%zmm0, %zmm1, %zmm1
	vmovapd	%zmm17, %zmm0
            subq     $64, %rsp
            vmovdqu64 %zmm4, (%rsp)
            subq     $64, %rsp
            vmovdqu64 %zmm5, (%rsp)
	call	mm_stiefel_Gs13_avx512
            vmovdqu64 (%rsp), %zmm5
            addq     $64, %rsp
            vmovdqu64 (%rsp), %zmm4
            addq     $64, %rsp
	vmovapd	-304(%rbp), %zmm17
	vmovapd	-240(%rbp), %zmm1
	vmulpd	%zmm5, %zmm17, %zmm0
	vmulpd	%zmm1, %zmm4, %zmm5
	kmovb	P512_MASK(%rax), %k1
	vfmadd231pd	%zmm1, %zmm7, %zmm0
	vmulpd	%zmm16, %zmm5, %zmm3
	vmovapd	%zmm8, %zmm1
	vfnmadd231pd	-176(%rbp), %zmm4, %zmm1
	vaddpd	%zmm6, %zmm0, %zmm0
	vdivpd	%zmm0, %zmm9, %zmm2
	vmulpd	%zmm17, %zmm4, %zmm0
	vmovapd	%zmm3, %zmm4
	vfnmadd132pd	Y, Y, %zmm4{%k1}{z}
	vmulpd	%zmm16, %zmm0, %zmm0
	vfmadd231pd	VY, %zmm1, %zmm4{%k1}{z}
	vmulpd	%zmm2, %zmm0, %zmm0
	vmulpd	%zmm5, %zmm2, %zmm2
	vmovapd	%zmm3, %zmm5
	vfnmadd132pd	X, X, %zmm5{%k1}{z}
	vfnmadd132pd	Z, Z, %zmm3{%k1}{z}
	vfnmadd132pd	%zmm2, VY, %zmm12{%k1}{z}
	vfmadd231pd	VX, %zmm1, %zmm5{%k1}{z}
	vfmadd132pd	VZ, %zmm3, %zmm1{%k1}{z}
	vfnmadd132pd	%zmm2, VX, %zmm14{%k1}{z}
	vfnmadd132pd	%zmm2, VZ, %zmm10{%k1}{z}
	vfnmadd132pd	%zmm0, %zmm12, %zmm11{%k1}{z}
	vfnmadd132pd	%zmm0, %zmm14, %zmm13{%k1}{z}
	vfnmadd132pd	%zmm15, %zmm10, %zmm0{%k1}{z}
	vmovapd	%zmm5, P512_X(%rax)
	vmovapd	%zmm13, P512_VX(%rax)
	vmovapd	%zmm4, P512_Y(%rax)
	vmovapd	%zmm11, P512_VY(%rax)
	vmovapd	%zmm1, P512_Z(%rax)
	vmovapd	%zmm0, P512_VZ(%rax)
	vzeroupper
	movq	-8(%rbp), %r10
	leave
	leaq	-8(%r10), %rsp
	ret

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
	.align 8
.LC7:
	.long	-416040929
	.long	1026222067
	.align 8
.LC9:
	.long	-416040929
	.long	1030416371
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
.LC37:      # 1.0
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
