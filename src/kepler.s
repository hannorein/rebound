.include "header.s"
.text
.p2align 4

.set GS0, %zmm25
.set GS1, %zmm26
.set GS2, %zmm27
.set GS3, %zmm28
.set BETA, %zmm29
.set XX, %zmm1

# High accuracy: (Gs1, Gs2, Gs3)
mm_stiefel_Gs13_avx512:
	vmulpd	XX, XX, %zmm2     # X^2
	vbroadcastsd	.LC1(%rip), %zmm3
	vbroadcastsd	.LC3(%rip), %zmm4
	vmulpd	%zmm2, BETA, %zmm0
	vfnmadd213pd	.LC5(%rip){1to8}, %zmm0, %zmm3
	vfnmadd213pd	.LC7(%rip){1to8}, %zmm0, %zmm4
	vfnmadd213pd	.LC9(%rip){1to8}, %zmm0, %zmm3
	vfnmadd213pd	.LC11(%rip){1to8}, %zmm0, %zmm4
	vfnmadd213pd	.LC13(%rip){1to8}, %zmm0, %zmm3
	vfnmadd213pd	.LC15(%rip){1to8}, %zmm0, %zmm4
	vfnmadd213pd	.LC17(%rip){1to8}, %zmm0, %zmm3
	vfnmadd213pd	.LC19(%rip){1to8}, %zmm0, %zmm4
	vfnmadd213pd	.LC21(%rip){1to8}, %zmm0, %zmm3
	vfnmadd213pd	.LC23(%rip){1to8}, %zmm0, %zmm4
	vfnmadd213pd	.LC25(%rip){1to8}, %zmm0, %zmm3
	vfnmadd213pd	.LC27(%rip){1to8}, %zmm0, %zmm4
	vfnmadd213pd	.LC29(%rip){1to8}, %zmm0, %zmm3
	vfnmadd213pd	.LC31(%rip){1to8}, %zmm0, %zmm4
	vfnmadd213pd	.LC33(%rip){1to8}, %zmm0, %zmm3
	vfnmadd213pd	.LC35(%rip){1to8}, %zmm0, %zmm4
	vmulpd	%zmm3, XX, %zmm3
	vfnmadd132pd	%zmm3, XX, %zmm0
	vmovapd	%zmm0, GS1  # TODO: combine with previous instruction
	vmulpd	%zmm3, %zmm2, GS3
	vmulpd	%zmm4, %zmm2, GS2
	ret

# Low accuracy: (Gs0, Gs1, Gs2, Gs3)
.p2align 4
mm_stiefel_Gs03_avx512:
	vmulpd	XX, XX, %zmm2
	vbroadcastsd	.LC17(%rip), %zmm3
	vbroadcastsd	.LC19(%rip), %zmm4
	vmulpd	%zmm2, BETA, %zmm0
	vfnmadd213pd	.LC21(%rip){1to8}, %zmm0, %zmm3
	vfnmadd213pd	.LC23(%rip){1to8}, %zmm0, %zmm4
	vfnmadd213pd	.LC25(%rip){1to8}, %zmm0, %zmm3
	vfnmadd213pd	.LC27(%rip){1to8}, %zmm0, %zmm4
	vfnmadd213pd	.LC29(%rip){1to8}, %zmm0, %zmm3
	vfnmadd213pd	.LC31(%rip){1to8}, %zmm0, %zmm4
	vfnmadd213pd	.LC33(%rip){1to8}, %zmm0, %zmm3
	vmovapd	%zmm4, GS0
	vfnmadd213pd	.LC35(%rip){1to8}, %zmm0, %zmm4
	vfnmadd213pd	.LC37(%rip){1to8}, %zmm0, GS0
	vmulpd	%zmm3, XX, %zmm3
	vfnmadd132pd	%zmm3, XX, %zmm0
	vmovapd	%zmm0, GS1  # TODO: combine with previous instruction
	vmulpd	%zmm3, %zmm2, GS3
	vmulpd	%zmm4, %zmm2, GS2
	ret

.p2align 4
.globl reb_whfast512_kepler_step
reb_whfast512_kepler_step:
.set X, %zmm13
.set Y, %zmm11
.set Z, %zmm15
.set VX, %zmm14
.set VY, %zmm12
.set VZ, %zmm10
.set M, %zmm22
.set DT, %zmm23
.set R, %zmm6
.set RI, %zmm16
.set ETA, %zmm8
.set ZETA, %zmm7
.set ONE, %zmm9
	vbroadcastsd	.LC39(%rip), %zmm21 # constants for Halley
	vbroadcastsd	.LC41(%rip), %zmm19
	vbroadcastsd	.LC43(%rip), %zmm20
	vmovapd	P512_X(%rdi), X
	vmovapd	P512_Y(%rdi), Y
	vmovapd	P512_Z(%rdi), Z
	vmovapd	P512_VX(%rdi), VX
	vmovapd	P512_VY(%rdi), VY
	vmovapd	P512_VZ(%rdi), VZ
	vmovapd	P512_DT(%rdi), DT           
	vmovapd	P512_M(%rdi), M
	kmovb	P512_MASK(%rdi), %k1
	vmulpd	X, X, %zmm0
	vfmadd231pd	Y, Y, %zmm0
	vfmadd231pd	Z, Z, %zmm0                 # r^2
	vsqrtpd	%zmm0, R                        # r
	vbroadcastsd	.LC37(%rip), ONE      # 1.0  #TODO: Keep in memory or try loading 8 doubles in one go
	vdivpd	R, ONE, RI                # 1/r
	vmulpd	VX, VX, %zmm0
	vfmadd231pd	VY, VY, %zmm0
	vfmadd231pd	VZ, VZ, %zmm0               # v^2
	vaddpd	M, M, BETA                      # 2*M
	vfmsub132pd	RI, %zmm0, BETA         # beta
	vmulpd	VX, X, ETA
	vfmadd231pd	VY, Y, ETA
	vfmadd231pd	VZ, Z, ETA                  # eta
	vmovapd	BETA, ZETA
	vfnmadd132pd	R, M, ZETA              # zeta
	vmulpd	RI, DT, %zmm5              # dt/r
	vmulpd	ETA, %zmm5, %zmm4             # eta*dt/r
	vmulpd	.LC35(%rip){1to8}, %zmm4, %zmm4     # 0.5*eta*dt/r
	vfnmadd132pd	RI, ONE, %zmm4        
	vmulpd	%zmm5, %zmm4, XX              # X (initial guess)
	
	call	mm_stiefel_Gs03_avx512

    # Halley
	vmovapd	ETA, %zmm5
	vfmadd132pd	GS1, R, %zmm5
	vmovapd	R, %zmm3
	vfmsub132pd	XX, DT, %zmm3
	vfmadd231pd	GS2, ETA, %zmm3
	vfmadd132pd	ZETA, %zmm5, GS2
	vfmadd231pd	GS3, ZETA, %zmm3        # f
	vmulpd	GS0, ETA, %zmm5
	vfmadd132pd	ZETA, %zmm5, GS1        # fpp
	vmulpd	GS2, GS2, %zmm4
	vmulpd	%zmm19, %zmm4, %zmm4
	vmulpd	%zmm3, GS1, %zmm0
	vfnmadd132pd	%zmm21, %zmm4, %zmm0
	vmulpd	%zmm20, %zmm3, %zmm3
	vsqrtpd	%zmm0, %zmm0
	vaddpd	%zmm0, %zmm2, %zmm2
	vfmsub132pd	%zmm2, %zmm3, XX
	vdivpd	%zmm2, XX, XX
	
    call	mm_stiefel_Gs03_avx512
	
    # Halley
	vmovapd	R, %zmm3
	vfmadd132pd	GS1, R, %zmm5
	vfmsub132pd	XX, DT, %zmm3
	vmovapd	GS2, %zmm2
	vfmadd231pd	GS2, ETA, %zmm3
	vfmadd132pd	ZETA, %zmm5, %zmm2
	vmulpd	GS0, ETA, %zmm5
	vfmadd231pd	GS3, ZETA, %zmm3
	vfmadd132pd	ZETA, %zmm5, GS1
	vmulpd	%zmm2, %zmm2, %zmm5
	vmulpd	GS1, %zmm3, %zmm0
	vmulpd	%zmm19, %zmm5, %zmm4
	vmulpd	%zmm20, %zmm3, %zmm3
	vfnmadd132pd	%zmm21, %zmm4, %zmm0
	vsqrtpd	%zmm0, %zmm0
	vaddpd	%zmm0, %zmm2, %zmm2
	vfmsub132pd	%zmm2, %zmm3, XX
	vdivpd	%zmm2, XX, XX

	call	mm_stiefel_Gs13_avx512
	
    # Newton
    vmulpd	GS1, ETA, %zmm2
	vfmadd231pd	GS2, ZETA, %zmm2
	vmulpd	%zmm2, %zmm3, %zmm3
	vmovapd	GS2, %zmm0
	vfnmadd132pd	ETA, %zmm3, %zmm0
	vaddpd	R, %zmm2, %zmm3
	vdivpd	%zmm3, ONE, %zmm3
	vfnmadd231pd	GS3, ZETA, %zmm0
	vaddpd	%zmm0, DT, %zmm0
	vmulpd	%zmm0, %zmm3, XX
	
    call	mm_stiefel_Gs13_avx512
	
    # Newton
    vmulpd	ETA, GS1, %zmm0
	vfmadd231pd	GS2, ZETA, %zmm0
	vaddpd	R, %zmm0, %zmm0
	vdivpd	%zmm0, ONE, %zmm4       # ri in C
	
	vmulpd	GS2, M, %zmm0
	vmulpd	RI, %zmm0, %zmm3        # negative f
	vmulpd	%zmm0, %zmm4, %zmm2     # negative gd
	vmovapd	DT, %zmm1
	vfnmadd231pd	GS3, M, %zmm1   # g 
	vmulpd	GS1, M, %zmm0
	vmulpd	RI, %zmm0, %zmm0
	vmulpd	%zmm4, %zmm0, %zmm0     # negative fd


	vmovapd	%zmm3, %zmm4
	vmovapd	%zmm3, %zmm5
    // Calculate new vx vy vz
	vfnmadd132pd	X, X, %zmm3{%k1}{z}
	vfnmadd132pd	Y, Y, %zmm4{%k1}{z}
	vfnmadd132pd	Z, Z, %zmm5{%k1}{z}
	vfmadd231pd	    VX, %zmm1, %zmm3{%k1}{z}
	vfmadd231pd	    VY, %zmm1, %zmm4{%k1}{z}
	vfmadd231pd	    VZ, %zmm1, %zmm5{%k1}{z}
    // Calculate new xyz
	vfnmadd132pd	%zmm2, VX, VX{%k1}{z}
	vfnmadd132pd	%zmm2, VY, VY{%k1}{z}
	vfnmadd132pd	%zmm2, VZ, VZ{%k1}{z}
	vfnmadd231pd	%zmm0, X, VX{%k1}{z}
	vfnmadd231pd	%zmm0, Y, VY{%k1}{z}
	vfnmadd231pd	%zmm0, Z, VZ{%k1}{z}
	vmovapd	%zmm3, P512_X(%rdi) # TODO: shuffle things around so that new X end up in X without copy
	vmovapd	%zmm4, P512_Y(%rdi)
	vmovapd	%zmm5, P512_Z(%rdi)
	vmovapd	VX, P512_VX(%rdi)  # TODO: Remove memory
	vmovapd	VY, P512_VY(%rdi)
	vmovapd	VZ, P512_VZ(%rdi)
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
