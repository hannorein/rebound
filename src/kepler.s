.include "header.s"

.text

# general purpose:
# zmm0
# zmm1
# zmm2
# zmm3
# zmm4
# zmm5
# reserved/persistant:
.set R, %zmm9
.set RI, %zmm13
.set ZETA, %zmm14
.set ETA, %zmm15
# zmm 16
.set XX, %zmm17
.set GS0, %zmm20
.set GS1, %zmm21
.set GS2, %zmm22
.set GS3, %zmm23
.set BETA, %zmm24
# zmm 27
# zmm 29
# zmm 30
# zmm 31

# High accuracy: (Gs1, Gs2, Gs3)
mm_stiefel_Gs13_avx512:
    vmulpd    XX, XX, %zmm2     # X^2
    vbroadcastsd    .IF19(%rip), %zmm3
    vbroadcastsd    .IF18(%rip), %zmm4
    vmulpd    %zmm2, BETA, %zmm0
    vfnmadd213pd    .IF17(%rip){1to8}, %zmm0, %zmm3
    vfnmadd213pd    .IF16(%rip){1to8}, %zmm0, %zmm4
    vfnmadd213pd    .IF15(%rip){1to8}, %zmm0, %zmm3
    vfnmadd213pd    .IF14(%rip){1to8}, %zmm0, %zmm4
    vfnmadd213pd    .IF13(%rip){1to8}, %zmm0, %zmm3
    vfnmadd213pd    .IF12(%rip){1to8}, %zmm0, %zmm4
    vfnmadd213pd    .IF11(%rip){1to8}, %zmm0, %zmm3
    vfnmadd213pd    .IF10(%rip){1to8}, %zmm0, %zmm4
    vfnmadd213pd    .IF9(%rip){1to8}, %zmm0, %zmm3
    vfnmadd213pd    .IF8(%rip){1to8}, %zmm0, %zmm4
    vfnmadd213pd    .IF7(%rip){1to8}, %zmm0, %zmm3
    vfnmadd213pd    .IF6(%rip){1to8}, %zmm0, %zmm4
    vfnmadd213pd    .IF5(%rip){1to8}, %zmm0, %zmm3
    vfnmadd213pd    .IF4(%rip){1to8}, %zmm0, %zmm4
    vfnmadd213pd    .IF3(%rip){1to8}, %zmm0, %zmm3
    vfnmadd213pd    .IF2(%rip){1to8}, %zmm0, %zmm4
    vmulpd    %zmm3, XX, %zmm3
    vfnmadd132pd    %zmm3, XX, %zmm0
    vmovapd    %zmm0, GS1  # TODO: combine with previous instruction
    vmulpd    %zmm3, %zmm2, GS3
    vmulpd    %zmm4, %zmm2, GS2
    ret

# Low accuracy: (Gs0, Gs1, Gs2, Gs3)
mm_stiefel_Gs03_avx512:
    vmulpd    XX, XX, %zmm2
    vbroadcastsd    .IF11(%rip), %zmm3
    vbroadcastsd    .IF10(%rip), %zmm4
    vmulpd    %zmm2, BETA, %zmm0
    vfnmadd213pd    .IF9(%rip){1to8}, %zmm0, %zmm3
    vfnmadd213pd    .IF8(%rip){1to8}, %zmm0, %zmm4
    vfnmadd213pd    .IF7(%rip){1to8}, %zmm0, %zmm3
    vfnmadd213pd    .IF6(%rip){1to8}, %zmm0, %zmm4
    vfnmadd213pd    .IF5(%rip){1to8}, %zmm0, %zmm3
    vfnmadd213pd    .IF4(%rip){1to8}, %zmm0, %zmm4
    vfnmadd213pd    .IF3(%rip){1to8}, %zmm0, %zmm3
    vmovapd    %zmm4, GS0
    vfnmadd213pd    .IF2(%rip){1to8}, %zmm0, %zmm4
    vfnmadd213pd    .IF1(%rip){1to8}, %zmm0, GS0
    vmulpd    %zmm3, XX, %zmm3
    vfnmadd132pd    %zmm3, XX, %zmm0
    vmovapd    %zmm0, GS1  # TODO: combine with previous instruction
    vmulpd    %zmm3, %zmm2, GS3
    vmulpd    %zmm4, %zmm2, GS2
    ret

halley:
    # In: GS0,GS1,GS2,GS3
    # Out: XX
    # No other registers used.
    vfmsub213pd     DT, ZETA, GS3
    vfmadd231pd     GS2, ETA, GS3
    vfmadd231pd     XX, R, GS3              # f

    vfmadd132pd     ZETA, R, GS2
    vfmadd231pd     ETA, GS1, GS2           # fp

    vmulpd          GS0, ETA, GS0
    vfmadd132pd     ZETA, GS0, GS1          # fpp

    vmulpd          GS1, GS3, GS1           # f*fpp
    # Next instruction uses exponent trick to speed up multiplication by 0.5 using integer sub
    vpsubq          HALF_MASK, GS1, GS1     # 0.5*f*fpp
    vfmsub231pd     GS2, GS2, GS1           # fp*fp-0.5*f*fpp
    vmulpd          GS3, GS2, GS3           # f*fp
    vdivpd          GS1, GS3, GS3
    vsubpd          GS3, XX, XX
    ret


.p2align 4
#.globl reb_whfast512_kepler_step #not used right now
.extern reb_whfast512_init_registers
.globl reb_whfast512_kepler_step_noinit
reb_whfast512_kepler_step:
    # Need to init registers here when not called after interaction step.
    # This will be required for synchronizing.  
    call reb_whfast512_init_registers

reb_whfast512_kepler_step_noinit:

    vmulpd    X, X, %zmm0
    vfmadd231pd    Y, Y, %zmm0
    vfmadd231pd    Z, Z, %zmm0                 # r^2
    vsqrtpd    %zmm0, R                        # r
    vdivpd    R, ONE, RI                # 1/r
    vmulpd    VX, VX, %zmm0
    vfmadd231pd    VY, VY, %zmm0
    vfmadd231pd    VZ, VZ, %zmm0               # v^2
    vaddpd    M, M, BETA                      # 2*M
    vfmsub132pd    RI, %zmm0, BETA         # beta
    vmulpd    VX, X, ETA
    vfmadd231pd    VY, Y, ETA
    vfmadd231pd    VZ, Z, ETA                  # eta
    vmovapd    BETA, ZETA
    vfnmadd132pd    R, M, ZETA              # zeta
    vmulpd    RI, DT, %zmm5              # dt/r
    vmulpd    ETA, %zmm5, %zmm4             # eta*dt/r
    vmulpd    .DOUBLE_HALF(%rip){1to8}, %zmm4, %zmm4     # 0.5*eta*dt/r
    vfnmadd132pd    RI, ONE, %zmm4        
    vmulpd    %zmm5, %zmm4, XX              # X (initial guess)
    
    call    mm_stiefel_Gs03_avx512
    call    halley

    call    mm_stiefel_Gs03_avx512
    call    halley
    

    call    mm_stiefel_Gs13_avx512
    
    # Newton
    vmulpd    GS1, ETA, %zmm2
    vfmadd231pd    GS2, ZETA, %zmm2
    vmulpd    %zmm2, XX, XX
    vmovapd    GS2, %zmm0
    vfnmadd132pd    ETA, XX, %zmm0
    vaddpd    R, %zmm2, XX
    vdivpd    XX, ONE, XX
    vfnmadd231pd    GS3, ZETA, %zmm0
    vaddpd    %zmm0, DT, %zmm0
    vmulpd    %zmm0, XX, XX
    
    call    mm_stiefel_Gs13_avx512
    
    # Newton (XX no longer needed)
    vmulpd    GS1, ETA, %zmm2
    vfmadd231pd    GS2, ZETA, %zmm2
    vaddpd    R, %zmm2, XX
    vdivpd  XX, ONE, %zmm4       # ri in C
    
    vmulpd    GS2, M, %zmm0
    vmulpd    RI, %zmm0, %zmm3        # negative f
    vmulpd    %zmm0, %zmm4, %zmm2     # negative gd
    vmovapd    DT, %zmm1
    vfnmadd231pd    GS3, M, %zmm1   # g 
    vmulpd    GS1, M, %zmm0
    vmulpd    RI, %zmm0, %zmm0
    vmulpd    %zmm4, %zmm0, %zmm0     # negative fd


    vmovapd    %zmm3, %zmm4
    vmovapd    %zmm3, %zmm5
    // Calculate new vx vy vz
    vfnmadd132pd    X, X, %zmm3{%k1}{z}
    vfnmadd132pd    Y, Y, %zmm4{%k1}{z}
    vfnmadd132pd    Z, Z, %zmm5{%k1}{z}
    vfmadd231pd        VX, %zmm1, %zmm3{%k1}{z}
    vfmadd231pd        VY, %zmm1, %zmm4{%k1}{z}
    vfmadd231pd        VZ, %zmm1, %zmm5{%k1}{z}
    // Calculate new xyz
    vfnmadd132pd    %zmm2, VX, VX{%k1}{z}
    vfnmadd132pd    %zmm2, VY, VY{%k1}{z}
    vfnmadd132pd    %zmm2, VZ, VZ{%k1}{z}
    vfnmadd231pd    %zmm0, X, VX{%k1}{z}
    vfnmadd231pd    %zmm0, Y, VY{%k1}{z}
    vfnmadd231pd    %zmm0, Z, VZ{%k1}{z}
    vmovapd    %zmm3, P512_X(%rdi) # TODO: shuffle things around so that new X end up in X without copy
    vmovapd    %zmm4, P512_Y(%rdi)
    vmovapd    %zmm5, P512_Z(%rdi)
    vmovapd    %zmm3, X # TODO: shuffle things around so that new X end up in X without copy
    vmovapd    %zmm4, Y
    vmovapd    %zmm5, Z
    ret

.section    .rodata
invfactorial:
.IF0:
    .double     1.0
.IF1:
    .double     1.0
.IF2:
.DOUBLE_HALF:
    .double     0.5
.IF3:
    .long    1431655765
    .long    1069897045
.IF4:
    .long    1431655765
    .long    1067799893
.IF5:
    .long    286331153
    .long    1065423121
.IF6:
    .long    381774871
    .long    1062650220
.IF7:
    .long    436314138
    .long    1059717536
.IF8:
    .long    436314138
    .long    1056571808
.IF9:
    .long    -1521039564
    .long    1053236707
.IF10:
    .long    -1216831652
    .long    1049787983
.IF11:
    .long    1744127204
    .long    1046144581
.IF12:
    .long    -268904296
    .long    1042411224
.IF13:
    .long    329805065
    .long    1038488134
.IF14:
    .long    -1463780195
    .long    1034500468
.IF15:
    .long    -416040929
    .long    1030416371
.IF16:
    .long    -416040929
    .long    1026222067
.IF17:
    .long    1882238282
    .long    1021924039
.IF18:
    .long    1673100695
    .long    1017545336
.IF19:
    .long    1182875991
    .long    1013118107
# Parts of inverse vactorial not used at the moment:    
#    .long    -1543372251
#    .long    1008620587
#    .long    -153291406
#    .long    1003953038
#    .long    1840773203
#    .long    999345721
#    .long    320223257
#    .long    994533812
#    .long    426964344
#    .long    989801712
#    .long    -585736279
#    .long    984871884
#    .long    -60141991
#    .long    979930757
#    .long    -1025716573
#    .long    974985905
#    .long    641009757
#    .long    969974154
#    .long    -1958520658
#    .long    964844025
#    .long    1346885134
#    .long    959681324
#    .long    -410782275
#    .long    954479826
#    .long    -410782275
#    .long    949236946
#    .long    1423773010
#    .long    943953938
#    .long    834731386
#    .long    938635522
    #


.section    .note.GNU-stack,"",@progbits
    
#.align 8
#.DOUBLE_FIVE:
#    .long    0
#    .long    1075052544
#.align 8
#.DOUBLE_SIXTEEN:
#    .long    0
#    .long    1076887552
#.align 8
#.DOUBLE_TWENTY:
#    .long    0
#    .long    1077149696
#.align 8
#
#.set FIVE, %zmm17
#.set SIXTEEN, %zmm18
#.set TWENTY, %zmm19
#    vbroadcastsd    .DOUBLE_FIVE(%rip), FIVE
#    vbroadcastsd    .DOUBLE_SIXTEEN(%rip), SIXTEEN
#    vbroadcastsd    .DOUBLE_TWENTY(%rip), TWENTY # constants for Halley

    # Old Halley
    # vmulpd    GS2, GS2, %zmm4
    # vmulpd    SIXTEEN, %zmm4, %zmm4
    # vmulpd    %zmm3, GS1, %zmm0
    # vfnmadd132pd    TWENTY, %zmm4, %zmm0
    # vmulpd    FIVE, %zmm3, %zmm3
    # vsqrtpd    %zmm0, %zmm0
    # vaddpd    %zmm0, %zmm2, %zmm2
    # vfmsub132pd    %zmm2, %zmm3, XX
    # vdivpd    %zmm2, XX, XX
