# file: integrator_whfast512.s
.section .text
.globl block1_gr
.globl block1_nogr
.globl reb_whfast512_kepler_step

#P512 Structure offsets
.set P512_M, 0
.set P512_DT, 64
.set P512_GR_PREFAC, 128
.set P512_GR_PREFAC2, 192
.set P512_m, 320
.set P512_X, 384
.set P512_Y, 448
.set P512_Z, 512
.set P512_VX, 576
.set P512_VY, 640
.set P512_VZ, 704
.set P512_MAT8_INERTIAL_TO_JACOBI, 768
.set P512_MAT8_JACOBI_TO_HELIOCENTRIC, 1280
.set P512_MAT8_INERTIAL_TO_JACOBI_T, 2304
.set P512_M0, 2816
.set P512_MASK, 2880

#####################################
# Register alias
#####################################

# Only used in Interaction step:
.set HVX, %zmm13
.set HVY, %zmm14
.set HVZ, %zmm15
.set HVXC, %zmm16
.set HVYC, %zmm17
.set HVZC, %zmm18
.set HX, %zmm19
.set HY, %zmm20
.set HZ, %zmm21

# Only used in Kepler step:
.set R, %zmm13
.set RI, %zmm14
.set ZETA, %zmm15
.set ETA, %zmm16
.set XX, %zmm17
.set GS0, %zmm18
.set GS1, %zmm0     # Note: reusing register zmm0
.set GS2, %zmm19
.set GS3, %zmm20
.set BETA, %zmm21

# Common register use
.set X, %zmm22
.set Y, %zmm23
.set Z, %zmm24
.set VX, %zmm25
.set VY, %zmm26
.set VZ, %zmm27
.set ONE, %zmm28
.set DT, %zmm29
.set HALF_MASK, %zmm30
.set M, %zmm31
.set M_DT, %zmm12   # Only used once per step.

.macro reb_whfast512_init_registers
    # Ignore exceptions
    subq $8, %rsp
    stmxcsr (%rsp)
    orl $0x8040, (%rsp)    # Set Global FTZ and DAZ
    ldmxcsr (%rsp)
    addq $8, %rsp

    # Load data
    kmovw           P512_MASK(%rdi), %k1
    vmovapd         P512_DT(%rdi), DT
    vmovapd         P512_M(%rdi), M
    vmulpd          DT, M, M_DT
    vbroadcastsd    .DOUBLE_ONE(%rip), ONE
    vpbroadcastq    .HALF_MASK(%rip), HALF_MASK
    vmovdqa64    .SIGN_MASK(%rip), %zmm10
    vmovapd      .EPS(%rip),     %zmm11
    
    vmovapd     P512_X(%rdi), X
    vmovapd     P512_Y(%rdi), Y
    vmovapd     P512_Z(%rdi), Z
    
    vmovapd     P512_VX(%rdi), VX
    vmovapd     P512_VY(%rdi), VY
    vmovapd     P512_VZ(%rdi), VZ
.endm  

.macro gravity_prefactor multiplier=ONE
    # Input:  zmm0=dx, zmm1=dy, zmm2=dy
    # Output: zmm6 = multiplier / r^3
    
    vmulpd      %zmm0, %zmm0, %zmm6     
    vfmadd231pd %zmm1, %zmm1, %zmm6      
    vfmadd231pd %zmm2, %zmm2, %zmm6     # zmm6 is now r^2
    
    vsqrtpd     %zmm6, %zmm7             
    vmulpd      %zmm6, %zmm7, %zmm6     # zmm6 is r^3
   
    vdivpd      %zmm6, \multiplier, %zmm6      
.endm


.macro REDUCE_ADD_AND_BROADCAST reg, temp_reg
    vshuff64x2 $0x4E, \reg, \reg, \temp_reg         # 01234567 -> 45670123
    vaddpd     \reg, \temp_reg, \temp_reg    
    vshufpd    $0x55, \temp_reg, \temp_reg, \reg    # 01234567 -> 10325476
    vaddpd     \reg, \temp_reg, \temp_reg
    vpermpd    $0x4E, \temp_reg, \reg               # 01234567 -> 23016745
    vaddpd     \reg, \temp_reg, \reg
.endm
        
.macro mat8_mul3 in0, in1, in2, out0, out1, out2
    # 8x8 matrix multiplied with 3 different 8 vectors
    # in: rax = vector to 64 matrix elements
    # zmm0, zmm1, zmm2  input and output vectors
    # uses: zmm3-zmm7
    # The idea is to use embedded broadcast loads
    # Note: matrix needs to be transposed.
    vmovupd \in0,   0(%rsp)
    vmovupd \in1,  64(%rsp)
    vmovupd \in2, 128(%rsp)

    # Keeping six independent FMA chains going
    vmovapd        (%rax), %zmm4
    vmovapd      64(%rax), %zmm3
    vmulpd        0(%rsp){1to8}, %zmm4, \out0
    vmulpd       64(%rsp){1to8}, %zmm4, \out1
    vmulpd      128(%rsp){1to8}, %zmm4, \out2

    vmulpd        8(%rsp){1to8}, %zmm3, %zmm5
    vmulpd       72(%rsp){1to8}, %zmm3, %zmm6
    vmulpd      136(%rsp){1to8}, %zmm3, %zmm7

    vmovapd     128(%rax), %zmm4
    vmovapd     192(%rax), %zmm3
    vfmadd231pd  16(%rsp){1to8}, %zmm4, \out0
    vfmadd231pd  80(%rsp){1to8}, %zmm4, \out1
    vfmadd231pd 144(%rsp){1to8}, %zmm4, \out2
    
    vfmadd231pd  24(%rsp){1to8}, %zmm3, %zmm5
    vfmadd231pd  88(%rsp){1to8}, %zmm3, %zmm6
    vfmadd231pd 152(%rsp){1to8}, %zmm3, %zmm7
    
    vmovapd     256(%rax), %zmm4
    vmovapd     320(%rax), %zmm3
    vfmadd231pd  32(%rsp){1to8}, %zmm4, \out0
    vfmadd231pd  96(%rsp){1to8}, %zmm4, \out1
    vfmadd231pd 160(%rsp){1to8}, %zmm4, \out2
    
    vfmadd231pd  40(%rsp){1to8}, %zmm3, %zmm5
    vfmadd231pd 104(%rsp){1to8}, %zmm3, %zmm6
    vfmadd231pd 168(%rsp){1to8}, %zmm3, %zmm7
    
    vmovapd     384(%rax), %zmm4
    vmovapd     448(%rax), %zmm3
    vfmadd231pd  48(%rsp){1to8}, %zmm4, \out0
    vfmadd231pd 112(%rsp){1to8}, %zmm4, \out1
    vfmadd231pd 176(%rsp){1to8}, %zmm4, \out2
    
    vfmadd231pd  56(%rsp){1to8}, %zmm3, %zmm5
    vfmadd231pd 120(%rsp){1to8}, %zmm3, %zmm6
    vfmadd231pd 184(%rsp){1to8}, %zmm3, %zmm7

    # Using two accumulators, adding at end
    vaddpd  %zmm5, \out0, \out0
    vaddpd  %zmm6, \out1, \out1
    vaddpd  %zmm7, \out2, \out2
   .endm

.macro vfnmadd_auto_inc reg1 reg2
    vfnmadd213pd    .IF0+(IF_offset*8)(%rip){1to8}, \reg1, \reg2
    .set IF_offset, IF_offset - 1
.endm

# High accuracy: (Gs1, Gs2, Gs3)
# Output: GS1==%zmm0, GS2, GS3
# numTerms must be an odd number
.macro mm_stiefel_Gs13_avx512 numTerms=19
    .set IF_offset, \numTerms
    vmulpd          XX, XX, %zmm2     # X^2
    vbroadcastsd    .IF0+(IF_offset*8)(%rip), %zmm3
    .set IF_offset, IF_offset - 1
    vbroadcastsd    .IF0+(IF_offset*8)(%rip), %zmm4
    .set IF_offset, IF_offset - 1
    vmulpd          %zmm2, BETA, %zmm0
    .set GS_iterations, (IF_offset -1)/2
    .rept GS_iterations
    vfnmadd_auto_inc %zmm0, %zmm3
    vfnmadd_auto_inc %zmm0, %zmm4
    .endr
    vmulpd          %zmm4, %zmm2, GS2
    vmulpd          %zmm3, XX, %zmm3
    vmulpd          %zmm3, %zmm2, GS3
    vfnmadd132pd    %zmm3, XX, %zmm0 # = GS1
.endm

# Low accuracy: (Gs0, Gs1, Gs2, Gs3)
# Output: GS0, GS1==%zmm0, GS2, GS3
# numTerms must be an odd number
.macro mm_stiefel_Gs03_avx512 numTerms=11
    .set IF_offset, \numTerms
    vmulpd          XX, XX, %zmm2     # X^2
    vbroadcastsd    .IF0+(IF_offset*8)(%rip), %zmm3
    .set IF_offset, IF_offset - 1
    vbroadcastsd    .IF0+(IF_offset*8)(%rip), %zmm4
    .set IF_offset, IF_offset - 1
    vmulpd          %zmm2, BETA, %zmm0
    .set GS_iterations, (IF_offset -1)/2 -1
    .rept GS_iterations
    vfnmadd_auto_inc %zmm0, %zmm3
    vfnmadd_auto_inc %zmm0, %zmm4
    .endr
    vfnmadd213pd    .IF0+(3*8)(%rip){1to8}, %zmm0, %zmm3
    vmovapd         %zmm4, GS0
    vfnmadd213pd    .IF0+(2*8)(%rip){1to8}, %zmm0, %zmm4
    vfnmadd213pd    .IF0+(1*8)(%rip){1to8}, %zmm0, GS0
    vmulpd          %zmm4, %zmm2, GS2
    vmulpd          %zmm3, XX, %zmm3
    vmulpd          %zmm3, %zmm2, GS3
    vfnmadd132pd    %zmm3, XX, %zmm0 # = GS1
.endm

.macro halley
    # In: GS0,GS1,GS2,GS3
    # Out: XX
    # No other registers used. Destroys input.
    vfmsub213pd     DT, ZETA, GS3
    vfmadd231pd     GS2, ETA, GS3
    vfmadd231pd     XX, R, GS3              # f

    vfmadd132pd     ZETA, R, GS2
    vfmadd231pd     ETA, GS1, GS2           # fp

    vmulpd          GS0, ETA, GS0
    vfmadd132pd     ZETA, GS0, GS1          # fpp

    vmulpd          GS1, GS3, GS1           # f*fpp
    # Next instruction uses exponent trick to speed up multiplication by 0.5 by using integer subtraction
    vpsubq          HALF_MASK, GS1, GS1     # 0.5*f*fpp
    vfmsub231pd     GS2, GS2, GS1           # fp*fp-0.5*f*fpp
    vmulpd          GS3, GS2, GS3           # f*fp
    vdivpd          GS1, GS3, GS3
    vsubpd          GS3, XX, XX
.endm

.macro newton
    # In: GS1,GS2,GS3
    # Out: XX
    # No other registers used. Destroys input.
    vmulpd          GS1, ETA, GS1
    vfmadd231pd     GS2, ZETA, GS1
    vmulpd          GS1, XX, XX
    vfnmadd132pd    ETA, XX, GS2
    vaddpd          R, GS1, XX
    vdivpd          XX, ONE, XX         # TODO: Hot spot
    vfnmadd231pd    GS3, ZETA, GS2
    vaddpd          GS2, DT, GS2
    vmulpd          GS2, XX, XX
.endm

###############################################################################
# Kepler Step
###############################################################################

.macro kepler_step grflag
    vmulpd          X, X, %zmm0
    vmulpd          VX, VX, %zmm1
    vfmadd231pd     Y, Y, %zmm0
    vfmadd231pd     VY, VY, %zmm1
    vfmadd231pd     Z, Z, %zmm0                 # r^2
    vfmadd231pd     VZ, VZ, %zmm1               # v^2
    vsqrtpd         %zmm0, R                    # r
    vdivpd          R, ONE, RI                  # 1/r
    vaddpd          M, M, BETA                  # 2*M
    vfmsub132pd     RI, %zmm1, BETA             # beta
    vmulpd          VX, X, ETA
    vfmadd231pd     VY, Y, ETA
    vfmadd231pd     VZ, Z, ETA                  # eta
    vmovapd         BETA, ZETA
    vfnmadd132pd    R, M, ZETA                  # zeta
    vmulpd          RI, DT, %zmm5               # dt/r
    vmulpd          ETA, %zmm5, %zmm4           # eta*dt/r
    vpsubq          HALF_MASK, %zmm4, %zmm4     # 0.5*eta*dt/r  (Note: integer sub trick)
    vfnmadd132pd    RI, ONE, %zmm4        
    vmulpd          %zmm5, %zmm4, XX            # X (second order initial guess)
   
    # Iterations to improve X
    # PYTHON REPLACE START
    mm_stiefel_Gs03_avx512
    halley
    #mm_stiefel_Gs03_avx512
    #halley
   
.NewtonLoop\grflag:    
    vmovapd         XX,     %zmm9
    mm_stiefel_Gs13_avx512
    newton 

    vsubpd      XX, %zmm9, %zmm9  # Delta XX
    vpandq      %zmm10, %zmm9, %zmm9  # abs(Delta XX)

    vcmppd     $25, %zmm11, %zmm9, %k2         # abs(Delta XX) < eps    25 = Not greater or equal, unordered (nans pass), quiet
    kmovb   %k2, %eax
    cmpb    $0xFF, %al
    jne .NewtonLoop\grflag

    mm_stiefel_Gs13_avx512
    # PYTHON REPLACE STOP
    
    # Calculate 1/r
    vmulpd          GS1, ETA, %zmm2
    vfmadd231pd     GS2, ZETA, %zmm2
    vaddpd          R, %zmm2, XX
    vdivpd          XX, ONE, %zmm4          # 1/r
    
    # Calculate f and g functions
    vmulpd          GS2, M, %zmm5
    vmulpd          RI, %zmm5, %zmm3        # negative f
    vmulpd          %zmm5, %zmm4, %zmm2     # negative gd
    vmovapd         DT, %zmm1
    vfnmadd231pd    GS3, M, %zmm1           # g 
    vmulpd          GS1, M, %zmm0
    vmulpd          RI, %zmm0, %zmm0
    vmulpd          %zmm4, %zmm0, %zmm0     # negative fd

    vmovapd         %zmm3, %zmm4
    vmovapd         %zmm3, %zmm5
    // Calculate new x y z
    vfnmadd132pd    X, X, %zmm3
    vfnmadd132pd    Y, Y, %zmm4
    vfnmadd132pd    Z, Z, %zmm5
    vfmadd231pd     VX, %zmm1, %zmm3{%k1}{z}
    vfmadd231pd     VY, %zmm1, %zmm4{%k1}{z}
    vfmadd231pd     VZ, %zmm1, %zmm5{%k1}{z}
    // Calculate new vx vy vz
    vfnmadd132pd    %zmm2, VX, VX
    vfnmadd132pd    %zmm2, VY, VY
    vfnmadd132pd    %zmm2, VZ, VZ
    vfnmadd231pd    %zmm0, X, VX{%k1}{z}
    vfnmadd231pd    %zmm0, Y, VY{%k1}{z}
    vfnmadd231pd    %zmm0, Z, VZ{%k1}{z}
    vmovapd    %zmm3, X 
    vmovapd    %zmm4, Y
    vmovapd    %zmm5, Z
.endm

###############################################################################
# Interaction Step
###############################################################################
.macro interaction_step grflag
    # TODO: Floating point error accumulation might be less if Jacobi and GR are added after P-P perturbations
    # Add Jacobi term in Jacobi coordinates
    vmulpd      X, X, %zmm4     
    vfmadd231pd Y, Y, %zmm4      
    vfmadd231pd Z, Z, %zmm4             # r^2
    vsqrtpd     %zmm4, %zmm5            # r 
    vmulpd      %zmm4, %zmm5, %zmm4     # r^3
  
    vdivpd      %zmm4, M_DT, %zmm6  # M*dt/r^3 (where M=(m0, m0+m1, m0+m1+m2,...)
    
    vfmadd231pd     X, %zmm6, VX{%k1}{z} 
    vfmadd231pd     Y, %zmm6, VY{%k1}{z} 
    vfmadd231pd     Z, %zmm6, VZ{%k1}{z} 
    
    leaq P512_MAT8_JACOBI_TO_HELIOCENTRIC(%rdi), %rax  # mat8_inertial_to_jacobi
    mat8_mul3 X, Y, Z, HX, HY, HZ
    
    # Calculating r, r^2, r^3 for Jacobi term and GR
    vmulpd      HX, HX, %zmm6
    vfmadd231pd HY, HY, %zmm6
    vfmadd231pd HZ, HZ, %zmm6               # r^2
    vsqrtpd     %zmm6, %zmm7                # r
        
    # Jacobi term
    vmovapd     P512_M0(%rdi), %zmm5        # -m0*dt
    vmulpd    %zmm6, %zmm7, %zmm7           # r^3    
    vdivpd    %zmm7, %zmm5, %zmm8{%k1}{z}   # -m0*dt/r^3 (jacobi term)

    .if \grflag == 1
        # GR term
        vmovapd     P512_GR_PREFAC(%rdi), %zmm3
        vmovapd     P512_GR_PREFAC2(%rdi), %zmm4

        vmulpd    %zmm6, %zmm6, %zmm5               # r^4
        vdivpd    %zmm5, %zmm3, %zmm7{%k1}{z}       # -dt*6*m0*m0/(c*c) /r^4

        vmulpd    %zmm7, HX, %zmm5                  # -x_j*dt*6*m0*m0/(c*c) /r^4
        vmulpd    %zmm7, HY, %zmm6
        vmulpd    %zmm7, HZ, %zmm7
        
        vmulpd    %zmm5, %zmm4, HVX{%k1}{z}         # x_j*dt*6*m0*m/(c*c) /r^4 
        vmulpd    %zmm6, %zmm4, HVY{%k1}{z}
        vmulpd    %zmm7, %zmm4, HVZ{%k1}{z}

        REDUCE_ADD_AND_BROADCAST HVX, %zmm4
        REDUCE_ADD_AND_BROADCAST HVY, %zmm4
        REDUCE_ADD_AND_BROADCAST HVZ, %zmm4

        vfmadd231pd  %zmm8, HX, HVX          # delta v_x due to Jacobi term + GR backreaction
        vfmadd231pd  %zmm8, HY, HVY
        vfmadd231pd  %zmm8, HZ, HVZ
        
        vaddpd    %zmm5, HVX, HVX            # delta v_x due to Jacobi term + GR backreaction + GR
        vaddpd    %zmm6, HVY, HVY
        vaddpd    %zmm7, HVZ, HVZ
    .else
        vmulpd    %zmm8, HX, HVX             # delta v_x due to Jacobi term, -x_j*m0*dt/r^3
        vmulpd    %zmm8, HY, HVY
        vmulpd    %zmm8, HZ, HVZ
    .endif

    #################################################################
    #// 0123 4567
    #// 3201 7645

    vmulpd  P512_m(%rdi), DT, %zmm3         # dt*m

    vpermpd $0x4B, HX, %zmm0                # 01234567 -> 32017645
    vpermpd $0x4B, HY, %zmm1
    vpermpd $0x4B, HZ, %zmm2
    vpermpd $0x4B, %zmm3, %zmm4

    vsubpd  %zmm0, HX, %zmm0                # d_x
    vsubpd  %zmm1, HY, %zmm1
    vsubpd  %zmm2, HZ, %zmm2
    
    gravity_prefactor                       # zmm6 is 1/r^3
    vmulpd      %zmm6, %zmm4, %zmm5         # dt*m/r^3
    
    vfnmadd231pd %zmm5, %zmm0,  HVX
    vfnmadd231pd %zmm5, %zmm1,  HVY
    vfnmadd231pd %zmm5, %zmm2,  HVZ

    vmulpd      %zmm6, %zmm3, %zmm5         # dt*m/r^3
    vpermpd $0x1E, %zmm0, %zmm0             # 32017645 -> 01234567
    vpermpd $0x1E, %zmm1, %zmm1
    vpermpd $0x1E, %zmm2, %zmm2
    vpermpd $0x1E, %zmm5, %zmm5

    #// 0123 4567
    #// 2310 6754
    
    vfmadd231pd %zmm5, %zmm0,  HVX
    vfmadd231pd %zmm5, %zmm1,  HVY
    vfmadd231pd %zmm5, %zmm2,  HVZ

    #################################################################
    #// 0123 4567
    #// 1032 5476
    
    vshufpd $0x55, HX, HX, %zmm0                # 01234567 -> 10325476
    vshufpd $0x55, HY, HY, %zmm1                # Using vshufpd (1 cycle) rather than vpermpd (3 cycles) 
    vshufpd $0x55, HZ, HZ, %zmm2
    vshufpd $0x55, %zmm3, %zmm3, %zmm4 

    vsubpd  %zmm0, HX, %zmm0                    # d_x
    vsubpd  %zmm1, HY, %zmm1
    vsubpd  %zmm2, HZ, %zmm2
    
    gravity_prefactor %zmm4                     # zmm6 is 1/r^3
    
    vfnmadd231pd %zmm6, %zmm0,  HVX
    vfnmadd231pd %zmm6, %zmm1,  HVY
    vfnmadd231pd %zmm6, %zmm2,  HVZ

    #################################################################
    #// 0123 4567
    #// 4567 1230
    
    vmovdqa64 b3idx(%rip), %zmm7

    vpermpd HX, %zmm7, %zmm0                    # 01234567 -> 45671230 
    vpermpd HY, %zmm7, %zmm1
    vpermpd HZ, %zmm7, %zmm2
    vpermpd %zmm3, %zmm7, %zmm4 

    vsubpd  %zmm0, HX, %zmm0                    # d_x
    vsubpd  %zmm1, HY, %zmm1
    vsubpd  %zmm2, HZ, %zmm2
    
    gravity_prefactor                           # zmm6 is 1/r^3
    vmulpd      %zmm6, %zmm4, %zmm5             # m/r^3
  
    vfnmadd231pd %zmm5, %zmm0,  HVX
    vfnmadd231pd %zmm5, %zmm1,  HVY
    vfnmadd231pd %zmm5, %zmm2,  HVZ

    vmulpd      %zmm6, %zmm3, %zmm5             # m/r^3
    
    #// 4567 1230
    #// 0123 4567
    vmulpd %zmm5, %zmm0,  HVXC
    vmulpd %zmm5, %zmm1,  HVYC
    vmulpd %zmm5, %zmm2,  HVZC

    #################################################################
    #// 0123 4567
    #// 5674 2301
    
    vmovdqa64 b4idx(%rip), %zmm7

    vpermpd HX, %zmm7, %zmm0                    # 01234567 -> 56742301  
    vpermpd HY, %zmm7, %zmm1                    # TODO: Make this an in-line shuffle be reusing block3 data
    vpermpd HZ, %zmm7, %zmm2
    vpermpd %zmm3, %zmm7, %zmm4 

    vsubpd  %zmm0, HX, %zmm0                    # d_x
    vsubpd  %zmm1, HY, %zmm1
    vsubpd  %zmm2, HZ, %zmm2
    
    gravity_prefactor                           # zmm6 is 1/r^3
    vmulpd      %zmm6, %zmm4, %zmm5             # m/r^3
  
    vfnmadd231pd %zmm5, %zmm0,  HVX
    vfnmadd231pd %zmm5, %zmm1,  HVY
    vfnmadd231pd %zmm5, %zmm2,  HVZ

    vmulpd      %zmm6, %zmm3, %zmm5             # m/r^3
    vpermpd $0x93, %zmm0, %zmm0                 # 5674 2301 -> 4567 1230
    vpermpd $0x93, %zmm1, %zmm1
    vpermpd $0x93, %zmm2, %zmm2
    vpermpd $0x93, %zmm5, %zmm5
    
    #// 4567 1230
    #// 3012 7456
    
    vfmadd231pd %zmm5, %zmm0,  HVXC
    vfmadd231pd %zmm5, %zmm1,  HVYC
    vfmadd231pd %zmm5, %zmm2,  HVZC
    
    #################################################################
    ## Final 256 bit lane crossing and add
    vmovdqa64 b34mergeidx(%rip), %zmm7

    vpermpd HVXC, %zmm7, %zmm0
    vpermpd HVYC, %zmm7, %zmm1
    vpermpd HVZC, %zmm7, %zmm2

    vaddpd %zmm0, HVX, %zmm0{%k1}{z}
    vaddpd %zmm1, HVY, %zmm1{%k1}{z}
    vaddpd %zmm2, HVZ, %zmm2{%k1}{z}

    # Convert accelerations (delta v) from heliocentric to Jacobi.
    leaq P512_MAT8_INERTIAL_TO_JACOBI(%rdi), %rax  # mat8_inertial_to_jacobi
   
    mat8_mul3 %zmm0, %zmm1, %zmm2, %zmm0, %zmm1, %zmm2

    # Update velocities
    # This could be combined with mat8_mul3.
    # However, that would increase floating point errors because sum(DVX) << VX
    vaddpd    VX, %zmm0, VX        
    vaddpd    VY, %zmm1, VY
    vaddpd    VZ, %zmm2, VZ
.endm 

###############################################################################
# Global functions
###############################################################################
reb_whfast512_kepler_step:
    # Need to init registers here when not called after interaction step.
    # This will be required for synchronizing.  
    reb_whfast512_init_registers
    kepler_step 2
    vmovapd    VX, P512_VX(%rdi)
    vmovapd    VY, P512_VY(%rdi)
    vmovapd    VZ, P512_VZ(%rdi)
    vmovapd    X, P512_X(%rdi)
    vmovapd    Y, P512_Y(%rdi)
    vmovapd    Z, P512_Z(%rdi)
    ret

# Macro creates two functions for branchless GR/no-GR
.macro BLOCK1 grflag
    # Input:   
    #           rdi = p512
    #           rsi = Number of steps (counting down)

    # Load constants
    reb_whfast512_init_registers
    # Allocate space on stack for matrix multiplications
    subq    $192, %rsp

    # Main loop
.LMainLoop\grflag:    
    kepler_step \grflag
    interaction_step \grflag
    subq    $1, %rsi
    jnz     .LMainLoop\grflag

    # Store final data in P512 structure
    vmovapd    VX, P512_VX(%rdi)
    vmovapd    VY, P512_VY(%rdi)
    vmovapd    VZ, P512_VZ(%rdi)
    vmovapd    X, P512_X(%rdi)
    vmovapd    Y, P512_Y(%rdi)
    vmovapd    Z, P512_Z(%rdi)

    addq    $192, %rsp
    ret
.endm

block1_gr: BLOCK1 1

block1_nogr: BLOCK1 0



.section    .rodata

.align 8
.HALF_MASK:  # Used for fast division by 2
    .quad   0x0010000000000000

# Shuffle Indicies
# Each is eight 64-bit integers
.align 64
b3idx:
    .quad 4,5,6,7,1,2,3,0   
b4idx:
    .quad 5,6,7,4,2,3,0,1
b34mergeidx:
    .quad 7,4,5,6,0,1,2,3

# Inverse factorial table
.align 64
.SIGN_MASK:
    .quad 0x7FFFFFFFFFFFFFFF
    .quad 0x7FFFFFFFFFFFFFFF
    .quad 0x7FFFFFFFFFFFFFFF
    .quad 0x7FFFFFFFFFFFFFFF
    .quad 0x7FFFFFFFFFFFFFFF
    .quad 0x7FFFFFFFFFFFFFFF
    .quad 0x7FFFFFFFFFFFFFFF
    .quad 0x7FFFFFFFFFFFFFFF
.align 64   
.EPS:
    .double 0.00000001
    .double 0.00000001
    .double 0.00000001
    .double 0.00000001
    .double 0.00000001
    .double 0.00000001
    .double 0.00000001
    .double 0.00000001
.align 64
.IF0:
.DOUBLE_ONE:
    .double     1.0
    .double     1.0
    .double     0.5
    .long    1431655765
    .long    1069897045
    .long    1431655765
    .long    1067799893
    .long    286331153
    .long    1065423121
    .long    381774871
    .long    1062650220
    .long    436314138
    .long    1059717536
    .long    436314138
    .long    1056571808
    .long    -1521039564
    .long    1053236707
    .long    -1216831652
    .long    1049787983
    .long    1744127204
    .long    1046144581
    .long    -268904296
    .long    1042411224
    .long    329805065
    .long    1038488134
    .long    -1463780195
    .long    1034500468
    .long    -416040929
    .long    1030416371
    .long    -416040929
    .long    1026222067
    .long    1882238282
    .long    1021924039
    .long    1673100695
    .long    1017545336
    .long    1182875991
    .long    1013118107
# Parts of inverse factorial not used at the moment:    
    .long    -1543372251
    .long    1008620587
    .long    -153291406
    .long    1003953038
    .long    1840773203
    .long    999345721
    .long    320223257
    .long    994533812
    .long    426964344
    .long    989801712
    .long    -585736279
    .long    984871884
    .long    -60141991
    .long    979930757
    .long    -1025716573
    .long    974985905
    .long    641009757
    .long    969974154
    .long    -1958520658
    .long    964844025
    .long    1346885134
    .long    959681324
    .long    -410782275
    .long    954479826
    .long    -410782275
    .long    949236946
    .long    1423773010
    .long    943953938
    .long    834731386
    .long    938635522

.section .note.GNU-stack,"",@progbits
