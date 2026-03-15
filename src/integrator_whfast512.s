# file: integrator_whfast512.s
.include "header.s"

.section .rodata

.align 8
.HALF_MASK:  # Used for fast division by 2
    .quad   0x0010000000000000

.align 64
.DOUBLE_ONE:
    .double 1.0

# Shuffle Indicies
# Each is eight 64-bit integers
.align 64
b3idx:
    .quad 4,5,6,7,1,2,3,0   

.align 64
b4idx:
    .quad 5,6,7,4,2,3,0,1

.align 64
b34mergeidx:
    .quad 7,4,5,6,0,1,2,3

.section .text
.globl block1_gr
.globl block1_nogr
.extern reb_whfast512_kepler_step_noinit
.extern local_reb_whfast512_kepler_step


.globl reb_whfast512_init_registers
reb_whfast512_init_registers:
    kmovw           P512_MASK(%rdi), %k1
    vmovapd         P512_DT(%rdi), DT
    vmovapd         P512_M(%rdi), M
    vbroadcastsd    .DOUBLE_ONE(%rip), ONE
    vpbroadcastq    .HALF_MASK(%rip), HALF_MASK
    
    vmovapd     P512_X(%rdi), X
    vmovapd     P512_Y(%rdi), Y
    vmovapd     P512_Z(%rdi), Z
    
    vmovapd     P512_VX(%rdi), VX
    vmovapd     P512_VY(%rdi), VY
    vmovapd     P512_VZ(%rdi), VZ
    
    ret


gravity_prefactor:
    # Input:  zmm0=dx, zmm1=dy, zmm2=dy
    
    vmulpd      %zmm0, %zmm0, %zmm6     
    vfmadd231pd %zmm1, %zmm1, %zmm6      
    vfmadd231pd %zmm2, %zmm2, %zmm6     # zmm6 is now r^2
    
    vsqrtpd     %zmm6, %zmm7             
    vmulpd      %zmm6, %zmm7, %zmm6     # zmm6 is r^3
   
    vdivpd      %zmm6, ONE, %zmm6      
    
    ret                                 # return 1 / r^3 in zmm0


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
    # uses: zmm3-zmm10
    # The idea is to use embedded broadcast loads
    # Note: matrix needs to be transposed.
    vmovupd \in0,   0(%rsp)
    vmovupd \in1,  64(%rsp)
    vmovupd \in2, 128(%rsp)

    # Preload all matrix elements into registers.
    # One register at a time might be just as fast?
    vmovapd     (%rax), %zmm3
    vmovapd     64(%rax), %zmm4
    vmovapd     128(%rax), %zmm5
    vmovapd     192(%rax), %zmm6
    vmovapd     256(%rax), %zmm7
    vmovapd     320(%rax), %zmm8
    vmovapd     384(%rax), %zmm9
    vmovapd     448(%rax), %zmm10

    # Keeping three independent FMA chains going
    vmulpd       0(%rsp){1to8}, %zmm3, \out0
    vmulpd      64(%rsp){1to8}, %zmm3, \out1
    vmulpd     128(%rsp){1to8}, %zmm3, \out2

    vfmadd231pd   8(%rsp){1to8}, %zmm4, \out0
    vfmadd231pd  72(%rsp){1to8}, %zmm4, \out1
    vfmadd231pd 136(%rsp){1to8}, %zmm4, \out2

    vfmadd231pd  16(%rsp){1to8}, %zmm5, \out0
    vfmadd231pd  80(%rsp){1to8}, %zmm5, \out1
    vfmadd231pd 144(%rsp){1to8}, %zmm5, \out2
    
    vfmadd231pd  24(%rsp){1to8}, %zmm6, \out0
    vfmadd231pd  88(%rsp){1to8}, %zmm6, \out1
    vfmadd231pd 152(%rsp){1to8}, %zmm6, \out2
    
    vfmadd231pd  32(%rsp){1to8}, %zmm7, \out0
    vfmadd231pd  96(%rsp){1to8}, %zmm7, \out1
    vfmadd231pd 160(%rsp){1to8}, %zmm7, \out2
    
    vfmadd231pd  40(%rsp){1to8}, %zmm8, \out0
    vfmadd231pd 104(%rsp){1to8}, %zmm8, \out1
    vfmadd231pd 168(%rsp){1to8}, %zmm8, \out2
    
    vfmadd231pd  48(%rsp){1to8}, %zmm9, \out0
    vfmadd231pd 112(%rsp){1to8}, %zmm9, \out1
    vfmadd231pd 176(%rsp){1to8}, %zmm9, \out2
    
    vfmadd231pd  56(%rsp){1to8}, %zmm10, \out0
    vfmadd231pd 120(%rsp){1to8}, %zmm10, \out1
    vfmadd231pd 184(%rsp){1to8}, %zmm10, \out2

   .endm


.macro BLOCK1 grflag
    # Input:   
    #           rdi = p512
    #           rsi = Number of steps (counting down)

    # Load Constant
    call reb_whfast512_init_registers
    # Allocate space on stack for matrix multiplications
    subq    $192, %rsp

####################################    
#   Start Main Loop
####################################
.LMainLoop\grflag:    

    call reb_whfast512_kepler_step_noinit

    # Interaction step:
    # Add Jacobi term in Jacobi coordinates
    vmulpd      X, X, %zmm4     
    vfmadd231pd Y, Y, %zmm4      
    vfmadd231pd Z, Z, %zmm4     # r^2
    vsqrtpd     %zmm4, %zmm5            # r 
    vmulpd      %zmm4, %zmm5, %zmm4     # r^3
  
    vdivpd      %zmm4, M, %zmm4  # 1/r^3*M (where M=(m0, m0+m1, m0+m1+m2,...)
    vmulpd      DT, %zmm4, %zmm6        # dt*1/r^3*M
    
    vfmadd231pd     X, %zmm6, VX{%k1}{z} 
    vfmadd231pd     Y, %zmm6, VY{%k1}{z} 
    vfmadd231pd     Z, %zmm6, VZ{%k1}{z} 
    
    leaq P512_MAT8_JACOBI_TO_HELIOCENTRIC(%rdi), %rax  # mat8_inertial_to_jacobi
    mat8_mul3 X, Y, Z, HX, HY, HZ
    
    # Calculating r, r^2, r^3 for Jacobi term and GR
    vmulpd      HX, HX, %zmm6
    vfmadd231pd HY, HY, %zmm6
    vfmadd231pd HZ, HZ, %zmm6         # r^2
    vsqrtpd     %zmm6, %zmm7                # r
        
    # Jacobi term
    vmovapd     P512_M0(%rdi), %zmm5   # -m0*dt
    vmulpd    %zmm6, %zmm7, %zmm7           # r^3    
    vdivpd    %zmm7, %zmm5, %zmm8{%k1}{z}   # -m0*dt/r^3 (jacobi term)

    .if \grflag == 1
        vmovapd     P512_GR_PREFAC(%rdi), %zmm3
        vmovapd     P512_GR_PREFAC2(%rdi), %zmm4

        vmulpd    %zmm6, %zmm6, %zmm5               # r^4
        vdivpd    %zmm5, %zmm3, %zmm7{%k1}{z}       # -dt*6*m0*m0/(c*c) /r^4

        vmulpd    %zmm7, HX, %zmm5               # -x_j*dt*6*m0*m0/(c*c) /r^4
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
        
        vaddpd    %zmm5, HVX, HVX               # delta v_x due to Jacobi term + GR backreaction + GR
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
    
    call gravity_prefactor                  # zmm6 is 1/r^3
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
    
    # TODO: Combine 1/r with multiplication
    call gravity_prefactor                      # zmm6 is 1/r^3
    vmulpd      %zmm6, %zmm4, %zmm5             # m/r^3
    
    vfnmadd231pd %zmm5, %zmm0,  HVX
    vfnmadd231pd %zmm5, %zmm1,  HVY
    vfnmadd231pd %zmm5, %zmm2,  HVZ

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
    
    call gravity_prefactor                      # zmm6 is 1/r^3
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
    
    call gravity_prefactor                      # zmm6 is 1/r^3
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

    subq    $1, %rsi
    jnz     .LMainLoop\grflag
####################################    
#   End  Main Loop
####################################


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
  

.section .note.GNU-stack,"",@progbits
