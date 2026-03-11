# file: integrator_whfast512.s

.section .rodata
    .align 64
.one: 
    .double 1.0

.align 64
b3idx:
    .quad 4,5,6,7,1,2,3,0   # Eight 64-bit integers

.align 64
b4idx:
    .quad 5,6,7,4,2,3,0,1

.align 64
b34mergeidx:
    .quad 7,4,5,6,0,1,2,3

matrixidx:
    .quad 0,0,0,0,0,0,0,0

.section .text
#.globl gravity_prefactor_avx512_one_test
.globl gravity_prefactor_avx512_one
.globl gravity_prefactor_avx512
.globl gr_potential
.globl block1
.globl mat8_mul3_avx512


gravity_prefactor_avx512_one:
    # Input:  zmm0=dx, zmm1=dy, zmm2=dy
    
    vmulpd      %zmm0, %zmm0, %zmm0     
    vfmadd231pd %zmm1, %zmm1, %zmm0      
    vfmadd231pd %zmm2, %zmm2, %zmm0     # zmm0 is now r^2
    
    vsqrtpd     %zmm0, %zmm1             
    vmulpd      %zmm0, %zmm1, %zmm0     # zmm0 is r^3
   
    vbroadcastsd .one(%rip), %zmm1      # Todo: keep 1 in a register at all times 
    vdivpd      %zmm0, %zmm1, %zmm0      
    
    ret                                 # return 1 / r^3 in zmm0

gravity_prefactor_avx512_one_zmm30:
    # Input:  zmm0=dx, zmm1=dy, zmm2=dy
    
    vmulpd      %zmm0, %zmm0, %zmm30     
    vfmadd231pd %zmm1, %zmm1, %zmm30      
    vfmadd231pd %zmm2, %zmm2, %zmm30     # zmm30 is now r^2
    
    vsqrtpd     %zmm30, %zmm31             
    vmulpd      %zmm30, %zmm31, %zmm30     # zmm30 is r^3
   
    vbroadcastsd .one(%rip), %zmm31      # Todo: keep 1 in a register at all times 
    vdivpd      %zmm30, %zmm31, %zmm30      
    
    ret                                 # return 1 / r^3 in zmm0
gravity_prefactor_avx512_one_test:
    # Input:  zmm0=dx, zmm1=dy, zmm2=dy
    
    vmulpd      %zmm0, %zmm0, %zmm0     
    vfmadd231pd %zmm1, %zmm1, %zmm0      
    vfmadd231pd %zmm2, %zmm2, %zmm0     # zmm0 is now r^2
    
    vsqrtpd     %zmm0, %zmm1             
    vmulpd      %zmm0, %zmm1, %zmm0     # zmm0 is r^3
   
    vbroadcastsd .one(%rip), %zmm1      # Todo: keep 1 in a register at all times 
    vdivpd      %zmm0, %zmm1, %zmm0      
    
    
    ret                                 # return 1 / r^3 in zmm0
    

gravity_prefactor_avx512:
    # Input:  zmm0=m, zmm1=dx, zmm2=dy, zmm3=dz
    
    vmulpd      %zmm1, %zmm1, %zmm1     
    vfmadd231pd %zmm2, %zmm2, %zmm1      
    vfmadd231pd %zmm3, %zmm3, %zmm1     # zmm1 is now r^2
    
    vsqrtpd     %zmm1, %zmm2             
    vmulpd      %zmm1, %zmm2, %zmm2      
    vdivpd      %zmm2, %zmm0, %zmm0      
    
    ret                                 # return m / r^3 in zmm0


.macro REDUCE_ADD_AND_BROADCAST reg, temp_reg
    vshuff64x2 $0x4E, \reg, \reg, \temp_reg         # 01234567 -> 45670123
    vaddpd     \reg, \temp_reg, \temp_reg    
    vshufpd    $0x55, \temp_reg, \temp_reg, \reg    # 01234567 -> 10325476
    vaddpd     \reg, \temp_reg, \temp_reg
    vpermpd    $0x4E, \temp_reg, \reg               # 01234567 -> 23016745
    vaddpd     \reg, \temp_reg, \reg
.endm
        
mat8_mul3_avx512:
    # 8x8 matrix multiplied with 3 different 8 vectors
    # in: rdi = vector to 64 matrix elements
    # zmm0, zmm1, zmm2  input vectors
    # rsi, rdx, rcx     output vectors 
	vmovapd	(%rdi), %zmm4
	vbroadcastsd	%xmm0, %zmm3
	movq	%rsi, %rax
	vmulpd	%zmm4, %zmm3, %zmm3
	vmovapd	%zmm3, (%rsi)
	vbroadcastsd	%xmm1, %zmm3
	movl	$1, %esi
	vmulpd	%zmm4, %zmm3, %zmm3
	vmovapd	%zmm3, (%rdx)
	vbroadcastsd	%xmm2, %zmm3
	vmulpd	%zmm4, %zmm3, %zmm3
	vmovapd	%zmm3, (%rcx)
	vpbroadcastq	%rsi, %zmm3
	vmovapd	64(%rdi), %zmm4
	movl	$2, %esi
	vpermpd	%zmm0, %zmm3, %zmm5
	vfmadd213pd	(%rax), %zmm4, %zmm5
	vmovapd	%zmm5, (%rax)
	vpermpd	%zmm1, %zmm3, %zmm5
	vfmadd213pd	(%rdx), %zmm4, %zmm5
	vpermpd	%zmm2, %zmm3, %zmm3
	vmovapd	%zmm5, (%rdx)
	vfmadd213pd	(%rcx), %zmm4, %zmm3
	vmovapd	%zmm3, (%rcx)
	vpbroadcastq	%rsi, %zmm3
	vmovapd	128(%rdi), %zmm4
	movl	$3, %esi
	vpermpd	%zmm0, %zmm3, %zmm5
	vfmadd213pd	(%rax), %zmm4, %zmm5
	vmovapd	%zmm5, (%rax)
	vpermpd	%zmm1, %zmm3, %zmm5
	vfmadd213pd	(%rdx), %zmm4, %zmm5
	vpermpd	%zmm2, %zmm3, %zmm3
	vmovapd	%zmm5, (%rdx)
	vfmadd213pd	(%rcx), %zmm4, %zmm3
	vmovapd	%zmm3, (%rcx)
	vpbroadcastq	%rsi, %zmm3
	vmovapd	192(%rdi), %zmm4
	movl	$4, %esi
	vpermpd	%zmm0, %zmm3, %zmm5
	vfmadd213pd	(%rax), %zmm4, %zmm5
	vmovapd	%zmm5, (%rax)
	vpermpd	%zmm1, %zmm3, %zmm5
	vfmadd213pd	(%rdx), %zmm4, %zmm5
	vpermpd	%zmm2, %zmm3, %zmm3
	vmovapd	%zmm5, (%rdx)
	vfmadd213pd	(%rcx), %zmm4, %zmm3
	vmovapd	%zmm3, (%rcx)
	vpbroadcastq	%rsi, %zmm3
	vmovapd	256(%rdi), %zmm4
	movl	$5, %esi
	vpermpd	%zmm0, %zmm3, %zmm5
	vfmadd213pd	(%rax), %zmm4, %zmm5
	vmovapd	%zmm5, (%rax)
	vpermpd	%zmm1, %zmm3, %zmm5
	vfmadd213pd	(%rdx), %zmm4, %zmm5
	vpermpd	%zmm2, %zmm3, %zmm3
	vmovapd	%zmm5, (%rdx)
	vfmadd213pd	(%rcx), %zmm4, %zmm3
	vmovapd	%zmm3, (%rcx)
	vpbroadcastq	%rsi, %zmm3
	vmovapd	320(%rdi), %zmm4
	movl	$6, %esi
	vpermpd	%zmm0, %zmm3, %zmm5
	vfmadd213pd	(%rax), %zmm4, %zmm5
	vmovapd	%zmm5, (%rax)
	vpermpd	%zmm1, %zmm3, %zmm5
	vfmadd213pd	(%rdx), %zmm4, %zmm5
	vpermpd	%zmm2, %zmm3, %zmm3
	vmovapd	%zmm5, (%rdx)
	vfmadd213pd	(%rcx), %zmm4, %zmm3
	vmovapd	%zmm3, (%rcx)
	vpbroadcastq	%rsi, %zmm3
	vmovapd	384(%rdi), %zmm4
	movl	$7, %esi
	vpermpd	%zmm0, %zmm3, %zmm5
	vfmadd213pd	(%rax), %zmm4, %zmm5
	vmovapd	%zmm5, (%rax)
	vpermpd	%zmm1, %zmm3, %zmm5
	vfmadd213pd	(%rdx), %zmm4, %zmm5
	vpermpd	%zmm2, %zmm3, %zmm3
	vmovapd	%zmm5, (%rdx)
	vfmadd213pd	(%rcx), %zmm4, %zmm3
	vmovapd	%zmm3, (%rcx)
	vpbroadcastq	%rsi, %zmm3
	vmovapd	448(%rdi), %zmm4
	vpermpd	%zmm0, %zmm3, %zmm0
	vfmadd213pd	(%rax), %zmm4, %zmm0
	vpermpd	%zmm1, %zmm3, %zmm1
	vpermpd	%zmm2, %zmm3, %zmm2
	vmovapd	%zmm0, (%rax)
	vfmadd213pd	(%rdx), %zmm4, %zmm1
	vmovapd	%zmm1, (%rdx)
	vfmadd213pd	(%rcx), %zmm4, %zmm2
	vmovapd	%zmm2, (%rcx)
    ret

mat8_mul3_avx512_nomem:
    # 8x8 matrix multiplied with 3 different 8 vectors
    # in: rdi = vector to 64 matrix elements
    # zmm0, zmm1, zmm2  input vectors
    # zmm10, zmm11, zmm12  output vectors
	vmovapd	(%rdi), %zmm4
	vbroadcastsd	%xmm0, %zmm3
	vmulpd	%zmm4, %zmm3, %zmm3
	vmovapd	%zmm3, %zmm10
	vbroadcastsd	%xmm1, %zmm3
	movl	$1, %esi
	vmulpd	%zmm4, %zmm3, %zmm3
	vmovapd	%zmm3, %zmm11
	vbroadcastsd	%xmm2, %zmm3
	vmulpd	%zmm4, %zmm3, %zmm3
	vmovapd	%zmm3, %zmm12
	vpbroadcastq	%rsi, %zmm3
	vmovapd	64(%rdi), %zmm4
	movl	$2, %esi
	vpermpd	%zmm0, %zmm3, %zmm5
	vfmadd213pd	%zmm10, %zmm4, %zmm5
	vmovapd	%zmm5, %zmm10  #TODO Remove mov
	vpermpd	%zmm1, %zmm3, %zmm5
	vfmadd213pd	%zmm11, %zmm4, %zmm5
	vpermpd	%zmm2, %zmm3, %zmm3
	vmovapd	%zmm5, %zmm11
	vfmadd213pd	%zmm12, %zmm4, %zmm3
	vmovapd	%zmm3, %zmm12
	vpbroadcastq	%rsi, %zmm3
	vmovapd	128(%rdi), %zmm4
	movl	$3, %esi
	vpermpd	%zmm0, %zmm3, %zmm5
	vfmadd213pd	%zmm10, %zmm4, %zmm5
	vmovapd	%zmm5, %zmm10
	vpermpd	%zmm1, %zmm3, %zmm5
	vfmadd213pd	%zmm11, %zmm4, %zmm5
	vpermpd	%zmm2, %zmm3, %zmm3
	vmovapd	%zmm5, %zmm11
	vfmadd213pd	%zmm12, %zmm4, %zmm3
	vmovapd	%zmm3, %zmm12
	vpbroadcastq	%rsi, %zmm3
	vmovapd	192(%rdi), %zmm4
	movl	$4, %esi
	vpermpd	%zmm0, %zmm3, %zmm5
	vfmadd213pd	%zmm10, %zmm4, %zmm5
	vmovapd	%zmm5, %zmm10
	vpermpd	%zmm1, %zmm3, %zmm5
	vfmadd213pd	%zmm11, %zmm4, %zmm5
	vpermpd	%zmm2, %zmm3, %zmm3
	vmovapd	%zmm5, %zmm11
	vfmadd213pd	%zmm12, %zmm4, %zmm3
	vmovapd	%zmm3, %zmm12
	vpbroadcastq	%rsi, %zmm3
	vmovapd	256(%rdi), %zmm4
	movl	$5, %esi
	vpermpd	%zmm0, %zmm3, %zmm5
	vfmadd213pd	%zmm10, %zmm4, %zmm5
	vmovapd	%zmm5, %zmm10
	vpermpd	%zmm1, %zmm3, %zmm5
	vfmadd213pd	%zmm11, %zmm4, %zmm5
	vpermpd	%zmm2, %zmm3, %zmm3
	vmovapd	%zmm5, %zmm11
	vfmadd213pd	%zmm12, %zmm4, %zmm3
	vmovapd	%zmm3, %zmm12
	vpbroadcastq	%rsi, %zmm3
	vmovapd	320(%rdi), %zmm4
	movl	$6, %esi
	vpermpd	%zmm0, %zmm3, %zmm5
	vfmadd213pd	%zmm10, %zmm4, %zmm5
	vmovapd	%zmm5, %zmm10
	vpermpd	%zmm1, %zmm3, %zmm5
	vfmadd213pd	%zmm11, %zmm4, %zmm5
	vpermpd	%zmm2, %zmm3, %zmm3
	vmovapd	%zmm5, %zmm11
	vfmadd213pd	%zmm12, %zmm4, %zmm3
	vmovapd	%zmm3, %zmm12
	vpbroadcastq	%rsi, %zmm3
	vmovapd	384(%rdi), %zmm4
	movl	$7, %esi
	vpermpd	%zmm0, %zmm3, %zmm5
	vfmadd213pd	%zmm10, %zmm4, %zmm5
	vmovapd	%zmm5, %zmm10
	vpermpd	%zmm1, %zmm3, %zmm5
	vfmadd213pd	%zmm11, %zmm4, %zmm5
	vpermpd	%zmm2, %zmm3, %zmm3
	vmovapd	%zmm5, %zmm11
	vfmadd213pd	%zmm12, %zmm4, %zmm3
	vmovapd	%zmm3, %zmm12
	vpbroadcastq	%rsi, %zmm3
	vmovapd	448(%rdi), %zmm4
	vpermpd	%zmm0, %zmm3, %zmm0
	vfmadd213pd	%zmm10, %zmm4, %zmm0
	vpermpd	%zmm1, %zmm3, %zmm1
	vpermpd	%zmm2, %zmm3, %zmm2
	vfmadd213pd	%zmm11, %zmm4, %zmm1
	vfmadd213pd	%zmm12, %zmm4, %zmm2
    ret


.set P512_DT, 64
.set P512_GR_PREFAC, 128
.set P512_GR_PREFAC2, 192
.set P512_M, 320
.set P512_HVX, 960
.set P512_HVY, 1024
.set P512_HVZ, 1088
.set P512_M0, 2688
.set P512_MASK, 2752

gr_potential:
    # Input:    zmm0=x_j, zmm1=y_j, zmm2=z_j
    #           zmm3 = gr_prefac, zmm4 = gr_prefac2
    #           zmm5 = -m0*dt
    #           edi = mask
    #           rsi=&hvx , rdx=&hvy, rcx=&hvz 
    kmovw   P512_MASK(%rdi), %k1             # mask
    vmovapd     P512_GR_PREFAC(%rdi), %zmm3
    vmovapd     P512_GR_PREFAC2(%rdi), %zmm4
    vmovapd     P512_M0(%rdi), %zmm5

    vmulpd    %zmm0, %zmm0, %zmm6
    vfmadd231pd    %zmm1, %zmm1, %zmm6
    vfmadd231pd    %zmm2, %zmm2, %zmm6     # r^2
    vsqrtpd    %zmm6, %zmm18               # r
    vmulpd    %zmm6, %zmm18, %zmm18       # r^3    
    vdivpd    %zmm18, %zmm5, %zmm8{%k1}{z}  # -m0*dt/r^3 (jacobi term)
    vmulpd    %zmm8, %zmm0, %zmm20          # -x_j*m0*dt/r^3
    vmulpd    %zmm8, %zmm1, %zmm21
    vmulpd    %zmm8, %zmm2, %zmm22

    vmulpd    %zmm6, %zmm6, %zmm7                 # r^4
    vdivpd    %zmm7, %zmm3, %zmm9{%k1}{z}         # -dt*6*m0*m0/(c*c) /r^4

    vmulpd    %zmm9, %zmm0, %zmm15                # -x_j*dt*6*m0*m0/(c*c) /r^4
    vmulpd    %zmm9, %zmm1, %zmm16
    vmulpd    %zmm9, %zmm2, %zmm17

    vmulpd    %zmm15, %zmm4, %zmm10{%k1}{z}       # x_j*dt*6*m0*m/(c*c) /r^4 
    vmulpd    %zmm16, %zmm4, %zmm11{%k1}{z}
    vmulpd    %zmm17, %zmm4, %zmm12{%k1}{z}

    REDUCE_ADD_AND_BROADCAST %zmm10, %zmm18   # sum
    REDUCE_ADD_AND_BROADCAST %zmm11, %zmm18
    REDUCE_ADD_AND_BROADCAST %zmm12, %zmm18

    vaddpd    %zmm10, %zmm15, %zmm10      # delta v_x due to gr TODO: Combine this with previous vmulpd
    vaddpd    %zmm11, %zmm16, %zmm11
    vaddpd    %zmm12, %zmm17, %zmm12
    
    vaddpd    %zmm10, %zmm20, %zmm10      # delta v_x due to gr + jacobi term
    vaddpd    %zmm11, %zmm21, %zmm11
    vaddpd    %zmm12, %zmm22, %zmm12

    vmovapd    %zmm10, P512_HVX(%rdi)              # TODO get rid of mov instruction
    vmovapd    %zmm11, P512_HVY(%rdi)
    vmovapd    %zmm12, P512_HVZ(%rdi)
    
    ret


block1:
    # Input:    zmm0=x_j, zmm1=y_j, zmm2=z_j
    #           rdi = p512
    #// 0123 4567
    #// 3201 7645
    
    kmovw   P512_MASK(%rdi), %k1             # mask

    vmovapd P512_DT(%rdi), %zmm3                 # dt
    vmulpd  P512_M(%rdi), %zmm3, %zmm3         # dt*m
    
    vmovapd %zmm0, %zmm23           # x
    vmovapd %zmm1, %zmm24
    vmovapd %zmm2, %zmm25

    vpermpd $0x4B, %zmm23, %zmm4               # 01234567 -> 32017645
    vpermpd $0x4B, %zmm24, %zmm5
    vpermpd $0x4B, %zmm25, %zmm6
    vpermpd $0x4B, %zmm3, %zmm15 

    vsubpd  %zmm4, %zmm23, %zmm0                # d_x
    vsubpd  %zmm5, %zmm24, %zmm1
    vsubpd  %zmm6, %zmm25, %zmm2
    
    call gravity_prefactor_avx512_one_zmm30  # zmm30 is 1/r^3
    vmulpd      %zmm30, %zmm15, %zmm14      # m/r^3
    
    vmovapd    960(%rdi),  %zmm10             # TODO get rid of mov instruction
    vmovapd    1024(%rdi),  %zmm11 
    vmovapd    1088(%rdi),  %zmm12
  
    vfnmadd231pd %zmm14, %zmm0,  %zmm10
    vfnmadd231pd %zmm14, %zmm1,  %zmm11
    vfnmadd231pd %zmm14, %zmm2,  %zmm12

    
    vmovapd    %zmm10, 960(%rdi)              # TODO get rid of mov instruction
    vmovapd    %zmm11, 1024(%rdi)
    vmovapd    %zmm12, 1088(%rdi)

    vpermpd $0x1E, %zmm0, %zmm0               # 32017645 -> 01234567
    vpermpd $0x1E, %zmm1, %zmm1
    vpermpd $0x1E, %zmm2, %zmm2
    vpermpd $0x1E, %zmm30, %zmm30
    vpermpd $0x1E, %zmm3, %zmm15                # 01234567 -> 32017645


    vmulpd      %zmm30, %zmm15, %zmm14      # m/r^3
    
    vmovapd    960(%rdi),  %zmm10             # TODO get rid of mov instruction
    vmovapd    1024(%rdi),  %zmm11 
    vmovapd    1088(%rdi),  %zmm12 

    #// 0123 4567
    #// 2310 6754
    
    vfmadd231pd %zmm14, %zmm0,  %zmm10
    vfmadd231pd %zmm14, %zmm1,  %zmm11
    vfmadd231pd %zmm14, %zmm2,  %zmm12

    vmovapd    %zmm10, 960(%rdi)              # TODO get rid of mov instruction
    vmovapd    %zmm11, 1024(%rdi)
    vmovapd    %zmm12, 1088(%rdi)

    #// 0123 4567
    #// 1032 5476
    
    vpermpd $0xB1, %zmm23, %zmm4               # 01234567 -> 10325476
    vpermpd $0xB1, %zmm24, %zmm5               # TODO: Use vshufpd instead here
    vpermpd $0xB1, %zmm25, %zmm6
    vpermpd $0xB1, %zmm3, %zmm15 

    vsubpd  %zmm4, %zmm23, %zmm0                # d_x
    vsubpd  %zmm5, %zmm24, %zmm1
    vsubpd  %zmm6, %zmm25, %zmm2
    
    call gravity_prefactor_avx512_one_zmm30  # zmm30 is 1/r^3
    vmulpd      %zmm30, %zmm15, %zmm14      # m/r^3
    
    vmovapd    960(%rdi),  %zmm10             # TODO get rid of mov instruction
    vmovapd    1024(%rdi),  %zmm11 
    vmovapd    1088(%rdi),  %zmm12
  
    vfnmadd231pd %zmm14, %zmm0,  %zmm10
    vfnmadd231pd %zmm14, %zmm1,  %zmm11
    vfnmadd231pd %zmm14, %zmm2,  %zmm12

    
    vmovapd    %zmm10, 960(%rdi)              # TODO get rid of mov instruction
    vmovapd    %zmm11, 1024(%rdi)
    vmovapd    %zmm12, 1088(%rdi)

    #// 0123 4567
    #// 4567 1230
    
    vmovdqa64 b3idx(%rip), %zmm18

    vpermpd %zmm23, %zmm18, %zmm4               # 01234567 -> 45671230 
    vpermpd %zmm24, %zmm18, %zmm5
    vpermpd %zmm25, %zmm18, %zmm6
    vpermpd %zmm3, %zmm18, %zmm15 

    vsubpd  %zmm4, %zmm23, %zmm0                # d_x
    vsubpd  %zmm5, %zmm24, %zmm1
    vsubpd  %zmm6, %zmm25, %zmm2
    

    call gravity_prefactor_avx512_one_zmm30  # zmm30 is 1/r^3
    vmulpd      %zmm30, %zmm15, %zmm14      # m/r^3
    
    vmovapd    960(%rdi),  %zmm10             # TODO get rid of mov instruction
    vmovapd    1024(%rdi),  %zmm11 
    vmovapd    1088(%rdi),  %zmm12
  
    vfnmadd231pd %zmm14, %zmm0,  %zmm10
    vfnmadd231pd %zmm14, %zmm1,  %zmm11
    vfnmadd231pd %zmm14, %zmm2,  %zmm12

    
    vmovapd    %zmm10, 960(%rdi)              # TODO get rid of mov instruction
    vmovapd    %zmm11, 1024(%rdi)
    vmovapd    %zmm12, 1088(%rdi)


    vmulpd      %zmm30, %zmm3, %zmm14      # m/r^3
    
    #// 4567 1230
    #// 0123 4567
    vmulpd %zmm14, %zmm0,  %zmm20
    vmulpd %zmm14, %zmm1,  %zmm21
    vmulpd %zmm14, %zmm2,  %zmm22


    #// 0123 4567
    #// 5674 2301
    
    vmovdqa64 b4idx(%rip), %zmm18

    vpermpd %zmm23, %zmm18, %zmm4               # 01234567 -> 56742301  
    vpermpd %zmm24, %zmm18, %zmm5                # TODO: Make this an in-line shuffle be reusing block3 data
    vpermpd %zmm25, %zmm18, %zmm6
    vpermpd %zmm3, %zmm18, %zmm15 

    vsubpd  %zmm4, %zmm23, %zmm0                # d_x
    vsubpd  %zmm5, %zmm24, %zmm1
    vsubpd  %zmm6, %zmm25, %zmm2
    

    call gravity_prefactor_avx512_one_zmm30  # zmm30 is 1/r^3
    vmulpd      %zmm30, %zmm15, %zmm14      # m/r^3
    
    vmovapd    960(%rdi),  %zmm10             # TODO get rid of mov instruction
    vmovapd    1024(%rdi),  %zmm11 
    vmovapd    1088(%rdi),  %zmm12
  
    vfnmadd231pd %zmm14, %zmm0,  %zmm10
    vfnmadd231pd %zmm14, %zmm1,  %zmm11
    vfnmadd231pd %zmm14, %zmm2,  %zmm12

    vmovapd    %zmm10, 960(%rdi)              # TODO get rid of mov instruction
    vmovapd    %zmm11, 1024(%rdi)
    vmovapd    %zmm12, 1088(%rdi)

    vpermpd $0x93, %zmm0, %zmm0               # 5674 2301 -> 4567 1230
    vpermpd $0x93, %zmm1, %zmm1
    vpermpd $0x93, %zmm2, %zmm2
    vpermpd $0x93, %zmm30, %zmm30
    vpermpd $0x93, %zmm3, %zmm15            


    vmulpd      %zmm30, %zmm15, %zmm14      # m/r^3
    
    #// 4567 1230
    #// 3012 7456
    
    vfmadd231pd %zmm14, %zmm0,  %zmm20
    vfmadd231pd %zmm14, %zmm1,  %zmm21
    vfmadd231pd %zmm14, %zmm2,  %zmm22
    
    ## Final 256 bit lane crossing and add
    vmovdqa64 b34mergeidx(%rip), %zmm18

    vpermpd %zmm20, %zmm18, %zmm10{%k1}{z}
    vpermpd %zmm21, %zmm18, %zmm11
    vpermpd %zmm22, %zmm18, %zmm12
    vmovapd    960(%rdi),  %zmm13             # TODO get rid of mov instruction
    vmovapd    1024(%rdi),  %zmm14 
    vmovapd    1088(%rdi),  %zmm15

    vaddpd %zmm10, %zmm13, %zmm0{%k1}{z}
    vaddpd %zmm11, %zmm14, %zmm1{%k1}{z}
    vaddpd %zmm12, %zmm15, %zmm2{%k1}{z}


    # Convert accelerations (delta v) from heliocentric to Jacobi.
	movq	%rdi, %rax
    leaq 1152(%rdi), %rdi  # mat8_inertial_to_jacobi
   
    call mat8_mul3_avx512_nomem

    # Update velocities
    vaddpd    576(%rax), %zmm0, %zmm3     #vx TODO get rid of memory
    vaddpd    640(%rax), %zmm1, %zmm4
    vaddpd    704(%rax), %zmm2, %zmm5
    

    # Add Jacobi term in Jacobi coordinates

    vmovapd    384(%rax),  %zmm0             # x  TODO get rid of mov instruction
    vmovapd    448(%rax),  %zmm1 
    vmovapd    512(%rax),  %zmm2

    lfence
    call gravity_prefactor_avx512_one_test
# The following is the same code but much slower than the function call.
# I do not understand why.    
# Also much slower when I change output register to zmm30???
#    vmulpd      %zmm0, %zmm0, %zmm0     
#    vfmadd231pd %zmm1, %zmm1, %zmm0      
#    vfmadd231pd %zmm2, %zmm2, %zmm0     # zmm0 is now r^2
#    
#    vsqrtpd     %zmm0, %zmm1             
#    vmulpd      %zmm0, %zmm1, %zmm0     # zmm0 is r^3
#   
#    vbroadcastsd .one(%rip), %zmm1      # Todo: keep 1 in a register at all times 
#    vdivpd      %zmm0, %zmm1, %zmm0      
    


    vmulpd  (%rax), %zmm0, %zmm0        # 1/r^3*M (where M=(m0, m0+m1, m0+m1+m2,...)
    vmulpd  64(%rax), %zmm0, %zmm7        # dt*1/r^3*M
    
    vmovapd    384(%rax),  %zmm0             # x  TODO get rid of mov instruction
    vmovapd    448(%rax),  %zmm1 
    vmovapd    512(%rax),  %zmm2

    vfmadd231pd     %zmm0, %zmm7, %zmm3{%k1}{z} 
    vfmadd231pd     %zmm1, %zmm7, %zmm4{%k1}{z} 
    vfmadd231pd     %zmm2, %zmm7, %zmm5{%k1}{z} 
    
    # Store final new velocities in Jacobi coordinates
    vmovapd    %zmm3, 576(%rax)              # vx TODO get rid of mov instruction
    vmovapd    %zmm4, 640(%rax)
    vmovapd    %zmm5, 704(%rax)
    
    ret


.section .note.GNU-stack,"",@progbits
