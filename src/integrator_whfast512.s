# file: integrator_whfast512.s

.section .rodata
    .align 64
.one: 
    .double 1.0
.t2: 
    .double 2.0
    .double 2.0
    .double 2.0
    .double 2.0
    .double 2.0
    .double 2.0
    .double 2.0
    .double 2.0
.t3: 
    .double 3.0
    .double 3.0
    .double 3.0
    .double 3.0
    .double 3.0
    .double 3.0
    .double 3.0
    .double 3.0
.t5: 
    .double 5.0
    .double 5.0
    .double 5.0
    .double 5.0
    .double 5.0
    .double 5.0
    .double 5.0
    .double 5.0



.section .text
.globl gravity_prefactor_avx512_one
.globl gravity_prefactor_avx512
.globl gr_potential
.globl block1


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
        

gr_potential:
    # Input:    zmm0=x_j, zmm1=y_j, zmm2=z_j
    #           zmm3 = gr_prefac, zmm4 = gr_prefac2
    #           zmm5 = -m0*dt
    #           edi = mask
    #           rsi=&hvx , rdx=&hvy, rcx=&hvz 
    kmovw	%edi, %k1

	vmulpd	%zmm0, %zmm0, %zmm6
	vfmadd231pd	%zmm1, %zmm1, %zmm6
	vfmadd231pd	%zmm2, %zmm2, %zmm6     # r^2
    vsqrtpd	%zmm6, %zmm18               # r
	vmulpd	%zmm6, %zmm18, %zmm18       # r^3    
	vdivpd	%zmm18, %zmm5, %zmm8{%k1}{z}  # -m0*dt/r^3 (jacobi term)
	vmulpd	%zmm8, %zmm0, %zmm20          # -x_j*m0*dt/r^3
	vmulpd	%zmm8, %zmm1, %zmm21
	vmulpd	%zmm8, %zmm2, %zmm22

	vmulpd	%zmm6, %zmm6, %zmm7                 # r^4
	vdivpd	%zmm7, %zmm3, %zmm9{%k1}{z}         # -dt*6*m0*m0/(c*c) /r^4

	vmulpd	%zmm9, %zmm0, %zmm15                # -x_j*dt*6*m0*m0/(c*c) /r^4
	vmulpd	%zmm9, %zmm1, %zmm16
	vmulpd	%zmm9, %zmm2, %zmm17

	vmulpd	%zmm15, %zmm4, %zmm10{%k1}{z}       # x_j*dt*6*m0*m/(c*c) /r^4 
	vmulpd	%zmm16, %zmm4, %zmm11{%k1}{z}
	vmulpd	%zmm17, %zmm4, %zmm12{%k1}{z}

    REDUCE_ADD_AND_BROADCAST %zmm10, %zmm18   # sum
    REDUCE_ADD_AND_BROADCAST %zmm11, %zmm18
    REDUCE_ADD_AND_BROADCAST %zmm12, %zmm18

	vaddpd	%zmm10, %zmm15, %zmm10      # delta v_x due to gr TODO: Combine this with previous vmulpd
    vaddpd	%zmm11, %zmm16, %zmm11
	vaddpd	%zmm12, %zmm17, %zmm12
	
    vaddpd	%zmm10, %zmm20, %zmm10      # delta v_x due to gr + jacobi term
    vaddpd	%zmm11, %zmm21, %zmm11
	vaddpd	%zmm12, %zmm22, %zmm12

	vmovapd	%zmm10, (%rsi)              # TODO get rid of mov instruction
	vmovapd	%zmm11, (%rdx)
	vmovapd	%zmm12, (%rcx)
	
    ret


block1:
    # Input:    zmm0=x_j, zmm1=y_j, zmm2=z_j
    #           zmm3=m_j*dt
    #           rsi=&hvx , rdx=&hvy, rcx=&hvz 
    #// 0123 4567
    #// 3201 7645
    
    vpermpd $0x4B, %zmm0, %zmm4               # 01234567 -> 32017645
    vpermpd $0x4B, %zmm1, %zmm5
    vpermpd $0x4B, %zmm2, %zmm6
    vpermpd $0x4B, %zmm3, %zmm15 

    vsubpd  %zmm4, %zmm0, %zmm7                # d_x
    vsubpd  %zmm5, %zmm1, %zmm8
    vsubpd  %zmm6, %zmm2, %zmm9
    

    # Prefactor calculation
    vmulpd      %zmm7, %zmm7, %zmm10     
    vfmadd231pd %zmm8, %zmm8, %zmm10      
    vfmadd231pd %zmm9, %zmm9, %zmm10     # zmm10 is now r^2
    
    vsqrtpd     %zmm10, %zmm11             
    vmulpd      %zmm10, %zmm11, %zmm10      # zmm10 is r^3
   
    vbroadcastsd .one(%rip), %zmm11         # Todo: keep 1 in a register at all times 
    vdivpd      %zmm10, %zmm11, %zmm17      # 1/r^3
    vmulpd      %zmm17, %zmm15, %zmm14      # m/r^3
    
    vmovapd	(%rdi),  %zmm10             # TODO get rid of mov instruction
	vmovapd	(%rsi),  %zmm11 
	vmovapd	(%rdx),  %zmm12
  
    vfnmadd231pd %zmm14, %zmm7,  %zmm10
    vfnmadd231pd %zmm14, %zmm8,  %zmm11
    vfnmadd231pd %zmm14, %zmm9,  %zmm12

	
    vmovapd	%zmm10, (%rdi)              # TODO get rid of mov instruction
	vmovapd	%zmm11, (%rsi)
	vmovapd	%zmm12, (%rdx)

    vpermpd $0x1E, %zmm7, %zmm7               # 32017645 -> 01234567
    vpermpd $0x1E, %zmm8, %zmm8
    vpermpd $0x1E, %zmm9, %zmm9
    vpermpd $0x1E, %zmm17, %zmm17
    vpermpd $0x1E, %zmm3, %zmm15                # 01234567 -> 32017645


    vmulpd      %zmm17, %zmm15, %zmm14      # m/r^3
	
    vmovapd	(%rdi),  %zmm10             # TODO get rid of mov instruction
	vmovapd	(%rsi),  %zmm11 
	vmovapd	(%rdx),  %zmm12 
    
    vfmadd231pd %zmm14, %zmm7,  %zmm10
    vfmadd231pd %zmm14, %zmm8,  %zmm11
    vfmadd231pd %zmm14, %zmm9,  %zmm12

	vmovapd	%zmm10, (%rdi)              # TODO get rid of mov instruction
	vmovapd	%zmm11, (%rsi)
	vmovapd	%zmm12, (%rdx)

    #// 0123 4567
    #// 2310 6754

    ret




.section .note.GNU-stack,"",@progbits
