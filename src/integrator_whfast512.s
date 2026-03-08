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
    #           zmm5 = mdt
    #           edi = mask
    #           rsi=&hvx , rdx=&hvy, rcx=&hvz 

	vmulpd	%zmm0, %zmm0, %zmm6
	vfmadd231pd	%zmm1, %zmm1, %zmm6
	vfmadd231pd	%zmm2, %zmm2, %zmm6     # zmm6 is r^2
	
    vsqrtpd	%zmm6, %zmm18
	vmulpd	%zmm6, %zmm18, %zmm18
	vdivpd	%zmm18, %zmm5, %zmm8{%k1}{z}  # zmm8 is  m0*dt/r^3

	vmulpd	%zmm6, %zmm6, %zmm7         # zmm7 is r^4
	vdivpd	%zmm7, %zmm3, %zmm9{%k1}{z}         

	vmulpd	%zmm9, %zmm0, %zmm15        # is x_j*gr_prefac/r^4 
	vmulpd	%zmm9, %zmm1, %zmm16
	vmulpd	%zmm9, %zmm2, %zmm17

    kmovw	%edi, %k1
	vmulpd	%zmm15, %zmm4, %zmm10{%k1}{z}
	vmulpd	%zmm16, %zmm4, %zmm11{%k1}{z}
	vmulpd	%zmm17, %zmm4, %zmm12{%k1}{z}

    REDUCE_ADD_AND_BROADCAST %zmm10, %zmm18   # sum_x
    REDUCE_ADD_AND_BROADCAST %zmm11, %zmm18
    REDUCE_ADD_AND_BROADCAST %zmm12, %zmm18

	vsubpd	%zmm10, %zmm15, %zmm10      # hvx = hvx - sum_x
	vsubpd	%zmm11, %zmm16, %zmm11
	vsubpd	%zmm12, %zmm17, %zmm12

    vbroadcastsd .t2(%rip), %zmm2
    vbroadcastsd .t3(%rip), %zmm3
    vbroadcastsd .t5(%rip), %zmm5
	vfnmadd213pd	%zmm2, %zmm3, %zmm5    # zmm5 = -(2*5 +3)
int3
    # 3                 
    # zmm0 = zmm10 - zmm8 * %zmm0 
	vfnmadd213pd	%zmm10, %zmm8, %zmm0    # hvx = hvx - x_j*prefact
	vmovapd	%zmm0, (%rsi)
	vfnmadd213pd	%zmm11, %zmm8, %zmm1
	vmovapd	%zmm1, (%rdx)
	vfnmadd213pd	%zmm12, %zmm8, %zmm2
	vmovapd	%zmm2, (%rcx)
	
    ret
