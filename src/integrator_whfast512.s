# file: integrator_whfast512.s

.section .rodata
    .align 64
.one: 
    .double 1.0


.section .text
.globl gravity_prefactor_avx512_one
.globl gravity_prefactor_avx512


gravity_prefactor_avx512_one:
    # Input:  zmm0=dx, zmm1=dy, zmm2=dy
    
    vmulpd      %zmm0, %zmm0, %zmm0     
    vfmadd231pd %zmm1, %zmm1, %zmm0      
    vfmadd231pd %zmm2, %zmm2, %zmm0      
    # zmm0 is now r^2
    
    vsqrtpd     %zmm0, %zmm1             
    vmulpd      %zmm0, %zmm1, %zmm0      # zmm0 is r^3
   
    vbroadcastsd .one(%rip), %zmm1       # Todo: keep 1 in a register at all times 
    vdivpd      %zmm0, %zmm1, %zmm0      
    
    # return 1 / r^3 in zmm0
    ret
    

gravity_prefactor_avx512:
    # Input:  zmm0=m, zmm1=dx, zmm2=dy, zmm3=dz
    
    vmulpd      %zmm1, %zmm1, %zmm1     
    vfmadd231pd %zmm2, %zmm2, %zmm1      
    vfmadd231pd %zmm3, %zmm3, %zmm1      
    # zmm1 is now r^2
    
    vsqrtpd     %zmm1, %zmm2             
    vmulpd      %zmm1, %zmm2, %zmm2      
    vdivpd      %zmm2, %zmm0, %zmm0      
    
    # return m / r^3 in zmm0
    ret

//static __m512d inline gravity_prefactor_avx512( __m512d m, __m512d dx, __m512d dy, __m512d dz) {
//    __m512d r2 = _mm512_mul_pd(dx, dx);
//    r2 = _mm512_fmadd_pd(dy,dy, r2);
//    r2 = _mm512_fmadd_pd(dz,dz, r2);
//    const __m512d r = _mm512_sqrt_pd(r2);
//    const __m512d r3 = _mm512_mul_pd(r, r2);
//    return _mm512_div_pd(m,r3);
//}
