# file: add_func.s
.section .text
.globl gravity_prefactor_avx512


gravity_prefactor_avx512:
    # Input:  zmm0=m, zmm1=dx, zmm2=dy, zmm3=dz
    
    # 1. r2 = dx*dx
    vmulpd      %zmm1, %zmm1, %zmm1     
    
    # 2. r2 = (dy*dy) + r2
    # vfmadd231pd: dest = (src1 * src2) + dest
    vfmadd231pd %zmm2, %zmm2, %zmm1      
    
    # 3. r2 = (dz*dz) + r2
    vfmadd231pd %zmm3, %zmm3, %zmm1      
    
    # 4. r = sqrt(r2)
    vsqrtpd     %zmm1, %zmm2             
    
    # 5. r3 = r * r2
    vmulpd      %zmm1, %zmm2, %zmm2      
    
    # 6. return m / r3
    # Result must be in zmm0
    vdivpd      %zmm2, %zmm0, %zmm0      

    ret

//static __m512d inline gravity_prefactor_avx512( __m512d m, __m512d dx, __m512d dy, __m512d dz) {
//    __m512d r2 = _mm512_mul_pd(dx, dx);
//    r2 = _mm512_fmadd_pd(dy,dy, r2);
//    r2 = _mm512_fmadd_pd(dz,dz, r2);
//    const __m512d r = _mm512_sqrt_pd(r2);
//    const __m512d r3 = _mm512_mul_pd(r, r2);
//    return _mm512_div_pd(m,r3);
//}
