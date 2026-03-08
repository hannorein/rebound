# file: add_func.s
.globl gravity_prefactor_avx512

gravity_prefactor_avx512:
    movq %rdi, %rax       # Move the first argument (n1) into the return register (rax)
    addq %rsi, %rax       # Add the second argument (n2) to rax
    ret                   # Return to the caller

//static __m512d inline gravity_prefactor_avx512( __m512d m, __m512d dx, __m512d dy, __m512d dz) {
//    __m512d r2 = _mm512_mul_pd(dx, dx);
//    r2 = _mm512_fmadd_pd(dy,dy, r2);
//    r2 = _mm512_fmadd_pd(dz,dz, r2);
//    const __m512d r = _mm512_sqrt_pd(r2);
//    const __m512d r3 = _mm512_mul_pd(r, r2);
//    return _mm512_div_pd(m,r3);
//}
