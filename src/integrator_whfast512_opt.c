/**
variants of WHFast512 functions for benchmarking
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rebound.h"
#include "integrator_whfast512.h"

#ifdef AVX512
#include <immintrin.h>

// newton with reciprocal approximation
// replace 1/fp division with rcp14 + newton


static inline __m512d fast_recip_pd(__m512d x) {
    __m512d t = _mm512_rcp14_pd(x);
    // One Newton-Raphson: t = t * (2 - x*t)
    t = _mm512_mul_pd(t, _mm512_fnmadd_pd(x, t, _mm512_set1_pd(2.0)));
    return t;
}

// same as above but with two newton iterations
static inline __m512d fast_recip_pd_2iter(__m512d x) {
    __m512d t = _mm512_rcp14_pd(x);
    // first newton iteration
    t = _mm512_mul_pd(t, _mm512_fnmadd_pd(x, t, _mm512_set1_pd(2.0)));
    // second newton iteration
    t = _mm512_mul_pd(t, _mm512_fnmadd_pd(x, t, _mm512_set1_pd(2.0)));
    return t;
}

// fast inverse square root for gravity
// replace sqrt + div with rsqrt14 + newton


static inline __m512d fast_gravity_prefactor(__m512d r2) {
    __m512d rsqrt = _mm512_rsqrt14_pd(r2);
    
    // newton for: y = y * (3 - x*y^2) / 2
    __m512d y2 = _mm512_mul_pd(rsqrt, rsqrt);
    __m512d three = _mm512_set1_pd(3.0);
    __m512d half = _mm512_set1_pd(0.5);
    rsqrt = _mm512_mul_pd(rsqrt, _mm512_mul_pd(half, _mm512_fnmadd_pd(r2, y2, three)));
    
    // second newton iteration
    y2 = _mm512_mul_pd(rsqrt, rsqrt);
    rsqrt = _mm512_mul_pd(rsqrt, _mm512_mul_pd(half, _mm512_fnmadd_pd(r2, y2, three)));
    
    // compute 1/r^3 = (1/sqrt(r2))^3 = rsqrt * rsqrt * rsqrt
    __m512d rsqrt2 = _mm512_mul_pd(rsqrt, rsqrt);
    return _mm512_mul_pd(rsqrt2, rsqrt);
}

// with mass multiplication: m / r^3
static inline __m512d fast_gravity_prefactor_with_mass(__m512d m, __m512d dx, __m512d dy, __m512d dz) {
    __m512d r2 = _mm512_mul_pd(dx, dx);
    r2 = _mm512_fmadd_pd(dy, dy, r2);
    r2 = _mm512_fmadd_pd(dz, dz, r2);
    
    __m512d inv_r3 = fast_gravity_prefactor(r2);
    return _mm512_mul_pd(m, inv_r3);
}

// simplely just compute 1/r^3 for a given position
static inline __m512d fast_gravity_prefactor_one(__m512d dx, __m512d dy, __m512d dz) {
    __m512d r2 = _mm512_mul_pd(dx, dx);
    r2 = _mm512_fmadd_pd(dy, dy, r2);
    r2 = _mm512_fmadd_pd(dz, dz, r2);
    return fast_gravity_prefactor(r2);
}

#endif
