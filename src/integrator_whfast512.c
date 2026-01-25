/**
 * @file    integrator_whfast.c
 * @brief   WHFAST512 integration scheme.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 *          Pejvak Javaheri <pejvak.javaheri@mail.utoronto.ca>
 * @details This file implements the WHFast512 integration scheme with 
 * optimizations for AVX512.
 * 
 * @section LICENSE
 * Copyright (c) 2023 Hanno Rein, Pejvak Javaheri
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "rebound.h"
#include "particle.h"
#include "tools.h"
#include "gravity.h"
#include "boundary.h"
#include "integrator.h"
#include "integrator_whfast.h"
#include "integrator_whfast512.h"

#ifdef AVX512
// Debug function to print vectors
void static inline printavx512(__m512d a) {
    double _nax[8];
    _mm512_storeu_pd(&_nax[0], a);
    printf("avx\n %.15e\n %.15e\n %.15e\n %.15e\n %.15e\n %.15e\n %.15e\n %.15e\n", _nax[0], _nax[1], _nax[2], _nax[3], _nax[4], _nax[5], _nax[6], _nax[7]);
}

// Print mask in binary format for debuggin
void static inline printmask8(__mmask8 mask) {
    for (int i = 7; i >= 0; i--) {
        printf("%d", (mask >> i) & 1);
        if (i == 4) printf(" "); 
    }
    printf("\n");
}
#endif

// Debug function to print 8x8 matrix
void static inline printmat8(double* a) {
    for (int i=0; i<8;i++){
        for (int j=0; j<8;j++){
            printf("%.16f ", a[i*8+j]);
        }
        printf("\n");
    }
}

#ifdef PROF
// Profiling counters
double walltime_interaction=0;;
double walltime_kepler=0;
double walltime_jump=0;
double walltime_com=0;
#endif

#ifdef AVX512
// 8x8 matrix multiplication using avx512
static inline __m512d mat8_mul_avx512(const double* matrix, const __m512d vector) {
    __m512d v_i = _mm512_set1_pd(vector[0]);
    __m512d col_i = _mm512_load_pd(matrix);
    __m512d res = _mm512_mul_pd(v_i, col_i);
    for (int i = 1; i < 8; i++) {
        __m512d v_i = _mm512_set1_pd(vector[i]);
        __m512d col_i = _mm512_load_pd(&matrix[i * 8]);
        res = _mm512_fmadd_pd(v_i, col_i, res);
    }
    return res;
}

// Three 8x8 matrix multiplications with the same matrix using avx512
// Used for coordinate transformations in x, y, and z
static inline void mat8_mul3_avx512(const double* matrix, const __m512d in1, const __m512d in2, const __m512d in3, __m512d* out1, __m512d* out2, __m512d* out3) {
    __m512d col_i = _mm512_load_pd(matrix);
    __m512d vin1 = _mm512_set1_pd(in1[0]);
    *out1 = _mm512_mul_pd(vin1, col_i);
    __m512d vin2 = _mm512_set1_pd(in2[0]);
    *out2 = _mm512_mul_pd(vin2, col_i);
    __m512d vin3 = _mm512_set1_pd(in3[0]);
    *out3 = _mm512_mul_pd(vin3, col_i);
    for (int i = 1; i < 8; i++) {
        __m512d col_i = _mm512_load_pd(&matrix[i * 8]);
        __m512d vin1 = _mm512_set1_pd(in1[i]);
        *out1 = _mm512_fmadd_pd(vin1, col_i, *out1);
        __m512d vin2 = _mm512_set1_pd(in2[i]);
        *out2 = _mm512_fmadd_pd(vin2, col_i, *out2);
        __m512d vin3 = _mm512_set1_pd(in3[i]);
        *out3 = _mm512_fmadd_pd(vin3, col_i, *out3);
    }
}
   
// Hepler function to load particle data into avx512 registers
static __m512d load_into_m512d(struct reb_simulation* r, size_t offset, const double* transformation, int N){
    struct reb_integrator_whfast512* const ri_whfast512 = &(r->ri_whfast512);
    struct reb_particle* particles = r->particles;
    const unsigned int N_systems = ri_whfast512->N_systems;
    const unsigned int p_per_system = 8/N_systems;
    const unsigned int N_per_system = r->N/N_systems;
    double tmp[8] = {0}; 
    for (int s=0; s<N_systems; s++){
        for (unsigned int i=1; i<N_per_system; i++){
            tmp[s*p_per_system+i-1] = *(double*)((char*)(&particles[s*N_per_system+i])+offset);
        }
    }
    __m512d tmp512 = _mm512_loadu_pd(tmp);
    if (transformation != NULL){
        return mat8_mul_avx512(transformation, tmp512);
    }else{
        return tmp512;
    }
}

// Hepler function to load particle data from avx512 registers
static void load_from_m512d(struct reb_simulation* r, size_t offset, const double* transformation, int N, __m512d vector){
    struct reb_integrator_whfast512* const ri_whfast512 = &(r->ri_whfast512);
    struct reb_particle* particles = r->particles;
    const unsigned int N_systems = ri_whfast512->N_systems;
    const unsigned int p_per_system = 8/N_systems;
    const unsigned int N_per_system = r->N/N_systems;
    double tmp[8] __attribute__((aligned(64)));; 
    __m512d tmp512;
    if (transformation != NULL){
        tmp512 = mat8_mul_avx512(transformation, vector);
    }else{
        tmp512 = vector;
    }
    _mm512_store_pd(tmp, tmp512);
    for (int s=0; s<N_systems; s++){
        for (unsigned int i=1; i<N_per_system; i++){
            *(double*)((char*)(&particles[s*N_per_system+i])+offset) = tmp[s*p_per_system+i-1]; 
        }
    }
}
#endif

// Performs one full center of mass step (H_0)
static void reb_whfast512_com_step(struct reb_simulation* r, const double _dt){
#ifdef PROF
    struct reb_timeval time_beginning;
    gettimeofday(&time_beginning,NULL);
#endif
    const unsigned int N_systems = r->ri_whfast512.N_systems;
    for (int s=0; s<N_systems; s++){
        r->ri_whfast512.p_jh0[s].x += _dt*r->ri_whfast512.p_jh0[s].vx;
        r->ri_whfast512.p_jh0[s].y += _dt*r->ri_whfast512.p_jh0[s].vy;
        r->ri_whfast512.p_jh0[s].z += _dt*r->ri_whfast512.p_jh0[s].vz;
    }
#ifdef PROF
    struct reb_timeval time_end;
    gettimeofday(&time_end,NULL);
    walltime_com += time_end.tv_sec-time_beginning.tv_sec+(time_end.tv_usec-time_beginning.tv_usec)/1e6;
#endif
}

// Convert democratic heliocentric coordinates to inertial coordinates
// Note: this is only called at the end. Speed is not a concern.
static void democraticheliocentric_to_inertial_posvel(struct reb_simulation* r){
    struct reb_integrator_whfast512* const ri_whfast512 = &(r->ri_whfast512);
    struct reb_particle* particles = r->particles;
    struct reb_particle_avx512* p512 = ri_whfast512->p512;
    const unsigned int N_systems = ri_whfast512->N_systems;
    const unsigned int p_per_system = 8/N_systems;
    const unsigned int N_per_system = r->N/N_systems;


#ifdef AVX512
    double m[8];
    double x[8];
    double y[8];
    double z[8];
    double vx[8];
    double vy[8];
    double vz[8];
    _mm512_storeu_pd(&m, p512->m);
    _mm512_storeu_pd(&x, p512->x);
    _mm512_storeu_pd(&y, p512->y);
    _mm512_storeu_pd(&z, p512->z);
    _mm512_storeu_pd(&vx, p512->vx);
    _mm512_storeu_pd(&vy, p512->vy);
    _mm512_storeu_pd(&vz, p512->vz);
#else // AVX512
      // Fallback for synchronization for when AVX512 is not available
    double* m = (double*)p512->m;
    double* x = (double*)p512->x;
    double* y = (double*)p512->y;
    double* z = (double*)p512->z;
    double* vx = p512->vx;
    double* vy = p512->vy;
    double* vz = p512->vz;
#endif // AVX512

    for (unsigned s=0;s<N_systems;s++){
        double x0s = 0;
        double y0s = 0;
        double z0s = 0;
        double vx0s = 0;
        double vy0s = 0;
        double vz0s = 0;
        for (unsigned int i=1; i<N_per_system; i++){
            x0s += x[s*p_per_system+(i-1)] * m[s*p_per_system+(i-1)];
            y0s += y[s*p_per_system+(i-1)] * m[s*p_per_system+(i-1)];
            z0s += z[s*p_per_system+(i-1)] * m[s*p_per_system+(i-1)];
            vx0s += vx[s*p_per_system+(i-1)] * m[s*p_per_system+(i-1)];
            vy0s += vy[s*p_per_system+(i-1)] * m[s*p_per_system+(i-1)];
            vz0s += vz[s*p_per_system+(i-1)] * m[s*p_per_system+(i-1)];
            particles[s*N_per_system+i].vx = vx[s*p_per_system+(i-1)] + ri_whfast512->p_jh0[s].vx;
            particles[s*N_per_system+i].vy = vy[s*p_per_system+(i-1)] + ri_whfast512->p_jh0[s].vy;
            particles[s*N_per_system+i].vz = vz[s*p_per_system+(i-1)] + ri_whfast512->p_jh0[s].vz;
        }
        x0s /= ri_whfast512->p_jh0[s].m;
        y0s /= ri_whfast512->p_jh0[s].m;
        z0s /= ri_whfast512->p_jh0[s].m;
        particles[s*N_per_system].x  = ri_whfast512->p_jh0[s].x - x0s;
        particles[s*N_per_system].y  = ri_whfast512->p_jh0[s].y - y0s;
        particles[s*N_per_system].z  = ri_whfast512->p_jh0[s].z - z0s;
        particles[s*N_per_system].vx = ri_whfast512->p_jh0[s].vx - vx0s;
        particles[s*N_per_system].vy = ri_whfast512->p_jh0[s].vy - vy0s;
        particles[s*N_per_system].vz = ri_whfast512->p_jh0[s].vz - vz0s;
        for (unsigned int i=1; i<N_per_system; i++){
            particles[s*N_per_system+i].x  = x[s*p_per_system+(i-1)] + particles[s*N_per_system].x;
            particles[s*N_per_system+i].y  = y[s*p_per_system+(i-1)] + particles[s*N_per_system].y;
            particles[s*N_per_system+i].z  = z[s*p_per_system+(i-1)] + particles[s*N_per_system].z;
        }
    }
}
                
// Speed is not a concern here 
static void jacobi_to_inertial_posvel(struct reb_simulation* r){
    struct reb_integrator_whfast512* const ri_whfast512 = &(r->ri_whfast512);
    const unsigned int N_systems = ri_whfast512->N_systems;
    const unsigned int N_per_system = r->N/N_systems;
#ifdef AVX512
    struct reb_particle* particles = r->particles;
    load_from_m512d(r, offsetof(struct reb_particle, x), ri_whfast512->mat8_jacobi_to_inertial, r->N, ri_whfast512->p512->x);
    load_from_m512d(r, offsetof(struct reb_particle, y), ri_whfast512->mat8_jacobi_to_inertial, r->N, ri_whfast512->p512->y);
    load_from_m512d(r, offsetof(struct reb_particle, z), ri_whfast512->mat8_jacobi_to_inertial, r->N, ri_whfast512->p512->z);
    load_from_m512d(r, offsetof(struct reb_particle, vx), ri_whfast512->mat8_jacobi_to_inertial, r->N, ri_whfast512->p512->vx);
    load_from_m512d(r, offsetof(struct reb_particle, vy), ri_whfast512->mat8_jacobi_to_inertial, r->N, ri_whfast512->p512->vy);
    load_from_m512d(r, offsetof(struct reb_particle, vz), ri_whfast512->mat8_jacobi_to_inertial, r->N, ri_whfast512->p512->vz);
    for (unsigned s=0;s<N_systems;s++){
        particles[s*N_per_system+0].x  = 0.0;
        particles[s*N_per_system+0].y  = 0.0;
        particles[s*N_per_system+0].z  = 0.0;
        particles[s*N_per_system+0].vx = 0.0;
        particles[s*N_per_system+0].vy = 0.0;
        particles[s*N_per_system+0].vz = 0.0;
        for (int i=1;i<N_per_system;i++){
            particles[s*N_per_system+0].x  -= particles[s*N_per_system+i].m/particles[s*N_per_system+0].m*particles[s*N_per_system+i].x;
            particles[s*N_per_system+0].y  -= particles[s*N_per_system+i].m/particles[s*N_per_system+0].m*particles[s*N_per_system+i].y;
            particles[s*N_per_system+0].z  -= particles[s*N_per_system+i].m/particles[s*N_per_system+0].m*particles[s*N_per_system+i].z;
            particles[s*N_per_system+0].vx -= particles[s*N_per_system+i].m/particles[s*N_per_system+0].m*particles[s*N_per_system+i].vx;
            particles[s*N_per_system+0].vy -= particles[s*N_per_system+i].m/particles[s*N_per_system+0].m*particles[s*N_per_system+i].vy;
            particles[s*N_per_system+0].vz -= particles[s*N_per_system+i].m/particles[s*N_per_system+0].m*particles[s*N_per_system+i].vz;
        }
        for (int i=0;i<N_per_system;i++){
            particles[s*N_per_system+i].x  += ri_whfast512->p_jh0[s].x;
            particles[s*N_per_system+i].y  += ri_whfast512->p_jh0[s].y;
            particles[s*N_per_system+i].z  += ri_whfast512->p_jh0[s].z;
            particles[s*N_per_system+i].vx += ri_whfast512->p_jh0[s].vx;
            particles[s*N_per_system+i].vy += ri_whfast512->p_jh0[s].vy;
            particles[s*N_per_system+i].vz += ri_whfast512->p_jh0[s].vz;
        }
    }
#else  // AVX512
    reb_simulation_error(r, "Fallback for Jacobi transformations in whfast512 not yet implemented.");
#endif // AVX512
}

#ifdef AVX512

// Fast inverse factorial lookup table
static const double invfactorial[35] = {1., 1., 1./2., 1./6., 1./24., 1./120., 1./720., 1./5040., 1./40320., 1./362880., 1./3628800., 1./39916800., 1./479001600., 1./6227020800., 1./87178291200., 1./1307674368000., 1./20922789888000., 1./355687428096000., 1./6402373705728000., 1./121645100408832000., 1./2432902008176640000., 1./51090942171709440000., 1./1124000727777607680000., 1./25852016738884976640000., 1./620448401733239439360000., 1./15511210043330985984000000., 1./403291461126605635584000000., 1./10888869450418352160768000000., 1./304888344611713860501504000000., 1./8841761993739701954543616000000., 1./265252859812191058636308480000000., 1./8222838654177922817725562880000000., 1./263130836933693530167218012160000000., 1./8683317618811886495518194401280000000., 1./295232799039604140847618609643520000000.};


// Stiefel function for Newton's method, returning Gs1, Gs2, and Gs3
static void inline mm_stiefel_Gs13_avx512(__m512d * Gs1, __m512d * Gs2, __m512d * Gs3, __m512d beta, __m512d X){
    __m512d X2 = _mm512_mul_pd(X,X); 
    __m512d z = _mm512_mul_pd(X2,beta); 

    // stumpff_cs. Note: assuming n = 0
    const int nmax = 19;
    *Gs3 = _mm512_set1_pd(invfactorial[nmax]); 
    *Gs2 = _mm512_set1_pd(invfactorial[nmax-1]); 

    for(int np=nmax-2;np>=3;np-=2){
        *Gs3 = _mm512_fnmadd_pd(z, *Gs3, _mm512_set1_pd(invfactorial[np]));
        *Gs2 = _mm512_fnmadd_pd(z, *Gs2, _mm512_set1_pd(invfactorial[np-1]));
    }
    *Gs3 = _mm512_mul_pd(*Gs3,X); 
    *Gs1 = _mm512_fnmadd_pd(z, *Gs3, X);
    *Gs3 = _mm512_mul_pd(*Gs3,X2); 
    *Gs2 = _mm512_mul_pd(*Gs2,X2); 
};

// Stiefel function for Halley's method, returning Gs0, Gs1, Gs2, and Gs3
static void inline mm_stiefel_Gs03_avx512(__m512d * Gs0, __m512d * Gs1, __m512d * Gs2, __m512d * Gs3, __m512d beta, __m512d X){
    __m512d X2 = _mm512_mul_pd(X,X); 
    __m512d z = _mm512_mul_pd(X2,beta); 

    // stumpff_cs. Note: assuming n = 0
    const int nmax = 11; // Note: reduced! needs to be improved with mm_stiefel_Gs13_avx512 on last step(s)
    *Gs3 = _mm512_set1_pd(invfactorial[nmax]); 
    *Gs2 = _mm512_set1_pd(invfactorial[nmax-1]); 

    for(int np=nmax-2;np>=3;np-=2){
        *Gs3 = _mm512_fnmadd_pd(z, *Gs3, _mm512_set1_pd(invfactorial[np]));
        *Gs2 = _mm512_fnmadd_pd(z, *Gs2, _mm512_set1_pd(invfactorial[np-1]));
    }
    *Gs0 = _mm512_fnmadd_pd(z, *Gs2, _mm512_set1_pd(1.0));
    *Gs3 = _mm512_mul_pd(*Gs3,X); 
    *Gs1 = _mm512_fnmadd_pd(z, *Gs3, X);
    *Gs3 = _mm512_mul_pd(*Gs3,X2); 
    *Gs2 = _mm512_mul_pd(*Gs2,X2); 
};

// Performs one full Kepler step
static void inline reb_whfast512_kepler_step(const struct reb_simulation* const r, const double dt){
#ifdef PROF
    struct reb_timeval time_beginning;
    gettimeofday(&time_beginning,NULL);
#endif
    struct reb_particle_avx512 * restrict p512  = r->ri_whfast512.p512;
    __m512d _dt = _mm512_set1_pd(dt); 

    __m512d r2 = _mm512_mul_pd(p512->x, p512->x);
    r2 = _mm512_fmadd_pd(p512->y, p512->y, r2);
    r2 = _mm512_fmadd_pd(p512->z, p512->z, r2);
    __m512d r0 = _mm512_sqrt_pd(r2);
    __m512d r0i = _mm512_div_pd(_mm512_set1_pd(1.0),r0);

    __m512d v2 = _mm512_mul_pd(p512->vx, p512->vx);
    v2 = _mm512_fmadd_pd(p512->vy, p512->vy, v2);
    v2 = _mm512_fmadd_pd(p512->vz, p512->vz, v2);

    __m512d beta = _mm512_mul_pd(_mm512_set1_pd(2.0), p512->M);
    beta = _mm512_fmsub_pd(beta, r0i, v2);

    __m512d eta0 = _mm512_mul_pd(p512->x, p512->vx);
    eta0 = _mm512_fmadd_pd(p512->y, p512->vy, eta0);
    eta0 = _mm512_fmadd_pd(p512->z, p512->vz, eta0);

    __m512d zeta0 = _mm512_fnmadd_pd(beta, r0, p512->M);

    __m512d Gs1;
    __m512d Gs2;
    __m512d Gs3;
    __m512d eta0Gs1zeta0Gs2; 
    __m512d ri; 

#define NEWTON_STEP() \
    mm_stiefel_Gs13_avx512(&Gs1, &Gs2, &Gs3, beta, X);\
    eta0Gs1zeta0Gs2 = _mm512_mul_pd(eta0, Gs1); \
    eta0Gs1zeta0Gs2 = _mm512_fmadd_pd(zeta0,Gs2, eta0Gs1zeta0Gs2); \
    ri = _mm512_add_pd(r0, eta0Gs1zeta0Gs2); \
    ri = _mm512_div_pd(_mm512_set1_pd(1.0), ri); \
    \
    X = _mm512_mul_pd(X, eta0Gs1zeta0Gs2);\
    X = _mm512_fnmadd_pd(eta0, Gs2, X);\
    X = _mm512_fnmadd_pd(zeta0, Gs3, X);\
    X = _mm512_add_pd(_dt, X);\
    X = _mm512_mul_pd(ri, X);


#define HALLEY_STEP() \
    mm_stiefel_Gs03_avx512(&Gs0, &Gs1, &Gs2, &Gs3, beta, X);\
    f = _mm512_fmsub_pd(r0,X,_dt);\
    f = _mm512_fmadd_pd(eta0, Gs2, f);\
    f = _mm512_fmadd_pd(zeta0, Gs3, f);\
    \
    fp = _mm512_fmadd_pd(eta0, Gs1, r0);\
    fp = _mm512_fmadd_pd(zeta0, Gs2, fp);\
    \
    fpp = _mm512_mul_pd(eta0, Gs0);\
    fpp = _mm512_fmadd_pd(zeta0, Gs1, fpp);\
    \
    denom = _mm512_mul_pd(fp,fp);\
    denom = _mm512_mul_pd(denom,_mm512_set1_pd(16.0));\
    \
    denom = _mm512_fnmadd_pd(_mm512_mul_pd(f,fpp),_mm512_set1_pd(20.0), denom);\
    /* not included: _mm512_abs_pd(denom) */;\
    denom = _mm512_sqrt_pd(denom);\
    denom = _mm512_add_pd(fp, denom);\
    \
    X = _mm512_fmsub_pd(X, denom, _mm512_mul_pd(f, _mm512_set1_pd(5.0)));\
    X = _mm512_div_pd(X, denom);
    
    // Division approximation by newton's method. Not much faster
    // __m512d t; 
    //t = _mm512_rcp14_pd(denom);
    //t = _mm512_mul_pd(t, _mm512_sub_pd(_mm512_set1_pd(2.0),_mm512_mul_pd(denom,t)));
    //X = _mm512_mul_pd(X, t);
    


    // Initial guess
    __m512d dtr0i = _mm512_mul_pd(_dt,r0i);
    __m512d X = _mm512_mul_pd(dtr0i,eta0);
    X = _mm512_mul_pd(X,_mm512_set1_pd(0.5));
    X = _mm512_fnmadd_pd(X,r0i,_mm512_set1_pd(1.0));
    X = _mm512_mul_pd(dtr0i,X);

    // Iterations
    __m512d f, fp, fpp, denom, Gs0;
    HALLEY_STEP();
    HALLEY_STEP();
    NEWTON_STEP();
    // +1 below

    // Final Newton step (note: X not needed after this) 
    mm_stiefel_Gs13_avx512(&Gs1, &Gs2, &Gs3, beta, X);
    eta0Gs1zeta0Gs2 = _mm512_mul_pd(eta0, Gs1); 
    eta0Gs1zeta0Gs2 = _mm512_fmadd_pd(zeta0,Gs2, eta0Gs1zeta0Gs2); 
    ri = _mm512_add_pd(r0, eta0Gs1zeta0Gs2); 
    ri = _mm512_div_pd(_mm512_set1_pd(1.0), ri); 

    // f and g function

    __m512d nf = _mm512_mul_pd(p512->M,Gs2); //negative f
    nf = _mm512_mul_pd(nf,r0i); 

    __m512d g = _mm512_fnmadd_pd(p512->M, Gs3, _dt);

    __m512d nfd = _mm512_mul_pd(p512->M, Gs1); // negative fd
    nfd = _mm512_mul_pd(nfd, r0i);
    nfd = _mm512_mul_pd(nfd, ri);

    __m512d ngd = _mm512_mul_pd(p512->M, Gs2); // negative gd
    ngd = _mm512_mul_pd(ngd, ri);

    __m512d nx = _mm512_maskz_fnmadd_pd(p512->mask, nf, p512->x, p512->x);
    nx = _mm512_maskz_fmadd_pd(p512->mask, g, p512->vx, nx);
    __m512d ny = _mm512_maskz_fnmadd_pd(p512->mask, nf, p512->y, p512->y);
    ny = _mm512_maskz_fmadd_pd(p512->mask, g, p512->vy, ny);
    __m512d nz = _mm512_maskz_fnmadd_pd(p512->mask, nf, p512->z, p512->z);
    nz = _mm512_maskz_fmadd_pd(p512->mask, g, p512->vz, nz);

    p512->vx = _mm512_maskz_fnmadd_pd(p512->mask, ngd, p512->vx, p512->vx);
    p512->vx = _mm512_maskz_fnmadd_pd(p512->mask, nfd, p512->x, p512->vx);
    p512->vy = _mm512_maskz_fnmadd_pd(p512->mask, ngd, p512->vy, p512->vy);
    p512->vy = _mm512_maskz_fnmadd_pd(p512->mask, nfd, p512->y, p512->vy);
    p512->vz = _mm512_maskz_fnmadd_pd(p512->mask, ngd, p512->vz, p512->vz);
    p512->vz = _mm512_maskz_fnmadd_pd(p512->mask, nfd, p512->z, p512->vz);
    

    p512->x = nx;
    p512->y = ny;
    p512->z = nz;
#ifdef PROF
    struct reb_timeval time_end;
    gettimeofday(&time_end,NULL);
    walltime_kepler += time_end.tv_sec-time_beginning.tv_sec+(time_end.tv_usec-time_beginning.tv_usec)/1e6;
#endif
}

// Helper functions for the interaction step
// Calculates 1/(dx**2+dy**2+dz**2)^(3/2)
static __m512d inline gravity_prefactor_avx512_one( __m512d dx, __m512d dy, __m512d dz) {
    __m512d r2 = _mm512_mul_pd(dx, dx);
    r2 = _mm512_fmadd_pd(dy,dy, r2);
    r2 = _mm512_fmadd_pd(dz,dz, r2);
    const __m512d r = _mm512_sqrt_pd(r2);
    const __m512d r3 = _mm512_mul_pd(r, r2);
    return _mm512_div_pd(_mm512_set1_pd(1.0),r3); 
}

static __m512d inline gravity_prefactor_avx512( __m512d m, __m512d dx, __m512d dy, __m512d dz) {
    __m512d r2 = _mm512_mul_pd(dx, dx);
    r2 = _mm512_fmadd_pd(dy,dy, r2);
    r2 = _mm512_fmadd_pd(dz,dz, r2);
    const __m512d r = _mm512_sqrt_pd(r2);
    const __m512d r3 = _mm512_mul_pd(r, r2);
    return _mm512_div_pd(m,r3);
}



// Performs one full interaction step
static void reb_whfast512_interaction_step_8planets_jacobi(struct reb_simulation * r){
#ifdef PROF
    struct reb_timeval time_beginning;
    gettimeofday(&time_beginning,NULL);
#endif
    struct reb_integrator_whfast512* const ri_whfast512 = &(r->ri_whfast512);
    struct reb_particle_avx512* restrict p512 = ri_whfast512->p512;

    // Convert position from Jacobi to heliocentric
    mat8_mul3_avx512(ri_whfast512->mat8_jacobi_to_heliocentric,
            p512->x, p512->y, p512->z,
            &p512->hx, &p512->hy, &p512->hz);
    __m512d x_j =  p512->hx;
    __m512d y_j =  p512->hy;
    __m512d z_j =  p512->hz;
    __m512d dt512 = _mm512_set1_pd(r->dt); 
    
    // General relativistic corrections
    if (ri_whfast512->gr_potential){
        __m512d r2 = _mm512_mul_pd(x_j, x_j);
        r2 = _mm512_fmadd_pd(y_j, y_j, r2);
        r2 = _mm512_fmadd_pd(z_j, z_j, r2);
        const __m512d r4 = _mm512_mul_pd(r2, r2);
        __m512d prefac = _mm512_div_pd(ri_whfast512->p512->gr_prefac,r4);
        __m512d dvx = _mm512_mul_pd(prefac, x_j); 
        __m512d dvy = _mm512_mul_pd(prefac, y_j); 
        __m512d dvz = _mm512_mul_pd(prefac, z_j); 
        p512->hvx  = dvx;
        p512->hvy  = dvy;
        p512->hvz  = dvz;

        // Calculate back reaction onto star and apply them to planets (heliocentric) 
        dvx = _mm512_maskz_mul_pd(p512->mask, ri_whfast512->p512->gr_prefac2, dvx); 
        dvy = _mm512_maskz_mul_pd(p512->mask, ri_whfast512->p512->gr_prefac2, dvy); 
        dvz = _mm512_maskz_mul_pd(p512->mask, ri_whfast512->p512->gr_prefac2, dvz); 

        double sum_x = _mm512_reduce_add_pd(dvx);
        double sum_y = _mm512_reduce_add_pd(dvy);
        double sum_z = _mm512_reduce_add_pd(dvz);
        p512->hvx = _mm512_maskz_sub_pd(p512->mask, p512->hvx, _mm512_set1_pd(sum_x));
        p512->hvy = _mm512_maskz_sub_pd(p512->mask, p512->hvy, _mm512_set1_pd(sum_y));
        p512->hvz = _mm512_maskz_sub_pd(p512->mask, p512->hvz, _mm512_set1_pd(sum_z));
    
        // Jacobi additions:
        // Need to add stellar term. Easy: already in heliocentric coordinates.
        // TODO: Should put a mask on particle 1 as +/- cancels.
        const __m512d r1 = _mm512_sqrt_pd(r2);
        const __m512d r3 = _mm512_mul_pd(r1, r2);
        __m512d prefact = _mm512_div_pd(_mm512_set1_pd(r->dt*r->particles[0].m),r3); // Note: using m0, m0, m0, ... for mass here
        p512->hvx = _mm512_maskz_fnmadd_pd(p512->mask, prefact, x_j, p512->hvx); 
        p512->hvy = _mm512_maskz_fnmadd_pd(p512->mask, prefact, y_j, p512->hvy); 
        p512->hvz = _mm512_maskz_fnmadd_pd(p512->mask, prefact, z_j, p512->hvz); 
    }else{
        // Jacobi additions:
        // TODO: Should put a mask on particle 1 as +/- cancels.
        // Need to add stellar term. Easy: already in heliocentric coordinates.
        __m512d prefact = gravity_prefactor_avx512_one(x_j, y_j, z_j);
        __m512d prefact1 = _mm512_mul_pd(prefact, _mm512_set1_pd(-r->particles[0].m)); // Note: using m0, m0, m0, ...
        prefact1 = _mm512_mul_pd(prefact1, dt512);
        p512->hvx = _mm512_maskz_mul_pd(p512->mask, prefact1, x_j); 
        p512->hvy = _mm512_maskz_mul_pd(p512->mask, prefact1, y_j); 
        p512->hvz = _mm512_maskz_mul_pd(p512->mask, prefact1, z_j); 
    }

    __m512d m_j = _mm512_mul_pd(p512->m, dt512);
    __m512d m_j_01234567 = m_j;

    {
        x_j = _mm512_permutex_pd(x_j, _MM_PERM_BACD); // within 256
        y_j = _mm512_permutex_pd(y_j, _MM_PERM_BACD);
        z_j = _mm512_permutex_pd(z_j, _MM_PERM_BACD);
        m_j = _mm512_permutex_pd(m_j, _MM_PERM_BACD);
        __m512d dx_j = _mm512_sub_pd(p512->hx, x_j);
        __m512d dy_j = _mm512_sub_pd(p512->hy, y_j);
        __m512d dz_j = _mm512_sub_pd(p512->hz, z_j);
        __m512d prefact = gravity_prefactor_avx512_one(dx_j, dy_j, dz_j);

        // 0123 4567
        // 3201 7645
        __m512d prefact1 = _mm512_mul_pd(prefact, m_j);
        p512->hvx = _mm512_fnmadd_pd(prefact1, dx_j, p512->hvx); 
        p512->hvy = _mm512_fnmadd_pd(prefact1, dy_j, p512->hvy); 
        p512->hvz = _mm512_fnmadd_pd(prefact1, dz_j, p512->hvz); 


        dx_j    = _mm512_permutex_pd(dx_j,    _MM_PERM_ABDC); // within 256
        dy_j    = _mm512_permutex_pd(dy_j,    _MM_PERM_ABDC);
        dz_j    = _mm512_permutex_pd(dz_j,    _MM_PERM_ABDC);
        prefact = _mm512_permutex_pd(prefact, _MM_PERM_ABDC);
        m_j     = _mm512_permute_pd(m_j,      0x55);    // within 128

        // 0123 4567
        // 2310 6754
        __m512d prefact2 = _mm512_mul_pd(prefact, m_j);
        p512->hvx = _mm512_fmadd_pd(prefact2, dx_j, p512->hvx); 
        p512->hvy = _mm512_fmadd_pd(prefact2, dy_j, p512->hvy); 
        p512->hvz = _mm512_fmadd_pd(prefact2, dz_j, p512->hvz); 
    }
    {
        x_j = _mm512_permutex_pd(x_j, _MM_PERM_BACD); // within 256
        y_j = _mm512_permutex_pd(y_j, _MM_PERM_BACD);
        z_j = _mm512_permutex_pd(z_j, _MM_PERM_BACD);
        m_j = _mm512_permutex_pd(m_j, _MM_PERM_ABDC); 

        const __m512d dx_j = _mm512_sub_pd(p512->hx, x_j);
        const __m512d dy_j = _mm512_sub_pd(p512->hy, y_j);
        const __m512d dz_j = _mm512_sub_pd(p512->hz, z_j);

        // 0123 4567
        // 1032 5476 
        const __m512d prefact = gravity_prefactor_avx512(m_j, dx_j, dy_j, dz_j);
        p512->hvx = _mm512_fnmadd_pd(prefact, dx_j, p512->hvx); 
        p512->hvy = _mm512_fnmadd_pd(prefact, dy_j, p512->hvy); 
        p512->hvz = _mm512_fnmadd_pd(prefact, dz_j, p512->hvz); 
    }

    // //////////////////////////////////////
    // 256 bit lane crossing
    // //////////////////////////////////////

    __m512d dvx; // delta vx for 4567 1230
    __m512d dvy;
    __m512d dvz;

    {
        x_j = _mm512_permutexvar_pd(_mm512_set_epi64(1,2,3,0,6,7,4,5), x_j); // accros 512
        y_j = _mm512_permutexvar_pd(_mm512_set_epi64(1,2,3,0,6,7,4,5), y_j);
        z_j = _mm512_permutexvar_pd(_mm512_set_epi64(1,2,3,0,6,7,4,5), z_j);
        m_j = _mm512_permutexvar_pd(_mm512_set_epi64(1,2,3,0,6,7,4,5), m_j);

        __m512d dx_j = _mm512_sub_pd(p512->hx, x_j);
        __m512d dy_j = _mm512_sub_pd(p512->hy, y_j);
        __m512d dz_j = _mm512_sub_pd(p512->hz, z_j);
        __m512d prefact = gravity_prefactor_avx512_one(dx_j, dy_j, dz_j);

        // 0123 4567
        // 4567 1230 
        __m512d prefact1 = _mm512_mul_pd(prefact, m_j);
        p512->hvx = _mm512_fnmadd_pd(prefact1, dx_j, p512->hvx); 
        p512->hvy = _mm512_fnmadd_pd(prefact1, dy_j, p512->hvy); 
        p512->hvz = _mm512_fnmadd_pd(prefact1, dz_j, p512->hvz); 


        // 4567 1230 
        // 0123 4567
        prefact = _mm512_mul_pd(prefact, m_j_01234567);
        dvx = _mm512_mul_pd(prefact, dx_j); 
        dvy = _mm512_mul_pd(prefact, dy_j); 
        dvz = _mm512_mul_pd(prefact, dz_j); 

    }

    {
        x_j = _mm512_permutex_pd(x_j, _MM_PERM_ADCB); // within 256
        y_j = _mm512_permutex_pd(y_j, _MM_PERM_ADCB);
        z_j = _mm512_permutex_pd(z_j, _MM_PERM_ADCB);
        m_j = _mm512_permutex_pd(m_j, _MM_PERM_ADCB);

        __m512d dx_j = _mm512_sub_pd(p512->hx, x_j);
        __m512d dy_j = _mm512_sub_pd(p512->hy, y_j);
        __m512d dz_j = _mm512_sub_pd(p512->hz, z_j);
        __m512d prefact = gravity_prefactor_avx512_one(dx_j, dy_j, dz_j);

        // 0123 4567
        // 5674 2301 
        __m512d prefact1 = _mm512_mul_pd(prefact, m_j);
        p512->hvx = _mm512_fnmadd_pd(prefact1, dx_j, p512->hvx); 
        p512->hvy = _mm512_fnmadd_pd(prefact1, dy_j, p512->hvy); 
        p512->hvz = _mm512_fnmadd_pd(prefact1, dz_j, p512->hvz); 



        dx_j         = _mm512_permutex_pd(dx_j,         _MM_PERM_CBAD); // within 256
        dy_j         = _mm512_permutex_pd(dy_j,         _MM_PERM_CBAD);
        dz_j         = _mm512_permutex_pd(dz_j,         _MM_PERM_CBAD);
        prefact      = _mm512_permutex_pd(prefact,      _MM_PERM_CBAD);
        m_j_01234567 = _mm512_permutex_pd(m_j_01234567, _MM_PERM_CBAD);

        // 4567 1230 
        // 3012 7456
        prefact = _mm512_mul_pd(prefact, m_j_01234567);
        dvx = _mm512_fmadd_pd(prefact, dx_j, dvx); 
        dvy = _mm512_fmadd_pd(prefact, dy_j, dvy); 
        dvz = _mm512_fmadd_pd(prefact, dz_j, dvz); 

    }

    // //////////////////////////////////////
    // 256 bit lane crossing for final add
    // //////////////////////////////////////

    {
        dvx = _mm512_permutexvar_pd(_mm512_set_epi64(3,2,1,0,6,5,4,7), dvx); //across 512
        dvy = _mm512_permutexvar_pd(_mm512_set_epi64(3,2,1,0,6,5,4,7), dvy);
        dvz = _mm512_permutexvar_pd(_mm512_set_epi64(3,2,1,0,6,5,4,7), dvz);

        p512->hvx = _mm512_maskz_add_pd(p512->mask, dvx, p512->hvx); 
        p512->hvy = _mm512_maskz_add_pd(p512->mask, dvy, p512->hvy); 
        p512->hvz = _mm512_maskz_add_pd(p512->mask, dvz, p512->hvz); 
    }
    
    // Convert accelerations (delta v) from heliocentric to Jacobi. Note: no difference between inertial and heliocentric here.
    mat8_mul3_avx512(ri_whfast512->mat8_inertial_to_jacobi, 
            p512->hvx, p512->hvy, p512->hvz,
            &dvx, &dvy, &dvz);
    ri_whfast512->p512->vx = _mm512_add_pd(dvx, ri_whfast512->p512->vx); 
    ri_whfast512->p512->vy = _mm512_add_pd(dvy, ri_whfast512->p512->vy); 
    ri_whfast512->p512->vz = _mm512_add_pd(dvz, ri_whfast512->p512->vz); 
    // Add Jacobi term using Jacobi coordinates
    __m512d prefact2 = gravity_prefactor_avx512_one(ri_whfast512->p512->x, ri_whfast512->p512->y, ri_whfast512->p512->z);
    __m512d prefact3 = _mm512_mul_pd(prefact2, ri_whfast512->p512->M); // Note: now M which is m0, m0+m1, m0+m1+m2, ...
    prefact3 = _mm512_mul_pd(prefact3, dt512);
    ri_whfast512->p512->vx = _mm512_maskz_fmadd_pd(p512->mask, prefact3, ri_whfast512->p512->x, ri_whfast512->p512->vx); 
    ri_whfast512->p512->vy = _mm512_maskz_fmadd_pd(p512->mask, prefact3, ri_whfast512->p512->y, ri_whfast512->p512->vy); 
    ri_whfast512->p512->vz = _mm512_maskz_fmadd_pd(p512->mask, prefact3, ri_whfast512->p512->z, ri_whfast512->p512->vz); 

#ifdef PROF
    struct reb_timeval time_end;
    gettimeofday(&time_end,NULL);
    walltime_interaction += time_end.tv_sec-time_beginning.tv_sec+(time_end.tv_usec-time_beginning.tv_usec)/1e6;
#endif
}
static void reb_whfast512_interaction_step_4planets_jacobi(struct reb_simulation * r){
#ifdef PROF
    struct reb_timeval time_beginning;
    gettimeofday(&time_beginning,NULL);
#endif
    struct reb_integrator_whfast512* const ri_whfast512 = &(r->ri_whfast512);
    struct reb_particle_avx512* restrict p512 = ri_whfast512->p512;

    // Convert position from Jacobi to heliocentric
    mat8_mul3_avx512(ri_whfast512->mat8_jacobi_to_heliocentric,
            p512->x, p512->y, p512->z,
            &p512->hx, &p512->hy, &p512->hz);
    __m512d x_j =  p512->hx;
    __m512d y_j =  p512->hy;
    __m512d z_j =  p512->hz;
    __m512d dt512 = _mm512_set1_pd(r->dt); 
    
    // General relativistic corrections
    if (ri_whfast512->gr_potential){
        __m512d r2 = _mm512_mul_pd(x_j, x_j);
        r2 = _mm512_fmadd_pd(y_j, y_j, r2);
        r2 = _mm512_fmadd_pd(z_j, z_j, r2);
        const __m512d r4 = _mm512_mul_pd(r2, r2);
        __m512d prefac = _mm512_div_pd(ri_whfast512->p512->gr_prefac,r4);
        __m512d dvx = _mm512_mul_pd(prefac, x_j); 
        __m512d dvy = _mm512_mul_pd(prefac, y_j); 
        __m512d dvz = _mm512_mul_pd(prefac, z_j); 
        p512->hvx  = dvx;
        p512->hvy  = dvy;
        p512->hvz  = dvz;
        
        // Calculate back reaction onto star and apply them to planets (heliocentric) 
        dvx = _mm512_maskz_mul_pd(p512->mask, ri_whfast512->p512->gr_prefac2, dvx); 
        dvy = _mm512_maskz_mul_pd(p512->mask, ri_whfast512->p512->gr_prefac2, dvy); 
        dvz = _mm512_maskz_mul_pd(p512->mask, ri_whfast512->p512->gr_prefac2, dvz); 

        dvx = _mm512_add_pd(_mm512_shuffle_pd(dvx, dvx, 0x55), dvx); // Swapping neighbouring elements
        dvx = _mm512_add_pd(_mm512_permutex_pd(dvx, _MM_PERM_ABCD), dvx);
        dvy = _mm512_add_pd(_mm512_shuffle_pd(dvy, dvy, 0x55), dvy);
        dvy = _mm512_add_pd(_mm512_permutex_pd(dvy, _MM_PERM_ABCD), dvy);
        dvz = _mm512_add_pd(_mm512_shuffle_pd(dvz, dvz, 0x55), dvz);
        dvz = _mm512_add_pd(_mm512_permutex_pd(dvz, _MM_PERM_ABCD), dvz);

        p512->hvx = _mm512_maskz_sub_pd(p512->mask, p512->hvx, dvx);
        p512->hvy = _mm512_maskz_sub_pd(p512->mask, p512->hvy, dvy);
        p512->hvz = _mm512_maskz_sub_pd(p512->mask, p512->hvz, dvz);
    
        // Jacobi additions:
        // Need to add stellar term. Easy: already in heliocentric coordinates.
        // TODO: Should put a mask on particle 1 as +/- cancels.
        const __m512d r1 = _mm512_sqrt_pd(r2);
        const __m512d r3 = _mm512_mul_pd(r1, r2);
        __m512d prefact = _mm512_div_pd(_mm512_set1_pd(r->dt*r->particles[0].m),r3); // Note: using m0, m0, m0, ... for mass here
        p512->hvx = _mm512_maskz_fnmadd_pd(p512->mask, prefact, x_j, p512->hvx); 
        p512->hvy = _mm512_maskz_fnmadd_pd(p512->mask, prefact, y_j, p512->hvy); 
        p512->hvz = _mm512_maskz_fnmadd_pd(p512->mask, prefact, z_j, p512->hvz); 
    }else{
        // Jacobi additions:
        // TODO: Should put a mask on particle 1 as +/- cancels.
        // Need to add stellar term. Easy: already in heliocentric coordinates.
        __m512d prefact = gravity_prefactor_avx512_one(x_j, y_j, z_j);
        __m512d prefact1 = _mm512_mul_pd(prefact, _mm512_set1_pd(-r->particles[0].m)); // Note: using m0, m0, m0, ...
        prefact1 = _mm512_mul_pd(prefact1, dt512);
        p512->hvx = _mm512_maskz_mul_pd(p512->mask, prefact1, x_j); 
        p512->hvy = _mm512_maskz_mul_pd(p512->mask, prefact1, y_j); 
        p512->hvz = _mm512_maskz_mul_pd(p512->mask, prefact1, z_j); 
    }

    __m512d dvx;
    __m512d dvy;
    __m512d dvz;
    __m512d m_j = _mm512_mul_pd(p512->m, dt512);

    {
        x_j = _mm512_permutex_pd(x_j, _MM_PERM_BACD); // within 256
        y_j = _mm512_permutex_pd(y_j, _MM_PERM_BACD);
        z_j = _mm512_permutex_pd(z_j, _MM_PERM_BACD);
        m_j = _mm512_permutex_pd(m_j, _MM_PERM_BACD);
        __m512d dx_j = _mm512_sub_pd(p512->hx, x_j);
        __m512d dy_j = _mm512_sub_pd(p512->hy, y_j);
        __m512d dz_j = _mm512_sub_pd(p512->hz, z_j);
        __m512d prefact = gravity_prefactor_avx512_one(dx_j, dy_j, dz_j);

        // 0123 4567
        // 3201 7645
        __m512d prefact1 = _mm512_mul_pd(prefact, m_j);
        p512->hvx = _mm512_maskz_fnmadd_pd(p512->mask, prefact1, dx_j, p512->hvx); 
        p512->hvy = _mm512_maskz_fnmadd_pd(p512->mask, prefact1, dy_j, p512->hvy); 
        p512->hvz = _mm512_maskz_fnmadd_pd(p512->mask, prefact1, dz_j, p512->hvz); 


        dx_j    = _mm512_permutex_pd(dx_j,    _MM_PERM_ABDC); // within 256
        dy_j    = _mm512_permutex_pd(dy_j,    _MM_PERM_ABDC);
        dz_j    = _mm512_permutex_pd(dz_j,    _MM_PERM_ABDC);
        prefact = _mm512_permutex_pd(prefact, _MM_PERM_ABDC);
        m_j     = _mm512_permute_pd(m_j,      0x55);    // within 128

        // 0123 4567
        // 2310 6754
        __m512d prefact2 = _mm512_mul_pd(prefact, m_j);
        p512->hvx = _mm512_maskz_fmadd_pd(p512->mask, prefact2, dx_j, p512->hvx); 
        p512->hvy = _mm512_maskz_fmadd_pd(p512->mask, prefact2, dy_j, p512->hvy); 
        p512->hvz = _mm512_maskz_fmadd_pd(p512->mask, prefact2, dz_j, p512->hvz); 
    }
    {
        x_j = _mm512_permutex_pd(x_j, _MM_PERM_BACD); // within 256
        y_j = _mm512_permutex_pd(y_j, _MM_PERM_BACD);
        z_j = _mm512_permutex_pd(z_j, _MM_PERM_BACD);
        m_j = _mm512_permutex_pd(m_j, _MM_PERM_ABDC); 

        const __m512d dx_j = _mm512_sub_pd(p512->hx, x_j);
        const __m512d dy_j = _mm512_sub_pd(p512->hy, y_j);
        const __m512d dz_j = _mm512_sub_pd(p512->hz, z_j);

        // 0123 4567
        // 1032 5476 
        const __m512d prefact = gravity_prefactor_avx512(m_j, dx_j, dy_j, dz_j);
        p512->hvx = _mm512_maskz_fnmadd_pd(p512->mask, prefact, dx_j, p512->hvx); 
        p512->hvy = _mm512_maskz_fnmadd_pd(p512->mask, prefact, dy_j, p512->hvy); 
        p512->hvz = _mm512_maskz_fnmadd_pd(p512->mask, prefact, dz_j, p512->hvz); 
    }
    
    // Convert accelerations (delta v) from heliocentric to Jacobi. Note: no difference between inertial and heliocentric here.
    mat8_mul3_avx512(ri_whfast512->mat8_inertial_to_jacobi, 
            p512->hvx, p512->hvy, p512->hvz,
            &dvx, &dvy, &dvz);
    ri_whfast512->p512->vx = _mm512_add_pd(dvx, ri_whfast512->p512->vx); 
    ri_whfast512->p512->vy = _mm512_add_pd(dvy, ri_whfast512->p512->vy); 
    ri_whfast512->p512->vz = _mm512_add_pd(dvz, ri_whfast512->p512->vz); 
    // Add Jacobi term using Jacobi coordinates
    __m512d prefact2 = gravity_prefactor_avx512_one(ri_whfast512->p512->x, ri_whfast512->p512->y, ri_whfast512->p512->z);
    __m512d prefact3 = _mm512_mul_pd(prefact2, ri_whfast512->p512->M); // Note: now M which is m0, m0+m1, m0+m1+m2, ...
    prefact3 = _mm512_mul_pd(prefact3, dt512);
    ri_whfast512->p512->vx = _mm512_maskz_fmadd_pd(p512->mask, prefact3, ri_whfast512->p512->x, ri_whfast512->p512->vx); 
    ri_whfast512->p512->vy = _mm512_maskz_fmadd_pd(p512->mask, prefact3, ri_whfast512->p512->y, ri_whfast512->p512->vy); 
    ri_whfast512->p512->vz = _mm512_maskz_fmadd_pd(p512->mask, prefact3, ri_whfast512->p512->z, ri_whfast512->p512->vz); 

#ifdef PROF
    struct reb_timeval time_end;
    gettimeofday(&time_end,NULL);
    walltime_interaction += time_end.tv_sec-time_beginning.tv_sec+(time_end.tv_usec-time_beginning.tv_usec)/1e6;
#endif
}
static void reb_whfast512_interaction_step_2planets_jacobi(struct reb_simulation * r){
#ifdef PROF
    struct reb_timeval time_beginning;
    gettimeofday(&time_beginning,NULL);
#endif
    struct reb_integrator_whfast512* const ri_whfast512 = &(r->ri_whfast512);
    struct reb_particle_avx512* restrict p512 = ri_whfast512->p512;

    // Convert position from Jacobi to heliocentric
    mat8_mul3_avx512(ri_whfast512->mat8_jacobi_to_heliocentric,
            p512->x, p512->y, p512->z,
            &p512->hx, &p512->hy, &p512->hz);
    __m512d x_j =  p512->hx;
    __m512d y_j =  p512->hy;
    __m512d z_j =  p512->hz;
    __m512d dt512 = _mm512_set1_pd(r->dt); 
    
    // General relativistic corrections
    if (ri_whfast512->gr_potential){
        __m512d r2 = _mm512_mul_pd(x_j, x_j);
        r2 = _mm512_fmadd_pd(y_j, y_j, r2);
        r2 = _mm512_fmadd_pd(z_j, z_j, r2);
        const __m512d r4 = _mm512_mul_pd(r2, r2);
        __m512d prefac = _mm512_div_pd(ri_whfast512->p512->gr_prefac,r4);
        __m512d dvx = _mm512_mul_pd(prefac, x_j); 
        __m512d dvy = _mm512_mul_pd(prefac, y_j); 
        __m512d dvz = _mm512_mul_pd(prefac, z_j); 
        p512->hvx  = dvx;
        p512->hvy  = dvy;
        p512->hvz  = dvz;

        // Calculate back reaction onto star and apply them to planets (heliocentric) 
        dvx = _mm512_maskz_mul_pd(p512->mask, ri_whfast512->p512->gr_prefac2, dvx); 
        dvy = _mm512_maskz_mul_pd(p512->mask, ri_whfast512->p512->gr_prefac2, dvy); 
        dvz = _mm512_maskz_mul_pd(p512->mask, ri_whfast512->p512->gr_prefac2, dvz); 
        
        dvx = _mm512_add_pd(_mm512_shuffle_pd(dvx, dvx, 0x55), dvx); // Swapping neighbouring elements
        dvy = _mm512_add_pd(_mm512_shuffle_pd(dvy, dvy, 0x55), dvy);
        dvz = _mm512_add_pd(_mm512_shuffle_pd(dvz, dvz, 0x55), dvz);
        
        p512->hvx  = _mm512_maskz_sub_pd(p512->mask, p512->hvx, dvx);
        p512->hvy  = _mm512_maskz_sub_pd(p512->mask, p512->hvy, dvy);
        p512->hvz  = _mm512_maskz_sub_pd(p512->mask, p512->hvz, dvz);
        
        // Jacobi additions:
        // Need to add stellar term. Easy: already in heliocentric coordinates.
        // TODO: Should put a mask on particle 1 as +/- cancels.
        const __m512d r1 = _mm512_sqrt_pd(r2);
        const __m512d r3 = _mm512_mul_pd(r1, r2);
        __m512d prefact = _mm512_div_pd(_mm512_set1_pd(r->dt*r->particles[0].m),r3); // Note: using m0, m0, m0, ... for mass here
        p512->hvx = _mm512_maskz_fnmadd_pd(p512->mask, prefact, x_j, p512->hvx); 
        p512->hvy = _mm512_maskz_fnmadd_pd(p512->mask, prefact, y_j, p512->hvy); 
        p512->hvz = _mm512_maskz_fnmadd_pd(p512->mask, prefact, z_j, p512->hvz); 
    }else{
        // Jacobi additions:
        // TODO: Should put a mask on particle 1 as +/- cancels.
        // Need to add stellar term. Easy: already in heliocentric coordinates.
        __m512d prefact = gravity_prefactor_avx512_one(x_j, y_j, z_j);
        __m512d prefact1 = _mm512_mul_pd(prefact, _mm512_set1_pd(-r->particles[0].m)); // Note: using m0, m0, m0, ...
        prefact1 = _mm512_mul_pd(prefact1, dt512);
        p512->hvx = _mm512_maskz_mul_pd(p512->mask, prefact1, x_j); 
        p512->hvy = _mm512_maskz_mul_pd(p512->mask, prefact1, y_j); 
        p512->hvz = _mm512_maskz_mul_pd(p512->mask, prefact1, z_j); 
    }

    __m512d m_j = _mm512_mul_pd(p512->m, dt512);

    {  
        __m512d dx_j = _mm512_sub_pd(p512->hx, _mm512_shuffle_pd(x_j, x_j, 0x55));
        __m512d dy_j = _mm512_sub_pd(p512->hy, _mm512_shuffle_pd(y_j, y_j, 0x55));
        __m512d dz_j = _mm512_sub_pd(p512->hz, _mm512_shuffle_pd(z_j, z_j, 0x55));
        m_j = _mm512_shuffle_pd(m_j, m_j, 0x55);

        __m512d prefact = gravity_prefactor_avx512_one(dx_j, dy_j, dz_j);

        __m512d prefact1 = _mm512_mul_pd(prefact, m_j);
        p512->hvx = _mm512_maskz_fnmadd_pd(p512->mask, prefact1, dx_j, p512->hvx); 
        p512->hvy = _mm512_maskz_fnmadd_pd(p512->mask, prefact1, dy_j, p512->hvy); 
        p512->hvz = _mm512_maskz_fnmadd_pd(p512->mask, prefact1, dz_j, p512->hvz); 
    }

    __m512d dvx;
    __m512d dvy;
    __m512d dvz;
    // Convert accelerations (delta v) from heliocentric to Jacobi. Note: no difference between inertial and heliocentric here.
    mat8_mul3_avx512(ri_whfast512->mat8_inertial_to_jacobi, 
            p512->hvx, p512->hvy, p512->hvz,
            &dvx, &dvy, &dvz);
    ri_whfast512->p512->vx = _mm512_add_pd(dvx, ri_whfast512->p512->vx); 
    ri_whfast512->p512->vy = _mm512_add_pd(dvy, ri_whfast512->p512->vy); 
    ri_whfast512->p512->vz = _mm512_add_pd(dvz, ri_whfast512->p512->vz); 
    // Add Jacobi term using Jacobi coordinates
    __m512d prefact2 = gravity_prefactor_avx512_one(ri_whfast512->p512->x, ri_whfast512->p512->y, ri_whfast512->p512->z);
    __m512d prefact3 = _mm512_mul_pd(prefact2, ri_whfast512->p512->M); // Note: now M which is m0, m0+m1, m0+m1+m2, ...
    prefact3 = _mm512_mul_pd(prefact3, dt512);
    ri_whfast512->p512->vx = _mm512_maskz_fmadd_pd(p512->mask, prefact3, ri_whfast512->p512->x, ri_whfast512->p512->vx); 
    ri_whfast512->p512->vy = _mm512_maskz_fmadd_pd(p512->mask, prefact3, ri_whfast512->p512->y, ri_whfast512->p512->vy); 
    ri_whfast512->p512->vz = _mm512_maskz_fmadd_pd(p512->mask, prefact3, ri_whfast512->p512->z, ri_whfast512->p512->vz); 

#ifdef PROF
    struct reb_timeval time_end;
    gettimeofday(&time_end,NULL);
    walltime_interaction += time_end.tv_sec-time_beginning.tv_sec+(time_end.tv_usec-time_beginning.tv_usec)/1e6;
#endif
}

static void reb_whfast512_interaction_step_8planets_democraticheliocentric(struct reb_simulation * r){
#ifdef PROF
    struct reb_timeval time_beginning;
    gettimeofday(&time_beginning,NULL);
#endif
    struct reb_integrator_whfast512* const ri_whfast512 = &(r->ri_whfast512);
    struct reb_particle_avx512* restrict p512 = ri_whfast512->p512;

    __m512d x_j =  p512->x;
    __m512d y_j =  p512->y;
    __m512d z_j =  p512->z;
    __m512d dt512 = _mm512_set1_pd(r->dt); 

    // General relativistic corrections
    if (ri_whfast512->gr_potential){
        __m512d r2 = _mm512_mul_pd(x_j, x_j);
        r2 = _mm512_fmadd_pd(y_j, y_j, r2);
        r2 = _mm512_fmadd_pd(z_j, z_j, r2);
        const __m512d r4 = _mm512_mul_pd(r2, r2);
        __m512d prefac = _mm512_div_pd(ri_whfast512->p512->gr_prefac,r4);
        __m512d dvx = _mm512_mul_pd(prefac, x_j); 
        __m512d dvy = _mm512_mul_pd(prefac, y_j); 
        __m512d dvz = _mm512_mul_pd(prefac, z_j); 
        p512->vx  = _mm512_sub_pd(p512->vx, dvx);
        p512->vy  = _mm512_sub_pd(p512->vy, dvy);
        p512->vz  = _mm512_sub_pd(p512->vz, dvz);

       // // Calculate back reaction onto star and apply them to planets (heliocentric) 
        dvx = _mm512_maskz_mul_pd(p512->mask, p512->gr_prefac2, dvx); 
        dvy = _mm512_maskz_mul_pd(p512->mask, p512->gr_prefac2, dvy); 
        dvz = _mm512_maskz_mul_pd(p512->mask, p512->gr_prefac2, dvz); 

        double sum_x = _mm512_reduce_add_pd(dvx);
        double sum_y = _mm512_reduce_add_pd(dvy);
        double sum_z = _mm512_reduce_add_pd(dvz);
        p512->vx = _mm512_maskz_sub_pd(p512->mask, p512->vx, _mm512_set1_pd(sum_x));
        p512->vy = _mm512_maskz_sub_pd(p512->mask, p512->vy, _mm512_set1_pd(sum_y));
        p512->vz = _mm512_maskz_sub_pd(p512->mask, p512->vz, _mm512_set1_pd(sum_z));
    }



    __m512d m_j = _mm512_mul_pd(p512->m, dt512);
    __m512d m_j_01234567 = m_j;

    {
        x_j = _mm512_permutex_pd(x_j, _MM_PERM_BACD); // within 256
        y_j = _mm512_permutex_pd(y_j, _MM_PERM_BACD);
        z_j = _mm512_permutex_pd(z_j, _MM_PERM_BACD);
        m_j = _mm512_permutex_pd(m_j, _MM_PERM_BACD);
        __m512d dx_j = _mm512_sub_pd(p512->x, x_j);
        __m512d dy_j = _mm512_sub_pd(p512->y, y_j);
        __m512d dz_j = _mm512_sub_pd(p512->z, z_j);
        __m512d prefact = gravity_prefactor_avx512_one(dx_j, dy_j, dz_j);

        // 0123 4567
        // 3201 7645
        __m512d prefact1 = _mm512_mul_pd(prefact, m_j);
        p512->vx = _mm512_maskz_fnmadd_pd(p512->mask, prefact1, dx_j, p512->vx); 
        p512->vy = _mm512_maskz_fnmadd_pd(p512->mask, prefact1, dy_j, p512->vy); 
        p512->vz = _mm512_maskz_fnmadd_pd(p512->mask, prefact1, dz_j, p512->vz); 


        dx_j    = _mm512_permutex_pd(dx_j,    _MM_PERM_ABDC); // within 256
        dy_j    = _mm512_permutex_pd(dy_j,    _MM_PERM_ABDC);
        dz_j    = _mm512_permutex_pd(dz_j,    _MM_PERM_ABDC);
        prefact = _mm512_permutex_pd(prefact, _MM_PERM_ABDC);
        m_j     = _mm512_permute_pd(m_j,      0x55);    // within 128

        // 0123 4567
        // 2310 6754
        __m512d prefact2 = _mm512_mul_pd(prefact, m_j);
        p512->vx = _mm512_maskz_fmadd_pd(p512->mask, prefact2, dx_j, p512->vx); 
        p512->vy = _mm512_maskz_fmadd_pd(p512->mask, prefact2, dy_j, p512->vy); 
        p512->vz = _mm512_maskz_fmadd_pd(p512->mask, prefact2, dz_j, p512->vz); 
    }
    {
        x_j = _mm512_permutex_pd(x_j, _MM_PERM_BACD); // within 256
        y_j = _mm512_permutex_pd(y_j, _MM_PERM_BACD);
        z_j = _mm512_permutex_pd(z_j, _MM_PERM_BACD);
        m_j = _mm512_permutex_pd(m_j, _MM_PERM_ABDC); 

        const __m512d dx_j = _mm512_sub_pd(p512->x, x_j);
        const __m512d dy_j = _mm512_sub_pd(p512->y, y_j);
        const __m512d dz_j = _mm512_sub_pd(p512->z, z_j);

        // 0123 4567
        // 1032 5476 
        const __m512d prefact = gravity_prefactor_avx512(m_j, dx_j, dy_j, dz_j);
        p512->vx = _mm512_maskz_fnmadd_pd(p512->mask, prefact, dx_j, p512->vx); 
        p512->vy = _mm512_maskz_fnmadd_pd(p512->mask, prefact, dy_j, p512->vy); 
        p512->vz = _mm512_maskz_fnmadd_pd(p512->mask, prefact, dz_j, p512->vz); 
    }

    // //////////////////////////////////////
    // 256 bit lane crossing
    // //////////////////////////////////////

    __m512d dvx; // delta vx for 4567 1230
    __m512d dvy;
    __m512d dvz;

    {
        x_j = _mm512_permutexvar_pd(_mm512_set_epi64(1,2,3,0,6,7,4,5), x_j); // accros 512
        y_j = _mm512_permutexvar_pd(_mm512_set_epi64(1,2,3,0,6,7,4,5), y_j);
        z_j = _mm512_permutexvar_pd(_mm512_set_epi64(1,2,3,0,6,7,4,5), z_j);
        m_j = _mm512_permutexvar_pd(_mm512_set_epi64(1,2,3,0,6,7,4,5), m_j);

        __m512d dx_j = _mm512_sub_pd(p512->x, x_j);
        __m512d dy_j = _mm512_sub_pd(p512->y, y_j);
        __m512d dz_j = _mm512_sub_pd(p512->z, z_j);
        __m512d prefact = gravity_prefactor_avx512_one(dx_j, dy_j, dz_j);

        // 0123 4567
        // 4567 1230 
        __m512d prefact1 = _mm512_mul_pd(prefact, m_j);
        p512->vx = _mm512_fnmadd_pd(prefact1, dx_j, p512->vx); 
        p512->vy = _mm512_fnmadd_pd(prefact1, dy_j, p512->vy); 
        p512->vz = _mm512_fnmadd_pd(prefact1, dz_j, p512->vz); 


        // 4567 1230 
        // 0123 4567
        prefact = _mm512_mul_pd(prefact, m_j_01234567);
        dvx = _mm512_mul_pd(prefact, dx_j); 
        dvy = _mm512_mul_pd(prefact, dy_j); 
        dvz = _mm512_mul_pd(prefact, dz_j); 

    }

    {
        x_j = _mm512_permutex_pd(x_j, _MM_PERM_ADCB); // within 256
        y_j = _mm512_permutex_pd(y_j, _MM_PERM_ADCB);
        z_j = _mm512_permutex_pd(z_j, _MM_PERM_ADCB);
        m_j = _mm512_permutex_pd(m_j, _MM_PERM_ADCB);

        __m512d dx_j = _mm512_sub_pd(p512->x, x_j);
        __m512d dy_j = _mm512_sub_pd(p512->y, y_j);
        __m512d dz_j = _mm512_sub_pd(p512->z, z_j);
        __m512d prefact = gravity_prefactor_avx512_one(dx_j, dy_j, dz_j);

        // 0123 4567
        // 5674 2301 
        __m512d prefact1 = _mm512_mul_pd(prefact, m_j);
        p512->vx = _mm512_maskz_fnmadd_pd(p512->mask, prefact1, dx_j, p512->vx); 
        p512->vy = _mm512_maskz_fnmadd_pd(p512->mask, prefact1, dy_j, p512->vy); 
        p512->vz = _mm512_maskz_fnmadd_pd(p512->mask, prefact1, dz_j, p512->vz); 



        dx_j         = _mm512_permutex_pd(dx_j,         _MM_PERM_CBAD); // within 256
        dy_j         = _mm512_permutex_pd(dy_j,         _MM_PERM_CBAD);
        dz_j         = _mm512_permutex_pd(dz_j,         _MM_PERM_CBAD);
        prefact      = _mm512_permutex_pd(prefact,      _MM_PERM_CBAD);
        m_j_01234567 = _mm512_permutex_pd(m_j_01234567, _MM_PERM_CBAD);

        // 4567 1230 
        // 3012 7456
        prefact = _mm512_mul_pd(prefact, m_j_01234567);
        dvx = _mm512_fmadd_pd(prefact, dx_j, dvx); 
        dvy = _mm512_fmadd_pd(prefact, dy_j, dvy); 
        dvz = _mm512_fmadd_pd(prefact, dz_j, dvz); 

    }

    // //////////////////////////////////////
    // 256 bit lane crossing for final add
    // //////////////////////////////////////

    {
        dvx = _mm512_permutexvar_pd(_mm512_set_epi64(3,2,1,0,6,5,4,7), dvx); //across 512
        dvy = _mm512_permutexvar_pd(_mm512_set_epi64(3,2,1,0,6,5,4,7), dvy);
        dvz = _mm512_permutexvar_pd(_mm512_set_epi64(3,2,1,0,6,5,4,7), dvz);

        p512->vx = _mm512_maskz_add_pd(p512->mask, dvx, p512->vx); 
        p512->vy = _mm512_maskz_add_pd(p512->mask, dvy, p512->vy); 
        p512->vz = _mm512_maskz_add_pd(p512->mask, dvz, p512->vz); 
    }

#ifdef PROF
    struct reb_timeval time_end;
    gettimeofday(&time_end,NULL);
    walltime_interaction += time_end.tv_sec-time_beginning.tv_sec+(time_end.tv_usec-time_beginning.tv_usec)/1e6;
#endif
}

// Performs one full interaction step
static void reb_whfast512_interaction_step_4planets_democraticheliocentric(struct reb_simulation * r){
#ifdef PROF
    struct reb_timeval time_beginning;
    gettimeofday(&time_beginning,NULL);
#endif
    struct reb_integrator_whfast512* const ri_whfast512 = &(r->ri_whfast512);
    struct reb_particle_avx512* restrict p512 = ri_whfast512->p512;

    __m512d x_j =  p512->x;
    __m512d y_j =  p512->y;
    __m512d z_j =  p512->z;
    __m512d dt512 = _mm512_set1_pd(r->dt); 

    // General relativistic corrections
    if (ri_whfast512->gr_potential){
        __m512d r2 = _mm512_mul_pd(x_j, x_j);
        r2 = _mm512_fmadd_pd(y_j, y_j, r2);
        r2 = _mm512_fmadd_pd(z_j, z_j, r2);
        const __m512d r4 = _mm512_mul_pd(r2, r2);
        __m512d prefac = _mm512_div_pd(ri_whfast512->p512->gr_prefac,r4);
        __m512d dvx = _mm512_mul_pd(prefac, x_j); 
        __m512d dvy = _mm512_mul_pd(prefac, y_j); 
        __m512d dvz = _mm512_mul_pd(prefac, z_j); 
        p512->vx  = _mm512_maskz_sub_pd(p512->mask, p512->vx, dvx);
        p512->vy  = _mm512_maskz_sub_pd(p512->mask, p512->vy, dvy);
        p512->vz  = _mm512_maskz_sub_pd(p512->mask, p512->vz, dvz);

        // Calculate back reaction onto star and apply them to planets (heliocentric) 
        dvx = _mm512_maskz_mul_pd(p512->mask, ri_whfast512->p512->gr_prefac2, dvx); 
        dvy = _mm512_maskz_mul_pd(p512->mask, ri_whfast512->p512->gr_prefac2, dvy); 
        dvz = _mm512_maskz_mul_pd(p512->mask, ri_whfast512->p512->gr_prefac2, dvz); 

        dvx = _mm512_add_pd(_mm512_shuffle_pd(dvx, dvx, 0x55), dvx); // Swapping neighbouring elements
        dvx = _mm512_add_pd(_mm512_permutex_pd(dvx, _MM_PERM_ABCD), dvx);
        dvy = _mm512_add_pd(_mm512_shuffle_pd(dvy, dvy, 0x55), dvy);
        dvy = _mm512_add_pd(_mm512_permutex_pd(dvy, _MM_PERM_ABCD), dvy);
        dvz = _mm512_add_pd(_mm512_shuffle_pd(dvz, dvz, 0x55), dvz);
        dvz = _mm512_add_pd(_mm512_permutex_pd(dvz, _MM_PERM_ABCD), dvz);

        p512->vx  = _mm512_maskz_sub_pd(p512->mask, p512->vx, dvx);
        p512->vy  = _mm512_maskz_sub_pd(p512->mask, p512->vy, dvy);
        p512->vz  = _mm512_maskz_sub_pd(p512->mask, p512->vz, dvz);
    }



    __m512d m_j = _mm512_mul_pd(p512->m, dt512);

    {
        x_j = _mm512_permutex_pd(x_j, _MM_PERM_BACD); // within 256
        y_j = _mm512_permutex_pd(y_j, _MM_PERM_BACD);
        z_j = _mm512_permutex_pd(z_j, _MM_PERM_BACD);
        m_j = _mm512_permutex_pd(m_j, _MM_PERM_BACD);
        __m512d dx_j = _mm512_sub_pd(p512->x, x_j);
        __m512d dy_j = _mm512_sub_pd(p512->y, y_j);
        __m512d dz_j = _mm512_sub_pd(p512->z, z_j);
        __m512d prefact = gravity_prefactor_avx512_one(dx_j, dy_j, dz_j);

        // 0123 4567
        // 3201 7645
        __m512d prefact1 = _mm512_mul_pd(prefact, m_j);
        p512->vx = _mm512_maskz_fnmadd_pd(p512->mask, prefact1, dx_j, p512->vx); 
        p512->vy = _mm512_maskz_fnmadd_pd(p512->mask, prefact1, dy_j, p512->vy); 
        p512->vz = _mm512_maskz_fnmadd_pd(p512->mask, prefact1, dz_j, p512->vz); 


        dx_j    = _mm512_permutex_pd(dx_j,    _MM_PERM_ABDC); // within 256
        dy_j    = _mm512_permutex_pd(dy_j,    _MM_PERM_ABDC);
        dz_j    = _mm512_permutex_pd(dz_j,    _MM_PERM_ABDC);
        prefact = _mm512_permutex_pd(prefact, _MM_PERM_ABDC);
        m_j     = _mm512_permute_pd(m_j,      0x55);    // within 128

        // 0123 4567
        // 2310 6754
        __m512d prefact2 = _mm512_mul_pd(prefact, m_j);
        p512->vx = _mm512_maskz_fmadd_pd(p512->mask, prefact2, dx_j, p512->vx); 
        p512->vy = _mm512_maskz_fmadd_pd(p512->mask, prefact2, dy_j, p512->vy); 
        p512->vz = _mm512_maskz_fmadd_pd(p512->mask, prefact2, dz_j, p512->vz); 
    }
    {
        x_j = _mm512_permutex_pd(x_j, _MM_PERM_BACD); // within 256
        y_j = _mm512_permutex_pd(y_j, _MM_PERM_BACD);
        z_j = _mm512_permutex_pd(z_j, _MM_PERM_BACD);
        m_j = _mm512_permutex_pd(m_j, _MM_PERM_ABDC); 

        const __m512d dx_j = _mm512_sub_pd(p512->x, x_j);
        const __m512d dy_j = _mm512_sub_pd(p512->y, y_j);
        const __m512d dz_j = _mm512_sub_pd(p512->z, z_j);

        // 0123 4567
        // 1032 5476 
        const __m512d prefact = gravity_prefactor_avx512(m_j, dx_j, dy_j, dz_j);
        p512->vx = _mm512_maskz_fnmadd_pd(p512->mask, prefact, dx_j, p512->vx); 
        p512->vy = _mm512_maskz_fnmadd_pd(p512->mask, prefact, dy_j, p512->vy); 
        p512->vz = _mm512_maskz_fnmadd_pd(p512->mask, prefact, dz_j, p512->vz); 
    }


#ifdef PROF
    struct reb_timeval time_end;
    gettimeofday(&time_end,NULL);
    walltime_interaction += time_end.tv_sec-time_beginning.tv_sec+(time_end.tv_usec-time_beginning.tv_usec)/1e6;
#endif
}

static void reb_whfast512_interaction_step_2planets_democraticheliocentric(struct reb_simulation * r){
#ifdef PROF
    struct reb_timeval time_beginning;
    gettimeofday(&time_beginning,NULL);
#endif
    struct reb_integrator_whfast512* const ri_whfast512 = &(r->ri_whfast512);
    struct reb_particle_avx512* restrict p512 = ri_whfast512->p512;

    __m512d x_j =  p512->x;
    __m512d y_j =  p512->y;
    __m512d z_j =  p512->z;
    __m512d dt512 = _mm512_set1_pd(r->dt); 

    // General relativistic corrections
    if (ri_whfast512->gr_potential){
        __m512d r2 = _mm512_mul_pd(x_j, x_j);
        r2 = _mm512_fmadd_pd(y_j, y_j, r2);
        r2 = _mm512_fmadd_pd(z_j, z_j, r2);
        const __m512d r4 = _mm512_mul_pd(r2, r2);
        __m512d prefac = _mm512_div_pd(ri_whfast512->p512->gr_prefac,r4);
        __m512d dvx = _mm512_mul_pd(prefac, x_j); 
        __m512d dvy = _mm512_mul_pd(prefac, y_j); 
        __m512d dvz = _mm512_mul_pd(prefac, z_j); 
        p512->vx  = _mm512_maskz_sub_pd(p512->mask, p512->vx, dvx);
        p512->vy  = _mm512_maskz_sub_pd(p512->mask, p512->vy, dvy);
        p512->vz  = _mm512_maskz_sub_pd(p512->mask, p512->vz, dvz);

        // Calculate back reaction onto star and apply them to planets (heliocentric) 
        dvx = _mm512_maskz_mul_pd(p512->mask, ri_whfast512->p512->gr_prefac2, dvx); 
        dvy = _mm512_maskz_mul_pd(p512->mask, ri_whfast512->p512->gr_prefac2, dvy); 
        dvz = _mm512_maskz_mul_pd(p512->mask, ri_whfast512->p512->gr_prefac2, dvz); 

        dvx = _mm512_add_pd(_mm512_shuffle_pd(dvx, dvx, 0x55), dvx); // Swapping neighbouring elements
        dvy = _mm512_add_pd(_mm512_shuffle_pd(dvy, dvy, 0x55), dvy);
        dvz = _mm512_add_pd(_mm512_shuffle_pd(dvz, dvz, 0x55), dvz);

        p512->vx  = _mm512_maskz_sub_pd(p512->mask, p512->vx, dvx);
        p512->vy  = _mm512_maskz_sub_pd(p512->mask, p512->vy, dvy);
        p512->vz  = _mm512_maskz_sub_pd(p512->mask, p512->vz, dvz);
    }



    __m512d m_j = _mm512_mul_pd(p512->m, dt512);

    {  
        __m512d dx_j = _mm512_sub_pd(p512->x, _mm512_shuffle_pd(x_j, x_j, 0x55));
        __m512d dy_j = _mm512_sub_pd(p512->y, _mm512_shuffle_pd(y_j, y_j, 0x55));
        __m512d dz_j = _mm512_sub_pd(p512->z, _mm512_shuffle_pd(z_j, z_j, 0x55));
        m_j = _mm512_shuffle_pd(m_j, m_j, 0x55);

        __m512d prefact = gravity_prefactor_avx512_one(dx_j, dy_j, dz_j);

        __m512d prefact1 = _mm512_mul_pd(prefact, m_j);
        p512->vx = _mm512_maskz_fnmadd_pd(p512->mask, prefact1, dx_j, p512->vx); 
        p512->vy = _mm512_maskz_fnmadd_pd(p512->mask, prefact1, dy_j, p512->vy); 
        p512->vz = _mm512_maskz_fnmadd_pd(p512->mask, prefact1, dz_j, p512->vz); 
    }


#ifdef PROF
    struct reb_timeval time_end;
    gettimeofday(&time_end,NULL);
    walltime_interaction += time_end.tv_sec-time_beginning.tv_sec+(time_end.tv_usec-time_beginning.tv_usec)/1e6;
#endif
}



// Convert inertial coordinates to democratic heliocentric coordinates
// Note: this is only called at the beginning. Speed is not a concern.
static void inertial_to_democratic_heliocentric_posvel(struct reb_simulation* r){
    struct reb_integrator_whfast512* const ri_whfast512 = &(r->ri_whfast512);
    struct reb_particle* particles = r->particles;
    const unsigned int N_systems = ri_whfast512->N_systems;
    const unsigned int p_per_system = 8/N_systems;
    const unsigned int N_per_system = r->N/N_systems;

    // Layout (2x 3 planet systems)
    //                    0    1  2  3   4    5  6  7
    // particles array    star p0 p1 p2  star p0 p1 p2
    // avx512 register    p0   p1 p2 NA  p0   p1 p2 NA

    // Layout (4x 2 planet systems)
    //                    0    1  2  3    4  5  6    7  8  9    10 11
    // particles array    star p0 p1 star p0 p1 star p0 p1 star p0 p1
    // avx512 register    p0   p1 p0 p1   p0 p1 p0   p1  

    double m[8];
    double x[8];
    double y[8];
    double z[8];
    double vx[8];
    double vy[8];
    double vz[8];
    for (unsigned s=0;s<N_systems;s++){
        double mtot = 0;
        double x0 = 0; // center of mass
        double y0 = 0;
        double z0 = 0;
        double vx0 = 0;
        double vy0 = 0;
        double vz0 = 0;
        for (unsigned int i=0; i<N_per_system; i++){
            mtot += particles[s*N_per_system+i].m;
            x0 += particles[s*N_per_system+i].x * particles[s*N_per_system+i].m;
            y0 += particles[s*N_per_system+i].y * particles[s*N_per_system+i].m;
            z0 += particles[s*N_per_system+i].z * particles[s*N_per_system+i].m;
            vx0 += particles[s*N_per_system+i].vx * particles[s*N_per_system+i].m;
            vy0 += particles[s*N_per_system+i].vy * particles[s*N_per_system+i].m;
            vz0 += particles[s*N_per_system+i].vz * particles[s*N_per_system+i].m;
        }
        for (unsigned int i=0; i<p_per_system; i++){
            m[s*p_per_system+i] = 0;
            x[s*p_per_system+i] = 0;
            y[s*p_per_system+i] = 0;
            z[s*p_per_system+i] = 0;
            vx[s*p_per_system+i] = 0.0;
            vy[s*p_per_system+i] = 0.0;
            vz[s*p_per_system+i] = 0.0;
        }
        ri_whfast512->p_jh0[s].m = mtot;
        ri_whfast512->p_jh0[s].x = x0/mtot;
        ri_whfast512->p_jh0[s].y = y0/mtot;
        ri_whfast512->p_jh0[s].z = z0/mtot;
        ri_whfast512->p_jh0[s].vx = vx0/mtot;
        ri_whfast512->p_jh0[s].vy = vy0/mtot;
        ri_whfast512->p_jh0[s].vz = vz0/mtot;
        for (unsigned int i=1; i<N_per_system; i++){
            m[s*p_per_system+(i-1)] = particles[s*N_per_system+i].m;
            x[s*p_per_system+(i-1)] = particles[s*N_per_system+i].x - particles[s*N_per_system].x; // heliocentric
            y[s*p_per_system+(i-1)] = particles[s*N_per_system+i].y - particles[s*N_per_system].y;
            z[s*p_per_system+(i-1)] = particles[s*N_per_system+i].z - particles[s*N_per_system].z;
            vx[s*p_per_system+(i-1)] = particles[s*N_per_system+i].vx - ri_whfast512->p_jh0[s].vx; // relative to com
            vy[s*p_per_system+(i-1)] = particles[s*N_per_system+i].vy - ri_whfast512->p_jh0[s].vy;
            vz[s*p_per_system+(i-1)] = particles[s*N_per_system+i].vz - ri_whfast512->p_jh0[s].vz;
        }
    }

    struct reb_particle_avx512* p512 = ri_whfast512->p512;
    p512->m = _mm512_loadu_pd(m);
    p512->x = _mm512_loadu_pd(x);
    p512->y = _mm512_loadu_pd(y);
    p512->z = _mm512_loadu_pd(z);
    p512->vx = _mm512_loadu_pd(vx);
    p512->vy = _mm512_loadu_pd(vy);
    p512->vz = _mm512_loadu_pd(vz);
}

static void inertial_to_jacobi_posvel(struct reb_simulation* r){
    struct reb_integrator_whfast512* const ri_whfast512 = &(r->ri_whfast512);
    const unsigned int N_systems = ri_whfast512->N_systems;
    const unsigned int N_per_system = r->N/N_systems;

    // Same layout as for democratic heliocentric
    for (int s=0; s<N_systems; s++){
        struct reb_particle com = reb_simulation_com_range(r,s*N_per_system, (s+1)*N_per_system);
        ri_whfast512->p_jh0[s] = com;
    }
    ri_whfast512->p512->x = load_into_m512d(r, offsetof(struct reb_particle,x),ri_whfast512->mat8_inertial_to_jacobi, r->N);
    ri_whfast512->p512->y = load_into_m512d(r, offsetof(struct reb_particle,y),ri_whfast512->mat8_inertial_to_jacobi, r->N);
    ri_whfast512->p512->z = load_into_m512d(r, offsetof(struct reb_particle,z),ri_whfast512->mat8_inertial_to_jacobi, r->N);
    ri_whfast512->p512->vx = load_into_m512d(r, offsetof(struct reb_particle,vx),ri_whfast512->mat8_inertial_to_jacobi, r->N);
    ri_whfast512->p512->vy = load_into_m512d(r, offsetof(struct reb_particle,vy),ri_whfast512->mat8_inertial_to_jacobi, r->N);
    ri_whfast512->p512->vz = load_into_m512d(r, offsetof(struct reb_particle,vz),ri_whfast512->mat8_inertial_to_jacobi, r->N);
    ri_whfast512->p512->m = load_into_m512d(r, offsetof(struct reb_particle,m),NULL, r->N);
}


// Performs one complete jump step
static void reb_whfast512_jump_step(struct reb_simulation* r, const double _dt){
#ifdef PROF
    struct reb_timeval time_beginning;
    gettimeofday(&time_beginning,NULL);
#endif
    struct reb_integrator_whfast512* ri_whfast512 = &(r->ri_whfast512);
    struct reb_particle_avx512* p512 = ri_whfast512->p512;
    double m0 = r->particles[0].m;

    __m512d pf512 = _mm512_set1_pd(_dt/m0);

    __m512d sumx = _mm512_mul_pd(p512->m, p512->vx);
    __m512d sumy = _mm512_mul_pd(p512->m, p512->vy);
    __m512d sumz = _mm512_mul_pd(p512->m, p512->vz);

    if (ri_whfast512->N_systems == 1){
        sumx = _mm512_add_pd(_mm512_shuffle_pd(sumx, sumx, 0x55), sumx); // Swapping neighbouring elements
        sumx = _mm512_add_pd(_mm512_permutex_pd(sumx, _MM_PERM_ABCD), sumx);
        sumx = _mm512_add_pd(_mm512_shuffle_f64x2(sumx,sumx, 78), sumx); // 78 is _MM_SHUFFLE(1,0,3,2), changed for icx

        sumy = _mm512_add_pd(_mm512_shuffle_pd(sumy, sumy, 0x55), sumy);
        sumy = _mm512_add_pd(_mm512_permutex_pd(sumy, _MM_PERM_ABCD), sumy);
        sumy = _mm512_add_pd(_mm512_shuffle_f64x2(sumy,sumy, 78), sumy);

        sumz = _mm512_add_pd(_mm512_shuffle_pd(sumz, sumz, 0x55), sumz);
        sumz = _mm512_add_pd(_mm512_permutex_pd(sumz, _MM_PERM_ABCD), sumz);
        sumz = _mm512_add_pd(_mm512_shuffle_f64x2(sumz,sumz, 78), sumz);
    }else if (ri_whfast512->N_systems == 2){
        sumx = _mm512_add_pd(_mm512_shuffle_pd(sumx, sumx, 0x55), sumx); // Swapping neighbouring elements
        sumx = _mm512_add_pd(_mm512_permutex_pd(sumx, _MM_PERM_ABCD), sumx);

        sumy = _mm512_add_pd(_mm512_shuffle_pd(sumy, sumy, 0x55), sumy);
        sumy = _mm512_add_pd(_mm512_permutex_pd(sumy, _MM_PERM_ABCD), sumy);

        sumz = _mm512_add_pd(_mm512_shuffle_pd(sumz, sumz, 0x55), sumz);
        sumz = _mm512_add_pd(_mm512_permutex_pd(sumz, _MM_PERM_ABCD), sumz);
    }else if (ri_whfast512->N_systems == 4){
        sumx = _mm512_add_pd(_mm512_shuffle_pd(sumx, sumx, 0x55), sumx); // Swapping neighbouring elements
        sumy = _mm512_add_pd(_mm512_shuffle_pd(sumy, sumy, 0x55), sumy);
        sumz = _mm512_add_pd(_mm512_shuffle_pd(sumz, sumz, 0x55), sumz);
    }

    p512->x = _mm512_maskz_fmadd_pd(p512->mask, sumx, pf512, p512->x); 
    p512->y = _mm512_maskz_fmadd_pd(p512->mask, sumy, pf512, p512->y); 
    p512->z = _mm512_maskz_fmadd_pd(p512->mask, sumz, pf512, p512->z); 

#ifdef PROF
    struct reb_timeval time_end;
    gettimeofday(&time_end,NULL);
    walltime_jump += time_end.tv_sec-time_beginning.tv_sec+(time_end.tv_usec-time_beginning.tv_usec)/1e6;
#endif
}

// Precalculate various constants and put them in 512 bit vectors.
void static recalculate_constants(struct reb_simulation* r){
    struct reb_integrator_whfast512* const ri_whfast512 = &(r->ri_whfast512);
    struct reb_particle* particles = r->particles;
    const unsigned int N_systems = ri_whfast512->N_systems;
    const unsigned int p_per_system = 8/N_systems;
    const unsigned int N_per_system = r->N/N_systems;
    double mat8_inertial_to_heliocentric[64];
    double M[8] = {0.0};
    switch (N_systems){
        case 1:
            ri_whfast512->p512->mask = (1 << (r->N -1)) - 1;
            break;
        case 2:
            if (N_per_system==5) ri_whfast512->p512->mask = 0xFF;
            if (N_per_system==4) ri_whfast512->p512->mask = 0x77;
            if (N_per_system==3) ri_whfast512->p512->mask = 0x33;
            if (N_per_system==2) ri_whfast512->p512->mask = 0x11;
            break;
        case 4:
            if (N_per_system==3) ri_whfast512->p512->mask = 0xFF;
            if (N_per_system==2) ri_whfast512->p512->mask = 0x55;
            break;
        default:
            reb_simulation_error(r,"Invalid value for N_systems.");
    }
    switch (ri_whfast512->coordinates){
        case REB_WHFAST512_COORDINATES_DEMOCRATICHELIOCENTRIC:
            for (int s=0; s<N_systems; s++){
                for (int p=0; p<p_per_system; p++){ // loop over all planets
                    M[s*p_per_system+p] = r->particles[s*N_per_system].m;
                }
            }
            break;
        case REB_WHFAST512_COORDINATES_JACOBI:
            // Allocate memory
            if (ri_whfast512->mat8_inertial_to_jacobi==NULL){
                ri_whfast512->mat8_inertial_to_jacobi = aligned_alloc(64, sizeof(double)*64);
            }
            if (ri_whfast512->mat8_jacobi_to_inertial==NULL){
                ri_whfast512->mat8_jacobi_to_inertial = aligned_alloc(64, sizeof(double)*64);
            }
            if (ri_whfast512->mat8_jacobi_to_heliocentric==NULL){
                ri_whfast512->mat8_jacobi_to_heliocentric = aligned_alloc(64, sizeof(double)*64);
            }
            // Zeroing.
            for (unsigned int i=1; i<9; i++){
                for (unsigned int j=1; j<9; j++){
                    ri_whfast512->mat8_inertial_to_jacobi[(i-1)+8*(j-1)] = 0.0;
                    ri_whfast512->mat8_jacobi_to_inertial[(i-1)+8*(j-1)] = 0.0;
                    ri_whfast512->mat8_jacobi_to_heliocentric[(i-1)+8*(j-1)] = 0.0;
                    mat8_inertial_to_heliocentric[(i-1)+8*(j-1)] = 0.0;
                }
            }
            
            // Filling vector
            for (int s=0; s<N_systems; s++){
                for (int i=0;i<N_per_system-1;i++){
                    for (int j=0;j<i+2;j++){
                        M[s*p_per_system+i] += r->particles[s*N_per_system+j].m;
                    }
                }
            }

            // Fill matricies
            for (int s=0; s<N_systems; s++){
                double ms = particles[s*N_per_system+0].m;
                for (unsigned int i=1; i<N_per_system; i++){
                    for (unsigned int j=i; j<N_per_system; j++){
                        ri_whfast512->mat8_inertial_to_jacobi[(s*p_per_system+i-1)+8*(s*p_per_system+j-1)] += particles[s*N_per_system+j].m/ms;
                    }
                    for (unsigned int j=1; j<N_per_system; j++){
                        mat8_inertial_to_heliocentric[(s*p_per_system+i-1)+8*(s*p_per_system+j-1)] += particles[s*N_per_system+j].m/particles[s*N_per_system+0].m;
                    }
                    mat8_inertial_to_heliocentric[(s*p_per_system+i-1)+8*(s*p_per_system+i-1)] += 1.0;
                    ri_whfast512->mat8_inertial_to_jacobi[(s*p_per_system+i-1)+8*(s*p_per_system+i-1)] += 1.0;
                    ri_whfast512->mat8_jacobi_to_inertial[(s*p_per_system+i-1)+8*(s*p_per_system+i-1)] += ms/(ms + particles[s*N_per_system+i].m);
                    ms += particles[s*N_per_system+i].m;
                    if (i<N_per_system-1){ 
                        for (unsigned int ii=i; ii>0; ii--){
                            int jj = i+1;
                            ri_whfast512->mat8_jacobi_to_inertial[(s*p_per_system+ii-1)+8*(s*p_per_system+jj-1)] -= particles[s*N_per_system+jj].m/(ms+particles[s*N_per_system+jj].m);
                        }
                    }
                }
            }

            // Might be numerically more stable to calculate this manually rather than do a matrix multiplication.
            for (unsigned int i=1; i<9; i++){
                for (unsigned int j=1; j<9; j++){
                    for (unsigned int k=1; k<9; k++){
                        ri_whfast512->mat8_jacobi_to_heliocentric[(i-1)+8*(j-1)] += mat8_inertial_to_heliocentric[(i-1)+8*(k-1)] * ri_whfast512->mat8_jacobi_to_inertial[(k-1)+8*(j-1)];
                    }
                }
            }
            break;
        default:
            reb_simulation_error(r,"Coordinate system not supported.");
    }

    ri_whfast512->p512->M = _mm512_loadu_pd(&M);

    // GR prefactors. Note: assumes units of AU, year/2pi.
    double c = 10065.32;
    double _gr_prefac[8];
    double _gr_prefac2[8];
    for(unsigned int i=0;i<8;i++){
        _gr_prefac[i] = 0; // for when N<8
        _gr_prefac2[i] = 0;
    }
    for (int s=0; s<N_systems; s++){
        double m0 = r->particles[s*N_per_system].m;
        for (int p=1; p<N_per_system; p++){
            switch (ri_whfast512->coordinates){
                case REB_WHFAST512_COORDINATES_JACOBI:
                    _gr_prefac[s*p_per_system+(p-1)] = -r->dt*6.*m0*m0/(c*c); // Minus sign!
                    _gr_prefac2[s*p_per_system+(p-1)] = -r->particles[s*N_per_system+p].m/m0;
                    break;
                case REB_WHFAST512_COORDINATES_DEMOCRATICHELIOCENTRIC:
                    _gr_prefac[s*p_per_system+(p-1)] = r->dt*6.*m0*m0/(c*c);
                    _gr_prefac2[s*p_per_system+(p-1)] = r->particles[s*N_per_system+p].m/m0;
                    break;
            }

        }
    }
    ri_whfast512->p512->gr_prefac = _mm512_loadu_pd(&_gr_prefac);
    ri_whfast512->p512->gr_prefac2 = _mm512_loadu_pd(&_gr_prefac2);

    ri_whfast512->dt_cached = r->dt;
    ri_whfast512->recalculate_constants = 0;
}

void static reb_integrator_whfast512_allocate(struct reb_simulation* const r){
    struct reb_integrator_whfast512* const ri_whfast512 = &(r->ri_whfast512);
    // Check if all assumptions are satisfied.
    // Note: These are not checked every timestep. 
    // So it is possible for the user to screw things up.
    if (r->dt<0.0){
        reb_simulation_error(r, "WHFast512 does not support negative timesteps. To integrate backwards, flip the sign of the velocities.");
        r->status = REB_STATUS_GENERIC_ERROR;
        return;
    }
    if (r->N_var!=0){
        reb_simulation_error(r, "WHFast512 does not support variational particles.");
        r->status = REB_STATUS_GENERIC_ERROR;
        return;
    }
    if (r->exact_finish_time!=0){
        reb_simulation_error(r, "WHFast512 requires exact_finish_time=0.");
        r->status = REB_STATUS_GENERIC_ERROR;
        return;
    }
    if (r->N>9 && ri_whfast512->N_systems == 1) {
        reb_simulation_error(r, "WHFast512 supports a maximum of 9 particles when N_systems is set to 1.");
        r->status = REB_STATUS_GENERIC_ERROR;
        return;
    }
    if (r->N>10 && ri_whfast512->N_systems == 2) {
        reb_simulation_error(r, "WHFast512 supports a maximum of 10 particles when N_systems is set to 2.");
        r->status = REB_STATUS_GENERIC_ERROR;
        return;
    }
    if (r->N>12 && ri_whfast512->N_systems == 4) {
        reb_simulation_error(r, "WHFast512 supports a maximum of 12 particles when N_systems is set to 4.");
        r->status = REB_STATUS_GENERIC_ERROR;
        return;
    }
    if (ri_whfast512->N_systems != 1 && ri_whfast512->N_systems !=2 && ri_whfast512->N_systems != 4){
        reb_simulation_error(r, "WHFast512 supports 1, 2, or 4 systems only.");
        r->status = REB_STATUS_GENERIC_ERROR;
        return;
    }
    if (r->N % ri_whfast512->N_systems != 0){
        reb_simulation_error(r, "Number of particles must be a multiple of ri_whfast512.N_systems.");
        r->status = REB_STATUS_GENERIC_ERROR;
        return;
    }
    if (r->G!=1.0){
        reb_simulation_error(r, "WHFast512 requires units in which G=1. Please rescale your system.");
        r->status = REB_STATUS_GENERIC_ERROR;
        return;
    }
    if (r->N_active!=-1 && r->N_active!=r->N){
        reb_simulation_error(r, "WHFast512 does not support test particles.");
        r->status = REB_STATUS_GENERIC_ERROR;
        return;
    }
    if (ri_whfast512->p512==NULL){
        ri_whfast512->p512 = aligned_alloc(64,sizeof(struct reb_particle_avx512));
        if (!ri_whfast512->p512){
            reb_simulation_error(r, "WHFast512 was not able to allocate memory.");
            r->status = REB_STATUS_GENERIC_ERROR;
            return;
        }
    }
    if (r->exit_min_distance || r->exit_max_distance){
        reb_simulation_warning(r, "You are using WHFast512 together with the flags exit_min_distance and/or exit_max_distance. With the current implementation, these flags will only check the last synchronized positions. In addition they might slow down WHFast512 significantly. If you need to use these flags, please open an issue on GitHub for further advice.");
    }
    ri_whfast512->N_allocated=1;
    r->gravity = REB_GRAVITY_NONE; // WHFast512 uses its own gravity routine.
}

// Main integration routine
void reb_integrator_whfast512_part1(struct reb_simulation* const r){
    struct reb_integrator_whfast512* const ri_whfast512 = &(r->ri_whfast512);
    const double dt = r->dt;

    if (ri_whfast512->N_allocated==0){
        reb_integrator_whfast512_allocate(r);
        ri_whfast512->recalculate_constants = 1;
    }
       
    if (ri_whfast512->recalculate_constants){
        recalculate_constants(r);
    }

    if (ri_whfast512->is_synchronized){
        ri_whfast512->time_since_last_synchronize = r->t;
        switch (ri_whfast512->coordinates){
            case REB_WHFAST512_COORDINATES_DEMOCRATICHELIOCENTRIC:
                inertial_to_democratic_heliocentric_posvel(r);
                break;
            case REB_WHFAST512_COORDINATES_JACOBI:
                inertial_to_jacobi_posvel(r);
                break;
            default:
                reb_simulation_error(r,"Coordinate system not supported.");
        }
    }

    if (ri_whfast512->is_synchronized){
        // First half DRIFT step
        reb_whfast512_kepler_step(r, dt/2.);    
    }else{
        // Combined DRIFT step
        reb_whfast512_kepler_step(r, dt);    // full timestep
    }

    if (ri_whfast512->coordinates==REB_WHFAST512_COORDINATES_DEMOCRATICHELIOCENTRIC){
        if (ri_whfast512->gr_potential){
            reb_whfast512_jump_step(r, dt/2.);
        }else{
            reb_whfast512_jump_step(r, dt);
        }
    }
    
    if (ri_whfast512->N_systems==1){
        if (ri_whfast512->coordinates==REB_WHFAST512_COORDINATES_DEMOCRATICHELIOCENTRIC){
            reb_whfast512_interaction_step_8planets_democraticheliocentric(r);
        }else{
            reb_whfast512_interaction_step_8planets_jacobi(r);
        }
    }else if (ri_whfast512->N_systems==2){
        if (ri_whfast512->coordinates==REB_WHFAST512_COORDINATES_DEMOCRATICHELIOCENTRIC){
            reb_whfast512_interaction_step_4planets_democraticheliocentric(r);
        }else{
            reb_whfast512_interaction_step_4planets_jacobi(r);
        }
    }else if (ri_whfast512->N_systems==4){
        if (ri_whfast512->coordinates==REB_WHFAST512_COORDINATES_DEMOCRATICHELIOCENTRIC){
            reb_whfast512_interaction_step_2planets_democraticheliocentric(r);
        }else{
            reb_whfast512_interaction_step_2planets_jacobi(r);
        }
    }

    if (ri_whfast512->coordinates==REB_WHFAST512_COORDINATES_DEMOCRATICHELIOCENTRIC){
        if (ri_whfast512->gr_potential){
            reb_whfast512_jump_step(r, dt/2.);
        }
    }

    ri_whfast512->is_synchronized = 0;

    r->t += dt;
    r->dt_last_done = dt;
}

#else // AVX512
      // Dummy function when AVX512 is not available
void reb_integrator_whfast512_part1(struct reb_simulation* const r){
    reb_simulation_error(r, "WHFast512 is not available. Please make sure your CPU supports AVX512 instructions, then recompile REBOUND with the AVX512 option turned on in the Makefile or set the AVX512 environment variable to 1 before running pip install.");
    r->status = REB_STATUS_GENERIC_ERROR;
}
#endif // AVX512

// Synchronization routine. Called every time an output is needed.
void reb_integrator_whfast512_synchronize(struct reb_simulation* const r){
#ifdef AVX512
    struct reb_integrator_whfast512* const ri_whfast512 = &(r->ri_whfast512);
    if (ri_whfast512->is_synchronized == 0){
        const unsigned int N_systems = ri_whfast512->N_systems;
        struct reb_particle_avx512* sync_pj = NULL;
        struct reb_particle sync_pj0[4];
        // Needed if no step has ever been done before (like SA)
        if (ri_whfast512->N_allocated==0){
            reb_integrator_whfast512_allocate(r);
            recalculate_constants(r);
        }
        if (ri_whfast512->keep_unsynchronized){
            sync_pj = aligned_alloc(64,sizeof(struct reb_particle_avx512));
            memcpy(sync_pj,ri_whfast512->p512, sizeof(struct reb_particle_avx512));
            for (int s=0; s<N_systems; s++){
                sync_pj0[s] = ri_whfast512->p_jh0[s];
            }
        }
        reb_whfast512_kepler_step(r, r->dt/2.);    
        reb_whfast512_com_step(r, r->t-ri_whfast512->time_since_last_synchronize);
        switch (ri_whfast512->coordinates){
            case REB_WHFAST512_COORDINATES_DEMOCRATICHELIOCENTRIC:
                democraticheliocentric_to_inertial_posvel(r);
                break;
            case REB_WHFAST512_COORDINATES_JACOBI:
                jacobi_to_inertial_posvel(r);
                break;
            default:
                reb_simulation_error(r,"Coordinate system not supported.");
        }

        if (ri_whfast512->keep_unsynchronized){
            memcpy(ri_whfast512->p512, sync_pj, sizeof(struct reb_particle_avx512));
            for (int s=0; s<N_systems; s++){
                ri_whfast512->p_jh0[s] = sync_pj0[s];
            }
            free(sync_pj);
        }else{
            ri_whfast512->is_synchronized = 1;
            ri_whfast512->time_since_last_synchronize = r->t;
        }
    }
#else 
    reb_integrator_whfast512_synchronize_fallback(r);
#endif // AVX512
}

void reb_integrator_whfast512_synchronize_fallback(struct reb_simulation* const r){
    // No AVX512 available
    // Using WHFast as a workaround.
    // Not bit-wise reproducible. 
    struct reb_integrator_whfast512* const ri_whfast512 = &(r->ri_whfast512);
    if (ri_whfast512->is_synchronized == 0){
        reb_simulation_warning(r, "WHFast512 is not available. Synchronization is provided using WHFast and is not bit-compatible to WHFast512.");
        const unsigned int N_systems = ri_whfast512->N_systems;
        const unsigned int p_per_system = 8/N_systems;
        const unsigned int N_per_system = r->N/N_systems;
        double dt = r->dt;
        for (int s=0; s<N_systems; s++){
            double m0 = r->particles[s*N_per_system].m;
            // 1/2 Kepler
            for (unsigned int i=1;i<N_per_system;i++){
                struct reb_particle p = {0};
                p.m = ri_whfast512->p512->m[s*p_per_system+i-1];
                p.x = ri_whfast512->p512->x[s*p_per_system+i-1];
                p.y = ri_whfast512->p512->y[s*p_per_system+i-1];
                p.z = ri_whfast512->p512->z[s*p_per_system+i-1];
                p.vx = ri_whfast512->p512->vx[s*p_per_system+i-1];
                p.vy = ri_whfast512->p512->vy[s*p_per_system+i-1];
                p.vz = ri_whfast512->p512->vz[s*p_per_system+i-1];
                reb_whfast_kepler_solver(r, &p, m0, 0, dt/2.0);
                ri_whfast512->p512->x[s*p_per_system+i-1]  = p.x;
                ri_whfast512->p512->y[s*p_per_system+i-1]  = p.y;
                ri_whfast512->p512->z[s*p_per_system+i-1]  = p.z;
                ri_whfast512->p512->vx[s*p_per_system+i-1] = p.vx;
                ri_whfast512->p512->vy[s*p_per_system+i-1] = p.vy;
                ri_whfast512->p512->vz[s*p_per_system+i-1] = p.vz;
            }
        }
        reb_whfast512_com_step(r, r->t-ri_whfast512->time_since_last_synchronize); // does not use AVX512
        democraticheliocentric_to_inertial_posvel(r);
        ri_whfast512->is_synchronized = 1;
    }
}

// Free memory and reset all constants.
// This needs to be called when the timestep, the number of particles, masses, etc are changed, 
void reb_integrator_whfast512_reset(struct reb_simulation* const r){
    struct reb_integrator_whfast512* const ri_whfast512 = &(r->ri_whfast512);
    if (ri_whfast512->N_allocated){
        free(ri_whfast512->p512);
    }
    free(ri_whfast512->mat8_inertial_to_jacobi);
    free(ri_whfast512->mat8_jacobi_to_inertial);
    free(ri_whfast512->mat8_jacobi_to_heliocentric);
    ri_whfast512->mat8_inertial_to_jacobi = NULL;
    ri_whfast512->mat8_jacobi_to_inertial = NULL;
    ri_whfast512->mat8_jacobi_to_heliocentric = NULL;
    ri_whfast512->p512 = NULL;
    ri_whfast512->N_allocated = 0;
    ri_whfast512->gr_potential = 0;
    ri_whfast512->is_synchronized = 1;
    ri_whfast512->keep_unsynchronized = 0;
}

// Everything is in part 1 for this integrator
void reb_integrator_whfast512_part2(struct reb_simulation* const r){
}
