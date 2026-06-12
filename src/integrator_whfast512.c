/**
 * integrator_whfast512.c: ASM version of WHFast512
 * 
 * Copyright (c) 2025 Hanno Rein
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
#include <sys/ioctl.h>
#include <string.h>
#if defined(__i386__) || defined(__x86_64__)
#include <immintrin.h>
#pragma GCC target("avx512f,avx512dq,avx512bw,avx512cd,avx512vl")
#else
typedef struct {
        double lanes[8];
} __attribute__((aligned(64))) __m512d;
typedef struct {
        uint64_t lanes[8];
} __attribute__((aligned(64))) __m512i;
typedef char  __attribute__((aligned(64))) __mmask8; 
#endif
#include "rebound.h"
#include "rebound_internal.h"
#include "particle.h"
#include "tools.h"
#include "gravity.h"
#include "boundary.h"
#include "integrator_whfast.h"
#include "integrator_whfast512.h"


#define DEBUG_AVX512 1

void reb_integrator_whfast512_free(void* state);		
void* reb_integrator_whfast512_create();		
void reb_integrator_whfast512_step(struct reb_simulation* r, void* state);
void reb_integrator_whfast512_synchronize(struct reb_simulation* r, void* state);
const struct reb_binarydata_field_descriptor reb_integrator_whfast512_field_descriptor_list[];

#define SIMD_DATA_MEMBERS X(M) X(dt) X(gr_prefac) X(m) X(x) X(y) X(z) X(vx) X(vy) X(vz) \
    X(mat8_inertial_to_jacobi) \
    X(mat8_jacobi_to_heliocentric) \
    X(M0) X(mask) \
    X(mat8_jacobi_to_inertial)\
    X(counter) 

struct simd_data{
    // Various constants
    __m512d M __attribute__ ((aligned (64)));                   //  Masses used in Kepler-Solver
    __m512d dt __attribute__ ((aligned (64)));                  //  Timestep
    __m512d gr_prefac __attribute__ ((aligned (64)));           //  Prefactor for GR
    __m512d m __attribute__ ((aligned (64)));
    __m512d x __attribute__ ((aligned (64)));
    __m512d y __attribute__ ((aligned (64)));
    __m512d z __attribute__ ((aligned (64)));
    __m512d vx __attribute__ ((aligned (64)));
    __m512d vy __attribute__ ((aligned (64)));
    __m512d vz __attribute__ ((aligned (64)));
    double mat8_inertial_to_jacobi[64] __attribute__ ((aligned (64))); // Coordinate transformation matricies. Can be recalculated from particle masses.
    double mat8_jacobi_to_heliocentric[64] __attribute__ ((aligned (64)));
    __m512d M0 __attribute__ ((aligned (64)));                   //  Masses used in Jacobi Term
    // Mask for cases with less than 8 planets
    __mmask8 mask __attribute__ ((aligned (64)));
    double mat8_jacobi_to_inertial[64] __attribute__ ((aligned (64)));
#ifdef DEBUG_AVX512
    __m512i counter __attribute__ ((aligned (64)));
#endif // DEBUG_AVX512
};

const struct reb_integrator reb_integrator_whfast512 = {
    .step = reb_integrator_whfast512_step,
    .create = reb_integrator_whfast512_create,
    .free = reb_integrator_whfast512_free,
    .synchronize = reb_integrator_whfast512_synchronize,
    .field_descriptor_list = reb_integrator_whfast512_field_descriptor_list,
};

const struct reb_binarydata_field_descriptor reb_integrator_whfast512_field_descriptor_list[] = {
    { "", REB_UINT,        "gr_potential",    offsetof(struct reb_integrator_whfast512_state, gr_potential), 0, 0, 0},
    { "", REB_UINT,        "corrector",       offsetof(struct reb_integrator_whfast512_state, corrector), 0, 0, 0},
    { "", REB_UINT,        "concatenate_steps", offsetof(struct reb_integrator_whfast512_state, concatenate_steps), 0, 0, 0},
    { "", REB_UINT,        "N_systems",       offsetof(struct reb_integrator_whfast512_state, N_systems), 0, 0, 0},
    { "", REB_POINTER_ALIGNED, "data",        offsetof(struct reb_integrator_whfast512_state, data), offsetof(struct reb_integrator_whfast512_state, N_allocated), sizeof(struct simd_data), 0},
    { "", REB_DOUBLE,      "last_synchronization", offsetof(struct reb_integrator_whfast512_state, last_synchronization), 0, 0, 0},
    { 0 }, // Null terminated list
};



void* reb_integrator_whfast512_create(){
    struct reb_integrator_whfast512_state* whfast512 = calloc(sizeof(struct reb_integrator_whfast512_state),1);
    whfast512->N_systems = 1;
    whfast512->gr_potential = 0;
    whfast512->concatenate_steps = 1;
    return whfast512;
}

void reb_integrator_whfast512_free(void* state){
    struct reb_integrator_whfast512_state* whfast512 = state;
    free(whfast512->data);
    free(whfast512);
}


#if defined(__i386__) || defined(__x86_64__)
#ifdef DEBUG_AVX512
uint64_t reb_whfast512_counter(struct reb_simulation* r, int test_p){
    struct reb_integrator_whfast512_state* whfast512 = r->integrator.state;
    struct simd_data* data = whfast512->data;
    uint64_t i[8];
    _mm512_store_epi64(&i[0], data->counter);
    return i[test_p];
}

// Debug function to print vectors
static inline void printavx512(__m512d a) {
    double _nax[8];
    _mm512_storeu_pd(&_nax[0], a);
    printf("avx = {%.17g, %.17g, %.17g, %.17g, %.17g, %.17g, %.17g, %.17g}\n", _nax[0], _nax[1], _nax[2], _nax[3], _nax[4], _nax[5], _nax[6], _nax[7]);
}

// Print mask in binary format for debuggin
static inline void printmask8(__mmask8 mask) {
    for (int i = 7; i >= 0; i--) {
        printf("%d", (mask >> i) & 1);
        if (i == 4) printf(" "); 
    }
    printf("\n");
}

// Debug function to print 8x8 matrix
static inline void printmat8(double* a) {
    for (int i=0; i<8;i++){
        for (int j=0; j<8;j++){
            printf("%.16f ", a[i*8+j]);
        }
        printf("\n");
    }
}
#endif // DEBUG_AVX512

// 8x8 matrix multiplication using avx512
__attribute__((target("avx512f,avx512vl,avx512bw,avx512dq")))
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
__attribute__((target("avx512f,avx512vl,avx512bw,avx512dq")))
void mat8_mul3_avx512(const double* matrix, const __m512d in1, const __m512d in2, const __m512d in3, __m512d* out1, __m512d* out2, __m512d* out3){
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
__attribute__((target("avx512f,avx512vl,avx512bw,avx512dq")))
static __m512d load_into_m512d(struct reb_simulation* r, size_t offset, const double* transformation, int N_systems){
    struct reb_particle* particles = r->particles;
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
__attribute__((target("avx512f,avx512vl,avx512bw,avx512dq")))
static void load_from_m512d(struct reb_simulation* r, size_t offset, const double* transformation, int N_systems, __m512d vector){
    struct reb_particle* particles = r->particles;
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

// Convert jacobi coordinates to inertial coordinates
// Also performs com step (assume original particles are unmodified)
// Note: Speed is not a concern here 
__attribute__((target("avx512f,avx512vl,avx512bw,avx512dq")))
static void jacobi_to_inertial_posvel_and_com(struct reb_simulation* r, struct simd_data* data, double dt_com, unsigned int N_systems){
    const unsigned int N_per_system = r->N/N_systems;
    struct reb_particle com[4];
    for (unsigned s=0;s<N_systems;s++){
        com[s] = reb_simulation_com_range(r,s*N_per_system, (s+1)*N_per_system); // original com
    }
    struct reb_particle* particles = r->particles;
    load_from_m512d(r, offsetof(struct reb_particle, x), data->mat8_jacobi_to_inertial, N_systems, data->x);
    load_from_m512d(r, offsetof(struct reb_particle, y), data->mat8_jacobi_to_inertial, N_systems, data->y);
    load_from_m512d(r, offsetof(struct reb_particle, z), data->mat8_jacobi_to_inertial, N_systems, data->z);
    load_from_m512d(r, offsetof(struct reb_particle, vx), data->mat8_jacobi_to_inertial, N_systems, data->vx);
    load_from_m512d(r, offsetof(struct reb_particle, vy), data->mat8_jacobi_to_inertial, N_systems, data->vy);
    load_from_m512d(r, offsetof(struct reb_particle, vz), data->mat8_jacobi_to_inertial, N_systems, data->vz);
    for (unsigned s=0;s<N_systems;s++){
        particles[s*N_per_system+0].x  = 0.0;
        particles[s*N_per_system+0].y  = 0.0;
        particles[s*N_per_system+0].z  = 0.0;
        particles[s*N_per_system+0].vx = 0.0;
        particles[s*N_per_system+0].vy = 0.0;
        particles[s*N_per_system+0].vz = 0.0;
        for (unsigned int i=1;i<N_per_system;i++){
            particles[s*N_per_system+0].x  -= particles[s*N_per_system+i].m/particles[s*N_per_system+0].m*particles[s*N_per_system+i].x;
            particles[s*N_per_system+0].y  -= particles[s*N_per_system+i].m/particles[s*N_per_system+0].m*particles[s*N_per_system+i].y;
            particles[s*N_per_system+0].z  -= particles[s*N_per_system+i].m/particles[s*N_per_system+0].m*particles[s*N_per_system+i].z;
            particles[s*N_per_system+0].vx -= particles[s*N_per_system+i].m/particles[s*N_per_system+0].m*particles[s*N_per_system+i].vx;
            particles[s*N_per_system+0].vy -= particles[s*N_per_system+i].m/particles[s*N_per_system+0].m*particles[s*N_per_system+i].vy;
            particles[s*N_per_system+0].vz -= particles[s*N_per_system+i].m/particles[s*N_per_system+0].m*particles[s*N_per_system+i].vz;
        }
        for (unsigned int i=0;i<N_per_system;i++){
            particles[s*N_per_system+i].x  += com[s].x + com[s].vx*dt_com;
            particles[s*N_per_system+i].y  += com[s].y + com[s].vy*dt_com;
            particles[s*N_per_system+i].z  += com[s].z + com[s].vz*dt_com;
            particles[s*N_per_system+i].vx += com[s].vx;
            particles[s*N_per_system+i].vy += com[s].vy;
            particles[s*N_per_system+i].vz += com[s].vz;
        }
    }
}


extern void reb_whfast512_full_steps_democraticheliocentric_gr(struct simd_data* data, long N_steps, int skip_first_kepler_step);
extern void reb_whfast512_full_steps_democraticheliocentric_nogr(struct simd_data* data, long N_steps, int skip_first_kepler_step);
extern void reb_whfast512_full_steps_jacobi_gr(struct simd_data* data, long N_steps, int skip_first_kepler_step);
extern void reb_whfast512_full_steps_jacobi_nogr(struct simd_data* data, long N_steps, int skip_first_kepler_step);
extern void reb_whfast512_corrector_step_gr(struct simd_data* data, double inv);
extern void reb_whfast512_corrector_step_nogr(struct simd_data* data, double inv);
extern void reb_whfast512_kepler_step(struct simd_data* data);

// _n2 = two systems of up to 4 planets, _n4 = four systems of 2 planets.
extern void reb_whfast512_full_steps_jacobi_gr_n2(struct simd_data* data, long N_steps, int skip_first_kepler_step);
extern void reb_whfast512_full_steps_jacobi_nogr_n2(struct simd_data* data, long N_steps, int skip_first_kepler_step);
extern void reb_whfast512_full_steps_jacobi_gr_n4(struct simd_data* data, long N_steps, int skip_first_kepler_step);
extern void reb_whfast512_full_steps_jacobi_nogr_n4(struct simd_data* data, long N_steps, int skip_first_kepler_step);
extern void reb_whfast512_corrector_step_gr_n2(struct simd_data* data, double inv);
extern void reb_whfast512_corrector_step_nogr_n2(struct simd_data* data, double inv);
extern void reb_whfast512_corrector_step_gr_n4(struct simd_data* data, double inv);
extern void reb_whfast512_corrector_step_nogr_n4(struct simd_data* data, double inv);

static void whfast512_full_steps(struct reb_integrator_whfast512_state* whfast512, long N_steps, int skip){
    struct simd_data* data = whfast512->data;
    if (whfast512->gr_potential){
        switch (whfast512->N_systems){
            case 2: reb_whfast512_full_steps_jacobi_gr_n2(data, N_steps, skip); break;
            case 4: reb_whfast512_full_steps_jacobi_gr_n4(data, N_steps, skip); break;
            default: reb_whfast512_full_steps_jacobi_gr(data, N_steps, skip); break;
        }
    }else{
        switch (whfast512->N_systems){
            case 2: reb_whfast512_full_steps_jacobi_nogr_n2(data, N_steps, skip); break;
            case 4: reb_whfast512_full_steps_jacobi_nogr_n4(data, N_steps, skip); break;
            default: reb_whfast512_full_steps_jacobi_nogr(data, N_steps, skip); break;
        }
    }
}

static void whfast512_corrector_step(struct reb_integrator_whfast512_state* whfast512, double inv){
    struct simd_data* data = whfast512->data;
    if (whfast512->gr_potential){
        switch (whfast512->N_systems){
            case 2: reb_whfast512_corrector_step_gr_n2(data, inv); break;
            case 4: reb_whfast512_corrector_step_gr_n4(data, inv); break;
            default: reb_whfast512_corrector_step_gr(data, inv); break;
        }
    }else{
        switch (whfast512->N_systems){
            case 2: reb_whfast512_corrector_step_nogr_n2(data, inv); break;
            case 4: reb_whfast512_corrector_step_nogr_n4(data, inv); break;
            default: reb_whfast512_corrector_step_nogr(data, inv); break;
        }
    }
}

__attribute__((target("avx512f,avx512vl,avx512bw,avx512dq")))
static void inertial_to_jacobi_posvel(struct reb_simulation* r, struct simd_data* data, unsigned int N_systems){
    const unsigned int N_per_system = r->N/N_systems;
    // Transformations assume system is in COM frame.
    struct reb_particle com[4];
    struct reb_particle* p_tmp = malloc(sizeof(struct reb_particle)*r->N);
    memcpy(p_tmp, r->particles, sizeof(struct reb_particle)*r->N);
    for (unsigned int s=0;s<N_systems;s++){
        com[s] = reb_simulation_com_range(r,s*N_per_system, (s+1)*N_per_system); // original com
    }
    for (unsigned int s=0;s<N_systems;s++){
        for (unsigned int i=0;i<N_per_system;i++){
            r->particles[s*N_per_system+i].x -= com[s].x;
            r->particles[s*N_per_system+i].y -= com[s].y;
            r->particles[s*N_per_system+i].z -= com[s].z;
            r->particles[s*N_per_system+i].vx -= com[s].vx;
            r->particles[s*N_per_system+i].vy -= com[s].vy;
            r->particles[s*N_per_system+i].vz -= com[s].vz;
        }
    }
    reb_simulation_move_to_com(r);
    // Same layout as for democratic heliocentric
    data->x = load_into_m512d(r, offsetof(struct reb_particle,x),data->mat8_inertial_to_jacobi, N_systems);
    data->y = load_into_m512d(r, offsetof(struct reb_particle,y),data->mat8_inertial_to_jacobi, N_systems);
    data->z = load_into_m512d(r, offsetof(struct reb_particle,z),data->mat8_inertial_to_jacobi, N_systems);
    data->vx = load_into_m512d(r, offsetof(struct reb_particle,vx),data->mat8_inertial_to_jacobi, N_systems);
    data->vy = load_into_m512d(r, offsetof(struct reb_particle,vy),data->mat8_inertial_to_jacobi, N_systems);
    data->vz = load_into_m512d(r, offsetof(struct reb_particle,vz),data->mat8_inertial_to_jacobi, N_systems);
    data->m = load_into_m512d(r, offsetof(struct reb_particle,m),NULL, N_systems);
    // Undo COM transformation. COM will be applied in jacobi_to_inertial_posvel_and_com().
    memcpy(r->particles, p_tmp, sizeof(struct reb_particle)*r->N);
    free(p_tmp);
}


// Precalculate various constants and put them in 512 bit vectors.
__attribute__((target("avx512f,avx512vl,avx512bw,avx512dq")))
static void recalculate_constants(struct reb_simulation* r, unsigned int N_systems){
    struct reb_integrator_whfast512_state* whfast512 = r->integrator.state;
    free(whfast512->data); // free in case previously allocated
    whfast512->data = aligned_alloc(64,sizeof(struct simd_data));
    whfast512->N_allocated=1;
    memset(whfast512->data, 0, sizeof(struct simd_data));
    if (!whfast512->data){
        reb_simulation_error(r, "WHFast512 was not able to allocate memory.");
        return;
    }
    struct simd_data* data = whfast512->data;
    struct reb_particle* particles = r->particles;
    const unsigned int p_per_system = 8/N_systems;
    const unsigned int N_per_system = r->N/N_systems;
    double mat8_inertial_to_heliocentric[64];
    double M[8] = {0.0};
    double M0[8] = {0.0};
    switch (N_systems){
        case 1:
            data->mask = (1 << (r->N -1)) - 1;
            break;
        case 2:
            if (N_per_system==5) data->mask = 0xFF;
            if (N_per_system==4) data->mask = 0x77;
            if (N_per_system==3) data->mask = 0x33;
            if (N_per_system==2) data->mask = 0x11;
            break;
        case 4:
            if (N_per_system==3) data->mask = 0xFF;
            if (N_per_system==2) data->mask = 0x55;
            break;
        default:
            reb_simulation_error(r,"Invalid value for N_systems.");
    }
    // Zeroing.
    for (unsigned int i=1; i<9; i++){
        for (unsigned int j=1; j<9; j++){
            data->mat8_inertial_to_jacobi[(i-1)+8*(j-1)] = 0.0;
            data->mat8_jacobi_to_inertial[(i-1)+8*(j-1)] = 0.0;
            data->mat8_jacobi_to_heliocentric[(i-1)+8*(j-1)] = 0.0;
            mat8_inertial_to_heliocentric[(i-1)+8*(j-1)] = 0.0;
        }
    }

    // Filling vector
    for (unsigned int s=0; s<N_systems; s++){
        for (unsigned int i=0;i<N_per_system-1;i++){
            for (unsigned int j=0;j<i+2;j++){
                M[s*p_per_system+i] += r->particles[s*N_per_system+j].m;
            }
        }
    }

    // Fill matricies
    for (unsigned int s=0; s<N_systems; s++){
        double ms = particles[s*N_per_system+0].m;
        for (unsigned int i=1; i<N_per_system; i++){
            for (unsigned int j=i; j<N_per_system; j++){
                data->mat8_inertial_to_jacobi[(s*p_per_system+i-1)+8*(s*p_per_system+j-1)] += particles[s*N_per_system+j].m/ms;
            }
            for (unsigned int j=1; j<N_per_system; j++){
                mat8_inertial_to_heliocentric[(s*p_per_system+i-1)+8*(s*p_per_system+j-1)] += particles[s*N_per_system+j].m/particles[s*N_per_system+0].m;
            }
            mat8_inertial_to_heliocentric[(s*p_per_system+i-1)+8*(s*p_per_system+i-1)] += 1.0;
            data->mat8_inertial_to_jacobi[(s*p_per_system+i-1)+8*(s*p_per_system+i-1)] += 1.0;
            data->mat8_jacobi_to_inertial[(s*p_per_system+i-1)+8*(s*p_per_system+i-1)] += ms/(ms + particles[s*N_per_system+i].m);
            ms += particles[s*N_per_system+i].m;
            if (i<N_per_system-1){ 
                for (unsigned int ii=i; ii>0; ii--){
                    int jj = i+1;
                    data->mat8_jacobi_to_inertial[(s*p_per_system+ii-1)+8*(s*p_per_system+jj-1)] -= particles[s*N_per_system+jj].m/(ms+particles[s*N_per_system+jj].m);
                }
            }
            M0[(s*p_per_system+i-1)] = particles[s*N_per_system+0].m;
        }
    }

    // Might be numerically more stable to calculate this manually rather than do a matrix multiplication.
    for (unsigned int i=1; i<9; i++){
        for (unsigned int j=1; j<9; j++){
            for (unsigned int k=1; k<9; k++){
                data->mat8_jacobi_to_heliocentric[(i-1)+8*(j-1)] += mat8_inertial_to_heliocentric[(i-1)+8*(k-1)] * data->mat8_jacobi_to_inertial[(k-1)+8*(j-1)];
            }
        }
    }

    data->M = _mm512_loadu_pd(&M);
    data->M0 = _mm512_loadu_pd(&M0); //  = particles[0].m 

    // GR prefactors. Note: assumes units of AU, year/2pi.
    double c = 10065.32;
    double _gr_prefac[8];
    for(unsigned int i=0;i<8;i++){
        _gr_prefac[i] = 0; // for when N<8
    }
    for (unsigned int s=0; s<N_systems; s++){
        double m0 = r->particles[s*N_per_system].m;
        for (unsigned int p=1; p<N_per_system; p++){
            _gr_prefac[s*p_per_system+(p-1)] = -6.*m0*m0/(c*c);
        }
    }
    data->gr_prefac = _mm512_loadu_pd(&_gr_prefac);
    data->dt = _mm512_set1_pd(r->dt); 
#define X(name) printf(".set P512_" #name ", %zu\n", offsetof(struct simd_data, name));
//    SIMD_DATA_MEMBERS
#undef X

}

__attribute__((target("avx512f,avx512vl,avx512bw,avx512dq")))
static int reb_integrator_whfast512_verify_setup(struct reb_simulation* const r){
    struct reb_integrator_whfast512_state* whfast512 = r->integrator.state;
    // Check if all assumptions are satisfied.
    // Note: These are not checked every timestep. 
    // So it is possible for the user to screw things up.
//    if (r->dt<0.0){
//        reb_simulation_error(r, "WHFast512 does not support negative timesteps. To integrate backwards, flip the sign of the velocities.");
//        return 1;
//    }
    if (!reb_avx512_available()){
        reb_simulation_error(r, "AVX512 is not supported by your CPU.");
        return 1;
    }
    if (r->N_var!=0){
        reb_simulation_error(r, "WHFast512 does not support variational particles.");
        return 1;
    }
    if (r->exact_finish_time!=0){
        reb_simulation_error(r, "WHFast512 requires exact_finish_time=0.");
        return 1;
    }
    if (r->N>9 && whfast512->N_systems == 1) {
        reb_simulation_error(r, "WHFast512 supports a maximum of 9 particles when N_systems is set to 1.");
        return 1;
    }
    if (r->N>10 && whfast512->N_systems == 2) {
        reb_simulation_error(r, "WHFast512 supports a maximum of 10 particles when N_systems is set to 2.");
        return 1;
    }
    if (r->N>12 && whfast512->N_systems == 4) {
        reb_simulation_error(r, "WHFast512 supports a maximum of 12 particles when N_systems is set to 4.");
        return 1;
    }
    if (whfast512->N_systems != 1 && whfast512->N_systems !=2 && whfast512->N_systems != 4){
        reb_simulation_error(r, "WHFast512 supports 1, 2, or 4 systems only.");
        return 1;
    }
    if (r->N % whfast512->N_systems != 0){
        reb_simulation_error(r, "Number of particles must be a multiple of whfast512.N_systems.");
        return 1;
    }
    if (r->G!=1.0){
        reb_simulation_error(r, "WHFast512 requires units in which G=1. Please rescale your system.");
        return 1;
    }
    if (r->N_active!=SIZE_MAX && r->N_active!=r->N){
        reb_simulation_error(r, "WHFast512 does not support test particles.");
        return 1;
    }
    if (r->exit_min_distance || r->exit_max_distance){
        reb_simulation_warning(r, "You are using WHFast512 together with the flags exit_min_distance and/or exit_max_distance. With the current implementation, these flags will only check the last synchronized positions. In addition they might slow down WHFast512 significantly. If you need to use these flags, please open an issue on GitHub for further advice.");
    }
    r->gravity = REB_GRAVITY_NONE; // WHFast512 uses its own gravity routine.
    return 0; // success
}

__attribute__((target("avx512f,avx512vl,avx512bw,avx512dq")))
void reb_integrator_whfast512_kepler_step(struct reb_simulation* const r, int N_steps){
    struct reb_integrator_whfast512_state* whfast512 = r->integrator.state;
    recalculate_constants(r, whfast512->N_systems);
    inertial_to_jacobi_posvel(r, whfast512->data, whfast512->N_systems);
    struct simd_data* data = whfast512->data;
    for (int i=0; i<N_steps; i++){
        reb_whfast512_kepler_step(data);    
    }
    jacobi_to_inertial_posvel_and_com(r, whfast512->data, 0.0, whfast512->N_systems);
    // This will leak memory. Cannot free because of counter access
    //free(whfast512->data);
    //whfast512->data = NULL;
}

// Optimized main loops allowing for concatenate_steps
__attribute__((target("avx512f,avx512vl,avx512bw,avx512dq")))
void reb_integrator_whfast512_step(struct reb_simulation* const r, void* state){
    struct reb_integrator_whfast512_state* whfast512 = state;
    const double dt = r->dt;
    const unsigned int N_steps = whfast512->concatenate_steps;

    if (reb_integrator_whfast512_verify_setup(r)){
        r->status = REB_STATUS_GENERIC_ERROR;
        return; // Error occured
    }

    int skip_first_kepler_step = 0;
    if (r->is_synchronized){
        recalculate_constants(r, whfast512->N_systems);
        struct simd_data* data = whfast512->data;
        // Use WHFast to apply the correctors.
        inertial_to_jacobi_posvel(r, data, whfast512->N_systems);
        whfast512->last_synchronization = r->t;
        if (whfast512->corrector){
            whfast512_corrector_step(whfast512, 1.0);
        }
        // First half DRIFT step.
        skip_first_kepler_step = 1;
        data->dt = _mm512_set1_pd(dt/2.0); 
        reb_whfast512_kepler_step(data);    
        data->dt = _mm512_set1_pd(dt); // Reset
    }

    whfast512_full_steps(whfast512, N_steps, skip_first_kepler_step);

    r->is_synchronized = 0;
    r->t += dt*N_steps;
    r->dt_last_done = dt;
}

// Synchronization routine. Called every time an output is needed.
__attribute__((target("avx512f,avx512vl,avx512bw,avx512dq")))
void reb_integrator_whfast512_synchronize(struct reb_simulation* const r, void* state){
    if (!reb_avx512_available()){
        reb_simulation_error(r, "AVX512 is not supported by your CPU.");
        return;
    }
    struct reb_integrator_whfast512_state* const whfast512 = state;
    if (!r->is_synchronized){
        struct simd_data * data = whfast512->data;
        if (!data){
            reb_simulation_error(r, "ASM512 is unable to synchronize. data is NULL.");
            return;
        }
        data->dt = _mm512_set1_pd(r->dt/2.0); 
        reb_whfast512_kepler_step(data);    
        data->dt = _mm512_set1_pd(r->dt); // Reset
                                         // TODO Add COM step
        if (whfast512->corrector){
            whfast512_corrector_step(whfast512, -1.0);
        }
        double dt_com = r->t - whfast512->last_synchronization;
        jacobi_to_inertial_posvel_and_com(r, data, dt_com, whfast512->N_systems);
        whfast512->last_synchronization = r->t;
        r->is_synchronized = 1;
        free(whfast512->data);
        whfast512->data = NULL;
    }
}

#else // defined(__i386__) || defined(__x86_64__)
void reb_integrator_whfast512_step(struct reb_simulation* r, void* state){
    (void)state;
    reb_simulation_error(r, "AVX512 is not supported on your platform");
}
void reb_integrator_whfast512_synchronize(struct reb_simulation* r, void* state){
    (void)state;
    reb_simulation_error(r, "AVX512 is not supported on your platform");
}
#endif // defined(__i386__) || defined(__x86_64__)
