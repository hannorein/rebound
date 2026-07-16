/**
 * integrator_whfast512.c: ASM version of WHFast512 in Jacobi Coordinates
 * 
 * Copyright (c) 2026 Rishit Dagli, Hanno Rein, Pejvak Javaheri
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
#include "rebound_internal.h"
#if !defined(_WIN32)
#define REB_ATTRIBUTE_ALIGNED_64 __attribute__((aligned(64)))
#else
#define REB_ATTRIBUTE_ALIGNED_64
#endif
#if (defined(__i386__) || defined(__x86_64__)) && !defined(_WIN32)
#include <immintrin.h>
#pragma GCC target("avx512f,avx512dq,avx512bw,avx512cd,avx512vl")
#else
typedef struct {
    double lanes[8];
} REB_ATTRIBUTE_ALIGNED_64 __m512d;
typedef struct {
    uint64_t lanes[8];
} REB_ATTRIBUTE_ALIGNED_64 __m512i;
typedef char  REB_ATTRIBUTE_ALIGNED_64 __mmask8; 
#endif
#include "particle.h"
#include "tools.h"
#include "gravity.h"
#include "boundary.h"
#include "integrator_whfast.h"
#include "integrator_whfast512.h"

//#define DEBUG_AVX512 1

void reb_integrator_whfast512_free(void* state);		
void* reb_integrator_whfast512_create();		
void reb_integrator_whfast512_step(struct reb_simulation* r, void* state);
void reb_integrator_whfast512_synchronize(struct reb_simulation* r, void* state);
const struct reb_binarydata_field_descriptor reb_integrator_whfast512_field_descriptor_list[];

// The main datasctructure. We pass this as a pointer to the assembly code.
struct simd_data{
    // Various constants
    __m512d M REB_ATTRIBUTE_ALIGNED_64;             //  Masses used in Kepler-Solver
    __m512d dt REB_ATTRIBUTE_ALIGNED_64;            //  Timestep
    __m512d gr_prefac REB_ATTRIBUTE_ALIGNED_64;     //  Prefactor for GR
    __m512d m REB_ATTRIBUTE_ALIGNED_64;
    __m512d x REB_ATTRIBUTE_ALIGNED_64;
    __m512d y REB_ATTRIBUTE_ALIGNED_64;
    __m512d z REB_ATTRIBUTE_ALIGNED_64;
    __m512d vx REB_ATTRIBUTE_ALIGNED_64;
    __m512d vy REB_ATTRIBUTE_ALIGNED_64;
    __m512d vz REB_ATTRIBUTE_ALIGNED_64;
    double mat8_inertial_to_jacobi[64] REB_ATTRIBUTE_ALIGNED_64; // Coordinate transformation matricies. Can be recalculated from particle masses.
    double mat8_jacobi_to_heliocentric[64] REB_ATTRIBUTE_ALIGNED_64;
    __m512d M0 REB_ATTRIBUTE_ALIGNED_64;            //  Masses used in Jacobi Term
    __mmask8 mask REB_ATTRIBUTE_ALIGNED_64;         // Mask for cases with less than 8 planets
    __m512d exit_max_distance REB_ATTRIBUTE_ALIGNED_64;
    __m512d exit_min_distance_r REB_ATTRIBUTE_ALIGNED_64; // Inverse of exit_min_distance
    double mat8_jacobi_to_inertial[64] REB_ATTRIBUTE_ALIGNED_64;
#ifdef DEBUG_AVX512
    __m512i counter REB_ATTRIBUTE_ALIGNED_64;
#endif // DEBUG_AVX512
};

const struct reb_integrator reb_integrator_whfast512 = {
    .documentation =
    "WHFast512 is a highly optimized implementation of the symplectic [Wisdom & Holman (1991)] integrator. " 
    "It supports simulations with up to 9 particles (8 planets + 1 central object). " 
    "Note that by default WHFast512 combines 1e6 timesteps to improve speed. "
    "You can set `concatenate_steps` to a lower value for more fine grained output. "
    "WHFast512 uses Jacobi coordinates and support symplectic correctors as well as general relativistic corrections. "
    "\n\n"
    "The algorithm is described in two papers [Javaheri et al. (2023)] and [Dagli & Rein (in prep)]. "
    "Note that in July 2026 significant changes have been made. "
    "See [Dagli & Rein (in prep)] for details on WHFast512 version 2. "
    "\n\n"
    "[Wisdom & Holman (1991)]: https://ui.adsabs.harvard.edu/abs/1991AJ....102.1528W/abstract\n"
    "[Javaheri et al. (2023)]: https://ui.adsabs.harvard.edu/abs/2023OJAp....6E..29J/abstract\n"
    "[Dagli & Rein (in prep)]: https://hanno-rein.de/publications.html\n"
    ,
    .step = reb_integrator_whfast512_step,
    .create = reb_integrator_whfast512_create,
    .free = reb_integrator_whfast512_free,
    .synchronize = reb_integrator_whfast512_synchronize,
    .field_descriptor_list = reb_integrator_whfast512_field_descriptor_list,
};

const struct reb_binarydata_field_descriptor reb_integrator_whfast512_field_descriptor_list[] = {
    { "If this flag is set to 1 (default is 0) then general relativistic corrections are included. "
        "The corrections are in the form of an additional potential term and reproduce the correct precession rate. "
        "The constants are hard coded for this effect and assume that the simulation is in units of G=1 and one length unit corresponds to one astronomical unit. ",
        REB_UINT,        "gr_potential",    offsetof(struct reb_integrator_whfast512_state, gr_potential), 0, 0, 0},
    { "If this flag is set to 17 (default is 0), then symplectic correctors are used.", 
        REB_UINT,        "corrector",       offsetof(struct reb_integrator_whfast512_state, corrector), 0, 0, 0},
    { "If this is set to a number other than 1 then timesteps are combined. "
        "By doing multiple timesteps in a row, WHFast512 can keep all simulation data in registers which significantly speeds up the calculation. "
        "This number should be as large as the output cadence allows. "
        "The default is 1e6. ",
        REB_UINT64,      "concatenate_steps", offsetof(struct reb_integrator_whfast512_state, concatenate_steps), 0, 0, 0},
    { "By default this value is set to 1, implying all 8 particles in the simulation correspond to one system. "
        "By setting N_systems to either 2 or 4, one can integrate multiple planetary systems with 2, 3, or 4 particles at the same time. "
        "See the example problems on how to setup the particles for this case. ",
        REB_UINT,        "N_systems",       offsetof(struct reb_integrator_whfast512_state, N_systems), 0, 0, 0},
    { "", REB_POINTER_ALIGNED, "data",        offsetof(struct reb_integrator_whfast512_state, data), SIZE_MAX, sizeof(struct simd_data), 0},
    { "", REB_DOUBLE,      "last_synchronization", offsetof(struct reb_integrator_whfast512_state, last_synchronization), 0, 0, 0},
    { 0 }, // Null terminated list
};

void* reb_integrator_whfast512_create(){
    struct reb_integrator_whfast512_state* whfast512 = calloc(sizeof(struct reb_integrator_whfast512_state),1);
    whfast512->N_systems = 1;
    whfast512->gr_potential = 0;
    whfast512->concatenate_steps = 1e6;
    return whfast512;
}

void reb_integrator_whfast512_free(void* state){
    struct reb_integrator_whfast512_state* whfast512 = state;
    free(whfast512->data);
    free(whfast512);
}

#if (defined(__i386__) || defined(__x86_64__)) && !defined(_WIN32)
// Helper macro to print out offsets in structure for assembly code
#define SIMD_DATA_MEMBERS X(M) X(dt) X(gr_prefac) X(m) X(x) X(y) X(z) X(vx) X(vy) X(vz) \
X(mat8_inertial_to_jacobi) \
X(mat8_jacobi_to_heliocentric) \
X(M0) X(mask) \
X(exit_max_distance) X(exit_min_distance_r) \
X(mat8_jacobi_to_inertial)\
X(counter) 

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
    if (transformation != NULL){
        double tmp2[8] = {0}; 
        for (int i=0; i<8; i++) {
            for (int j=0; j<8; j++) {
                tmp2[i] += tmp[j]*transformation[8*j+i];
            }
        }
        return _mm512_loadu_pd(tmp2);
    }else{
        return _mm512_loadu_pd(tmp);
    }
}

// Hepler function to load particle data from avx512 registers
__attribute__((target("avx512f,avx512vl,avx512bw,avx512dq")))
static void load_from_m512d(struct reb_simulation* r, size_t offset, const double* transformation, int N_systems, __m512d vector){
    struct reb_particle* particles = r->particles;
    const unsigned int p_per_system = 8/N_systems;
    const unsigned int N_per_system = r->N/N_systems;
    double tmp[8];
    _mm512_storeu_pd(tmp, vector);
    double tmp2[8] = {0}; 
    for (int i=0; i<8; i++) {
        for (int j=0; j<8; j++) {
            tmp2[i] += tmp[j]*transformation[8*j+i];
        }
    }
    for (int s=0; s<N_systems; s++){
        for (unsigned int i=1; i<N_per_system; i++){
            *(double*)((char*)(&particles[s*N_per_system+i])+offset) = tmp2[s*p_per_system+i-1];
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

// External functions. Implemented in integrator_whfast512.s.
// _n2 = two systems of up to 4 planets, _n4 = four systems of 2 planets.

#define MACRO_4(gr, sys, mind, maxd) \
    extern enum REB_STATUS reb_whfast512_full_steps_##gr##_n##sys##_##encounter##_##escape(struct simd_data* data, uint64_t* N_steps, int skip_first_kepler_step, volatile sig_atomic_t* sigint);

#define MACRO_3(gr, sys, encounter) \
    MACRO_4(gr, sys, encounter, escape) \
    MACRO_4(gr, sys, encounter, noescape)

#define MACRO_2(gr, sys) \
    MACRO_3(gr, nsys, encounter) \
    MACRO_3(gr, nsys, noencounter) 

#define MACRO_1(gr) \
    MACRO_2(gr, 1) \
    MACRO_2(gr, 2) \
    MACRO_2(gr, 4) 

#define MACRO_0() \
    MACRO_1(0) \
    MACRO_1(1)

MACRO_0()


extern enum REB_STATUS reb_whfast512_kepler_step(struct simd_data* data);
extern void reb_whfast512_corrector_step_gr_n1(struct simd_data* data, double inv);
extern void reb_whfast512_corrector_step_nogr_n1(struct simd_data* data, double inv);
extern void reb_whfast512_corrector_step_gr_n2(struct simd_data* data, double inv);
extern void reb_whfast512_corrector_step_nogr_n2(struct simd_data* data, double inv);
extern void reb_whfast512_corrector_step_gr_n4(struct simd_data* data, double inv);
extern void reb_whfast512_corrector_step_nogr_n4(struct simd_data* data, double inv);

static enum REB_STATUS whfast512_full_steps(struct reb_integrator_whfast512_state* whfast512, uint64_t* N_steps, int skip){
    struct simd_data* data = whfast512->data;
    if (whfast512->gr_potential){
        switch (whfast512->N_systems){
            case 2: return reb_whfast512_full_steps_gr_n2(data, N_steps, skip, &reb_sigint); break;
            case 4: return reb_whfast512_full_steps_gr_n4(data, N_steps, skip, &reb_sigint); break;
            default: return reb_whfast512_full_steps_gr_n1(data, N_steps, skip, &reb_sigint); break;
        }
    }else{
        switch (whfast512->N_systems){
            case 2: return reb_whfast512_full_steps_nogr_n2(data, N_steps, skip, &reb_sigint); break;
            case 4: return reb_whfast512_full_steps_nogr_n4(data, N_steps, skip, &reb_sigint); break;
            default: return reb_whfast512_full_steps_nogr_n1(data, N_steps, skip, &reb_sigint); break;
        }
    }
}

static void whfast512_corrector_step(struct reb_integrator_whfast512_state* whfast512, double inv){
    struct simd_data* data = whfast512->data;
    if (whfast512->gr_potential){
        switch (whfast512->N_systems){
            case 2: reb_whfast512_corrector_step_gr_n2(data, inv); break;
            case 4: reb_whfast512_corrector_step_gr_n4(data, inv); break;
            default: reb_whfast512_corrector_step_gr_n1(data, inv); break;
        }
    }else{
        switch (whfast512->N_systems){
            case 2: reb_whfast512_corrector_step_nogr_n2(data, inv); break;
            case 4: reb_whfast512_corrector_step_nogr_n4(data, inv); break;
            default: reb_whfast512_corrector_step_nogr_n1(data, inv); break;
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
    data->exit_max_distance = _mm512_set1_pd(r->exit_max_distance);
    data->exit_min_distance_r = _mm512_set1_pd(1.0/r->exit_min_distance);
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
    r->gravity = REB_GRAVITY_NONE; // WHFast512 uses its own gravity routine.
    return 0; // success
}

// Optimized main loops allowing for concatenate_steps
__attribute__((target("avx512f,avx512vl,avx512bw,avx512dq")))
void reb_integrator_whfast512_step(struct reb_simulation* const r, void* state){
    struct reb_integrator_whfast512_state* whfast512 = state;
    const double dt = r->dt;
    uint64_t N_steps = whfast512->concatenate_steps;

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

    r->status = whfast512_full_steps(whfast512, &N_steps, skip_first_kepler_step);

    r->is_synchronized = 0;
    r->t += dt*N_steps;     // Note: N_steps might have been changed by whfast512_full_steps.
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

#else // (defined(__i386__) || defined(__x86_64__)) && !defined(_WIN32)
void reb_integrator_whfast512_step(struct reb_simulation* r, void* state){
    (void)state;
    reb_simulation_error(r, "AVX512 is not supported on your platform");
}
void reb_integrator_whfast512_synchronize(struct reb_simulation* r, void* state){
    (void)state;
    reb_simulation_error(r, "AVX512 is not supported on your platform");
}
#endif // (defined(__i386__) || defined(__x86_64__)) && !defined(_WIN32)
