/**
 * integrator_asm512.c: ASM version of WHFast512
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
#include <linux/perf_event.h>
#include <asm/unistd.h>
#include "rebound.h"
#include "particle.h"
#include "tools.h"
#include "gravity.h"
#include "boundary.h"
#include "integrator_whfast.h"
#include "integrator_asm512.h"

#ifdef AVX512
// Debug function to print vectors
static inline void printavx512(__m512d a) {
    double _nax[8];
    _mm512_storeu_pd(&_nax[0], a);
    printf("avx = {%.17g, %.17g, %.17g, %.17g, %.17g, %.17g, %.17g, %.17g\n", _nax[0], _nax[1], _nax[2], _nax[3], _nax[4], _nax[5], _nax[6], _nax[7]);
}

// Print mask in binary format for debuggin
static inline void printmask8(__mmask8 mask) {
    for (int i = 7; i >= 0; i--) {
        printf("%d", (mask >> i) & 1);
        if (i == 4) printf(" "); 
    }
    printf("\n");
}
const int reb_integrator_asm512_available = 1;   // Let python check if REBOUND was compiled with AVX512 enabled
#else //AVX512
const int reb_integrator_asm512_available = 0;   // Let python check if REBOUND was compiled with AVX512 enabled
#endif //AVX512

// Debug function to print 8x8 matrix
static inline void printmat8(double* a) {
    for (int i=0; i<8;i++){
        for (int j=0; j<8;j++){
            printf("%.16f ", a[i*8+j]);
        }
        printf("\n");
    }
}

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
#endif

                
// Convert jacobi coordinates to inertial coordinates
// Also performs com step (assume original particles are unmodified)
// Note: Speed is not a concern here 
static void jacobi_to_inertial_posvel_and_com(struct reb_simulation* r, double dt_com){
    struct reb_integrator_asm512* const ri_whfast512 = &(r->ri_whfast512);
    const unsigned int N_systems = ri_whfast512->N_systems;
    const unsigned int N_per_system = r->N/N_systems;
    struct reb_particle com[4];
    for (unsigned s=0;s<N_systems;s++){
        com[s] = reb_simulation_com_range(r,s*N_per_system, (s+1)*N_per_system); // original com
    }
#ifdef AVX512
    struct reb_particle* particles = r->particles;
    load_from_m512d(r, offsetof(struct reb_particle, x), ri_whfast512->p512->mat8_jacobi_to_inertial, N_systems, ri_whfast512->p512->x);
    load_from_m512d(r, offsetof(struct reb_particle, y), ri_whfast512->p512->mat8_jacobi_to_inertial, N_systems, ri_whfast512->p512->y);
    load_from_m512d(r, offsetof(struct reb_particle, z), ri_whfast512->p512->mat8_jacobi_to_inertial, N_systems, ri_whfast512->p512->z);
    load_from_m512d(r, offsetof(struct reb_particle, vx), ri_whfast512->p512->mat8_jacobi_to_inertial, N_systems, ri_whfast512->p512->vx);
    load_from_m512d(r, offsetof(struct reb_particle, vy), ri_whfast512->p512->mat8_jacobi_to_inertial, N_systems, ri_whfast512->p512->vy);
    load_from_m512d(r, offsetof(struct reb_particle, vz), ri_whfast512->p512->mat8_jacobi_to_inertial, N_systems, ri_whfast512->p512->vz);
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
            particles[s*N_per_system+i].x  += com[s].x + com[s].vx*dt_com;
            particles[s*N_per_system+i].y  += com[s].y + com[s].vy*dt_com;
            particles[s*N_per_system+i].z  += com[s].z + com[s].vz*dt_com;
            particles[s*N_per_system+i].vx += com[s].vx;
            particles[s*N_per_system+i].vy += com[s].vy;
            particles[s*N_per_system+i].vz += com[s].vz;
        }
    }
#else  // AVX512
    reb_simulation_error(r, "Fallback for Jacobi transformations in whfast512 not yet implemented.");
#endif // AVX512
}

#ifdef AVX512

extern void block1_gr(struct reb_particle_avx512* p512, long N_steps);
extern void block1_nogr(struct reb_particle_avx512* p512, long N_steps);
extern void reb_whfast512_kepler_step(struct reb_particle_avx512* p512);

static void inertial_to_jacobi_posvel(struct reb_simulation* r){
    struct reb_integrator_asm512* const ri_whfast512 = &(r->ri_whfast512);
    const unsigned int N_systems = ri_whfast512->N_systems;
    const unsigned int N_per_system = r->N/N_systems;
    // Transformations assume system is in COM frame.
    struct reb_particle com[4];
    struct reb_particle* p_tmp = malloc(sizeof(struct reb_particle)*r->N);
    memcpy(p_tmp, r->particles, sizeof(struct reb_particle)*r->N);
    for (unsigned s=0;s<N_systems;s++){
        com[s] = reb_simulation_com_range(r,s*N_per_system, (s+1)*N_per_system); // original com
    }
    for (unsigned s=0;s<N_systems;s++){
        for (int i=0;i<N_per_system;i++){
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
    ri_whfast512->p512->x = load_into_m512d(r, offsetof(struct reb_particle,x),ri_whfast512->p512->mat8_inertial_to_jacobi, N_systems);
    ri_whfast512->p512->y = load_into_m512d(r, offsetof(struct reb_particle,y),ri_whfast512->p512->mat8_inertial_to_jacobi, N_systems);
    ri_whfast512->p512->z = load_into_m512d(r, offsetof(struct reb_particle,z),ri_whfast512->p512->mat8_inertial_to_jacobi, N_systems);
    ri_whfast512->p512->vx = load_into_m512d(r, offsetof(struct reb_particle,vx),ri_whfast512->p512->mat8_inertial_to_jacobi, N_systems);
    ri_whfast512->p512->vy = load_into_m512d(r, offsetof(struct reb_particle,vy),ri_whfast512->p512->mat8_inertial_to_jacobi, N_systems);
    ri_whfast512->p512->vz = load_into_m512d(r, offsetof(struct reb_particle,vz),ri_whfast512->p512->mat8_inertial_to_jacobi, N_systems);
    ri_whfast512->p512->m = load_into_m512d(r, offsetof(struct reb_particle,m),NULL, N_systems);
    // Undo COM transformation. COM will be applied in jacobi_to_inertial_posvel_and_com().
    memcpy(r->particles, p_tmp, sizeof(struct reb_particle)*r->N);
    free(p_tmp);
}


// Precalculate various constants and put them in 512 bit vectors.
void static recalculate_constants(struct reb_simulation* r){
    struct reb_integrator_asm512* const ri_whfast512 = &(r->ri_whfast512);
    struct reb_particle* particles = r->particles;
    const unsigned int N_systems = ri_whfast512->N_systems;
    const unsigned int p_per_system = 8/N_systems;
    const unsigned int N_per_system = r->N/N_systems;
    double mat8_inertial_to_heliocentric[64];
    double M[8] = {0.0};
    double M0[8] = {0.0};
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
            // Zeroing.
            for (unsigned int i=1; i<9; i++){
                for (unsigned int j=1; j<9; j++){
                    ri_whfast512->p512->mat8_inertial_to_jacobi[(i-1)+8*(j-1)] = 0.0;
                    ri_whfast512->p512->mat8_jacobi_to_inertial[(i-1)+8*(j-1)] = 0.0;
                    ri_whfast512->p512->mat8_jacobi_to_heliocentric[(i-1)+8*(j-1)] = 0.0;
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
                        ri_whfast512->p512->mat8_inertial_to_jacobi[(s*p_per_system+i-1)+8*(s*p_per_system+j-1)] += particles[s*N_per_system+j].m/ms;
                    }
                    for (unsigned int j=1; j<N_per_system; j++){
                        mat8_inertial_to_heliocentric[(s*p_per_system+i-1)+8*(s*p_per_system+j-1)] += particles[s*N_per_system+j].m/particles[s*N_per_system+0].m;
                    }
                    mat8_inertial_to_heliocentric[(s*p_per_system+i-1)+8*(s*p_per_system+i-1)] += 1.0;
                    ri_whfast512->p512->mat8_inertial_to_jacobi[(s*p_per_system+i-1)+8*(s*p_per_system+i-1)] += 1.0;
                    ri_whfast512->p512->mat8_jacobi_to_inertial[(s*p_per_system+i-1)+8*(s*p_per_system+i-1)] += ms/(ms + particles[s*N_per_system+i].m);
                    ms += particles[s*N_per_system+i].m;
                    if (i<N_per_system-1){ 
                        for (unsigned int ii=i; ii>0; ii--){
                            int jj = i+1;
                            ri_whfast512->p512->mat8_jacobi_to_inertial[(s*p_per_system+ii-1)+8*(s*p_per_system+jj-1)] -= particles[s*N_per_system+jj].m/(ms+particles[s*N_per_system+jj].m);
                        }
                    }
                    M0[(s*p_per_system+i-1)] = -particles[s*N_per_system+0].m*r->dt;
                }
            }

            // Might be numerically more stable to calculate this manually rather than do a matrix multiplication.
            for (unsigned int i=1; i<9; i++){
                for (unsigned int j=1; j<9; j++){
                    for (unsigned int k=1; k<9; k++){
                        ri_whfast512->p512->mat8_jacobi_to_heliocentric[(i-1)+8*(j-1)] += mat8_inertial_to_heliocentric[(i-1)+8*(k-1)] * ri_whfast512->p512->mat8_jacobi_to_inertial[(k-1)+8*(j-1)];
                    }
                }
            }
            //// Quick transpose test
            //// TODO Cleanup.
            //double tmp[64];
            //for (unsigned int i=0; i<8; i++){
            //    for (unsigned int j=0; j<8; j++){
            //        tmp[i+8*j] = ri_whfast512->p512->mat8_jacobi_to_heliocentric[i+8*j];
            //    }
            //}
            //for (unsigned int i=0; i<8; i++){
            //    for (unsigned int j=0; j<8; j++){
            //        ri_whfast512->p512->mat8_jacobi_to_heliocentric[i+8*j] = tmp[j+8*i];
            //        ri_whfast512->p512->mat8_inertial_to_jacobi_T[i+8*j] = ri_whfast512->p512->mat8_inertial_to_jacobi[j+8*i];
            //    }
            //}

    ri_whfast512->p512->M = _mm512_loadu_pd(&M);
    ri_whfast512->p512->M0 = _mm512_loadu_pd(&M0); //  = -particles[0].m * dt

    // GR and jump prefactors. Note: assumes units of AU, year/2pi.
    double c = 10065.32;
    double _gr_prefac[8];
    double _gr_prefac2[8];
    double _jump_prefac[8];
    for(unsigned int i=0;i<8;i++){
        _gr_prefac[i] = 0; // for when N<8
        _gr_prefac2[i] = 0;
        _jump_prefac[i] = 0;
    }
    for (int s=0; s<N_systems; s++){
        double m0 = r->particles[s*N_per_system].m;
        for (int p=1; p<N_per_system; p++){
                    _gr_prefac[s*p_per_system+(p-1)] = -r->dt*6.*m0*m0/(c*c);
                    _gr_prefac2[s*p_per_system+(p-1)] = -r->particles[s*N_per_system+p].m / m0;
        }
    }
    ri_whfast512->p512->gr_prefac = _mm512_loadu_pd(&_gr_prefac);
    ri_whfast512->p512->gr_prefac2 = _mm512_loadu_pd(&_gr_prefac2);
    ri_whfast512->p512->dt = _mm512_set1_pd(r->dt); 

    ri_whfast512->recalculate_constants = 0;
}

static int reb_integrator_whfast512_allocate(struct reb_simulation* const r){
    struct reb_integrator_asm512* const ri_whfast512 = &(r->ri_whfast512);
    // Check if all assumptions are satisfied.
    // Note: These are not checked every timestep. 
    // So it is possible for the user to screw things up.
    if (r->dt<0.0){
        reb_simulation_error(r, "WHFast512 does not support negative timesteps. To integrate backwards, flip the sign of the velocities.");
        r->status = REB_STATUS_GENERIC_ERROR;
        return 1;
    }
    if (r->N_var!=0){
        reb_simulation_error(r, "WHFast512 does not support variational particles.");
        r->status = REB_STATUS_GENERIC_ERROR;
        return 1;
    }
    if (ri_whfast512->corrector && ri_whfast512->coordinates != REB_WHFAST512_COORDINATES_JACOBI){
        reb_simulation_error(r, "Symplectic correctors require Jacobi coordinates.");
        r->status = REB_STATUS_GENERIC_ERROR;
        return 1;
    }
    if (r->exact_finish_time!=0){
        reb_simulation_error(r, "WHFast512 requires exact_finish_time=0.");
        r->status = REB_STATUS_GENERIC_ERROR;
        return 1;
    }
    if (r->N>9 && ri_whfast512->N_systems == 1) {
        reb_simulation_error(r, "WHFast512 supports a maximum of 9 particles when N_systems is set to 1.");
        r->status = REB_STATUS_GENERIC_ERROR;
        return 1;
    }
    if (r->N>10 && ri_whfast512->N_systems == 2) {
        reb_simulation_error(r, "WHFast512 supports a maximum of 10 particles when N_systems is set to 2.");
        r->status = REB_STATUS_GENERIC_ERROR;
        return 1;
    }
    if (r->N>12 && ri_whfast512->N_systems == 4) {
        reb_simulation_error(r, "WHFast512 supports a maximum of 12 particles when N_systems is set to 4.");
        r->status = REB_STATUS_GENERIC_ERROR;
        return 1;
    }
    if (ri_whfast512->N_systems != 1 && ri_whfast512->N_systems !=2 && ri_whfast512->N_systems != 4){
        reb_simulation_error(r, "WHFast512 supports 1, 2, or 4 systems only.");
        r->status = REB_STATUS_GENERIC_ERROR;
        return 1;
    }
    if (r->N % ri_whfast512->N_systems != 0){
        reb_simulation_error(r, "Number of particles must be a multiple of ri_whfast512.N_systems.");
        r->status = REB_STATUS_GENERIC_ERROR;
        return 1;
    }
    if (r->G!=1.0){
        reb_simulation_error(r, "WHFast512 requires units in which G=1. Please rescale your system.");
        r->status = REB_STATUS_GENERIC_ERROR;
        return 1;
    }
    if (r->N_active!=-1 && r->N_active!=r->N){
        reb_simulation_error(r, "WHFast512 does not support test particles.");
        r->status = REB_STATUS_GENERIC_ERROR;
        return 1;
    }
    if (ri_whfast512->p512==NULL){
        ri_whfast512->p512 = aligned_alloc(64,sizeof(struct reb_particle_avx512));
        if (!ri_whfast512->p512){
            reb_simulation_error(r, "WHFast512 was not able to allocate memory.");
            r->status = REB_STATUS_GENERIC_ERROR;
            return 1;
        }
    }
    if (r->exit_min_distance || r->exit_max_distance){
        reb_simulation_warning(r, "You are using WHFast512 together with the flags exit_min_distance and/or exit_max_distance. With the current implementation, these flags will only check the last synchronized positions. In addition they might slow down WHFast512 significantly. If you need to use these flags, please open an issue on GitHub for further advice.");
    }
    ri_whfast512->N_allocated=1;
    r->gravity = REB_GRAVITY_NONE; // WHFast512 uses its own gravity routine.
    return 0; // success
}

// Implementation of the GR force for WHFast correctors.
// (WHFast512 comes with built-in support) 
static void gr_force(struct reb_simulation* r){
    double C2 = 10065.32 * 10065.32;
    struct reb_particle* particles = r->particles;
    const struct reb_particle source = particles[0];
    const double prefac1 = 6.*(r->G*source.m)*(r->G*source.m)/C2;
    for (int i=1; i<r->N; i++){
        const struct reb_particle p = particles[i];
        const double dx = p.x - source.x;
        const double dy = p.y - source.y;
        const double dz = p.z - source.z;
        const double r2 = dx*dx + dy*dy + dz*dz;
        const double prefac = prefac1/(r2*r2);

        particles[i].ax -= prefac*dx;
        particles[i].ay -= prefac*dy;
        particles[i].az -= prefac*dz;
        particles[0].ax += p.m/source.m*prefac*dx;
        particles[0].ay += p.m/source.m*prefac*dy;
        particles[0].az += p.m/source.m*prefac*dz;
    }
}

// Creates a new simulation just for doing the correctors.
// Slow, but convenient. Only called when a simulation is synchronized.
static void apply_corrector(struct reb_simulation* r, double direction){
    struct reb_integrator_asm512* restrict const ri_whfast512 = &(r->ri_whfast512);
    if (ri_whfast512->corrector){
        const unsigned int N_systems = ri_whfast512->N_systems;
        const unsigned int N_per_system = r->N/N_systems;
        for (int s=0; s<N_systems; s++){
            struct reb_simulation* rt = reb_simulation_create();
            rt->dt = r->dt;
            rt->G = r->G;
            if (ri_whfast512->gr_potential){
                rt->additional_forces = gr_force;
            }
            for (int i=0;i<N_per_system;i++){
                reb_simulation_add(rt, r->particles[s*N_per_system+i]);
            }
            reb_integrator_whfast_init(rt);
            reb_integrator_whfast_from_inertial(rt);
            reb_whfast_apply_corrector(rt, direction, ri_whfast512->corrector);
            reb_integrator_whfast_to_inertial(rt);
            for (int i=0;i<N_per_system;i++){
                r->particles[s*N_per_system+i] = rt->particles[i];
                r->particles[s*N_per_system+i].sim = r;
            }
            reb_simulation_free(rt);
        }
    }
}

// Optimized main loops allowing for concatenate_steps
void reb_integrator_whfast512_part1(struct reb_simulation* const r){
    struct reb_integrator_asm512* restrict const ri_whfast512 = &(r->ri_whfast512);
    const double dt = r->dt;
    const unsigned int N_steps = ri_whfast512->concatenate_steps;

    if (ri_whfast512->N_allocated==0){
        if (reb_integrator_whfast512_allocate(r)){
            return; // Error occured
        }
        ri_whfast512->recalculate_constants = 1;
    }

    if (ri_whfast512->recalculate_constants){
        recalculate_constants(r);
    }

    // Reset previously synchronized particle data when keep_unsynchronized is turned on.
    if (ri_whfast512->keep_unsynchronized && ri_whfast512->particles_keep_unsynchronized){
        memcpy(r->particles, ri_whfast512->particles_keep_unsynchronized, sizeof(struct reb_particle)*r->N);
        free(ri_whfast512->particles_keep_unsynchronized);
        ri_whfast512->particles_keep_unsynchronized = NULL;
        ri_whfast512->N_allocated_particles_keep_unsynchronized = 0;
    }
    // Normal initial synchronize step (not trigger when synchronize was kalled with keep_unsynchronized=1)
    if (ri_whfast512->is_synchronized){
        ri_whfast512->time_of_last_synchronize = r->t;
        // Use WHFast to apply the correctors.
        apply_corrector(r, 1.0);
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
        // First half DRIFT step. Note negative sign. We will do a full step below.
        ri_whfast512->p512->dt = _mm512_set1_pd(-dt/2.0); 
        local_reb_whfast512_kepler_step(r->ri_whfast512.p512);    
        ri_whfast512->p512->dt = _mm512_set1_pd(dt); // Reset
    }
        
#ifdef PERF
    int fd_insn   = setup_counter(PERF_TYPE_HARDWARE, PERF_COUNT_HW_INSTRUCTIONS);
    int fd_cycles = setup_counter(PERF_TYPE_HARDWARE, PERF_COUNT_HW_CPU_CYCLES);
    int fd_l1d    = setup_counter(PERF_TYPE_HW_CACHE,
                        PERF_COUNT_HW_CACHE_L1D |
                        (PERF_COUNT_HW_CACHE_OP_READ << 8) |
                        (PERF_COUNT_HW_CACHE_RESULT_ACCESS << 16));

    // start all
    ioctl(fd_insn,   PERF_EVENT_IOC_RESET,  0);
    ioctl(fd_cycles, PERF_EVENT_IOC_RESET,  0);
    ioctl(fd_l1d,    PERF_EVENT_IOC_RESET,  0);
    ioctl(fd_insn,   PERF_EVENT_IOC_ENABLE, 0);
    ioctl(fd_cycles, PERF_EVENT_IOC_ENABLE, 0);
    ioctl(fd_l1d,    PERF_EVENT_IOC_ENABLE, 0);

#endif // PERF

    // Tight inner loops for speed
    if (ri_whfast512->coordinates==REB_WHFAST512_COORDINATES_DEMOCRATICHELIOCENTRIC){
        if (ri_whfast512->N_systems==1){
            if (ri_whfast512->gr_potential){
                for (unsigned int i=0;i<N_steps;i++){
                    local_reb_whfast512_kepler_step(r->ri_whfast512.p512);       // full timestep
                    reb_whfast512_jump_step(r);         // half timestep
                    reb_whfast512_interaction_step_8planets_democraticheliocentric(r);
                    reb_whfast512_jump_step(r);         // half timestep
                }
            }else{
                for (unsigned int i=0;i<N_steps;i++){
                    local_reb_whfast512_kepler_step(r->ri_whfast512.p512);       // full timestep
                    reb_whfast512_jump_step(r);         // full timestep
                    reb_whfast512_interaction_step_8planets_democraticheliocentric(r);
                }
            }
        }else if (ri_whfast512->N_systems==2){
            if (ri_whfast512->gr_potential){
                for (unsigned int i=0;i<N_steps;i++){
                    local_reb_whfast512_kepler_step(r->ri_whfast512.p512);       // full timestep
                    reb_whfast512_jump_step(r);         // half timstep
                    reb_whfast512_interaction_step_4planets_democraticheliocentric(r);
                    reb_whfast512_jump_step(r);         // half timestep
                }
            }else{
                for (unsigned int i=0;i<N_steps;i++){
                    local_reb_whfast512_kepler_step(r->ri_whfast512.p512);       // full timestep
                    reb_whfast512_jump_step(r);         // full timestep
                    reb_whfast512_interaction_step_4planets_democraticheliocentric(r);
                }
            }
        }else if (ri_whfast512->N_systems==4){
            if (ri_whfast512->gr_potential){
                for (unsigned int i=0;i<N_steps;i++){
                    local_reb_whfast512_kepler_step(r->ri_whfast512.p512);       // full timestep
                    reb_whfast512_jump_step(r);         // half timestep
                    reb_whfast512_interaction_step_2planets_democraticheliocentric(r);
                    reb_whfast512_jump_step(r);         // half timestep
                }
            }else{
                for (unsigned int i=0;i<N_steps;i++){
                    local_reb_whfast512_kepler_step(r->ri_whfast512.p512);       // full timestep
                    reb_whfast512_jump_step(r);         // full timestep
                    reb_whfast512_interaction_step_2planets_democraticheliocentric(r);
                }
            }
        }
    }else{ //JACOBI
        if (ri_whfast512->N_systems==1){
            if (ri_whfast512->gr_potential){
                //for (unsigned int i=0;i<N_steps;i++){
                    block1_gr(ri_whfast512->p512, N_steps);
                //}
            }else{
                //for (unsigned int i=0;i<N_steps;i++){
                    block1_nogr(ri_whfast512->p512, N_steps);
                //}
            }
        }else if (ri_whfast512->N_systems==2){
            for (unsigned int i=0;i<N_steps;i++){
                local_reb_whfast512_kepler_step(r->ri_whfast512.p512);    // full timestep
                reb_whfast512_interaction_step_4planets_jacobi(r);
            }
        }else if (ri_whfast512->N_systems==4){
            for (unsigned int i=0;i<N_steps;i++){
                local_reb_whfast512_kepler_step(r->ri_whfast512.p512);    // full timestep
                reb_whfast512_interaction_step_2planets_jacobi(r);
            }
        }
    }
#ifdef PERF
     // stop all
    ioctl(fd_insn,   PERF_EVENT_IOC_DISABLE, 0);
    ioctl(fd_cycles, PERF_EVENT_IOC_DISABLE, 0);
    ioctl(fd_l1d,    PERF_EVENT_IOC_DISABLE, 0);

    long long insn, cycles, l1d;
    read(fd_insn,   &insn,   sizeof(insn));
    read(fd_cycles, &cycles, sizeof(cycles));
    read(fd_l1d,    &l1d,    sizeof(l1d));

    printf("Instructions:    %lld\n", insn);
    printf("Cycles:          %lld\n", cycles);
    printf("IPC:             %.3f\n", (double)insn / cycles);
    printf("L1D accesses:    %lld\n", l1d);
    printf("Mem ops/insn:    %.3f\n", (double)l1d / insn);

    close(fd_insn);
    close(fd_cycles);
    close(fd_l1d);

#endif // PERF

    ri_whfast512->is_synchronized = 0;
    r->t += dt*N_steps;
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
    struct reb_integrator_asm512* const ri_whfast512 = &(r->ri_whfast512);
    if (ri_whfast512->is_synchronized == 0){
        struct reb_particle_avx512* sync_pj = NULL;
        // Needed if no step has ever been done before (like SA)
        if (ri_whfast512->N_allocated==0){
            if (reb_integrator_whfast512_allocate(r)){
                return; // error occured
            }
            ri_whfast512->recalculate_constants = 1;
        }
        if (ri_whfast512->recalculate_constants){
            recalculate_constants(r);
        }
        if (ri_whfast512->keep_unsynchronized){
            free(ri_whfast512->particles_keep_unsynchronized);
            ri_whfast512->particles_keep_unsynchronized = malloc(sizeof(struct reb_particle)*r->N);
            memcpy(ri_whfast512->particles_keep_unsynchronized, r->particles, sizeof(struct reb_particle)*r->N);
            ri_whfast512->N_allocated_particles_keep_unsynchronized = r->N;
            sync_pj = aligned_alloc(64,sizeof(struct reb_particle_avx512));
            memcpy(sync_pj,ri_whfast512->p512, sizeof(struct reb_particle_avx512));
        }
        ri_whfast512->p512->dt = _mm512_set1_pd(r->dt/2.0); 
        reb_whfast512_kepler_step(r->ri_whfast512.p512);    
        ri_whfast512->p512->dt = _mm512_set1_pd(r->dt); // Reset
        switch (ri_whfast512->coordinates){
            case REB_WHFAST512_COORDINATES_DEMOCRATICHELIOCENTRIC:
                democraticheliocentric_to_inertial_posvel_and_com(r, r->t-ri_whfast512->time_of_last_synchronize);
                break;
            case REB_WHFAST512_COORDINATES_JACOBI:
                jacobi_to_inertial_posvel_and_com(r, r->t-ri_whfast512->time_of_last_synchronize);
                break;
            default:
                reb_simulation_error(r,"Coordinate system not supported.");
        }
        // Use WHFast to applyt the correctors
        apply_corrector(r, -1.0);
        if (ri_whfast512->keep_unsynchronized){
            memcpy(ri_whfast512->p512, sync_pj, sizeof(struct reb_particle_avx512));
            free(sync_pj);
        }else{
            ri_whfast512->is_synchronized = 1;
            ri_whfast512->time_of_last_synchronize = r->t;
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
    struct reb_integrator_asm512* const ri_whfast512 = &(r->ri_whfast512);
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
        democraticheliocentric_to_inertial_posvel_and_com(r, r->t-ri_whfast512->time_of_last_synchronize);
        ri_whfast512->is_synchronized = 1;
    }
}

// Free memory and reset all constants.
// This needs to be called when the timestep, the number of particles, masses, etc are changed, 
void reb_integrator_whfast512_reset(struct reb_simulation* const r){
    struct reb_integrator_asm512* const ri_whfast512 = &(r->ri_whfast512);
    free(ri_whfast512->p512);
    free(ri_whfast512->particles_keep_unsynchronized);
    ri_whfast512->particles_keep_unsynchronized = NULL;
    ri_whfast512->p512 = NULL;
    ri_whfast512->N_allocated = 0;
    ri_whfast512->N_allocated_particles_keep_unsynchronized = 0;
    ri_whfast512->N_systems = 1;
    ri_whfast512->gr_potential = 0;
    ri_whfast512->corrector = 0;
    ri_whfast512->coordinates = REB_WHFAST512_COORDINATES_JACOBI;
    ri_whfast512->is_synchronized = 1;
    ri_whfast512->time_of_last_synchronize = 0.0;
    ri_whfast512->keep_unsynchronized = 0;
    ri_whfast512->concatenate_steps = 1;
}

// Everything is in part 1 for this integrator
void reb_integrator_whfast512_part2(struct reb_simulation* const r){
}
