#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "rebound.h"
#include "integrator_asm512.h"
#include "integrator_asm512_opt.h"

#ifdef AVX512
#include <immintrin.h>

// Mirrors struct simd_data in integrator_asm512.c. Must stay in sync with
// the layout that integrator_asm512.s / integrator_asm512_opt.s expect.
struct simd_data {
    __m512d M  __attribute__ ((aligned (64)));
    __m512d dt __attribute__ ((aligned (64)));
    __m512d gr_prefac  __attribute__ ((aligned (64)));
    __m512d gr_prefac2 __attribute__ ((aligned (64)));
    __m512d m  __attribute__ ((aligned (64)));
    __m512d x  __attribute__ ((aligned (64)));
    __m512d y  __attribute__ ((aligned (64)));
    __m512d z  __attribute__ ((aligned (64)));
    __m512d vx __attribute__ ((aligned (64)));
    __m512d vy __attribute__ ((aligned (64)));
    __m512d vz __attribute__ ((aligned (64)));
    double  mat8_inertial_to_jacobi[64]      __attribute__ ((aligned (64)));
    double  mat8_jacobi_to_heliocentric[64]  __attribute__ ((aligned (64)));
    __m512d M0  __attribute__ ((aligned (64)));
    __mmask8 mask __attribute__ ((aligned (64)));
    double  mat8_jacobi_to_inertial[64] __attribute__ ((aligned (64)));
    __m512i counter __attribute__ ((aligned (64)));
};

static inline __m512d mat8_mul_avx512(const double* matrix, const __m512d vector) {
    __m512d v_i = _mm512_set1_pd(vector[0]);
    __m512d col_i = _mm512_load_pd(matrix);
    __m512d res = _mm512_mul_pd(v_i, col_i);
    for (int i = 1; i < 8; i++) {
        __m512d v_i_in = _mm512_set1_pd(vector[i]);
        __m512d col_i_in = _mm512_load_pd(&matrix[i * 8]);
        res = _mm512_fmadd_pd(v_i_in, col_i_in, res);
    }
    return res;
}

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

static void load_from_m512d(struct reb_simulation* r, size_t offset, const double* transformation, int N_systems, __m512d vector){
    struct reb_particle* particles = r->particles;
    const unsigned int p_per_system = 8/N_systems;
    const unsigned int N_per_system = r->N/N_systems;
    double tmp[8] __attribute__((aligned(64)));
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

static void jacobi_to_inertial_posvel_and_com(struct reb_simulation* r, struct simd_data* data, double dt_com, unsigned int N_systems){
    const unsigned int N_per_system = r->N/N_systems;
    struct reb_particle com[4];
    for (unsigned s=0;s<N_systems;s++){
        com[s] = reb_simulation_com_range(r,s*N_per_system, (s+1)*N_per_system);
    }
    struct reb_particle* particles = r->particles;
    load_from_m512d(r, offsetof(struct reb_particle, x),  data->mat8_jacobi_to_inertial, N_systems, data->x);
    load_from_m512d(r, offsetof(struct reb_particle, y),  data->mat8_jacobi_to_inertial, N_systems, data->y);
    load_from_m512d(r, offsetof(struct reb_particle, z),  data->mat8_jacobi_to_inertial, N_systems, data->z);
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

static void inertial_to_jacobi_posvel(struct reb_simulation* r, struct simd_data* data, unsigned int N_systems){
    const unsigned int N_per_system = r->N/N_systems;
    struct reb_particle com[4];
    struct reb_particle* p_tmp = malloc(sizeof(struct reb_particle)*r->N);
    memcpy(p_tmp, r->particles, sizeof(struct reb_particle)*r->N);
    for (unsigned int s=0;s<N_systems;s++){
        com[s] = reb_simulation_com_range(r,s*N_per_system, (s+1)*N_per_system);
    }
    for (unsigned int s=0;s<N_systems;s++){
        for (unsigned int i=0;i<N_per_system;i++){
            r->particles[s*N_per_system+i].x  -= com[s].x;
            r->particles[s*N_per_system+i].y  -= com[s].y;
            r->particles[s*N_per_system+i].z  -= com[s].z;
            r->particles[s*N_per_system+i].vx -= com[s].vx;
            r->particles[s*N_per_system+i].vy -= com[s].vy;
            r->particles[s*N_per_system+i].vz -= com[s].vz;
        }
    }
    reb_simulation_move_to_com(r);
    data->x  = load_into_m512d(r, offsetof(struct reb_particle,x), data->mat8_inertial_to_jacobi, N_systems);
    data->y  = load_into_m512d(r, offsetof(struct reb_particle,y), data->mat8_inertial_to_jacobi, N_systems);
    data->z  = load_into_m512d(r, offsetof(struct reb_particle,z), data->mat8_inertial_to_jacobi, N_systems);
    data->vx = load_into_m512d(r, offsetof(struct reb_particle,vx), data->mat8_inertial_to_jacobi, N_systems);
    data->vy = load_into_m512d(r, offsetof(struct reb_particle,vy), data->mat8_inertial_to_jacobi, N_systems);
    data->vz = load_into_m512d(r, offsetof(struct reb_particle,vz), data->mat8_inertial_to_jacobi, N_systems);
    data->m  = load_into_m512d(r, offsetof(struct reb_particle,m),  NULL, N_systems);
    memcpy(r->particles, p_tmp, sizeof(struct reb_particle)*r->N);
    free(p_tmp);
}

static void recalculate_constants(struct reb_simulation* r, struct simd_data* data, unsigned int N_systems){
    struct reb_particle* particles = r->particles;
    const unsigned int p_per_system = 8/N_systems;
    const unsigned int N_per_system = r->N/N_systems;
    double mat8_inertial_to_heliocentric[64];
    double M[8]  = {0.0};
    double M0[8] = {0.0};
    switch (N_systems){
        case 1: data->mask = (1 << (r->N -1)) - 1; break;
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
    for (unsigned int i=1; i<9; i++){
        for (unsigned int j=1; j<9; j++){
            data->mat8_inertial_to_jacobi[(i-1)+8*(j-1)] = 0.0;
            data->mat8_jacobi_to_inertial[(i-1)+8*(j-1)] = 0.0;
            data->mat8_jacobi_to_heliocentric[(i-1)+8*(j-1)] = 0.0;
            mat8_inertial_to_heliocentric[(i-1)+8*(j-1)] = 0.0;
        }
    }
    for (unsigned int s=0; s<N_systems; s++){
        for (unsigned int i=0;i<N_per_system-1;i++){
            for (unsigned int j=0;j<i+2;j++){
                M[s*p_per_system+i] += r->particles[s*N_per_system+j].m;
            }
        }
    }
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
    for (unsigned int i=1; i<9; i++){
        for (unsigned int j=1; j<9; j++){
            for (unsigned int k=1; k<9; k++){
                data->mat8_jacobi_to_heliocentric[(i-1)+8*(j-1)] += mat8_inertial_to_heliocentric[(i-1)+8*(k-1)] * data->mat8_jacobi_to_inertial[(k-1)+8*(j-1)];
            }
        }
    }
    data->M  = _mm512_loadu_pd(&M);
    data->M0 = _mm512_loadu_pd(&M0);

    double c = 10065.32;
    double _gr_prefac[8]  = {0};
    double _gr_prefac2[8] = {0};
    for (unsigned int s=0; s<N_systems; s++){
        double m0 = r->particles[s*N_per_system].m;
        for (unsigned int p=1; p<N_per_system; p++){
            _gr_prefac[s*p_per_system+(p-1)]  = -6.*m0*m0/(c*c);
            _gr_prefac2[s*p_per_system+(p-1)] = -r->particles[s*N_per_system+p].m / m0;
        }
    }
    data->gr_prefac  = _mm512_loadu_pd(&_gr_prefac);
    data->gr_prefac2 = _mm512_loadu_pd(&_gr_prefac2);
    data->dt         = _mm512_set1_pd(r->dt);
}

extern void reb_asm512_opt_full_steps_gr(struct simd_data* data, long N_steps, int skip_first_kepler_step);
extern void reb_asm512_opt_full_steps_nogr(struct simd_data* data, long N_steps, int skip_first_kepler_step);
extern void reb_asm512_mom_full_steps_gr(struct simd_data* data, long N_steps, int skip_first_kepler_step);
extern void reb_asm512_mom_full_steps_nogr(struct simd_data* data, long N_steps, int skip_first_kepler_step);
extern void reb_asm512_fused_full_steps_gr(struct simd_data* data, long N_steps, int skip_first_kepler_step);
extern void reb_asm512_fused_full_steps_nogr(struct simd_data* data, long N_steps, int skip_first_kepler_step);
extern void reb_asm512_opt_corrector_step_gr(struct simd_data* data, double inv);
extern void reb_asm512_opt_corrector_step_nogr(struct simd_data* data, double inv);
extern void reb_asm512_opt_kepler_step(struct simd_data* data);

typedef void (*asm512_full_steps_fn)(struct simd_data*, long, int);

static int verify_setup_opt(struct reb_simulation* const r){
    struct reb_integrator_asm512_state* asm512 = r->integrator.state;
    if (r->N_var != 0
        || r->exact_finish_time != 0
        || (r->N > 9 && asm512->N_systems == 1)
        || asm512->N_systems != 1
        || r->G != 1.0
        || (r->N_active != SIZE_MAX && r->N_active != r->N)){
        reb_simulation_error(r, "asm512_opt setup invalid (mirror asm512 constraints).");
        r->status = REB_STATUS_GENERIC_ERROR;
        return 1;
    }
    asm512->N_allocated = 1;
    r->gravity = REB_GRAVITY_NONE;
    return 0;
}

static void asm512_variant_step(struct reb_simulation* const r, void* state,
                                asm512_full_steps_fn fn_gr,
                                asm512_full_steps_fn fn_nogr){
    struct reb_integrator_asm512_state* asm512 = state;
    const double dt = r->dt;
    const unsigned int N_steps = asm512->concatenate_steps;
    if (verify_setup_opt(r)){
        return;
    }
    if (asm512->recalculate_constants){
        recalculate_constants(r, asm512->data, asm512->N_systems);
        asm512->recalculate_constants = 0;
    }
    struct simd_data* data = asm512->data;
    int skip_first_kepler_step = 0;
    if (r->is_synchronized){
        inertial_to_jacobi_posvel(r, data, asm512->N_systems);
        if (asm512->corrector){
            if (asm512->gr_potential){
                reb_asm512_opt_corrector_step_gr(data, 1.0);
            }else{
                reb_asm512_opt_corrector_step_nogr(data, 1.0);
            }
        }
        skip_first_kepler_step = 1;
        data->dt = _mm512_set1_pd(dt/2.0);
        reb_asm512_opt_kepler_step(data);
        data->dt = _mm512_set1_pd(dt);
    }
    if (asm512->gr_potential){
        fn_gr(asm512->data, N_steps, skip_first_kepler_step);
    }else{
        fn_nogr(asm512->data, N_steps, skip_first_kepler_step);
    }
    r->is_synchronized = 0;
    r->t += dt*N_steps;
    r->dt_last_done = dt;
}

static void integrator_asm512_opt_step(struct reb_simulation* const r, void* state){
    asm512_variant_step(r, state, reb_asm512_opt_full_steps_gr, reb_asm512_opt_full_steps_nogr);
}
static void integrator_asm512_mom_step(struct reb_simulation* const r, void* state){
    asm512_variant_step(r, state, reb_asm512_mom_full_steps_gr, reb_asm512_mom_full_steps_nogr);
}
static void integrator_asm512_fused_step(struct reb_simulation* const r, void* state){
    asm512_variant_step(r, state, reb_asm512_fused_full_steps_gr, reb_asm512_fused_full_steps_nogr);
}

static void integrator_asm512_opt_synchronize(struct reb_simulation* const r, void* state){
    struct reb_integrator_asm512_state* const asm512 = state;
    if (r->is_synchronized){
        return;
    }
    struct simd_data* data = asm512->data;
    data->dt = _mm512_set1_pd(r->dt/2.0);
    reb_asm512_opt_kepler_step(data);
    data->dt = _mm512_set1_pd(r->dt);
    if (asm512->gr_potential){
        reb_asm512_opt_corrector_step_gr(data, -1.0);
    }else{
        reb_asm512_opt_corrector_step_nogr(data, -1.0);
    }
    jacobi_to_inertial_posvel_and_com(r, data, 0.0, asm512->N_systems);
    r->is_synchronized = 1;
}

#else

static void integrator_asm512_opt_step(struct reb_simulation* const r, void* state){
    (void)state;
    reb_simulation_error(r, "asm512_opt requires AVX512.");
    r->status = REB_STATUS_GENERIC_ERROR;
}
static void integrator_asm512_mom_step(struct reb_simulation* const r, void* state){
    (void)state;
    reb_simulation_error(r, "asm512_mom requires AVX512.");
    r->status = REB_STATUS_GENERIC_ERROR;
}
static void integrator_asm512_fused_step(struct reb_simulation* const r, void* state){
    (void)state;
    reb_simulation_error(r, "asm512_fused requires AVX512.");
    r->status = REB_STATUS_GENERIC_ERROR;
}

static void integrator_asm512_opt_synchronize(struct reb_simulation* const r, void* state){
    (void)r; (void)state;
}

#endif

extern void* reb_integrator_asm512_create(void);
extern void  reb_integrator_asm512_free(void* state);

static const struct reb_binarydata_field_descriptor asm512_opt_fields[] = {
    { 0 },
};

const struct reb_integrator reb_integrator_asm512_opt = {
    .step                  = integrator_asm512_opt_step,
    .create                = reb_integrator_asm512_create,
    .free                  = reb_integrator_asm512_free,
    .synchronize           = integrator_asm512_opt_synchronize,
    .field_descriptor_list = asm512_opt_fields,
};

const struct reb_integrator reb_integrator_asm512_mom = {
    .step                  = integrator_asm512_mom_step,
    .create                = reb_integrator_asm512_create,
    .free                  = reb_integrator_asm512_free,
    .synchronize           = integrator_asm512_opt_synchronize,
    .field_descriptor_list = asm512_opt_fields,
};

const struct reb_integrator reb_integrator_asm512_fused = {
    .step                  = integrator_asm512_fused_step,
    .create                = reb_integrator_asm512_create,
    .free                  = reb_integrator_asm512_free,
    .synchronize           = integrator_asm512_opt_synchronize,
    .field_descriptor_list = asm512_opt_fields,
};
