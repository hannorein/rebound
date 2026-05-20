
#include "rebound.h"
#include "rebound_internal.h"
#include <string.h>
#include <math.h>
#include "transformations.h"
#include "particle.h"
#include "tools.h"
#include "gravity.h"
#include "boundary.h"
#include "integrator_whfast_HJ.h"

#define MAX(a, b) ((a) < (b) ? (b) : (a))   ///< Returns the maximum of a and b
#define MIN(a, b) ((a) > (b) ? (b) : (a))   ///< Returns the minimum of a and b

void* reb_integrator_whfast_hj_create();
void reb_integrator_whfast_hj_free(void* state);
void reb_integrator_whfast_hj_step(struct reb_simulation* const r, void* state);

const struct reb_integrator reb_integrator_whfast_hj = {
    .step = reb_integrator_whfast_hj_step,
    .create = reb_integrator_whfast_hj_create,
    .free = reb_integrator_whfast_hj_free,
};


void* reb_integrator_whfast_hj_create(){
    // Allocate memory and set default parameters.
    struct reb_integrator_whfast_hj_state* whfast = calloc(sizeof(struct reb_integrator_whfast_hj_state),1);
    return whfast;
}

void reb_integrator_whfast_hj_free(void* p){
    struct reb_integrator_whfast_hj_state* whfast = p;
    free(whfast->p_jh);
    free(whfast);
}


/***************************** 
 * Interaction Hamiltonian  */
void reb_integrator_whfast_hj_interaction_step(struct reb_simulation* const r, struct reb_particle* p_jh, const double _dt){
    const size_t N = r->N;
    const size_t N_active = (r->N_active==SIZE_MAX || r->testparticle_type ==1)?N:r->N_active;
    const double G = r->G;
    struct reb_particle* particles = r->particles;
    const double m0 = particles[0].m;
    reb_transformations_inertial_to_jacobi_acc(particles, p_jh, particles, N, N_active);
    double eta = m0;
    for (size_t i=1;i<N;i++){
        // Eq 132
        const struct reb_particle pji = p_jh[i];
        if (i<N_active){
            eta += pji.m;
        }
        // ax was calculate by update_acceleration O(N^2), last term from the right
        p_jh[i].vx += _dt * pji.ax;
        p_jh[i].vy += _dt * pji.ay;
        p_jh[i].vz += _dt * pji.az;
        
        // Additional Jacobi terms (second term from the right)
        const double rj2i = 1./(pji.x*pji.x + pji.y*pji.y + pji.z*pji.z);
        const double rji  = sqrt(rj2i);
        const double rj3iM = rji*rj2i*G*eta;
        const double prefac1 = _dt*rj3iM;
        p_jh[i].vx += prefac1*pji.x;
        p_jh[i].vy += prefac1*pji.y;
        p_jh[i].vz += prefac1*pji.z;
    }
}

/***************************** 
 * DKD Scheme                */

void reb_integrator_whfast_hj_kepler_step(const struct reb_simulation* const r, struct reb_particle* p_jh, const double _dt){
    const double m0 = r->particles[0].m;
    const double G = r->G;
    const size_t N = r->N;
    const size_t N_active = (r->N_active==SIZE_MAX || r->testparticle_type ==1)?N:r->N_active;
    double eta = m0;
    for (size_t i=1;i<N;i++){
        if (i<N_active){
            eta += p_jh[i].m;
        }
        reb_integrator_whfast_kepler_solver(&p_jh[i], eta*G, _dt, r);
    }
}

void reb_integrator_whfast_hj_com_step(const struct reb_simulation* const r, struct reb_particle* p_jh, const double _dt){
    (void)r;
    p_jh[0].x += _dt*p_jh[0].vx;
    p_jh[0].y += _dt*p_jh[0].vy;
    p_jh[0].z += _dt*p_jh[0].vz;
}

int reb_integrator_whfast_hj_init(struct reb_simulation* const r, struct reb_integrator_whfast_hj_state* whfast){
    const size_t N = r->N;
    if (whfast->N_allocated != N){
        whfast->N_allocated = N;
        whfast->p_jh = realloc(whfast->p_jh,sizeof(struct reb_particle)*N);
        r->did_modify_particles = 1;
    }
    return 0;
}

void reb_integrator_whfast_hj_from_inertial(struct reb_simulation* const r, struct reb_particle* p_jh){
    struct reb_particle* restrict const particles = r->particles;
    const size_t N = r->N;
    const size_t N_active = (r->N_active==SIZE_MAX || r->testparticle_type==1)?N:r->N_active;
    reb_transformations_inertial_to_jacobi_posvel(particles, p_jh, particles, N, N_active);
}

void reb_integrator_whfast_hj_to_inertial(struct reb_simulation* const r, struct reb_particle* p_jh){
    struct reb_particle* restrict const particles = r->particles;
    const size_t N = r->N;
    const size_t N_active = (r->N_active==SIZE_MAX || r->testparticle_type==1)?N:r->N_active;
    reb_transformations_jacobi_to_inertial_posvel(particles, p_jh, particles, N, N_active);
}

void reb_integrator_whfast_hj_step(struct reb_simulation* const r, void* state){
    struct reb_integrator_whfast_hj_state* whfast = state;
    const double dt = r->dt;
    if (reb_integrator_whfast_hj_init(r, whfast)){
        // Non recoverable error occurred.
        return;
    }

    reb_integrator_whfast_hj_from_inertial(r, whfast->p_jh);
    
    reb_integrator_whfast_hj_kepler_step(r, whfast->p_jh, r->dt/2.);    
    reb_integrator_whfast_hj_com_step(r, whfast->p_jh, r->dt/2.);

    reb_integrator_whfast_hj_to_inertial(r, whfast->p_jh);
    reb_simulation_update_acceleration(r);
    //   0 1 2 3 4 5
    // 0   x x x x x
    // 1 x   x x x x
    // 2 x x   x x x 
    // 3 x x x   x x
    // 4 x x x x   x
    // 5 x x x x x
    reb_integrator_whfast_hj_interaction_step(r, whfast->p_jh, dt);

    reb_integrator_whfast_hj_kepler_step(r, whfast->p_jh, r->dt/2.);    
    reb_integrator_whfast_hj_com_step(r, whfast->p_jh, r->dt/2.);
    reb_integrator_whfast_hj_to_inertial(r, whfast->p_jh);

    r->t+=r->dt;
    r->dt_last_done = r->dt;
}
