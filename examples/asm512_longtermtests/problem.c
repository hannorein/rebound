/**
 * Unit tests for WHFast512
 *
 * This file contains units tests for WHFast512.
 * Note that these are not run automatically 
 * because GitHub's CI does not support AVX5212.
 */

#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <sched.h>
#include <stdbool.h>

void gr_force(struct reb_simulation* r){
    double C2 = 10065.32 * 10065.32;
    struct reb_particle* particles = r->particles;
    const struct reb_particle source = particles[0];
    const double prefac1 = 6.*(r->G*source.m)*(r->G*source.m)/C2;
    for (size_t i=1; i<r->N; i++){
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
double gr_potential(struct reb_simulation* const r){
    const struct reb_particle* const particles = r->particles;
    double C2 = 10065.32 * 10065.32;
    const struct reb_particle source = particles[0];
	const double mu = r->G*source.m;
    const double prefac = 3.*mu*mu/C2;
    double H = 0.;

	for (size_t i=1;i<r->N;i++){
		struct reb_particle pi = particles[i];
        double dx = pi.x - source.x;
        double dy = pi.y - source.y;
        double dz = pi.z - source.z;
        double r2 = dx*dx + dy*dy + dz*dz;
        H -= prefac*pi.m/r2;
    }		
    return H;
}
struct reb_simulation* setup_sim(){
    struct reb_simulation* r = reb_simulation_create();
    // Setup constants
    r->dt = 5.0/365.25*2*M_PI;
    r->G = 1.;
    r->exact_finish_time = 0;
   
    reb_simulation_add_fmt(r, "solarsystem");

    return r;
}

extern void reb_integrator_asm512_kepler_step(struct reb_simulation* const r, int N_steps);
extern uint64_t reb_asm512_counter(struct reb_simulation* r, int test_p);


int main(int argc, char* argv[]) {
    struct reb_simulation* r_asm = setup_sim();
    struct reb_simulation* r_old = setup_sim();
    struct reb_simulation* r_whf = setup_sim();
    double E0 = reb_simulation_energy(r_asm)+gr_potential(r_asm);
    reb_simulation_set_integrator(r_asm, "asm512");
    struct reb_integrator_asm512_state* asm512 = r_asm->integrator.state;
    asm512->gr_potential = 1;
    asm512->concatenate_steps = 1e6;
    asm512->corrector = 17;
    reb_simulation_set_integrator(r_old, "whfast512");
    struct reb_integrator_whfast512_state* whfast512 = r_old->integrator.state;
    whfast512->gr_potential = 1;
    reb_simulation_set_integrator(r_whf, "whfast");
    struct reb_integrator_whfast_state* whfast = r_whf->integrator.state;
    whfast->safe_mode = 0;
    whfast->corrector = 17;
    r_whf->additional_forces = gr_force;
    for (double dT=1e1; r_asm->t<1e7*2*M_PI; dT=dT*1.05){
        reb_simulation_integrate(r_asm, r_asm->t+dT);
        reb_simulation_integrate(r_old, r_asm->t);
        reb_simulation_integrate(r_whf, r_asm->t);
        double E1_asm = reb_simulation_energy(r_asm)+gr_potential(r_asm);
        double E1_old = reb_simulation_energy(r_old)+gr_potential(r_old);
        double E1_whf = reb_simulation_energy(r_whf)+gr_potential(r_whf);
        printf("%e %e %e %e\n", r_asm->t, fabs((E0-E1_asm)/E0), fabs((E0-E1_old)/E0), fabs((E0-E1_whf)/E0));
    }
    reb_simulation_free(r_asm);
    reb_simulation_free(r_old);
    reb_simulation_free(r_whf);
    return 1;
}
