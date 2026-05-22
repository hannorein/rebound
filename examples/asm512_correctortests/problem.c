#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <sched.h>
#include <stdbool.h>
#include <string.h>

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

void setup_sim(char* integrator, int corrector){
    printf("integrator=%s, corrector=%d: \t", integrator, corrector);
    struct reb_simulation* r;
    for(int gr = 0; gr<2; gr++){
        r = reb_simulation_create();
        // Setup constants
        r->dt = 6.0/365.25*2*M_PI;
        r->G = 1.;
        r->exact_finish_time = 0;
        reb_simulation_add_fmt(r, "solarsystem");
        reb_simulation_set_integrator(r, integrator);
        if (strcmp(integrator,"asm512")==0){
            struct reb_integrator_asm512_state* asm512 = r->integrator.state;
            asm512->gr_potential = 0;
            asm512->concatenate_steps = 1e4;
            asm512->corrector = corrector;
            asm512->gr_potential = gr;
        }
        if (strcmp(integrator,"whfast512")==0){
            struct reb_integrator_whfast512_state* whfast512 = r->integrator.state;
            whfast512->gr_potential = gr;
        }
        if (strcmp(integrator,"whfast")==0){
            struct reb_integrator_whfast_state* whfast = r->integrator.state;
            whfast->safe_mode = 0;
            whfast->corrector = corrector;
            if (gr) r->additional_forces = gr_force;
        }
        double E0 = reb_simulation_energy(r);
        if (gr) E0+= gr_potential(r);
        reb_simulation_integrate(r, 1e2*M_PI*2);
        double E1 = reb_simulation_energy(r);
        if (gr) E1+= gr_potential(r);
        if (!gr) printf("no-");
        printf("gr= %.16e ", fabs((E0-E1)/E0));
        reb_simulation_free(r);
    }
    printf("\n");
}



extern void reb_integrator_asm512_kepler_step(struct reb_simulation* const r, int N_steps);
extern uint64_t reb_asm512_counter(struct reb_simulation* r, int test_p);

int main(int argc, char* argv[]) {
    setup_sim("asm512", 0);
    setup_sim("asm512", 17);
    setup_sim("whfast", 0);
    setup_sim("whfast", 17);
    setup_sim("whfast512", 0);

}
