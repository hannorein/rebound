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
struct reb_simulation* setup_sim(int i){
    struct reb_simulation* r = reb_simulation_create();
    // Setup constants
    r->dt = 6.0/365.25*2*M_PI;
    r->G = 1.;
    r->exact_finish_time = 0;
   
    reb_simulation_add_fmt(r, "solarsystem");
    r->particles[3].x += 1e-15*i;
    reb_simulation_move_to_com(r);

    return r;
}

int main(int argc, char* argv[]) {
    int id = 0;
    if (argc>1){
        id = atoi(argv[1]);
    }
    struct reb_simulation* r_asm = setup_sim(id);
    double E0 = reb_simulation_energy(r_asm);
    reb_simulation_set_integrator(r_asm, "whfast512");
    struct reb_integrator_whfast512_state* whfast512 = r_asm->integrator.state;
    whfast512->gr_potential = 1;
    whfast512->concatenate_steps = 1e8;
    whfast512->corrector = 17;
    if (whfast512->gr_potential){
        E0 += gr_potential(r_asm);
    }
    char filename[1024];
    sprintf(filename, "/scratch/rein/whfast512_tests/out_rishit2_%03d.txt", id);
    for (double dT=1e1; r_asm->t<2*5e9*2*M_PI && r_asm->status <=0; dT=dT*1.05){
        reb_simulation_integrate(r_asm, r_asm->t+dT);
        double E1_asm = reb_simulation_energy(r_asm);
        char* mode = "a";
        if (dT==1e1) mode = "w";
        FILE* f;
        if (whfast512->gr_potential){
           E1_asm += gr_potential(r_asm);
        }
        f = fopen(filename,mode);
        fprintf(f,"%e ", r_asm->t);
        fprintf(f,"%e ", fabs((E0-E1_asm)/E0));
        fprintf(f,"\n");
        fclose(f);
    }
    reb_simulation_free(r_asm);
    return 1;
}
