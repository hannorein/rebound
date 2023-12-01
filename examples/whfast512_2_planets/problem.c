/**
 * Integrating a two planet system with WHFast512
 *
 * This example integrates four two-planets systems
 * using the WHFast512 integrator. Note that you need 
 * a CPU which support AVX512 instructions to run 
 * this example.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sched.h>
#include <stdbool.h>
#include "rebound.h"

// Initial conditions for the Sun, Mercury, and Venus
// from NASA horizons
double all_ss_pos[3][3] = {
    {-0.008816286905115728, -0.0010954664916791675, 0.0002143249385447027},
    {-0.05942272929227954, -0.46308699693348293, -0.032897989948949075},
    {-0.7276101005375593, 0.006575003332463933, 0.041795901908847084},
};

double all_ss_vel[3][3] = {
    {0.00014315746073017681, -0.0004912441820893999, 8.127678560998346e-07},
    {1.2978664284760637, -0.09524541469911743, -0.12677574364801253},
    {-0.019239782390457125, -1.1813975672919448, -0.01509392594251431},
};

double all_ss_mass[3] = {
    0.9999999999950272,
    1.6601208254808336e-07,
    2.447838287784771e-06,
};

struct reb_simulation* setup_single(){
    struct reb_simulation* r = reb_simulation_create();
    for (int i = 0; i < 3; i++){
        struct reb_particle p = {
            .m = all_ss_mass[i],
            .x = all_ss_pos[i][0], .y = all_ss_pos[i][1], .z = all_ss_pos[i][2],
            .vx = all_ss_vel[i][0], .vy = all_ss_vel[i][1], .vz = all_ss_vel[i][2]
        };
        reb_simulation_add(r, p);
    }
    reb_simulation_move_to_com(r);
    return r;
}

double run(int use_whfast512){
    struct timeval time_beginning;
    struct timeval time_end;
    double tmax = 2.*M_PI*1e5; // 100 kyr
    
    // We integrate four 2 planet systems in parallel. 
    // To do that, simply add all the planet to one simulation
    // so that the order of particles is:
    // Star 1
    //     Planet 1
    //     Planet 2 
    // Star 2
    //     Planet 1
    //     Planet 2 
    // Star 3
    //     Planet 1
    //     Planet 2 
    // Star 4
    //     Planet 1
    //     Planet 2 

    if (use_whfast512){ 
        gettimeofday(&time_beginning,NULL);
        // One simulation with all 4x2 = 8 planets.
        struct reb_simulation* r = reb_simulation_create();
        r->exact_finish_time = 0;
        r->dt = 5.0/365.25*2*M_PI; // 5 days
        r->G = 1.;
        r->force_is_velocity_dependent = 0; 
        // Tell WHFast512 how many systems we are integrating in parallel.
        // This parameter can be either 1, 2, or 4.
        r->ri_whfast512.N_systems = 4;
        for (int s = 0; s < 4; s++){
            struct reb_simulation* r_single = setup_single();
            // We're adding a small perturbation to each simulation so they are
            // not all exactly the same. In principle the simulations can be
            // completely different, the only thing that needs to be the same
            // is the timestep.
            r_single->particles[1].x += 1e-14*s; 
            for (int i=0; i<r_single->N; i++){
                reb_simulation_add(r, r_single->particles[i]);
            }
            reb_simulation_free(r_single);
        }
        r->integrator = REB_INTEGRATOR_WHFAST512;
        int err = reb_simulation_integrate(r,  tmax);
        if (err>0){
            printf("An error occured during the integration.\n");
            exit(EXIT_FAILURE);
        }
        reb_simulation_free(r);
        gettimeofday(&time_end,NULL);
    }else{
        // Without WHFast512 we need to integrate 4 simulations one after the other
        gettimeofday(&time_beginning,NULL);
        for (int s = 0; s < 4; s++){
            struct reb_simulation* r = setup_single();
            r->exact_finish_time = 0;
            r->dt = 5.0/365.25*2*M_PI; // 5 days
            r->G = 1.;
            r->force_is_velocity_dependent = 0; 
            r->integrator = REB_INTEGRATOR_WHFAST;
            r->ri_whfast.coordinates = REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC;
            r->particles[1].x += 1e-14*s; 
            int err = reb_simulation_integrate(r,  tmax);
            if (err>0){
                printf("An error occured during the integration.\n");
                exit(EXIT_FAILURE);
            }
            reb_simulation_free(r);
        }
        gettimeofday(&time_end,NULL);
    }
   
    double walltime = time_end.tv_sec-time_beginning.tv_sec+(time_end.tv_usec-time_beginning.tv_usec)/1e6;
    double gypday = 1e-9*(tmax/M_PI/2.)/walltime*86400;
    printf("walltime= %.2fs  (time required to integrate to 5 Gyr= %.2fdays)\n", walltime, 5./gypday);
    return walltime;
}

int main(int argc, char* argv[]) {
    printf("Integrating for 100 kyr with WHFast512:\n");
    double w1= run(1);
    printf("Integrating for 100 kyr with WHFast:\n");
    double w0= run(0);
    printf("\nSpeedup: %.2fx\n", w0/w1);
    return EXIT_SUCCESS;
}
