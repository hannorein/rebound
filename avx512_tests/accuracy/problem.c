#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sched.h>
#include <stdbool.h>
#include "rebound.h"

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_simulation_create();
    r->dt = 6.0/365.25*2*M_PI; // 6 days
    r->integrator = REB_INTEGRATOR_WHFAST512;
    r->exact_finish_time = 0;
    r->ri_whfast512.gr_potential = 0;
    r->ri_whfast512.coordinates = REB_WHFAST512_COORDINATES_JACOBI;
    r->ri_whfast512.concatenate_steps = 1234567;
    
    reb_simulation_add_fmt(r, "m", 1.0);
    const double a_init = 0.387098;
    for (int i=0;i<8;i++){
        reb_simulation_add_fmt(r, "m a e", 0.0, a_init, 0.1*(double)i);
    }

    double min_walltime = 1e300; 

    for (int i=0; i<10; i++){
        struct timeval time_beginning;
        struct timeval time_end;
        gettimeofday(&time_beginning,NULL);
        reb_simulation_steps(r,1);
        gettimeofday(&time_end,NULL);
        double walltime = time_end.tv_sec-time_beginning.tv_sec+(time_end.tv_usec-time_beginning.tv_usec)/1e6;
        if (walltime < min_walltime){
            min_walltime = walltime;
        }
    }
    reb_simulation_synchronize(r);

    const double n_init = 1.0/sqrt(a_init*a_init*a_init);
    const double P_init = 2.0*M_PI*sqrt(a_init*a_init*a_init);
    const double periods = fmod(r->t, P_init);
    //printf("%g %g %g\n", r->t, P_init, periods);
    printf("%.8f\t", min_walltime);
    for (int i=0;i<8;i++){
        struct reb_orbit o = reb_orbit_from_particle(r->G, r->particles[i+1], r->particles[0]);
        //printf("%.8f\n", o.M-n_init*periods);
        printf("%.8e\t%.8e\t%.8e\t%.8e\t", o.e-0.1*i, o.a-a_init, o.pomega, o.M-n_init*periods);
    }
    printf("\n");
    
    return EXIT_SUCCESS;
}
