/**
 * Apophis
 * 
 * This problem uses TES to integrate the Apophis encounter with the Earth in 2029. 
 * Toy model consisting of Sun, Earth and Apophis is used.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

double pos[3][3] = 
{
{2.046127869986152144e-06,	2.237089773829956975e-06,	-1.666185482284776807e-10},
{-6.812497776258548132e-01,	-7.448297505282599484e-01,	5.547495373765781696e-05},
{7.547693442365438488e-01,	-1.898579598772653920e-03,	1.866781385938447516e-02},
};

double vel[3][3] = 
{
    {-3.729016715508550764e-08,	3.509207245432022034e-08,	-2.679097518437730580e-12}, 
    {1.241560630431407311e-02,	-1.168375980143308152e-02,	8.919943936329344472e-07}, 
    {1.722330612025966146e-03,	2.138934116281841075e-02,	-1.089248121544092191e-03},
};

double mass[3] =
{
    2.959122082855910945e-04, 8.887697821383033199e-10, 4.016599676151094958e-24,
};

double final_pos[3][3] = 
{
    {0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00},
    {-9.9604746697684021e-01,  3.5311915613020404e-03,  -1.2054180564475472e-06},
    {-8.1089946547081804e-01, -5.4094730893500276e-01,  6.8972157890442951e-03},
};


void heartbeat(struct reb_simulation* r);
double e_init;
double e_final;
double tmax;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    // Setup constants
    double years = 100;
    r->dt             = 365.24141698304345/100;
    tmax            = years*365.24141698304345;  
    r->G            = 1.4881806877180788e-34;        // in AU^3 / kg / day^2.
    r->integrator        = REB_INTEGRATOR_TES;
    r->heartbeat        = heartbeat;
    r->exact_finish_time = 1; // Finish exactly at tmax in reb_integrate(). Default is already 1.
    r->ri_tes.dq_max = 1e-3;
    r->ri_tes.recti_per_orbit = 1.61803398875;
    r->ri_tes.epsilon = 1e-6;

    // Initial conditions
    for (int i=0;i<3;i++){
        struct reb_particle p = {0};
        p.x  = pos[i][0];         p.y  = pos[i][1];         p.z  = pos[i][2];
        p.vx = vel[i][0];         p.vy = vel[i][1];         p.vz = vel[i][2];
        p.m  = mass[i]/r->G;
        reb_add(r, p); 
    }
    reb_move_to_com(r);
    e_init = reb_tools_energy(r);
    system("rm -f energy.txt");
    reb_integrate(r, tmax);
    e_final = reb_tools_energy(r);

    printf("\n\nRelative energy error: %.2E \n", fabs((e_final-e_init)/e_init));

    double max_pos_err = 0;
    for(uint32_t i = 0; i < 3; i ++)
    {
        double dx = fabs(r->particles[i].x-final_pos[i][0]);
        double dy = fabs(r->particles[i].y-final_pos[i][1]);
        double dz = fabs(r->particles[i].z-final_pos[i][2]);
        double err = dx+dy+dz;

        max_pos_err = err > max_pos_err ?  err : max_pos_err;
    }

    printf("\nMax position error: %.3e AU", max_pos_err);

    reb_free_simulation(r);

}

void heartbeat(struct reb_simulation* r){
    if (reb_output_check(r, 365)){
        reb_output_timing(r, tmax);
        reb_integrator_synchronize(r);
        FILE* f = fopen("energy.txt","a");
        double e = reb_tools_energy(r);
        fprintf(f,"%e %e\n",r->t, fabs((e-e_init)/e_init));
        fclose(f);
    }
}

