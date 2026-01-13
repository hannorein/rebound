/**
 * Solar System with WHFast512
 *
 * This example integrates the Solar System using
 * the WHFast512 integrator. Note that you need a
 * CPU which support AVX512 instructions to run 
 * this example.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sched.h>
#include <stdbool.h>
#include "rebound.h"

// Initial conditions for the Solar System
// from NASA horizons
double all_ss_pos[9][3] = {
    {-0.008816286905115728, -0.0010954664916791675, 0.0002143249385447027},
    {-0.05942272929227954, -0.46308699693348293, -0.032897989948949075},
    {-0.7276101005375593, 0.006575003332463933, 0.041795901908847084},
    {-0.5789530452882667, -0.8361530119313055, 0.0002611520181901174},
    {-1.45202492400084, 0.827519404194876, 0.052981833432457694},
    {4.492983939852296, 2.0661626247490354, -0.10909246996001629},
    {8.4974210980544, -4.8620394993693585, -0.2537835862373596},
    {12.959111916929283, 14.760785302864473, -0.1130656917933948},
    {29.787987348666505, -2.51460654509393, -0.6347108842010732}
};

double all_ss_vel[9][3] = {
    {0.00014315746073017681, -0.0004912441820893999, 8.127678560998346e-07},
    {1.2978664284760637, -0.09524541469911743, -0.12677574364801253},
    {-0.019239782390457125, -1.1813975672919448, -0.01509392594251431},
    {0.8098712561282222, -0.5682496529341624, 2.6169897281383047e-05},
    {-0.37436417754222295, -0.6365841544564991, -0.004143932260467942},
    {-0.18818907783656452, 0.41919544951404614, 0.0024710497024424977},
    {0.14292308496870448, 0.2808676923735748, -0.010574288572728728},
    {-0.1734971049470612, 0.14019515029516152, 0.0027683484887051457},
    {0.014142947617173336, 0.18292110872737416, -0.004092845767710294}
};

double all_ss_mass[9] = {
    1.0,
    1.6601208254808336e-07,
    2.447838287784771e-06,
    3.0404326489511185e-06,
    3.2271560828978514e-07,
    0.0009547919099366768,
    0.0002858856700231729,
    4.366249613200406e-05,
    5.151383772628957e-05
};

// Implementation of the GR force for WHFast.
// (WHFast512 comes with built-in support) 
void gr_force(struct reb_simulation* r){
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


struct reb_simulation* run(int use_whfast512, int coordinates){
    
    struct reb_simulation* r = reb_simulation_create();
    // Setup constants
    r->dt = 5.0/365.25*2*M_PI; // 5 days
    r->G = 1.;
    r->exact_finish_time = 0;
    r->force_is_velocity_dependent = 0; 
    if (use_whfast512){ 
        r->integrator = REB_INTEGRATOR_WHFAST512;
        r->ri_whfast512.gr_potential = 1;
        r->ri_whfast512.coordinates = coordinates;
    }else{
        r->integrator = REB_INTEGRATOR_WHFAST;
        r->ri_whfast.coordinates = coordinates;
        r->gravity = REB_GRAVITY_JACOBI;
        r->ri_whfast.safe_mode = 0;
        r->additional_forces = gr_force;
    }
    
    // Initial conditions
    for (int i = 0; i < 9; i++) {
        struct reb_particle p = {
            .m = all_ss_mass[i],
            .x = all_ss_pos[i][0], .y = all_ss_pos[i][1], .z = all_ss_pos[i][2],
            .vx = all_ss_vel[i][0], .vy = all_ss_vel[i][1], .vz = all_ss_vel[i][2]
        };
        reb_simulation_add(r, p);
    }



    reb_simulation_move_to_com(r);
    reb_simulation_steps(r,1);

    struct timeval time_beginning;
    struct timeval time_end;
    gettimeofday(&time_beginning,NULL);
    int Nsteps = 1000000;
    reb_simulation_steps(r,Nsteps);

    gettimeofday(&time_end,NULL);
    double walltime = time_end.tv_sec-time_beginning.tv_sec+(time_end.tv_usec-time_beginning.tv_usec)/1e6;
    double tmax = Nsteps * r->dt;
    double gypday = 1e-9*(tmax/M_PI/2.)/walltime*86400;
    r->walltime = walltime;
    printf("walltime= %.2fs  (time required to integrate to 5 Gyr= %.2fdays)\n", walltime, 5./gypday);
    reb_simulation_synchronize(r);
    return r;
}

int main(int argc, char* argv[]) {
    //printf("Integrating for 1 Myr with WHFast512:\n");
    //double w1= run(1);
    struct reb_simulation* r1 = run(1,0);
    struct reb_simulation* r0 = run(0,0);
    struct reb_simulation* r2 = run(1,1); //DHC
    for (int i = 0; i < 9; i++) {
        printf("%d: %22.15e %22.15e %22.15e\n", i, r0->particles[i].x, r1->particles[i].x, (r0->particles[i].x - r1->particles[i].x)/r1->particles[i].x);
    }
    printf("------------------------\n");
    printf("Speedup (jac) %.4f\n", r0->walltime/r1->walltime);
    printf("Speedup (dhc) %.4f\n", r0->walltime/r2->walltime);

    
    return EXIT_SUCCESS;
}
