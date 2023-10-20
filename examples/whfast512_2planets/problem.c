/**
 * 4 Planets with WHFast512
 *
 * This example integrates two 4 planet systems using
 * the WHFast512 integrator. Note that you need a
 * CPU which support AVX512 instructions to run 
 * this example.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sched.h>
#include <stdbool.h>
#include "rebound.h"
#include "integrator_whfast512.h"

// Initial conditions for the Solar System
// from NASA horizons

double all_ss_pos[9][3] = {
    {-0.008816286905115728, -0.0010954664916791675, 0.0002143249385447027},
    {-0.05942272929227954, -0.46308699693348293, -0.032897989948949075},
    {-0.7276101005375593, 0.006575003332463933, 0.041795901908847084},
    {-0.5789530452882667, -0.8361530119313055, 0.0002611520181901174},
    {-1.45202492400084, 0.827519404194876, 0.052981833432457694},
    {-0.05942272929227954, -0.46308699693348293, -0.032897989948949075},
    {-0.7276101005375593, 0.006575003332463933, 0.041795901908847084},
    {-0.5789530452882667, -0.8361530119313055, 0.0002611520181901174},
    {-1.45202492400084, 0.827519404194876, 0.052981833432457694}
};

double all_ss_vel[9][3] = {
    {0.00014315746073017681, -0.0004912441820893999, 8.127678560998346e-07},
    {1.2978664284760637, -0.09524541469911743, -0.12677574364801253},
    {-0.019239782390457125, -1.1813975672919448, -0.01509392594251431},
    {0.8098712561282222, -0.5682496529341624, 2.6169897281383047e-05},
    {-0.37436417754222295, -0.6365841544564991, -0.004143932260467942},
    {1.2978664284760637, -0.09524541469911743, -0.12677574364801253},
    {-0.019239782390457125, -1.1813975672919448, -0.01509392594251431},
    {0.8098712561282222, -0.5682496529341624, 2.6169897281383047e-05},
    {-0.37436417754222295, -0.6365841544564991, -0.004143932260467942}
};

double all_ss_mass[9] = {
    0.9999999999950272,
    1.6601208254808336e-07,
    2.447838287784771e-06,
    3.0404326489511185e-06,
    3.2271560828978514e-07,
    1.6601208254808336e-07,
    2.447838287784771e-06,
    3.0404326489511185e-06,
    3.2271560828978514e-07
};


double all_2_pos1[3][3] = {
    {-0.008816286905115728, -0.0010954664916791675, 0.0002143249385447027},
    {-0.05942272929227954, -0.46308699693348293, -0.032897989948949075},
    {-0.7276101005375593, 0.006575003332463933, 0.041795901908847084}
};

double all_2_vel1[3][3] = {
    {0.00014315746073017681, -0.0004912441820893999, 8.127678560998346e-07},
    {1.2978664284760637, -0.09524541469911743, -0.12677574364801253},
    {-0.019239782390457125, -1.1813975672919448, -0.01509392594251431}
};

double all_2_mass1[3] = {
    0.9999999999950272,
    1.6601208254808336e-07,
    2.447838287784771e-06
};

double all_2_pos2[3][3] = {
    {-0.008816286905115728, -0.0010954664916791675, 0.0002143249385447027},
    {-0.5789530452882667, -0.8361530119313055, 0.0002611520181901174},
    {-1.45202492400084, 0.827519404194876, 0.052981833432457694}
};

double all_2_vel2[3][3] = {
    {0.00014315746073017681, -0.0004912441820893999, 8.127678560998346e-07},
    {0.8098712561282222, -0.5682496529341624, 2.6169897281383047e-05},
    {-0.37436417754222295, -0.6365841544564991, -0.004143932260467942}
};

double all_2_mass2[3] = {
    0.9999999999950272,
    3.0404326489511185e-06,
    3.2271560828978514e-07
};

double all_2_pos3[3][3] = {
    {-0.008816286905115728, -0.0010954664916791675, 0.0002143249385447027},
    {4.492983939852296, 2.0661626247490354, -0.10909246996001629},
    {8.4974210980544, -4.8620394993693585, -0.2537835862373596}
};

double all_2_vel3[3][3] = {
    {0.00014315746073017681, -0.0004912441820893999, 8.127678560998346e-07},
    {-0.18818907783656452, 0.41919544951404614, 0.0024710497024424977},
    {0.14292308496870448, 0.2808676923735748, -0.010574288572728728}
};

double all_2_mass3[3] = {
    0.9999999999950272,
    0.0009547919099366768,
    0.0002858856700231729
};

double all_2_pos4[3][3] = {
    {-0.008816286905115728, -0.0010954664916791675, 0.0002143249385447027},
    {12.959111916929283, 14.760785302864473, -0.1130656917933948},
    {29.787987348666505, -2.51460654509393, -0.6347108842010732}
};

double all_2_vel4[3][3] = {
    {0.00014315746073017681, -0.0004912441820893999, 8.127678560998346e-07},
    {-0.1734971049470612, 0.14019515029516152, 0.0027683484887051457},
    {0.014142947617173336, 0.18292110872737416, -0.004092845767710294}
};

double all_2_mass4[3] = {
    0.9999999999950272,
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


double run(int use_whfast512){
    struct timeval time_beginning;
    struct timeval time_end;

    struct reb_simulation* r4 = reb_create_simulation();
    // Setup constants
    r4->dt = 5.0/365.25*2*M_PI; // 5 days
    r4->G = 1.;
    r4->exact_finish_time = 0;
    r4->force_is_velocity_dependent = 0; 
    r4->N = 0;
    if (use_whfast512){ 
        r4->integrator = REB_INTEGRATOR_WHFAST512;
        r4->ri_whfast512.gr_potential = 1;
    }else{
        r4->integrator = REB_INTEGRATOR_WHFAST;
        r4->ri_whfast.coordinates = REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC;
        r4->ri_whfast.safe_mode = 0;
        r4->additional_forces = gr_force;
    }
    
    // Initial conditions
    for (int i = 0; i < 3; i++) {
        struct reb_particle p = {
            .m = all_2_mass4[i],
            .x = all_2_pos4[i][0], .y = all_2_pos4[i][1], .z = all_2_pos4[i][2],
            .vx = all_2_vel4[i][0], .vy = all_2_vel4[i][1], .vz = all_2_vel4[i][2]
        };
        reb_add(r4, p);
    }

    reb_move_to_com(r4);
    
    struct reb_simulation* r3 = reb_create_simulation();
    // Setup constants
    r3->dt = 5.0/365.25*2*M_PI; // 5 days
    r3->G = 1.;
    r3->exact_finish_time = 0;
    r3->force_is_velocity_dependent = 0; 
    r3->N = 0;
    if (use_whfast512){ 
        r3->integrator = REB_INTEGRATOR_WHFAST512;
        r3->ri_whfast512.gr_potential = 1;
    }else{
        r3->integrator = REB_INTEGRATOR_WHFAST;
        r3->ri_whfast.coordinates = REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC;
        r3->ri_whfast.safe_mode = 0;
        r3->additional_forces = gr_force;
    }
    
    // Initial conditions
    for (int i = 0; i < 3; i++) {
        struct reb_particle p = {
            .m = all_2_mass3[i],
            .x = all_2_pos3[i][0], .y = all_2_pos3[i][1], .z = all_2_pos3[i][2],
            .vx = all_2_vel3[i][0], .vy = all_2_vel3[i][1], .vz = all_2_vel3[i][2]
        };
        reb_add(r3, p);
    }

    reb_move_to_com(r3);


    struct reb_simulation* r2 = reb_create_simulation();
    // Setup constants
    r2->dt = 5.0/365.25*2*M_PI; // 5 days
    r2->G = 1.;
    r2->exact_finish_time = 0;
    r2->force_is_velocity_dependent = 0; 
    r2->N = 0;
    if (use_whfast512){ 
        r2->integrator = REB_INTEGRATOR_WHFAST512;
        r2->ri_whfast512.gr_potential = 1;
    }else{
        r2->integrator = REB_INTEGRATOR_WHFAST;
        r2->ri_whfast.coordinates = REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC;
        r2->ri_whfast.safe_mode = 0;
        r2->additional_forces = gr_force;
    }
    
    // Initial conditions
    for (int i = 0; i < 3; i++) {
        struct reb_particle p = {
            .m = all_2_mass2[i],
            .x = all_2_pos2[i][0], .y = all_2_pos2[i][1], .z = all_2_pos2[i][2],
            .vx = all_2_vel2[i][0], .vy = all_2_vel2[i][1], .vz = all_2_vel2[i][2]
        };
        reb_add(r2, p);
    }

    reb_move_to_com(r2);

    struct reb_simulation* r1 = reb_create_simulation();
    // Setup constants
    r1->dt = 5.0/365.25*2*M_PI; // 5 days
    r1->G = 1.;
    r1->exact_finish_time = 0;
    r1->force_is_velocity_dependent = 0; 
    r1->N = 0;
    if (use_whfast512){ 
        r1->integrator = REB_INTEGRATOR_WHFAST512;
        r1->ri_whfast512.gr_potential = 1;
    }else{
        r1->integrator = REB_INTEGRATOR_WHFAST;
        r1->ri_whfast.coordinates = REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC;
        r1->ri_whfast.safe_mode = 0;
        r1->additional_forces = gr_force;
    }
    
    // Initial conditions
    for (int i = 0; i < 3; i++) {
        struct reb_particle p = {
            .m = all_2_mass1[i],
            .x = all_2_pos1[i][0], .y = all_2_pos1[i][1], .z = all_2_pos1[i][2],
            .vx = all_2_vel1[i][0], .vy = all_2_vel1[i][1], .vz = all_2_vel1[i][2]
        };
        reb_add(r1, p);
    }

    reb_move_to_com(r1);

    struct reb_simulation *r = combine_simulations(r1, r2, r3, r4);

    gettimeofday(&time_beginning,NULL);
    double tmax = 3000000*2*M_PI; // 15 days
    int err = 0;
    if (use_whfast512){
        err = reb_integrate(r,  tmax);
        for (int i=0; i<r->N; i++){
        struct reb_particle p = r->particles[i];
        printf("%f %f %f\n", p.x, p.y, p.z);
        }
    }else{
        err = reb_integrate(r1, tmax);
        if (err>0){
            printf("An error occured during the integration.\n");
            exit(EXIT_FAILURE);
        }
        for (int i=0; i<r1->N; i++){
        struct reb_particle p = r1->particles[i];
        printf("%f %f %f\n", p.x, p.y, p.z);
        }
        err = reb_integrate(r2, tmax);
        if (err>0){
            printf("An error occured during the integration.\n");
            exit(EXIT_FAILURE);
        }
        for (int i=0; i<r2->N; i++){
        struct reb_particle p = r2->particles[i];
        printf("%f %f %f\n", p.x, p.y, p.z);
        }
        err = reb_integrate(r3, tmax);
        if (err>0){
            printf("An error occured during the integration.\n");
            exit(EXIT_FAILURE);
        }
        for (int i=0; i<r3->N; i++){
        struct reb_particle p = r3->particles[i];
        printf("%f %f %f\n", p.x, p.y, p.z);
        }

        err = reb_integrate(r4, tmax);
        for (int i=0; i<r4->N; i++){
        struct reb_particle p = r4->particles[i];
        printf("%f %f %f\n", p.x, p.y, p.z);
        }
    }
    if (err>0){
        printf("An error occured during the integration.\n");
        exit(EXIT_FAILURE);
    }
    gettimeofday(&time_end,NULL);

    double walltime = time_end.tv_sec-time_beginning.tv_sec+(time_end.tv_usec-time_beginning.tv_usec)/1e6;
    double gypday = 1e-9*(tmax/M_PI/2.)/walltime*86400;
    printf("walltime= %.2fs  (time required to integrate to 5 Gyr= %.2fdays)\n", walltime, 5./gypday);
    return walltime;
}

int main(int argc, char* argv[]) {
    printf("Integrating for 1 Myrs with WHFast512:\n");
    double w1= run(1);
    printf("Integrating for 1 Myrs with WHFast:\n");
    double w0= run(0);
    printf("\nSpeedup: %.2fx\n", w0/w1);
    return EXIT_SUCCESS;
}
