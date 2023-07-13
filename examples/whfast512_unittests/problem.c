/**
 * Unit tests for WHFast512
 *
 * This file contains units tests for WHFast512.
 * Note that these are not run automatically 
 * because GitHub's CI does not support AVX5212.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include <sched.h>
#include <stdbool.h>
#include "rebound.h"

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
    0.9999999999950272,
    1.6601208254808336e-07,
    2.447838287784771e-06,
    3.0404326489511185e-06,
    3.2271560828978514e-07,
    0.0009547919099366768,
    0.0002858856700231729,
    4.366249613200406e-05,
    5.151383772628957e-05
};

struct reb_simulation* setup_sim(){
    struct reb_simulation* r = reb_create_simulation();
    // Setup constants
    r->dt = 4.0/365.25*2*M_PI; //6 days
    r->G = 1.;
    r->exact_finish_time = 0;
    r->force_is_velocity_dependent = 0; 
    
    // Initial conditions
    for (int i = 0; i < 9; i++) {
        struct reb_particle p = {0};
        p.x = all_ss_pos[i][0];
        p.y = all_ss_pos[i][1];
        p.z = all_ss_pos[i][2];
        p.vx = all_ss_vel[i][0];
        p.vy = all_ss_vel[i][1];
        p.vz = all_ss_vel[i][2];
        p.m = all_ss_mass[i];
        reb_add(r, p);
    }
    reb_move_to_com(r);
    return r;
}

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

int test_basic(){
    struct reb_simulation* r = setup_sim();
    struct reb_simulation* r512 = reb_copy_simulation(r);
     
    r512->integrator = REB_INTEGRATOR_WHFAST512;
    r512->ri_whfast512.gr_potential = 0;
    r->integrator = REB_INTEGRATOR_WHFAST;
    r->ri_whfast.coordinates = REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC;
    r->ri_whfast.safe_mode = 0;

    double tmax = 1e2;
    if (reb_integrate(r, tmax)>0) return 0;
    if (reb_integrate(r512, tmax)>0) return 0;

    for (int i=0;i<r->N;i++){
        if (fabs(r->particles[i].x - r512->particles[i].x)>1e-11){
            printf("Accuracy not met in GR test.\n");
            printf("%.16e\n",fabs(r->particles[i].x - r512->particles[i].x));
            return 0;
        }
    }

    reb_free_simulation(r);
    reb_free_simulation(r512);
    return 1;
}


int test_com(){
    struct reb_simulation* r512 = setup_sim();
     
    r512->integrator = REB_INTEGRATOR_WHFAST512;
    r512->ri_whfast512.gr_potential = 0;

    double tmax = 1e5;
    if (reb_integrate(r512, tmax)>0) return 0;
    struct reb_particle com = reb_get_com(r512);
    assert(fabs(com.x)<1e-14);
    assert(fabs(com.y)<1e-14);
    assert(fabs(com.z)<1e-14);

    reb_free_simulation(r512);
    return 1;
}

int test_twobody(){
    struct reb_simulation* r512 = reb_create_simulation();
    r512->exact_finish_time = 0;
    r512->integrator = REB_INTEGRATOR_WHFAST512;
    r512->ri_whfast512.gr_potential = 0;
    reb_add_fmt(r512, "m", 1.0);
    reb_add_fmt(r512, "a", 1.0);

    double tmax = 10.*M_PI*2.;
    r512->dt = tmax / 128; 
    if (reb_integrate(r512, tmax)>0) return 0;
    assert(fabs(r512->particles[1].x-1.0)<2e-15);
    assert(fabs(r512->particles[1].y)<2e-13);
    assert(fabs(r512->particles[1].z)==0.0);

    reb_free_simulation(r512);
    return 1;
}

int test_gr(){
    struct reb_simulation* r = setup_sim();
    struct reb_simulation* r512 = reb_copy_simulation(r);
     
    r512->integrator = REB_INTEGRATOR_WHFAST512;
    r512->ri_whfast512.gr_potential = 1;
    r->integrator = REB_INTEGRATOR_WHFAST;
    r->ri_whfast.coordinates = REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC;
    r->ri_whfast.safe_mode = 0;
    r->additional_forces = gr_force;

    double tmax = 1e2;
    if (reb_integrate(r, tmax)>0) return 0;
    if (reb_integrate(r512, tmax)>0) return 0;

    for (int i=0;i<r->N;i++){
        if (fabs(r->particles[i].x - r512->particles[i].x)>1e-9){
            printf("Accuracy not met in GR test.\n");
            printf("%.16e\n",fabs(r->particles[i].x - r512->particles[i].x));
            return 0;
        }
    }

    reb_free_simulation(r);
    reb_free_simulation(r512);
    return 1;
}

int test_restart(){
    struct reb_simulation* r512 = setup_sim();
    r512->integrator = REB_INTEGRATOR_WHFAST512;
    r512->exact_finish_time = 0;
    r512->ri_whfast512.gr_potential = 1;
    r512->ri_whfast512.keep_unsynchronized = 1;
    
    double tmax = 1e2;
    double tmaxfinal = 4.*tmax;

    struct reb_simulation* r512c = reb_copy_simulation(r512);
    if (reb_integrate(r512c, tmaxfinal)>0) return 0;

    if (reb_integrate(r512, tmax)>0) return 0;
    if (reb_integrate(r512, 2.*tmax)>0) return 0;
    reb_output_binary(r512, "test.bin");
    if (reb_integrate(r512, 3.*tmax)>0) return 0;
    if (reb_integrate(r512, tmaxfinal)>0) return 0;
    
    struct reb_simulation* r512c2 = reb_create_simulation_from_binary("test.bin");
    if (r512c2 == NULL) return 0;
    if (reb_integrate(r512c2, 3.*tmax)>0) return 0;
    if (reb_integrate(r512c2, tmaxfinal)>0) return 0;

    for (int i=0;i<r512->N;i++){
        assert(r512->t == r512c->t);
        assert(r512->particles[i].m == r512c->particles[i].m);
        assert(r512->particles[i].x == r512c->particles[i].x);
        assert(r512->particles[i].vx == r512c->particles[i].vx);
        assert(r512->t == r512c2->t);
        assert(r512->particles[i].x == r512c2->particles[i].x);
        assert(r512->particles[i].m == r512c2->particles[i].m);
        assert(r512->particles[i].vx == r512c2->particles[i].vx);
    }

    reb_free_simulation(r512);
    reb_free_simulation(r512c);
    reb_free_simulation(r512c2);
    return 1;
}

int main(int argc, char* argv[]) {
    assert(test_basic());
    assert(test_restart());
    assert(test_com());
    assert(test_twobody());
    assert(test_gr());
    printf("All tests passed.\n");
}
