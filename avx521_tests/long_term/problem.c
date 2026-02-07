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
const double all_ss_pos[9][3] = {
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

const double all_ss_vel[9][3] = {
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

const double all_ss_mass[9] = {
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

const double C2 = 10065.32 * 10065.32;
static double energy(struct reb_simulation* const sim){
    const struct reb_particle* const particles = sim->particles;
	const double G = sim->G;
    const struct reb_particle source = particles[0];
	const double mu = G*source.m;
    const double prefac = 3.*mu*mu/C2;
    double H = 0.;

	for (int i=1;i<sim->N;i++){
		struct reb_particle pi = particles[i];
        double dx = pi.x - source.x;
        double dy = pi.y - source.y;
        double dz = pi.z - source.z;
        double r2 = dx*dx + dy*dy + dz*dz;
        H -= prefac*pi.m/r2;
    }		
	
    return reb_simulation_energy(sim)+H;
}

int main(int argc, char* argv[]) {
    if (argc!=2){
        printf("Usage: rebound [id]\n");
        exit(0);
    }
    int id = atoi(argv[1]);
    
    struct reb_simulation* r = reb_simulation_create();
    r->dt = 1.0/365.25*2*M_PI;
    r->G = 1.;
    r->exact_finish_time = 0;
    r->integrator = REB_INTEGRATOR_WHFAST512;
    r->ri_whfast512.gr_potential = 1;
    r->ri_whfast512.keep_unsynchronized = 1;
    r->ri_whfast512.coordinates = REB_WHFAST512_COORDINATES_JACOBI;
    for (int i = 0; i < 9; i++) {
        struct reb_particle p = {
            .m = all_ss_mass[i],
            .x = all_ss_pos[i][0], .y = all_ss_pos[i][1], .z = all_ss_pos[i][2],
            .vx = all_ss_vel[i][0], .vy = all_ss_vel[i][1], .vz = all_ss_vel[i][2]
        };
        reb_simulation_add(r, p);
    }
    r->particles[1].x += 1e-11*(double)id; // ~meter
    reb_simulation_move_to_com(r);

    double E0 = energy(r);
    struct timeval time_beginning;
    struct timeval time_end;
    gettimeofday(&time_beginning,NULL);
    r->ri_whfast512.concatenate_steps = 1e5;

    double next_steps = 1;
    char filename_log[1024];
    sprintf(filename_log, "/scratch/rein/whfast512_tests/log_%05d.txt", id);
    remove(filename_log);
    char filename_archive[1024];
    sprintf(filename_archive, "/scratch/rein/whfast512_tests/archive_%05d.bin", id);
    remove(filename_archive);
        
    reb_simulation_save_to_file(r, filename_archive);
    while (r->t/2.0/M_PI < 5e9){
        reb_simulation_steps(r,(int)next_steps);
        reb_simulation_synchronize(r);
        double maxe = 0.0;
        for (int i=1;i<r->N;i++){
            struct reb_orbit o = reb_orbit_from_particle(r->G, r->particles[i], r->particles[0]);
            if (o.e>maxe){
                maxe = o.e;
            }
        }
        //r->t += r->dt * 1e5 * (int)next_steps;

        reb_simulation_save_to_file(r, filename_archive);

        gettimeofday(&time_end,NULL);
        double walltime = time_end.tv_sec-time_beginning.tv_sec+(time_end.tv_usec-time_beginning.tv_usec)/1e6;
        double E1 = energy(r);
        double E = fabs((E1-E0)/E0);

        FILE* f = fopen(filename_log, "a");
        fprintf(f,"%.16e\t%.16e\t%.16e\t%.16e\n", r->t, E, walltime, maxe);
        fclose(f);
        
        next_steps *= 1.005;
    }

    reb_simulation_free(r);
    
    return EXIT_SUCCESS;
}
