#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "../minunit.h"

double ss_pos[6][3] =
    {
     {-4.06428567034226e-3, -6.08813756435987e-3, -1.66162304225834e-6}, // Sun
     {+3.40546614227466e+0, +3.62978190075864e+0, +3.42386261766577e-2}, // Jupiter
     {+6.60801554403466e+0, +6.38084674585064e+0, -1.36145963724542e-1}, // Saturn
     {+1.11636331405597e+1, +1.60373479057256e+1, +3.61783279369958e-1}, // Uranus
     {-3.01777243405203e+1, +1.91155314998064e+0, -1.53887595621042e-1}, // Neptune
     {-2.13858977531573e+1, +3.20719104739886e+1, +2.49245689556096e+0}  // Pluto
};
double ss_vel[6][3] =
    {
     {+6.69048890636161e-6, -6.33922479583593e-6, -3.13202145590767e-9}, // Sun
     {-5.59797969310664e-3, +5.51815399480116e-3, -2.66711392865591e-6}, // Jupiter
     {-4.17354020307064e-3, +3.99723751748116e-3, +1.67206320571441e-5}, // Saturn
     {-3.25884806151064e-3, +2.06438412905916e-3, -2.17699042180559e-5}, // Uranus
     {-2.17471785045538e-4, -3.11361111025884e-3, +3.58344705491441e-5}, // Neptune
     {-1.76936577252484e-3, -2.06720938381724e-3, +6.58091931493844e-4}  // Pluto
};

double ss_mass[6] =
    {
     1.00000597682, // Sun + inner planets
     1. / 1047.355, // Jupiter
     1. / 3501.6,   // Saturn
     1. / 22869.,   // Uranus
     1. / 19314.,   // Neptune
     7.4074074e-09  // Pluto
};

double tmax = 7.3e5;

int tests_run = 0;


static char * test_oss(){
	struct reb_simulation* r = reb_create_simulation();
	// Setup constants
	const double k = 0.01720209895; // Gaussian constant
	r->dt = 40;			// in days
	r->G = k * k;			// These are the same units as used by the mercury6 code.
	r->ri_whfast.safe_mode = 0;     // Turn of safe mode. Need to call integrator_synchronize() before outputs.
	r->ri_whfast.corrector = 11;    // Turn on symplectic correctors (11th order).

	// Setup callbacks:
	r->force_is_velocity_dependent = 0; // Force only depends on positions.
	r->integrator = REB_INTEGRATOR_WHFAST;
	//r->integrator	= REB_INTEGRATOR_IAS15;
	//r->integrator	= REB_INTEGRATOR_WH;

	// Initial conditions
	for (int i = 0; i < 6; i++) {
		struct reb_particle p = {0};
		p.x = ss_pos[i][0];
		p.y = ss_pos[i][1];
		p.z = ss_pos[i][2];
		p.vx = ss_vel[i][0];
		p.vy = ss_vel[i][1];
		p.vz = ss_vel[i][2];
		p.m = ss_mass[i];
		reb_add(r, p);
	}
	if (r->integrator == REB_INTEGRATOR_WH) {
		struct reb_particle* const particles = r->particles;
		// Move to heliocentric frame (required by WH integrator)
		for (int i = 1; i < r->N; i++) {
			particles[i].x -= particles[0].x;
			particles[i].y -= particles[0].y;
			particles[i].z -= particles[0].z;
			particles[i].vx -= particles[0].vx;
			particles[i].vy -= particles[0].vy;
			particles[i].vz -= particles[0].vz;
		}
		particles[0].x = 0;
		particles[0].y = 0;
		particles[0].z = 0;
		particles[0].vx = 0;
		particles[0].vy = 0;
		particles[0].vz = 0;
	} else {
		reb_move_to_com(r);
	}

	r->N_active = r->N - 1; // Pluto is treated as a test-particle.

	double e_initial = reb_tools_energy(r);

	reb_integrate(r, tmax);

	double e_final = reb_tools_energy(r);
    
    // Relative energy error: 
    double e_err = fabs((e_final - e_initial) / e_initial);
    mu_assert_almost_equal("Energy error too large.",e_err,0.,1e-8);
    return 0;
}

static char * all_tests() {
    mu_run_test(test_oss);
    return 0;
}

int main(int argc, char* argv[]) {
    char *result = all_tests();
    if (result !=0) {
        printf("%s\n", result);
    }else{
        printf("ALL TESTS PASSED\n");
    }
    printf("Tests run: %d\n", tests_run);

    return result != 0;
}

