/**
 * Sun-Earth-Apophis 3-body run loaded from a saved snapshot.
 * This simulation stop at close encounter
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

double tmax = 7.3e8;
double dmin = 1e100;
double t_dmin = 0.0;
double encounter_threshold = 0.01; // AU
int encounter_reported = 0;
int stop_on_encounter = 1;

void heartbeat(struct reb_simulation *const r);

int main(int argc, char *argv[])
{
    // struct reb_simulation *r = reb_simulation_create();
    struct reb_simulation *r = reb_simulation_create_from_file("find_close_encounter/apophis_1yr_before.bin", -1);
    if (r == NULL){
        fprintf(stderr, "Error: could not load snapshot file find_close_encounter/apophis_1yr_before.bin\n");
        return 1;
    }

    // This allows you to connect to the simulation using
    // a web browser by pointing it to http://localhost:1234
    reb_simulation_start_server(r, 1234);

    // Setup constants
    const double k = 0.01720209895; // Gaussian constant
    r->dt = 1.0;                    // in days
    r->G = k * k;                   // same units as used by the mercury6 code.
    // r->ri_whfast.safe_mode = 1;     // Keep synchronized for reliable visualization/output.
    // r->ri_whfast.corrector = 11;    // Turn on symplectic correctors (11th order).
    // r->ri_whfast.coordinates = REB_WHFAST_COORDINATES_JACOBI;

    // Setup callbacks:
    r->heartbeat = heartbeat;
    r->force_is_velocity_dependent = 0; // Force only depends on positions.
    // r->integrator = REB_INTEGRATOR_WHFAST;
    r->integrator    = REB_INTEGRATOR_IAS15;

    // Initial conditions are loaded from snapshot file.
    if (r->N < 3){
        fprintf(stderr, "Error: snapshot must contain at least 3 particles, found N=%u\n", r->N);
        reb_simulation_free(r);
        return 1;
    }
    r->N_active = 2;        // Sun and Earth are active.
    r->particles[2].m = 0.; // Apophis as test particle.

    double e_initial = reb_simulation_energy(r);

    // Start integration
    reb_simulation_integrate(r, INFINITY); // Runs forever
    // Do not overwrite the seed snapshot in this run-mode.

    // reb_simulation_integrate(r, tmax);      // Integrates only to tmax

    double e_final = reb_simulation_energy(r);
    printf("\nDone. Final time: %.4f. Relative energy error: %e\n", r->t, fabs((e_final - e_initial) / e_initial));

    // Cleanup
    reb_simulation_free(r);
}

void heartbeat(struct reb_simulation *const r)
{
    const struct reb_particle* p = r->particles;
    const double dx = p[2].x - p[1].x; // Apophis - Earth
    const double dy = p[2].y - p[1].y;
    const double dz = p[2].z - p[1].z;
    const double d = sqrt(dx * dx + dy * dy + dz * dz);

    if (d < dmin){
        dmin = d;
        t_dmin = r->t;
        printf("New dmin(Earth-Apophis) = %.10g AU at t = %.10g days\n", dmin, t_dmin);
    }

    if (!encounter_reported && d < encounter_threshold){
        encounter_reported = 1;
        printf("Close encounter detected: d = %.10g AU at t = %.10g days (threshold = %.6g AU)\n",
               d, r->t, encounter_threshold);
        if (stop_on_encounter){
            ((struct reb_simulation*)r)->status = REB_STATUS_USER;
            printf("Stopping integration at first encounter.\n");
        }
    }

    if (reb_simulation_output_check(r, 100.))
    {
        printf("t = %.3f d, d(Earth-Apophis) = %.10g AU, dmin = %.10g AU at t_dmin = %.10g d\n",
               r->t, d, dmin, t_dmin);
    }
}
