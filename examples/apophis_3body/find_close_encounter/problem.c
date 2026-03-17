/**
 * Sun Earth apophis 3 body System
 * This example use Jacobi coordinate with WHFAST integrator
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

// Approximate Apophis mass in solar-mass units.
// static const double apophis_mass = 1.358e-20;
static const double apophis_mass = 3.0034896149156e-6;

double ss_pos[3][3] =
    {
        {0.0, 0.0, 0.0}, // Sun
        {1.0, 0.0, 0.0}, // Earth
        {0.9224, 0.0, 0.0}, // Apophis
    };
double ss_vel[3][3] =
    {
        {0.0, 0.0, 0.0}, // Sun
        {0.0, +0.01720209895, 0.0}, // Earth
        {0.0, +0.0185, 0.0}, // Apophis 
    };

double ss_mass[3] =
    {
        1.0, // Sun
        3.0034896149156e-6, // Earth
        apophis_mass, // Apophis
    };

double tmax = 7.3e8;
double dmin = 1e100;
double t_dmin = 0.0;
double encounter_threshold = 0.01; // AU
int encounter_reported = 0;
int stop_on_encounter = 1;

void heartbeat(struct reb_simulation *const r);

int main(int argc, char *argv[])
{
    struct reb_simulation *r = reb_simulation_create();

    // Setup constants
    const double k = 0.01720209895; // Gaussian constant
    r->dt = 1.0;                    // in days
    r->G = k * k;                   // same units as used by the mercury6 code.
    
    // Setup callbacks:
    r->heartbeat = heartbeat;
    r->force_is_velocity_dependent = 0; // Force only depends on positions.
    r->integrator    = REB_INTEGRATOR_IAS15;

    // Initial conditions
    for (int i = 0; i < 3; i++)
    {
        struct reb_particle p = {0};
        p.x = ss_pos[i][0];
        p.y = ss_pos[i][1];
        p.z = ss_pos[i][2];
        p.vx = ss_vel[i][0];
        p.vy = ss_vel[i][1];
        p.vz = ss_vel[i][2];
        p.m = ss_mass[i];
        reb_simulation_add(r, p);
    }

    reb_simulation_move_to_com(r);

    r->N_active = 3; // Sun, Earth, and Apophis all contribute to gravity.

    double e_initial = reb_simulation_energy(r);

    // Start integration

    // reb_simulation_integrate(r, 23009.9522-365.25); // Runs to 1 year before close encounter

    reb_simulation_integrate(r, 54207.6304 - 365.25); // big astoid

    // Save a clean restart state without live server/runtime pointers attached.
    reb_simulation_save_to_file(r, "apophis_1yr_before.bin"); // save it

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
