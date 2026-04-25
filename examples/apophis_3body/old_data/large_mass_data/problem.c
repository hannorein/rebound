/**
 * Sun-Earth-Apophis 3-body run loaded from a saved snapshot.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <inttypes.h>
#include "rebound.h"

static FILE* diag_fp = NULL;
static double e_initial = 0.0;

void heartbeat(struct reb_simulation* r){
    const double e = reb_simulation_energy(r);
    const double rel_e = (e - e_initial) / e_initial;
    fprintf(diag_fp, "%" PRIu64 ",%.16e,%.16e,%.16e\n",
            r->steps_done, r->t, r->dt_last_done, rel_e);
}

int main(int argc, char *argv[])
{
    const char* integrator_name = (argc > 1) ? argv[1] : "ias15";
    const char* output_csv = (argc > 2) ? argv[2] : "diag.csv";

    // struct reb_simulation *r = reb_simulation_create();
    struct reb_simulation *r = reb_simulation_create_from_file("find_close_encounter/apophis_1yr_before.bin", -1);
    if (r == NULL){
        fprintf(stderr, "Error: could not load snapshot file find_close_encounter/apophis_1yr_before.bin\n");
        return 1;
    }

    const double tmax = r->t + 365.25 * 2.0; // Run for 2 years from the saved snapshot.

    // This allows you to connect to the simulation using
    // a web browser by pointing it to http://localhost:1234
    reb_simulation_start_server(r, 1234);

    // Keep the saved IAS15 state unless we explicitly switch integrators.
    if (strcmp(integrator_name, "whfast") == 0){
        const double k = 0.01720209895; // Gaussian constant
        r->dt = 10.0;                   // in days
        r->G = k * k;                   // same units as used by the mercury6 code.
        r->force_is_velocity_dependent = 0; // Force only depends on positions.
        r->integrator = REB_INTEGRATOR_WHFAST;
        r->ri_whfast.safe_mode = 1;     // Keep synchronized for reliable output.
        r->ri_whfast.corrector = 11;    // Turn on symplectic corrector.
        r->ri_whfast.coordinates = REB_WHFAST_COORDINATES_JACOBI;
    }else if (strcmp(integrator_name, "ias15") == 0){
        r->integrator = REB_INTEGRATOR_IAS15;
    }else{
        fprintf(stderr, "Error: unknown integrator '%s'. Use 'ias15' or 'whfast'.\n", integrator_name);
        reb_simulation_free(r);
        return 1;
    }

    // Initial conditions are loaded from snapshot file.
    if (r->N < 3){
        fprintf(stderr, "Error: snapshot must contain at least 3 particles, found N=%u\n", r->N);
        reb_simulation_free(r);
        return 1;
    }
    r->N_active = 3;                // Sun, Earth, and Apophis are active.

    diag_fp = fopen(output_csv, "w");
    if (diag_fp == NULL){
        fprintf(stderr, "Error: could not open output file '%s'.\n", output_csv);
        reb_simulation_free(r);
        return 1;
    }
    fprintf(diag_fp, "step,t,dt,relE\n");
    e_initial = reb_simulation_energy(r);
    r->heartbeat = heartbeat;

    // Start integration
    // reb_simulation_integrate(r, INFINITY); // Runs forever

    reb_simulation_integrate(r, tmax);      // Integrates only to tmax

    double e_final = reb_simulation_energy(r);
    printf("\nDone (%s). Final time: %.4f. Relative energy error: %e\n", integrator_name, r->t, fabs((e_final - e_initial) / e_initial));
    fclose(diag_fp);

    // Cleanup
    reb_simulation_free(r);
}
