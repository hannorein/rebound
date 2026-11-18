/**
 * Solar System
 *
 * This example integrates all planets of the Solar
 * System. The data comes from the NASA HORIZONS system. 
 */
#include "rebound.h"
#include <math.h>

#define MAX_PARTICLES 1024

void heartbeat(struct reb_simulation* r);
int num_particles;
double e_init;
double tmax;
int particle_index = 0;

// NOTE: If visualizations are not enabled, then the orbit color logic will not affect anything.
// To change the strength of the glow, multiply the color by a scalar less than 1.
    // When rendering the particles, they are looped over in a list.
    // Adding orbit colors adds to a list that is looped over at the same time as the particles being rendered.
    // Therefore, adding one particle, and one color will have that newly added particle render its trail with that new color.
    // More Colors Case: particles = [p0, p1, p2]; colors = [c0, c1, c2, c3];
        // p0 renders with c0
        // p1 renders with c1
        // p2 renders with c2
        // c3 is left alone.
    // More Particles Case: particles = [p0, p1, p2, p3]; colors = [c0, c1, c2]
        // p0 renders with c0
        // p1 renders with c1
        // p2 renders with c2
        // p3 renders with default white color
int main(int argc, char* argv[]){
    const struct reb_vec3df yellow = (struct reb_vec3df){.25f, .25f, 0.f};
    const struct reb_vec3df cyan = (struct reb_vec3df){0.f, .125f, .25f};
    const struct reb_vec3df cyan_strong = (struct reb_vec3df){0.f, 1.f, 1.f};
    struct reb_simulation* r = reb_simulation_create();

    // Starting the REBOUND visualization server. This
    // allows you to visualize the simulation by pointing 
    // your web browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);

    // Setup constants
    r->dt                   = 0.1;       // in units of years/(2*pi)
    tmax                    = 1e6*M_PI*2.0;        // in units of years/(2*pi)

    reb_simulation_synchronize(r);
    reb_simulation_set_integrator(r, "whfast");
    struct reb_integrator_whfast_state* whfast = r->integrator.state;
    whfast->safe_mode  = 1;             // Turn off safe mode. Need to call reb_simulation_synchronize() before outputs. 
    whfast->corrector  = 11;            // 11th order symplectic corrector
    r->exact_finish_time    = 1;     // Finish exactly at tmax in reb_simulation_integrate(). Default is already 1.

    // Initial conditions
    reb_simulation_add_fmt(r, "solar system"); // Built in dataset for testing. (8 bodies)
    struct reb_particle sun = r->particles[0]; // The built in dataset has the first particle added being the sun.
    for(int i = 0; i < 8; i++) {
        reb_display_add_orbit_colors(r, cyan);
    }

    // Orbital elements taken from NASA Horizons on November 10, masses were sourced from wikipedia
    struct reb_particle pluto = reb_particle_from_orbit(
        1.0,
        sun,
        6.515e-9,   // Mass sourced from wikipedia
        39.2362261921,
        0.2444,
        0.299289060132,
        1.92350466542,
        2.00568749997,
        13.9423174709
    );

    struct reb_particle makemake = reb_particle_from_orbit(
        1.0,
        sun,
        1.345e-9,   // Mass sourced from wikipedia
        45.60667833331847,
        0.1579631260528772,
        0.506566438697,
        1.38432832145,
        5.18471276384,
        3.01890826923
    );

    struct reb_particle eris = reb_particle_from_orbit(
        1.0,
        sun,
        8.233e-9,   // Mass sourced from wikipedia
        67.89489239788317,
        0.4389806488279149,
        0.767303993189,
        0.628147076228,
        2.63277495655,
        3.38772513067
    );

    struct reb_particle haumea = reb_particle_from_orbit(
        1.0,
        sun,
        1.976e-9,
        43.09677939130120,
        0.1936255936047158,
        0.492330940402,
        2.12555476604,
        4.19805238537,
        3.68127001789
    );

    // The cool new guy (data from wikipedia)
    // Mass was estimated to be on the same order of magnitude of the other bodies (no data on mass at time of writing example)
    struct reb_particle ed138 = reb_particle_from_orbit(
        1.0,
        sun,
        1e-9,
        43.09,
        0.1387,
        0.426383936262,
        5.63060669986,
        2.94437044811,
        2.446386
    );

    // Add all the dwarf planets and their orbit colors (ensuring order so that colors are correctly applied)
    for(int i = 0; i < 4; i++) {
        reb_display_add_orbit_colors(r, yellow);
    }

    reb_simulation_add(r, pluto);
    reb_simulation_add(r, makemake);
    reb_simulation_add(r, haumea);
    reb_simulation_add(r, eris);

    // Add the new body with the stronger cyan color for greater glow
    reb_display_add_orbit_colors(r, cyan_strong);
    reb_simulation_add(r, ed138);

    reb_simulation_move_to_com(r);
    reb_simulation_integrate(r, tmax);
}