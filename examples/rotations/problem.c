/**
 * Using rotations in REBOUND
 * 
 * A simple example showing how to use the built-in rotations framework in REBOUND. See also Rotations.ipynb for more details
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
    // REBOUND provides a reb_rotation struct.
    // Internally, it is implemented using quaternions but
    // you don't need to understand how quaternions work!
    // Simply put: this struct can rotate a vector or an 
    // entire simulation.
    
    // The following example rotates the vector v around the 
    // x axis by 90 degrees:
    struct reb_vec3d axis = {.x = 1, .y = 0, .z = 0};
    struct reb_rotation r1 = reb_rotation_init_angle_axis(M_PI/2.0, axis);
    
    struct reb_vec3d v = {.x = 1, .y = 2, .z = 3};
    struct reb_vec3d v_rotated = reb_vec3d_rotate(v, r1);
    printf("v_rotated = %.5f %.5f %.5f\n", v_rotated.x, v_rotated.y, v_rotated.z);


    // You can rotate a particle (its position and velocity)
    struct reb_particle p = {.m=1, .x=1, .vy=1};
    reb_particle_irotate(&p, r1); // irotate means rotate in place

    // You can also rotate all the particles in a simulation:
    struct reb_simulation* r = reb_simulation_create();
    reb_simulation_add_fmt(r, "m", 1.);                // Central object
    reb_simulation_add_fmt(r, "m a e", 1e-3, 1., 0.1); // Jupiter mass planet
    reb_simulation_irotate(r, r1);
    reb_simulation_free(r);

    // You can chain rotations by multiplying them together.
    // Note that the order of rotations matters, just like
    // the order matters when multiplying together two matricies.
    struct reb_vec3d axis_y = {.y = 1,};
    struct reb_rotation r2 = reb_rotation_init_angle_axis(M_PI/2.0, axis_y);
    struct reb_rotation r_combined = reb_rotation_mul(r2, r1); 
    v_rotated = reb_vec3d_rotate(v, r_combined); // equal to applying r1 first, then r2
   
    // You can easily calculate the inverse of rotations.
    struct reb_rotation r_inverse = reb_rotation_inverse(r_combined);
    v_rotated = reb_vec3d_rotate(v, r_inverse);
   
    // For celestial mechanics, we provide a special init method that
    // uses the ascending node, inclination and longitude of periastron.
    // Applying this method to an orbit in the xy plane with the
    // pericenter on the x axis gives the same result as initializing an 
    // particle with orbital parameters the "normal way". 
    double Omega = 0.12;
    double inc = 0.223;
    double omega = 0.345;
    struct reb_rotation r_orbit = reb_rotation_init_orbit(Omega, inc, omega);
    
    r = reb_simulation_create();
    reb_simulation_add_fmt(r, "m", 1.);          // Central object
    reb_simulation_add_fmt(r, "a e", 1., 0.001); // orbit in the xy plane 
    reb_simulation_add_fmt(r, "a e Omega inc omega ", 1., 0.001, Omega, inc, omega); // 3d orbit
    reb_particle_irotate(&r->particles[1], r_orbit);  
    printf("particle[1] = %.5f %.5f %.5f   %.5f %.5f %.5f\n", r->particles[1].x, r->particles[1].y, r->particles[1].z, r->particles[1].vx, r->particles[1].vy, r->particles[1].vz);
    printf("particle[2] = %.5f %.5f %.5f   %.5f %.5f %.5f\n", r->particles[2].x, r->particles[2].y, r->particles[2].z, r->particles[2].vx, r->particles[2].vy, r->particles[2].vz);
    reb_simulation_free(r);

    
    // REBOUND also comes with a built-in constructor that generates a rotation
    // which rotates a given vector to a new vector. For example:
    struct reb_vec3d v1 = {.x = 1, .y = 0, .z = 0};
    struct reb_vec3d v2 = {.x = 4, .y = 5, .z = 6};

    struct reb_rotation r3 = reb_rotation_init_from_to(v1, v2);
    v_rotated = reb_vec3d_rotate(v1, r3);
    
    v2 = reb_vec3d_normalize(v2); // for easy comparison 
    printf("v2        = %.5f %.5f %.5f\n", v2.x, v2.y, v2.z);
    printf("v_rotated = %.5f %.5f %.5f\n", v_rotated.x, v_rotated.y, v_rotated.z);
}

