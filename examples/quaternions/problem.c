/**
 * Using Quaternions in REBOUND
 * 
 * A simple example showing various use cases of quaternions
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_create_simulation();

    reb_add_fmt(r, "m", 1.);                // Central object
    reb_add_fmt(r, "m a e", 1e-3, 1., 0.1); // Jupiter mass planet


    struct reb_quat q1 = reb_quat_identity();
    struct reb_vec3d v = {.x = 1, .y = 2, .z = 3};

    printf("v = %.5f %.5f %.5f\n", v.x, v.y, v.z);
    v = reb_quat_act(q1,v);
    printf("v = %.5f %.5f %.5f\n", v.x, v.y, v.z);

    struct reb_vec3d axis = {.x = 1, .y = 0, .z = 0};
    q1 = reb_quat_from_angle_axis(M_PI/2.12, axis);
    v = reb_quat_act(q1, v);
    printf("v = %.5f %.5f %.5f\n", v.x, v.y, v.z);
    
    axis.x = 0; axis.y = 1; axis.z = 0;
    struct reb_quat q2 = reb_quat_from_angle_axis(M_PI/2.12, axis);
    v = reb_quat_act(q2, v);
    printf("v = %.5f %.5f %.5f\n", v.x, v.y, v.z);
    
    struct reb_quat q3 = reb_quat_mul(q2, q1);
    struct reb_vec3d v3 = {.x = 1, .y = 2, .z = 3};
    v3 = reb_quat_act(q3, v3);
    printf("v3 = %.5f %.5f %.5f\n", v3.x, v3.y, v3.z);

    reb_free_simulation(r);
}

