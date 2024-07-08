#include <stdio.h>
#include <stdlib.h>
#include <ejecta_transform.h>

typedef struct {
    double ID;
    double x, y, z;
    double vx, vy, vz;
    double mass, density;
} ReadParticle;

void transform(ReadParticle *p, double x_translation) {
    // Rotate frame 180 degrees around y-axis
    double new_x = -p->x;
    double new_z = -p->z;
    double new_vx = -p->vx;
    double new_vz = -p->vz;

    // Translate frame along -x direction
    new_x += x_translation;

    // Update particle's position and velocity
    p->x = new_x;
    p->z = new_z;
    p->vx = new_vx;
    p->vz = new_vz;
}

