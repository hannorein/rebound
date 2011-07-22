#ifndef _MAIN_H
#define _MAIN_H

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

extern double softening;
extern double G;
extern double t;
extern double tmax;
extern double dt;
extern int N;
extern int N_active_first;
extern int N_active_last;
extern double timing_initial;
extern int root_nx;
extern int root_ny;
extern int root_nz;
extern int root_n;
extern double boxsize;
extern double boxsize_x;
extern double boxsize_y;
extern double boxsize_z;
extern double boxsize_max;

void init_box();
void iterate();

#endif
