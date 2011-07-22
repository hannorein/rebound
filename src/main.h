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

void iterate();

#endif
