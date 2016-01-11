//
//  functions.h
//  
//
//  Created by Ari Silburt on 2015-06-12.
//
//

#ifndef ____functions__
#define ____functions__

#include <stdio.h>
#include "../../src/rebound.h"

//FUNCTIONS******************************
void legend(char* planetdir, char* legenddir, char* xyz_print, char* CEprint, struct reb_simulation* r, double tmax, double m_planetesimal, double total_planetesimal_mass, int N_planetesimals, double inner, double outer, double powerlaw, double Ms, double drh, double epsilon, int seed, int HYBRID_ON);

void output_to_mercury_swifter(struct reb_simulation* r, double HSR, double tmax, int n_output, int movie_output_interval);

double calc_dt(struct reb_simulation* r, double mp, double Ms, double a, double dRHill, double dt_prev);

void calc_Hill2(struct reb_simulation* r);

double calc_Etot(struct reb_simulation* a, double soft, double dKE);

//void calc_ELtot(double* Etot, double* Ktot, double* Utot, double* Ltot, double planetesimal_mass, struct reb_simulation* r);

void calc_ae(double* a, double* e, double* d, struct reb_simulation* r, int i, double t);

void planetesimal_forces_global(struct reb_simulation *a);

void planetesimal_forces_mini(struct reb_simulation *a);

void check_for_encounter(struct reb_simulation* r, struct reb_simulation* s, int* N_encounters, int N_encounters_previous, double* minimum_r, double* maximum_val, char* removeddir, int* output_it, double* dKE, double soft, double ejection_distance2);

void ini_mini(struct reb_simulation* const r, struct reb_simulation* s, double ias_epsilon, int turn_planetesimal_forces_on, double ias_timestep, double timestep);

void update_global(struct reb_simulation* const s, struct reb_simulation* const r, int N_encounters_previous);

void add_or_subtract_particles(struct reb_simulation* r, struct reb_simulation* s, int N_encounters,int N_encounters_previous, char* CEprint, double soft, double dE_collison, double E0);

void update_previous_global_positions(struct reb_simulation* r, int N_encounters);

void update_encounter_indices(int* N_encounters, int* N_encounters_previous);

void output_frame_per_body(struct reb_particle* particles, char* dir, int N, double t);

void output_frame_per_time(struct reb_particle* particles, char* name, int N, double t, int* movie_counter);

time_t clock_start();

void clock_finish(clock_t t_ini, int N_encounters, int N, char* legenddir);

void global_free();


//EXTERNAL VARIABLES******************************
extern double planetesimal_mass;
extern int* encounter_index;
extern int* previous_encounter_index;
extern double* Hill2;
extern double* x_prev; extern double* y_prev; extern double* z_prev;
extern double t_prev; 
extern int N_encounters_tot; extern int N_encounters_previous;
extern int N_tot;
extern struct reb_simulation* r;
extern double soft;


#endif /* defined(____functions__) */
