/**
 * @file    output.c
 * @brief   Output routines.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section     LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include "particle.h"
#include "rebound.h"
#include "tools.h"
#include "output.h"
#include "integrator.h"
#include "integrator_sei.h"

#include "input.h"
#ifdef MPI
#include "communication_mpi.h"
#include "mpi.h"
#endif // MPI


// List of REBOUND parameters to be written to a file.
// Modify this list if you wish to input/output additional fields in the reb_simulation structure.
const struct reb_binary_field_descriptor reb_binary_field_descriptor_list[]= {
    { 0,  REB_DOUBLE,       "t",                            offsetof(struct reb_simulation, t), 0, 0},
    { 1,  REB_DOUBLE,       "G",                            offsetof(struct reb_simulation, G), 0, 0},
    { 2,  REB_DOUBLE,       "softening",                    offsetof(struct reb_simulation, softening), 0, 0},
    { 3,  REB_DOUBLE,       "dt",                           offsetof(struct reb_simulation, dt), 0, 0},
    { 4,  REB_UINT,         "N",                            offsetof(struct reb_simulation, N), 0, 0},
    { 5,  REB_INT,          "N_var",                        offsetof(struct reb_simulation, N_var), 0, 0},
    // 6 Used to be varconfig
    { 7,  REB_INT,          "N_active",                     offsetof(struct reb_simulation, N_active), 0, 0},
    { 8,  REB_INT,          "testparticle_type",            offsetof(struct reb_simulation, testparticle_type), 0, 0},
    { 9,  REB_INT,          "hash_ctr",                     offsetof(struct reb_simulation, hash_ctr), 0, 0},
    { 10, REB_DOUBLE,       "opening_angle2",               offsetof(struct reb_simulation, opening_angle2), 0, 0},
    { 11, REB_INT,          "status",                       offsetof(struct reb_simulation, status), 0, 0},
    { 12, REB_INT,          "exact_finish_time",            offsetof(struct reb_simulation, exact_finish_time), 0, 0},
    { 13, REB_UINT,         "force_is_velocity_dependent",  offsetof(struct reb_simulation, force_is_velocity_dependent), 0, 0},
    { 14, REB_UINT,         "gravity_ignore_terms",         offsetof(struct reb_simulation, gravity_ignore_terms), 0, 0},
    { 15, REB_DOUBLE,       "output_timing_last",           offsetof(struct reb_simulation, output_timing_last), 0, 0},
    { 16, REB_INT,          "save_messages",                offsetof(struct reb_simulation, save_messages), 0, 0},
    { 17, REB_DOUBLE,       "exit_max_distance",            offsetof(struct reb_simulation, exit_max_distance), 0, 0},
    { 18, REB_DOUBLE,       "exit_min_distance",            offsetof(struct reb_simulation, exit_min_distance), 0, 0},
    { 19, REB_DOUBLE,       "usleep",                       offsetof(struct reb_simulation, usleep), 0, 0},
    { 20, REB_INT,          "track_energ_yoffset",          offsetof(struct reb_simulation, track_energy_offset), 0, 0},
    { 21, REB_DOUBLE,       "energy_offset",                offsetof(struct reb_simulation, energy_offset), 0, 0},
    { 22, REB_VEC3D,        "boxsize",                      offsetof(struct reb_simulation, boxsize), 0, 0},
    { 23, REB_DOUBLE,       "boxsize_max",                  offsetof(struct reb_simulation, boxsize_max), 0, 0},
    { 24, REB_DOUBLE,       "root_size",                    offsetof(struct reb_simulation, root_size), 0, 0},
    { 25, REB_INT,          "root_n",                       offsetof(struct reb_simulation, root_n), 0, 0},
    { 26, REB_INT,          "root_nx",                      offsetof(struct reb_simulation, root_nx), 0, 0},
    { 27, REB_INT,          "root_ny",                      offsetof(struct reb_simulation, root_ny), 0, 0},
    { 28, REB_INT,          "root_nz",                      offsetof(struct reb_simulation, root_nz), 0, 0},
    { 29, REB_INT,          "nghostx",                      offsetof(struct reb_simulation, nghostx), 0, 0},
    { 30, REB_INT,          "nghosty",                      offsetof(struct reb_simulation, nghosty), 0, 0},
    { 31, REB_INT,          "nghostz",                      offsetof(struct reb_simulation, nghostz), 0, 0},
    { 32, REB_INT,          "collision_resolve_keep_sorted",offsetof(struct reb_simulation, collision_resolve_keep_sorted), 0, 0},
    { 33, REB_DOUBLE,       "minimum_collision_velocity",   offsetof(struct reb_simulation, minimum_collision_velocity), 0, 0},
    { 34, REB_DOUBLE,       "collisions_plog",              offsetof(struct reb_simulation, collisions_plog), 0, 0},
    { 36, REB_LONG,         "collisions_Nlog",              offsetof(struct reb_simulation, collisions_Nlog), 0, 0},
    { 37, REB_INT,          "calculate_megno",              offsetof(struct reb_simulation, calculate_megno), 0, 0},
    { 38, REB_DOUBLE,       "megno_Ys",                     offsetof(struct reb_simulation, megno_Ys), 0, 0},
    { 39, REB_DOUBLE,       "megno_Yss",                    offsetof(struct reb_simulation, megno_Yss), 0, 0},
    { 40, REB_DOUBLE,       "megno_cov_Yt",                 offsetof(struct reb_simulation, megno_cov_Yt), 0, 0},
    { 41, REB_DOUBLE,       "megno_var_t",                  offsetof(struct reb_simulation, megno_var_t), 0, 0},
    { 42, REB_DOUBLE,       "megno_mean_t",                 offsetof(struct reb_simulation, megno_mean_t), 0, 0},
    { 43, REB_DOUBLE,       "megno_mean_Y",                 offsetof(struct reb_simulation, megno_mean_Y), 0, 0},
    { 44, REB_LONG,         "megno_n",                      offsetof(struct reb_simulation, megno_n), 0, 0},
    { 45, REB_OTHER,        "simulationarchive_size_first", offsetof(struct reb_simulation, simulationarchive_size_first), 0, 0}, // Manually calculated
    { 46, REB_LONG,         "simulationarchive_size_snapshot", offsetof(struct reb_simulation, simulationarchive_size_snapshot), 0, 0},
    { 47, REB_DOUBLE,       "simulationarchive_auto_interval", offsetof(struct reb_simulation, simulationarchive_auto_interval), 0, 0},
    { 102, REB_DOUBLE,      "simulationarchive_auto_walltime", offsetof(struct reb_simulation, simulationarchive_auto_walltime), 0, 0},
    { 48, REB_DOUBLE,       "simulationarchive_next",       offsetof(struct reb_simulation, simulationarchive_next), 0, 0},
    { 50, REB_INT,          "collision",                    offsetof(struct reb_simulation, collision), 0, 0},
    { 51, REB_INT,          "integrator",                   offsetof(struct reb_simulation, integrator), 0, 0},
    { 52, REB_INT,          "boundary",                     offsetof(struct reb_simulation, boundary), 0, 0},
    { 53, REB_INT,          "gravity",                      offsetof(struct reb_simulation, gravity), 0, 0},
    { 54, REB_DOUBLE,       "ri_sei.OMEGA",                 offsetof(struct reb_simulation, ri_sei.OMEGA), 0, 0},
    { 55, REB_DOUBLE,       "ri_sei.OMEGAZ",                offsetof(struct reb_simulation, ri_sei.OMEGAZ), 0, 0},
    { 56, REB_DOUBLE,       "ri_sei.lastdt",                offsetof(struct reb_simulation, ri_sei.lastdt), 0, 0},
    { 57, REB_DOUBLE,       "ri_sei.sindt",                 offsetof(struct reb_simulation, ri_sei.sindt), 0, 0},
    { 58, REB_DOUBLE,       "ri_sei.tandt",                 offsetof(struct reb_simulation, ri_sei.tandt), 0, 0},
    { 59, REB_DOUBLE,       "ri_sei.sindtz",                offsetof(struct reb_simulation, ri_sei.sindtz), 0, 0},
    { 60, REB_DOUBLE,       "ri_sei.tandtz",                offsetof(struct reb_simulation, ri_sei.tandtz), 0, 0},
    { 61, REB_UINT,         "ri_whfast.corrector",          offsetof(struct reb_simulation, ri_whfast.corrector), 0, 0},
    { 62, REB_UINT,         "ri_whfast.recalculate_coordinates_this_timestep", offsetof(struct reb_simulation, ri_whfast.recalculate_coordinates_this_timestep), 0, 0},
    { 63, REB_UINT,         "ri_whfast.safe_mode",          offsetof(struct reb_simulation, ri_whfast.safe_mode), 0, 0},
    { 64, REB_UINT,         "ri_whfast.keep_unsynchronized",offsetof(struct reb_simulation, ri_whfast.keep_unsynchronized), 0, 0},
    { 65, REB_UINT,         "ri_whfast.is_synchronized",    offsetof(struct reb_simulation, ri_whfast.is_synchronized), 0, 0},
    { 66, REB_UINT,         "ri_whfast.timestep_warnning",  offsetof(struct reb_simulation, ri_whfast.timestep_warning), 0, 0},
    { 69, REB_DOUBLE,       "ri_ias15.epsilon",             offsetof(struct reb_simulation, ri_ias15.epsilon), 0, 0},
    { 70, REB_DOUBLE,       "ri_ias15.min_dt",              offsetof(struct reb_simulation, ri_ias15.min_dt), 0, 0},
    { 71, REB_UINT,         "ri_ias15.adaptive_mode",       offsetof(struct reb_simulation, ri_ias15.adaptive_mode), 0, 0},
    { 72, REB_ULONG,        "ri_ias15.iterations_max_exceeded", offsetof(struct reb_simulation, ri_ias15.iterations_max_exceeded), 0, 0},
    { 85, REB_POINTER,      "particles",                    offsetof(struct reb_simulation, particles), offsetof(struct reb_simulation, N), sizeof(struct reb_particle)},
    { 86, REB_POINTER,      "var_config",                   offsetof(struct reb_simulation, var_config), offsetof(struct reb_simulation, var_config_N), sizeof(struct reb_variational_configuration)},
    { 87, REB_OTHER,        "functionpointers", 0, 0, 0},
    { 89, REB_POINTER,      "ri_ias15.at",                  offsetof(struct reb_simulation, ri_ias15.at), offsetof(struct reb_simulation, ri_ias15.allocated_N), sizeof(double)},
    { 90, REB_POINTER,      "ri_ias15.x0",                  offsetof(struct reb_simulation, ri_ias15.x0), offsetof(struct reb_simulation, ri_ias15.allocated_N), sizeof(double)},
    { 91, REB_POINTER,      "ri_ias15.v0",                  offsetof(struct reb_simulation, ri_ias15.v0), offsetof(struct reb_simulation, ri_ias15.allocated_N), sizeof(double)},
    { 92, REB_POINTER,      "ri_ias15.a0",                  offsetof(struct reb_simulation, ri_ias15.a0), offsetof(struct reb_simulation, ri_ias15.allocated_N), sizeof(double)},
    { 93, REB_POINTER,      "ri_ias15.csx",                 offsetof(struct reb_simulation, ri_ias15.csx), offsetof(struct reb_simulation, ri_ias15.allocated_N), sizeof(double)},
    { 94, REB_POINTER,      "ri_ias15.csv",                 offsetof(struct reb_simulation, ri_ias15.csv), offsetof(struct reb_simulation, ri_ias15.allocated_N), sizeof(double)},
    { 95, REB_POINTER,      "ri_ias15.csa0",                offsetof(struct reb_simulation, ri_ias15.csa0), offsetof(struct reb_simulation, ri_ias15.allocated_N), sizeof(double)},
    { 96, REB_DP7,          "ri_ias15.g",                   offsetof(struct reb_simulation, ri_ias15.g), offsetof(struct reb_simulation, ri_ias15.allocated_N), 7*sizeof(double)},
    { 97, REB_DP7,          "ri_ias15.b",                   offsetof(struct reb_simulation, ri_ias15.b), offsetof(struct reb_simulation, ri_ias15.allocated_N), 7*sizeof(double)},
    { 98, REB_DP7,          "ri_ias15.csb",                 offsetof(struct reb_simulation, ri_ias15.csb), offsetof(struct reb_simulation, ri_ias15.allocated_N), 7*sizeof(double)},
    { 99, REB_DP7,          "ri_ias15.e",                   offsetof(struct reb_simulation, ri_ias15.e), offsetof(struct reb_simulation, ri_ias15.allocated_N), 7*sizeof(double)},
    { 100, REB_DP7,         "ri_ias15.br",                  offsetof(struct reb_simulation, ri_ias15.br), offsetof(struct reb_simulation, ri_ias15.allocated_N), 7*sizeof(double)},
    { 101, REB_DP7,         "ri_ias15.er",                  offsetof(struct reb_simulation, ri_ias15.er), offsetof(struct reb_simulation, ri_ias15.allocated_N), 7*sizeof(double)},
    { 104, REB_POINTER,     "ri_whfast.p_jh",               offsetof(struct reb_simulation, ri_whfast.p_jh), offsetof(struct reb_simulation, ri_whfast.allocated_N), sizeof(struct reb_particle)},
    { 107, REB_INT,         "visualization",                offsetof(struct reb_simulation, visualization), 0, 0},
    { 112, REB_POINTER,     "ri_janus.p_int",               offsetof(struct reb_simulation, ri_janus.p_int), offsetof(struct reb_simulation, ri_janus.allocated_N), sizeof(struct reb_particle_int)},
    { 113, REB_DOUBLE,      "ri_janus.scale_pos",           offsetof(struct reb_simulation, ri_janus.scale_pos), 0, 0},
    { 114, REB_DOUBLE,      "ri_janus.scale_vel",           offsetof(struct reb_simulation, ri_janus.scale_vel), 0, 0},
    { 115, REB_UINT,        "ri_janus.order",               offsetof(struct reb_simulation, ri_janus.order), 0, 0},
    { 116, REB_UINT,        "ri_janus.recalculate_integer_coordinates_this_timestep", offsetof(struct reb_simulation, ri_janus.recalculate_integer_coordinates_this_timestep), 0, 0},
    { 117, REB_INT,         "ri_whfast.coordinates",        offsetof(struct reb_simulation, ri_whfast.coordinates), 0, 0},
    { 118, REB_DOUBLE,      "ri_mercurius.hillfac",         offsetof(struct reb_simulation, ri_mercurius.hillfac), 0, 0},
    { 119, REB_UINT,        "ri_mercurius.safe_mode",       offsetof(struct reb_simulation, ri_mercurius.safe_mode), 0, 0},
    { 120, REB_UINT,        "ri_mercurius.is_synchronized", offsetof(struct reb_simulation, ri_mercurius.is_synchronized), 0, 0},
    { 122, REB_POINTER,     "ri_mercurius.dcrit",           offsetof(struct reb_simulation, ri_mercurius.dcrit), offsetof(struct reb_simulation, ri_mercurius.dcrit_allocated_N), sizeof(double)},
    { 123, REB_UINT,        "ri_mercurius.recalculate_coordinates_this_timestep", offsetof(struct reb_simulation, ri_mercurius.recalculate_coordinates_this_timestep), 0, 0},
    { 125, REB_INT,         "simulationarchive_version",    offsetof(struct reb_simulation, simulationarchive_version), 0, 0},
    { 126, REB_DOUBLE,      "walltime",                     offsetof(struct reb_simulation, walltime), 0, 0},
    { 130, REB_UINT32,      "python_unit_l",                offsetof(struct reb_simulation, python_unit_l), 0, 0},
    { 131, REB_UINT32,      "python_unit_m",                offsetof(struct reb_simulation, python_unit_m), 0, 0},
    { 132, REB_UINT32,      "python_unit_t",                offsetof(struct reb_simulation, python_unit_t), 0, 0},
    { 133, REB_VEC3D,       "ri_mercurius.com_pos",         offsetof(struct reb_simulation, ri_mercurius.com_pos), 0, 0},
    { 134, REB_VEC3D,       "ri_mercurius.com_vel",         offsetof(struct reb_simulation, ri_mercurius.com_vel), 0, 0},
    { 135, REB_ULONGLONG,   "simulationarchive_auto_step",  offsetof(struct reb_simulation, simulationarchive_auto_step), 0, 0},
    { 136, REB_ULONGLONG,   "simulationarchive_next_step",  offsetof(struct reb_simulation, simulationarchive_next_step), 0, 0},
    { 137, REB_ULONGLONG,   "steps_done",                   offsetof(struct reb_simulation, steps_done), 0, 0},
    { 140, REB_UINT,        "ri_saba.safe_mode",            offsetof(struct reb_simulation, ri_saba.safe_mode), 0, 0},
    { 141, REB_UINT,        "ri_saba.is_synchronized",      offsetof(struct reb_simulation, ri_saba.is_synchronized), 0, 0},
    { 143, REB_UINT,        "ri_whfast.corrector2",         offsetof(struct reb_simulation, ri_whfast.corrector2), 0, 0},
    { 144, REB_INT,         "ri_whfast.kernel",             offsetof(struct reb_simulation, ri_whfast.kernel), 0, 0},
    { 145, REB_DOUBLE,      "dt_last_done",                 offsetof(struct reb_simulation, dt_last_done), 0, 0},
    { 146, REB_INT,         "ri_saba.type",                 offsetof(struct reb_simulation, ri_saba.type), 0, 0},
    { 147, REB_UINT,        "ri_saba.keep_unsynchronized",  offsetof(struct reb_simulation, ri_saba.keep_unsynchronized), 0, 0},
    { 148, REB_INT,         "ri_eos.phi0",                  offsetof(struct reb_simulation, ri_eos.phi0), 0, 0},
    { 149, REB_INT,         "ri_eos.phi1",                  offsetof(struct reb_simulation, ri_eos.phi1), 0, 0},
    { 150, REB_UINT,        "ri_eos.n",                     offsetof(struct reb_simulation, ri_eos.n), 0, 0},
    { 151, REB_UINT,        "ri_eos.safe_mode",             offsetof(struct reb_simulation, ri_eos.safe_mode), 0, 0},
    { 152, REB_UINT,        "ri_eos.is_synchronized",       offsetof(struct reb_simulation, ri_eos.is_synchronized), 0, 0},
    { 154, REB_UINT,        "rand_seed",                    offsetof(struct reb_simulation, rand_seed), 0, 0},
    { 155, REB_INT,         "testparticle_hidewarnings",    offsetof(struct reb_simulation, testparticle_hidewarnings), 0, 0},
    { 156, REB_DOUBLE,      "ri_bs.eps_abs",                offsetof(struct reb_simulation, ri_bs.eps_abs), 0, 0},
    { 157, REB_DOUBLE,      "ri_bs.eps_rel",                offsetof(struct reb_simulation, ri_bs.eps_rel), 0, 0},
    { 158, REB_DOUBLE,      "ri_bs.min_dt",                 offsetof(struct reb_simulation, ri_bs.min_dt), 0, 0},
    { 159, REB_DOUBLE,      "ri_bs.max_dt",                 offsetof(struct reb_simulation, ri_bs.max_dt), 0, 0},
    { 160, REB_INT,         "ri_bs.firstOrLastStep",        offsetof(struct reb_simulation, ri_bs.firstOrLastStep), 0, 0},
    { 161, REB_INT,         "ri_bs.previousRejected",       offsetof(struct reb_simulation, ri_bs.previousRejected), 0, 0},
    { 162, REB_INT,         "ri_bs.targetIter",             offsetof(struct reb_simulation, ri_bs.targetIter), 0, 0},
//    { 163, REB_INT,         "var_rescale_warning", offsetof(struct reb_simulation, var_rescale_warning), 0, 0},
    // TES Variables used to have ids 300 - 388. Do not reuse. 
    { 390, REB_UINT,        "ri_whfast512.keep_unsynchronized", offsetof(struct reb_simulation, ri_whfast512.keep_unsynchronized), 0, 0},
    { 391, REB_UINT,        "ri_whfast512.is_synchronized", offsetof(struct reb_simulation, ri_whfast512.is_synchronized), 0, 0},
    { 392, REB_UINT,        "ri_whfast512.gr_potential",    offsetof(struct reb_simulation, ri_whfast512.gr_potential), 0, 0},
    { 394, REB_POINTER_ALIGNED, "ri_whfast512.pjh",         offsetof(struct reb_simulation, ri_whfast512.p_jh), offsetof(struct reb_simulation, ri_whfast512.allocated_N), sizeof(struct reb_particle_avx512)},
    { 396, REB_DOUBLE,      "max_radius0",                  offsetof(struct reb_simulation, max_radius0), 0, 0},
    { 397, REB_DOUBLE,      "max_radius1",                  offsetof(struct reb_simulation, max_radius1), 0, 0},
    { 398, REB_UINT,        "ri_whfast512.systems_N",       offsetof(struct reb_simulation, ri_whfast512.systems_N), 0, 0},
    { 399, REB_PARTICLE4,   "ri_whfast512.pjh0",            offsetof(struct reb_simulation, ri_whfast512.p_jh0), 0, 0},
    { 1329743186, REB_OTHER,"header", 0, 0, 0},
    { 9998, REB_OTHER,      "sablob", 0, 0, 0},
    { 9999, REB_FIELD_END,  "end", 0, 0, 0}
};

// required for python pickling
void reb_output_free_stream(char* buf){
    free(buf);
}

/** 
 * @brief Replacement for open_memstream
 */
void reb_output_stream_write(char** bufp, size_t* allocatedsize, size_t* sizep, void* restrict data, size_t size){
    // Increase size
    int increased = 0;
    while (*allocatedsize==0 || (*sizep)+size>(*allocatedsize)){
        increased = 1;
	    *allocatedsize = (*allocatedsize) ? (*allocatedsize) * 2 : 32;
    }
    if (increased){
        *bufp = realloc(*bufp,*allocatedsize);
    }
    // Copy data to buffer
    memcpy((*bufp)+(*sizep),data,size);
    *sizep += size;
}

/**
 * @brief Same as reb_output_check but with a phase argument
 */
int reb_output_check_phase(struct reb_simulation* r, double interval,double phase){
    double shift = r->t+interval*phase;
    if (floor(shift/interval)!=floor((shift-r->dt)/interval)){
        return 1;
    }
    // Output at beginning 
    if (r->t==0){
        return 1;
    }
    return 0;
}

int reb_output_check(struct reb_simulation* r, double interval){
    return reb_output_check_phase(r, interval,0);
}


#ifdef PROFILING
#warning PROFILING enabled. Rebound is NOT thread-safe.
double profiling_time_sum[PROFILING_CAT_NUM];
double profiling_time_initial   = 0;
double profiling_timing_initial = 0;
double profiling_time_final     = 0;
void profiling_start(void){
    struct timeval tim;
    gettimeofday(&tim, NULL);
    profiling_time_initial = tim.tv_sec+(tim.tv_usec/1000000.0);
}
void profiling_stop(int cat){
    struct timeval tim;
    gettimeofday(&tim, NULL);
    profiling_time_final = tim.tv_sec+(tim.tv_usec/1000000.0);
    profiling_time_sum[cat] += profiling_time_final - profiling_time_initial;
}
#endif // PROFILING

void reb_output_timing(struct reb_simulation* r, const double tmax){
    const int N = r->N;
#ifdef MPI
    int N_tot = 0;
    MPI_Reduce(&N, &N_tot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
    if (r->mpi_id!=0) return;
#else
    int N_tot = N;
#endif
    struct timeval tim;
    gettimeofday(&tim, NULL);
    double temp = tim.tv_sec+(tim.tv_usec/1000000.0);
    if (r->output_timing_last==-1){
        r->output_timing_last = temp;
    }else{
        printf("\r");
#ifdef PROFILING
        fputs("\033[A\033[2K",stdout);
        for (int i=0;i<=PROFILING_CAT_NUM;i++){
            fputs("\033[A\033[2K",stdout);
        }
#endif // PROFILING
    }
    printf("N_tot= %- 9d  ",N_tot);
    if (r->integrator==REB_INTEGRATOR_SEI){
        printf("t= %- 9f [orb]  ",r->t*r->ri_sei.OMEGA/2./M_PI);
    }else{
        printf("t= %- 9f  ",r->t);
    }
    printf("dt= %- 9f  ",r->dt);
    printf("cpu= %- 9f [s]  ",temp-r->output_timing_last);
    if (tmax>0){
        printf("t/tmax= %5.2f%%",r->t/tmax*100.0);
    }
#ifdef PROFILING
    if (profiling_timing_initial==0){
        struct timeval tim;
        gettimeofday(&tim, NULL);
        profiling_timing_initial = tim.tv_sec+(tim.tv_usec/1000000.0);
    }
    printf("\nCATEGORY       TIME \n");
    double _sum = 0;
    for (int i=0;i<=PROFILING_CAT_NUM;i++){
        switch (i){
            case PROFILING_CAT_INTEGRATOR:
                printf("Integrator     ");
                break;
            case PROFILING_CAT_BOUNDARY:
                printf("Boundary check ");
                break;
            case PROFILING_CAT_GRAVITY:
                printf("Gravity/Forces ");
                break;
            case PROFILING_CAT_COLLISION:
                printf("Collisions     ");
                break;
#ifdef OPENGL
            case PROFILING_CAT_VISUALIZATION:
                printf("Visualization  ");
                break;
#endif // OPENGL
            case PROFILING_CAT_NUM:
                printf("Other          ");
                break;
        }
        if (i==PROFILING_CAT_NUM){
            printf("%5.2f%%",(1.-_sum/(profiling_time_final - profiling_timing_initial))*100.);
        }else{
            printf("%5.2f%%\n",profiling_time_sum[i]/(profiling_time_final - profiling_timing_initial)*100.);
            _sum += profiling_time_sum[i];
        }
    }
#endif // PROFILING
    fflush(stdout);
    r->output_timing_last = temp;
}


void reb_output_ascii(struct reb_simulation* r, char* filename){
    const int N = r->N;
#ifdef MPI
    char filename_mpi[1024];
    sprintf(filename_mpi,"%s_%d",filename,r->mpi_id);
    FILE* of = fopen(filename_mpi,"ab"); 
#else // MPI
    FILE* of = fopen(filename,"ab"); 
#endif // MPI
    if (of==NULL){
        reb_error(r, "Can not open file.");
        return;
    }
    for (int i=0;i<N;i++){
        struct reb_particle p = r->particles[i];
        fprintf(of,"%e\t%e\t%e\t%e\t%e\t%e\n",p.x,p.y,p.z,p.vx,p.vy,p.vz);
    }
    fclose(of);
}

void reb_output_orbits(struct reb_simulation* r, char* filename){
    const int N = r->N;
#ifdef MPI
    char filename_mpi[1024];
    sprintf(filename_mpi,"%s_%d",filename,r->mpi_id);
    FILE* of = fopen(filename_mpi,"ab"); 
#else // MPI
    FILE* of = fopen(filename,"ab"); 
#endif // MPI
    if (of==NULL){
        reb_error(r, "Can not open file.");
        return;
    }
    struct reb_particle com = r->particles[0];
    for (int i=1;i<N;i++){
        struct reb_orbit o = reb_tools_particle_to_orbit(r->G, r->particles[i],com);
        fprintf(of,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",r->t,o.a,o.e,o.inc,o.Omega,o.omega,o.l,o.P,o.f);
        com = reb_get_com_of_pair(com,r->particles[i]);
    }
    fclose(of);
}

// Macro to write a single field to a binary file.
// Memset forces padding to be set to 0 (not necessary but
// helps when comparing binary files)
#define WRITE_FIELD_TYPE(typen, value, length) {\
        struct reb_binary_field field;\
        memset(&field,0,sizeof(struct reb_binary_field));\
        field.type = typen;\
        field.size = (length);\
        reb_output_stream_write(bufp, &allocatedsize, sizep, &field,sizeof(struct reb_binary_field));\
        reb_output_stream_write(bufp, &allocatedsize, sizep, value,field.size);\
    }


void reb_output_binary_to_stream(struct reb_simulation* r, char** bufp, size_t* sizep){
    size_t allocatedsize = 0;
    *bufp = NULL;
    *sizep = 0;
    // Init integrators. This helps with bit-by-bit reproducibility.
    reb_integrator_init(r);

    // Output header.
    char header[64] = "\0";
    int cwritten = sprintf(header,"REBOUND Binary File. Version: %s",reb_version_str);
    snprintf(header+cwritten+1,64-cwritten-1,"%s",reb_githash_str);
    reb_output_stream_write(bufp, &allocatedsize, sizep, header,sizeof(char)*64);

    // Compress data if possible
    // This does not affect future calculation, but might trigger a realloc.
    if (r->ri_ias15.allocated_N > 3*r->N){
        r->ri_ias15.allocated_N = 3*r->N;
    }
    /// Output all fields
    int i=0;
    while (reb_binary_field_descriptor_list[i].dtype!=REB_FIELD_END){
        int dtype = reb_binary_field_descriptor_list[i].dtype;
        // Simple data types:
        if (dtype == REB_DOUBLE || dtype == REB_INT || dtype == REB_UINT || dtype == REB_UINT32 ||
                dtype == REB_LONG || dtype == REB_ULONG || dtype == REB_ULONGLONG ||
                dtype == REB_PARTICLE || dtype == REB_PARTICLE4 || dtype == REB_VEC3D ){
            struct reb_binary_field field;
            memset(&field,0,sizeof(struct reb_binary_field));
            field.type = reb_binary_field_descriptor_list[i].type;
            switch (dtype){
                case REB_DOUBLE: 
                    field.size = sizeof(double);
                    break;
                case REB_INT: 
                    field.size = sizeof(int);
                    break;
                case REB_UINT: 
                    field.size = sizeof(unsigned int);
                    break;
                case REB_UINT32: 
                    field.size = sizeof(uint32_t);
                    break;
                case REB_LONG:
                    field.size = sizeof(long);
                    break;
                case REB_ULONG:
                    field.size = sizeof(unsigned long);
                    break;
                case REB_ULONGLONG:
                    field.size = sizeof(unsigned long long);
                    break;
                case REB_VEC3D:
                    field.size = sizeof(struct reb_vec3d);
                    break;
                case REB_PARTICLE:
                    field.size = sizeof(struct reb_particle);
                    break;
                case REB_PARTICLE4:
                    field.size = 4*sizeof(struct reb_particle);
                    break;
            }
            reb_output_stream_write(bufp, &allocatedsize, sizep, &field, sizeof(struct reb_binary_field));
            char* pointer = (char*)r + reb_binary_field_descriptor_list[i].offset;
            reb_output_stream_write(bufp, &allocatedsize, sizep, pointer, field.size);
        }
        // Pointer data types
        if (dtype == REB_POINTER || dtype == REB_POINTER_ALIGNED ){
            struct reb_binary_field field;
            memset(&field,0,sizeof(struct reb_binary_field));
            field.type = reb_binary_field_descriptor_list[i].type;
            unsigned int* pointer_N = (unsigned int*)((char*)r + reb_binary_field_descriptor_list[i].offset_N);
            field.size = (*pointer_N) * reb_binary_field_descriptor_list[i].element_size;
                
            if (field.size){
                reb_output_stream_write(bufp, &allocatedsize, sizep, &field, sizeof(struct reb_binary_field));
                char* pointer = (char*)r + reb_binary_field_descriptor_list[i].offset;
                pointer = *(char**)pointer;
                reb_output_stream_write(bufp, &allocatedsize, sizep, pointer, field.size);
            }
        }
        // Special datatype for IAS15. Similar to POINTER
        if (dtype == REB_DP7 ){
            struct reb_binary_field field;
            memset(&field,0,sizeof(struct reb_binary_field));
            field.type = reb_binary_field_descriptor_list[i].type;
            unsigned int* pointer_N = (unsigned int*)((char*)r + reb_binary_field_descriptor_list[i].offset_N);
            field.size = (*pointer_N) * reb_binary_field_descriptor_list[i].element_size;
                
            if (field.size){
                reb_output_stream_write(bufp, &allocatedsize, sizep, &field, sizeof(struct reb_binary_field));
                char* pointer = (char*)r + reb_binary_field_descriptor_list[i].offset;
                struct reb_dp7* dp7 = (struct reb_dp7*)pointer;
                reb_output_stream_write(bufp, &allocatedsize, sizep, dp7->p0,field.size/7);
                reb_output_stream_write(bufp, &allocatedsize, sizep, dp7->p1,field.size/7);
                reb_output_stream_write(bufp, &allocatedsize, sizep, dp7->p2,field.size/7);
                reb_output_stream_write(bufp, &allocatedsize, sizep, dp7->p3,field.size/7);
                reb_output_stream_write(bufp, &allocatedsize, sizep, dp7->p4,field.size/7);
                reb_output_stream_write(bufp, &allocatedsize, sizep, dp7->p5,field.size/7);
                reb_output_stream_write(bufp, &allocatedsize, sizep, dp7->p6,field.size/7);
            }
        }
        i++;
    }

    // Write function pointer warning flag
    int functionpointersused = 0;
    if (r->coefficient_of_restitution ||
        r->collision_resolve ||
        r->additional_forces ||
        r->heartbeat ||
        r->post_timestep_modifications ||
        r->free_particle_ap){
        functionpointersused = 1;
    }

    struct reb_binary_field field_functionp;
    memset(&field_functionp,0,sizeof(struct reb_binary_field));
    field_functionp.type = 87; // TODO do not hardcode. 
    field_functionp.size = sizeof(int);
    reb_output_stream_write(bufp, &allocatedsize, sizep, &field_functionp, sizeof(struct reb_binary_field));
    reb_output_stream_write(bufp, &allocatedsize, sizep, &functionpointersused, field_functionp.size);

    // To output size of binary file, need to calculate it first. 
    if (r->simulationarchive_version<3){ // to be removed in a future release
        r->simulationarchive_size_first = (*sizep)+sizeof(struct reb_binary_field)*2+sizeof(long)+sizeof(struct reb_simulationarchive_blob16);
    }else{
        r->simulationarchive_size_first = (*sizep)+sizeof(struct reb_binary_field)*2+sizeof(long)+sizeof(struct reb_simulationarchive_blob);
    }
    WRITE_FIELD_TYPE( 45 ,        &r->simulationarchive_size_first,   sizeof(long));
    int end_null = 0;
    
    struct reb_binary_field_descriptor fd_end = reb_binary_field_descriptor_for_name("end");
    WRITE_FIELD_TYPE(fd_end.type, &end_null, 0);
    if (r->simulationarchive_version<3){ // to be removed in a future release
        struct reb_simulationarchive_blob16 blob = {0};
        reb_output_stream_write(bufp, &allocatedsize, sizep, &blob, sizeof(struct reb_simulationarchive_blob16));
    }else{
        struct reb_simulationarchive_blob blob = {0};
        reb_output_stream_write(bufp, &allocatedsize, sizep, &blob, sizeof(struct reb_simulationarchive_blob));
    }
}

void reb_output_binary(struct reb_simulation* r, const char* filename){
#ifdef MPI
    char filename_mpi[1024];
    sprintf(filename_mpi,"%s_%d",filename,r->mpi_id);
    FILE* of = fopen(filename_mpi,"wb"); 
#else // MPI
    FILE* of = fopen(filename,"wb"); 
#endif // MPI
    if (of==NULL){
        reb_error(r, "Can not open file.");
        return;
    }
    char* bufp;
    size_t sizep;
    reb_output_binary_to_stream(r, &bufp,&sizep);
    fwrite(bufp,sizep,1,of);
    free(bufp);
    fclose(of);
}

void reb_output_binary_positions(struct reb_simulation* r, const char* filename){
    const int N = r->N;
#ifdef MPI
    char filename_mpi[1024];
    sprintf(filename_mpi,"%s_%d",filename,r->mpi_id);
    FILE* of = fopen(filename_mpi,"wb"); 
#else // MPI
    FILE* of = fopen(filename,"wb"); 
#endif // MPI
    if (of==NULL){
        reb_error(r, "Can not open file.");
        return;
    }
    for (int i=0;i<N;i++){
        struct reb_vec3d v;
        v.x = r->particles[i].x;
        v.y = r->particles[i].y;
        v.z = r->particles[i].z;
        fwrite(&(v),sizeof(struct reb_vec3d),1,of);
    }
    fclose(of);
}

void reb_output_velocity_dispersion(struct reb_simulation* r, char* filename){
    const int N = r->N;
    // Algorithm with reduced roundoff errors (see wikipedia)
    struct reb_vec3d A = {.x=0, .y=0, .z=0};
    struct reb_vec3d Q = {.x=0, .y=0, .z=0};
    for (int i=0;i<N;i++){
        struct reb_vec3d Aim1 = A;
        struct reb_particle p = r->particles[i];
        A.x = A.x + (p.vx-A.x)/(double)(i+1);
        if (r->integrator==REB_INTEGRATOR_SEI){
            A.y = A.y + (p.vy+1.5*r->ri_sei.OMEGA*p.x-A.y)/(double)(i+1);
        }else{
            A.y = A.y + (p.vy-A.y)/(double)(i+1);
        }
        A.z = A.z + (p.vz-A.z)/(double)(i+1);
        Q.x = Q.x + (p.vx-Aim1.x)*(p.vx-A.x);
        if (r->integrator==REB_INTEGRATOR_SEI){
            Q.y = Q.y + (p.vy+1.5*r->ri_sei.OMEGA*p.x-Aim1.y)*(p.vy+1.5*r->ri_sei.OMEGA*p.x-A.y);
        }else{
            Q.y = Q.y + (p.vy-Aim1.y)*(p.vy-A.y);
        }
        Q.z = Q.z + (p.vz-Aim1.z)*(p.vz-A.z);
    }
#ifdef MPI
    int N_tot = 0;
    struct reb_vec3d A_tot = {.x=0, .y=0, .z=0};
    struct reb_vec3d Q_tot = {.x=0, .y=0, .z=0};
    MPI_Reduce(&N, &N_tot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
    MPI_Reduce(&A, &A_tot, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
    MPI_Reduce(&Q, &Q_tot, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
    if (r->mpi_id!=0) return;
#else
    int N_tot = N;
    struct reb_vec3d A_tot = A;
    struct reb_vec3d Q_tot = Q;
#endif
    Q_tot.x = sqrt(Q_tot.x/(double)N_tot);
    Q_tot.y = sqrt(Q_tot.y/(double)N_tot);
    Q_tot.z = sqrt(Q_tot.z/(double)N_tot);
    FILE* of = fopen(filename,"ab"); 
    if (of==NULL){
        reb_error(r, "Can not open file.");
        return;
    }
    fprintf(of,"%e\t%e\t%e\t%e\t%e\t%e\t%e\n",r->t,A_tot.x,A_tot.y,A_tot.z,Q_tot.x,Q_tot.y,Q_tot.z);
    fclose(of);
}

    
