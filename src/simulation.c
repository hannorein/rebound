/**
 * @file 	simulation.c
 * @brief 	Functions which create, free, and operate on simulations.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2026 Hanno Rein
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
#include "rebound.h"
#include "rebound_internal.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifdef MPI
#include <mpi.h>
#include "communication_mpi.h"
#endif // MPI
#ifdef OPENMP
#include <omp.h>
#endif // OPENMP
#include "tree.h"
#include "particle.h"
#include "simulation.h"
#include "simulationarchive.h"
#include "binarydata.h"
#include "output.h"
#include "server.h"
#include "display.h"
#include "boundary.h"
#include "binarydata.h"
#include "fmemopen.h" // own implementation of fmemopen
#include "integrator.h"
#include "collision.h"
#include "tools.h"
#include "integrator_saba.h"
#include "integrator_whfast.h"
#include "integrator_ias15.h"
#include "integrator_mercurius.h"
#include "integrator_trace.h"
#include "integrator_bs.h"

#define MAX(a, b) ((a) < (b) ? (b) : (a))       ///< Returns the maximum of a and b

struct reb_thread_info {
    struct reb_simulation* r;
    double tmax;
};

// For unit tests to check if python struct has same size as c struct.
size_t reb_simulation_struct_size(){
    return sizeof(struct reb_simulation);
}

// User triggered stop.
void reb_simulation_stop(struct reb_simulation* const r){
    r->status = REB_STATUS_USER;
}

// Print or save warning message. r can be NULL.
void reb_simulation_warning(struct reb_simulation* const r, const char* const msg){
    int save_messages = 0;
    if (r != NULL && r->save_messages) save_messages = 1;
    reb_message(&r->messages, save_messages, REB_MESSAGE_TYPE_WARNING, msg);
}

// Print or save error message. r can be NULL.
void reb_simulation_error(struct reb_simulation* const r, const char* const msg){
    int save_messages = 0;
    if (r != NULL && r->save_messages) save_messages = 1;
    reb_message(&r->messages, save_messages, REB_MESSAGE_TYPE_ERROR, msg);
}

// Print or save info message. r can be NULL.
void reb_simulation_info(struct reb_simulation* const r, const char* const msg){
    int save_messages = 0;
    if (r != NULL && r->save_messages) save_messages = 1;
    reb_message(&r->messages, save_messages, REB_MESSAGE_TYPE_INFO, msg);
}

// Workaround for setting a python callback function. 
void reb_simulation_set_collision_resolve(struct reb_simulation* r, enum REB_COLLISION_RESOLVE_OUTCOME (*resolve) (struct reb_simulation* const r, struct reb_collision c)){
    r->collision_resolve = resolve;
}

// One bare timestep. Without heartbeat.
static void reb_simulation_step(struct reb_simulation* const r);

struct reb_simulation* reb_simulation_create(){
    struct reb_simulation* r = calloc(1,sizeof(struct reb_simulation));
    reb_simulation_init(r);
    return r;
}

void reb_simulation_free(struct reb_simulation* const r){
    reb_simulation_free_pointers(r);
    free(r);
}

// Temporary fix for REBOUNDx. 
void reb_clear_pre_post_pointers(struct reb_simulation* const r){
    r->pre_timestep_modifications  = NULL;
    r->post_timestep_modifications  = NULL;
}

// Heartbeat wrapper. Runs actual heartbeat and does exit checks.
static void run_heartbeat(struct reb_simulation* const r){
    if (r->heartbeat){ r->heartbeat(r); }               // Heartbeat
    if (r->exit_max_distance){
        // Check for escaping particles
        const double max2 = r->exit_max_distance * r->exit_max_distance;
        const struct reb_particle* const particles = r->particles;
        const size_t N = r->N;
        for (size_t i=0;i<N;i++){
            struct reb_particle p = particles[i];
            double r2 = p.x*p.x + p.y*p.y + p.z*p.z;
            if (r2>max2){
                r->status = REB_STATUS_ESCAPE;
            }
        }
    }
    if (r->exit_min_distance){
        // Check for close encounters
        const double min2 = r->exit_min_distance * r->exit_min_distance;
        const struct reb_particle* const particles = r->particles;
        const size_t N = r->N;
        for (size_t i=0;i<N;i++){
            struct reb_particle pi = particles[i];
            for (size_t j=0;j<i;j++){
                struct reb_particle pj = particles[j];
                const double x = pi.x-pj.x;
                const double y = pi.y-pj.y;
                const double z = pi.z-pj.z;
                const double r2 = x*x + y*y + z*z;
                if (r2<min2){
                    r->status = REB_STATUS_ENCOUNTER;
                }
            }
        }
    }
}

static int error_message_waiting(struct reb_simulation* const r){
    if (r->messages){
        for (size_t i=0;i<reb_messages_max_N;i++){
            if (r->messages[i]!=NULL){
                if (r->messages[i][0]=='e'){
                    return 1;
                }
            }
        }
    }
    return 0;
}


int reb_check_exit(struct reb_simulation* const r, const double tmax, double* last_full_dt){
    if(r->status <= REB_STATUS_SINGLE_STEP){
        if(r->status == REB_STATUS_SINGLE_STEP){
            r->status = REB_STATUS_PAUSED;
        }else{
            // This allows an arbitrary number of steps before the simulation is paused
            r->status++;
        }
    }
    while(r->status == REB_STATUS_PAUSED || r->status == REB_STATUS_SCREENSHOT){
        // Wait for user to disable paused simulation
#ifdef __EMSCRIPTEN__
        emscripten_sleep(100);
#else
        usleep(1000);
#endif 
        if (reb_sigint){ // cancel while paused
            r->status = REB_STATUS_SIGINT;
        }
    }
    const double dtsign = copysign(1.,r->dt);   // Used to determine integration direction
    if (error_message_waiting(r)){
        r->status = REB_STATUS_GENERIC_ERROR;
    }
    if (r->status>=0){
        // Exit now.
    }else if(tmax!=INFINITY){
        if(r->exact_finish_time==1){
            if ((r->t+r->dt)*dtsign>=tmax*dtsign){  // Next step would overshoot
                if (r->t==tmax){
                    r->status = REB_STATUS_SUCCESS;
                }else if(r->status == REB_STATUS_LAST_STEP){
                    double tscale = 1e-12*fabs(tmax);   // Find order of magnitude for time
                    if (tscale<1e-200){     // Failsafe if tmax==0.
                        tscale = 1e-12;
                    }
                    if (fabs(r->t-tmax)<tscale){
                        r->status = REB_STATUS_SUCCESS;
                    }else{
                        // not there yet, do another step.
                        reb_simulation_synchronize(r);
                        r->dt = tmax-r->t;
                    }
                }else{
                    r->status = REB_STATUS_LAST_STEP; // Do one small step, then exit.
                    reb_simulation_synchronize(r);
                    if (r->dt_last_done!=0.){   // If first timestep is also last, do not use dt_last_done (which would be 0.)
                        *last_full_dt = r->dt_last_done; // store last full dt before decreasing the timestep to match finish time
                    }
                    r->dt = tmax-r->t;
                }
            }else{
                if (r->status == REB_STATUS_LAST_STEP){
                    // This will get executed if an adaptive integrator reduces
                    // the timestep in what was supposed to be the last timestep.
                    r->status = REB_STATUS_RUNNING;
                }
            }
        }else{
            if (r->t*dtsign>=tmax*dtsign){  // Past the integration time
                r->status = REB_STATUS_SUCCESS; // Exit now.
            }
        }
    }
#ifndef MPI
    if (!r->N){
        if (!r->N_odes){
            reb_simulation_warning(r,"No particles found. Will exit.");
            r->status = REB_STATUS_NO_PARTICLES; // Exit now.
        }else{
            if (r->integrator != REB_INTEGRATOR_BS){
                reb_simulation_warning(r,"No particles found. Will exit. Use BS integrator to integrate user-defined ODEs without any particles present.");
                r->status = REB_STATUS_NO_PARTICLES; // Exit now.
            }
        }
    }
#else
    int status_max = 0;
    MPI_Allreduce(&(r->status), &status_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 
    if (status_max>=0){
        r->status = status_max;
    }

#endif // MPI
    return r->status;
}



void reb_simulation_free_pointers(struct reb_simulation* const r){
    if (r->simulationarchive_filename){
        free(r->simulationarchive_filename);
    }
    if(r->display_settings){
        free(r->display_settings);
    }
#ifdef OPENGL
    if(r->display_data){
        // Waiting for visualization to shut down.
        if (r->display_data->window){ // Not needed under normal circumstances
            usleep(100);
        }
        if (r->display_data->window){ // still running?
            reb_simulation_info(NULL, "Waiting for OpenGL visualization to shut down...");
            while(r->display_data->window){
                usleep(100);
            }
        }
        pthread_mutex_destroy(&(r->display_data->mutex));
        if (r->display_data->r_copy){
            reb_simulation_free(r->display_data->r_copy);
            r->display_data->r_copy = NULL;
        }
        if (r->display_data->particle_data){
            free(r->display_data->particle_data);
            r->display_data->particle_data = NULL;
        }
        if (r->display_data->orbit_data){
            free(r->display_data->orbit_data);
            r->display_data->orbit_data = NULL;
        }
        free(r->display_data);
        r->display_data = NULL;
    }
#endif //OPENGL
#ifdef SERVER
    reb_simulation_stop_server(r);
#endif // SERVER
    free(r->gravity_cs);
    free(r->collisions);
    reb_integrator_whfast_reset(r);
    reb_integrator_ias15_reset(r);
    reb_integrator_mercurius_reset(r);
    reb_integrator_trace_reset(r);
    reb_integrator_bs_reset(r);
    free(r->tree_root);
    r->tree_root = NULL;
    if (r->ri_custom.reset){
        r->ri_custom.reset(r);
    }
    if(r->free_particle_ap){
        for(size_t i=0; i<r->N; i++){
            r->free_particle_ap(&r->particles[i]);
        }
    }
    free(r->particles);
    free(r->particles_var);
    for (size_t i=0;i<r->N_name_list;i++){
        free(r->name_list[i]);
    }
    free(r->name_list);
    free(r->name_hash_table);
    if (r->messages){
        for (size_t i=0;i<reb_messages_max_N;i++){
            free(r->messages[i]);
        }
    }
    if (r->messages){
        free(r->messages);
    }
    if (r->extras_cleanup){
        r->extras_cleanup(r);
    }
    if (r->var_config){
        free(r->var_config);
    }
    for (size_t s=0; s<r->N_odes; s++){
        r->odes[s]->r = NULL;
    }
    free(r->odes);
}

static void* reb_simulation_integrate_raw(void* args){
    reb_sigint = 0;
    signal(SIGINT, reb_sigint_handler);
    struct reb_thread_info* thread_info = (struct reb_thread_info*)args;
    struct reb_simulation* const r = thread_info->r;

    if (thread_info->tmax != r->t){
        int dt_sign = (thread_info->tmax > r->t) ? 1.0 : -1.0; // determine integration direction
        r->dt = copysign(r->dt, dt_sign);
    }

    double last_full_dt = r->dt; // need to store r->dt in case timestep gets artificially shrunk to meet exact_finish_time=1
    r->dt_last_done = 0.; // Reset in case first timestep attempt will fail

    if (r->testparticle_hidewarnings==0 && reb_particle_check_testparticles(r)){
        reb_simulation_warning(r,"At least one test particle (type 0) has finite mass. This might lead to unexpected behaviour. Set testparticle_hidewarnings=1 to hide this warning.");
    }
    if (r->status != REB_STATUS_PAUSED && r->status != REB_STATUS_SCREENSHOT){ // Allow simulation to be paused initially
        r->status = REB_STATUS_RUNNING;
    }
    run_heartbeat(r);
#ifdef __EMSCRIPTEN__
    double t0 = emscripten_performance_now();
#endif
    while(reb_check_exit(r,thread_info->tmax,&last_full_dt)<0){
#ifdef __EMSCRIPTEN__
        double t1 = emscripten_performance_now();
        if (t1-t0>1000./120.){ // max framerate 120Hz
            t0 = t1;
            emscripten_sleep(0); // allow drawing and event handling
        }

#endif 
#ifdef OPENGL
        if (r->display_data){
            // Note: Mutex is not FIFO.
            // Allow time for mutex to lock in display.c before it is relocked here.
            while (r->display_data->need_copy == 1){
                usleep(10);
            }
            pthread_mutex_lock(&(r->display_data->mutex)); 
        }
#endif //OPENGL
#ifdef SERVER
        if (r->server_data){
            // Note: Mutex is not FIFO.
            // Allow time for mutex to lock in display.c before it is relocked here.
            while (r->server_data->need_copy == 1){
                usleep(10);
            }
#ifdef _WIN32
            WaitForSingleObject(r->server_data->mutex, INFINITE);
#else // _WIN32
            pthread_mutex_lock(&(r->server_data->mutex)); 
#endif // _WIN32
            r->server_data->mutex_locked_by_integrate = 1;
        }
#endif //SERVER
        if (r->simulationarchive_filename){ reb_simulationarchive_heartbeat(r);}
        reb_simulation_step(r); 
        run_heartbeat(r);
        if (reb_sigint){
            r->status = REB_STATUS_SIGINT;
        }
#ifdef OPENGL
        if (r->display_data){
            pthread_mutex_unlock(&(r->display_data->mutex));
        }
#endif //OPENGL
#ifdef SERVER
        if (r->server_data){
#ifdef _WIN32
            ReleaseMutex(r->server_data->mutex);
#else // _WIN32
            pthread_mutex_unlock(&(r->server_data->mutex));
#endif // _WIN32
            r->server_data->mutex_locked_by_integrate = 0;
        }
#endif //SERVER
        if (r->usleep > 0){
            usleep(r->usleep);
        }
    }
    reb_simulation_synchronize(r);
    if(r->exact_finish_time==1){ // if finish_time = 1, r->dt could have been shrunk, so set to the last full timestep
        r->dt = last_full_dt; 
    }
    if (r->simulationarchive_filename){ reb_simulationarchive_heartbeat(r);}

    return NULL;
}


enum REB_STATUS reb_simulation_integrate(struct reb_simulation* const r, double tmax){
    struct reb_thread_info thread_info = {
        .r = r,
        .tmax = tmax, 
    };

#ifdef OPENGL
#ifdef __EMSCRIPTEN__
    if (r->display_data==NULL){
        r->display_data = calloc(sizeof(struct reb_display_data),1);
        r->display_data->r = r;
        // Setup windows, compile shaders, etc.
        reb_display_init(r); // Will return. Display routines running in animation_loop.
    }
    reb_simulation_integrate_raw(&thread_info);
#else // __EMSCRIPTEN__
    if (r->display_data==NULL){
        r->display_data = calloc(sizeof(struct reb_display_data),1);
        r->display_data->r = r;
        if (pthread_mutex_init(&(r->display_data->mutex), NULL)){
            reb_simulation_error(r,"Mutex creation failed.");
        }
    }

    if (pthread_create(&r->display_data->compute_thread,NULL,reb_simulation_integrate_raw,&thread_info)){
        reb_simulation_error(r, "Error creating compute thread.");
    }

    reb_display_init(r); // Display routines need to run on main thread. Will not return until r->status<0.
    if (pthread_join(r->display_data->compute_thread,NULL)){
        reb_simulation_error(r, "Error joining compute thread.");
    }
#endif // __EMSCRIPTEN__
#else // OPENGL
    reb_simulation_integrate_raw(&thread_info);
#endif // OPENGL
    return r->status;
}

static void reset_function_pointers(struct reb_simulation* const r){
    r->coefficient_of_restitution   = NULL;
    r->collision_resolve        = NULL;
    r->additional_forces        = NULL;
    r->heartbeat            = NULL;
    r->pre_timestep_modifications  = NULL;
    r->post_timestep_modifications  = NULL;
    r->free_particle_ap = NULL;
    r->ri_custom.step = NULL;
    r->ri_custom.synchronize = NULL;
    r->ri_custom.reset = NULL;
    r->extras_cleanup = NULL;
}

void reb_simulation_steps(struct reb_simulation* const r, unsigned int N_steps){
    run_heartbeat(r);
    for (unsigned int i=0;i<N_steps;i++){
        reb_simulation_step(r);
        run_heartbeat(r);
    }
}
static void reb_simulation_step(struct reb_simulation* const r){
    // Update walltime
    struct reb_timeval time_beginning;
    gettimeofday(&time_beginning,NULL);

    // A 'DKD'-like integrator will do the first 'D' part.
    if (r->pre_timestep_modifications){
        reb_simulation_synchronize(r);
        r->pre_timestep_modifications(r);
        r->ri_whfast.recalculate_coordinates_this_timestep = 1;
        r->ri_mercurius.recalculate_coordinates_this_timestep = 1;
    }

    PROFILING_START();
    reb_integrator_step(r);
    PROFILING_STOP(PROFILING_CAT_INTEGRATOR);

    if (r->post_timestep_modifications){
        reb_simulation_synchronize(r);
        r->post_timestep_modifications(r);
        r->ri_whfast.recalculate_coordinates_this_timestep = 1;
        r->ri_mercurius.recalculate_coordinates_this_timestep = 1;
    }

    if (r->N_var){
        reb_simulation_rescale_var(r);
    }

    PROFILING_START();
    reb_boundary_check(r);     
    PROFILING_STOP(PROFILING_CAT_BOUNDARY);

    if (r->collision != REB_COLLISION_NONE){
#ifdef MPI
        reb_communication_mpi_distribute_particles(r);
#endif // MPI

        PROFILING_START();
        reb_collision_search(r);
        PROFILING_STOP(PROFILING_CAT_COLLISION);
    }

    // Update walltime
    struct reb_timeval time_end;
    gettimeofday(&time_end,NULL);
    r->walltime_last_step = time_end.tv_sec-time_beginning.tv_sec+(time_end.tv_usec-time_beginning.tv_usec)/1e6;
    r->walltime_last_steps_sum += r->walltime_last_step;
    r->walltime_last_steps_N +=1;
    if (r->walltime_last_steps_sum > 0.1){
        r->walltime_last_steps = r->walltime_last_steps_sum/r->walltime_last_steps_N;
        r->walltime_last_steps_sum =0;
        r->walltime_last_steps_N = 0;
    }
    r->walltime += r->walltime_last_step;
    // Update step counter
    r->steps_done++; // This also counts failed IAS15 steps
}

void reb_simulation_copy_with_messages(struct reb_simulation* r_copy,  struct reb_simulation* r, enum REB_BINARYDATA_ERROR_CODE* warnings){
    char* bufp;
    size_t sizep;
    reb_binarydata_simulation_to_stream(r, &bufp,&sizep);

    reb_simulation_free_pointers(r_copy);
    memset(r_copy, 0, sizeof(struct reb_simulation));
    reb_simulation_init(r_copy);

    FILE* fin = reb_fmemopen(bufp, sizep, "r");
    reb_binarydata_input_fields(r_copy, fin, warnings);
    fclose(fin);

    free(bufp);
}

char* reb_simulation_diff_char(struct reb_simulation* r1, struct reb_simulation* r2){
    char* bufp1;
    char* bufp2;
    char* bufp;
    size_t sizep1, sizep2, size;
    reb_binarydata_simulation_to_stream(r1, &bufp1,&sizep1);
    reb_binarydata_simulation_to_stream(r2, &bufp2,&sizep2);

    reb_binarydata_diff(bufp1, sizep1, bufp2, sizep2, &bufp, &size, 3);

    free(bufp1);
    free(bufp2);
    return bufp;
}

int reb_simulation_diff(struct reb_simulation* r1, struct reb_simulation* r2, int output_option){
    if (output_option!=1 && output_option!=2){
        // Not implemented
        return -1;
    }
    char* bufp1;
    char* bufp2;
    size_t sizep1, sizep2;
    reb_binarydata_simulation_to_stream(r1, &bufp1,&sizep1);
    reb_binarydata_simulation_to_stream(r2, &bufp2,&sizep2);

    int ret = reb_binarydata_diff(bufp1, sizep1, bufp2, sizep2, NULL, NULL, output_option);

    free(bufp1);
    free(bufp2);
    return ret;
}

struct reb_simulation* reb_simulation_copy(struct reb_simulation* r){
    struct reb_simulation* r_copy = reb_simulation_create();
    enum REB_BINARYDATA_ERROR_CODE warnings = REB_BINARYDATA_WARNING_NONE;
    reb_simulation_copy_with_messages(r_copy,r,&warnings);
    r = reb_binarydata_process_warnings(r, warnings);
    return r_copy;
}

void reb_simulation_init(struct reb_simulation* r){
    memset(r, 0, sizeof(struct reb_simulation));
    r->rand_seed = reb_tools_get_rand_seed();
    reset_function_pointers(r);
    r->t        = 0; 
    r->G        = 1;
    r->softening    = 0;
    r->dt       = 0.001;
    r->dt_last_done = 0.;
    r->steps_done = 0;
    r->root_size    = -1;
    r->N_root_x  = 1;
    r->N_root_y  = 1;
    r->N_root_z  = 1;
    r->N_ghost_x  = 0;
    r->N_ghost_y  = 0;
    r->N_ghost_z  = 0;
    r->N        = 0;    
    r->N_allocated   = 0;    
    r->N_active     = SIZE_MAX; 
    r->var_rescale_warning   = 0;   
    r->name_list = NULL;
    r->name_hash_table = NULL;
    r->N_name_list = 0;
    r->testparticle_type = 0;   
    r->testparticle_hidewarnings = 0;
    r->N_var    = 0;    
    r->N_var_config = 0;    
    r->var_config   = NULL;     
    r->exit_min_distance    = 0;    
    r->exit_max_distance    = 0;    
    r->status       = REB_STATUS_SUCCESS;
    r->exact_finish_time    = 1;
    r->force_is_velocity_dependent = 0;
    r->gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_NONE;
    r->calculate_megno  = 0;
    r->output_timing_last   = -1;
    r->save_messages = 0;
    r->track_energy_offset = 0;
    r->server_data = NULL;
    r->display_data = NULL;
    r->display_settings = NULL;
    r->walltime = 0;

    r->minimum_collision_velocity = 0;
    r->collisions_plog  = 0;
    r->collisions_log_n  = 0;    
    r->collisions_N  = 0;    
    r->collision_resolve_keep_sorted   = 0;    

    r->simulationarchive_version       = 3;
    r->simulationarchive_auto_interval = 0.;    
    r->simulationarchive_auto_walltime = 0.;    
    r->simulationarchive_auto_step     = 0;    
    r->simulationarchive_next          = 0.;    
    r->simulationarchive_next_step     = 0;    
    r->simulationarchive_filename      = NULL;    

    // Default modules
    r->integrator   = REB_INTEGRATOR_IAS15;
    r->boundary     = REB_BOUNDARY_NONE;
    r->gravity      = REB_GRAVITY_BASIC;
    r->collision    = REB_COLLISION_NONE;


    // Integrators  
    // ********** WHFAST
    // the defaults below are chosen to safeguard the user against spurious results, but
    // will be slower and less accurate
    r->ri_whfast.corrector = 0;
    r->ri_whfast.corrector2 = 0;
    r->ri_whfast.kernel = 0;
    r->ri_whfast.coordinates = REB_WHFAST_COORDINATES_JACOBI;
    r->ri_whfast.safe_mode = 1;
    r->ri_whfast.recalculate_coordinates_this_timestep = 0;
    r->ri_whfast.is_synchronized = 1;
    r->ri_whfast.timestep_warning = 0;
    r->ri_whfast.recalculate_coordinates_but_not_synchronized_warning = 0;

    // ********** WHFAST512
    r->ri_whfast512.is_synchronized = 1;
    r->ri_whfast512.gr_potential = 0;
    r->ri_whfast512.keep_unsynchronized = 0;
    r->ri_whfast512.recalculate_constants = 1;
    r->ri_whfast512.N_systems = 1;

    // ********** SABA
    r->ri_saba.type = REB_SABA_10_6_4;
    r->ri_saba.safe_mode = 1;
    r->ri_saba.is_synchronized = 1;

    // ********** IAS15
    r->ri_ias15.epsilon         = 1e-9;
    r->ri_ias15.min_dt      = 0;
    r->ri_ias15.adaptive_mode = REB_IAS15_PRS23; // new default since January 2024
    r->ri_ias15.iterations_max_exceeded = 0;    

    // ********** SEI
    r->ri_sei.OMEGA     = 1;
    r->ri_sei.OMEGAZ    = -1;
    r->ri_sei.lastdt    = 0;

    // ********** CUSTOM
    r->ri_custom.data = NULL;
    r->ri_custom.data_size = 0;

    // ********** LEAPFROG
    r->ri_leapfrog.order = 2;

    // ********** MERCURIUS
    r->ri_mercurius.mode = 0;
    r->ri_mercurius.safe_mode = 1;
    r->ri_mercurius.recalculate_coordinates_this_timestep = 0;
    r->ri_mercurius.recalculate_r_crit_this_timestep = 0;
    r->ri_mercurius.is_synchronized = 1;
    r->ri_mercurius.encounter_N = 0;
    r->ri_mercurius.r_crit_hill = 3;

    // ********** JANUS
    r->ri_janus.recalculate_integer_coordinates_this_timestep = 0;
    r->ri_janus.order = 6;
    r->ri_janus.scale_pos = 1e-16;
    r->ri_janus.scale_vel = 1e-16;

    // ********** TRACE
    r->ri_trace.mode = REB_TRACE_MODE_NONE;
    r->ri_trace.peri_mode = REB_TRACE_PERI_FULL_BS;
    r->ri_trace.encounter_N = 0;
    r->ri_trace.r_crit_hill = 3.;
    r->ri_trace.peri_crit_eta = 1.0;
    r->ri_trace.force_accept = 0;

    // ********** EOS
    r->ri_eos.n = 2;
    r->ri_eos.phi0 = REB_EOS_LF;
    r->ri_eos.phi1 = REB_EOS_LF;
    r->ri_eos.safe_mode = 1;
    r->ri_eos.is_synchronized = 1;


    // ********** BS
    reb_integrator_bs_reset(r);

    // Tree parameters. Will not be used unless gravity or collision search makes use of tree.
    r->tree_root        = NULL;
    r->opening_angle2   = 0.25;

#ifdef MPI
    r->mpi_id = 0;                            
    r->mpi_num = 0;                           
    r->particles_send = NULL;  
    r->N_particles_send = 0;                  
    r->N_particles_send_max = 0;               
    r->particles_recv = NULL;     
    r->N_particles_recv = 0;                  
    r->N_particles_recv_max = 0;               

    r->tree_essential_send = NULL;
    r->N_tree_essential_send = 0;             
    r->N_tree_essential_send_max = 0;          
    r->tree_essential_recv = NULL;
    r->N_tree_essential_recv = 0;             
    r->N_tree_essential_recv_max = 0;          
#endif // MPI
#ifdef OPENMP
    char msg[1024];
    sprintf(msg, "Using OpenMP with %d threads per node.\n", omp_get_max_threads());
    reb_simulation_info(r, msg);
#endif // OPENMP
}

// Finds the two largest particles in the simulation. *p1 and *p2 will be set to the indices of the particles.
void reb_simulation_two_largest_particles(struct reb_simulation* r, size_t* p1, size_t* p2) {
    struct reb_particle* particles = r->particles;
    *p1 = SIZE_MAX;
    *p2 = SIZE_MAX;
    double largest1 = -1.0;
    double largest2 = -1.0;
#ifdef OPENMP
    int num_threads;
    // A struct to hold the two largest values found by each thread
    struct two_max {
        double largest1;
        double largest2;
        size_t p1;
        size_t p2;
    };

    // Array to store the two largest values from each thread
    struct two_max *thread_max;
#pragma omp parallel
    {
        num_threads = omp_get_num_threads();
#pragma omp master
        {
            thread_max = (struct two_max *)malloc(num_threads * sizeof(struct two_max));
        }

#pragma omp barrier
        int thread_id = omp_get_thread_num();
        thread_max[thread_id].largest1 = -1.0;
        thread_max[thread_id].largest2 = -1.0;
        thread_max[thread_id].p1 = SIZE_MAX;
        thread_max[thread_id].p2 = SIZE_MAX;

#pragma omp for
        for (size_t i=0; i<r->N; i++) {
            if (particles[i].r > thread_max[thread_id].largest1) {
                thread_max[thread_id].largest2 = thread_max[thread_id].largest1;
                thread_max[thread_id].p2 = thread_max[thread_id].p1;
                thread_max[thread_id].largest1 = particles[i].r;
                thread_max[thread_id].p1 = i;
            } else if (particles[i].r > thread_max[thread_id].largest2) {
                thread_max[thread_id].largest2 = particles[i].r;
                thread_max[thread_id].p2 = i;
            }
        }
    }

    // Reduce the results from all threads
    for (int i=0; i<num_threads; i++) {
        if (thread_max[i].largest1 > largest1) {
            largest2 = largest1;
            *p2 = *p1;
            largest1 = thread_max[i].largest1;
            *p1 = thread_max[i].p1;
        } else if (thread_max[i].largest1 > largest2) {
            largest2 = thread_max[i].largest1;
            *p2 = thread_max[i].p1;
        }

        if (thread_max[i].largest2 > largest2) {
            largest2 = thread_max[i].largest2;
            *p2 = thread_max[i].p2;
        }
    }

    free(thread_max);
#else // OPENMP
    for (size_t i=0; i<r->N; i++) {
        if (particles[i].r > largest1) {
            largest2 = largest1;
            *p2 = *p1;
            largest1 = particles[i].r;
            *p1 = i;
        }else{
            if (particles[i].r > largest2) {
                largest2 = particles[i].r;
                *p2 = i;
            }
        }
    }
#endif // OPENMP
}

void reb_simulation_get_serialized_particle_data(struct reb_simulation* r, double* m, double* radius, double (*xyz)[3], double (*vxvyvz)[3], double (*xyzvxvyvz)[6]){
    const size_t N = r->N;
    struct reb_particle* restrict const particles = r->particles;
    for (size_t i=0;i<N;i++){
        if (m){
            m[i] = particles[i].m;
        }
        if (radius){
            radius[i] = particles[i].r;
        }
        if (xyz){
            xyz[i][0] = particles[i].x;
            xyz[i][1] = particles[i].y;
            xyz[i][2] = particles[i].z;
        }
        if (vxvyvz){
            vxvyvz[i][0] = particles[i].vx;
            vxvyvz[i][1] = particles[i].vy;
            vxvyvz[i][2] = particles[i].vz;
        }
        if (xyzvxvyvz){
            xyzvxvyvz[i][0] = particles[i].x;
            xyzvxvyvz[i][1] = particles[i].y;
            xyzvxvyvz[i][2] = particles[i].z;
            xyzvxvyvz[i][3] = particles[i].vx;
            xyzvxvyvz[i][4] = particles[i].vy;
            xyzvxvyvz[i][5] = particles[i].vz;
        }
    }
}

void reb_simulation_set_serialized_particle_data(struct reb_simulation* r, double* m, double* radius, double (*xyz)[3], double (*vxvyvz)[3], double (*xyzvxvyvz)[6]){
    const size_t N = r->N;
    struct reb_particle* restrict const particles = r->particles;
    for (size_t i=0;i<N;i++){
        if (m){
            particles[i].m = m[i];
        }
        if (radius){
            particles[i].r = radius[i] ;
        }
        if (xyz){
            particles[i].x = xyz[i][0];
            particles[i].y = xyz[i][1];
            particles[i].z = xyz[i][2];
        }
        if (vxvyvz){
            particles[i].vx = vxvyvz[i][0];
            particles[i].vy = vxvyvz[i][1];
            particles[i].vz = vxvyvz[i][2];
        }
        if (xyzvxvyvz){
            particles[i].x = xyzvxvyvz[i][0];
            particles[i].y = xyzvxvyvz[i][1];
            particles[i].z = xyzvxvyvz[i][2];
            particles[i].vx = xyzvxvyvz[i][3];
            particles[i].vy = xyzvxvyvz[i][4];
            particles[i].vz = xyzvxvyvz[i][5];
        }
    }
}
