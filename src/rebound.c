/**
 * @file    rebound.c
 * @brief   Main REBOUND control structures and routine, iteration loop.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section LICENSE
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
#define _NO_CRT_STDIO_INLINE // WIN32 to use _vsprintf_s
#if defined(_WIN32) && defined(_MSC_VER)
#pragma comment(lib, "legacy_stdio_definitions.lib")
#endif
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h> // for offsetof()
#include <sys/types.h>
#include <errno.h>
#include <string.h>
#include <fcntl.h>
#include "rebound.h"
#include "fmemopen.h" // own implementation of fmemopen
#include "integrator.h"
#include "integrator_saba.h"
#include "integrator_whfast.h"
#include "integrator_ias15.h"
#include "integrator_mercurius.h"
#include "integrator_trace.h"
#include "integrator_bs.h"
#include "boundary.h"
#include "gravity.h"
#include "collision.h"
#include "tree.h"
#include "output.h"
#include "tools.h"
#include "particle.h"
#include "input.h"
#include "binarydiff.h"
#include "simulationarchive.h"
#include "server.h"
#ifdef MPI
#include "communication_mpi.h"
#endif
#include "display.h"
#ifdef OPENMP
#include <omp.h>
#endif
#define MAX(a, b) ((a) < (b) ? (b) : (a))       ///< Returns the maximum of a and b
#define STRINGIFY(s) str(s)
#define str(s) #s
#ifdef _WIN32
void usleep(__int64 usec);
#endif // _WIN32
const int reb_max_messages_length = 1024;   // needs to be constant expression for array size
const int reb_N_max_messages = 10;
const char* reb_build_str = __DATE__ " " __TIME__;  // Date and time build string. 
const char* reb_version_str = "4.4.8";         // **VERSIONLINE** This line gets updated automatically. Do not edit manually.
const char* reb_githash_str = STRINGIFY(GITHASH);             // This line gets updated automatically. Do not edit manually.

static int reb_simulation_error_message_waiting(struct reb_simulation* const r);

void reb_simulation_steps(struct reb_simulation* const r, unsigned int N_steps){
    for (unsigned int i=0;i<N_steps;i++){
        reb_simulation_step(r);
    }
}
void reb_simulation_step(struct reb_simulation* const r){
    // Update walltime
    struct reb_timeval time_beginning;
    gettimeofday(&time_beginning,NULL);

    // A 'DKD'-like integrator will do the first 'D' part.
    PROFILING_START()
    if (r->pre_timestep_modifications){
        reb_simulation_synchronize(r);
        r->pre_timestep_modifications(r);
        r->ri_whfast.recalculate_coordinates_this_timestep = 1;
        r->ri_mercurius.recalculate_coordinates_this_timestep = 1;
    }
    reb_integrator_part1(r);
    PROFILING_STOP(PROFILING_CAT_INTEGRATOR)

    // Update and simplify tree. 
    // Prepare particles for distribution to other nodes. 
    // This function also creates the tree if called for the first time.
    if (r->tree_needs_update || r->gravity==REB_GRAVITY_TREE || r->collision==REB_COLLISION_TREE || r->collision==REB_COLLISION_LINETREE){
        // Check for root crossings.
        PROFILING_START()
        reb_boundary_check(r);     
        PROFILING_STOP(PROFILING_CAT_BOUNDARY)

        // Update tree (this will remove particles which left the box)
        PROFILING_START()
        reb_simulation_update_tree(r);          
        PROFILING_STOP(PROFILING_CAT_GRAVITY)
    }

    PROFILING_START()
#ifdef MPI
    // Distribute particles and add newly received particles to tree.
    reb_communication_mpi_distribute_particles(r);
#endif // MPI

    if (r->tree_root!=NULL && r->gravity==REB_GRAVITY_TREE){
        // Update center of mass and quadrupole moments in tree in preparation of force calculation.
        reb_simulation_update_tree_gravity_data(r); 
#ifdef MPI
        // Prepare essential tree (and particles close to the boundary needed for collisions) for distribution to other nodes.
        reb_tree_prepare_essential_tree_for_gravity(r);

        // Transfer essential tree and particles needed for collisions.
        reb_communication_mpi_distribute_essential_tree_for_gravity(r);
#endif // MPI
    }

    // Calculate accelerations. 
    reb_calculate_acceleration(r);
    if (r->N_var){
        reb_calculate_acceleration_var(r);
    }
    // Calculate non-gravity accelerations. 
    if (r->additional_forces) r->additional_forces(r);
    PROFILING_STOP(PROFILING_CAT_GRAVITY)

    // A 'DKD'-like integrator will do the 'KD' part.
    PROFILING_START()
    reb_integrator_part2(r);
    
    if (r->post_timestep_modifications){
        reb_simulation_synchronize(r);
        r->post_timestep_modifications(r);
        r->ri_whfast.recalculate_coordinates_this_timestep = 1;
        r->ri_mercurius.recalculate_coordinates_this_timestep = 1;
    }
    
    if (r->N_var){
        reb_simulation_rescale_var(r);
    }
    PROFILING_STOP(PROFILING_CAT_INTEGRATOR)

    // Do collisions here. We need both the positions and velocities at the same time.
    // Check for root crossings.
    PROFILING_START()
    reb_boundary_check(r);     
    if (r->tree_needs_update){
        // Update tree (this will remove particles which left the box)
        reb_simulation_update_tree(r);          
    }
    PROFILING_STOP(PROFILING_CAT_BOUNDARY)

    // Search for collisions using local and essential tree.
    PROFILING_START()
    reb_collision_search(r);
    PROFILING_STOP(PROFILING_CAT_COLLISION)
    
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

void reb_exit(const char* const msg){
    // This function should also kill all children. 
    // Not implemented as pid is not easy to get to.
    // kill(pid, SIGKILL);
    fprintf(stderr,"\n\033[1mFatal error! Exiting now.\033[0m %s\n",msg);
    exit(EXIT_FAILURE);
}

void reb_message(struct reb_simulation* const r, char type, const char* const msg){
    int save_messages = 0;
    if (r != NULL){
        save_messages = r->save_messages;
    }
    if (!save_messages || strlen(msg)>=reb_max_messages_length){
        if (type=='w'){
            fprintf(stderr,"\n\033[1mWarning!\033[0m %s\n",msg);
        }else if (type=='e'){
            fprintf(stderr,"\n\033[1mError!\033[0m %s\n",msg);
        }
    }else{
        // TODO: Should be protected by MUTEX
        if (r->messages==NULL){
            r->messages = calloc(reb_N_max_messages,sizeof(char*));
        }
        int n = 0;
        for (;n<reb_N_max_messages;n++){
            if (r->messages[n]==NULL){
                break;
            }
        }
        if (n==reb_N_max_messages){
            free(r->messages[0]);
            for (int i=0;i<reb_N_max_messages-1;i++){
                r->messages[i] = r->messages[i+1];
            }
            r->messages[reb_N_max_messages-1] = NULL;
            n= reb_N_max_messages-1;
        }
        r->messages[n] = malloc(sizeof(char*)*reb_max_messages_length);
        r->messages[n][0] = type;
        strcpy(r->messages[n]+1, msg);
    }
}

void reb_simulation_warning(struct reb_simulation* const r, const char* const msg){
    reb_message(r, 'w', msg);
}

void reb_simulation_error(struct reb_simulation* const r, const char* const msg){
    reb_message(r, 'e', msg);
}

void reb_simulation_stop(struct reb_simulation* const r){
    r->status = REB_STATUS_USER;
}

int reb_simulation_get_next_message(struct reb_simulation* const r, char* const buf){
    if (r->messages){
        char* w0 = r->messages[0];
        if (w0){
            for(int i=0;i<reb_N_max_messages-1;i++){
                r->messages[i] = r->messages[i+1];
            }
            r->messages[reb_N_max_messages-1] = NULL;
            strcpy(buf,w0);
            free(w0);
            return 1;
        }
    }
    return 0;
}

static int reb_simulation_error_message_waiting(struct reb_simulation* const r){
    if (r->messages){
        for (int i=0;i<reb_N_max_messages;i++){
            if (r->messages[i]!=NULL){
                if (r->messages[i][0]=='e'){
                    return 1;
                }
            }
        }
    }
    return 0;
}


void reb_simulation_configure_box(struct reb_simulation* const r, const double root_size, const int N_root_x, const int N_root_y, const int N_root_z){
    r->root_size = root_size;
    r->N_root_x = N_root_x;
    r->N_root_y = N_root_y;
    r->N_root_z = N_root_z;
    // Setup box sizes
    r->boxsize.x = r->root_size *(double)r->N_root_x;
    r->boxsize.y = r->root_size *(double)r->N_root_y;
    r->boxsize.z = r->root_size *(double)r->N_root_z;
    r->N_root = r->N_root_x*r->N_root_y*r->N_root_z;
    r->boxsize_max = MAX(r->boxsize.x, MAX(r->boxsize.y, r->boxsize.z));
    if (r->N_root_x <=0 || r->N_root_y <=0 || r->N_root_z <= 0){
        reb_exit("Number of root boxes must be greater or equal to 1 in each direction.");
    }
}
#ifdef MPI
void reb_mpi_init(struct reb_simulation* const r){
    reb_communication_mpi_init(r,0,NULL);
    // Make sure domain can be decomposed into equal number of root boxes per node.
    if ((r->N_root/r->mpi_num)*r->mpi_num != r->N_root){
        if (r->mpi_id==0) fprintf(stderr,"ERROR: Number of root boxes (%d) not a multiple of mpi nodes (%d).\n",r->N_root,r->mpi_num);
        exit(-1);
    }
    printf("MPI-node: %d. Process id: %d.\n",r->mpi_id, getpid());
}

void reb_mpi_finalize(struct reb_simulation* const r){
    r->mpi_id = 0;
    r->mpi_num = 0;
    MPI_Finalize();
}
#endif // MPI

void reb_simulation_free(struct reb_simulation* const r){
    reb_simulation_free_pointers(r);
    free(r);
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
            printf("Waiting for OpenGL visualization to shut down...\n");
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
    reb_tree_delete(r);
    if (r->gravity_cs){
        free(r->gravity_cs  );
    }
    if (r->collisions){
        free(r->collisions  );
    }
    reb_integrator_whfast_reset(r);
    reb_integrator_ias15_reset(r);
    reb_integrator_mercurius_reset(r);
    reb_integrator_trace_reset(r);
    reb_integrator_bs_reset(r);
    if(r->free_particle_ap){
        for(unsigned int i=0; i<r->N; i++){
            r->free_particle_ap(&r->particles[i]);
        }
    }
    if (r->particles){
        free(r->particles   );
    }
    if (r->particle_lookup_table){
        free(r->particle_lookup_table);
    }
    if (r->messages){
        for (int i=0;i<reb_N_max_messages;i++){
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
    for (int s=0; s<r->N_odes; s++){
        r->odes[s]->r = NULL;
    }
    free(r->odes);
}

int reb_simulation_reset_function_pointers(struct reb_simulation* const r){
    int wasnotnull = 0;
    if (r->coefficient_of_restitution ||
        r->collision_resolve ||
        r->additional_forces ||
        r->heartbeat ||
        r->pre_timestep_modifications ||
        r->post_timestep_modifications ||
        r->free_particle_ap ||
        r->extras_cleanup){
      wasnotnull = 1;
    }
    r->coefficient_of_restitution   = NULL;
    r->collision_resolve        = NULL;
    r->additional_forces        = NULL;
    r->heartbeat            = NULL;
    r->pre_timestep_modifications  = NULL;
    r->post_timestep_modifications  = NULL;
    r->free_particle_ap = NULL;
    r->extras_cleanup = NULL;
    return wasnotnull;
}

struct reb_simulation* reb_simulation_create(){
    struct reb_simulation* r = calloc(1,sizeof(struct reb_simulation));
    reb_simulation_init(r);
    return r;
}


void reb_simulation_copy_with_messages(struct reb_simulation* r_copy,  struct reb_simulation* r, enum reb_simulation_binary_error_codes* warnings){
    char* bufp;
    size_t sizep;
    reb_simulation_save_to_stream(r, &bufp,&sizep);
    
    reb_simulation_free_pointers(r_copy);
    memset(r_copy, 0, sizeof(struct reb_simulation));
    reb_simulation_init(r_copy);

    FILE* fin = reb_fmemopen(bufp, sizep, "r");
    reb_input_fields(r_copy, fin, warnings);
    fclose(fin);

    free(bufp);
}

char* reb_simulation_diff_char(struct reb_simulation* r1, struct reb_simulation* r2){
    char* bufp1;
    char* bufp2;
    char* bufp;
    size_t sizep1, sizep2, size;
    reb_simulation_save_to_stream(r1, &bufp1,&sizep1);
    reb_simulation_save_to_stream(r2, &bufp2,&sizep2);

    reb_binary_diff(bufp1, sizep1, bufp2, sizep2, &bufp, &size, 3);
    
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
    reb_simulation_save_to_stream(r1, &bufp1,&sizep1);
    reb_simulation_save_to_stream(r2, &bufp2,&sizep2);

    int ret = reb_binary_diff(bufp1, sizep1, bufp2, sizep2, NULL, NULL, output_option);
    
    free(bufp1);
    free(bufp2);
    return ret;
}

struct reb_simulation* reb_simulation_copy(struct reb_simulation* r){
    struct reb_simulation* r_copy = reb_simulation_create();
    enum reb_simulation_binary_error_codes warnings = REB_SIMULATION_BINARY_WARNING_NONE;
    reb_simulation_copy_with_messages(r_copy,r,&warnings);
    r = reb_input_process_warnings(r, warnings);
    return r_copy;
}

void reb_clear_pre_post_pointers(struct reb_simulation* const r){
    // Temporary fix for REBOUNDx. 
    r->pre_timestep_modifications  = NULL;
    r->post_timestep_modifications  = NULL;
}

void reb_simulation_init(struct reb_simulation* r){
    memset(r, 0, sizeof(struct reb_simulation));
    r->rand_seed = reb_tools_get_rand_seed();
    reb_simulation_reset_function_pointers(r);
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
    r->N_root   = 1;
    r->N_ghost_x  = 0;
    r->N_ghost_y  = 0;
    r->N_ghost_z  = 0;
    r->N        = 0;    
    r->N_allocated   = 0;    
    r->N_active     = -1;   
    r->var_rescale_warning   = 0;   
    r->particle_lookup_table = NULL;
    r->hash_ctr = 0;
    r->N_lookup = 0;
    r->N_allocated_lookup = 0;
    r->testparticle_type = 0;   
    r->testparticle_hidewarnings = 0;
    r->N_var    = 0;    
    r->N_var_config = 0;    
    r->var_config   = NULL;     
    r->exit_min_distance    = 0;    
    r->exit_max_distance    = 0;    
    r->max_radius0    = 0.;   
    r->max_radius1    = 0.;   
    r->status       = REB_STATUS_SUCCESS;
    r->exact_finish_time    = 1;
    r->force_is_velocity_dependent = 0;
    r->gravity_ignore_terms    = 0;
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
    r->ri_ias15.adaptive_mode = 2; // new default since January 2024
    r->ri_ias15.iterations_max_exceeded = 0;    
    
    // ********** SEI
    r->ri_sei.OMEGA     = 1;
    r->ri_sei.OMEGAZ    = -1;
    r->ri_sei.lastdt    = 0;
    
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
    r->tree_needs_update= 0;
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
    printf("Using OpenMP with %d threads per node.\n",omp_get_max_threads());
#endif // OPENMP
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
    if (reb_simulation_error_message_waiting(r)){
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


void reb_run_heartbeat(struct reb_simulation* const r){
    if (r->heartbeat){ r->heartbeat(r); }               // Heartbeat
    if (r->exit_max_distance){
        // Check for escaping particles
        const double max2 = r->exit_max_distance * r->exit_max_distance;
        const struct reb_particle* const particles = r->particles;
        const int N = r->N - r->N_var;
        for (int i=0;i<N;i++){
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
        const int N = r->N - r->N_var;
        for (int i=0;i<N;i++){
            struct reb_particle pi = particles[i];
            for (int j=0;j<i;j++){
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

////////////////////////////////////////////////////
///  Integrate functions and visualization stuff

volatile sig_atomic_t reb_sigint;

void reb_sigint_handler(int signum) {
    // Handles graceful shutdown for interrupts
    if (signum == SIGINT){
        reb_sigint += 1;
    }
}

struct reb_thread_info {
    struct reb_simulation* r;
    double tmax;
};

static void* reb_simulation_integrate_raw(void* args){
    reb_sigint = 0;
    signal(SIGINT, reb_sigint_handler);
    struct reb_thread_info* thread_info = (struct reb_thread_info*)args;
	struct reb_simulation* const r = thread_info->r;
#ifdef MPI
    // Distribute particles
    reb_communication_mpi_distribute_particles(r);
#endif // MPI

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
    reb_run_heartbeat(r);
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
        reb_run_heartbeat(r);
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

size_t reb_simulation_struct_size(){
    // For unit tests to check if python struct has same size
    return sizeof(struct reb_simulation);
}
int reb_check_fp_contract(){
    // Checks if floating point contractions are on. 
    // If so, this will prevent unit tests from passing
    // and bit-wise reproducibility will fail.
    double a = 1.2382309285234567;
    double b = 2.123478623874623234567;
    double c = 6.0284234234234567;

    double r1 = a*b+c;
    double ab = a*b;
    double r2 = ab+c;

    return r1 != r2;
}

// Wrapper to free pointers from python.
void reb_free(void* p){
    free(p);
}

#ifdef _WIN32

void PyInit_librebound() {};

// Source: https://codebrowser.dev/glibc/glibc/stdlib/rand_r.c.html 
int rand_r(unsigned int *seed) {
     unsigned int next = *seed;
     int result;

     next *= 1103515245;
     next += 12345;
     result = (unsigned int) (next / 65536) % 2048;

     next *= 1103515245;
     next += 12345;
     result <<= 10;
     result ^= (unsigned int) (next / 65536) % 1024;

     next *= 1103515245;
     next += 12345;
     result <<= 10;
     result ^= (unsigned int) (next / 65536) % 1024;

     *seed = next;

     return result;
}


// Source: https://stackoverflow.com/a/40160038/115102
#include <stdarg.h>
int vasprintf(char **strp, const char *fmt, va_list ap) {
    // _vscprintf tells you how big the buffer needs to be
    int len = _vscprintf(fmt, ap);
    if (len == -1) {
        return -1;
    }
    size_t size = (size_t)len + 1;
    char *str = malloc(size);
    if (!str) {
        return -1;
    }
    // _vsprintf_s is the "secure" version of vsprintf
    int r = vsprintf_s(str, len + 1, fmt, ap);
    if (r == -1) {
        free(str);
        return -1;
    }
    *strp = str;
    return r;
}
int asprintf(char **strp, const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    int r = vasprintf(strp, fmt, ap);
    va_end(ap);
    return r;
}

// Source: https://stackoverflow.com/a/26085827/115102
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <stdint.h> // portable: uint64_t   MSVC: __int64
int gettimeofday(struct reb_timeval * tp, struct timezone * tzp)
{
    // Note: some broken versions only have 8 trailing zero's, the correct epoch has 9 trailing zero's
    // This magic number is the number of 100 nanosecond intervals since January 1, 1601 (UTC)
    // until 00:00:00 January 1, 1970
    static const uint64_t EPOCH = ((uint64_t) 116444736000000000ULL);

    SYSTEMTIME  system_time;
    FILETIME    file_time;
    uint64_t    time;

    GetSystemTime( &system_time );
    SystemTimeToFileTime( &system_time, &file_time );
    time =  ((uint64_t)file_time.dwLowDateTime )      ;
    time += ((uint64_t)file_time.dwHighDateTime) << 32;

    tp->tv_sec  = (int64_t) ((time - EPOCH) / 10000000L);
    tp->tv_usec = (int64_t) (system_time.wMilliseconds * 1000);
    return 0;
}

// Source: https://stackoverflow.com/a/17283549/115102
void usleep(__int64 usec)
{
    HANDLE timer;
    LARGE_INTEGER ft;

    ft.QuadPart = -(10*usec); // Convert to 100 nanosecond interval, negative value indicates relative time

    timer = CreateWaitableTimer(NULL, TRUE, NULL);
    SetWaitableTimer(timer, &ft, 0, NULL, NULL, 0);
    WaitForSingleObject(timer, INFINITE);
    CloseHandle(timer);
}
#endif // _WIN32
  
#ifdef OPENMP
void reb_omp_set_num_threads(int num_threads){
    omp_set_num_threads(num_threads);
}
#endif // OPENMP

const char* reb_logo[26] = {
"          _                           _  ",
"         | |                         | | ",
" _ __ ___| |__   ___  _   _ _ __   __| | ",
"| '__/ _ \\ '_ \\ / _ \\| | | | '_ \\ / _` | ",
"| | |  __/ |_) | (_) | |_| | | | | (_| | ",
"|_|  \\___|_.__/ \\___/ \\__,_|_| |_|\\__,_| ",
"                                         ",
"              `-:://::.`                 ",
"          `/oshhoo+++oossso+:`           ",
"       `/ssooys++++++ossssssyyo:`        ",
"     `+do++oho+++osssso++++++++sy/`      ",
"    :yoh+++ho++oys+++++++++++++++ss.     ",
"   /y++hooyyooshooo+++++++++++++++oh-    ",
"  -dsssdssdsssdssssssssssooo+++++++oh`   ",
"  ho++ys+oy+++ho++++++++oosssssooo++so   ",
" .d++oy++ys+++oh+++++++++++++++oosssod   ",
" -h+oh+++yo++++oyo+++++++++++++++++oom   ",
" `d+ho+++ys+++++oys++++++++++++++++++d   ",
"  yys++++oy+++++++oys+++++++++++++++s+   ",
"  .m++++++h+++++++++oys++++++++++++oy`   ",
"   -yo++++ss++++++++++oyso++++++++oy.    ",
"    .ss++++ho+++++++++++osys+++++yo`     ",
"      :ss+++ho+++++++++++++osssss-       ",
"        -ossoys++++++++++++osso.         ",
"          `-/oyyyssosssyso+/.            ",
"                ``....`                  ",
};

const unsigned char reb_favicon_png[] = {
  0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a, 0x00, 0x00, 0x00, 0x0d,
  0x49, 0x48, 0x44, 0x52, 0x00, 0x00, 0x00, 0x10, 0x00, 0x00, 0x00, 0x10,
  0x08, 0x06, 0x00, 0x00, 0x00, 0x1f, 0xf3, 0xff, 0x61, 0x00, 0x00, 0x00,
  0x01, 0x73, 0x52, 0x47, 0x42, 0x00, 0xae, 0xce, 0x1c, 0xe9, 0x00, 0x00,
  0x00, 0x44, 0x65, 0x58, 0x49, 0x66, 0x4d, 0x4d, 0x00, 0x2a, 0x00, 0x00,
  0x00, 0x08, 0x00, 0x01, 0x87, 0x69, 0x00, 0x04, 0x00, 0x00, 0x00, 0x01,
  0x00, 0x00, 0x00, 0x1a, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0xa0, 0x01,
  0x00, 0x03, 0x00, 0x00, 0x00, 0x01, 0x00, 0x01, 0x00, 0x00, 0xa0, 0x02,
  0x00, 0x04, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x10, 0xa0, 0x03,
  0x00, 0x04, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x10, 0x00, 0x00,
  0x00, 0x00, 0x34, 0x55, 0x71, 0xf2, 0x00, 0x00, 0x01, 0xaf, 0x49, 0x44,
  0x41, 0x54, 0x38, 0x11, 0x6d, 0xd3, 0x3f, 0x48, 0x55, 0x51, 0x1c, 0xc0,
  0xf1, 0xf7, 0xb4, 0x87, 0x29, 0x1a, 0xa4, 0x90, 0x6d, 0x89, 0xa2, 0x0e,
  0x49, 0x69, 0xbd, 0x10, 0x1a, 0x2a, 0x74, 0x14, 0x57, 0x77, 0x09, 0x5d,
  0x0c, 0x23, 0x52, 0x70, 0x11, 0x6c, 0x70, 0x73, 0x30, 0x78, 0x98, 0x36,
  0x94, 0x20, 0x51, 0xad, 0xd2, 0x5f, 0x88, 0x1c, 0x55, 0x5c, 0x24, 0x0c,
  0x27, 0x51, 0xb4, 0xa9, 0x10, 0x41, 0x6a, 0xe8, 0x9f, 0xf8, 0xfd, 0x5e,
  0xee, 0x79, 0x1e, 0xd4, 0x1f, 0x7c, 0xee, 0xf9, 0x77, 0xcf, 0xb9, 0xf7,
  0xfc, 0xee, 0xb9, 0xd9, 0xcc, 0xc9, 0xb8, 0x42, 0x57, 0x1f, 0x6e, 0xe1,
  0x07, 0xfe, 0xe2, 0x00, 0x1b, 0x78, 0x82, 0x75, 0x9c, 0x1a, 0x65, 0xf4,
  0x3e, 0x46, 0x01, 0x8d, 0xe8, 0xc5, 0x03, 0x18, 0x25, 0x18, 0xc6, 0x3e,
  0x26, 0xe1, 0xbd, 0x49, 0x38, 0x60, 0x94, 0xe3, 0x43, 0x5a, 0xce, 0x51,
  0x6e, 0x63, 0x16, 0xbe, 0xcd, 0x23, 0xbc, 0x83, 0x71, 0x09, 0xf3, 0x78,
  0x81, 0xb3, 0xc8, 0x64, 0xbd, 0x10, 0x53, 0x58, 0xc1, 0x6f, 0x34, 0xa1,
  0x15, 0x7b, 0xd8, 0x41, 0x1b, 0x7c, 0xd0, 0x2f, 0x78, 0xcf, 0x22, 0x9c,
  0xdc, 0x8d, 0xfb, 0x67, 0xb8, 0x5c, 0x83, 0x13, 0x9f, 0x23, 0x8e, 0x09,
  0x1a, 0x79, 0x98, 0x87, 0x87, 0xf8, 0x83, 0x1b, 0xe8, 0x80, 0xf9, 0x69,
  0xc1, 0x33, 0xdf, 0x60, 0x06, 0xde, 0xec, 0x42, 0x9d, 0xf8, 0x87, 0x2d,
  0x54, 0x63, 0x04, 0x97, 0xe1, 0xb8, 0x39, 0x58, 0x83, 0xe1, 0xbc, 0x21,
  0x34, 0xd8, 0xf8, 0x8a, 0xcf, 0x30, 0xf3, 0x35, 0xb8, 0x03, 0x33, 0x5f,
  0x8b, 0x10, 0xe7, 0xa9, 0xbc, 0x85, 0xdb, 0x8b, 0xe3, 0xbd, 0x8d, 0x65,
  0x54, 0x46, 0xbd, 0xd3, 0xd4, 0x07, 0x31, 0x16, 0xf5, 0x59, 0xbd, 0x00,
  0x13, 0x5d, 0x61, 0x23, 0x8d, 0x79, 0x93, 0xe3, 0x1e, 0x7f, 0xa6, 0x1d,
  0x17, 0x29, 0xcf, 0xa1, 0x80, 0xf6, 0xb4, 0x2f, 0x14, 0xdf, 0xa9, 0xb8,
  0x95, 0xd1, 0xd0, 0x41, 0x99, 0x75, 0x81, 0xd2, 0xa8, 0xa3, 0x87, 0xfa,
  0x4b, 0x78, 0x70, 0x5c, 0xb4, 0x0a, 0x71, 0x7c, 0xa2, 0xe1, 0xbe, 0x7d,
  0x90, 0x91, 0x73, 0x81, 0x6f, 0xa8, 0x83, 0x61, 0x76, 0x17, 0x92, 0x5a,
  0x26, 0xf3, 0x85, 0xf2, 0x6a, 0x5a, 0x8f, 0x0b, 0x93, 0x7e, 0x17, 0xf5,
  0xd8, 0x74, 0xc0, 0x6f, 0xee, 0xe9, 0x32, 0x4c, 0x54, 0x08, 0xbf, 0xf3,
  0xbd, 0xd0, 0x88, 0x4a, 0xdf, 0xd8, 0x5c, 0x38, 0xa7, 0xd5, 0x37, 0x58,
  0x45, 0x0e, 0xb7, 0x11, 0x6f, 0xc7, 0x43, 0x93, 0xc7, 0xf1, 0xf8, 0x4f,
  0xc7, 0x3e, 0xc2, 0xdc, 0x64, 0xdc, 0xb3, 0xfd, 0x0a, 0x1f, 0x93, 0xd6,
  0xd1, 0xc5, 0x27, 0xc5, 0x8b, 0x3a, 0xd2, 0x81, 0xd7, 0x48, 0x8e, 0x72,
  0x18, 0x74, 0xd5, 0x37, 0xb8, 0x8e, 0x2e, 0xf8, 0xe7, 0xed, 0xc2, 0x2f,
  0x72, 0x13, 0x4b, 0x68, 0xc6, 0x38, 0x1a, 0x30, 0x00, 0x4f, 0x6f, 0xf1,
  0x5f, 0xb0, 0x1e, 0xc2, 0x1f, 0xa8, 0x1f, 0x75, 0xf0, 0xc4, 0xb9, 0x35,
  0xcf, 0x8a, 0x07, 0xee, 0x29, 0xd6, 0x50, 0x8c, 0x43, 0x2d, 0xa4, 0x52,
  0x79, 0x8a, 0xe5, 0x13, 0x77, 0x00, 0x00, 0x00, 0x00, 0x49, 0x45, 0x4e,
  0x44, 0xae, 0x42, 0x60, 0x82
};
const unsigned int reb_favicon_len = 581;
