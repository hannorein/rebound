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
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <string.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/wait.h>
#include <pthread.h>
#include <fcntl.h>
#include "rebound.h"
#include "integrator.h"
#include "integrator_saba.h"
#include "integrator_whfast.h"
#include "integrator_ias15.h"
#include "integrator_mercurius.h"
#include "integrator_bs.h"
#include "integrator_tes.h"
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

const int reb_max_messages_length = 1024;   // needs to be constant expression for array size
const int reb_max_messages_N = 10;
const char* reb_build_str = __DATE__ " " __TIME__;  // Date and time build string. 
const char* reb_version_str = "3.26.3";         // **VERSIONLINE** This line gets updated automatically. Do not edit manually.
const char* reb_githash_str = STRINGIFY(GITHASH);             // This line gets updated automatically. Do not edit manually.

static int reb_error_message_waiting(struct reb_simulation* const r);

void reb_steps(struct reb_simulation* const r, unsigned int N_steps){
    for (unsigned int i=0;i<N_steps;i++){
        reb_step(r);
    }
}
void reb_step(struct reb_simulation* const r){
    // Update walltime
    struct timeval time_beginning;
    gettimeofday(&time_beginning,NULL);

    // A 'DKD'-like integrator will do the first 'D' part.
    PROFILING_START()
    if (r->pre_timestep_modifications){
        reb_integrator_synchronize(r);
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
        reb_tree_update(r);          
        PROFILING_STOP(PROFILING_CAT_GRAVITY)
    }

    PROFILING_START()
#ifdef MPI
    // Distribute particles and add newly received particles to tree.
    reb_communication_mpi_distribute_particles(r);
#endif // MPI

    if (r->tree_root!=NULL && r->gravity==REB_GRAVITY_TREE){
        // Update center of mass and quadrupole moments in tree in preparation of force calculation.
        reb_tree_update_gravity_data(r); 
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
        reb_integrator_synchronize(r);
        r->post_timestep_modifications(r);
        r->ri_whfast.recalculate_coordinates_this_timestep = 1;
        r->ri_mercurius.recalculate_coordinates_this_timestep = 1;
    }
    
    if (r->N_var){
        reb_var_rescale(r);
    }
    PROFILING_STOP(PROFILING_CAT_INTEGRATOR)

    // Do collisions here. We need both the positions and velocities at the same time.
    // Check for root crossings.
    PROFILING_START()
    reb_boundary_check(r);     
    if (r->tree_needs_update){
        // Update tree (this will remove particles which left the box)
        reb_tree_update(r);          
    }
    PROFILING_STOP(PROFILING_CAT_BOUNDARY)

    // Search for collisions using local and essential tree.
    PROFILING_START()
    reb_collision_search(r);
    PROFILING_STOP(PROFILING_CAT_COLLISION)
    
    // Update walltime
    struct timeval time_end;
    gettimeofday(&time_end,NULL);
    r->walltime += time_end.tv_sec-time_beginning.tv_sec+(time_end.tv_usec-time_beginning.tv_usec)/1e6;
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
        if (r->messages==NULL){
            r->messages = calloc(reb_max_messages_N,sizeof(char*));
        }
        int n = 0;
        for (;n<reb_max_messages_N;n++){
            if (r->messages[n]==NULL){
                break;
            }
        }
        if (n==reb_max_messages_N){
            free(r->messages[0]);
            for (int i=0;i<reb_max_messages_N-1;i++){
                r->messages[i] = r->messages[i+1];
            }
            r->messages[reb_max_messages_N-1] = NULL;
            n= reb_max_messages_N-1;
        }
        r->messages[n] = malloc(sizeof(char*)*reb_max_messages_length);
        r->messages[n][0] = type;
        strcpy(r->messages[n]+1, msg);
    }
}

void reb_warning(struct reb_simulation* const r, const char* const msg){
    reb_message(r, 'w', msg);
}

void reb_error(struct reb_simulation* const r, const char* const msg){
    reb_message(r, 'e', msg);
}

void reb_stop(struct reb_simulation* const r){
    r->status = REB_EXIT_USER;
}

int reb_get_next_message(struct reb_simulation* const r, char* const buf){
    if (r->messages){
        char* w0 = r->messages[0];
        if (w0){
            for(int i=0;i<reb_max_messages_N-1;i++){
                r->messages[i] = r->messages[i+1];
            }
            r->messages[reb_max_messages_N-1] = NULL;
            strcpy(buf,w0);
            free(w0);
            return 1;
        }
    }
    return 0;
}

static int reb_error_message_waiting(struct reb_simulation* const r){
    if (r->messages){
        for (int i=0;i<reb_max_messages_N;i++){
            if (r->messages[i]!=NULL){
                if (r->messages[i][0]=='e'){
                    return 1;
                }
            }
        }
    }
    return 0;
}


void reb_configure_box(struct reb_simulation* const r, const double root_size, const int root_nx, const int root_ny, const int root_nz){
    r->root_size = root_size;
    r->root_nx = root_nx;
    r->root_ny = root_ny;
    r->root_nz = root_nz;
    // Setup box sizes
    r->boxsize.x = r->root_size *(double)r->root_nx;
    r->boxsize.y = r->root_size *(double)r->root_ny;
    r->boxsize.z = r->root_size *(double)r->root_nz;
    r->root_n = r->root_nx*r->root_ny*r->root_nz;
    r->boxsize_max = MAX(r->boxsize.x, MAX(r->boxsize.y, r->boxsize.z));
    if (r->root_nx <=0 || r->root_ny <=0 || r->root_nz <= 0){
        reb_exit("Number of root boxes must be greater or equal to 1 in each direction.");
    }
}
#ifdef MPI
void reb_mpi_init(struct reb_simulation* const r){
    reb_communication_mpi_init(r,0,NULL);
    // Make sure domain can be decomposed into equal number of root boxes per node.
    if ((r->root_n/r->mpi_num)*r->mpi_num != r->root_n){
        if (r->mpi_id==0) fprintf(stderr,"ERROR: Number of root boxes (%d) not a multiple of mpi nodes (%d).\n",r->root_n,r->mpi_num);
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

static void set_dp7_null(struct reb_dp7 * dp){
    dp->p0 = NULL;
    dp->p1 = NULL;
    dp->p2 = NULL;
    dp->p3 = NULL;
    dp->p4 = NULL;
    dp->p5 = NULL;
    dp->p6 = NULL;
}

void reb_free_simulation(struct reb_simulation* const r){
    reb_free_pointers(r);
    free(r);
}

void reb_free_pointers(struct reb_simulation* const r){
    if (r->simulationarchive_filename){
        free(r->simulationarchive_filename);
    }
    reb_tree_delete(r);
    if(r->display_data){
        pthread_mutex_destroy(&(r->display_data->mutex));
        free(r->display_data->r_copy);
        free(r->display_data->particles_copy);
        free(r->display_data->p_jh_copy);
        free(r->display_data->particle_data);
        free(r->display_data->orbit_data);
        free(r->display_data); // TODO: Free other pointers in display_data
    }
    if (r->gravity_cs){
        free(r->gravity_cs  );
    }
    if (r->collisions){
        free(r->collisions  );
    }
    reb_integrator_whfast_reset(r);
    reb_integrator_ias15_reset(r);
    reb_integrator_mercurius_reset(r);
    reb_integrator_bs_reset(r);
    reb_integrator_tes_reset(r);
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
        for (int i=0;i<reb_max_messages_N;i++){
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
    for (int s=0; s<r->odes_N; s++){
        r->odes[s]->r = NULL;
    }
}

void reb_reset_temporary_pointers(struct reb_simulation* const r){
    // Note: this will not clear the particle array.
    r->gravity_cs_allocatedN    = 0;
    r->gravity_cs           = NULL;
    r->collisions_allocatedN    = 0;
    r->collisions           = NULL;
    r->extras               = NULL;
    r->messages             = NULL;
    // ********** Lookup Table
    r->particle_lookup_table = NULL;
    r->N_lookup = 0;
    r->allocatedN_lookup = 0;
    // ********** WHFAST
    r->ri_whfast.allocated_N    = 0;
    r->ri_whfast.allocated_Ntemp= 0;
    r->ri_whfast.p_jh           = NULL;
    r->ri_whfast.p_temp         = NULL;
    r->ri_whfast.keep_unsynchronized = 0;
    // ********** IAS15
    r->ri_ias15.allocatedN      = 0;
    set_dp7_null(&(r->ri_ias15.g));
    set_dp7_null(&(r->ri_ias15.b));
    set_dp7_null(&(r->ri_ias15.csb));
    set_dp7_null(&(r->ri_ias15.e));
    set_dp7_null(&(r->ri_ias15.br));
    set_dp7_null(&(r->ri_ias15.er));
    r->ri_ias15.at          = NULL;
    r->ri_ias15.x0          = NULL;
    r->ri_ias15.v0          = NULL;
    r->ri_ias15.a0          = NULL;
    r->ri_ias15.csx         = NULL;
    r->ri_ias15.csv         = NULL;
    r->ri_ias15.csa0        = NULL;
    r->ri_ias15.at          = NULL;
    r->ri_ias15.map_allocated_N      = 0;
    r->ri_ias15.map         = NULL;
    // ********** MERCURIUS
    r->ri_mercurius.allocatedN = 0;
    r->ri_mercurius.allocatedN_additionalforces = 0;
    r->ri_mercurius.dcrit_allocatedN = 0;
    r->ri_mercurius.dcrit = NULL;
    r->ri_mercurius.particles_backup = NULL;
    r->ri_mercurius.particles_backup_additionalforces = NULL;
    r->ri_mercurius.encounter_map = NULL;
    // ********** JANUS
    r->ri_janus.allocated_N = 0;
    r->ri_janus.p_int = NULL;
    r->ri_janus.recalculate_integer_coordinates_this_timestep = 0;
    r->ri_janus.order = 6;
    r->ri_janus.scale_pos = 1e-16;
    r->ri_janus.scale_vel = 1e-16;
    // ********** ODEs
    r->odes = NULL;
    r->odes_N = 0;
    r->odes_allocatedN = 0;
    // ********** TES
    r->ri_tes.particles_dh = NULL;
}

int reb_reset_function_pointers(struct reb_simulation* const r){
    int wasnotnull = 0;
    if (r->coefficient_of_restitution ||
        r->collision_resolve ||
        r->additional_forces ||
        r->heartbeat ||
        r->display_heartbeat ||
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
    r->display_heartbeat    = NULL;
    r->pre_timestep_modifications  = NULL;
    r->post_timestep_modifications  = NULL;
    r->free_particle_ap = NULL;
    r->extras_cleanup = NULL;
    return wasnotnull;
}

struct reb_simulation* reb_create_simulation(){
    struct reb_simulation* r = calloc(1,sizeof(struct reb_simulation));
    reb_init_simulation(r);
    return r;
}


void reb_copy_simulation_with_messages(struct reb_simulation* r_copy,  struct reb_simulation* r, enum reb_input_binary_messages* warnings){
    char* bufp;
    size_t sizep;
    reb_output_binary_to_stream(r, &bufp,&sizep);
    
    reb_reset_temporary_pointers(r_copy);
    reb_reset_function_pointers(r_copy);
    r_copy->simulationarchive_filename = NULL;
    
    // Set to old version by default. Will be overwritten if new version was used.
    r_copy->simulationarchive_version = 0;

    char* bufp_beginning = bufp; // bufp will be changed
    while(reb_input_field(r_copy, NULL, warnings, &bufp)){ }
    free(bufp_beginning);
    
}

int reb_diff_simulations(struct reb_simulation* r1, struct reb_simulation* r2, int output_option){
    if (output_option!=1 && output_option!=2){
        // Not implemented
        return -1;
    }
    char* bufp1;
    char* bufp2;
    size_t sizep1, sizep2;
    reb_output_binary_to_stream(r1, &bufp1,&sizep1);
    reb_output_binary_to_stream(r2, &bufp2,&sizep2);

    int ret = reb_binary_diff_with_options(bufp1, sizep1, bufp2, sizep2, NULL, NULL, output_option);
    
    free(bufp1);
    free(bufp2);
    return ret;
}

struct reb_simulation* reb_copy_simulation(struct reb_simulation* r){
    struct reb_simulation* r_copy = reb_create_simulation();
    enum reb_input_binary_messages warnings = REB_INPUT_BINARY_WARNING_NONE;
    reb_copy_simulation_with_messages(r_copy,r,&warnings);
    r = reb_input_process_warnings(r, warnings);
    return r_copy;
}

void reb_clear_pre_post_pointers(struct reb_simulation* const r){
    // Temporary fix for REBOUNDx. 
    r->pre_timestep_modifications  = NULL;
    r->post_timestep_modifications  = NULL;
}

void reb_init_simulation(struct reb_simulation* r){
    reb_tools_init_srand(r);
    reb_reset_temporary_pointers(r);
    reb_reset_function_pointers(r);
    r->t        = 0; 
    r->G        = 1;
    r->softening    = 0;
    r->dt       = 0.001;
    r->dt_last_done = 0.;
    r->steps_done = 0;
    r->root_size    = -1;
    r->root_nx  = 1;
    r->root_ny  = 1;
    r->root_nz  = 1;
    r->root_n   = 1;
    r->nghostx  = 0;
    r->nghosty  = 0;
    r->nghostz  = 0;
    r->N        = 0;    
    r->allocatedN   = 0;    
    r->N_active     = -1;   
    r->var_rescale_warning   = 0;   
    r->particle_lookup_table = NULL;
    r->hash_ctr = 0;
    r->N_lookup = 0;
    r->allocatedN_lookup = 0;
    r->testparticle_type = 0;   
    r->testparticle_hidewarnings = 0;
    r->N_var    = 0;    
    r->var_config_N = 0;    
    r->var_config   = NULL;     
    r->exit_min_distance    = 0;    
    r->exit_max_distance    = 0;    
    r->max_radius[0]    = 0.;   
    r->max_radius[1]    = 0.;   
    r->status       = REB_RUNNING;
    r->exact_finish_time    = 1;
    r->force_is_velocity_dependent = 0;
    r->gravity_ignore_terms    = 0;
    r->calculate_megno  = 0;
    r->output_timing_last   = -1;
    r->save_messages = 0;
    r->track_energy_offset = 0;
    r->display_data = NULL;
    r->walltime = 0;

    r->minimum_collision_velocity = 0;
    r->collisions_plog  = 0;
    r->collisions_Nlog  = 0;    
    r->collision_resolve_keep_sorted   = 0;    
    
    r->simulationarchive_size_first    = 0;    
    r->simulationarchive_size_snapshot = 0;    
    r->simulationarchive_version       = 3;
    r->simulationarchive_auto_interval = 0.;    
    r->simulationarchive_auto_walltime = 0.;    
    r->simulationarchive_auto_step     = 0;    
    r->simulationarchive_next          = 0.;    
    r->simulationarchive_next_step     = 0;    
    r->simulationarchive_filename      = NULL;    
    
    // Default modules
#ifdef OPENGL
    r->visualization= REB_VISUALIZATION_OPENGL;
#else // OPENGL
    r->visualization= REB_VISUALIZATION_NONE;
#endif // OPENGL
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
    
    // ********** SABA
    r->ri_saba.type = REB_SABA_10_6_4;
    r->ri_saba.safe_mode = 1;
    r->ri_saba.is_synchronized = 1;
    
    // ********** IAS15
    r->ri_ias15.epsilon         = 1e-9;
    r->ri_ias15.min_dt      = 0;
    r->ri_ias15.epsilon_global  = 1;
    r->ri_ias15.iterations_max_exceeded = 0;    
    
    // ********** SEI
    r->ri_sei.OMEGA     = 1;
    r->ri_sei.OMEGAZ    = -1;
    r->ri_sei.lastdt    = 0;
    
    // ********** MERCURIUS
    r->ri_mercurius.mode = 0;
    r->ri_mercurius.safe_mode = 1;
    r->ri_mercurius.recalculate_coordinates_this_timestep = 0;
    r->ri_mercurius.recalculate_dcrit_this_timestep = 0;
    r->ri_mercurius.is_synchronized = 1;
    r->ri_mercurius.encounterN = 0;
    r->ri_mercurius.hillfac = 3;
    
    // ********** EOS
    r->ri_eos.n = 2;
    r->ri_eos.phi0 = REB_EOS_LF;
    r->ri_eos.phi1 = REB_EOS_LF;
    r->ri_eos.safe_mode = 1;
    r->ri_eos.is_synchronized = 1;
    
    
    // ********** NS
    reb_integrator_bs_reset(r);

    // ********** TES
    r->ri_tes.dq_max = 1e-2;                   // good fall back value and should be OK for reasonable planet masses.
    r->ri_tes.recti_per_orbit = 1.61803398875; // golden ratio 
    r->ri_tes.epsilon = 1e-6;
    r->ri_tes.allocated_N = 0;

    // Tree parameters. Will not be used unless gravity or collision search makes use of tree.
    r->tree_needs_update= 0;
    r->tree_root        = NULL;
    r->opening_angle2   = 0.25;

#ifdef MPI
    r->mpi_id = 0;                            
    r->mpi_num = 0;                           
    r->particles_send = NULL;  
    r->particles_send_N = 0;                  
    r->particles_send_Nmax = 0;               
    r->particles_recv = NULL;     
    r->particles_recv_N = 0;                  
    r->particles_recv_Nmax = 0;               
    
    r->tree_essential_send = NULL;
    r->tree_essential_send_N = 0;             
    r->tree_essential_send_Nmax = 0;          
    r->tree_essential_recv = NULL;
    r->tree_essential_recv_N = 0;             
    r->tree_essential_recv_Nmax = 0;          

#else // MPI
#ifndef LIBREBOUND
    printf("Process id: %d.\n", getpid());
#endif // LIBREBOUND
#endif // MPI
#ifdef OPENMP
    printf("Using OpenMP with %d threads per node.\n",omp_get_max_threads());
#endif // OPENMP
}

int reb_check_exit(struct reb_simulation* const r, const double tmax, double* last_full_dt){
    while(r->status == REB_RUNNING_PAUSED){
        // Wait for user to disable paused simulation
        usleep(1000);
    }
    const double dtsign = copysign(1.,r->dt);   // Used to determine integration direction
    if (reb_error_message_waiting(r)){
        r->status = REB_EXIT_ERROR;
    }
    if (r->status>=0){
        // Exit now.
    }else if(tmax!=INFINITY){
        if(r->exact_finish_time==1){
            if ((r->t+r->dt)*dtsign>=tmax*dtsign){  // Next step would overshoot
                if (r->t==tmax){
                    r->status = REB_EXIT_SUCCESS;
                }else if(r->status == REB_RUNNING_LAST_STEP){
                    double tscale = 1e-12*fabs(tmax);   // Find order of magnitude for time
                    if (tscale<1e-200){     // Failsafe if tmax==0.
                        tscale = 1e-12;
                    }
                    if (fabs(r->t-tmax)<tscale){
                        r->status = REB_EXIT_SUCCESS;
                    }else{
                        // not there yet, do another step.
                        reb_integrator_synchronize(r);
                        r->dt = tmax-r->t;
                    }
                }else{
                    r->status = REB_RUNNING_LAST_STEP; // Do one small step, then exit.
                    reb_integrator_synchronize(r);
                    if (r->dt_last_done!=0.){   // If first timestep is also last, do not use dt_last_done (which would be 0.)
                        *last_full_dt = r->dt_last_done; // store last full dt before decreasing the timestep to match finish time
                    }
                    r->dt = tmax-r->t;
                }
            }else{
                if (r->status == REB_RUNNING_LAST_STEP){
                    // This will get executed if an adaptive integrator reduces
                    // the timestep in what was supposed to be the last timestep.
                    r->status = REB_RUNNING;
                }
            }
        }else{
            if (r->t*dtsign>=tmax*dtsign){  // Past the integration time
                r->status = REB_EXIT_SUCCESS; // Exit now.
            }
        }
    }
#ifndef MPI
    if (!r->N){
        if (!r->odes_N){
            reb_warning(r,"No particles found. Will exit.");
            r->status = REB_EXIT_NOPARTICLES; // Exit now.
        }else{
            if (r->integrator != REB_INTEGRATOR_BS){
                reb_warning(r,"No particles found. Will exit. Use BS integrator to integrate user-defined ODEs without any particles present.");
                r->status = REB_EXIT_NOPARTICLES; // Exit now.
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
    if (r->display_heartbeat){ reb_check_for_display_heartbeat(r); } 
    if (r->exit_max_distance){
        // Check for escaping particles
        const double max2 = r->exit_max_distance * r->exit_max_distance;
        const struct reb_particle* const particles = r->particles;
        const int N = r->N - r->N_var;
        for (int i=0;i<N;i++){
            struct reb_particle p = particles[i];
            double r2 = p.x*p.x + p.y*p.y + p.z*p.z;
            if (r2>max2){
                r->status = REB_EXIT_ESCAPE;
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
                    r->status = REB_EXIT_ENCOUNTER;
                }
            }
        }
    }
}

////////////////////////////////////////////////////
///  Integrate functions and visualization stuff

struct reb_thread_info {
    struct reb_simulation* r;
    double tmax;
};

volatile sig_atomic_t reb_sigint;

void reb_sigint_handler(int signum) {
    // Handles graceful shutdown for interrupts
    if (signum == SIGINT){
        reb_sigint = 1;
    }
}

static void* reb_integrate_raw(void* args){
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
        reb_warning(r,"At least one test particle (type 0) has finite mass. This might lead to unexpected behaviour. Set testparticle_hidewarnings=1 to hide this warning.");
    }

    r->status = REB_RUNNING;
    reb_run_heartbeat(r);
    while(reb_check_exit(r,thread_info->tmax,&last_full_dt)<0){
#ifdef OPENGL
        if (r->display_data){
            if (r->display_data->opengl_enabled){ pthread_mutex_lock(&(r->display_data->mutex)); }
        }
#endif // OPENGL
        if (r->simulationarchive_filename){ reb_simulationarchive_heartbeat(r);}
        reb_step(r); 
        reb_run_heartbeat(r);
        if (reb_sigint== 1){
            r->status = REB_EXIT_SIGINT;
        }
#ifdef OPENGL
        if (r->display_data){
            if (r->display_data->opengl_enabled){ pthread_mutex_unlock(&(r->display_data->mutex)); }
        }
#endif // OPENGL
        if (r->usleep > 0){
            usleep(r->usleep);
        }
    }

    reb_integrator_synchronize(r);
    if (r->display_heartbeat){                          // Display Heartbeat
        r->display_heartbeat(r); 
    }
    if(r->exact_finish_time==1){ // if finish_time = 1, r->dt could have been shrunk, so set to the last full timestep
        r->dt = last_full_dt; 
    }
    if (r->simulationarchive_filename){ reb_simulationarchive_heartbeat(r);}

    return NULL;
}

enum REB_STATUS reb_integrate(struct reb_simulation* const r, double tmax){
    struct reb_thread_info thread_info = {
        .r = r,
        .tmax = tmax, 
    };
    switch (r->visualization){
        case REB_VISUALIZATION_NONE:
            {
                if (r->display_data){
                    r->display_data->opengl_enabled = 0;
                }
                reb_integrate_raw(&thread_info);
            }
            break;
        case REB_VISUALIZATION_OPENGL:
            {
#ifdef OPENGL
                reb_display_init_data(r);
                r->display_data->opengl_enabled = 1;

                pthread_t compute_thread;
                if (pthread_create(&compute_thread,NULL,reb_integrate_raw,&thread_info)){
                    reb_error(r, "Error creating display thread.");
                }
                
                reb_display_init(r); // Display routines running on main thread.

                if (pthread_join(compute_thread,NULL)){
                    reb_error(r, "Error joining display thread.");
                }
#else // OPENGL
                reb_error(r,"REBOUND was not compiled/linked with OPENGL libraries.");
                return REB_EXIT_ERROR; 
#endif // OPENGL
            }
            break;
        case REB_VISUALIZATION_WEBGL:
            {
                reb_display_init_data(r);
                reb_integrate_raw(&thread_info);
            }
            break;
    }
    return r->status;
}

  
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


// Recreate following with 
// grep REB_BINARY_FIELD_TYPE_  src/rebound.h | awk -F 'REB_BINARY_FIELD_TYPE_|=|,' '{print $3, $2}' | xarg
// Add space at beginning and end.
const char* reb_binary_field_type_reverse = 
" 0 T 1 G 2 SOFTENING 3 DT 4 N 5 NVAR 6 VARCONFIGN 7 NACTIVE 8 TESTPARTICLETYPE 9 HASHCTR 10 OPENINGANGLE2 11 STATUS 12 EXACTFINISHTIME 13 FORCEISVELOCITYDEP 14 GRAVITYIGNORETERMS 15 OUTPUTTIMINGLAST 16 SAVEMESSAGES 17 EXITMAXDISTANCE 18 EXITMINDISTANCE 19 USLEEP 20 TRACKENERGYOFFSET 21 ENERGYOFFSET 22 BOXSIZE 23 BOXSIZEMAX 24 ROOTSIZE 25 ROOTN 26 ROOTNX 27 ROOTNY 28 ROOTNZ 29 NGHOSTX 30 NGHOSTY 31 NGHOSTZ 32 COLLISIONRESOLVEKEEPSORTED 33 MINIMUMCOLLISIONVELOCITY 34 COLLISIONSPLOG 35 MAXRADIUS 36 COLLISIONSNLOG 37 CALCULATEMEGNO 38 MEGNOYS 39 MEGNOYSS 40 MEGNOCOVYT 41 MEGNOVART 42 MEGNOMEANT 43 MEGNOMEANY 44 MEGNON 45 SASIZEFIRST 46 SASIZESNAPSHOT 47 SAAUTOINTERVAL 102 SAAUTOWALLTIME 48 SANEXT 50 COLLISION 51 INTEGRATOR 52 BOUNDARY 53 GRAVITY 54 SEI_OMEGA 55 SEI_OMEGAZ 56 SEI_LASTDT 57 SEI_SINDT 58 SEI_TANDT 59 SEI_SINDTZ 60 SEI_TANDTZ 61 WHFAST_CORRECTOR 62 WHFAST_RECALCJAC 63 WHFAST_SAFEMODE 64 WHFAST_KEEPUNSYNC 65 WHFAST_ISSYNCHRON 66 WHFAST_TIMESTEPWARN 69 IAS15_EPSILON 70 IAS15_MINDT 71 IAS15_EPSILONGLOBAL 72 IAS15_ITERATIONSMAX 85 PARTICLES 86 VARCONFIG 87 FUNCTIONPOINTERS 88 IAS15_ALLOCATEDN 89 IAS15_AT 90 IAS15_X0 91 IAS15_V0 92 IAS15_A0 93 IAS15_CSX 94 IAS15_CSV 95 IAS15_CSA0 96 IAS15_G 97 IAS15_B 98 IAS15_CSB 99 IAS15_E 100 IAS15_BR 101 IAS15_ER 104 WHFAST_PJ 107 VISUALIZATION 110 JANUS_ALLOCATEDN 112 JANUS_PINT 113 JANUS_SCALEPOS 114 JANUS_SCALEVEL 115 JANUS_ORDER 116 JANUS_RECALC 117 WHFAST_COORDINATES 118 MERCURIUS_HILLFAC 119 MERCURIUS_SAFEMODE 120 MERCURIUS_ISSYNCHRON 122 MERCURIUS_DCRIT 123 MERCURIUS_RECALCULATE_COORD 125 SAVERSION 126 WALLTIME 130 PYTHON_UNIT_L 131 PYTHON_UNIT_M 132 PYTHON_UNIT_T 133 MERCURIUS_COMPOS 134 MERCURIUS_COMVEL 135 SAAUTOSTEP 136 SANEXTSTEP 137 STEPSDONE 140 SABA_SAFEMODE 141 SABA_ISSYNCHRON 143 WHFAST_CORRECTOR2 144 WHFAST_KERNEL 145 DTLASTDONE 146 SABA_TYPE 147 SABA_KEEPUNSYNC 148 EOS_PHI0 149 EOS_PHI1 150 EOS_N 151 EOS_SAFEMODE 152 EOS_ISSYNCHRON 154 RAND_SEED 155 TESTPARTICLEHIDEWARNINGS 156 BS_EPSABS 157 BS_EPSREL 158 BS_MINDT 159 BS_MAXDT 160 BS_FIRSTORLASTSTEP 161 BS_PREVIOUSREJECTED 162 BS_TARGETITER 163 VARRESCALEWARNING 300 TES_DQ_MAX 301 TES_RECTI_PER_ORBIT 302 TES_EPSILON 303 TES_PERIOD 304 TES_ALLOCATED_N 305 TES_PARTICLES_DH 306 TES_SV_LEN 307 TES_SV_SIZE 308 TES_CV_LEN 309 TES_CV_SIZE 310 TES_MASS 311 TES_X_DH 312 TES_COM 313 TES_COM_DOT 314 TES_MASS_STAR_LAST 320 TES_UVARS_SV_SIZE 321 TES_UVARS_CV_SIZE 322 TES_UVARS_T0 323 TES_UVARS_TLAST 324 TES_UVARS_CSQ 325 TES_UVARS_CSP 326 TES_UVARS_CSV 327 TES_UVARS_Q0 328 TES_UVARS_V0 329 TES_UVARS_P0 330 TES_UVARS_Q1 331 TES_UVARS_V1 332 TES_UVARS_P1 333 TES_UVARS_X 334 TES_UVARS_Q0_NORM 335 TES_UVARS_BETA 336 TES_UVARS_ETA 337 TES_UVARS_ZETA 338 TES_UVARS_PERIOD 339 TES_UVARS_XPERIOD 340 TES_UVARS_STUMPF_C0 341 TES_UVARS_STUMPF_C1 342 TES_UVARS_STUMPF_C2 343 TES_UVARS_STUMPF_C3 344 TES_UVARS_MU 350 TES_RADAU_DX 351 TES_RADAU_XOUT 352 TES_RADAU_RECTI_ARRAY 353 TES_RADAU_PREDICTORS 354 TES_RADAU_DSTATE0 355 TES_RADAU_DDSTATE0 356 TES_RADAU_DSTATE 357 TES_RADAU_DDSTATE 358 TES_RADAU_CS_DSTATE0 359 TES_RADAU_CS_DDSTATE0 360 TES_RADAU_CS_DSTATE 361 TES_RADAU_CS_DDSTATE 362 TES_RADAU_CS_DX 363 TES_RADAU_FCALLS 364 TES_RADAU_RECTIS 365 TES_RADAU_ITERS 366 TES_RADAU_B6 367 TES_RADAU_B 368 TES_RADAU_BLAST 369 TES_RADAU_B_1ST 370 TES_RADAU_BLAST_1ST 371 TES_RADAU_CS_B 372 TES_RADAU_CS_B_1ST 373 TES_RADAU_G 374 TES_RADAU_G_1ST 380 TES_DHEM_XOSC_STORE 381 TES_DHEM_XOSC_PRED_STORE 382 TES_DHEM_XOSC_CS_STORE 383 TES_DHEM_XOSC_DOT_STORE 384 TES_DHEM_X 385 TES_DHEM_M_INV 386 TES_DHEM_M_TOTAL 387 TES_DHEM_RECTI_TIME 388 TES_DHEM_RECTI_PERIOD 390 WHFAST512_KEEPUNSYNC 391 WHFAST512_ISSYNCHRON 392 WHFAST512_GRPOTENTIAL 393 WHFAST512_ALLOCATEDN 394 WHFAST512_PJH 395 WHFAST512_PJH0 1329743186 HEADER 9998 SABLOB 9999 END ";
