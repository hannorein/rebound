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
#include <signal.h>
#include <string.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/wait.h>
#include <pthread.h>
#include <fcntl.h>
#include "rebound.h"
#include "integrator.h"
#include "integrator_whfast.h"
#include "integrator_whfasthelio.h"
#include "integrator_ias15.h"
#include "integrator_hermes.h"
#include "boundary.h"
#include "gravity.h"
#include "collision.h"
#include "tree.h"
#include "output.h"
#include "tools.h"
#include "particle.h"
#include "simulationarchive.h"
#ifdef MPI
#include "communication_mpi.h"
#endif
#ifdef OPENGL
#include "display.h"
#endif // OPENGL
#ifdef OPENMP
#include <omp.h>
#endif
#define MAX(a, b) ((a) < (b) ? (b) : (a))       ///< Returns the maximum of a and b
#define STRINGIFY(s) str(s)
#define str(s) #s

const int reb_max_messages_length = 1024;   // needs to be constant expression for array size
const int reb_max_messages_N = 10;
const char* reb_build_str = __DATE__ " " __TIME__;  // Date and time build string. 
const char* reb_version_str = "3.1.1";         // **VERSIONLINE** This line gets updated automatically. Do not edit manually.
const char* reb_githash_str = STRINGIFY(GITHASH);             // This line gets updated automatically. Do not edit manually.

void reb_step(struct reb_simulation* const r){
    // A 'DKD'-like integrator will do the first 'D' part.
    PROFILING_START()
    reb_integrator_part1(r);
    PROFILING_STOP(PROFILING_CAT_INTEGRATOR)

    // Update and simplify tree. 
    // Prepare particles for distribution to other nodes. 
    // This function also creates the tree if called for the first time.
    if (r->tree_needs_update || r->gravity==REB_GRAVITY_TREE || r->collision==REB_COLLISION_TREE){
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
        r->ri_whfast.recalculate_jacobi_this_timestep = 1;
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
    if (r->integrator!=REB_INTEGRATOR_HERMES){ //Hybrid integrator will search for collisions in mini simulation.
        reb_collision_search(r);
    }
    PROFILING_STOP(PROFILING_CAT_COLLISION)
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
    reb_tree_delete(r);
    free(r->gravity_cs  );
    free(r->collisions  );
    reb_integrator_whfast_reset(r);
    reb_integrator_whfasthelio_reset(r);
    reb_integrator_ias15_reset(r);
    free(r->particles   );
    free(r->particle_lookup_table);
    if (r->messages){
        for (int i=0;i<reb_max_messages_N;i++){
            free(r->messages[i]);
        }
    }
    free(r->messages);
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
    r->ri_whfast.eta            = NULL;
    r->ri_whfast.p_j            = NULL;
    r->ri_whfast.keep_unsynchronized = 0;
    // ********** WHFASTHELIO
    r->ri_whfasthelio.allocated_N  = 0;
    r->ri_whfasthelio.p_h          = NULL;
    r->ri_whfasthelio.keep_unsynchronized = 0;
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
    // ********** HERMES
    r->ri_hermes.mini      = NULL;
    r->ri_hermes.global    = NULL;
    r->ri_hermes.global_index_from_mini_index = NULL;
    r->ri_hermes.global_index_from_mini_index_N = 0;
    r->ri_hermes.global_index_from_mini_index_Nmax = 0;
    r->ri_hermes.is_in_mini = NULL;
    r->ri_hermes.is_in_mini_Nmax = 0;
    r->ri_hermes.a_Nmax = 0;
    r->ri_hermes.a_i = NULL;
    r->ri_hermes.a_f = NULL;
}

int reb_reset_function_pointers(struct reb_simulation* const r){
    int wasnotnull = 0;
    if (r->coefficient_of_restitution ||
        r->collision_resolve ||
        r->additional_forces ||
        r->heartbeat ||
        r->post_timestep_modifications ||
        r->free_particle_ap){
      wasnotnull = 1;
    }
    r->coefficient_of_restitution   = NULL;
    r->collision_resolve        = NULL;
    r->additional_forces        = NULL;
    r->heartbeat            = NULL;
    r->post_timestep_modifications  = NULL;
    r->free_particle_ap = NULL;
    return wasnotnull;
}

struct reb_simulation* reb_create_simulation(){
    struct reb_simulation* r = calloc(1,sizeof(struct reb_simulation));
    reb_init_simulation(r);
    return r;
}

void reb_init_simulation(struct reb_simulation* r){
    reb_tools_init_srand();
    reb_reset_temporary_pointers(r);
    reb_reset_function_pointers(r);
    r->t        = 0; 
    r->G        = 1;
    r->softening    = 0;
    r->dt       = 0.001;
    r->dt_last_done = 0.;
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
    r->particle_lookup_table = NULL;
    r->hash_ctr = 0;
    r->N_lookup = 0;
    r->allocatedN_lookup = 0;
    r->testparticle_type = 0;   
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

    r->minimum_collision_velocity = 0;
    r->collisions_plog  = 0;
    r->collisions_Nlog  = 0;    
    r->collision_resolve_keep_sorted  = 0;    
    
    r->simulationarchive_size_first  = 0;    
    r->simulationarchive_size_snapshot   = 0;    
    r->simulationarchive_interval    = 0.;    
    r->simulationarchive_interval_walltime = 0.;    
    r->simulationarchive_walltime    = 0.;    
    r->simulationarchive_next        = 0.;    
    r->simulationarchive_filename    = NULL;    
    
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
    r->ri_whfast.safe_mode = 1;
    r->ri_whfast.recalculate_jacobi_this_timestep = 0;
    r->ri_whfast.is_synchronized = 1;
    r->ri_whfast.timestep_warning = 0;
    r->ri_whfast.recalculate_jacobi_but_not_synchronized_warning = 0;
    // ********** WHFASTHELIO
    r->ri_whfasthelio.corrector = 0;
    r->ri_whfasthelio.safe_mode = 1;
    r->ri_whfasthelio.is_synchronized = 1;
    r->ri_whfasthelio.recalculate_heliocentric_this_timestep = 0;
    r->ri_whfasthelio.recalculate_heliocentric_but_not_synchronized_warning = 0;
    
    // ********** IAS15
    r->ri_ias15.epsilon         = 1e-9;
    r->ri_ias15.min_dt      = 0;
    r->ri_ias15.epsilon_global  = 1;
    r->ri_ias15.iterations_max_exceeded = 0;    
    
    // ********** SEI
    r->ri_sei.OMEGA     = 1;
    r->ri_sei.OMEGAZ    = -1;
    r->ri_sei.lastdt    = 0;
    
    // ********** HERMES
    r->ri_hermes.mini_active = 0;
    r->ri_hermes.collision_this_global_dt = 0;
    r->ri_hermes.steps = 0;
    r->ri_hermes.steps_miniactive = 0;
    r->ri_hermes.steps_miniN = 0;
    r->ri_hermes.timestep_too_large_warning = 0;
    r->ri_hermes.solar_switch_factor = 15.;
    r->ri_hermes.hill_switch_factor = 3.;            
    r->ri_hermes.adaptive_hill_switch_factor = 1;    
    r->ri_hermes.current_hill_switch_factor = 3.;     //Internal 
    
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
    if (r->status>=0){
        // Exit now.
    }else if(tmax!=INFINITY){
        if(r->exact_finish_time==1){
            if ((r->t+r->dt)*dtsign>=tmax*dtsign){  // Next step would overshoot
                double tscale = 1e-12*fabs(tmax);   // Find order of magnitude for time
                if (tscale<1e-200){     // Failsafe if tmax==0.
                    tscale = 1e-12;
                }
                if (r->t==tmax){
                    r->status = REB_EXIT_SUCCESS;
                }else if(r->status == REB_RUNNING_LAST_STEP){
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
    if (r->N<=0){
        reb_warning(r,"No particles found. Will exit.");
        r->status = REB_EXIT_NOPARTICLES; // Exit now.
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
    if (r->simulationarchive_filename){ reb_simulationarchive_heartbeat(r);}
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
    if (r->usleep > 0){
        usleep(r->usleep);
    }
}

////////////////////////////////////////////////////
///  Integrate functions and visualization stuff


void * reb_integrate_without_visualization(void* args){
#ifdef MPI
    // Distribute particles
    reb_communication_mpi_distribute_particles(r);
#endif // MPI
    struct reb_display_data* data = (struct reb_display_data*)args;
    struct reb_simulation* r = data->r;
    double tmax = data->tmax;
#ifdef OPENGL
    pthread_mutex_t* display_mutex = data->mutex;
    unsigned int opengl_enabled = data->opengl_enabled;
#endif //: OPENGL

#ifndef LIBREBOUND
    struct timeval tim;
    gettimeofday(&tim, NULL);
    double timing_initial = tim.tv_sec+(tim.tv_usec/1000000.0);
#endif // LIBREBOUND

    double last_full_dt = r->dt; // need to store r->dt in case timestep gets artificially shrunk to meet exact_finish_time=1

    r->status = REB_RUNNING;
    reb_run_heartbeat(r);
    while(reb_check_exit(r,tmax,&last_full_dt)<0){
#ifdef OPENGL
        if (opengl_enabled){ pthread_mutex_lock(display_mutex); }
#endif // OPENGL
        reb_step(r); 
        reb_run_heartbeat(r);
#ifdef OPENGL
        if (opengl_enabled){ pthread_mutex_unlock(display_mutex); }
#endif // OPENGL
    }

    reb_integrator_synchronize(r);
    if(r->exact_finish_time==1){ // if finish_time = 1, r->dt could have been shrunk, so set to the last full timestep
        r->dt = last_full_dt; 
    }


#ifndef LIBREBOUND
    gettimeofday(&tim, NULL);
    double timing_final = tim.tv_sec+(tim.tv_usec/1000000.0);
    printf("\nComputation finished. Total runtime: %f s\n",timing_final-timing_initial);
#endif // LIBREBOUND
    data->return_status = r->status;
    return NULL;
}

enum REB_STATUS reb_integrate(struct reb_simulation* const r, double tmax){
#ifdef OPENGL
    int opengl_enabled = (r->usleep<0)?0:1;
    if (!opengl_enabled){
#endif // OPENGL
        struct reb_display_data data = { 
            .r = r, 
#ifdef OPENGL
            .opengl_enabled = opengl_enabled, 
            .mutex = NULL, 
#endif // OPENGL
            .tmax = tmax, 
            .return_status = REB_RUNNING,
        };
        reb_integrate_without_visualization(&data);
        return data.return_status;
#ifdef OPENGL
    }else{
        // Need root_size for visualization. Creating one. 
        if (r->root_size==-1){  
            reb_warning(r,"Configuring box automatically for vizualization based on particle positions.");
            const struct reb_particle* p = r->particles;
            double max_r = 0;
            for (int i=0;i<r->N;i++){
                const double _r = sqrt(p[i].x*p[i].x+p[i].y*p[i].y+p[i].z*p[i].z);
                max_r = MAX(max_r, _r);
            }
            reb_configure_box(r, max_r*2.3,MAX(1.,r->root_nx),MAX(1.,r->root_ny),MAX(1.,r->root_nz));
        }

        pthread_mutex_t mutex;
        if (pthread_mutex_init(&mutex, NULL)){
            reb_error(r,"Mutex creation failed.");
        }
        
        struct reb_display_data data = { 
            .r = r, 
            .opengl_enabled = opengl_enabled, 
            .mutex = &mutex, 
            .tmax = tmax, 
            .return_status = REB_RUNNING,
        };
        
        pthread_t compute_thread;
        if (pthread_create(&compute_thread,NULL,reb_integrate_without_visualization,&data)){
            reb_error(r, "Error creating display thread.");
        }

        reb_display_init(&data);

        if (pthread_join(compute_thread,NULL)){
            reb_error(r, "Error joining display thread.");
        }
        pthread_mutex_destroy(&mutex);
        return data.return_status;
    }
#endif //OPENGL
}

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
