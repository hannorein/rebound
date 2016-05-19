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
#include <semaphore.h>
#include <fcntl.h>
#include "rebound.h"
#include "integrator.h"
#include "integrator_wh.h"
#include "integrator_whfast.h"
#include "integrator_ias15.h"
#include "boundary.h"
#include "gravity.h"
#include "collision.h"
#include "tree.h"
#include "output.h"
#include "tools.h"
#include "particle.h"
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

#ifndef LIBREBOUND
static const char* logo[];              /**< Logo of rebound. */
#endif // LIBREBOUND
const char* reb_build_str = __DATE__ " " __TIME__;  // Date and time build string. 
const char* reb_version_str = "2.16.3";         // **VERSIONLINE** This line gets updated automatically. Do not edit manually.


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
    reb_collision_search(r);
    PROFILING_STOP(PROFILING_CAT_COLLISION)
}

void reb_exit(const char* const msg){
    // This function should also kill all children. 
    // Not implemented as pid is not easy to get to.
    // kill(pid, SIGKILL);
    fprintf(stderr,"\n\033[1mError!\033[0m %s\n",msg);
    exit(EXIT_FAILURE);
}

void reb_warning(const char* const msg){
    fprintf(stderr,"\n\033[1mWarning!\033[0m %s\n",msg);
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
    reb_integrator_wh_reset(r);
    reb_integrator_whfast_reset(r);
    reb_integrator_ias15_reset(r);
    free(r->particles   );
}

void reb_reset_temporary_pointers(struct reb_simulation* const r){
    // Note: this will not clear the particle array.
    r->gravity_cs_allocatedN    = 0;
    r->gravity_cs           = NULL;
    r->collisions_allocatedN    = 0;
    r->collisions           = NULL;
    // ********** WHFAST
    r->ri_whfast.allocated_N    = 0;
    r->ri_whfast.eta        = NULL;
    r->ri_whfast.p_j        = NULL;
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
    // ********** WH
    r->ri_wh.allocatedN         = 0;
    r->ri_wh.eta            = NULL;
}

void reb_reset_function_pointers(struct reb_simulation* const r){
    r->coefficient_of_restitution   = NULL;
    r->collision_resolve        = NULL;
    r->additional_forces        = NULL;
    r->heartbeat            = NULL;
    r->post_timestep_modifications  = NULL;
}

struct reb_simulation* reb_create_simulation(){
    struct reb_simulation* r = calloc(1,sizeof(struct reb_simulation));
    reb_init_simulation(r);
    return r;
}

void reb_init_simulation(struct reb_simulation* r){
#ifndef LIBREBOUND
    int i =0;
    while (logo[i]!=NULL){ printf("%s",logo[i++]); }
    printf("Built: %s\n\n",reb_build_str);
#endif // LIBREBOUND
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
    r->gravity_ignore_10    = 0;
    r->calculate_megno  = 0;
    r->output_timing_last   = -1;

    r->minimum_collision_velocity = 0;
    r->collisions_plog  = 0;
    r->collisions_Nlog  = 0;    
    
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
    
    // ********** IAS15
    r->ri_ias15.epsilon         = 1e-9;
    r->ri_ias15.min_dt      = 0;
    r->ri_ias15.epsilon_global  = 1;
    r->ri_ias15.iterations_max_exceeded = 0;    
    
    // ********** SEI
    r->ri_sei.OMEGA     = 1;
    r->ri_sei.OMEGAZ    = -1;
    r->ri_sei.lastdt    = 0;
    
    r->ri_hybrid.switch_ratio = 8; // Default of 8 mutual Hill radii
    r->ri_hybrid.mode = SYMPLECTIC;

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
        reb_warning("No particles found. Will exit.");
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
enum REB_STATUS reb_integrate(struct reb_simulation* const r_user, double tmax){
#ifdef MPI
    // Distribute particles
    reb_communication_mpi_distribute_particles(r_user);
#endif // MPI
#ifdef OPENGL
    int opengl_enabled = 1;
    if (r_user->usleep<0){
        opengl_enabled = 0;
    }

    // Copy and share simulation struct 
    struct reb_simulation* r = NULL;
    char sem_name[256] = "reb_display";
    sem_t* display_mutex = NULL;
   
    if (opengl_enabled){
        r = (struct reb_simulation*)mmap(r_user, sizeof(struct reb_simulation), PROT_READ|PROT_WRITE, MAP_ANON|MAP_SHARED, -1, 0);
        memcpy(r, r_user, sizeof(struct reb_simulation));
        // Copy and share particle array
        r->particles = (struct reb_particle*)mmap(NULL, r->allocatedN*sizeof(struct reb_particle), PROT_READ|PROT_WRITE, MAP_ANON|MAP_SHARED, -1, 0);
        memcpy(r->particles, r_user->particles, r->allocatedN*sizeof(struct reb_particle));
        for(int i=0; i<r->N; i++){
            r->particles[i].sim = r;
        }
        // Create Semaphore
#ifdef MPI
        sprintf(sem_name,"/reb_display_%d;",r_user->mpi_id);
#endif // MPI
        sem_unlink(sem_name); // unlink first
        if ((display_mutex = sem_open(sem_name, O_CREAT | O_EXCL, 0666, 1))==SEM_FAILED){
            perror("sem_open");
            exit(EXIT_FAILURE);
        }
        // Need root_size for visualization. Creating one. 
        if (r->root_size==-1){  
            reb_warning("Configuring box automatically for vizualization based on particle positions.");
            const struct reb_particle* p = r->particles;
            double max_r = 0;
            for (int i=0;i<r->N;i++){
                const double _r = sqrt(p[i].x*p[i].x+p[i].y*p[i].y+p[i].z*p[i].z);
                max_r = MAX(max_r, _r);
            }
            reb_configure_box(r, max_r*2.3,MAX(1.,r->root_nx),MAX(1.,r->root_ny),MAX(1.,r->root_nz));
        }
    }else{
        r = r_user;
    }

#else // OPENGL
    struct reb_simulation* const r = r_user;
#endif // OPENGL



#ifndef LIBREBOUND
    struct timeval tim;
    gettimeofday(&tim, NULL);
    double timing_initial = tim.tv_sec+(tim.tv_usec/1000000.0);
#endif // LIBREBOUND

    double last_full_dt = r->dt; // need to store r->dt in case timestep gets artificially shrunk to meet exact_finish_time=1

    r->status = REB_RUNNING;
    reb_run_heartbeat(r);

    
#ifdef OPENGL
    if (opengl_enabled){
        pid_t   childpid;
        if((childpid = fork()) == -1) {
                perror("fork");
                exit(EXIT_FAILURE);
        }
        if(childpid == 0) {  // Child (vizualization)
            reb_display_init(0,NULL,r, display_mutex);
            exit(EXIT_SUCCESS); // NEVER REACHED
        } else {        // Parent (computation)
            PROFILING_START()
            while(reb_check_exit(r,tmax,&last_full_dt)<0){
                sem_wait(display_mutex);    
                PROFILING_STOP(PROFILING_CAT_VISUALIZATION)
                reb_step(r);
                reb_run_heartbeat(r);
                PROFILING_START()
                sem_post(display_mutex);
            }
            PROFILING_STOP(PROFILING_CAT_VISUALIZATION)
        }
    }else{
#endif // OPENGL
        while(reb_check_exit(r,tmax,&last_full_dt)<0){
            reb_step(r); 
            reb_run_heartbeat(r);
        }
#ifdef OPENGL
    }
#endif // OPENGL

    reb_integrator_synchronize(r);
    if(r->exact_finish_time==1){ // if finish_time = 1, r->dt could have been shrunk, so set to the last full timestep
        r->dt = last_full_dt; 
    }

#ifdef OPENGL
    if (opengl_enabled){
        int status;
        wait(&status);
        struct reb_particle* const particles_user_loc = r_user->particles;
        memcpy(r_user, r, sizeof(struct reb_simulation));
        r_user->particles = particles_user_loc;
        memcpy(r_user->particles, r->particles, r->N*sizeof(struct reb_particle));
        for(int i=0; i<r_user->N; i++){
            r_user->particles[i].sim = r_user;
        }
        sem_unlink(sem_name);
        sem_close(display_mutex);
    }
#endif //OPENGL

#ifndef LIBREBOUND
    gettimeofday(&tim, NULL);
    double timing_final = tim.tv_sec+(tim.tv_usec/1000000.0);
    printf("\nComputation finished. Total runtime: %f s\n",timing_final-timing_initial);
#endif // LIBREBOUND
    return r->status;
}

#ifndef LIBREBOUND
static const char* logo[] = {
"          _                           _  \n",   
"         | |                         | | \n",  
" _ __ ___| |__   ___  _   _ _ __   __| | \n", 
"| '__/ _ \\ '_ \\ / _ \\| | | | '_ \\ / _` | \n", 
"| | |  __/ |_) | (_) | |_| | | | | (_| | \n", 
"|_|  \\___|_.__/ \\___/ \\__,_|_| |_|\\__,_| \n", 
"                                         \n",   
"              `-:://::.`                 \n",
"          `/oshhoo+++oossso+:`           \n", 
"       `/ssooys++++++ossssssyyo:`        \n", 
"     `+do++oho+++osssso++++++++sy/`      \n", 
"    :yoh+++ho++oys+++++++++++++++ss.     \n", 
"   /y++hooyyooshooo+++++++++++++++oh-    \n", 
"  -dsssdssdsssdssssssssssooo+++++++oh`   \n", 
"  ho++ys+oy+++ho++++++++oosssssooo++so   \n", 
" .d++oy++ys+++oh+++++++++++++++oosssod   \n", 
" -h+oh+++yo++++oyo+++++++++++++++++oom   \n", 
" `d+ho+++ys+++++oys++++++++++++++++++d   \n", 
"  yys++++oy+++++++oys+++++++++++++++s+   \n", 
"  .m++++++h+++++++++oys++++++++++++oy`   \n", 
"   -yo++++ss++++++++++oyso++++++++oy.    \n", 
"    .ss++++ho+++++++++++osys+++++yo`     \n", 
"      :ss+++ho+++++++++++++osssss-       \n", 
"        -ossoys++++++++++++osso.         \n", 
"          `-/oyyyssosssyso+/.            \n", 
"                ``....`                  \n", 
"                                         \n",   
"Written by Hanno Rein, Shangfei Liu,     \n",
"David Spiegel, Daniel Tamayo and many    \n",
"other. REBOUND project website:          \n",  
"http://github.com/hannorein/rebound/     \n",    
"                                         \n", 
NULL};
#endif // LIBREBOUND
