/**
 * @file 	particle.c
 * @brief 	reb_particle structure and main particle routines.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
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
#include <stdint.h>
#include <string.h>
#include "rebound.h"
#include "tree.h"
#include "boundary.h"
#include "particle.h"
#include "integrator_ias15.h"
#include "integrator_bs.h"
#include "integrator_mercurius.h"
#include "integrator_trace.h"
#include "collision.h"
#ifdef OPENMP
#include <omp.h>
#endif
#ifdef MPI
#include "communication_mpi.h"
#endif // MPI

#ifdef GRAVITY_GRAPE
#warning Fix this. 
extern double gravity_minimum_mass;
#endif // GRAVITY_GRAPE

static void reb_simulation_add_local(struct reb_simulation* const r, struct reb_particle pt){
    if (reb_boundary_particle_is_in_box(r, pt)==0){
        if (r->boxsize.x==0 && r->boxsize.y==0 && r->boxsize.z==0){ 
            reb_simulation_error(r,"Cannot add particle because simulation box not initialized. Call reb_simulation_configure_box() before adding particles.");
        }else{
            // reb_particle has left the box. Do not add.
            reb_simulation_error(r,"Particle outside of box boundaries. Did not add particle.");
        }
        return;
    }
    while (r->N_allocated<=r->N){
        unsigned int old_N_allocated = r->N_allocated;
        r->N_allocated = r->N_allocated ? r->N_allocated * 2 : 128;
        r->particles = realloc(r->particles,sizeof(struct reb_particle)*r->N_allocated);
        memset(r->particles + old_N_allocated, 0, (r->N_allocated - old_N_allocated) * sizeof(struct reb_particle));
    }

    r->particles[r->N] = pt;
    r->particles[r->N].sim = r;
    if (r->gravity==REB_GRAVITY_TREE || r->collision==REB_COLLISION_TREE || r->collision==REB_COLLISION_LINETREE){
        if (r->root_size==-1){
            reb_simulation_error(r,"root_size is -1. Make sure you call reb_simulation_configure_box() before using a tree based gravity or collision solver.");
            return;
        }
        if(fabs(pt.x)>r->boxsize.x/2. || fabs(pt.y)>r->boxsize.y/2. || fabs(pt.z)>r->boxsize.z/2.){
            reb_simulation_error(r,"Cannot add particle outside of simulation box.");
            return;
        }
        reb_tree_add_particle_to_tree(r, r->N);
    }
    (r->N)++;
    // Particle was added successfully. Do other work now.
    if (pt.name){
        r->particles[r->N-1].name = reb_simulation_register_name(r,pt.name);
    }
    if (r->integrator == REB_INTEGRATOR_MERCURIUS){
        struct reb_integrator_mercurius* rim = &(r->ri_mercurius);
        if (r->ri_mercurius.mode==0){ //WHFast part
            rim->recalculate_r_crit_this_timestep       = 1;
            rim->recalculate_coordinates_this_timestep = 1;
        }else{  // IAS15 part
            reb_integrator_ias15_reset(r);
            if (rim->N_allocated_dcrit<r->N){
                rim->dcrit              = realloc(rim->dcrit, sizeof(double)*r->N);
                rim->N_allocated_dcrit = r->N;
            }
            rim->dcrit[r->N-1] = reb_integrator_mercurius_calculate_dcrit_for_particle(r,r->N-1);
            if (rim->N_allocated<r->N){
                rim->particles_backup   = realloc(rim->particles_backup,sizeof(struct reb_particle)*r->N);
                rim->encounter_map      = realloc(rim->encounter_map,sizeof(int)*r->N);
                rim->N_allocated = r->N;
            }
            rim->encounter_map[rim->encounter_N] = r->N-1;
            rim->encounter_N++;
            if (r->N_active==-1){ 
                // If global N_active is not set, then all particles are active, so the new one as well.
                // Otherwise, assume we're adding non active particle. 
                rim->encounter_N_active++;
            }
        }
    }

    // TRACE can add particles mid-timestep now
    if (r->integrator == REB_INTEGRATOR_TRACE){
        struct reb_integrator_trace* ri_trace = &(r->ri_trace);
        if (r->ri_trace.mode==1 || r->ri_trace.mode==3){ // BS part
            const int old_N = r->N-1;
            if (ri_trace->N_allocated < r->N){
                ri_trace->current_Ks    = realloc(ri_trace->current_Ks, sizeof(int)*r->N*r->N);
                ri_trace->particles_backup = realloc(ri_trace->particles_backup, sizeof(struct reb_particle)*r->N);
                ri_trace->particles_backup_kepler = realloc(ri_trace->particles_backup_kepler, sizeof(struct reb_particle)*r->N);
                ri_trace->current_Ks    = realloc(ri_trace->current_Ks, sizeof(int)*r->N*r->N);
                ri_trace->encounter_map = realloc(ri_trace->encounter_map, sizeof(int)*r->N);
                ri_trace->N_allocated   = r->N;
            }

            // First reshuffle existing Ks
            for (int i = old_N-1; i >= 0; i--){
                for (int j = old_N-1; j >= 0; j--){ 
                    ri_trace->current_Ks[i*old_N+j+i] = ri_trace->current_Ks[i*old_N+j];
                }
            }

            // add in new particle, we want it to interact with all currently interacting particles
            // exclude star
            for (int i = 1; i < ri_trace->encounter_N; i++){
                ri_trace->current_Ks[ri_trace->encounter_map[i]*r->N+old_N] = 1;
            }

            ri_trace->encounter_map[ri_trace->encounter_N] = old_N;
            ri_trace->encounter_N++;

            if (r->N_active==-1){ 
                // If global N_active is not set, then all particles are active, so the new one as well.
                // Otherwise, assume we're adding non active particle. 
                ri_trace->encounter_N_active++;
            }

        }
    }
}

void reb_simulation_add(struct reb_simulation* const r, struct reb_particle pt){
#ifdef GRAVITY_GRAPE
    if (pt.m<gravity_minimum_mass){
        gravity_minimum_mass = pt.m;
    }
#endif // GRAVITY_GRAPE
#ifdef MPI
    int rootbox = reb_get_rootbox_for_particle(r, pt);
    int N_root_per_node = r->N_root/r->mpi_num;
    int proc_id = rootbox/N_root_per_node;
    const unsigned int N_active = (r->N_active==-1)?r->N: (unsigned int)r->N_active;
    if (proc_id != r->mpi_id && r->N >= N_active){
        // Add particle to array and send them to proc_id later. 
        reb_communication_mpi_add_particle_to_send_queue(r,pt,proc_id);
        return;
    }
#endif // MPI
       // Add particle to local particle array.
    reb_simulation_add_local(r, pt);
}

int reb_particle_check_testparticles(struct reb_simulation* const r){
    if (r->N_active == (int)r->N || r->N_active == -1){
        return 0;
    }
    // Check if testparticle of type 0 has mass!=0
    if (r->testparticle_type == 0){
        int found_issue = 0;
        const int N_real = r->N - r->N_var;
#pragma omp parallel for
        for (int i=r->N_active; i<N_real; i++){
            if (r->particles[i].m!=0.){
                found_issue = 1;
            }
        }
        if (found_issue){
            return 1;
        }
    }
    return 0;
}

// Finds the two largest particles in the simulation. *p1 and *p2 will be set to the indices of the particles.
void reb_simulation_two_largest_particles(struct reb_simulation* r, int* p1, int* p2) {
    struct reb_particle* particles = r->particles;
    *p1 = -1;
    *p2 = -1;
    double largest1 = -1.0;
    double largest2 = -1.0;
#ifdef OPENMP
    int num_threads;
    // A struct to hold the two largest values found by each thread
    struct two_max {
        double largest1;
        double largest2;
        int p1;
        int p2;
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
        thread_max[thread_id].p1 = -1;
        thread_max[thread_id].p2 = -1;

#pragma omp for
        for (int i=0; i<r->N; i++) {
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
    for (int i=0; i<r->N; i++) {
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


int reb_get_rootbox_for_particle(const struct reb_simulation* const r, struct reb_particle pt){
    if (r->root_size==-1) return 0;
    int i = ((int)floor((pt.x + r->boxsize.x/2.)/r->root_size)+r->N_root_x)%r->N_root_x;
    int j = ((int)floor((pt.y + r->boxsize.y/2.)/r->root_size)+r->N_root_y)%r->N_root_y;
    int k = ((int)floor((pt.z + r->boxsize.z/2.)/r->root_size)+r->N_root_z)%r->N_root_z;
    int index = (k*r->N_root_y+j)*r->N_root_x+i;
    return index;
}

int reb_simulation_particle_index(struct reb_particle* p){
    struct reb_simulation* r = p->sim;
    int i = 0;
    const int N = r->N;
    while(&r->particles[i] != p){
        i++;
        if(i>=N){
            return -1;	// p not in simulation.  Shouldn't happen unless you mess with p.sim after creating the particle
        }	
    }
    return i;
}


const char* reb_simulation_get_registered_name(struct reb_simulation* r, const char* const name){
#ifdef MPI
    return name; // Does not register.
#else // MPI
    if (name==NULL) return NULL; // NULL string not allowed.
    for (int i=0; i<r->N_name_list; i++){
        if (strcmp(name,r->name_list[i])==0){
            return r->name_list[i];
        }
    }
    return NULL; // Not found
#endif // MPI
}
const char* reb_simulation_register_name(struct reb_simulation* r, const char* const name){
    const char* registered_name = reb_simulation_get_registered_name(r,name);
    if (registered_name) return registered_name;
    registered_name = strdup(name);
    r->N_name_list++;
    r->name_list = realloc(r->name_list,sizeof(char*)*r->N_name_list);
    r->name_list[r->N_name_list-1] = (char*)registered_name;
    return registered_name;
}
struct reb_particle* reb_simulation_get_particle_by_name(struct reb_simulation* r, const char* const name){
    for (int i=0; i<r->N; i++){
        const char* p_name = r->particles[i].name;
        if (p_name){
            if (strcmp(p_name,name)==0){
                return &(r->particles[i]);
            }
        }
    }
    return NULL; // Not found
}


#ifdef MPI
struct reb_particle* reb_simulation_get_particle_by_id(struct reb_simulation* r, int id){
    for (int i=0; i<r->N; i++){
        int p_id = (int)(uintptr_t)(r->particles[i].name);
        if (p_id){
            if (p_id == id){
                return &(r->particles[i]);
            }
        }
    }
    return NULL; // Not found
}

struct reb_particle reb_simulation_particle_by_id_mpi(struct reb_simulation* const r, int id){
    struct reb_particle* p = reb_simulation_get_particle_by_id(r, id);
    int found = (p==NULL)?0:1;
    MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (found == 0){
        return reb_particle_nan();
    }
    if (found > 1){
        reb_simulation_error(r, "Multiple particles with same id found.");
        return reb_particle_nan();
    }
    struct reb_particle ph = {0};
    if (p!=NULL){
        ph = *p;
        ph.sim = NULL;
    }
    int root = (p==NULL) ? 0 : r->mpi_id;
    MPI_Allreduce(MPI_IN_PLACE, &root, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Bcast(&ph, sizeof(struct reb_particle), MPI_CHAR, root, MPI_COMM_WORLD);
    return ph;
}
#endif // MPI

int reb_simulation_remove_particle_by_name(struct reb_simulation* r, const char* const name, int keep_sorted){
    struct reb_particle* p = reb_simulation_get_particle_by_name(r, name);
    if (!p){
        reb_simulation_error(r, "Particle not found.");
        return 1; // Not found.
    }
    size_t index = p - r->particles;
    return !reb_simulation_remove_particle(r, index, keep_sorted); // TODO: return value is different between the two functions. 
}

void reb_particle_set_name(struct reb_particle* p, const char* const name){
    if (name==NULL){
        // Delete name.
        p->name = NULL;
        return;
    }
    struct reb_simulation* r = p->sim;
    if (!r){
        reb_simulation_error(NULL,"Cannot set particle name using_reb_particle_set_name() as the particle is not part of a simulation. You can set the name manually.");
        return;
    }
    p->name = reb_simulation_register_name(r,name);
}


void reb_simulation_remove_all_particles(struct reb_simulation* const r){
    r->N 		= 0;
    r->N_allocated 	= 0;
    r->N_active 	= -1;
    r->N_var 	= 0;
    free(r->particles);
    r->particles 	= NULL;
}

int reb_simulation_remove_particle(struct reb_simulation* const r, int index, int keep_sorted){
    if (r->integrator == REB_INTEGRATOR_MERCURIUS){
        keep_sorted = 1; // Force keep_sorted for hybrid integrator
        struct reb_integrator_mercurius* rim = &(r->ri_mercurius);
        if (rim->N_allocated_dcrit>0 && index<(int)rim->N_allocated_dcrit){
            for (unsigned int i=0;i<r->N-1;i++){
                if ((int)i>=index){
                    rim->dcrit[i] = rim->dcrit[i+1];
                }
            }
        }
        reb_integrator_ias15_reset(r);
        if (r->ri_mercurius.mode==1){
            struct reb_integrator_mercurius* rim = &(r->ri_mercurius);
            int after_to_be_removed_particle = 0;
            int encounter_index = -1;
            for (unsigned int i=0;i<rim->encounter_N;i++){
                if (after_to_be_removed_particle == 1){
                    rim->encounter_map[i-1] = rim->encounter_map[i] - 1; 
                }
                if (rim->encounter_map[i]==index){
                    encounter_index = i;
                    after_to_be_removed_particle = 1;
                }
            }
            if (encounter_index<(int)rim->encounter_N_active){
                rim->encounter_N_active--;
            }
            rim->encounter_N--;
        }
    }

    if (r->integrator == REB_INTEGRATOR_TRACE){
        keep_sorted = 1; // Force keepSorted for hybrid integrator
        struct reb_integrator_trace* ri_trace = &(r->ri_trace);
        reb_integrator_bs_reset(r);
        if (r->ri_trace.mode==1 || r->ri_trace.mode==3){
            // Only removed mid-timestep if collision - BS Step!
            int after_to_be_removed_particle = 0;
            int encounter_index = -1;
            for (int i=0;i<ri_trace->encounter_N;i++){
                if (after_to_be_removed_particle == 1){
                    ri_trace->encounter_map[i-1] = ri_trace->encounter_map[i] - 1;
                }
                if (ri_trace->encounter_map[i]==index){
                    encounter_index = i;
                    after_to_be_removed_particle = 1;
                }
            }

            // reshuffle current_Ks
            unsigned int counter = 0;
            const int new_N = r->N-1;
            for (unsigned int i = 0; i < new_N; i++){
                if (i == index) counter += r->N;
                for (unsigned int j = 0; j < new_N; j++){
                    if (j == index) counter++;
                    ri_trace->current_Ks[i*new_N+j] = ri_trace->current_Ks[i*new_N+j+counter];
                }
            }
            if (encounter_index<ri_trace->encounter_N_active){
                ri_trace->encounter_N_active--;
            }
            ri_trace->encounter_N--;
        }
    }

    if (r->N==1){
        r->N = 0;
        if(r->free_particle_ap){
            r->free_particle_ap(&r->particles[index]);
        }
        reb_simulation_warning(r, "Last particle removed.");
        return 1;
    }
    if (index >= (int)r->N || index < 0){
        char warning[1024];
        sprintf(warning, "Index %d passed to particles_remove was out of range (N=%d).  Did not remove particle.", index, r->N);
        reb_simulation_error(r, warning);
        return 0;
    }
    if (r->N_var){
        reb_simulation_error(r, "Removing particles not supported when calculating MEGNO.  Did not remove particle.");
        return 0;
    }
    if(keep_sorted){
        r->N--;
        if(r->free_particle_ap){
            r->free_particle_ap(&r->particles[index]);
        }
        if(index<r->N_active){
            r->N_active--;
        }
        for(unsigned int j=index; j<r->N; j++){
            r->particles[j] = r->particles[j+1];
        }
        if (r->tree_root){
            reb_simulation_error(r, "REBOUND cannot remove a particle a tree and keep the particles sorted. Did not remove particle.");
            return 0;
        }
    }else{
        if (r->tree_root){
            // Just flag particle, will be removed in update_tree.
            r->particles[index].y = nan("");
            if(r->free_particle_ap){
                r->free_particle_ap(&r->particles[index]);
            }
        }else{
            r->N--;
            if(r->free_particle_ap){
                r->free_particle_ap(&r->particles[index]);
            }
            r->particles[index] = r->particles[r->N];
        }
    }

    return 1;
}

void reb_particle_isub(struct reb_particle* p1, struct reb_particle* p2){
    p1->x -= p2->x;
    p1->y -= p2->y;
    p1->z -= p2->z;
    p1->vx -= p2->vx;
    p1->vy -= p2->vy;
    p1->vz -= p2->vz;
    p1->m -= p2->m;
}

void reb_particle_iadd(struct reb_particle* p1, struct reb_particle* p2){
    p1->x += p2->x;
    p1->y += p2->y;
    p1->z += p2->z;
    p1->vx += p2->vx;
    p1->vy += p2->vy;
    p1->vz += p2->vz;
    p1->m += p2->m;
}

void reb_particle_imul(struct reb_particle* p1, double value){
    p1->x *= value;
    p1->y *= value;
    p1->z *= value;
    p1->vx *= value;
    p1->vy *= value;
    p1->vz *= value;
    p1->m *= value;
}

double reb_particle_distance(struct reb_particle* p1, struct reb_particle* p2){
    double dx = p1->x - p2->x;
    double dy = p1->y - p2->y;
    double dz = p1->z - p2->z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

struct reb_particle reb_particle_nan(void){
    struct reb_particle p = { 0 };
    p.x = nan("");
    p.y = nan("");
    p.z = nan("");
    p.vx = nan("");
    p.vy = nan("");
    p.vz = nan("");
    p.ax = nan("");
    p.ay = nan("");
    p.az = nan("");
    p.m = nan("");
    p.r = nan("");
    p.last_collision = nan("");
    return p;
}
