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
#ifndef COLLISIONS_NONE
#include "collision.h"
#endif // COLLISIONS_NONE
#ifdef MPI
#include "communication_mpi.h"
#endif // MPI

#ifdef GRAVITY_GRAPE
#warning Fix this. 
extern double gravity_minimum_mass;
#endif // GRAVITY_GRAPE

static void reb_simulation_add_local(struct reb_simulation* const r, struct reb_particle pt){
	if (reb_boundary_particle_is_in_box(r, pt)==0){
		// reb_particle has left the box. Do not add.
		reb_simulation_error(r,"Particle outside of box boundaries. Did not add particle.");
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
#ifndef COLLISIONS_NONE
	if (pt.r>=r->max_radius0){
		r->max_radius1 = r->max_radius0;
		r->max_radius0 = pt.r;
	}else{
		if (pt.r>=r->max_radius1){
			r->max_radius1 = pt.r;
		}
	}
#endif 	// COLLISIONS_NONE
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
	// Add particle to local partical array.
	reb_simulation_add_local(r, pt);
}

int reb_particle_check_testparticles(struct reb_simulation* const r){
    if (r->N_active == (int)r->N || r->N_active == -1){
        return 0;
    }
    // Check if testparticle of type 0 has mass!=0
    if (r->testparticle_type == 0){
        const int N_real = r->N - r->N_var;
        for (int i=r->N_active; i<N_real; i++){
            if (r->particles[i].m!=0.){
                return 1;
            }
        }
    }
    return 0;
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

static struct reb_particle* reb_search_lookup_table(struct reb_simulation* const r, uint32_t hash){
    const struct reb_hash_pointer_pair* const lookup = r->particle_lookup_table;
    if (lookup == NULL){
        return NULL;
    }

    int left = 0;
    int right = r->N_lookup-1;
    int middle;
    while(left <= right){
        middle = (left + right)/2;
        uint32_t lookuphash = lookup[middle].hash;
        if(lookuphash < hash){
            left = middle+1;
        }
        else if(lookuphash > hash){
            right = middle-1;
        }
        else if(lookuphash == hash){
            if(lookup[middle].index < (int)r->N){
                return &r->particles[lookup[middle].index];
            }
            else{ // found lookup table entry pointing beyond r->N in particles array. Needs update
                return NULL;
            }
        }
    }
    return NULL;
}

static int compare_hash(const void* a, const void* b){
    struct reb_hash_pointer_pair* ia = (struct reb_hash_pointer_pair*)a; 
    struct reb_hash_pointer_pair* ib = (struct reb_hash_pointer_pair*)b;
    return (ia->hash > ib->hash) - (ia->hash < ib->hash); // to avoid overflow possibilities
}

static void reb_update_particle_lookup_table(struct reb_simulation* const r){
    const struct reb_particle* const particles = r->particles;
    int N_hash = 0;
    int zerohash = -1;
    for(unsigned int i=0; i<r->N; i++){
        if(N_hash >= r->N_allocated_lookup){
            r->N_allocated_lookup = r->N_allocated_lookup ? r->N_allocated_lookup * 2 : 128;
            r->particle_lookup_table = realloc(r->particle_lookup_table, sizeof(struct reb_hash_pointer_pair)*r->N_allocated_lookup);
        }
        if(particles[i].hash == 0){ // default hash (0) special case
            if (zerohash == -1){    // first zero hash
                zerohash = i;
                r->particle_lookup_table[zerohash].hash = particles[i].hash;
                r->particle_lookup_table[zerohash].index = i;
                N_hash++;
            }
            else{                   // update zero hash entry in lookup without incrementing N_hash
                r->particle_lookup_table[zerohash].index = i;
            }
        }
        else{                   
            r->particle_lookup_table[N_hash].hash = particles[i].hash;
            r->particle_lookup_table[N_hash].index = i;
            N_hash++;
        }
    }
    r->N_lookup = N_hash;
    qsort(r->particle_lookup_table, r->N_lookup, sizeof(*r->particle_lookup_table), compare_hash); // only sort the first N_lookup entries that are initialized.
}

struct reb_particle* reb_simulation_particle_by_hash(struct reb_simulation* const r, uint32_t hash){
    struct reb_particle* p; 
    p = reb_search_lookup_table(r, hash);
    if (p == NULL){
        reb_update_particle_lookup_table(r);
        p = reb_search_lookup_table(r, hash);
    }
    else{
        if (p->hash != hash){
            reb_update_particle_lookup_table(r);
            p = reb_search_lookup_table(r, hash);
        }
    }
    return p;
}

struct reb_particle reb_simulation_particle_by_hash_mpi(struct reb_simulation* const r, uint32_t hash){
#ifdef MPI
    struct reb_particle* p = reb_simulation_particle_by_hash(r, hash);
    int found = p==0?0:1;
    int found_sum = 0;
    MPI_Allreduce(&found, &found_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (found_sum == 0){
        return reb_particle_nan();
    }
    if (found_sum > 1){
        reb_simulation_error(r, "Multiple particles with same hash found.");
    }
    int root;
    int mayberoot = found ? r->mpi_id : 0;
    MPI_Allreduce(&mayberoot, &root, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    struct reb_particle ph = {0};
    if (found){
        ph = *p;
        ph.sim = NULL;
    }
    MPI_Bcast(&ph, sizeof(struct reb_particle), MPI_CHAR, root, MPI_COMM_WORLD);
    return ph;
#else // MPI
    struct reb_particle* p = reb_simulation_particle_by_hash(r, hash);
    if (p==0){
        return reb_particle_nan();
    }else{
       return *p;
    }
#endif // MPI
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

int reb_simulation_remove_particle_by_hash(struct reb_simulation* const r, uint32_t hash, int keep_sorted){
    struct reb_particle* p = reb_simulation_particle_by_hash(r, hash);
    if(p == NULL){
		reb_simulation_error(r,"Particle to be removed not found in simulation.  Did not remove particle.");
        return 0;
    }
    else{
        int index = reb_simulation_particle_index(p);
        return reb_simulation_remove_particle(r, index, keep_sorted);
    }
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
    struct reb_particle p;
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
    p.c = NULL;
    p.hash = 0;
    p.ap = NULL;
    p.sim = NULL;

    return p;
}
