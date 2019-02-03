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
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include "rebound.h"
#include "tree.h"
#include "boundary.h"
#include "particle.h"
#include "integrator_ias15.h"
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

static void reb_add_local(struct reb_simulation* const r, struct reb_particle pt){
	if (reb_boundary_particle_is_in_box(r, pt)==0){
		// reb_particle has left the box. Do not add.
		reb_error(r,"Particle outside of box boundaries. Did not add particle.");
		return;
	}
	while (r->allocatedN<=r->N){
		r->allocatedN = r->allocatedN ? r->allocatedN * 2 : 128;
		r->particles = realloc(r->particles,sizeof(struct reb_particle)*r->allocatedN);
	}

	r->particles[r->N] = pt;
	r->particles[r->N].sim = r;
	if (r->gravity==REB_GRAVITY_TREE || r->collision==REB_COLLISION_TREE){
        if (r->root_size==-1){
            reb_error(r,"root_size is -1. Make sure you call reb_configure_box() before using a tree based gravity or collision solver.");
            return;
        }
		reb_tree_add_particle_to_tree(r, r->N);
	}
	(r->N)++;
}

void reb_add(struct reb_simulation* const r, struct reb_particle pt){
#ifndef COLLISIONS_NONE
	if (pt.r>=r->max_radius[0]){
		r->max_radius[1] = r->max_radius[0];
		r->max_radius[0] = pt.r;
	}else{
		if (pt.r>=r->max_radius[1]){
			r->max_radius[1] = pt.r;
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
	int root_n_per_node = r->root_n/r->mpi_num;
	int proc_id = rootbox/root_n_per_node;
	if (proc_id != r->mpi_id && r->N >= r->N_active){
		// Add particle to array and send them to proc_id later. 
		reb_communication_mpi_add_particle_to_send_queue(r,pt,proc_id);
		return;
	}
#endif // MPI
	// Add particle to local partical array.
	reb_add_local(r, pt);
}

int reb_get_rootbox_for_particle(const struct reb_simulation* const r, struct reb_particle pt){
	if (r->root_size==-1) return 0;
	int i = ((int)floor((pt.x + r->boxsize.x/2.)/r->root_size)+r->root_nx)%r->root_nx;
	int j = ((int)floor((pt.y + r->boxsize.y/2.)/r->root_size)+r->root_ny)%r->root_ny;
	int k = ((int)floor((pt.z + r->boxsize.z/2.)/r->root_size)+r->root_nz)%r->root_nz;
	int index = (k*r->root_ny+j)*r->root_nx+i;
	return index;
}

int reb_get_particle_index(struct reb_particle* p){
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
            if(lookup[middle].index < r->N){
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
    for(int i=0; i<r->N; i++){
        if(N_hash >= r->allocatedN_lookup){
            r->allocatedN_lookup = r->allocatedN_lookup ? r->allocatedN_lookup * 2 : 128;
            r->particle_lookup_table = realloc(r->particle_lookup_table, sizeof(struct reb_hash_pointer_pair)*r->allocatedN_lookup);
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

struct reb_particle* reb_get_particle_by_hash(struct reb_simulation* const r, uint32_t hash){
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

void reb_remove_all(struct reb_simulation* const r){
	r->N 		= 0;
	r->allocatedN 	= 0;
	r->N_active 	= -1;
	r->N_var 	= 0;
	free(r->particles);
	r->particles 	= NULL;
}

int reb_remove(struct reb_simulation* const r, int index, int keepSorted){
    if (r->integrator == REB_INTEGRATOR_MERCURIUS){
        keepSorted = 1; // Force keepSorted for hybrid integrator
        struct reb_simulation_integrator_mercurius* rim = &(r->ri_mercurius);
        for (int i=0;i<r->N-1;i++){
            if (i>=index){
                rim->dcrit[i] = rim->dcrit[i+1];
            }
        }
        reb_integrator_ias15_reset(r);
        if (r->ri_mercurius.mode==1){
            struct reb_simulation_integrator_mercurius* rim = &(r->ri_mercurius);
            int after_to_be_removed_particle = 0;
            for (int i=0;i<rim->encounterN;i++){
                if (after_to_be_removed_particle == 1){
                    rim->encounter_map[i-1] = rim->encounter_map[i] - 1; 
                }
                if (rim->encounter_map[i]==index){
                    after_to_be_removed_particle = 1;
                }
            }
            if (index<rim->encounterNactive){
                rim->encounterNactive--;
            }
            rim->encounterN--;
        }
    }
	if (r->N==1){
	    r->N = 0;
        if(r->free_particle_ap){
            r->free_particle_ap(&r->particles[index]);
        }
		reb_warning(r, "Last particle removed.");
		return 1;
	}
	if (index >= r->N){
		char warning[1024];
        sprintf(warning, "Index %d passed to particles_remove was out of range (N=%d).  Did not remove particle.", index, r->N);
		reb_error(r, warning);
		return 0;
	}
	if (r->N_var){
		reb_error(r, "Removing particles not supported when calculating MEGNO.  Did not remove particle.");
		return 0;
	}
	if(keepSorted){
	    r->N--;
        if(r->free_particle_ap){
            r->free_particle_ap(&r->particles[index]);
        }
        if(index<r->N_active){
            r->N_active--;
        }
		for(int j=index; j<r->N; j++){
			r->particles[j] = r->particles[j+1];
		}
        if (r->tree_root){
		    reb_error(r, "REBOUND cannot remove a particle a tree and keep the particles sorted. Did not remove particle.");
		    return 0;
        }
	}else{
        if (r->tree_root){
            // Just flag particle, will be removed in tree_update.
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

int reb_remove_by_hash(struct reb_simulation* const r, uint32_t hash, int keepSorted){
    struct reb_particle* p = reb_get_particle_by_hash(r, hash);
    if(p == NULL){
		reb_error(r,"Particle to be removed not found in simulation.  Did not remove particle.");
        return 0;
    }
    else{
        int index = reb_get_particle_index(p);
        return reb_remove(r, index, keepSorted);
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
    p.lastcollision = nan("");
    p.c = NULL;
    p.hash = 0;
    p.ap = NULL;
    p.sim = NULL;

    return p;
}
