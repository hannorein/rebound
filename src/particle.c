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
#include "rebound.h"
#include "rebound_internal.h"
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
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

// Default size for hash table. Can be increased for very large simulations.
#define REB_NAME_HASH_TABLE_SIZE 1024



void reb_simulation_add(struct reb_simulation* const r, struct reb_particle pt){
    if (reb_boundary_particle_is_in_box(r, pt)==0){
        reb_simulation_error(r,"Particle outside of box boundaries. Did not add particle.");
        return;
    }
    // Allocate memory if needed.
    if (r->N_allocated<=r->N){
        unsigned int old_N_allocated = r->N_allocated;
        r->N_allocated = r->N_allocated ? r->N_allocated * 2 : 8;
        r->particles = realloc(r->particles,sizeof(struct reb_particle)*r->N_allocated);
        memset(r->particles + old_N_allocated, 0, (r->N_allocated - old_N_allocated) * sizeof(struct reb_particle));
    }

    r->particles[r->N] = pt;
    r->particles[r->N].sim = r;
    (r->N)++;
    // Particle was added successfully. Do other work now.
#ifndef MPI
    if (pt.name){
        reb_particle_set_name(&r->particles[r->N-1], pt.name);
    }
#endif // MPI

    // Check if any integrators need to do extra work
    if (r->integrator.did_add_particle){
        r->integrator.did_add_particle(r);
    }
}

// Compares two particles. Return 0 if identical.
int reb_particle_cmp(struct reb_particle p1, struct reb_particle p2){
    int differ = 0;
    differ = differ || (p1.x != p2.x);
    differ = differ || (p1.y != p2.y);
    differ = differ || (p1.z != p2.z);
    differ = differ || (p1.vx != p2.vx);
    differ = differ || (p1.vy != p2.vy);
    differ = differ || (p1.vz != p2.vz);
    differ = differ || (p1.ax != p2.ax);
    differ = differ || (p1.ay != p2.ay);
    differ = differ || (p1.az != p2.az);
    differ = differ || (p1.m != p2.m);
    differ = differ || (p1.r != p2.r);
    differ = differ || (p1.name != p2.name);
    return differ;
}


int reb_particle_check_testparticles(struct reb_simulation* const r){
    if (r->N_active == r->N || r->N_active == SIZE_MAX){
        return 0;
    }
    // Check if testparticle of type 0 has mass!=0
    if (r->testparticle_type == 0){
        int found_issue = 0;
        const int N = r->N;
#pragma omp parallel for
        for (int i=r->N_active; i<N; i++){
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


int reb_get_rootbox_for_particle(const struct reb_simulation* const r, struct reb_particle pt){
    if (r->root_size==-1) return 0;
    int i = ((int)floor((pt.x + r->root_size*(double)r->N_root_x/2.)/r->root_size)+r->N_root_x)%r->N_root_x;
    int j = ((int)floor((pt.y + r->root_size*(double)r->N_root_y/2.)/r->root_size)+r->N_root_y)%r->N_root_y;
    int k = ((int)floor((pt.z + r->root_size*(double)r->N_root_z/2.)/r->root_size)+r->N_root_z)%r->N_root_z;
    int index = (k*r->N_root_y+j)*r->N_root_x+i;
    return index;
}

int reb_simulation_particle_index(struct reb_particle* p){
    if (!p) return -1;
    struct reb_simulation* r = p->sim;
    if (!r) return -1;
    if (p<r->particles) return -1;
    size_t index = p-r->particles;
    if (index<r->N) return index;
    return -1; // Not found.
}
int reb_simulation_particle_var_index(struct reb_particle* p){
    if (!p) return -1;
    struct reb_simulation* r = p->sim;
    if (!r) return -1;
    if (p<r->particles_var) return -1;
    size_t index = p-r->particles_var;
    if (index<r->N_var) return index;
    return -1; // Not found.
}


static const char* get_registered_name(struct reb_simulation* r, const char* const name){
#ifdef MPI
    (void)name; // not used
    reb_simulation_error(r, "Particle names are not supported with MPI. Use integer ids instead.\n");
    return NULL; // Does not register.
#else // MPI
    if (name==NULL) return NULL; // NULL string not allowed.
    if (r->name_hash_table){
        uint32_t hash = reb_hash(name)%REB_NAME_HASH_TABLE_SIZE;
        struct reb_name_hash_item* item = &r->name_hash_table[hash];
        if (item->index>0){ // Entry exists
            do {            // Loop over linked list
                if (item->index < r->N+1){
                    struct reb_particle* p = &r->particles[item->index-1];
                    const char* p_name = p->name;
                    if (p_name){
                        if (strcmp(p_name,name)==0){
                            return p_name;
                        }
                    }
                }
                item = item->next;
            } while(item);
        }
    }
    // If not found yet. Go through entire name list.
    for (size_t i=0; i<r->N_name_list; i++){
        if (strcmp(name,r->name_list[i])==0){
            return r->name_list[i];
        }
    }
    return NULL; // Not found
#endif // MPI
}

const char* reb_simulation_register_name(struct reb_simulation* r, const char* const name){
    const char* registered_name = get_registered_name(r,name);
    if (registered_name) return registered_name;
    registered_name = strdup(name);
    r->N_name_list++;
    r->name_list = realloc(r->name_list,sizeof(char*)*r->N_name_list);
    r->name_list[r->N_name_list-1] = (char*)registered_name;
    return registered_name;
}

static void add_to_name_hash_table(struct reb_simulation* r, int index, const char* name){
    uint32_t hash = reb_hash(name)%REB_NAME_HASH_TABLE_SIZE;
    if (!r->name_hash_table){
        r->name_hash_table = calloc(REB_NAME_HASH_TABLE_SIZE,sizeof(struct reb_name_hash_item));
    }
    // Check for collision
    struct reb_name_hash_item* item = &r->name_hash_table[hash];
    if (!item->index){ //empty bucket
        item->index = index + 1;
    }else{
        // Find end of linked list
        while (item->next){
            item = item->next;
        }
        item->next = calloc(1, sizeof(struct reb_name_hash_item));
        item->next->index = index +1;
    }
}

struct reb_particle* reb_simulation_get_particle_by_name(struct reb_simulation* r, const char* const name){
    // Try hash table lookup first
    if (r->name_hash_table){
        int table_needs_repair = 0;
        uint32_t hash = reb_hash(name)%REB_NAME_HASH_TABLE_SIZE;
        struct reb_name_hash_item* item = &r->name_hash_table[hash];
        if (item->index>0){ // Entry exists
            do {            // Loop over linked list
                if (item->index > r->N){
                    table_needs_repair = 1; // Out of bounds. Particle got removed?
                }else{
                    struct reb_particle* p = &r->particles[item->index-1];
                    const char* p_name = p->name;
                    if (p_name){
                        if (strcmp(p_name,name)==0){
                            return p;
                        }
                    }else{
                        table_needs_repair = 1; // Expected a name.
                    }
                }
                item = item->next;
            } while(item);
        }
        if (table_needs_repair){
            reb_simulation_warning(r, "Name hash table needs repairs.");
        }
    }
    // If not found loop over all particles
    for (size_t i=0; i<r->N; i++){
        const char* p_name = r->particles[i].name;
        if (p_name){
            if (strcmp(p_name,name)==0){
                // Update hash table
                add_to_name_hash_table(r, i, name);
                return &(r->particles[i]);
            }
        }
    }
    return NULL; // Not found
}


#ifdef MPI
struct reb_particle reb_simulation_particle_by_id(struct reb_simulation* const r, size_t id){
    struct reb_particle* p = NULL;
    for (size_t i=0; i<r->N; i++){
        size_t p_id = (size_t)(r->particles[i].name);
        if (p_id == id){
            p = &(r->particles[i]);
            break;
        }
    }

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

int reb_simulation_remove_particle_by_name(struct reb_simulation* r, const char* const name){
    struct reb_particle* p = reb_simulation_get_particle_by_name(r, name);
    if (!p){
        reb_simulation_error(r, "Particle not found.");
        return 1; // Not found.
    }
    size_t index = p - r->particles;
    return reb_simulation_remove_particle(r, index);
}


void reb_particle_set_name(struct reb_particle* p, const char* const name){
    if (name==NULL){
        // Delete name.
        p->name = NULL;
        return;
    }
    struct reb_simulation* r = p->sim;
#ifndef MPI
    if (!r){
        reb_simulation_error(NULL,"Cannot set particle name using_reb_particle_set_name() as the particle is not part of a simulation. You can set the name manually.");
        return;
    }
    p->name = reb_simulation_register_name(r,name);

    uint32_t index = p - r->particles;
    add_to_name_hash_table(r, index, name);
#else //  MPI
    (void)name; // not used
    reb_simulation_error(r, "Particle names are not supported with MPI. Use integer ids instead.\n");
    __asm__ ("int3");
#endif //  MPI
}


void reb_simulation_remove_all_particles(struct reb_simulation* const r){
    r->N = 0;
    r->N_allocated = 0;
    r->N_active = SIZE_MAX;
    r->N_var = 0;
    free(r->particles);
    r->particles = NULL;
    free(r->particles_var);
    r->particles_var = NULL;
}

int reb_simulation_remove_particle(struct reb_simulation* const r, size_t index){
    if (r->N_var){
        reb_simulation_error(r, "Removing particles not supported when variational particles are in use. Did not remove particle.");
        return 1;
    }
    // Check if any integrators need to do work before removing particle
    if (r->integrator.will_remove_particle){
        r->integrator.will_remove_particle(r, index);
    }

    if (r->N==1){
        r->N = 0;
        if(r->free_particle_ap){
            r->free_particle_ap(&r->particles[index]);
        }
#ifndef MPI // There might be empty nodes when MPI is used.
        reb_simulation_warning(r, "Last particle removed.");
#endif // MPI
        return 0;
    }
    if (index >= r->N){
        char warning[1024];
        sprintf(warning, "Index %zu passed to particles_remove was out of range (N=%zu).  Did not remove particle.", index, r->N);
        reb_simulation_error(r, warning);
        return 1;
    }
    r->N--;
    if(r->free_particle_ap){
        r->free_particle_ap(&r->particles[index]);
    }
    if(index<r->N_active && r->N_active!=SIZE_MAX){
        r->N_active--;
    }
    for(size_t j=index; j<r->N; j++){
        r->particles[j] = r->particles[j+1];
    }

    return 0; // Success
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
    return p;
}
