/**
 * @file        collision.c
 * @brief       Collision search routine.
 * @author      Hanno Rein <hanno@hanno-rein.de>
 *
 * @details     A collision is defined as an overlap between two particles. This
 * is only an approximation and works only if the timestep is small
 * enough. More precisely, dt << v / Rp, where v is the typical velocity
 * and Rp the radius of a particle. Furthermore, particles must be 
 * approaching each other at the time when they overlap. 
 * 
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
#include "particle.h"
#include "collision.h"
#include "rebound.h"
#include "boundary.h"
#include "tree.h"
#ifdef MPI
#include "communication_mpi.h"
#endif // MPI
#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< Returns the minimum of a and b
#define MAX(a, b) ((a) > (b) ? (a) : (b))    ///< Returns the maximum of a and b

static void reb_tree_get_nearest_neighbour_in_cell(struct reb_simulation* const r, int* collisions_N, struct reb_ghostbox gb, struct reb_ghostbox gbunmod, int ri, double p1_r,  double* nearest_r2, struct reb_collision* collision_nearest, struct reb_treecell* c);
static void reb_tree_check_for_overlapping_trajectories_in_cell(struct reb_simulation* const r, int* collisions_N, struct reb_ghostbox gb, struct reb_ghostbox gbunmod, int ri, double p1_r, double p1_r_plus_dtv, struct reb_collision* collision_nearest, struct reb_treecell* c, double maxdrift);

void reb_collision_search(struct reb_simulation* const r){
    int N = r->N - r->N_var;
    int Ninner = N;
    int* mercurius_map = NULL;
    if (r->integrator==REB_INTEGRATOR_MERCURIUS){
        if (r->ri_mercurius.mode==0){
            // After jump step, only collisions with star might occur.
            // All other collisions in encounter step/
            Ninner = 1;
        }else{
            N = r->ri_mercurius.encounterN;
            Ninner = N;
            mercurius_map = r->ri_mercurius.encounter_map;
        }
    }
    int collisions_N = 0;
    const struct reb_particle* const particles = r->particles;
    switch (r->collision){
        case REB_COLLISION_NONE:
        break;
        case REB_COLLISION_DIRECT:
        {
            // Loop over ghost boxes, but only the inner most ring.
            int nghostxcol = (r->nghostx>1?1:r->nghostx);
            int nghostycol = (r->nghosty>1?1:r->nghosty);
            int nghostzcol = (r->nghostz>1?1:r->nghostz);
            for (int gbx=-nghostxcol; gbx<=nghostxcol; gbx++){
            for (int gby=-nghostycol; gby<=nghostycol; gby++){
            for (int gbz=-nghostzcol; gbz<=nghostzcol; gbz++){
                // Loop over all particles
                for (int i=0;i<N;i++){
#ifndef OPENMP
                    if (reb_sigint) return;
#endif // OPENMP
                    int ip = i;
                    if (mercurius_map){
                        ip = mercurius_map[i];
                    }
                    struct reb_particle p1 = particles[ip];
                    struct reb_ghostbox gborig = reb_boundary_get_ghostbox(r, gbx,gby,gbz);
                    struct reb_ghostbox gb = gborig;
                    // Precalculate shifted position 
                    gb.shiftx += p1.x;
                    gb.shifty += p1.y;
                    gb.shiftz += p1.z;
                    gb.shiftvx += p1.vx;
                    gb.shiftvy += p1.vy;
                    gb.shiftvz += p1.vz;
                    // Loop over all particles again
                    for (int j=0;j<Ninner;j++){
                        // Do not collide particle with itself.
                        if (i==j) continue;
                        int jp = j;
                        if (mercurius_map){
                            jp = mercurius_map[j];
                        }
                        struct reb_particle p2 = particles[jp];
                        double dx = gb.shiftx - p2.x; 
                        double dy = gb.shifty - p2.y; 
                        double dz = gb.shiftz - p2.z; 
                        double sr = p1.r + p2.r; 
                        double r2 = dx*dx+dy*dy+dz*dz;
                        // Check if particles are overlapping 
                        if (r2>sr*sr) continue;    
                        double dvx = gb.shiftvx - p2.vx; 
                        double dvy = gb.shiftvy - p2.vy; 
                        double dvz = gb.shiftvz - p2.vz; 
                        // Check if particles are approaching each other
                        if (dvx*dx + dvy*dy + dvz*dz >0) continue; 
                        // Add particles to collision array.
                        if (r->collisions_allocatedN<=collisions_N){
                            // Allocate memory if there is no space in array.
                            // Init to 32 if no space has been allocated yet, otherwise double it.
                            r->collisions_allocatedN = r->collisions_allocatedN ? r->collisions_allocatedN * 2 : 32;
                            r->collisions = realloc(r->collisions,sizeof(struct reb_collision)*r->collisions_allocatedN);
                        }
                        r->collisions[collisions_N].p1 = ip;
                        r->collisions[collisions_N].p2 = jp;
                        r->collisions[collisions_N].gb = gborig;
                        collisions_N++;
                    }
                }
            }
            }
            }
        }
        break;
        case REB_COLLISION_LINE:
        {
            double dt_last_done = r->dt_last_done;
            // Loop over ghost boxes, but only the inner most ring.
            int nghostxcol = (r->nghostx>1?1:r->nghostx);
            int nghostycol = (r->nghosty>1?1:r->nghosty);
            int nghostzcol = (r->nghostz>1?1:r->nghostz);
            for (int gbx=-nghostxcol; gbx<=nghostxcol; gbx++){
            for (int gby=-nghostycol; gby<=nghostycol; gby++){
            for (int gbz=-nghostzcol; gbz<=nghostzcol; gbz++){
                // Loop over all particles
                for (int i=0;i<N;i++){
#ifndef OPENMP
                    if (reb_sigint) return;
#endif // OPENMP
                    struct reb_particle p1 = particles[i];
                    struct reb_ghostbox gborig = reb_boundary_get_ghostbox(r, gbx,gby,gbz);
                    struct reb_ghostbox gb = gborig;
                    // Precalculate shifted position 
                    gb.shiftx += p1.x;
                    gb.shifty += p1.y;
                    gb.shiftz += p1.z;
                    gb.shiftvx += p1.vx;
                    gb.shiftvy += p1.vy;
                    gb.shiftvz += p1.vz;
                    // Loop over all particles again
                    for (int j=i+1;j<N;j++){
                        struct reb_particle p2 = particles[j];
                        const double dx1 = gb.shiftx - p2.x; // distance at end
                        const double dy1 = gb.shifty - p2.y;
                        const double dz1 = gb.shiftz - p2.z;
                        const double r1 = (dx1*dx1 + dy1*dy1 + dz1*dz1);
                        const double dvx1 = gb.shiftvx - p2.vx; 
                        const double dvy1 = gb.shiftvy - p2.vy;
                        const double dvz1 = gb.shiftvz - p2.vz;
                        const double dx2 = dx1 -dt_last_done*dvx1; // distance at beginning
                        const double dy2 = dy1 -dt_last_done*dvy1;
                        const double dz2 = dz1 -dt_last_done*dvz1;
                        const double r2 = (dx2*dx2 + dy2*dy2 + dz2*dz2);
                        const double t_closest = (dx1*dvx1 + dy1*dvy1 + dz1*dvz1)/(dvx1*dvx1 + dvy1*dvy1 + dvz1*dvz1);

                        double rmin2_ab = MIN(r1,r2);
                        if (t_closest/dt_last_done>=0. && t_closest/dt_last_done<=1.){
                            const double dx3 = dx1-t_closest*dvx1; // closest approach
                            const double dy3 = dy1-t_closest*dvy1;
                            const double dz3 = dz1-t_closest*dvz1;
                            const double r3 = (dx3*dx3 + dy3*dy3 + dz3*dz3);
                            rmin2_ab = MIN(rmin2_ab, r3);
                        }
                        double rsum = p1.r + p2.r;
                        if (rmin2_ab>rsum*rsum) continue;

                        // Add particles to collision array.
                        if (r->collisions_allocatedN<=collisions_N){
                            // Allocate memory if there is no space in array.
                            // Init to 32 if no space has been allocated yet, otherwise double it.
                            r->collisions_allocatedN = r->collisions_allocatedN ? r->collisions_allocatedN * 2 : 32;
                            r->collisions = realloc(r->collisions,sizeof(struct reb_collision)*r->collisions_allocatedN);
                        }
                        r->collisions[collisions_N].p1 = i;
                        r->collisions[collisions_N].p2 = j;
                        r->collisions[collisions_N].gb = gborig;
                        collisions_N++;
                    }
                }
            }
            }
            }
        }
        break;
        case REB_COLLISION_TREE:
        {
            // Update and simplify tree. 
            // Prepare particles for distribution to other nodes. 
            reb_tree_update(r);          

#ifdef MPI
            // Distribute particles and add newly received particles to tree.
            reb_communication_mpi_distribute_particles(r);
            
            // Prepare essential tree (and particles close to the boundary needed for collisions) for distribution to other nodes.
            reb_tree_prepare_essential_tree_for_collisions(r);

            // Transfer essential tree and particles needed for collisions.
            reb_communication_mpi_distribute_essential_tree_for_collisions(r);
#endif // MPI

            // Loop over ghost boxes, but only the inner most ring.
            int nghostxcol = (r->nghostx>1?1:r->nghostx);
            int nghostycol = (r->nghosty>1?1:r->nghosty);
            int nghostzcol = (r->nghostz>1?1:r->nghostz);
            const struct reb_particle* const particles = r->particles;
            const int N = r->N - r->N_var;
            // Loop over all particles
#pragma omp parallel for schedule(guided)
            for (int i=0;i<N;i++){
#ifndef OPENMP
                if (reb_sigint) return;
#endif // OPENMP
                struct reb_particle p1 = particles[i];
                struct reb_collision collision_nearest;
                collision_nearest.p1 = i;
                collision_nearest.p2 = -1;
                double p1_r = p1.r;
                double nearest_r2 = r->boxsize_max*r->boxsize_max/4.;
                // Loop over ghost boxes.
                for (int gbx=-nghostxcol; gbx<=nghostxcol; gbx++){
                for (int gby=-nghostycol; gby<=nghostycol; gby++){
                for (int gbz=-nghostzcol; gbz<=nghostzcol; gbz++){
                    // Calculated shifted position (for speedup). 
                    struct reb_ghostbox gb = reb_boundary_get_ghostbox(r, gbx,gby,gbz);
                    struct reb_ghostbox gbunmod = gb;
                    gb.shiftx += p1.x; 
                    gb.shifty += p1.y; 
                    gb.shiftz += p1.z; 
                    gb.shiftvx += p1.vx; 
                    gb.shiftvy += p1.vy; 
                    gb.shiftvz += p1.vz; 
                    // Loop over all root boxes.
                    for (int ri=0;ri<r->root_n;ri++){
                        struct reb_treecell* rootcell = r->tree_root[ri];
                        if (rootcell!=NULL){
                            reb_tree_get_nearest_neighbour_in_cell(r, &collisions_N, gb, gbunmod,ri,p1_r,&nearest_r2,&collision_nearest,rootcell);
                        }
                    }
                }
                }
                }
                // Continue if no collision was found
                if (collision_nearest.p2==-1) continue;
            }
        }
        break;
        case REB_COLLISION_LINETREE:
        {
            // Calculate max drift (can also be stored in tree for further speedup)
            double vmax2 = 0.;
            for (int i=0;i<N;i++){
                struct reb_particle p1 = particles[i];
                vmax2 = MAX(vmax2, p1.vx*p1.vx + p1.vy*p1.vy + p1.vz*p1.vz);
            }
            double maxdrift = r->dt_last_done*sqrt(vmax2);
            // Update and simplify tree. 
            // Prepare particles for distribution to other nodes. 
            reb_tree_update(r);          

            // Loop over ghost boxes, but only the inner most ring.
            int nghostxcol = (r->nghostx>1?1:r->nghostx);
            int nghostycol = (r->nghosty>1?1:r->nghosty);
            int nghostzcol = (r->nghostz>1?1:r->nghostz);
            const struct reb_particle* const particles = r->particles;
            const int N = r->N - r->N_var;
            // Loop over all particles
#pragma omp parallel for schedule(guided)
            for (int i=0;i<N;i++){
#ifndef OPENMP
                if (reb_sigint) return;
#endif // OPENMP
                struct reb_particle p1 = particles[i];
                struct reb_collision collision_nearest;
                collision_nearest.p1 = i;
                collision_nearest.p2 = -1;
                double p1_r = p1.r;
                // Add drift during last timestep
                double p1_r_plus_dtv = p1_r + r->dt_last_done*sqrt(p1.vx*p1.vx + p1.vy*p1.vy + p1.vz*p1.vz);
                // Loop over ghost boxes.
                for (int gbx=-nghostxcol; gbx<=nghostxcol; gbx++){
                for (int gby=-nghostycol; gby<=nghostycol; gby++){
                for (int gbz=-nghostzcol; gbz<=nghostzcol; gbz++){
                    // Calculated shifted position (for speedup). 
                    struct reb_ghostbox gb = reb_boundary_get_ghostbox(r, gbx,gby,gbz);
                    struct reb_ghostbox gbunmod = gb;
                    gb.shiftx += p1.x; 
                    gb.shifty += p1.y; 
                    gb.shiftz += p1.z; 
                    gb.shiftvx += p1.vx; 
                    gb.shiftvy += p1.vy; 
                    gb.shiftvz += p1.vz; 
                    // Loop over all root boxes.
                    for (int ri=0;ri<r->root_n;ri++){
                        struct reb_treecell* rootcell = r->tree_root[ri];
                        if (rootcell!=NULL){
                            reb_tree_check_for_overlapping_trajectories_in_cell(r, &collisions_N, gb, gbunmod,ri,p1_r,p1_r_plus_dtv,&collision_nearest,rootcell,maxdrift);
                        }
                    }
                }
                }
                }
                // Continue if no collision was found
                if (collision_nearest.p2==-1) continue;
            }
        }
        break;
        default:
            reb_exit("Collision routine not implemented.");
    }

    // randomize
    for (int i=0;i<collisions_N;i++){
        int new = rand_r(&(r->rand_seed))%collisions_N;
        struct reb_collision c1 = r->collisions[i];
        r->collisions[i] = r->collisions[new];
        r->collisions[new] = c1;
    }
    // Loop over all collisions previously found in reb_collision_search().
    
    int (*resolve) (struct reb_simulation* const r, struct reb_collision c) = r->collision_resolve;
    if (resolve==NULL){
        // Default is to throw an exception
        resolve = reb_collision_resolve_halt;
    }
    unsigned int collision_resolve_keep_sorted = r->collision_resolve_keep_sorted;
    if (r->integrator == REB_INTEGRATOR_MERCURIUS){
        collision_resolve_keep_sorted = 1; // Force keep_sorted for hybrid integrator
    }

    for (int i=0;i<collisions_N;i++){
        
        struct reb_collision c = r->collisions[i];
        if (c.p1 != -1 && c.p2 != -1){
            // Resolve collision
            int outcome = resolve(r, c);
            
            // Remove particles
            if (outcome & 1){
                // Remove p1
                int removedp1 = reb_remove(r,c.p1,collision_resolve_keep_sorted);
                if (removedp1){
                    if (r->tree_root){ // In a tree, particles get removed later. 
                        for (int j=i+1;j<collisions_N;j++){ // Update other collisions
                            struct reb_collision* cp = &(r->collisions[j]);
                            // Skip collisions which involved the removed particle
                            if (cp->p1==c.p1 || cp->p2==c.p1){
                                cp->p1 = -1;
                                cp->p2 = -1;
                            }
                        }
                    }else{ // Not in a tree, particles get removed immediately 
                        // Update p2 of current collision
                        if (collision_resolve_keep_sorted){
                            if (c.p2 > c.p1){
                                c.p2--;
                            }
                        }else{
                            if (c.p2 == r->N-r->N_var){
                                c.p2 = c.p1;
                            }
                        }
                        for (int j=i+1;j<collisions_N;j++){ // Update other collisions
                            struct reb_collision* cp = &(r->collisions[j]);
                            // Skip collisions which involve the removed particle
                            if (cp->p1==c.p1 || cp->p2==c.p1){
                                cp->p1 = -1;
                                cp->p2 = -1;
                            }
                            // Adjust collisions
                            if (collision_resolve_keep_sorted){
                                if (cp->p1 > c.p1){
                                    cp->p1--;
                                }
                                if (cp->p2 > c.p1){
                                    cp->p2--;
                                }
                            }else{
                                if (cp->p1 == r->N-r->N_var){
                                    cp->p1 = c.p1;
                                }
                                if (cp->p2 == r->N-r->N_var){
                                    cp->p2 = c.p1;
                                }
                            }
                        }
                    }
                }
            }
            if (outcome & 2){
                // Remove p1
                int removedp2 = reb_remove(r,c.p2,collision_resolve_keep_sorted);
                if (removedp2){ // Update other collisions
                    if (r->tree_root){ // In a tree, particles get removed later. 
                        for (int j=i+1;j<collisions_N;j++){ // Update other collisions
                            struct reb_collision* cp = &(r->collisions[j]);
                            // Skip collisions which involved the removed particle
                            if (cp->p1==c.p2 || cp->p2==c.p2){
                                cp->p1 = -1;
                                cp->p2 = -1;
                            }
                        }
                    }else{ // Not in a tree, particles get removed immediately 
                        for (int j=i+1;j<collisions_N;j++){
                            struct reb_collision* cp = &(r->collisions[j]);
                            // Skip collisions which involve the removed particle
                            if (cp->p1==c.p2 || cp->p2==c.p2){
                                cp->p1 = -1;
                                cp->p2 = -1;
                            }
                            // Adjust collisions
                            if (collision_resolve_keep_sorted){
                                if (cp->p1 > c.p2){
                                    cp->p1--;
                                }
                                if (cp->p2 > c.p2){
                                    cp->p2--;
                                }
                            }else{
                                if (cp->p1 == r->N-r->N_var){
                                    cp->p1 = c.p2;
                                }
                                if (cp->p2 == r->N-r->N_var){
                                    cp->p2 = c.p2;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/**
 * @brief Workaround for python setters.
 **/
void reb_set_collision_resolve(struct reb_simulation* r, int (*resolve) (struct reb_simulation* const r, struct reb_collision c)){
    r->collision_resolve = resolve;
}

/**
 * @brief Find the nearest neighbour in a cell or its daughters.
 * @details The function only returns a positive result if the particles
 * are overlapping. Thus, the name nearest neighbour is not
 * exactly true.
 * @param r REBOUND simulation to work on.
 * @param gb (Shifted) position and velocity of the particle.
 * @param ri Index of the root box currently being searched in.
 * @param p1_r Radius of the particle (this is not in gb).
 * @param nearest_r2 Pointer to the nearest neighbour found so far.
 * @param collision_nearest Pointer to the nearest collision found so far.
 * @param c Pointer to the cell currently being searched in.
 * @param collisions_N Pointer to current number of collisions
 * @param gbunmod Ghostbox unmodified
 */
static void reb_tree_get_nearest_neighbour_in_cell(struct reb_simulation* const r, int* collisions_N, struct reb_ghostbox gb, struct reb_ghostbox gbunmod, int ri, double p1_r, double* nearest_r2, struct reb_collision* collision_nearest, struct reb_treecell* c){
    const struct reb_particle* const particles = r->particles;
    if (c->pt>=0){     
        // c is a leaf node
        int condition     = 1;
#ifdef MPI
        int isloc    = 1 ;
        isloc = reb_communication_mpi_rootbox_is_local(r, ri);
        if (isloc==1){
#endif // MPI
            /**
             * If this is a local cell, make sure particle is not colliding with itself.
             * If this is a remote cell, the particle number might be the same, even for 
             * different particles. 
             * TODO: This can probably be written in a cleaner way.
             */
            condition = (c->pt != collision_nearest->p1);
#ifdef MPI
        }
#endif // MPI
        if (condition){
            struct reb_particle p2;
#ifdef MPI
            if (isloc==1){
#endif // MPI
                p2 = particles[c->pt];
#ifdef MPI
            }else{
                int root_n_per_node = r->root_n/r->mpi_num;
                int proc_id = ri/root_n_per_node;
                p2 = r->particles_recv[proc_id][c->pt];
            }
#endif // MPI

            double dx = gb.shiftx - p2.x;
            double dy = gb.shifty - p2.y;
            double dz = gb.shiftz - p2.z;
            double r2 = dx*dx+dy*dy+dz*dz;
            // A closer neighbour has already been found 
            //if (r2 > *nearest_r2) return;
            double rp = p1_r+p2.r;
            // reb_particles are not overlapping 
            if (r2 > rp*rp) return;
            double dvx = gb.shiftvx - p2.vx;
            double dvy = gb.shiftvy - p2.vy;
            double dvz = gb.shiftvz - p2.vz;
            // reb_particles are not approaching each other
            if (dvx*dx + dvy*dy + dvz*dz >0) return;
            // Found a new nearest neighbour. Save it for later.
            *nearest_r2 = r2;
            collision_nearest->ri = ri;
            collision_nearest->p2 = c->pt;
            collision_nearest->gb = gbunmod;
            // Save collision in collisions array.
#pragma omp critical
            {
                if (r->collisions_allocatedN<=(*collisions_N)){
                    // Init to 32 if no space has been allocated yet, otherwise double it.
                    r->collisions_allocatedN = r->collisions_allocatedN ? r->collisions_allocatedN * 2 : 32;
                    r->collisions = realloc(r->collisions,sizeof(struct reb_collision)*r->collisions_allocatedN);
                }
                r->collisions[(*collisions_N)] = *collision_nearest;
                (*collisions_N)++;
            }
        }
    }else{        
        // c is not a leaf node
        double dx = gb.shiftx - c->x;
        double dy = gb.shifty - c->y;
        double dz = gb.shiftz - c->z;
        double r2 = dx*dx + dy*dy + dz*dz;
        double rp  = p1_r + r->max_radius[1] + 0.86602540378443*c->w;
        // Check if we need to decent into daughter cells
        if (r2 < rp*rp ){
            for (int o=0;o<8;o++){
                struct reb_treecell* d = c->oct[o];
                if (d!=NULL){
                    reb_tree_get_nearest_neighbour_in_cell(r, collisions_N, gb,gbunmod,ri,p1_r,nearest_r2,collision_nearest,d);
                }
            }
        }
    }
}


static void reb_tree_check_for_overlapping_trajectories_in_cell(struct reb_simulation* const r, int* collisions_N, struct reb_ghostbox gb, struct reb_ghostbox gbunmod, int ri, double p1_r, double p1_r_plus_dtv, struct reb_collision* collision_nearest, struct reb_treecell* c, double maxdrift){
    const struct reb_particle* const particles = r->particles;
    if (c->pt>=0){     
        // c is a leaf node
        if (c->pt != collision_nearest->p1){
            struct reb_particle p2 = particles[c->pt];
            double dt_done_last = r->dt_last_done;
            const double dx1 = gb.shiftx - p2.x; // distance at beginning
            const double dy1 = gb.shifty - p2.y;
            const double dz1 = gb.shiftz - p2.z;
            const double r1 = (dx1*dx1 + dy1*dy1 + dz1*dz1);
            const double dvx1 = gb.shiftvx - p2.vx; 
            const double dvy1 = gb.shiftvy - p2.vy;
            const double dvz1 = gb.shiftvz - p2.vz;
            const double dx2 = dx1 -dt_done_last*dvx1; // distance at end
            const double dy2 = dy1 -dt_done_last*dvy1;
            const double dz2 = dz1 -dt_done_last*dvz1;
            const double r2 = (dx2*dx2 + dy2*dy2 + dz2*dz2);
            const double t_closest = (dx1*dvx1 + dy1*dvy1 + dz1*dvz1)/(dvx1*dvx1 + dvy1*dvy1 + dvz1*dvz1);

            double rmin2_ab = MIN(r1,r2);
            if (t_closest/dt_done_last>=0. && t_closest/dt_done_last<=1.){
                const double dx3 = dx1-t_closest*dvx1; // closest approach
                const double dy3 = dy1-t_closest*dvy1;
                const double dz3 = dz1-t_closest*dvz1;
                const double r3 = (dx3*dx3 + dy3*dy3 + dz3*dz3);
                rmin2_ab = MIN(rmin2_ab, r3);
            }
            double rsum = p1_r + p2.r;
            if (rmin2_ab>rsum*rsum) return;
            collision_nearest->ri = ri;
            collision_nearest->p2 = c->pt;
            collision_nearest->gb = gbunmod;
            // Save collision in collisions array.
#pragma omp critical
            {
                if (r->collisions_allocatedN<=(*collisions_N)){
                    // Init to 32 if no space has been allocated yet, otherwise double it.
                    r->collisions_allocatedN = r->collisions_allocatedN ? r->collisions_allocatedN * 2 : 32;
                    r->collisions = realloc(r->collisions,sizeof(struct reb_collision)*r->collisions_allocatedN);
                }
                r->collisions[(*collisions_N)] = *collision_nearest;
                (*collisions_N)++;
            }
        }
    }else{        
        // c is not a leaf node
        double dx = gb.shiftx - c->x;
        double dy = gb.shifty - c->y;
        double dz = gb.shiftz - c->z;
        double r2 = dx*dx + dy*dy + dz*dz;
        double rp  = p1_r_plus_dtv + maxdrift + 0.86602540378443*c->w;
        // Check if we need to decent into daughter cells
        if (r2 < rp*rp ){
            for (int o=0;o<8;o++){
                struct reb_treecell* d = c->oct[o];
                if (d!=NULL){
                    reb_tree_check_for_overlapping_trajectories_in_cell(r, collisions_N, gb,gbunmod,ri,p1_r,p1_r_plus_dtv,collision_nearest,d,maxdrift);
                }
            }
        }
    }
}




int reb_collision_resolve_hardsphere(struct reb_simulation* const r, struct reb_collision c){
    struct reb_particle* const particles = r->particles;
    struct reb_particle p1 = particles[c.p1];
    struct reb_particle p2;
#ifdef MPI
    int isloc = reb_communication_mpi_rootbox_is_local(r, c.ri);
    if (isloc==1){
#endif // MPI
        p2 = particles[c.p2];
#ifdef MPI
    }else{
        int root_n_per_node = r->root_n/r->mpi_num;
        int proc_id = c.ri/root_n_per_node;
        p2 = r->particles_recv[proc_id][c.p2];
    }
#endif // MPI
//    if (p1.lastcollision==t || p2.lastcollision==t) return;
    struct reb_ghostbox gb = c.gb;
    double x21  = p1.x + gb.shiftx  - p2.x; 
    double y21  = p1.y + gb.shifty  - p2.y; 
    double z21  = p1.z + gb.shiftz  - p2.z; 
    double rp   = p1.r+p2.r;
    double oldvyouter;
    if (x21>0){
         oldvyouter = p1.vy;
    }else{
        oldvyouter = p2.vy;
    }
    if (rp*rp < x21*x21 + y21*y21 + z21*z21) return 0;
    double vx21 = p1.vx + gb.shiftvx - p2.vx; 
    double vy21 = p1.vy + gb.shiftvy - p2.vy; 
    double vz21 = p1.vz + gb.shiftvz - p2.vz; 
    if (vx21*x21 + vy21*y21 + vz21*z21 >0) return 0; // not approaching
    // Bring the to balls in the xy plane.
    // NOTE: this could probabely be an atan (which is faster than atan2)
    double theta = atan2(z21,y21);
    double stheta = sin(theta);
    double ctheta = cos(theta);
    double vy21n = ctheta * vy21 + stheta * vz21;    
    double y21n = ctheta * y21 + stheta * z21;    
    
    // Bring the two balls onto the positive x axis.
    double phi = atan2(y21n,x21);
    double cphi = cos(phi);
    double sphi = sin(phi);
    double vx21nn = cphi * vx21  + sphi * vy21n;        

    // Coefficient of restitution
    double eps= 1; // perfect bouncing by default 
    if (r->coefficient_of_restitution){
        eps = r->coefficient_of_restitution(r, vx21nn);
    }
    double dvx2 = -(1.0+eps)*vx21nn;
    double minr = (p1.r>p2.r)?p2.r:p1.r;
    double maxr = (p1.r<p2.r)?p2.r:p1.r;
    double mindv= minr*r->minimum_collision_velocity;
    double _r = sqrt(x21*x21 + y21*y21 + z21*z21);
    mindv *= 1.-(_r - maxr)/minr;
    if (mindv>maxr*r->minimum_collision_velocity)mindv = maxr*r->minimum_collision_velocity;
    if (dvx2<mindv) dvx2 = mindv;
    // Now we are rotating backwards
    double dvx2n = cphi * dvx2;        
    double dvy2n = sphi * dvx2;        
    double dvy2nn = ctheta * dvy2n;    
    double dvz2nn = stheta * dvy2n;    


    // Applying the changes to the particles.
#ifdef MPI
    if (isloc==1){
#endif // MPI
    const double p2pf = p1.m/(p1.m+p2.m);
    particles[c.p2].vx -=    p2pf*dvx2n;
    particles[c.p2].vy -=    p2pf*dvy2nn;
    particles[c.p2].vz -=    p2pf*dvz2nn;
    particles[c.p2].lastcollision = r->t;
#ifdef MPI
    }
#endif // MPI
    const double p1pf = p2.m/(p1.m+p2.m);
    particles[c.p1].vx +=    p1pf*dvx2n; 
    particles[c.p1].vy +=    p1pf*dvy2nn; 
    particles[c.p1].vz +=    p1pf*dvz2nn; 
    particles[c.p1].lastcollision = r->t;
        
    // Return y-momentum change
    if (x21>0){
        r->collisions_plog += -fabs(x21)*(oldvyouter-particles[c.p1].vy) * p1.m;
        r->collisions_Nlog ++;
    }else{
        r->collisions_plog += -fabs(x21)*(oldvyouter-particles[c.p2].vy) * p2.m;
        r->collisions_Nlog ++;
    }
    return 0;
}

int reb_collision_resolve_halt(struct reb_simulation* const r, struct reb_collision c){
    r->status = REB_EXIT_COLLISION;
    r->particles[c.p1].lastcollision = r->t;
    r->particles[c.p2].lastcollision = r->t;
    return 0; // don't remove either particle
}

int reb_collision_resolve_merge(struct reb_simulation* const r, struct reb_collision c){
    if (r->particles[c.p1].lastcollision==r->t || r->particles[c.p2].lastcollision==r->t) return 0;

    // Every collision will cause two callbacks (with p1/p2 interchanged).
    // Always remove particle with larger index and merge into lower index particle.
    // This will keep N_active meaningful even after mergers.
    int swap = 0;
    int i = c.p1;
    int j = c.p2;   //want j to be removed particle
    if (j<i){
        swap = 1;
        i = c.p2;
        j = c.p1;
    }

    struct reb_particle* pi = &(r->particles[i]);
    struct reb_particle* pj = &(r->particles[j]);
                
    double invmass = 1.0/(pi->m + pj->m);
    
    //Scale out energy from collision - initial energy
    double Ei=0, Ef=0;
    if(r->track_energy_offset){
        {
            double vx = pi->vx;
            double vy = pi->vy;
            double vz = pi->vz;
            // Calculate energy difference in inertial frame
            if (r->integrator == REB_INTEGRATOR_MERCURIUS && r->ri_mercurius.mode==1){
                vx += r->ri_mercurius.com_vel.x;
                vy += r->ri_mercurius.com_vel.y;
                vz += r->ri_mercurius.com_vel.z;
            }

            Ei += 0.5*pi->m*(vx*vx + vy*vy + vz*vz);
        }
        {
            double vx = pj->vx;
            double vy = pj->vy;
            double vz = pj->vz;
            if (r->integrator == REB_INTEGRATOR_MERCURIUS && r->ri_mercurius.mode==1){
                vx += r->ri_mercurius.com_vel.x;
                vy += r->ri_mercurius.com_vel.y;
                vz += r->ri_mercurius.com_vel.z;
            }

            Ei += 0.5*pj->m*(vx*vx + vy*vy + vz*vz);
        }
        const int N_active = ((r->N_active==-1)?r->N-r->N_var:r->N_active);
        // No potential energy between test particles
        if (i<N_active || j<N_active){
            double x = pi->x - pj->x;
            double y = pi->y - pj->y;
            double z = pi->z - pj->z;
            double _r = sqrt(x*x + y*y + z*z);

            Ei += - r->G*pi->m*pj->m/_r;
        }
    }
    
    // Merge by conserving mass, volume and momentum
    pi->vx = (pi->vx*pi->m + pj->vx*pj->m)*invmass;
    pi->vy = (pi->vy*pi->m + pj->vy*pj->m)*invmass;
    pi->vz = (pi->vz*pi->m + pj->vz*pj->m)*invmass;
    pi->x  = (pi->x*pi->m + pj->x*pj->m)*invmass;
    pi->y  = (pi->y*pi->m + pj->y*pj->m)*invmass;
    pi->z  = (pi->z*pi->m + pj->z*pj->m)*invmass;
    pi->m  = pi->m + pj->m;
    pi->r  = pow(pow(pi->r,3.)+pow(pj->r,3.),1./3.);
    pi->lastcollision = r->t;
    

    // Keeping track of energy offst
    if(r->track_energy_offset){
        {
            double vx = pi->vx;
            double vy = pi->vy;
            double vz = pi->vz;
            if (r->integrator == REB_INTEGRATOR_MERCURIUS && r->ri_mercurius.mode==1){
                vx += r->ri_mercurius.com_vel.x;
                vy += r->ri_mercurius.com_vel.y;
                vz += r->ri_mercurius.com_vel.z;
            }

            Ef += 0.5*pi->m*(vx*vx + vy*vy + vz*vz);
        }
        r->energy_offset += Ei - Ef;
    }
    
    return swap?1:2; // Remove particle p2 from simulation
}
