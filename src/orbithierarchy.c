/**
 * @file 	orbithierarchy.c
 * @brief 	Tools for creating and working with a tree that represents the orbital architecture.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2025 Hanno Rein
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
#include "rebound.h"

struct reb_orbit reb_orbit_from_orbithierarchy(double G, struct reb_orbithierarchy* oh){
    if (oh->primary==NULL || oh->secondary==NULL){
        printf("Orbithierarchy is a single particle, not an orbit.");
        return reb_orbit_nan();
    }
    return reb_orbit_from_particle(G, *(oh->secondary->com), *(oh->primary->com));
}
    
double reb_orbithierarchy_eccentricity(double G, struct reb_orbithierarchy* p1, struct reb_orbithierarchy* p2){
    double mu,dx,dy,dz,dvx,dvy,dvz,vsquared,vcircsquared,d;
    double ex, ey, ez, vr, rvr, vdiffsquared, muinv;
    mu = G*(p1->com->m+p2->com->m);
    dx = p1->com->x - p2->com->x;
    dy = p1->com->y - p2->com->y;
    dz = p1->com->z - p2->com->z;
    dvx = p1->com->vx - p2->com->vx;
    dvy = p1->com->vy - p2->com->vy;
    dvz = p1->com->vz - p2->com->vz;
    d = sqrt ( dx*dx + dy*dy + dz*dz );
    vsquared = dvx*dvx + dvy*dvy + dvz*dvz;
    vcircsquared = mu/d;	

    vdiffsquared = vsquared - vcircsquared;	
    vr = (dx*dvx + dy*dvy + dz*dvz)/d;	
    rvr = d*vr;
    muinv = 1./mu;

    ex = muinv*( vdiffsquared*dx - rvr*dvx );
    ey = muinv*( vdiffsquared*dy - rvr*dvy );
    ez = muinv*( vdiffsquared*dz - rvr*dvz );
    return sqrt( ex*ex + ey*ey + ez*ez );		// eccentricity
}

void reb_orbithierarchy_free(struct reb_orbithierarchy* oh) {
    if (!oh) return;
    if (oh->primary || oh->secondary) free(oh->com); // not a leaf node
    reb_orbithierarchy_free(oh->primary);
    reb_orbithierarchy_free(oh->secondary);
    free(oh);
}

struct reb_orbithierarchy* reb_simulation_create_orbithierarchy(struct reb_simulation* r) {
    int N = r->N; // number of objects to be sorted
    struct reb_orbithierarchy** to_be_sorted = calloc(N,sizeof(struct reb_orbithierarchy*));
    for(int i=0;i<N;i++){
        to_be_sorted[i] = calloc(1,sizeof(struct reb_orbithierarchy));
        to_be_sorted[i]->com = &(r->particles[i]);
    }
    while (N>1){
        int imin = -1;
        int jmin = -1;
        double Pmin = INFINITY;
        for(int i=1;i<N;i++){
            for(int j=0;j<i;j++){
                double Pij = reb_orbit_from_particle(r->G, *(to_be_sorted[i]->com), *(to_be_sorted[j]->com)).P;
                if (Pij>0 && Pij<Pmin){
                    Pmin = Pij;
                    imin = i; jmin=j;
                }
            }
        }
        if (Pmin < INFINITY){
            struct reb_orbithierarchy* oh = calloc(1,sizeof(struct reb_orbithierarchy));
            if (to_be_sorted[imin]->com->m>to_be_sorted[jmin]->com->m){ 
                oh->primary = to_be_sorted[imin]; 
                oh->secondary = to_be_sorted[jmin]; 
            }else{
                oh->primary = to_be_sorted[jmin]; 
                oh->secondary = to_be_sorted[imin]; 
            }
            oh->com = malloc(sizeof(struct reb_particle));
            *(oh->com) = reb_particle_com_of_pair(*(oh->primary->com), *(oh->secondary->com));
            to_be_sorted[imin] = oh;
            to_be_sorted[jmin] = to_be_sorted[N-1];
            N--;
        }else{
            printf("Error: Cannot find bound orbit.\n");
            for(int i=0;i<N;i++){
                reb_orbithierarchy_free(to_be_sorted[i]);
            }
            free(to_be_sorted);
            return NULL;
        }
    }
    struct reb_orbithierarchy* root = to_be_sorted[0];
    free(to_be_sorted);
    return root;
}

void reb_orbithierarchy_print(struct reb_orbithierarchy* oh, struct reb_simulation* r, int level) {
    if (!oh) {
        printf("NULL Orbit Node\n");
        return;
    }
    for (int i=0;i<level;i++){
        printf("  ");
    }
    if (oh->primary && oh->secondary){
        struct reb_orbit o = reb_orbit_from_orbithierarchy(r->G, oh);
        printf("binary   m= %f    P = %f   e = %f",oh->com->m, o.P, o.e);
    }else{
        int pid = -1;
        for (int i=0; i<r->N;i++){
            if (oh->com == &(r->particles[i])){
                pid = i;
                break;
            }
        }
        printf("particle id=%d   m= %f    x = %f y = %f", pid, oh->com->m, oh->com->x, oh->com->y);
    }
    printf("\n");
    if (oh->primary){
        reb_orbithierarchy_print(oh->primary, r, level+1);
    }
    if (oh->secondary){
        reb_orbithierarchy_print(oh->secondary, r, level+1);
    }
}

int reb_orbithierarchy_is_jacobi(struct reb_orbithierarchy* oh) {
    if (!oh) return 0;
    int is_particle = oh->primary==NULL && oh->secondary==NULL;
    if (is_particle) return 1;
    int primary_is_particle = oh->primary->primary==NULL && oh->primary->secondary==NULL;
    int secondary_is_particle = oh->secondary->primary==NULL && oh->secondary->secondary==NULL;
    if (primary_is_particle && secondary_is_particle) return 1;
    if (!secondary_is_particle) return 0;
    return reb_orbithierarchy_is_jacobi(oh->primary);
}

int reb_orbithierarchy_is_jacobi_ordered(struct reb_simulation* r, struct reb_orbithierarchy* oh) {
    if (reb_orbithierarchy_is_jacobi(oh)==0) return 0;
    for (int i=r->N-1;i>=1;i--){
        if (oh->secondary->com != &(r->particles[i])) return 0;
        oh = oh->primary;
    }
    if (oh->com != &(r->particles[0])) return 0;
    return 1;
}
