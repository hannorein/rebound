#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>


struct reb_orbit_node {
    struct reb_orbit_node* primary;
    struct reb_orbit_node* secondary;
    struct reb_particle com;
    double P;
};

double period(double G, struct reb_orbit_node* p1, struct reb_orbit_node* p2){
    double mu,dx,dy,dz,dvx,dvy,dvz,vsquared,vcircsquared,n,a,d;
    mu = G*(p1->com.m+p2->com.m);
    dx = p1->com.x - p2->com.x;
    dy = p1->com.y - p2->com.y;
    dz = p1->com.z - p2->com.z;
    dvx = p1->com.vx - p2->com.vx;
    dvy = p1->com.vy - p2->com.vy;
    dvz = p1->com.vz - p2->com.vz;
    d = sqrt ( dx*dx + dy*dy + dz*dz );
    vsquared = dvx*dvx + dvy*dvy + dvz*dvz;
    vcircsquared = mu/d;	
    a = -mu/( vsquared - 2.*vcircsquared );	// semi major axis
    n = a/fabs(a)*sqrt(fabs(mu/(a*a*a)));	// mean motion (negative if hyperbolic)
    return 2*M_PI/n;						// period (negative if hyperbolic)
}
    
double eccentricity(double G, struct reb_orbit_node* p1, struct reb_orbit_node* p2){
    double mu,dx,dy,dz,dvx,dvy,dvz,vsquared,vcircsquared,d;
    double ex, ey, ez, vr, rvr, vdiffsquared, muinv;
    mu = G*(p1->com.m+p2->com.m);
    dx = p1->com.x - p2->com.x;
    dy = p1->com.y - p2->com.y;
    dz = p1->com.z - p2->com.z;
    dvx = p1->com.vx - p2->com.vx;
    dvy = p1->com.vy - p2->com.vy;
    dvz = p1->com.vz - p2->com.vz;
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

void free_orbit_node(struct reb_orbit_node* on) {
    if (!on) return;
    free_orbit_node(on->primary);
    free_orbit_node(on->secondary);
    free(on);
}

struct reb_orbit_node* hierarchy(struct reb_simulation* r) {
    int N = r->N; // number of objects to be sorted
    struct reb_orbit_node** to_be_sorted = calloc(N,sizeof(struct reb_orbit_node*));
    for(int i=0;i<N;i++){
        to_be_sorted[i] = calloc(1,sizeof(struct reb_orbit_node));
        to_be_sorted[i]->com = r->particles[i];
    }
    while (N>1){
        int imin = -1;
        int jmin = -1;
        double Pmin = INFINITY;
        for(int i=1;i<N;i++){
            for(int j=0;j<i;j++){
                double Pij = period(r->G, to_be_sorted[i], to_be_sorted[j]);
                if (Pij>0 && Pij<Pmin){
                    Pmin = Pij;
                    imin = i; jmin=j;
                }
            }
        }
        if (Pmin < INFINITY){
            struct reb_orbit_node* on = calloc(1,sizeof(struct reb_orbit_node));
            if (to_be_sorted[imin]->com.m>to_be_sorted[jmin]->com.m){ 
                on->primary = to_be_sorted[imin]; 
                on->secondary = to_be_sorted[jmin]; 
            }else{
                on->primary = to_be_sorted[jmin]; 
                on->secondary = to_be_sorted[imin]; 
            }
            on->P = Pmin;
            on->com = reb_particle_com_of_pair(on->primary->com, on->secondary->com); 
            to_be_sorted[imin] = on;
            to_be_sorted[jmin] = to_be_sorted[N-1];
            N--;
        }else{
            printf("Error: Cannot find bound orbit.\n");
            for(int i=0;i<N;i++){
                free_orbit_node(to_be_sorted[i]);
            }
            free(to_be_sorted);
            return NULL;
        }
    }
    struct reb_orbit_node* root = to_be_sorted[0];
    free(to_be_sorted);
    return root;
}

void print_orbit_node(struct reb_simulation* r, struct reb_orbit_node* on, int level) {
    if (!on) {
        printf("NULL Orbit Node\n");
        return;
    }
    for (int i=0;i<level;i++){
        printf("  ");
    }
    if (on->primary && on->secondary){
        double e = eccentricity(r->G, on->primary, on->secondary);
        printf("binary   m= %f    P = %f   e = %f",on->com.m, on->P, e);
    }else{
        int pid = -1;
        for (int i=0; i<r->N;i++){
            if (!reb_particle_diff(on->com, r->particles[i])){
                pid = i;
                break;
            }
        }
        printf("particle id=%d   m= %f    x = %f y = %f", pid, on->com.m, on->com.x, on->com.y);
    }
    printf("\n");
    if (on->primary){
        print_orbit_node(r, on->primary, level+1);
    }
    if (on->secondary){
        print_orbit_node(r, on->secondary, level+1);
    }
}

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_simulation_create();
    
    reb_simulation_add_fmt(r, "m", 1.);
    reb_simulation_add_fmt(r, "m P primary", 0.01, 2.0, r->particles[0]);
    reb_simulation_add_fmt(r, "m P primary", 0.01, 1.0, r->particles[0]);
    reb_simulation_add_fmt(r, "m a e", 0.01, -1.0, 9.3);

//    reb_simulation_add_fmt(r, "m", 1.);
//    reb_simulation_add_fmt(r, "m P e", 1e-3, 9., 0.1);
//    reb_simulation_add_fmt(r, "m P e", 1e-3, 9., 0.1);
//    reb_simulation_add_fmt(r, "P e", 10.0, 0.1);
//    reb_simulation_add_fmt(r, "m P e primary", 1e-6, 4.0, 0.1, r->particles[1]);

    struct reb_orbit_node* on = hierarchy(r);
    print_orbit_node(r, on,0);
    free_orbit_node(on);


    reb_simulation_free(r);
}

