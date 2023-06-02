/**
 * MPI Unittests
 *
 * This file contains a couple of unit tests for REBOUND
 * using MPI.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <assert.h>
#include "rebound.h"
#include "tools.h"


void test_twobody(){
    struct reb_simulation* const r = reb_create_simulation();
    r->integrator       = REB_INTEGRATOR_LEAPFROG;
    r->gravity          = REB_GRAVITY_TREE;
    r->boundary         = REB_BOUNDARY_OPEN;
    r->opening_angle2   = 1.5; 
    r->G                = 1;        
    r->dt               = 0.1;
    reb_configure_box(r,10,2,2,1);

    printf("MPI init...\n");
    reb_mpi_init(r);
    if (r->mpi_id==0){
        reb_add_fmt(r, "m hash", 2., reb_hash("star1"));
    }
    struct reb_particle com = reb_get_com(r); // Need to call this on all machines. 
    if (r->mpi_id==0){
        reb_add_fmt(r, "m a e primary hash", 1., 1., 0.1, com, reb_hash("star2"));
    }
    printf("Init done. (%d)   N = %d\n", r->mpi_id, r->N);

    printf("Moving to com...\n"); // Will also distribute particles
    reb_move_to_com(r);
    
    printf("Checking com...\n");
    com = reb_get_com(r);
    assert(fabs(com.x)<1e-15);
    assert(fabs(com.y)<1e-15);
    assert(fabs(com.z)<1e-15);
    assert(fabs(com.vx)<1e-15);
    assert(fabs(com.vy)<1e-15);
    assert(fabs(com.vz)<1e-15);

    printf("Starting the integration...\n");
    reb_integrate(r, 10.);
    
    printf("Checking conservation of orbital elements...\n");
    struct reb_particle star1 = reb_get_remote_particle_by_hash(r, reb_hash("star1"));
    struct reb_particle star2 = reb_get_remote_particle_by_hash(r, reb_hash("star2"));
    struct reb_orbit o = reb_tools_particle_to_orbit(r->G, star2, star1);
    
    assert(fabs(o.a-1.)<1e-3);
    assert(fabs(o.e-0.1)<1e-2);

    printf("Cleanup...\n");
    reb_mpi_finalize(r);
    reb_free_simulation(r); 


}

int main(int argc, char* argv[]){
    test_twobody();
}

