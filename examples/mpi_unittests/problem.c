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
        reb_add_fmt(r, "m", 1.);
    }
    struct reb_particle com = reb_get_com(r); // Need to call this on all machines. 
    if (r->mpi_id==0){
        reb_add_fmt(r, "m a e primary", 1., 1., 0.1, com);
    }
    printf("Init done. (%d)   N = %d\n", r->mpi_id, r->N);

    printf("Moving to com...\n");
    reb_move_to_com(r);
    
    printf("Checking com...\n");
    com = reb_get_com(r);
    assert(com.x==0);
    assert(com.y==0);
    assert(com.z==0);
    assert(com.vx==0);
    assert(com.vy==0);
    assert(com.vz==0);

    printf("Starting the integration...\n");
    reb_integrate(r, 10.);

    printf("Cleanup...\n");
    reb_mpi_finalize(r);
    reb_free_simulation(r); 


}

int main(int argc, char* argv[]){
    test_twobody();
}

