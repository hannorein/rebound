/**
 * MPI Unittests
 *
 * This file contains a couple of unit tests for REBOUND
 * using MPI.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <assert.h>
#include "rebound.h"
#include "tools.h"


void test_twobody(){
    struct reb_simulation* const r = reb_simulation_create();
    r->integrator       = REB_INTEGRATOR_LEAPFROG;
    r->gravity          = REB_GRAVITY_TREE;
    r->boundary         = REB_BOUNDARY_OPEN;
    r->opening_angle2   = 1.5; 
    r->G                = 1;        
    r->dt               = 0.1;
    reb_simulation_configure_box(r,10,2,2,1);

    printf("MPI init...\n");
    reb_mpi_init(r);
    if (r->mpi_id==0){
        reb_simulation_add_fmt(r, "m hash", 2., reb_hash("star1"));
    }
    struct reb_particle com = reb_simulation_com(r); // Need to call this on all machines. 
    if (r->mpi_id==0){
        reb_simulation_add_fmt(r, "m a e primary hash", 1., 1., 0.1, com, reb_hash("star2"));
    }
    printf("Init done. (%d)   N = %d\n", r->mpi_id, r->N);

    printf("Moving to com...\n"); // Will also distribute particles
    reb_simulation_move_to_com(r);
    
    printf("Checking com...\n");
    com = reb_simulation_com(r);
    assert(fabs(com.x)<1e-15);
    assert(fabs(com.y)<1e-15);
    assert(fabs(com.z)<1e-15);
    assert(fabs(com.vx)<1e-15);
    assert(fabs(com.vy)<1e-15);
    assert(fabs(com.vz)<1e-15);

    printf("Starting the integration...\n");
    reb_simulation_integrate(r, 10.);
    
    printf("Checking conservation of orbital elements...\n");
    struct reb_particle star1 = reb_simulation_particle_by_hash_mpi(r, reb_hash("star1"));
    struct reb_particle star2 = reb_simulation_particle_by_hash_mpi(r, reb_hash("star2"));
    struct reb_orbit o = reb_orbit_from_particle(r->G, star2, star1);
    
    assert(fabs(o.a-1.)<1e-3);
    assert(fabs(o.e-0.1)<1e-2);

    printf("Checking input/output...\n");

    com = reb_simulation_com(r); // Need to call this on all machines. 
    if (r->mpi_id==0){
        for (int i=0; i<10; i++){
            reb_simulation_add_fmt(r, "m a primary hash", 0.01, 2.0+0.1*i, com, i);
        }
    }
    reb_simulation_steps(r, 1);
    {
        // Delete any previous files
        char filename[1024];
        sprintf(filename, "out.bin_%d", r->mpi_id);
        remove(filename);
    }
    reb_simulation_save_to_file(r, "out.bin");
    reb_simulation_steps(r, 1);
    reb_simulation_save_to_file(r, "out.bin");
    reb_simulation_steps(r, 10);

    struct reb_simulationarchive* sa = reb_simulationarchive_create_from_file("out.bin");
    assert(sa->nblobs==2);
    struct reb_simulation* r2 = reb_simulation_create_from_simulationarchive(sa,-1);
    reb_simulation_steps(r2, 10);

    assert(r->N == r2->N);
    assert(r->t == r2->t);

    // Order of particles will be different. Need to compare them by hash
    for(int i=0; i<10; i++){
        struct reb_particle p1 = reb_simulation_particle_by_hash_mpi(r, i);
        struct reb_particle p2 = reb_simulation_particle_by_hash_mpi(r2, i);
        assert(p1.x==p2.x);
        assert(p1.y==p2.y);
        assert(p1.z==p2.z);
    }




    printf("Cleanup...\n");
    reb_mpi_finalize(r);
    reb_simulation_free(r); 
    reb_simulation_free(r2); 


}

int main(int argc, char* argv[]){
    test_twobody();
}

