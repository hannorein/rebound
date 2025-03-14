/**
 * Close encounters with simple MERCURIUS (Hernandez 2023)
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

int test_collision(struct reb_simulation* const r, struct reb_collision c){
    // Simple collision. Two particles become three, conserving mass and angular momentum
    // remove second particle
    unsigned int i = c.p1;
    unsigned int j = c.p2;   //want j to be removed particle

    struct reb_particle* pi = &(r->particles[i]);
    struct reb_particle* pj = &(r->particles[j]);
    struct reb_particle com = reb_particle_com_of_pair(*pi, *pj);

    // total mass
    double mtot = pi->m + pj->m;
    double rtot = pi->r + pj->r;
    double mfrag = mtot/3.;

    // angular momentum
    double ptot_x = pi->m * pi->vx + pj->m * pj->vx;
    double ptot_y = pi->m * pi->vy + pj->m * pj->vy;
    double ptot_z = pi->m * pi->vz + pj->m * pj->vz;

    // define vector
    struct reb_vec3d ptot = {.x=ptot_x, .y=ptot_y, .z=ptot_z};
    double p = sqrt(reb_vec3d_length_squared(ptot));
    double pfrag = p/3.;
    struct reb_vec3d phat = reb_vec3d_normalize(ptot);

    double vfrag_mag = pfrag/mfrag;

    // first particle in direction of momentum
    struct reb_particle p1 = {0};
    p1.m = mfrag;
    p1.x = com.x + rtot * phat.x;
    p1.y = com.y + rtot * phat.y;
    p1.z = com.y + rtot * phat.z;
    p1.vx = vfrag_mag * phat.x;
    p1.vy = vfrag_mag * phat.y;
    p1.vz = vfrag_mag * phat.z;

    // other two perpendicular
    struct reb_rotation r2 = {.ix=0,.iy=0,.iz=1,.r=M_PI/2.};
    struct reb_vec3d p2_hat = reb_vec3d_rotate(phat, r2);
    struct reb_particle p2 = {0};
    p2.m = mfrag;
    p2.x = com.x + rtot * p2_hat.x;
    p2.y = com.y + rtot * p2_hat.y;
    p2.z = com.y + rtot * p2_hat.z;
    p2.vx = vfrag_mag * p2_hat.x;
    p2.vy = vfrag_mag * p2_hat.y;
    p2.vz = vfrag_mag * p2_hat.z;
    
    struct reb_rotation r3 = {.ix=0,.iy=0,.iz=1,.r=-1.*M_PI/2.};
    struct reb_vec3d p3_hat = reb_vec3d_rotate(phat, r3);
    struct reb_particle p3 = {0};
    p3.m = mfrag;
    p3.x = com.x + rtot * p3_hat.x;
    p3.y = com.y + rtot * p3_hat.y;
    p3.z = com.y + rtot * p3_hat.z;
    p3.vx = vfrag_mag * p3_hat.x;
    p3.vy = vfrag_mag * p3_hat.y;
    p3.vz = vfrag_mag * p3_hat.z;
    /*
    // Scale out energy from collision
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

            if (r->integrator == REB_INTEGRATOR_TRACE && r->ri_trace.mode==REB_TRACE_MODE_KEPLER){
                vx += r->ri_trace.com_vel.x;
                vy += r->ri_trace.com_vel.y;
                vz += r->ri_trace.com_vel.z;
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

            if (r->integrator == REB_INTEGRATOR_TRACE && r->ri_trace.mode==REB_TRACE_MODE_KEPLER){
                vx += r->ri_trace.com_vel.x;
                vy += r->ri_trace.com_vel.y;
                vz += r->ri_trace.com_vel.z;
            }

            Ei += 0.5*pj->m*(vx*vx + vy*vy + vz*vz);
        }
        
	double x = pi->x - pj->x;
        double y = pi->y - pj->y;
        double z = pi->z - pj->z;
        double _r = sqrt(x*x + y*y + z*z);
        Ei += - r->G*pi->m*pj->m/_r;
	
	// Final energy
	{
            double vx = p1.vx;
            double vy = p1.vy;
            double vz = p1.vz;
            if (r->integrator == REB_INTEGRATOR_MERCURIUS && r->ri_mercurius.mode==1){
                vx += r->ri_mercurius.com_vel.x;
                vy += r->ri_mercurius.com_vel.y;
                vz += r->ri_mercurius.com_vel.z;
            }

            if (r->integrator == REB_INTEGRATOR_TRACE && r->ri_trace.mode==REB_TRACE_MODE_KEPLER){
                vx += r->ri_trace.com_vel.x;
                vy += r->ri_trace.com_vel.y;
                vz += r->ri_trace.com_vel.z;
            }

            Ef += 0.5*p1.m*(vx*vx + vy*vy + vz*vz);
        }
	{
            double vx = p2.vx;
            double vy = p2.vy;
            double vz = p2.vz;
            if (r->integrator == REB_INTEGRATOR_MERCURIUS && r->ri_mercurius.mode==1){
                vx += r->ri_mercurius.com_vel.x;
                vy += r->ri_mercurius.com_vel.y;
                vz += r->ri_mercurius.com_vel.z;
            }

            if (r->integrator == REB_INTEGRATOR_TRACE && r->ri_trace.mode==REB_TRACE_MODE_KEPLER){
                vx += r->ri_trace.com_vel.x;
                vy += r->ri_trace.com_vel.y;
                vz += r->ri_trace.com_vel.z;
            }

            Ef += 0.5*p2.m*(vx*vx + vy*vy + vz*vz);
        }
	{
            double vx = p3.vx;
            double vy = p3.vy;
            double vz = p3.vz;
            if (r->integrator == REB_INTEGRATOR_MERCURIUS && r->ri_mercurius.mode==1){
                vx += r->ri_mercurius.com_vel.x;
                vy += r->ri_mercurius.com_vel.y;
                vz += r->ri_mercurius.com_vel.z;
            }

            if (r->integrator == REB_INTEGRATOR_TRACE && r->ri_trace.mode==REB_TRACE_MODE_KEPLER){
                vx += r->ri_trace.com_vel.x;
                vy += r->ri_trace.com_vel.y;
                vz += r->ri_trace.com_vel.z;
            }

            Ef += 0.5*p3.m*(vx*vx + vy*vy + vz*vz);
        }
        
	double x12 = p1.x - p2.x;
        double y12 = p1.y - p2.y;
        double z12 = p1.z - p2.z;
	double x23 = p2.x - p3.x;
        double y23 = p2.y - p3.y;
        double z23 = p2.z - p3.z;
	double x13 = p1.x - p3.x;
        double y13 = p1.y - p3.y;
        double z13 = p1.z - p3.z;
        double _r12 = sqrt(x12*x12 + y12*y12 + z12*z12);
        double _r23 = sqrt(x23*x23 + y23*y23 + z23*z23);
        double _r13 = sqrt(x13*x13 + y13*y13 + z13*z13);
        Ef += - r->G*p1.m*p2.m/_r12;
        Ef += - r->G*p2.m*p3.m/_r23;
        Ef += - r->G*p1.m*p3.m/_r13;
        r->energy_offset += Ei - Ef;
    }
    */
    
    double E0 = reb_simulation_energy(r);
    
    reb_simulation_add(r, p1);
    reb_simulation_add(r, p2);
    reb_simulation_add(r, p3);
    reb_simulation_remove_particle(r, 2, 1);
    reb_simulation_remove_particle(r, 1, 1);
    
    r->collisions[1].p1 = -1;
    r->collisions[1].p2 = -1;

    if (r->track_energy_offset){
        r->energy_offset += E0 - reb_simulation_energy(r);
    }
    return 0;
}

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
	
    struct reb_particle p = {0};
    p.m = 1.;
    reb_simulation_add(r, p);
    reb_simulation_add_fmt(r, "m r a e", 1e-5, 1.6e-4, 0.5, 0.1);
    reb_simulation_add_fmt(r, "m r a e f", 1e-8, 4e-5, 0.55, 0.4, -0.94);

    r->integrator = REB_INTEGRATOR_MERCURIUS;
    r->dt = 0.01;
    r->collision = REB_COLLISION_DIRECT;
    //r->collision_resolve = reb_collision_resolve_merge;
    r->collision_resolve = test_collision;
    r->track_energy_offset = 1;
    //r->heartbeat  = heartbeat;

    double e0 = reb_simulation_energy(r);
    printf("%d\n", r->N); 
    reb_simulation_integrate(r, 1.);
    printf("%d %e\n", r->N, fabs((reb_simulation_energy(r) - e0) / e0));
    reb_simulation_free(r);
}

