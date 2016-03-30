/**
 * A very simple test problem
 * 
 * We first create a REBOUND simulation, then we add 
 * two particles and integrate the system for 100 time 
 * units.
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>


int main(int argc, char* argv[]) {
    struct reb_simulation* sim = reb_create_simulation();

    struct reb_particle p = {0.};
    p.m = 1;
    reb_add(sim,p);
    p = reb_tools_orbit_to_particle(1, sim->particles[0],1e-2,
            1., //a
            0.2940098, //e
            0.80004, // i
            0.3020, // Omega
            0.1, // omega
            0.); //f
    reb_add(sim,p);

    double a, lambda, k, h, ix, iy;

    reb_integrate(sim,10.);
   
    p = sim->particles[1]; 
    reb_particle_to_pal(1.,p,sim->particles[0],&a,&lambda,&k,&h,&ix,&iy);
    struct reb_particle pp = reb_pal_to_particle(1.,sim->particles[0],1e-2,a,lambda,k,h,ix,iy);
    printf("%f\n",a);
    printf("%f %f %f %f %f %f\n",p.x,p.y,p.z, p.vx,p.vy,p.vz);
    printf("%f %f %f %f %f %f\n",pp.x,pp.y,pp.z, pp.vx,pp.vy,pp.vz);
    
    { 
    struct reb_particle pv = reb_vary_pal_lambda(1.,sim->particles[1],sim->particles[0]);
    double Delta = 1e-8;
    struct reb_particle ppp = reb_pal_to_particle(1.,sim->particles[0],1e-2,a,lambda+Delta,k,h,ix,iy);
    printf("\n");
    printf("%f %f %f %f %f %f\n",pv.x,pv.y,pv.z, pv.vx,pv.vy,pv.vz);
    printf("%f %f %f %f %f %f\n",(ppp.x-pp.x)/Delta,(ppp.y-pp.y)/Delta,(ppp.z-pp.z)/Delta, (ppp.vx-pp.vx)/Delta,(ppp.vy-pp.vy)/Delta,(ppp.vz-pp.vz)/Delta);
    }

    {
    struct reb_particle pv = reb_vary_pal_h(1.,sim->particles[1],sim->particles[0]);
    double Delta = 1e-8;
    struct reb_particle ppp = reb_pal_to_particle(1.,sim->particles[0],1e-2,a,lambda,k,h+Delta,ix,iy);
    printf("\n");
    printf("%f %f %f %f %f %f\n",pv.x,pv.y,pv.z, pv.vx,pv.vy,pv.vz);
    printf("%f %f %f %f %f %f\n",(ppp.x-pp.x)/Delta,(ppp.y-pp.y)/Delta,(ppp.z-pp.z)/Delta, (ppp.vx-pp.vx)/Delta,(ppp.vy-pp.vy)/Delta,(ppp.vz-pp.vz)/Delta);
    }

    {
    struct reb_particle pv = reb_vary_pal_k(1.,sim->particles[1],sim->particles[0]);
    double Delta = 1e-8;
    struct reb_particle ppp = reb_pal_to_particle(1.,sim->particles[0],1e-2,a,lambda,k+Delta,h,ix,iy);
    printf("\n");
    printf("%f %f %f %f %f %f\n",pv.x,pv.y,pv.z, pv.vx,pv.vy,pv.vz);
    printf("%f %f %f %f %f %f\n",(ppp.x-pp.x)/Delta,(ppp.y-pp.y)/Delta,(ppp.z-pp.z)/Delta, (ppp.vx-pp.vx)/Delta,(ppp.vy-pp.vy)/Delta,(ppp.vz-pp.vz)/Delta);
    }

    {
    struct reb_particle pv = reb_vary_pal_kk(1.,sim->particles[1],sim->particles[0]);
    double Delta = 1e-4;
    struct reb_particle ppp = reb_pal_to_particle(1.,sim->particles[0],1e-2,a,lambda,k+Delta,h,ix,iy);
    struct reb_particle ppm = reb_pal_to_particle(1.,sim->particles[0],1e-2,a,lambda,k-Delta,h,ix,iy);
    printf("\n");
    printf("%f %f %f %f %f %f\n",pv.x,pv.y,pv.z, pv.vx,pv.vy,pv.vz);
    printf("%f %f %f %f %f %f\n",(ppp.x-2.*pp.x+ppm.x)/(Delta*Delta),(ppp.y-2.*pp.y+ppm.y)/(Delta*Delta),(ppp.z-2.*pp.z+ppm.z)/(Delta*Delta), 
    	(ppp.vx-2.*pp.vx+ppm.vx)/(Delta*Delta),(ppp.vy-2.*pp.vy+ppm.vy)/(Delta*Delta),(ppp.vz-2.*pp.vz+ppm.vz)/(Delta*Delta));
    }

    {
    struct reb_particle pv = reb_vary_pal_hh(1.,sim->particles[1],sim->particles[0]);
    double Delta = 1e-4;
    struct reb_particle ppp = reb_pal_to_particle(1.,sim->particles[0],1e-2,a,lambda,k,h+Delta,ix,iy);
    struct reb_particle ppm = reb_pal_to_particle(1.,sim->particles[0],1e-2,a,lambda,k,h-Delta,ix,iy);
    printf("\n");
    printf("%f %f %f %f %f %f\n",pv.x,pv.y,pv.z, pv.vx,pv.vy,pv.vz);
    printf("%f %f %f %f %f %f\n",(ppp.x-2.*pp.x+ppm.x)/(Delta*Delta),(ppp.y-2.*pp.y+ppm.y)/(Delta*Delta),(ppp.z-2.*pp.z+ppm.z)/(Delta*Delta), 
    	(ppp.vx-2.*pp.vx+ppm.vx)/(Delta*Delta),(ppp.vy-2.*pp.vy+ppm.vy)/(Delta*Delta),(ppp.vz-2.*pp.vz+ppm.vz)/(Delta*Delta));
    }

    {
    struct reb_particle pv = reb_vary_pal_lambdalambda(1.,sim->particles[1],sim->particles[0]);
    double Delta = 1e-4;
    struct reb_particle ppp = reb_pal_to_particle(1.,sim->particles[0],1e-2,a,lambda+Delta,k,h,ix,iy);
    struct reb_particle ppm = reb_pal_to_particle(1.,sim->particles[0],1e-2,a,lambda-Delta,k,h,ix,iy);
    printf("\n");
    printf("%f %f %f %f %f %f\n",pv.x,pv.y,pv.z, pv.vx,pv.vy,pv.vz);
    printf("%f %f %f %f %f %f\n",(ppp.x-2.*pp.x+ppm.x)/(Delta*Delta),(ppp.y-2.*pp.y+ppm.y)/(Delta*Delta),(ppp.z-2.*pp.z+ppm.z)/(Delta*Delta), 
    	(ppp.vx-2.*pp.vx+ppm.vx)/(Delta*Delta),(ppp.vy-2.*pp.vy+ppm.vy)/(Delta*Delta),(ppp.vz-2.*pp.vz+ppm.vz)/(Delta*Delta));
    }

    {
    struct reb_particle pv = reb_vary_pal_klambda(1.,sim->particles[1],sim->particles[0]);
    double Delta = 1e-4;
    struct reb_particle pppp = reb_pal_to_particle(1.,sim->particles[0],1e-2,a,lambda+Delta,k+Delta,h,ix,iy);
    struct reb_particle ppmp = reb_pal_to_particle(1.,sim->particles[0],1e-2,a,lambda-Delta,k+Delta,h,ix,iy);
    struct reb_particle pppm = reb_pal_to_particle(1.,sim->particles[0],1e-2,a,lambda+Delta,k-Delta,h,ix,iy);
    struct reb_particle ppmm = reb_pal_to_particle(1.,sim->particles[0],1e-2,a,lambda-Delta,k-Delta,h,ix,iy);
    printf("\n");
    printf("%f %f %f %f %f %f\n",pv.x,pv.y,pv.z, pv.vx,pv.vy,pv.vz);
    printf("%f %f %f %f %f %f\n",(pppp.x+ppmm.x-pppm.x-ppmp.x)/(4.*Delta*Delta),(pppp.y+ppmm.y-pppm.y-ppmp.y)/(4.*Delta*Delta),(pppp.z+ppmm.z-pppm.z-ppmp.z)/(4.*Delta*Delta), 
    	(pppp.vx+ppmm.vx-pppm.vx-ppmp.vx)/(4.*Delta*Delta),(pppp.vy+ppmm.vy-pppm.vy-ppmp.vy)/(4.*Delta*Delta),(pppp.vz+ppmm.vz-pppm.vz-ppmp.vz)/(4.*Delta*Delta));
    }

    {
    struct reb_particle pv = reb_vary_pal_hlambda(1.,sim->particles[1],sim->particles[0]);
    double Delta = 1e-4;
    struct reb_particle pppp = reb_pal_to_particle(1.,sim->particles[0],1e-2,a,lambda+Delta,k,h+Delta,ix,iy);
    struct reb_particle ppmp = reb_pal_to_particle(1.,sim->particles[0],1e-2,a,lambda-Delta,k,h+Delta,ix,iy);
    struct reb_particle pppm = reb_pal_to_particle(1.,sim->particles[0],1e-2,a,lambda+Delta,k,h-Delta,ix,iy);
    struct reb_particle ppmm = reb_pal_to_particle(1.,sim->particles[0],1e-2,a,lambda-Delta,k,h-Delta,ix,iy);
    printf("\n");
    printf("%f %f %f %f %f %f\n",pv.x,pv.y,pv.z, pv.vx,pv.vy,pv.vz);
    printf("%f %f %f %f %f %f\n",(pppp.x+ppmm.x-pppm.x-ppmp.x)/(4.*Delta*Delta),(pppp.y+ppmm.y-pppm.y-ppmp.y)/(4.*Delta*Delta),(pppp.z+ppmm.z-pppm.z-ppmp.z)/(4.*Delta*Delta), 
    	(pppp.vx+ppmm.vx-pppm.vx-ppmp.vx)/(4.*Delta*Delta),(pppp.vy+ppmm.vy-pppm.vy-ppmp.vy)/(4.*Delta*Delta),(pppp.vz+ppmm.vz-pppm.vz-ppmp.vz)/(4.*Delta*Delta));
    }

    {
    struct reb_particle pv = reb_vary_pal_kh(1.,sim->particles[1],sim->particles[0]);
    double Delta = 1e-4;
    struct reb_particle pppp = reb_pal_to_particle(1.,sim->particles[0],1e-2,a,lambda,k+Delta,h+Delta,ix,iy);
    struct reb_particle ppmp = reb_pal_to_particle(1.,sim->particles[0],1e-2,a,lambda,k-Delta,h+Delta,ix,iy);
    struct reb_particle pppm = reb_pal_to_particle(1.,sim->particles[0],1e-2,a,lambda,k+Delta,h-Delta,ix,iy);
    struct reb_particle ppmm = reb_pal_to_particle(1.,sim->particles[0],1e-2,a,lambda,k-Delta,h-Delta,ix,iy);
    printf("\n");
    printf("%f %f %f %f %f %f\n",pv.x,pv.y,pv.z, pv.vx,pv.vy,pv.vz);
    printf("%f %f %f %f %f %f\n",(pppp.x+ppmm.x-pppm.x-ppmp.x)/(4.*Delta*Delta),(pppp.y+ppmm.y-pppm.y-ppmp.y)/(4.*Delta*Delta),(pppp.z+ppmm.z-pppm.z-ppmp.z)/(4.*Delta*Delta), 
    	(pppp.vx+ppmm.vx-pppm.vx-ppmp.vx)/(4.*Delta*Delta),(pppp.vy+ppmm.vy-pppm.vy-ppmp.vy)/(4.*Delta*Delta),(pppp.vz+ppmm.vz-pppm.vz-ppmp.vz)/(4.*Delta*Delta));
    }
}

