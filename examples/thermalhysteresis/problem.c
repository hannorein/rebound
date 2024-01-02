/**
 * Thermal Hysteresis
 *
 * This example can be used as a starting point to reproduce 
 * the results of Larue, Latter, and Rein (2022).
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"


// Structure to store simulation parameters and output data
// Having this structure and setting it as r->extras allows 
// us to avoid having any global variables.
struct collisions_log {
    int Nslices;
    int Nsamples;   // Samples to avergae over (since last output)
    double lastsample; 
    double twarmup; 
    int isHot;      // Used during warmup
    double tau;
    double* plog;
    double* Elog;
    long* Nlog;
    double* T;
    double* qNL;
    double* qL;
    double* nuT;
    double* nuC;
};

// The "realistic" coefficient of restitution. See Eq (2). Contains code for warmup.
double eps_realistic(const struct reb_simulation* const r, double v, double x){
    struct collisions_log* log = (struct collisions_log* )r->extras; 
    v = fabs(v);
    double vc = 5.;
    double offset = 0;
    if (v<vc){
        return offset; 
    }
    double b = 1;
    double _x = (v-vc)/b;
    double epsmax = 0.923;
    double eps = offset+epsmax*1.624*_x/(1.+pow(_x,1.234));
    if (r->t<log->twarmup){
        if (!log->isHot){
            eps *= r->t/log->twarmup;
        }
    }
    if (eps>1) eps=1;
    if (eps<0) eps=0;
    return eps;
}

// A custom collision resolve routine. Needed because coefficient of resitution 
// depends on position and because of extra logging.
int collision_resolve(struct reb_simulation* const r, struct reb_collision c){
	struct reb_particle* const particles = r->particles;
	struct reb_particle p1 = particles[c.p1];
	struct reb_particle p2 = particles[c.p2];
	struct reb_vec6d gb = c.gb;
	double x21  = p1.x + gb.x  - p2.x; 
	double y21  = p1.y + gb.y  - p2.y; 
	double z21  = p1.z + gb.z  - p2.z; 
	double rp   = p1.r+p2.r;
	double oldvyouter;
    struct reb_particle old1 = p1; 
    struct reb_particle old2 = p2; 
	if (x21>0){
	 	oldvyouter = p1.vy;
	}else{
		oldvyouter = p2.vy;
	}
	if (rp*rp < x21*x21 + y21*y21 + z21*z21) return 0;
	double vx21 = p1.vx + gb.vx - p2.vx; 
	double vy21 = p1.vy + gb.vy - p2.vy; 
	double vz21 = p1.vz + gb.vz - p2.vz; 
	if (vx21*x21 + vy21*y21 + vz21*z21 >0) return 0; // not approaching
	// Bring the to balls in the xy plane.
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

	// Determine coefficient of restitution
	double eps = eps_realistic(r, vx21nn, (p1.x + gb.x  + p2.x)/2.);

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
	const double p2pf = p1.m/(p1.m+p2.m);
	particles[c.p2].vx -=	p2pf*dvx2n;
	particles[c.p2].vy -=	p2pf*dvy2nn;
	particles[c.p2].vz -=	p2pf*dvz2nn;
	particles[c.p2].last_collision = r->t;
	const double p1pf = p2.m/(p1.m+p2.m);
	particles[c.p1].vx +=	p1pf*dvx2n; 
	particles[c.p1].vy +=	p1pf*dvy2nn; 
	particles[c.p1].vz +=	p1pf*dvz2nn; 
	particles[c.p1].last_collision = r->t;
		
    struct reb_particle new1 = particles[c.p1];
    struct reb_particle new2 = particles[c.p2];
    new1.vy += 1.5*r->ri_sei.OMEGA*new1.x;
    new2.vy += 1.5*r->ri_sei.OMEGA*new2.x;
    old1.vy += 1.5*r->ri_sei.OMEGA*old1.x;
    old2.vy += 1.5*r->ri_sei.OMEGA*old2.x;

    // Logging
    struct collisions_log* log = (struct collisions_log* )r->extras; 
    double xmid = (p1.x+p2.x)/2.;
    int i = ((int)floor((xmid/r->boxsize.x+0.5)*log->Nslices))%log->Nslices;
    double E1 = 0.5*(old1.vx*old1.vx + old1.vy*old1.vy  + old1.vz*old1.vz); 
    double E2 = 0.5*(old2.vx*old2.vx + old2.vy*old2.vy  + old2.vz*old2.vz); 
    double E1p = 0.5*(new1.vx*new1.vx + new1.vy*new1.vy  + new1.vz*new1.vz); 
    double E2p = 0.5*(new2.vx*new2.vx + new2.vy*new2.vy  + new2.vz*new2.vz); 
    double E1s = E1 - 0.5*(E1+E2); 
    double E2s = E2 - 0.5*(E1+E2); 
    double E1sp = E1p - 0.5*(E1p+E2p); 
    double E2sp = E2p - 0.5*(E1p+E2p); 
    double dE1s = E1sp - E1s;
    double dE2s = E2sp - E2s;
	if (x21>0){
		log->Elog[i] += fabs(x21)*dE1s;
        log->plog[i] += -fabs(x21)*(oldvyouter-particles[c.p1].vy) * p1.m;
	}else{
		log->Elog[i] += fabs(x21)*dE2s;
        log->plog[i] += -fabs(x21)*(oldvyouter-particles[c.p2].vy) * p2.m;
	}
	log->Nlog[i]++;
    return 0;
}

double midplane_fillingfactor(const struct reb_simulation* const r){
    double area = 0.;
    for (int i=0;i<r->N;i++){
        struct reb_particle p = r->particles[i];
        double R2 = p.r*p.r-p.z*p.z;
        if (R2>0.){
            area += M_PI*R2;
        }
    }
    return area/(r->boxsize.x*r->boxsize.y);
}

struct reb_vec3d velocity_dispersion(const struct reb_simulation* const r, double xmin, double xmax){
    // Algorithm with reduced roundoff errors (see wikipedia)
    struct reb_vec3d A = {.x=0, .y=0, .z=0};
    struct reb_vec3d W = {.x=0, .y=0, .z=0};
    int Ncounted = 0;
    for (int i=0;i<r->N;i++){
        struct reb_vec3d Aim1 = A;
        struct reb_particle p = r->particles[i];
        if (p.x>xmin && p.x<xmax){
            A.x = A.x + (p.vx-A.x)/(double)(Ncounted+1);
            A.y = A.y + (p.vy+1.5*r->ri_sei.OMEGA*p.x-A.y)/(double)(Ncounted+1);
            A.z = A.z + (p.vz-A.z)/(double)(Ncounted+1);
            W.x = W.x + (p.vx-Aim1.x)*(p.vx-A.x);
            W.y = W.y + (p.vy+1.5*r->ri_sei.OMEGA*p.x-Aim1.y)*(p.vy+1.5*r->ri_sei.OMEGA*p.x-A.y);
            W.z = W.z + (p.vz-Aim1.z)*(p.vz-A.z);
            Ncounted++;
        }
    }
    W.x = sqrt(W.x/(double)Ncounted);
    W.y = sqrt(W.y/(double)Ncounted);
    W.z = sqrt(W.z/(double)Ncounted);

    // Return velocity dispersion in xx, yy, zz
    return W; 
}

void heartbeat(struct reb_simulation* const r){
    if (reb_simulation_output_check(r, 1e-3*2.*M_PI/r->ri_sei.OMEGA)){
        reb_simulation_output_timing(r, 0);
        //reb_output_append_velocity_dispersion("veldisp.txt");
    }
    return;
    struct collisions_log* log = (struct collisions_log* )r->extras; 

    // Calculate quantities for each slice
    int Nslices = log->Nslices;
    for (int i=0;i<Nslices;i++){
        double xmin = -r->boxsize.x/2. + r->boxsize.x * i /Nslices; 
        double xmax = xmin + r->boxsize.x /Nslices; 
        double sigma = 1./((xmax - xmin) * r->boxsize.y);
        struct reb_vec3d W = velocity_dispersion(r,xmin,xmax);
        double _T = 1./3.*(W.x*W.x+ W.y*W.y+ W.z*W.z)/(r->ri_sei.OMEGA*r->ri_sei.OMEGA);
        log->T[i] += _T;

        // qL
        double u_x = 0; 
        double u_y = 0; 
        double u_z = 0; 
        double Wxy = 0;
        int _N=0;
        for (int j=0;j<r->N;j++){
            struct reb_particle p = r->particles[j];
            if (p.x>xmin && p.x<xmax){
                double vx = p.vx;
                double vy = p.vy+1.5*r->ri_sei.OMEGA*p.x;
                double vz = p.vz;
                u_x += vx;
                u_y += vy;
                u_z += vz;
                Wxy += vx*vy;
                _N ++;
            }
        }
        sigma *= _N;
        u_x /= _N;
        u_y /= _N;
        u_z /= _N;
        Wxy /= _N;
        double _qL=0;
        for (int j=0;j<r->N;j++){
            struct reb_particle p = r->particles[j];
            if (p.x>xmin && p.x<xmax){
                double cx = p.vx - u_x;
                double cy = p.vy+1.5*r->ri_sei.OMEGA*p.x - u_y;
                double cz = p.vz - u_z;
                double c2 = cx*cx + cy*cy + cz*cz;
                _qL += 0.5*c2*cx;
            }
        }
        log->qL[i] += sigma * _qL / _N;
    
        // qNL
        double dt = r->t-log->lastsample;
        if (dt>0.5e-4*2.*M_PI/r->ri_sei.OMEGA){
            log->qNL[i] += sigma*log->Elog[i] /(dt*_N);
        }
        log->Elog[i] = 0;

        // nuT
        log->nuT[i] += 2./3. * Wxy / r->ri_sei.OMEGA;

        // nuC
        if (dt>0.5e-4*2.*M_PI/r->ri_sei.OMEGA){
            log->nuC[i] += 2.*log->plog[i] /(3.0*r->ri_sei.OMEGA*_N*dt);
        }
        log->plog[i] = 0;

    }
    log->lastsample = r->t;
    log->Nsamples ++;

    // Save output 10 times per orbit
    if (reb_simulation_output_check(r,0.1*2.*M_PI/r->ri_sei.OMEGA)){
        char buf[256];
        sprintf(buf,"out_tau%.1f_hot%d/out.txt",log->tau,log->isHot);
        FILE* f = fopen(buf,"a+");
        fprintf(f, "%5.3f\t",r->t/(2.*M_PI/r->ri_sei.OMEGA));     // 0
        double FF = midplane_fillingfactor(r);
        fprintf(f, "%5.7f\t",FF);                                 // 1
        
        for (int i=0;i<Nslices;i++){
            fprintf(f, "%5.3f\t", log->T[i]/log->Nsamples);       // 2  (c^2)
            fprintf(f, "%5.3f\t", log->qL[i]/log->Nsamples);      // 3
            fprintf(f, "%5.3f\t", log->qNL[i]/log->Nsamples);     // 4
            fprintf(f, "%5.3f\t", log->nuT[i]/log->Nsamples);     // 5
            fprintf(f, "%5.3f\t", log->nuC[i]/log->Nsamples);     // 6
            log->T[i] = 0;
            log->qL[i] = 0;
            log->qNL[i] = 0;
            log->nuT[i] = 0;
            log->nuC[i] = 0;
        }
        log->Nsamples = 0;
        fprintf(f, "\n");
        fclose(f);
    }

    // On screen update every 10 orbit
    if (reb_simulation_output_check(r,10.*2.*M_PI/r->ri_sei.OMEGA)){
        printf("tau = %.3f\t t = %.2f\n", log->tau, r->t/(2.*M_PI/r->ri_sei.OMEGA));
    }
}

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_simulation_create();
    
    // Starting the REBOUND visualization server. This
    // allows you to visualize the simulation by pointing 
    // your web browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);

    const double OMEGA    = 1;    
    r->opening_angle2     = .5;
    r->integrator         = REB_INTEGRATOR_SEI;
    r->boundary           = REB_BOUNDARY_SHEAR;
    r->gravity            = REB_GRAVITY_NONE;
    r->collision          = REB_COLLISION_LINETREE;
    r->collision_resolve  = collision_resolve;
    r->ri_sei.OMEGA       = OMEGA;
    r->ri_sei.OMEGAZ      = OMEGA;
    r->dt                 = 1e-2*2.*M_PI/OMEGA;
    r->heartbeat          = heartbeat;
    double boxsize        = 200; 
    reb_simulation_configure_box(r, boxsize, 8, 1, 1);
    r->N_ghost_x = 1;
    r->N_ghost_y = 1;
    r->N_ghost_z = 0;
    
    r->minimum_collision_velocity = OMEGA*0.001;  // small fraction of the shear accross a particle

    // Setup memory for logging.
    struct collisions_log* log= malloc(sizeof(struct collisions_log));

    // Read in command line arguments to overwrite defaults: tau, hot, tmax
    log->tau = 0.1;
    if (argc>1){
        log->tau = atof(argv[1]);
    }
    log->isHot = 0;
    if (argc>2){
        log->isHot = atoi(argv[2]);
    }
    double tmax = 200;  // in orbits
    if (argc>3){
        tmax = atof(argv[3]);
    }
    
    
    
    // Add all ring paricles
    double area = 0.;
    while (log->tau> area/(r->boxsize.x*r->boxsize.y)){
        struct reb_particle pt = {0};
        double fac = 1;
        pt.x     = reb_random_uniform(r, -r->boxsize.x/2.,r->boxsize.x/2.);
        if (log->isHot && pt.x > 0){
            fac = 20;
        }
        pt.y     = reb_random_uniform(r, -r->boxsize.y/2.,r->boxsize.y/2.);
        pt.vx    = fac*reb_random_normal(r, 1.)*OMEGA;
        pt.vy    = -1.5*pt.x*OMEGA+fac*reb_random_normal(r, 1.)*OMEGA;
        double a = fac*0.1*reb_random_normal(r, 1.)*OMEGA;
        double f = reb_random_uniform(r, 0,2.*M_PI);
        pt.z     = a*cos(f);
        pt.vz    = -a*sin(f);
        pt.r     = 1.;
        pt.m     = 1.;
        reb_simulation_add(r, pt);
        area += M_PI*pt.r*pt.r;
    }

    r->extras = log;
    log->Nslices = 1;
    log->lastsample = 0;
    log->Nsamples = 0;
    log->twarmup = 20 * 2.*M_PI/OMEGA;
    log->T = malloc(sizeof(double)*log->Nslices);
    log->qNL = malloc(sizeof(double)*log->Nslices);
    log->qL = malloc(sizeof(double)*log->Nslices);
    log->nuT = malloc(sizeof(double)*log->Nslices);
    log->nuC= malloc(sizeof(double)*log->Nslices);
    log->plog = malloc(sizeof(double)*log->Nslices);
    log->Elog = malloc(sizeof(double)*log->Nslices);
    log->Nlog = malloc(sizeof(long)*log->Nslices);
    for (int i=0;i<log->Nslices;i++){
        log->T[i] = 0;
        log->qNL[i] = 0;
        log->qL[i] = 0;
        log->nuT[i] = 0;
        log->nuC[i] = 0;
        log->plog[i] = 0;
        log->Elog[i] = 0;
        log->Nlog[i] = 0;
    }

    // Prepare output directories
    char buf[256];
    sprintf(buf,"rm -fr out_tau%.1f_hot%d",log->tau,log->isHot);
    system(buf);
    sprintf(buf,"mkdir out_tau%.1f_hot%d",log->tau,log->isHot);
    system(buf);

    // Integrate
    reb_simulation_integrate(r, tmax*2.*M_PI/OMEGA);
    
    // Final output
    sprintf(buf,"out_tau%.1f_hot%d/final.txt",log->tau,log->isHot);
    reb_simulation_output_ascii(r, buf);

    // Cleanup
    free(log->T);
    free(log->qNL);
    free(log->qL);
    free(log->nuT);
    free(log->nuC);
    free(log->plog);
    free(log->Elog);
    free(log->Nlog);
    free(log);
    reb_simulation_free(r);
    
}

