/**
 * @file 	tools.c
 * @brief 	Tools for creating distributions.
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
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "main.h"
#include "tools.h"


// Random number helper routines
double	tools_normaldistribution2_rsq;		/**< Used for speedup**/ 
double 	tools_normaldistribution2_v2;		/**< Used for speedup**/
int 	tools_normaldistribution2_ready = 0;	/**< Used for speedup**/

double tools_uniform(double min, double max){
	return ((double)rand())/((double)(RAND_MAX))*(max-min)+min;
}


double tools_powerlaw(double min, double max, double slope){
	double y = tools_uniform(0., 1.);
	return pow( (pow(max,slope+1.)-pow(min,slope+1.))*y+pow(min,slope+1.), 1./(slope+1.));
}

double tools_normal(double variance){
	if (tools_normaldistribution2_ready==1){
		tools_normaldistribution2_ready=0;
		return tools_normaldistribution2_v2*sqrt(-2.*log(tools_normaldistribution2_rsq)/tools_normaldistribution2_rsq*variance);
	}
	double v1,v2,rsq=1.;
	while(rsq>=1. || rsq<1.0e-12){
		v1=2.*((double)rand())/((double)(RAND_MAX))-1.0;
		v2=2.*((double)rand())/((double)(RAND_MAX))-1.0;
		rsq=v1*v1+v2*v2;
	}
	tools_normaldistribution2_ready = 1;
	tools_normaldistribution2_rsq   = rsq;
	tools_normaldistribution2_v2	  = v2;
	return 	v1*sqrt(-2.*log(rsq)/rsq*variance);
}

double tools_rayleigh(double sigma){
	double y = tools_uniform(0.,1.);
	return sigma*sqrt(-2*log(y));
}



/// Other helper routines
double tools_energy(void){
	double e_kin = 0.;
	double e_pot = 0.;
	for (int i=0;i<N-N_megno;i++){
		struct particle pi = particles[i];
		e_kin += 0.5 * pi.m * (pi.vx*pi.vx + pi.vy*pi.vy + pi.vz*pi.vz);
		for (int j=i+1;j<N-N_megno;j++){
			struct particle pj = particles[j];
			double dx = pi.x - pj.x;
			double dy = pi.y - pj.y;
			double dz = pi.z - pj.z;
			e_pot -= G*pj.m*pi.m/sqrt(dx*dx + dy*dy + dz*dz);
		}
	}
	return e_kin +e_pot;
}

void tools_move_to_center_of_momentum(void){
	double m = 0;
	double x = 0;
	double y = 0;
	double z = 0;
	double vx = 0;
	double vy = 0;
	double vz = 0;
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
		m  += p.m;
		x  += p.x*p.m;
		y  += p.y*p.m;
		z  += p.z*p.m;
		vx += p.vx*p.m;
		vy += p.vy*p.m;
		vz += p.vz*p.m;
	}
	x /= m;
	y /= m;
	z /= m;
	vx /= m;
	vy /= m;
	vz /= m;
	for (int i=0;i<N;i++){
		particles[i].x  -= x;
		particles[i].y  -= y;
		particles[i].z  -= z;
		particles[i].vx -= vx;
		particles[i].vy -= vy;
		particles[i].vz -= vz;
	}
}

struct particle tools_get_center_of_mass(struct particle p1, struct particle p2){
	p1.x   = p1.x*p1.m + p2.x*p2.m;		
	p1.y   = p1.y*p1.m + p2.y*p2.m;
	p1.z   = p1.z*p1.m + p2.z*p2.m;
	p1.vx  = p1.vx*p1.m + p2.vx*p2.m;
	p1.vy  = p1.vy*p1.m + p2.vy*p2.m;
	p1.vz  = p1.vz*p1.m + p2.vz*p2.m;
	p1.m  += p2.m;
	if (p1.m>0.){
		p1.x  /= p1.m;
		p1.y  /= p1.m;
		p1.z  /= p1.m;
		p1.vx /= p1.m;
		p1.vy /= p1.m;
		p1.vz /= p1.m;
	}
	return p1;
}

void tools_init_plummer(int _N, double M, double R) {
	// Algorithm from:	
	// http://adsabs.harvard.edu/abs/1974A%26A....37..183A
	
	double E = 3./64.*M_PI*M*M/R;
	for (int i=0;i<_N;i++){
		struct particle star;
		double r = pow(pow(tools_uniform(0,1),-2./3.)-1.,-1./2.);
		double x2 = tools_uniform(0,1);
		double x3 = tools_uniform(0,2.*M_PI);
		star.z = (1.-2.*x2)*r;
		star.x = sqrt(r*r-star.z*star.z)*cos(x3);
		star.y = sqrt(r*r-star.z*star.z)*sin(x3);
		double x5,g,q;
		do{
			x5 = tools_uniform(0.,1.);
			q = tools_uniform(0.,1.);
			g = q*q*pow(1.-q*q,7./2.);
		}while(0.1*x5>g);
		double ve = pow(2.,1./2.)*pow(1.+r*r,-1./4.);
		double v = q*ve;
		double x6 = tools_uniform(0.,1.);
		double x7 = tools_uniform(0.,2.*M_PI);
		star.vz = (1.-2.*x6)*v;
		star.vx = sqrt(v*v-star.vz*star.vz)*cos(x7);
		star.vy = sqrt(v*v-star.vz*star.vz)*sin(x7);
		
		star.x *= 3.*M_PI/64.*M*M/E;
		star.y *= 3.*M_PI/64.*M*M/E;
		star.z *= 3.*M_PI/64.*M*M/E;
		
		star.vx *= sqrt(E*64./3./M_PI/M);
		star.vy *= sqrt(E*64./3./M_PI/M);
		star.vz *= sqrt(E*64./3./M_PI/M);

		star.m = M/(double)_N;


		particles_add(star);
	}
}

struct particle tools_init_orbit2d(double M, double m, double a, double e, double omega, double f){
	struct particle p;
	p.m = m;
	double r = a*(1.-e*e)/(1.+e*cos(f));
	double n = sqrt(G*(m+M)/(a*a*a));
	double tx = r*cos(f);
	double ty = r*sin(f);
	p.z = 0;
	double tvx = -n*a/sqrt(1.-e*e)*sin(f);
	double tvy =  n*a/sqrt(1.-e*e)*(e+cos(f));
	p.vz = 0;
	p.ax = 0; p.ay = 0; p.az = 0;
	
	p.x  =  cos(omega)*tx  - sin(omega)*ty;
	p.y  =  sin(omega)*tx  + cos(omega)*ty;
	
	p.vx =  cos(omega)*tvx - sin(omega)*tvy;
	p.vy =  sin(omega)*tvx + cos(omega)*tvy;

	return p;
}

struct particle tools_init_orbit3d(double M, double m, double a, double e, double i, double Omega, double omega, double f){
	struct particle p;
	p.m = m;
	double r = a*(1-e*e)/(1 + e*cos(f));

	// Murray & Dermott Eq 2.122
	p.x  = r*(cos(Omega)*cos(omega+f) - sin(Omega)*sin(omega+f)*cos(i));
	p.y  = r*(sin(Omega)*cos(omega+f) + cos(Omega)*sin(omega+f)*cos(i));
	p.z  = r*sin(omega+f)*sin(i);

	double n = sqrt(G*(m+M)/(a*a*a));

	// Murray & Dermott Eq. 2.36 after applying the 3 rotation matrices from Sec. 2.8 to the velocities in the orbital plane
	p.vx = (n*a/sqrt(1-e*e))*((e+cos(f))*(-cos(i)*cos(omega)*sin(Omega) - cos(Omega)*sin(omega)) - sin(f)*(cos(omega)*cos(Omega) - cos(i)*sin(omega)*sin(Omega)));
	p.vy = (n*a/sqrt(1-e*e))*((e+cos(f))*(cos(i)*cos(omega)*cos(Omega) - sin(omega)*sin(Omega)) - sin(f)*(cos(omega)*sin(Omega) + cos(i)*cos(Omega)*sin(omega)));
	p.vz = (n*a/sqrt(1-e*e))*((e+cos(f))*cos(omega)*sin(i) - sin(f)*sin(i)*sin(omega));

	p.ax = 0; 	p.ay = 0; 	p.az = 0;

	return p;
}

#define TINY 1.0e-12
struct orbit tools_p2orbit(struct particle p, struct particle star){
	struct orbit o;
	double h0,h1,h2,e0,e1,e2,n0,n1,n,er,vr,mu,ea;
	mu = G*(p.m+star.m);
	p.x -= star.x;
	p.y -= star.y;
	p.z -= star.z;
	p.vx -= star.vx;
	p.vy -= star.vy;
	p.vz -= star.vz;
	h0 = (p.y*p.vz - p.z*p.vy); 			//angular momentum vector
	h1 = (p.z*p.vx - p.x*p.vz);
	h2 = (p.x*p.vy - p.y*p.vx);
	o.h = sqrt ( h0*h0 + h1*h1 + h2*h2 );		// abs value of angular moment 
	double v = sqrt ( p.vx*p.vx + p.vy*p.vy + p.vz*p.vz );
	o.r = sqrt ( p.x*p.x + p.y*p.y + p.z*p.z );
	vr = (p.x*p.vx + p.y*p.vy + p.z*p.vz)/o.r;
	e0 = 1./mu*( (v*v-mu/o.r)*p.x - o.r*vr*p.vx );
	e1 = 1./mu*( (v*v-mu/o.r)*p.y - o.r*vr*p.vy );
	e2 = 1./mu*( (v*v-mu/o.r)*p.z - o.r*vr*p.vz );
 	o.e = sqrt( e0*e0 + e1*e1 + e2*e2 );		// eccentricity
	o.a = -mu/( v*v - 2.*mu/o.r );			// semi major axis
	o.P = 2.*M_PI*sqrt( o.a*o.a*o.a/mu );		// period
	o.inc = acos( h2/o.h ) ;				// inclination (wrt xy-plane)   -  Note if pi/2 < i < pi then the orbit is retrograde
	n0 = -h1;					// vector of nodes lies in xy plane => no z component
	n1 =  h0;		
	n = sqrt( n0*n0 + n1*n1 );
	er = p.x*e0 + p.y*e1 + p.z*e2;
	if (n<=1.e-30||o.inc<=1.e-30){			// we are in the xy plane
		o.Omega=0.;
		if (e1>=0.) { o.omega=acos(e0/o.e); }else{ o.omega = 2.*M_PI-acos(e0/o.e); }
	}else{
		if (e2>=0.) { o.omega=acos(( n0*e0 + n1*e1 )/(n*o.e)); }else{ o.omega=2.*M_PI-acos(( n0*e0 + n1*e1 )/(n*o.e)); }// pericenter = 0 if pericenter = ascending node
		if (n1>=0.) { o.Omega = acos(n0/n); }else{  o.Omega=2.*M_PI-acos(n0/n);} 					// longitude of ascending node in xy plane, measured from x axis
	//	if (isnan(o.Omega)||isinf(o.Omega)) o.Omega=0.;
	}
	o.f = er/(o.e*o.r);
	ea = (1.-o.r/o.a)/o.e;
	if (o.f>1.||o.f<-1.){				// failsafe
		o.f = M_PI - M_PI * o.f;
		ea  = M_PI - M_PI * ea;
	}else{
		o.f = acos(o.f);			// true anomaly = 0 if planet at pericenter
		ea  = acos(ea);				// eccentric anomaly
	}
	
	if (vr<0.) { 
		o.f=2.*M_PI-o.f;	
		ea =2.*M_PI-ea;
	}
	o.l = ea -o.e*sin(ea)+o.omega + o.Omega;	// mean longitude
	if (o.e<=1.e-10){ 				//circular orbit
		o.omega=0.;
		o.f=0.; 				// f has no meaning
		o.l=0.;
	}

	return o;
}

/**************************
 * MEGNO Routines         */

double tools_megno_Ys;
double tools_megno_Yss;
double tools_megno_cov_Yt;	// covariance of <Y> and t
double tools_megno_var_t;  	// variance of t 
double tools_megno_mean_t; 	// mean of t
double tools_megno_mean_Y; 	// mean of Y
double tools_megno_delta0; 	// initial scale of delta (for one particle)
long   tools_megno_n; 		// number of covariance updates
void tools_megno_init(double delta){
	int _N_megno = N;
	tools_megno_Ys = 0.;
	tools_megno_Yss = 0.;
	tools_megno_cov_Yt = 0.;
	tools_megno_var_t = 0.;
	tools_megno_n = 0;
	tools_megno_mean_Y = 0;
	tools_megno_mean_t = 0;
	tools_megno_delta0 = delta;
        for (int i=0;i<_N_megno;i++){ 
                struct particle megno = {
			.m  = particles[i].m,
			.x  = tools_normal(1.),
			.y  = tools_normal(1.),
			.z  = tools_normal(1.),
			.vx = tools_normal(1.),
			.vy = tools_normal(1.),
			.vz = tools_normal(1.) };
		double deltad = delta/sqrt(megno.x*megno.x + megno.y*megno.y + megno.z*megno.z + megno.vx*megno.vx + megno.vy*megno.vy + megno.vz*megno.vz); // rescale
		megno.x *= deltad;
		megno.y *= deltad;
		megno.z *= deltad;
		megno.vx *= deltad;
		megno.vy *= deltad;
		megno.vz *= deltad;

                particles_add(megno);
        }
	N_megno = _N_megno;
}
double tools_megno(void){ // Returns the MEGNO <Y>
	if (t==0.) return 0.;
	return tools_megno_Yss/t;
}
double tools_lyapunov(void){ // Returns the largest Lyapunov characteristic number (LCN), or maximal Lyapunov exponent
	if (t==0.) return 0.;
	return tools_megno_cov_Yt/tools_megno_var_t;
}
double tools_megno_deltad_delta(void){
        double deltad = 0;
        double delta2 = 0;
        for (int i=N-N_megno;i<N;i++){
                deltad += particles[i].vx * particles[i].x; 
                deltad += particles[i].vy * particles[i].y; 
                deltad += particles[i].vz * particles[i].z; 
                deltad += particles[i].ax * particles[i].vx; 
                deltad += particles[i].ay * particles[i].vy; 
                deltad += particles[i].az * particles[i].vz; 
                delta2 += particles[i].x  * particles[i].x; 
                delta2 += particles[i].y  * particles[i].y;
                delta2 += particles[i].z  * particles[i].z;
                delta2 += particles[i].vx * particles[i].vx; 
                delta2 += particles[i].vy * particles[i].vy;
                delta2 += particles[i].vz * particles[i].vz;
        }
        return deltad/delta2;
}

// Automatically rescale megno particles if their size grew by more then 10 orders of magnitude.
//void tools_megno_rescale(void){
//	int rescale_flag = 0;
//	double delta2_limit = 1e100*tools_megno_delta0*tools_megno_delta0;
//        for (int i=N-N_megno;i<N;i++){
//		double delta2 = 0;
//		delta2 += particles[i].x  * particles[i].x; 
//		delta2 += particles[i].y  * particles[i].y;
//		delta2 += particles[i].z  * particles[i].z;
//		delta2 += particles[i].vx * particles[i].vx; 
//		delta2 += particles[i].vy * particles[i].vy;
//		delta2 += particles[i].vz * particles[i].vz;
//		if (delta2>delta2_limit){
//			rescale_flag = 1;
//		}
//        }
//	if (rescale_flag){
//		for (int i=N-N_megno;i<N;i++){
//			//particles[i].x  = tools_normal(1.);
//			//particles[i].y  = tools_normal(1.);
//			//particles[i].z  = tools_normal(1.);
//			//particles[i].vx = tools_normal(1.);
//			//particles[i].vy = tools_normal(1.);
//			//particles[i].vz = tools_normal(1.);
//			double deltad = tools_megno_delta0/sqrt(particles[i].x*particles[i].x + particles[i].y*particles[i].y + particles[i].z*particles[i].z + particles[i].vx*particles[i].vx + particles[i].vy*particles[i].vy + particles[i].vz*particles[i].vz); // rescale factor
//			particles[i].x *= deltad;
//			particles[i].y *= deltad;
//			particles[i].z *= deltad;
//			particles[i].vx *= deltad;
//			particles[i].vy *= deltad;
//			particles[i].vz *= deltad;
//		}
//	}
//}

void tools_megno_update(double dY){
	// Calculate running Y(t)
	tools_megno_Ys += dY;
	double Y = tools_megno_Ys/t;
	// Calculate averge <Y> 
	tools_megno_Yss += Y * dt;
	// Update covariance of (Y,t) and variance of t
	tools_megno_n++;
	double _d_t = t - tools_megno_mean_t;
	tools_megno_mean_t += _d_t/(double)tools_megno_n;
	double _d_Y = tools_megno() - tools_megno_mean_Y;
	tools_megno_mean_Y += _d_Y/(double)tools_megno_n;
	tools_megno_cov_Yt += ((double)tools_megno_n-1.)/(double)tools_megno_n 
					*(t-tools_megno_mean_t)
					*(tools_megno()-tools_megno_mean_Y);
	tools_megno_var_t  += ((double)tools_megno_n-1.)/(double)tools_megno_n 
					*(t-tools_megno_mean_t)
					*(t-tools_megno_mean_t);
//	// Check to see if we need to rescale
//	if (tools_megno_delta0!=0){
//		tools_megno_rescale();
//	}
}
