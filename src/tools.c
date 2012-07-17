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

void tools_move_to_center_of_momentum(){
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

void tools_init_plummer(int _N, double mlow, double rfrac, int quiet, double scale, double* shift) {
	struct particle* _particles = calloc(_N,sizeof(struct particle));
	double scalefactor = (scale < 0 ?  16.0 / (3.0 * M_PI)  : scale);
	double inv_scalefactor = 1.0 / scalefactor;
	double sqrt_scalefactor = sqrt( scalefactor );
	rfrac *= scalefactor;          /* from VIRIAL to STRUCTURAL units */
	double mfrac = rfrac*rfrac*rfrac / pow(1.0 + rfrac*rfrac, 1.5);

	// Setup particles
	for (int i = 0; i < _N; i++) {
		_particles[i].m = 1.0/ (double) _N;
		double radius = 0;
		if (quiet==0){
			radius = 1.0 / sqrt( pow (tools_uniform(mlow,mfrac), -2.0/3.0) - 1.0);
		} else if (quiet==1) {
			double m_min = (i * mfrac)/(double)_N;
			double m_max = ((i+1) * mfrac)/(double)_N;
			radius = 1.0 / sqrt( pow (tools_uniform(m_min,m_max), -2.0/3.0) - 1.0);
		} else if (quiet==2) {
			double m_med = ((i+0.5) * mfrac)/(double)_N;
			radius = 1.0 / sqrt( pow (m_med, -2.0/3.0) - 1.0);
		} 	
		double theta = acos(tools_uniform(-1.0, 1.0));
		double phi = tools_uniform(0.0, 2.*M_PI);
		_particles[i].x = radius * sin( theta ) * cos( phi );
		_particles[i].y = radius * sin( theta ) * sin( phi );
		_particles[i].z = radius * cos( theta );
		double x = 0.0;
		double y = 0.1;
		while (y > x*x*pow( 1.0 - x*x, 3.5)) {
			x = tools_uniform(0.0,1.0);
			y = tools_uniform(0.0,0.1);
		}
		double velocity = x * sqrt(2.0) * pow( 1.0 + radius*radius, -0.25);
		theta = acos(tools_uniform(-1.0, 1.0));
		phi = tools_uniform(0.0,2.*M_PI);
		_particles[i].vx = velocity * sin( theta ) * cos( phi );
		_particles[i].vy = velocity * sin( theta ) * sin( phi );
		_particles[i].vz = velocity * cos( theta );

	}

	// Scale model and calculate center of mass.
	double w_x  = 0, w_y  = 0, w_z  = 0;
	double w_vx = 0, w_vy = 0, w_vz = 0;
	for (int i = 0; i < _N; i++) {
		_particles[i].x  *= inv_scalefactor;
		_particles[i].y  *= inv_scalefactor;
		_particles[i].z  *= inv_scalefactor;
		_particles[i].vx *= sqrt_scalefactor;
		_particles[i].vy *= sqrt_scalefactor;
		_particles[i].vz *= sqrt_scalefactor;
		w_x  += _particles[i].x;
		w_y  += _particles[i].y;
		w_z  += _particles[i].z;
		w_vx += _particles[i].vx;
		w_vy += _particles[i].vy;
		w_vz += _particles[i].vz;
	}
	double w_ins = 1./(double)_N;
	w_x *= w_ins;
	w_y *= w_ins;
	w_z *= w_ins;
	w_x -= shift[0];
	w_y -= shift[1];
	w_z -= shift[2];
	w_vx *= w_ins;
	w_vy *= w_ins;
	w_vz *= w_ins;
	w_vx -= shift[3];
	w_vy -= shift[4];
	w_vz -= shift[5];
	for (int i = 0; i < _N; i++) {
		_particles[i].x -= w_x;
		_particles[i].y -= w_y;
		_particles[i].z -= w_z;
		_particles[i].vx -= w_vx;
		_particles[i].vy -= w_vy;
		_particles[i].vz -= w_vz;
		particles_add(_particles[i]);
	}
	free(_particles);
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
	o.l = ea -o.e*sin(ea)+o.omega;			// mean longitude
	if (o.e<=1.e-10){ 				//circular orbit
		o.omega=0.;
		o.f=0.; 				// f has no meaning
		o.l=0.;
	}

	return o;
}
