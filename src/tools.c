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
#define TWOPI 6.283185307179586

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

struct orbit tools_p2orbit(struct particle p, struct particle star){
	// Algorithm from: mercury6 
	
	struct orbit o;
	double hx, hy, hz, h2, h, v2, r, rv, s, tru, mu;
	double ci, to, temp, tmp2, bige, f, cf, ce;
	mu = G*(p.m+star.m);
	
	p.x -= star.x;
	p.y -= star.y;
	p.z -= star.z;
	p.vx -= star.vx;
	p.vy -= star.vy;
	p.vz -= star.vz;

	hx = p.y*p.vz - p.z*p.vy; // angular momentum vector
	hy = p.z*p.vx - p.x*p.vz;
	hz = p.x*p.vy - p.y*p.vx;
	h2 = hx*hx + hy*hy + hz*hz;
	
	v2 = p.vx*p.vx + p.vy*p.vy + p.vz*p.vz;
	rv = p.x*p.vx + p.y*p.vy + p.z*p.vz;
	
	r = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
	h = sqrt(h2);
	s = h2/mu;

	// Inclination and ascending node
	ci = hz/h;
	if (fabs(ci) < 1){
		o.inc = acos(ci);
		o.Omega = atan2(hx, -hy);
		if (o.Omega < 0) o.Omega = o.Omega + TWOPI;
	}
	else{
		if (ci > 0) o.inc = 0;
		if (ci < 0) o.inc = M_PI;
		o.Omega = 0;
	}

	// Eccentricity and pericentre distance
	temp = 1 + s*(v2/mu - 2/r);
	if (temp<0){
		o.e = 0;
	}
	else{
		o.e = sqrt(temp);
	}
	o.q = s/(1 + o.e);

	// True longitude
	if (hy != 0){
		to = -hx/hy;
		temp = (1 - ci)*to;
		tmp2 = to*to;
		tru = atan2((p.y*(1 + tmp2*ci) - p.x*temp), (p.x*(tmp2 + ci) - p.y*temp));
	}
	else{
		tru = atan2(p.y*ci, p.x);
	}
	if (ci<0) tru = tru + M_PI;

	if (o.e < 3.e-8){
		o.p = 0.0;
		o.l = tru; // mean longitude if "circular"
	}
	else{
		ce = (v2*r - mu)/(o.e*mu);
		
		// Mean anomaly for ellipse
		if (o.e < 1){
			if (fabs(ce) > 1) ce = copysign(1.0, ce);
			bige = acos(ce);
			if (rv<0) bige = TWOPI - bige;
			o.l = bige - o.e*sin(bige);
		}
		else{ // Mean anomaly for hyperbola
			if (ce < 1.0) ce = 1.0;
			bige = log(ce + sqrt(ce*ce - 1));
			if (rv < 0) bige = -bige;
			o.l = o.e*sinh(bige) - bige;
		}

		// Longitude of pericentre 
		cf = (s-r)/(o.e*r);
		if (fabs(cf)>1) cf = copysign(1.0, cf);
		f = acos(cf);
		if (rv<0) f = TWOPI - f;
		o.p = tru - f;
		o.p = fmod(o.p + TWOPI + TWOPI, TWOPI);
	}

	o.omega = o.p - o.Omega;
	if (o.l < 0 && o.e < 1) o.l = o.l + TWOPI;
	if (o.l > TWOPI && o.e < 1) o.l = fmod(o.l, TWOPI);
	
	return o;
}
