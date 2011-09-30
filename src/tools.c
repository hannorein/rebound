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

#define TINY 1.0e-12
struct orbit tools_p2orbit(struct particle p, double cmass){
	// Compute the angular momentum H, and thereby the inclination INC.
	double Omega, a, e, M, E=0, f=0, omega, inc; 	// orbital paramaers

	double u; int ialpha;			// internals
	
	double hx = p.y*p.vz - p.z*p.vy;
	double hy = p.z*p.vx - p.x*p.vz;
	double hz = p.x*p.vy - p.y*p.vx;
	double h2 = hx*hx + hy*hy +hz*hz;
	double h  = sqrt(h2);
	inc = acos(hz/h);

	// Compute longitude of ascending node CAPOM and the argument of latitude u.
	double fac = sqrt(hx*hx + hy*hy)/h;
	
	if(fac < TINY ){
	  	Omega = 0.;
	  	u = atan2(p.y,p.x);
	  	if(fabs(inc - M_PI) < 10.*TINY){
			u = -u;
	  	}
	}else{
	  	Omega = atan2(hx,-hy); 
	  	u = atan2 ( p.z/sin(inc) , p.x*cos(Omega) + p.y*sin(Omega));
	}

	while(Omega < 0.) Omega = Omega + 2.*M_PI;
	
	while(u < 0.) u = u + 2.*M_PI;
	

	//  Compute the radius R and velocity squared V2, and the dot product RDOTV, the energy per unit mass ENERGY .

	double r = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
	double v2 = p.vx*p.vx + p.vy*p.vy + p.vz*p.vz;
	double vdotr = p.x*p.vx + p.y*p.vy + p.z*p.vz;
	double energy = 0.5*v2 - G*cmass/r;

	//  Determine type of conic section and label it via IALPHA
	if(fabs(energy*r/(G*cmass)) < sqrt(TINY)){
		ialpha = 0;
	}else{
	   	if(energy < 0.) ialpha = -1;
	   	if(energy > 0.) ialpha = +1;
	}

	// Depending on the conic type, determine the remaining elements

	// ELLIPSE :
	if(ialpha == -1){
		a = -0.5*G*cmass/energy;
	  	fac = 1. - h2/(G*cmass*a);
		double cape, w;
		if (fac > TINY){
			e = sqrt ( fac );
             		double face =(a-r)/(a*e);

			//... Apr. 16/93 : watch for case where face is slightly outside unity
             		if ( face > 1.){
                		cape = 0.;
             		}else{
                		if ( face > -1.){
                   			cape = acos( face );
				}else{
                   			cape = M_PI;
				}
			}
			
            		if ( vdotr < 0. ) cape = 2.*M_PI - cape;
			double cw, sw;
	    		cw = (cos( cape) -e)/(1. - e*cos(cape));
	    		sw = sqrt(1. - e*e)*sin(cape)/(1. - e*cos(cape));
	    		w = atan2(sw,cw);
	    		while(w < 0.) w = w + 2.*M_PI;
	  	}else{
	    		e = 0.;
	    		w = u;
	    		cape = u;
		}
		f = w;
		E = cape;
	  	M = cape - e*sin (cape);
	  	omega = u - w;
	  	while(omega < 0.) omega = omega + 2.*M_PI;
	  	omega = omega - floor(omega/(2.*M_PI))*2.*M_PI;
	}
	// HYPERBOLA :
	if(ialpha == 1){
	  	a = 0.5*G*cmass/energy;
	  	fac = h2/(G*cmass*a);
		double w, capf;
          	if (fac > TINY){
 	    		e = sqrt ( 1. + fac );
	    		double tmpf = (a+r)/(a*e);
            		if (tmpf < 1.0){
              			 tmpf = 1.0;
			}
	    		capf = log(tmpf + sqrt(tmpf*tmpf -1.));
	    		if ( vdotr < 0. ) capf = - capf;
			double cw,sw;
	    		cw = (e - cosh(capf))/(e*cosh(capf) - 1. );
	    		sw = sqrt(e*e - 1.)*sinh(capf)/(e*cosh(capf) - 1. );
	    		w = atan2(sw,cw);
	    		if(w < 0.) w = w + 2.*M_PI;
	  	}else{
	// we only get here if a hyperbola is essentially a parabola so we calculate e and w accordingly to avoid singularities
	    		e = 1.;
	    		double tmpf = 0.5*h2/(G*cmass);
	    		w = acos(2.*tmpf/r -1.);
	    		if ( vdotr < 0.) w = 2.*M_PI - w;
	    		tmpf = (a+r)/(a*e);
	    		capf = log(tmpf + sqrt(tmpf*tmpf -1.));
	  	}

	  	M = e * sinh(capf) - capf;
	  	omega = u - w;
	  	if(omega < 0.) omega = omega + 2.*M_PI;
		omega = omega - floor(omega/(2.*M_PI))*2.*M_PI;
	}

	// PARABOLA : ( NOTE - in this case we use "a" to mean pericentric distance)

	if(ialpha == 0){
		double w;
		a =  0.5*h2/(G*cmass);
	  	e = 1.;
	  	w = acos(2.*a/r -1.);
	  	if ( vdotr < 0.) w = 2.*M_PI - w;
	  	double tmpf = tan(0.5 * w);
	  	M = tmpf* (1. + tmpf*tmpf/3.);
	  	omega = u - w;
	  	if(omega < 0.) omega = omega + 2.*M_PI;
	  	omega = omega - floor(omega/(2.*M_PI))*2.*M_PI; 	 
	}
	
	struct orbit o; 
	o.Omega 	= Omega;
	o.omega 	= omega;
	o.M		= M;
	o.r		= r;
	o.f		= f;
	o.E		= E;
	o.inc		= inc;
	o.h		= h;
	o.e		= e;
	o.a		= a;
	return o;
}
