/**
 * @file    derivatives.c
 * @brief   Functions to calculate derivatives of Keplerian orbits.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section     LICENSE
 * Copyright (c) 2016 Hanno Rein, Dan Tamayp, Rejean Leblanc 
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
#include <sys/time.h>
#include "rebound.h"
#include "tools.h"
#include "derivatives.h"



struct reb_particle reb_derivatives_lambda(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};

    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);
    double dq_dlambda = -p/(1.-q);
    double dp_dlambda = q/(1.-q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dlambda = -1./(1.-q)*slp;
    double dslp_dlambda = 1./(1.-q)*clp;
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dxi_dlambda = a*(dclp_dlambda + dp_dlambda/(2.-l)*h);
    double deta_dlambda = a*(dslp_dlambda - dp_dlambda/(2.-l)*k);

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dlambda = deta_dlambda*ix-dxi_dlambda*iy;

    np.x = dxi_dlambda+0.5*iy*dW_dlambda;
    np.y = deta_dlambda-0.5*ix*dW_dlambda;
    np.z = 0.5*iz*dW_dlambda;

    double an = sqrt(G*(po.m+primary.m)/a);
    double ddxi_dlambda  = an/((1.-q)*(1.-q))*dq_dlambda*(-slp+q/(2.-l)*h)
    + an/(1.-q)*(-dslp_dlambda+dq_dlambda/(2.-l)*h);
    double ddeta_dlambda = an/((1.-q)*(1.-q))*dq_dlambda*(+clp-q/(2.-l)*k)
    + an/(1.-q)*(dclp_dlambda-dq_dlambda/(2.-l)*k);
    double ddW_dlambda = ddeta_dlambda*ix-ddxi_dlambda*iy;
    np.vx = ddxi_dlambda+0.5*iy*ddW_dlambda;
    np.vy = ddeta_dlambda-0.5*ix*ddW_dlambda;
    np.vz = 0.5*iz*ddW_dlambda;

    return np;
}

struct reb_particle reb_derivatives_h(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dh = -1./(1.-q)*(-slp*clp);
    double dslp_dh = -1./(1.-q)*(clp*clp);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dl_dh = 1./sqrt(1.-h*h-k*k)*h;
    double dp_dh = 1./(1.-q)*(-clp);
    double dxi_dh = a*(dclp_dh + dp_dh/(2.-l)*h + p/(2.-l) + p/((2.-l)*(2.-l))*dl_dh*h);
    double deta_dh = a*(dslp_dh - dp_dh/(2.-l)*k - p/((2.-l)*(2.-l))*k*dl_dh -1);

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dh = deta_dh*ix-dxi_dh*iy;

    np.x = dxi_dh+0.5*iy*dW_dh;
    np.y = deta_dh-0.5*ix*dW_dh;
    np.z = 0.5*iz*dW_dh;

    double dq_dh = 1./(1.-q)*(slp-h);

    double an = sqrt(G*(po.m+primary.m)/a);
    double ddxi_dh  = dq_dh*an/((1.-q)*(1.-q))*(-slp+q/(2.-l)*h)
                + an/(1.-q) * (-dslp_dh+dq_dh/(2.-l)*h+dl_dh*q/((2.-l)*(2.-l))*h+q/(2.-l));
    double ddeta_dh = dq_dh*an/((1.-q)*(1.-q))*(+clp-q/(2.-l)*k)
                + an/(1.-q) * (+dclp_dh-dq_dh/(2.-l)*k-dl_dh*q/((2.-l)*(2.-l))*k);
    double ddW_dh = ddeta_dh*ix-ddxi_dh*iy;

    np.vx = ddxi_dh+0.5*iy*ddW_dh;
    np.vy = ddeta_dh-0.5*ix*ddW_dh;
    np.vz = 0.5*iz*ddW_dh;

    return np;
}

struct reb_particle reb_derivatives_k(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dk = -1./(1.-q)*(slp*slp);
    double dslp_dk = -1./(1.-q)*(-slp*clp);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dl_dk = 1./sqrt(1.-h*h-k*k)*k;
    double dp_dk = 1./(1.-q)*(slp);
    double dxi_dk = a*(dclp_dk + dp_dk/(2.-l)*h + p/((2.-l)*(2.-l))*dl_dk*h -1);
    double deta_dk = a*(dslp_dk - dp_dk/(2.-l)*k - p/(2.-l) - p/((2.-l)*(2.-l))*dl_dk*k);

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dk = deta_dk*ix-dxi_dk*iy;

    np.x = dxi_dk+0.5*iy*dW_dk;
    np.y = deta_dk-0.5*ix*dW_dk;
    np.z = 0.5*iz*dW_dk;

    double dq_dk = 1./(1.-q)*(clp-k);

    double an = sqrt(G*(po.m+primary.m)/a);
    double ddxi_dk  = dq_dk*an/((1.-q)*(1.-q))*(-slp+q/(2.-l)*h)
                + an/(1.-q) * (-dslp_dk+dq_dk/(2.-l)*h+dl_dk*q/((2.-l)*(2.-l))*h);
    double ddeta_dk = dq_dk*an/((1.-q)*(1.-q))*(+clp-q/(2.-l)*k)
                + an/(1.-q) * (+dclp_dk-dq_dk/(2.-l)*k-dl_dk*q/((2.-l)*(2.-l))*k-q/(2.-l));
    double ddW_dk = ddeta_dk*ix-ddxi_dk*iy;

    np.vx = ddxi_dk+0.5*iy*ddW_dk;
    np.vy = ddeta_dk-0.5*ix*ddW_dk;
    np.vz = 0.5*iz*ddW_dk;

    return np;
}

struct reb_particle reb_derivatives_k_k(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dk = -1./(1.-q)*(slp*slp);
    double dslp_dk = -1./(1.-q)*(-slp*clp);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dl_dk = 1./sqrt(1.-h*h-k*k)*k;
    double dl_dkk = 1./sqrt(1.-h*h-k*k) + (k*k)/(sqrt(1.-h*h-k*k)*sqrt(1.-h*h-k*k)*sqrt(1.-h*h-k*k));
    double dp_dk = 1./(1.-q)*(slp);
    double dq_dk = 1./(1.-q)*(clp-k);
    double dp_dkk = dq_dk/((1.-q)*(1.-q))*(slp) + 1./(1.-q)*(dslp_dk);
    double dq_dkk = dq_dk/((1.-q)*(1.-q))*(clp-k) + 1./(1.-q)*(dclp_dk -1.);
    double dclp_dkk = -dq_dk/((1.-q)*(1.-q))*(slp*slp) -2./(1.-q)*slp*dslp_dk;
    double dslp_dkk = -dq_dk/((1.-q)*(1.-q))*(-slp*clp) -1./(1.-q)*-slp*dclp_dk -1./(1.-q)*-dslp_dk*clp;

    double dxi_dkk = a*(dclp_dkk + dp_dkk/(2.-l)*h + dl_dk*dp_dk/((2.-l)*(2.-l))*h + dp_dk/((2.-l)*(2.-l))*dl_dk*h + 2.*dl_dk*p/((2.-l)*(2.-l)*(2.-l))*dl_dk*h + p/((2.-l)*(2.-l))*dl_dkk*h);
    double deta_dkk = a*(dslp_dkk - dp_dkk/(2.-l)*k - dl_dk*dp_dk/((2.-l)*(2.-l))*k - dp_dk/(2.-l) - dp_dk/(2.-l) - dl_dk*p/((2.-l)*(2.-l)) 
                - dp_dk/((2.-l)*(2.-l))*dl_dk*k - 2.*dl_dk*p/((2.-l)*(2.-l)*(2.-l))*dl_dk*k - p/((2.-l)*(2.-l))*dl_dkk*k - p/((2.-l)*(2.-l))*dl_dk);

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dkk = deta_dkk*ix-dxi_dkk*iy;

    np.x = dxi_dkk+0.5*iy*dW_dkk;
    np.y = deta_dkk-0.5*ix*dW_dkk;
    np.z = 0.5*iz*dW_dkk;

    double an = sqrt(G*(po.m+primary.m)/a);
    double ddxi_dkk  = dq_dkk*an/((1.-q)*(1.-q))*(-slp+q/(2.-l)*h) + 2.*dq_dk*dq_dk*an/((1.-q)*(1.-q)*(1.-q))*(-slp+q/(2.-l)*h)
                + dq_dk*an/((1.-q)*(1.-q))*(-dslp_dk+dq_dk/(2.-l)*h+dl_dk*q/((2.-l)*(2.-l))*h)
                + dq_dk*an/((1.-q)*(1.-q))*(-dslp_dk+dq_dk/(2.-l)*h+dl_dk*q/((2.-l)*(2.-l))*h)
                + an/(1.-q)*(-dslp_dkk + dq_dkk/(2.-l)*h + dl_dk*dq_dk/((2.-l)*(2.-l))*h 
                + dl_dkk*q/((2.-l)*(2.-l))*h + dl_dk*dq_dk/((2.-l)*(2.-l))*h + 2.*dl_dk*dl_dk*q/((2.-l)*(2.-l)*(2.-l))*h );
    double ddeta_dkk = dq_dkk*an/((1.-q)*(1.-q))*(+clp-q/(2.-l)*k) + 2.*dq_dk*dq_dk*an/((1.-q)*(1.-q)*(1.-q))*(+clp-q/(2.-l)*k) 
                + dq_dk*an/((1.-q)*(1.-q))*(+dclp_dk-dq_dk/(2.-l)*k-dl_dk*q/((2.-l)*(2.-l))*k-q/(2.-l))
                + dq_dk*an/((1.-q)*(1.-q))*(+dclp_dk-dq_dk/(2.-l)*k-dl_dk*q/((2.-l)*(2.-l))*k-q/(2.-l))
                + an/(1.-q)*(+dclp_dkk - dq_dkk/(2.-l)*k - dq_dk*dl_dk/((2.-l)*(2.-l))*k - dq_dk/(2.-l) 
                - dl_dkk*q/((2.-l)*(2.-l))*k - dl_dk*dq_dk/((2.-l)*(2.-l))*k - 2.*dl_dk*dl_dk*q/((2.-l)*(2.-l)*(2.-l))*k - dl_dk*q/((2.-l)*(2.-l)) - dq_dk/(2.-l) - dl_dk*q/((2.-l)*(2.-l)) );
    double ddW_dkk = ddeta_dkk*ix-ddxi_dkk*iy;

    np.vx = ddxi_dkk+0.5*iy*ddW_dkk;
    np.vy = ddeta_dkk-0.5*ix*ddW_dkk;
    np.vz = 0.5*iz*ddW_dkk;

    return np;
}

struct reb_particle reb_derivatives_h_h(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dh = -1./(1.-q)*(-slp*clp);
    double dslp_dh = -1./(1.-q)*(clp*clp);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dl_dh = 1./sqrt(1.-h*h-k*k)*h;
    double dl_dhh = 1./sqrt(1.-h*h-k*k) + (h*h)/(sqrt(1.-h*h-k*k)*sqrt(1.-h*h-k*k)*sqrt(1.-h*h-k*k));
    double dp_dh = 1./(1.-q)*(-clp);
    double dq_dh = 1./(1.-q)*(slp-h);
    double dq_dhh = 1./((1.-q)*(1.-q))*dq_dh*(slp-h) + 1./(1.-q)*(dslp_dh-1);
    double dp_dhh = 1./((1.-q)*(1.-q))*dq_dh*(-clp) + 1./(1.-q)*(-dclp_dh);
    double dclp_dhh = -1./((1.-q)*(1.-q))*dq_dh*(-slp*clp) - 1./(1.-q)*(-dslp_dh*clp) - 1./(1.-q)*(-slp*dclp_dh);
    double dslp_dhh = -1./((1.-q)*(1.-q))*dq_dh*(clp*clp) - 2./(1.-q)*(clp*dclp_dh);

    double dxi_dhh = a*(dclp_dhh + (dp_dhh/(2.-l)*h + dl_dh*dp_dh/((2.-l)*(2.-l))*h + dp_dh/(2.-l)) + (dp_dh/(2.-l)+ dl_dh*p/((2.-l)*(2.-l)))
        + (dp_dh/((2.-l)*(2.-l))*dl_dh*h + 2.*p/((2.-l)*(2.-l)*(2.-l))*dl_dh*dl_dh*h + p/((2.-l)*(2.-l))*dl_dhh*h + p/((2.-l)*(2.-l))*dl_dh));
    double deta_dhh = a*(dslp_dhh + (-dp_dhh/(2.-l)*k - dl_dh*dp_dh/((2.-l)*(2.-l))*k) +(- dp_dh/((2.-l)*(2.-l))*k*dl_dh - 2.*p/((2.-l)*(2.-l)*(2.-l))*k*dl_dh*dl_dh- p/((2.-l)*(2.-l))*k*dl_dhh ));

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dhh = deta_dhh*ix-dxi_dhh*iy;

    np.x = dxi_dhh+0.5*iy*dW_dhh;
    np.y = deta_dhh-0.5*ix*dW_dhh;
    np.z = 0.5*iz*dW_dhh;

    double an = sqrt(G*(po.m+primary.m)/a);
    double ddxi_dhh  = dq_dhh*an/((1.-q)*(1.-q))*(-slp+q/(2.-l)*h) + 2.*dq_dh*dq_dh*an/((1.-q)*(1.-q)*(1.-q))*(-slp+q/(2.-l)*h) 
                + dq_dh*an/((1.-q)*(1.-q))*(-dslp_dh+dq_dh/(2.-l)*h+dl_dh*q/((2.-l)*(2.-l))*h+q/(2.-l))
                + dq_dh*an/((1.-q)*(1.-q))*(-dslp_dh+dq_dh/(2.-l)*h+dl_dh*q/((2.-l)*(2.-l))*h+q/(2.-l)) 
                + an/(1.-q)*(-dslp_dhh + (dq_dhh/(2.-l)*h+dl_dh*dq_dh/((2.-l)*(2.-l))*h+dq_dh/(2.-l)) 
                + (dl_dhh*q/((2.-l)*(2.-l))*h+dl_dh*dq_dh/((2.-l)*(2.-l))*h+2.*dl_dh*dl_dh*q/((2.-l)*(2.-l)*(2.-l))*h+dl_dh*q/((2.-l)*(2.-l))) + (dq_dh/(2.-l)+dl_dh*q/((2.-l)*(2.-l)))  );
    double ddeta_dhh = dq_dhh*an/((1.-q)*(1.-q))*(+clp-q/(2.-l)*k) + 2.*dq_dh*dq_dh*an/((1.-q)*(1.-q)*(1.-q))*(+clp-q/(2.-l)*k) 
                + dq_dh*an/((1.-q)*(1.-q))*(+dclp_dh-dq_dh/(2.-l)*k-dl_dh*q/((2.-l)*(2.-l))*k)
                + dq_dh*an/((1.-q)*(1.-q))*(+dclp_dh-dq_dh/(2.-l)*k-dl_dh*q/((2.-l)*(2.-l))*k)
                + an/(1.-q)*(+dclp_dhh - dq_dhh/(2.-l)*k - dl_dh*dq_dh/((2.-l)*(2.-l))*k 
                - dl_dhh*q/((2.-l)*(2.-l))*k - dl_dh*dq_dh/((2.-l)*(2.-l))*k - 2.*dl_dh*dl_dh*q/((2.-l)*(2.-l)*(2.-l))*k );

    double ddW_dhh = ddeta_dhh*ix-ddxi_dhh*iy;

    np.vx = ddxi_dhh+0.5*iy*ddW_dhh;
    np.vy = ddeta_dhh-0.5*ix*ddW_dhh;
    np.vz = 0.5*iz*ddW_dhh;

    return np;
}

struct reb_particle reb_derivatives_lambda_lambda(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};

    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);
    double dq_dlambda = -p/(1.-q);
    double dp_dlambda = q/(1.-q);
    double dq_dlambdalambda = -dp_dlambda/(1.-q) - p/((1.-q)*(1.-q))*dq_dlambda ;
    double dp_dlambdalambda = dq_dlambda/(1.-q) + q/((1.-q)*(1.-q))*dq_dlambda ;

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dlambda = -1./(1.-q)*slp;
    double dslp_dlambda = 1./(1.-q)*clp;
    double dclp_dlambdalambda = -1./((1.-q)*(1.-q))*dq_dlambda*slp -1./(1.-q)*dslp_dlambda;
    double dslp_dlambdalambda = 1./((1.-q)*(1.-q))*dq_dlambda*clp + 1./(1.-q)*dclp_dlambda;    
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dxi_dlambdalambda = a*(dclp_dlambdalambda + dp_dlambdalambda/(2.-l)*h);
    double deta_dlambdalambda = a*(dslp_dlambdalambda - dp_dlambdalambda/(2.-l)*k);

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dlambdalambda = deta_dlambdalambda*ix-dxi_dlambdalambda*iy;

    np.x = dxi_dlambdalambda+0.5*iy*dW_dlambdalambda;
    np.y = deta_dlambdalambda-0.5*ix*dW_dlambdalambda;
    np.z = 0.5*iz*dW_dlambdalambda;

    double an = sqrt(G*(po.m+primary.m)/a);
    double ddxi_dlambdalambda  = 2.*an/((1.-q)*(1.-q)*(1.-q))*dq_dlambda*dq_dlambda*(-slp+q/(2.-l)*h) 
                + an/((1.-q)*(1.-q))*dq_dlambdalambda*(-slp+q/(2.-l)*h) + an/((1.-q)*(1.-q))*dq_dlambda*(-dslp_dlambda+dq_dlambda/(2.-l)*h)
                + an/((1.-q)*(1.-q))*dq_dlambda*(-dslp_dlambda+dq_dlambda/(2.-l)*h) + an/(1.-q)*(-dslp_dlambdalambda+dq_dlambdalambda/(2.-l)*h);
    double ddeta_dlambdalambda = 2.*an/((1.-q)*(1.-q)*(1.-q))*dq_dlambda*dq_dlambda*(+clp-q/(2.-l)*k) 
                + an/((1.-q)*(1.-q))*dq_dlambdalambda*(+clp-q/(2.-l)*k) + an/((1.-q)*(1.-q))*dq_dlambda*(dclp_dlambda-dq_dlambda/(2.-l)*k)
                + an/((1.-q)*(1.-q))*dq_dlambda*(dclp_dlambda-dq_dlambda/(2.-l)*k) + an/(1.-q)*(dclp_dlambdalambda-dq_dlambdalambda/(2.-l)*k);

    double ddW_dlambdalambda = ddeta_dlambdalambda*ix-ddxi_dlambdalambda*iy;
    np.vx = ddxi_dlambdalambda+0.5*iy*ddW_dlambdalambda;
    np.vy = ddeta_dlambdalambda-0.5*ix*ddW_dlambdalambda;
    np.vz = 0.5*iz*ddW_dlambdalambda;

    return np;
}

struct reb_particle reb_derivatives_k_lambda(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dk = -1./(1.-q)*(slp*slp);
    double dslp_dk = -1./(1.-q)*(-slp*clp);
    double dclp_dlambda = -1./(1.-q)*slp;
    double dslp_dlambda = 1./(1.-q)*clp;
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dl_dk = 1./sqrt(1.-h*h-k*k)*k;
    double dp_dk = 1./(1.-q)*(slp);
    double dq_dk = 1./(1.-q)*(clp-k);
    double dq_dlambda = -p/(1.-q);
    double dp_dlambda = q/(1.-q);
    double dq_dklambda = -dp_dk/(1.-q) -p/((1.-q)*(1.-q))*dq_dk;
    double dp_dklambda = dq_dk/(1.-q) + q/((1.-q)*(1.-q))*dq_dk;
    double dclp_dklambda = -1./(1.-q)*dslp_dk -1./((1.-q)*(1.-q))*dq_dk*slp;
    double dslp_dklambda = 1./(1.-q)*dclp_dk + 1./((1.-q)*(1.-q))*dq_dk*clp;


    double dxi_dklambda = a*(dclp_dklambda + dp_dklambda/(2.-l)*h + dp_dlambda/((2.-l)*(2.-l))*dl_dk*h);
    double deta_dklambda = a*(dslp_dklambda - dp_dklambda/(2.-l)*k - dp_dlambda/(2.-l) - dp_dlambda/((2.-l)*(2.-l))*dl_dk*k);

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dklambda = deta_dklambda*ix-dxi_dklambda*iy;

    np.x = dxi_dklambda+0.5*iy*dW_dklambda;
    np.y = deta_dklambda-0.5*ix*dW_dklambda;
    np.z = 0.5*iz*dW_dklambda;

    double an = sqrt(G*(po.m+primary.m)/a);
    double ddxi_dklambda  = dq_dklambda*an/((1.-q)*(1.-q))*(-slp+q/(2.-l)*h) + 2.*dq_dk*dq_dlambda*an/((1.-q)*(1.-q)*(1.-q))*(-slp+q/(2.-l)*h) 
                + dq_dk*an/((1.-q)*(1.-q))*(-dslp_dlambda+dq_dlambda/(2.-l)*h)
                + dq_dlambda*an/((1.-q)*(1.-q))*(-dslp_dk+dq_dk/(2.-l)*h+dl_dk*q/((2.-l)*(2.-l))*h)
                + an/(1.-q)*(-dslp_dklambda+dq_dklambda/(2.-l)*h+dl_dk*dq_dlambda/((2.-l)*(2.-l))*h);
    double ddeta_dklambda = dq_dklambda*an/((1.-q)*(1.-q))*(+clp-q/(2.-l)*k) + 2.*dq_dk*dq_dlambda*an/((1.-q)*(1.-q)*(1.-q))*(+clp-q/(2.-l)*k) 
                + dq_dk*an/((1.-q)*(1.-q))*(+dclp_dlambda-dq_dlambda/(2.-l)*k)
                + dq_dlambda*an/((1.-q)*(1.-q))*(+dclp_dk-dq_dk/(2.-l)*k-dl_dk*q/((2.-l)*(2.-l))*k-q/(2.-l)) 
                + an/(1.-q)*(+dclp_dklambda-dq_dklambda/(2.-l)*k-dl_dk*dq_dlambda/((2.-l)*(2.-l))*k-dq_dlambda/(2.-l));
    double ddW_dklambda = ddeta_dklambda*ix-ddxi_dklambda*iy;

    np.vx = ddxi_dklambda+0.5*iy*ddW_dklambda;
    np.vy = ddeta_dklambda-0.5*ix*ddW_dklambda;
    np.vz = 0.5*iz*ddW_dklambda;

    return np;
}

struct reb_particle reb_derivatives_h_lambda(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dh = -1./(1.-q)*(-slp*clp);
    double dslp_dh = -1./(1.-q)*(clp*clp);
    double dclp_dlambda = -1./(1.-q)*slp;
    double dslp_dlambda = 1./(1.-q)*clp;
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dl_dh = 1./sqrt(1.-h*h-k*k)*h;
    double dp_dh = 1./(1.-q)*(-clp);
    double dq_dh = 1./(1.-q)*(slp-h);
    double dq_dlambda = -p/(1.-q);
    double dp_dlambda = q/(1.-q);
    double dq_dhlambda = -dp_dh/(1.-q) -p/((1.-q)*(1.-q))*dq_dh;
    double dp_dhlambda = dq_dh/(1.-q) +q/((1.-q)*(1.-q))*dq_dh;
    double dclp_dhlambda = -1./((1.-q)*(1.-q))*(-slp*clp)*dq_dlambda -1./(1.-q)*(-dslp_dlambda*clp) -1./(1.-q)*(-slp*dclp_dlambda);
    double dslp_dhlambda = -1./((1.-q)*(1.-q))*(clp*clp)*dq_dlambda -2./(1.-q)*(clp*dclp_dlambda);

    double dxi_dhlambda = a*(dclp_dhlambda + dp_dhlambda/(2.-l)*h + dp_dlambda/(2.-l) + dp_dlambda/((2.-l)*(2.-l))*dl_dh*h);
    double deta_dhlambda = a*(dslp_dhlambda - dp_dhlambda/(2.-l)*k - dp_dlambda/((2.-l)*(2.-l))*k*dl_dh);

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dhlambda = deta_dhlambda*ix-dxi_dhlambda*iy;

    np.x = dxi_dhlambda+0.5*iy*dW_dhlambda;
    np.y = deta_dhlambda-0.5*ix*dW_dhlambda;
    np.z = 0.5*iz*dW_dhlambda;

    double an = sqrt(G*(po.m+primary.m)/a);
    double ddxi_dhlambda  = dq_dhlambda*an/((1.-q)*(1.-q))*(-slp+q/(2.-l)*h) + 2.*dq_dlambda*dq_dh*an/((1.-q)*(1.-q)*(1.-q))*(-slp+q/(2.-l)*h) 
                + dq_dh*an/((1.-q)*(1.-q))*(-dslp_dlambda+dq_dlambda/(2.-l)*h)
                + dq_dlambda*an/((1.-q)*(1.-q))*(-dslp_dh+dq_dh/(2.-l)*h+dl_dh*q/((2.-l)*(2.-l))*h+q/(2.-l)) 
                + an/(1.-q)*(-dslp_dhlambda+dq_dhlambda/(2.-l)*h+dl_dh*dq_dlambda/((2.-l)*(2.-l))*h+dq_dlambda/(2.-l));
    double ddeta_dhlambda = dq_dhlambda*an/((1.-q)*(1.-q))*(+clp-q/(2.-l)*k) + 2.*dq_dh*dq_dlambda*an/((1.-q)*(1.-q)*(1.-q))*(+clp-q/(2.-l)*k) 
                + dq_dh*an/((1.-q)*(1.-q))*(+dclp_dlambda-dq_dlambda/(2.-l)*k)
                + dq_dlambda*an/((1.-q)*(1.-q))*(+dclp_dh-dq_dh/(2.-l)*k-dl_dh*q/((2.-l)*(2.-l))*k) 
                + an/(1.-q)*(+dclp_dhlambda-dq_dhlambda/(2.-l)*k-dl_dh*dq_dlambda/((2.-l)*(2.-l))*k);
    double ddW_dhlambda = ddeta_dhlambda*ix-ddxi_dhlambda*iy;

    np.vx = ddxi_dhlambda+0.5*iy*ddW_dhlambda;
    np.vy = ddeta_dhlambda-0.5*ix*ddW_dhlambda;
    np.vz = 0.5*iz*ddW_dhlambda;

    return np;
}

struct reb_particle reb_derivatives_k_h(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dk = -1./(1.-q)*(slp*slp);
    double dslp_dk = -1./(1.-q)*(-slp*clp);
    double dclp_dh = -1./(1.-q)*(-slp*clp);
    double dslp_dh = -1./(1.-q)*(clp*clp);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dl_dh = 1./sqrt(1.-h*h-k*k)*h;
    double dp_dh = 1./(1.-q)*(-clp);
    double dq_dh = 1./(1.-q)*(slp-h);
    double dl_dk = 1./sqrt(1.-h*h-k*k)*k;
    double dp_dk = 1./(1.-q)*(slp);
    double dq_dk = 1./(1.-q)*(clp-k);
    double dl_dkh = k*h/(sqrt(1.-h*h-k*k)*sqrt(1.-h*h-k*k)*sqrt(1.-h*h-k*k));
    double dp_dkh = 1./((1.-q)*(1.-q))*dq_dh*(slp) + 1./(1.-q)*(dslp_dh);
    double dq_dkh = 1./((1.-q)*(1.-q))*dq_dh*(clp-k) + 1./(1.-q)*(dclp_dh);
    double dclp_dkh = -1./((1.-q)*(1.-q))*dq_dh*(slp*slp) -2./(1.-q)*(slp*dslp_dh);
    double dslp_dkh = -1./((1.-q)*(1.-q))*dq_dh*(-slp*clp) -1./(1.-q)*(-dslp_dh*clp) -1./(1.-q)*(-slp*dclp_dh);

    double dxi_dkh = a*(dclp_dkh + dp_dkh/(2.-l)*h + dl_dh*dp_dk/((2.-l)*(2.-l))*h + dp_dk/(2.-l) 
                + dp_dh/((2.-l)*(2.-l))*dl_dk*h + 2.*p/((2.-l)*(2.-l)*(2.-l))*dl_dk*dl_dh*h + p/((2.-l)*(2.-l))*dl_dkh*h + p/((2.-l)*(2.-l))*dl_dk);
    double deta_dkh = a*(dslp_dkh - dp_dkh/(2.-l)*k - dl_dh*dp_dk/((2.-l)*(2.-l))*k - dp_dh/(2.-l)- dl_dh*p/((2.-l)*(2.-l)) 
                - dp_dh/((2.-l)*(2.-l))*dl_dk*k - p/((2.-l)*(2.-l))*dl_dkh*k - 2.*p/((2.-l)*(2.-l)*(2.-l))*dl_dk*dl_dh*k);

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dkh = deta_dkh*ix-dxi_dkh*iy;

    np.x = dxi_dkh+0.5*iy*dW_dkh;
    np.y = deta_dkh-0.5*ix*dW_dkh;
    np.z = 0.5*iz*dW_dkh;

    double an = sqrt(G*(po.m+primary.m)/a);
    double ddxi_dkh = dq_dkh*an/((1.-q)*(1.-q))*(-slp+q/(2.-l)*h) + 2.*dq_dh*dq_dk*an/((1.-q)*(1.-q)*(1.-q))*(-slp+q/(2.-l)*h) 
                + dq_dk*an/((1.-q)*(1.-q))*(-dslp_dh+dq_dh/(2.-l)*h + dl_dh*q/((2.-l)*(2.-l))*h + q/(2.-l))
                + dq_dh*an/((1.-q)*(1.-q))*(-dslp_dk+dq_dk/(2.-l)*h+dl_dk*q/((2.-l)*(2.-l))*h) 
                + an/(1.-q)*(-dslp_dkh+(dq_dkh/(2.-l)*h+dl_dh*dq_dk/((2.-l)*(2.-l))*h+dq_dk/(2.-l)) 
                + dl_dkh*q/((2.-l)*(2.-l))*h + dl_dk*dq_dh/((2.-l)*(2.-l))*h + 2.*dl_dh*dl_dk*q/((2.-l)*(2.-l)*(2.-l))*h + dl_dk*q/((2.-l)*(2.-l)));
    double ddeta_dkh = dq_dkh*an/((1.-q)*(1.-q))*(+clp-q/(2.-l)*k) + 2.*dq_dh*dq_dk*an/((1.-q)*(1.-q)*(1.-q))*(+clp-q/(2.-l)*k) 
                + dq_dk*an/((1.-q)*(1.-q))*(+dclp_dh-dq_dh/(2.-l)*k-dl_dh*q/((2.-l)*(2.-l))*k)
                + dq_dh*an/((1.-q)*(1.-q))*(+dclp_dk-dq_dk/(2.-l)*k-dl_dk*q/((2.-l)*(2.-l))*k-q/(2.-l))
                + an/(1.-q)*(+dclp_dkh-dq_dkh/(2.-l)*k-dl_dh*dq_dk/((2.-l)*(2.-l))*k 
                -dl_dkh*q/((2.-l)*(2.-l))*k -dl_dk*dq_dh/((2.-l)*(2.-l))*k -2.*dl_dk*dl_dh*q/((2.-l)*(2.-l)*(2.-l))*k -dq_dh/(2.-l)-dl_dh*q/((2.-l)*(2.-l)) );
    double ddW_dkh = ddeta_dkh*ix-ddxi_dkh*iy;

    np.vx = ddxi_dkh+0.5*iy*ddW_dkh;
    np.vy = ddeta_dkh-0.5*ix*ddW_dkh;
    np.vz = 0.5*iz*ddW_dkh;

    return np;
}

struct reb_particle reb_derivatives_a(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};

    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dxi_da = clp + p/(2.-l)*h -k;
    double deta_da = slp - p/(2.-l)*k -h;

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_da = deta_da*ix-dxi_da*iy;

    np.x = dxi_da+0.5*iy*dW_da;
    np.y = deta_da-0.5*ix*dW_da;
    np.z = 0.5*iz*dW_da;

    double dan_da = -0.5*sqrt(G*(po.m+primary.m)/(a*a*a));
    double ddxi_da  = dan_da/(1.-q)*(-slp+q/(2.-l)*h);
    double ddeta_da = dan_da/(1.-q)*(+clp-q/(2.-l)*k);

    double ddW_da = ddeta_da*ix-ddxi_da*iy;
    np.vx = ddxi_da+0.5*iy*ddW_da;
    np.vy = ddeta_da-0.5*ix*ddW_da;
    np.vz = 0.5*iz*ddW_da;

    return np;
}

struct reb_particle reb_derivatives_a_a(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};

    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dxi_daa = 0.0;//clp + p/(2.-l)*h -k;
    double deta_daa = 0.0;//slp - p/(2.-l)*k -h;

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_daa = deta_daa*ix-dxi_daa*iy;

    np.x = dxi_daa+0.5*iy*dW_daa;
    np.y = deta_daa-0.5*ix*dW_daa;
    np.z = 0.5*iz*dW_daa;

    double dan_daa = 0.75*sqrt(G*(po.m+primary.m)/(a*a*a*a*a));
    double ddxi_daa  = dan_daa/(1.-q)*(-slp+q/(2.-l)*h);
    double ddeta_daa = dan_daa/(1.-q)*(+clp-q/(2.-l)*k);

    double ddW_daa = ddeta_daa*ix-ddxi_daa*iy;
    np.vx = ddxi_daa+0.5*iy*ddW_daa;
    np.vy = ddeta_daa-0.5*ix*ddW_daa;
    np.vz = 0.5*iz*ddW_daa;

    return np;
}

struct reb_particle reb_derivatives_ix(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double xi = a*(clp + p/(2.-l)*h -k);
    double eta = a*(slp - p/(2.-l)*k -h);

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double diz_dix = -ix/sqrt(fabs(4.-ix*ix-iy*iy));
    double W = eta*ix-xi*iy;
    double dW_dix = eta;

    np.x = 0.5*iy*dW_dix;
    np.y = -0.5*W-0.5*ix*dW_dix;
    np.z = 0.5*diz_dix*W + 0.5*iz*dW_dix;

    double an = sqrt(G*(po.m+primary.m)/a);
    double dxi  = an/(1.-q)*(-slp+q/(2.-l)*h);
    double deta = an/(1.-q)*(+clp-q/(2.-l)*k);
    double dW = deta*ix-dxi*iy;
    double ddW_dix = deta;

    np.vx = 0.5*iy*ddW_dix;
    np.vy = -0.5*dW-0.5*ix*ddW_dix;
    np.vz = 0.5*diz_dix*dW + 0.5*iz*ddW_dix;

    return np;
}

struct reb_particle reb_derivatives_ix_ix(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double xi = a*(clp + p/(2.-l)*h -k);
    double eta = a*(slp - p/(2.-l)*k -h);

    double iz = sqrt(4.-ix*ix-iy*iy);
    double diz_dix = -ix/iz;
    double diz_dixix = -1./iz - ix*ix/(iz*iz*iz);
    double W = eta*ix-xi*iy;
    double dW_dix = eta;
    double dW_dixix = 0.0;

    np.x = 0.5*iy*dW_dixix;
    np.y = -dW_dix-0.5*ix*dW_dixix;
    np.z = 0.5*diz_dixix*W+diz_dix*dW_dix;

    double an = sqrt(G*(po.m+primary.m)/a);
    double dxi  = an/(1.-q)*(-slp+q/(2.-l)*h);
    double deta = an/(1.-q)*(+clp-q/(2.-l)*k);
    double dW = deta*ix-dxi*iy;
    double ddW_dix = deta;
    double ddW_dixix = 0.0;

    np.vx = 0.5*iy*ddW_dixix;
    np.vy = -ddW_dix-0.5*ix*ddW_dixix;
    np.vz = 0.5*diz_dixix*dW+diz_dix*ddW_dix;

    return np;
}

struct reb_particle reb_derivatives_iy(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double xi = a*(clp + p/(2.-l)*h -k);
    double eta = a*(slp - p/(2.-l)*k -h);

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double diz_diy = -iy/sqrt(fabs(4.-ix*ix-iy*iy));
    double W = eta*ix-xi*iy;
    double dW_diy = -xi;

    np.x = 0.5*W+0.5*iy*dW_diy;
    np.y = -0.5*ix*dW_diy;
    np.z = 0.5*diz_diy*W + 0.5*iz*dW_diy;

    double an = sqrt(G*(po.m+primary.m)/a);
    double dxi  = an/(1.-q)*(-slp+q/(2.-l)*h);
    double deta = an/(1.-q)*(+clp-q/(2.-l)*k);
    double dW = deta*ix-dxi*iy;
    double ddW_diy = -dxi;

    np.vx = 0.5*dW+0.5*iy*ddW_diy;
    np.vy = -0.5*ix*ddW_diy;
    np.vz = 0.5*diz_diy*dW + 0.5*iz*ddW_diy;

    return np;
}

struct reb_particle reb_derivatives_iy_iy(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double xi = a*(clp + p/(2.-l)*h -k);
    double eta = a*(slp - p/(2.-l)*k -h);

    //double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double diz_diy = -iy/sqrt(fabs(4.-ix*ix-iy*iy));
    double diz_diyiy = -1./sqrt(fabs(4.-ix*ix-iy*iy)) -iy*iy/( sqrt(fabs(4.-ix*ix-iy*iy))*sqrt(fabs(4.-ix*ix-iy*iy))*sqrt(fabs(4.-ix*ix-iy*iy)) );
    double W = eta*ix-xi*iy;
    double dW_diy = -xi;
    double dW_diyiy = 0.0;

    np.x = dW_diy+0.5*iy*dW_diyiy;
    np.y = -0.5*ix*dW_diyiy;
    np.z = 0.5*diz_diyiy*W + diz_diy*dW_diy;

    double an = sqrt(G*(po.m+primary.m)/a);
    double dxi  = an/(1.-q)*(-slp+q/(2.-l)*h);
    double deta = an/(1.-q)*(+clp-q/(2.-l)*k);
    double dW = deta*ix-dxi*iy;
    double ddW_diy = -dxi;
    double ddW_diyiy = 0.0;

    np.vx = ddW_diy-0.5*iy*ddW_diyiy;
    np.vy = -0.5*ix*ddW_diyiy;
    np.vz = 0.5*diz_diyiy*dW + diz_diy*ddW_diy;

    return np;
}

struct reb_particle reb_derivatives_k_ix(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dk = -1./(1.-q)*(slp*slp);
    double dslp_dk = -1./(1.-q)*(-slp*clp);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dl_dk = 1./sqrt(1.-h*h-k*k)*k;
    double dp_dk = 1./(1.-q)*(slp);
    double dxi_dk = a*(dclp_dk + dp_dk/(2.-l)*h + p/((2.-l)*(2.-l))*dl_dk*h -1);
    double deta_dk = a*(dslp_dk - dp_dk/(2.-l)*k - p/(2.-l) - p/((2.-l)*(2.-l))*dl_dk*k);

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double diz_dix = -ix/sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dk = deta_dk*ix-dxi_dk*iy;
    double dW_dkix = deta_dk;

    np.x = 0.5*iy*dW_dkix;
    np.y = -0.5*dW_dk-0.5*ix*dW_dkix;
    np.z = 0.5*diz_dix*dW_dk+0.5*iz*dW_dkix;

    double dq_dk = 1./(1.-q)*(clp-k);

    double an = sqrt(G*(po.m+primary.m)/a);
    double ddxi_dk  = dq_dk*an/((1.-q)*(1.-q))*(-slp+q/(2.-l)*h)
                + an/(1.-q) * (-dslp_dk+dq_dk/(2.-l)*h+dl_dk*q/((2.-l)*(2.-l))*h);
    double ddeta_dk = dq_dk*an/((1.-q)*(1.-q))*(+clp-q/(2.-l)*k)
                + an/(1.-q) * (+dclp_dk-dq_dk/(2.-l)*k-dl_dk*q/((2.-l)*(2.-l))*k-q/(2.-l));
    double ddW_dk = ddeta_dk*ix-ddxi_dk*iy;
    double ddW_dkix = ddeta_dk;

    np.vx = 0.5*iy*ddW_dkix;
    np.vy = -0.5*ddW_dk-0.5*ix*ddW_dkix;
    np.vz = 0.5*diz_dix*ddW_dk+0.5*iz*ddW_dkix;

    return np;
}

struct reb_particle reb_derivatives_h_ix(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dh = -1./(1.-q)*(-slp*clp);
    double dslp_dh = -1./(1.-q)*(clp*clp);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dl_dh = 1./sqrt(1.-h*h-k*k)*h;
    double dp_dh = 1./(1.-q)*(-clp);
    double dxi_dh = a*(dclp_dh + dp_dh/(2.-l)*h + p/(2.-l) + p/((2.-l)*(2.-l))*dl_dh*h);
    double deta_dh = a*(dslp_dh - dp_dh/(2.-l)*k - p/((2.-l)*(2.-l))*k*dl_dh -1);

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double diz_dix = -ix/sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dh = deta_dh*ix-dxi_dh*iy;
    double dW_dhix = deta_dh;

    np.x = 0.5*iy*dW_dhix;
    np.y = -0.5*dW_dh-0.5*ix*dW_dhix;
    np.z = 0.5*diz_dix*dW_dh+0.5*iz*dW_dhix;

    double dq_dh = 1./(1.-q)*(slp-h);

    double an = sqrt(G*(po.m+primary.m)/a);
    double ddxi_dh  = dq_dh*an/((1.-q)*(1.-q))*(-slp+q/(2.-l)*h)
                + an/(1.-q) * (-dslp_dh+dq_dh/(2.-l)*h+dl_dh*q/((2.-l)*(2.-l))*h+q/(2.-l));
    double ddeta_dh = dq_dh*an/((1.-q)*(1.-q))*(+clp-q/(2.-l)*k)
                + an/(1.-q) * (+dclp_dh-dq_dh/(2.-l)*k-dl_dh*q/((2.-l)*(2.-l))*k);
    double ddW_dh = ddeta_dh*ix-ddxi_dh*iy;
    double ddW_dhix = ddeta_dh;

    np.vx = 0.5*iy*ddW_dhix;
    np.vy = -0.5*ddW_dh-0.5*ix*ddW_dhix;
    np.vz = 0.5*diz_dix*ddW_dh+0.5*iz*ddW_dhix;

    return np;
}

struct reb_particle reb_derivatives_lambda_ix(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};

    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);
    double dq_dlambda = -p/(1.-q);
    double dp_dlambda = q/(1.-q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dlambda = -1./(1.-q)*slp;
    double dslp_dlambda = 1./(1.-q)*clp;
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dxi_dlambda = a*(dclp_dlambda + dp_dlambda/(2.-l)*h);
    double deta_dlambda = a*(dslp_dlambda - dp_dlambda/(2.-l)*k);

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double diz_dix = -ix/sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dlambda = deta_dlambda*ix-dxi_dlambda*iy;
    double dW_dlambdaix = deta_dlambda;

    np.x = 0.5*iy*dW_dlambdaix;
    np.y = -0.5*dW_dlambda-0.5*ix*dW_dlambdaix;
    np.z = 0.5*diz_dix*dW_dlambda+0.5*iz*dW_dlambdaix;

    double an = sqrt(G*(po.m+primary.m)/a);
    double ddxi_dlambda  = an/((1.-q)*(1.-q))*dq_dlambda*(-slp+q/(2.-l)*h)
    + an/(1.-q)*(-dslp_dlambda+dq_dlambda/(2.-l)*h);
    double ddeta_dlambda = an/((1.-q)*(1.-q))*dq_dlambda*(+clp-q/(2.-l)*k)
    + an/(1.-q)*(dclp_dlambda-dq_dlambda/(2.-l)*k);
    double ddW_dlambda = ddeta_dlambda*ix-ddxi_dlambda*iy;
    double ddW_dlambdaix = ddeta_dlambda;
    np.vx = 0.5*iy*ddW_dlambdaix;
    np.vy = -0.5*ddW_dlambda-0.5*ix*ddW_dlambdaix;
    np.vz = 0.5*diz_dix*ddW_dlambda+0.5*iz*ddW_dlambdaix;

    return np;
}

struct reb_particle reb_derivatives_lambda_iy(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};

    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);
    double dq_dlambda = -p/(1.-q);
    double dp_dlambda = q/(1.-q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dlambda = -1./(1.-q)*slp;
    double dslp_dlambda = 1./(1.-q)*clp;
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dxi_dlambda = a*(dclp_dlambda + dp_dlambda/(2.-l)*h);
    double deta_dlambda = a*(dslp_dlambda - dp_dlambda/(2.-l)*k);

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double diz_diy = -iy/sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dlambda = deta_dlambda*ix-dxi_dlambda*iy;
    double dW_dlambdaiy = -dxi_dlambda;
    np.x = 0.5*dW_dlambda+0.5*iy*dW_dlambdaiy;
    np.y = -0.5*ix*dW_dlambdaiy;
    np.z = 0.5*diz_diy*dW_dlambda+0.5*iz*dW_dlambdaiy;

    double an = sqrt(G*(po.m+primary.m)/a);
    double ddxi_dlambda  = an/((1.-q)*(1.-q))*dq_dlambda*(-slp+q/(2.-l)*h)
    + an/(1.-q)*(-dslp_dlambda+dq_dlambda/(2.-l)*h);
    double ddeta_dlambda = an/((1.-q)*(1.-q))*dq_dlambda*(+clp-q/(2.-l)*k)
    + an/(1.-q)*(dclp_dlambda-dq_dlambda/(2.-l)*k);
    double ddW_dlambda = ddeta_dlambda*ix-ddxi_dlambda*iy;
    double ddW_dlambdaiy = -ddxi_dlambda;
    np.vx = 0.5*ddW_dlambda+0.5*iy*ddW_dlambdaiy;
    np.vy = -0.5*ix*ddW_dlambdaiy;
    np.vz = 0.5*diz_diy*ddW_dlambda+0.5*iz*ddW_dlambdaiy;

    return np;
}

struct reb_particle reb_derivatives_h_iy(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dh = -1./(1.-q)*(-slp*clp);
    double dslp_dh = -1./(1.-q)*(clp*clp);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dl_dh = 1./sqrt(1.-h*h-k*k)*h;
    double dp_dh = 1./(1.-q)*(-clp);
    double dxi_dh = a*(dclp_dh + dp_dh/(2.-l)*h + p/(2.-l) + p/((2.-l)*(2.-l))*dl_dh*h);
    double deta_dh = a*(dslp_dh - dp_dh/(2.-l)*k - p/((2.-l)*(2.-l))*k*dl_dh -1);

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double diz_diy = -iy/sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dh = deta_dh*ix-dxi_dh*iy;
    double dW_dhiy = -dxi_dh;
    np.x = 0.5*dW_dh+0.5*iy*dW_dhiy;
    np.y = -0.5*ix*dW_dhiy;
    np.z = 0.5*diz_diy*dW_dh+0.5*iz*dW_dhiy;

    double dq_dh = 1./(1.-q)*(slp-h);

    double an = sqrt(G*(po.m+primary.m)/a);
    double ddxi_dh  = dq_dh*an/((1.-q)*(1.-q))*(-slp+q/(2.-l)*h)
                + an/(1.-q) * (-dslp_dh+dq_dh/(2.-l)*h+dl_dh*q/((2.-l)*(2.-l))*h+q/(2.-l));
    double ddeta_dh = dq_dh*an/((1.-q)*(1.-q))*(+clp-q/(2.-l)*k)
                + an/(1.-q) * (+dclp_dh-dq_dh/(2.-l)*k-dl_dh*q/((2.-l)*(2.-l))*k);
    double ddW_dh = ddeta_dh*ix-ddxi_dh*iy;
    double ddW_dhiy = -ddxi_dh;
    np.vx = 0.5*ddW_dh+0.5*iy*ddW_dhiy;
    np.vy = -0.5*ix*ddW_dhiy;
    np.vz = 0.5*diz_diy*ddW_dh+0.5*iz*ddW_dhiy;

    return np;
}

struct reb_particle reb_derivatives_k_iy(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dk = -1./(1.-q)*(slp*slp);
    double dslp_dk = -1./(1.-q)*(-slp*clp);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dl_dk = 1./sqrt(1.-h*h-k*k)*k;
    double dp_dk = 1./(1.-q)*(slp);
    double dxi_dk = a*(dclp_dk + dp_dk/(2.-l)*h + p/((2.-l)*(2.-l))*dl_dk*h -1);
    double deta_dk = a*(dslp_dk - dp_dk/(2.-l)*k - p/(2.-l) - p/((2.-l)*(2.-l))*dl_dk*k);

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double diz_diy = -iy/sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dk = deta_dk*ix-dxi_dk*iy;
    double dW_dkiy = -dxi_dk;
    np.x = 0.5*dW_dk+0.5*iy*dW_dkiy;
    np.y = -0.5*ix*dW_dkiy;
    np.z = 0.5*diz_diy*dW_dk+0.5*iz*dW_dkiy;

    double dq_dk = 1./(1.-q)*(clp-k);

    double an = sqrt(G*(po.m+primary.m)/a);
    double ddxi_dk  = dq_dk*an/((1.-q)*(1.-q))*(-slp+q/(2.-l)*h)
                + an/(1.-q) * (-dslp_dk+dq_dk/(2.-l)*h+dl_dk*q/((2.-l)*(2.-l))*h);
    double ddeta_dk = dq_dk*an/((1.-q)*(1.-q))*(+clp-q/(2.-l)*k)
                + an/(1.-q) * (+dclp_dk-dq_dk/(2.-l)*k-dl_dk*q/((2.-l)*(2.-l))*k-q/(2.-l));
    double ddW_dk = ddeta_dk*ix-ddxi_dk*iy;
    double ddW_dkiy = -ddxi_dk;
    np.vx = 0.5*ddW_dk+0.5*iy*ddW_dkiy;
    np.vy = -0.5*ix*ddW_dkiy;
    np.vz = 0.5*diz_diy*ddW_dk+0.5*iz*ddW_dkiy;

    return np;
}

struct reb_particle reb_derivatives_ix_iy(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double xi = a*(clp + p/(2.-l)*h -k);
    double eta = a*(slp - p/(2.-l)*k -h);

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double diz_dix = -ix/sqrt(fabs(4.-ix*ix-iy*iy));
    double diz_diy = -iy/sqrt(fabs(4.-ix*ix-iy*iy));
    double diz_dixiy = -ix*iy/(sqrt(fabs(4.-ix*ix-iy*iy))*sqrt(fabs(4.-ix*ix-iy*iy))*sqrt(fabs(4.-ix*ix-iy*iy)));
    double W = eta*ix-xi*iy;
    double dW_dix = eta;
    double dW_diy = -xi;
    double dW_dixiy = 0.0;
    np.x = 0.5*dW_dix+0.5*iy*dW_dixiy;
    np.y = -0.5*dW_diy-0.5*ix*dW_dixiy;
    np.z = 0.5*diz_dixiy*W+0.5*diz_dix*dW_diy + 0.5*diz_diy*dW_dix+0.5*iz*dW_dixiy;

    double an = sqrt(G*(po.m+primary.m)/a);
    double dxi  = an/(1.-q)*(-slp+q/(2.-l)*h);
    double deta = an/(1.-q)*(+clp-q/(2.-l)*k);
    double dW = deta*ix-dxi*iy;
    double ddW_dix = deta;
    double ddW_diy = -dxi;
    double ddW_dixiy = 0.0;

    np.vx = 0.5*ddW_dix+0.5*iy*ddW_dixiy;
    np.vy = -0.5*ddW_diy-0.5*ix*ddW_dixiy;
    np.vz = 0.5*diz_dixiy*dW+0.5*diz_dix*ddW_diy + 0.5*diz_diy*ddW_dix+0.5*iz*ddW_dixiy;

    return np;
}

struct reb_particle reb_derivatives_a_ix(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double deta_da = (slp - p/(2.-l)*k -h);
    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double diz_dix = -ix/sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_daix = deta_da;
    double dxi_da = clp + p/(2.-l)*h -k;
    double dW_da = deta_da*ix-dxi_da*iy;

    np.x = 0.5*iy*dW_daix;
    np.y = -0.5*dW_da-0.5*ix*dW_daix;
    np.z = 0.5*diz_dix*dW_da + 0.5*iz*dW_daix;

    double dan_da = -0.5*sqrt(G*(po.m+primary.m)/(a*a*a));
    double ddeta_da = dan_da/(1.-q)*(+clp-q/(2.-l)*k);
    double ddW_daix = ddeta_da;
    double ddxi_da  = dan_da/(1.-q)*(-slp+q/(2.-l)*h);
    double ddW_da = ddeta_da*ix-ddxi_da*iy;

    np.vx = 0.5*iy*ddW_daix;
    np.vy = -0.5*ddW_da-0.5*ix*ddW_daix;
    np.vz = 0.5*diz_dix*ddW_da + 0.5*iz*ddW_daix;

    return np;
}

struct reb_particle reb_derivatives_a_iy(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    
    double l = 1.-sqrt(1.-h*h-k*k);

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double diz_diy = -iy/sqrt(fabs(4.-ix*ix-iy*iy));
    double dxi_da = clp + p/(2.-l)*h -k;
    double deta_da = slp - p/(2.-l)*k -h;
    double dW_da = deta_da*ix-dxi_da*iy;
    double dW_daiy = -dxi_da;
    np.x = 0.5*dW_da+0.5*iy*dW_daiy;
    np.y = -0.5*ix*dW_daiy;
    np.z = 0.5*diz_diy*dW_da + 0.5*iz*dW_daiy;

    double dan_da = -0.5*sqrt(G*(po.m+primary.m)/(a*a*a));
    double ddxi_da  = dan_da/(1.-q)*(-slp+q/(2.-l)*h);
    double ddW_daiy = -ddxi_da;
    double ddeta_da = dan_da/(1.-q)*(+clp-q/(2.-l)*k);
    double ddW_da = ddeta_da*ix-ddxi_da*iy;

    np.vx = 0.5*ddW_da+0.5*iy*ddW_daiy;
    np.vy = -0.5*ix*ddW_daiy;
    np.vz = 0.5*diz_diy*ddW_da + 0.5*iz*ddW_daiy;

    return np;
}


struct reb_particle reb_derivatives_a_lambda(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};

    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);
    double dq_dlambda = -p/(1.-q);
    double dp_dlambda = q/(1.-q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dlambda = -1./(1.-q)*slp;
    double dslp_dlambda = 1./(1.-q)*clp;
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dxi_dalambda = (dclp_dlambda + dp_dlambda/(2.-l)*h);
    double deta_dalambda = (dslp_dlambda - dp_dlambda/(2.-l)*k);

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dalambda = deta_dalambda*ix-dxi_dalambda*iy;

    np.x = dxi_dalambda+0.5*iy*dW_dalambda;
    np.y = deta_dalambda-0.5*ix*dW_dalambda;
    np.z = 0.5*iz*dW_dalambda;

    double dan_da = -0.5*sqrt(G*(po.m+primary.m)/(a*a*a));
    double ddxi_dalambda  = dan_da/((1.-q)*(1.-q))*dq_dlambda*(-slp+q/(2.-l)*h)
    + dan_da/(1.-q)*(-dslp_dlambda+dq_dlambda/(2.-l)*h);
    double ddeta_dalambda = dan_da/((1.-q)*(1.-q))*dq_dlambda*(+clp-q/(2.-l)*k)
    + dan_da/(1.-q)*(dclp_dlambda-dq_dlambda/(2.-l)*k);
    double ddW_dalambda = ddeta_dalambda*ix-ddxi_dalambda*iy;
    np.vx = ddxi_dalambda+0.5*iy*ddW_dalambda;
    np.vy = ddeta_dalambda-0.5*ix*ddW_dalambda;
    np.vz = 0.5*iz*ddW_dalambda;

    return np;
}

struct reb_particle reb_derivatives_a_h(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dh = -1./(1.-q)*(-slp*clp);
    double dslp_dh = -1./(1.-q)*(clp*clp);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dl_dh = 1./sqrt(1.-h*h-k*k)*h;
    double dp_dh = 1./(1.-q)*(-clp);
    double dxi_dah = (dclp_dh + dp_dh/(2.-l)*h + p/(2.-l) + p/((2.-l)*(2.-l))*dl_dh*h);
    double deta_dah = (dslp_dh - dp_dh/(2.-l)*k - p/((2.-l)*(2.-l))*k*dl_dh -1);

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dah = deta_dah*ix-dxi_dah*iy;

    np.x = dxi_dah+0.5*iy*dW_dah;
    np.y = deta_dah-0.5*ix*dW_dah;
    np.z = 0.5*iz*dW_dah;

    double dq_dh = 1./(1.-q)*(slp-h);

    double dan_da = -0.5*sqrt(G*(po.m+primary.m)/(a*a*a));
    double ddxi_dah  = dq_dh*dan_da/((1.-q)*(1.-q))*(-slp+q/(2.-l)*h)
                + dan_da/(1.-q) * (-dslp_dh+dq_dh/(2.-l)*h+dl_dh*q/((2.-l)*(2.-l))*h+q/(2.-l));
    double ddeta_dah = dq_dh*dan_da/((1.-q)*(1.-q))*(+clp-q/(2.-l)*k)
                + dan_da/(1.-q) * (+dclp_dh-dq_dh/(2.-l)*k-dl_dh*q/((2.-l)*(2.-l))*k);
    double ddW_dah = ddeta_dah*ix-ddxi_dah*iy;

    np.vx = ddxi_dah+0.5*iy*ddW_dah;
    np.vy = ddeta_dah-0.5*ix*ddW_dah;
    np.vz = 0.5*iz*ddW_dah;

    return np;
}

struct reb_particle reb_derivatives_a_k(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dk = -1./(1.-q)*(slp*slp);
    double dslp_dk = -1./(1.-q)*(-slp*clp);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dl_dk = 1./sqrt(1.-h*h-k*k)*k;
    double dp_dk = 1./(1.-q)*(slp);
    double dxi_dak = (dclp_dk + dp_dk/(2.-l)*h + p/((2.-l)*(2.-l))*dl_dk*h -1);
    double deta_dak = (dslp_dk - dp_dk/(2.-l)*k - p/(2.-l) - p/((2.-l)*(2.-l))*dl_dk*k);

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double dW_dak = deta_dak*ix-dxi_dak*iy;

    np.x = dxi_dak+0.5*iy*dW_dak;
    np.y = deta_dak-0.5*ix*dW_dak;
    np.z = 0.5*iz*dW_dak;

    double dq_dk = 1./(1.-q)*(clp-k);

    double dan_da = -0.5*sqrt(G*(po.m+primary.m)/(a*a*a));
    double ddxi_dak  = dq_dk*dan_da/((1.-q)*(1.-q))*(-slp+q/(2.-l)*h)
                + dan_da/(1.-q) * (-dslp_dk+dq_dk/(2.-l)*h+dl_dk*q/((2.-l)*(2.-l))*h);
    double ddeta_dak = dq_dk*dan_da/((1.-q)*(1.-q))*(+clp-q/(2.-l)*k)
                + dan_da/(1.-q) * (+dclp_dk-dq_dk/(2.-l)*k-dl_dk*q/((2.-l)*(2.-l))*k-q/(2.-l));
    double ddW_dak = ddeta_dak*ix-ddxi_dak*iy;

    np.vx = ddxi_dak+0.5*iy*ddW_dak;
    np.vy = ddeta_dak-0.5*ix*ddW_dak;
    np.vz = 0.5*iz*ddW_dak;

    return np;
}

struct reb_particle reb_derivatives_m(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    np.m = 1.;
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    np.x = 0.0;
    np.y = 0.0;
    np.z = 0.0;

    double dan_dm = 0.5*sqrt(G/(a*(po.m+primary.m)));
    double ddxi_dm  = dan_dm/(1.-q)*(-slp+q/(2.-l)*h);
    double ddeta_dm = dan_dm/(1.-q)*(+clp-q/(2.-l)*k);

    double ddW_dm = ddeta_dm*ix-ddxi_dm*iy;
    np.vx = ddxi_dm+0.5*iy*ddW_dm;
    np.vy = ddeta_dm-0.5*ix*ddW_dm;
    np.vz = 0.5*iz*ddW_dm;

    return np;
}

struct reb_particle reb_derivatives_m_a(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double l = 1.-sqrt(1.-h*h-k*k);

    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    np.x = 0.0;
    np.y = 0.0;
    np.z = 0.0;

    double dan_dma = -0.5*0.5*sqrt(G/(a*a*a*(po.m+primary.m)));
    double ddxi_dma  = dan_dma/(1.-q)*(-slp+q/(2.-l)*h);
    double ddeta_dma = dan_dma/(1.-q)*(+clp-q/(2.-l)*k);

    double ddW_dma = ddeta_dma*ix-ddxi_dma*iy;
    np.vx = ddxi_dma+0.5*iy*ddW_dma;
    np.vy = ddeta_dma-0.5*ix*ddW_dma;
    np.vz = 0.5*iz*ddW_dma;

    return np;
}

struct reb_particle reb_derivatives_m_lambda(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);
    double dq_dlambda = -p/(1.-q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dlambda = -1./(1.-q)*slp;
    double dslp_dlambda = 1./(1.-q)*clp;
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double iz = sqrt(fabs(4.-ix*ix-iy*iy));

    np.x = 0.0;
    np.y = 0.0;
    np.z = 0.0;

    double dan_dm = 0.5*sqrt(G/(a*(po.m+primary.m)));
    double ddxi_dmlambda  = dan_dm/((1.-q)*(1.-q))*dq_dlambda*(-slp+q/(2.-l)*h)
    + dan_dm/(1.-q)*(-dslp_dlambda+dq_dlambda/(2.-l)*h);
    double ddeta_dmlambda = dan_dm/((1.-q)*(1.-q))*dq_dlambda*(+clp-q/(2.-l)*k)
    + dan_dm/(1.-q)*(dclp_dlambda-dq_dlambda/(2.-l)*k);
    double ddW_dmlambda = ddeta_dmlambda*ix-ddxi_dmlambda*iy;
    np.vx = ddxi_dmlambda+0.5*iy*ddW_dmlambda;
    np.vy = ddeta_dmlambda-0.5*ix*ddW_dmlambda;
    np.vz = 0.5*iz*ddW_dmlambda;

    return np;
}

struct reb_particle reb_derivatives_m_h(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);
    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dh = -1./(1.-q)*(-slp*clp);
    double dslp_dh = -1./(1.-q)*(clp*clp);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dl_dh = 1./sqrt(1.-h*h-k*k)*h;
    double iz = sqrt(fabs(4.-ix*ix-iy*iy));

    np.x = 0.0;
    np.y = 0.0;
    np.z = 0.0;

    double dq_dh = 1./(1.-q)*(slp-h);

    double dan_dm = 0.5*sqrt(G/(a*(po.m+primary.m)));
    double ddxi_dmh  = dq_dh*dan_dm/((1.-q)*(1.-q))*(-slp+q/(2.-l)*h)
                + dan_dm/(1.-q) * (-dslp_dh+dq_dh/(2.-l)*h+dl_dh*q/((2.-l)*(2.-l))*h+q/(2.-l));
    double ddeta_dmh = dq_dh*dan_dm/((1.-q)*(1.-q))*(+clp-q/(2.-l)*k)
                + dan_dm/(1.-q) * (+dclp_dh-dq_dh/(2.-l)*k-dl_dh*q/((2.-l)*(2.-l))*k);
    double ddW_dmh = ddeta_dmh*ix-ddxi_dmh*iy;

    np.vx = ddxi_dmh+0.5*iy*ddW_dmh;
    np.vy = ddeta_dmh-0.5*ix*ddW_dmh;
    np.vz = 0.5*iz*ddW_dmh;

    return np;
}

struct reb_particle reb_derivatives_m_k(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);
    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    double dclp_dk = -1./(1.-q)*(slp*slp);
    double dslp_dk = -1./(1.-q)*(-slp*clp);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double dl_dk = 1./sqrt(1.-h*h-k*k)*k;
    double iz = sqrt(fabs(4.-ix*ix-iy*iy));

    np.x = 0.0;
    np.y = 0.0;
    np.z = 0.0;

    double dq_dk = 1./(1.-q)*(clp-k);

    double dan_dm = 0.5*sqrt(G/(a*(po.m+primary.m)));
    double ddxi_dmk  = dq_dk*dan_dm/((1.-q)*(1.-q))*(-slp+q/(2.-l)*h)
                + dan_dm/(1.-q) * (-dslp_dk+dq_dk/(2.-l)*h+dl_dk*q/((2.-l)*(2.-l))*h);
    double ddeta_dmk = dq_dk*dan_dm/((1.-q)*(1.-q))*(+clp-q/(2.-l)*k)
                + dan_dm/(1.-q) * (+dclp_dk-dq_dk/(2.-l)*k-dl_dk*q/((2.-l)*(2.-l))*k-q/(2.-l));
    double ddW_dmk = ddeta_dmk*ix-ddxi_dmk*iy;

    np.vx = ddxi_dmk+0.5*iy*ddW_dmk;
    np.vy = ddeta_dmk-0.5*ix*ddW_dmk;
    np.vz = 0.5*iz*ddW_dmk;

    return np;
}

struct reb_particle reb_derivatives_m_ix(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);
    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double diz_dix = -ix/sqrt(fabs(4.-ix*ix-iy*iy));

    np.x = 0.0;
    np.y = 0.0;
    np.z = 0.0;

    double dan_dm = 0.5*sqrt(G/(a*(po.m+primary.m)));
    double ddxi_dm  = dan_dm/(1.-q)*(-slp+q/(2.-l)*h);
    double ddeta_dm = dan_dm/(1.-q)*(+clp-q/(2.-l)*k);
    double ddW_dm = ddeta_dm*ix-ddxi_dm*iy;
    double ddW_dmix = ddeta_dm;

    np.vx = 0.5*iy*ddW_dmix;
    np.vy = -0.5*ddW_dm-0.5*ix*ddW_dmix;
    np.vz = 0.5*diz_dix*ddW_dm + 0.5*iz*ddW_dmix;

    return np;
}

struct reb_particle reb_derivatives_m_iy(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);
    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    double diz_diy = -iy/sqrt(fabs(4.-ix*ix-iy*iy));

    np.x = 0.0;
    np.y = 0.0;
    np.z = 0.0;

    double dan_dm = 0.5*sqrt(G/(a*(po.m+primary.m)));
    double ddxi_dm  = dan_dm/(1.-q)*(-slp+q/(2.-l)*h);
    double ddeta_dm = dan_dm/(1.-q)*(+clp-q/(2.-l)*k);
    double ddW_dm = ddeta_dm*ix-ddxi_dm*iy;
    double ddW_dmiy = -ddxi_dm;

    np.vx = 0.5*ddW_dm+0.5*iy*ddW_dmiy;
    np.vy = -0.5*ix*ddW_dmiy;
    np.vz = 0.5*diz_diy*ddW_dm + 0.5*iz*ddW_dmiy;

    return np;
}

struct reb_particle reb_derivatives_m_m(double G, struct reb_particle primary, struct reb_particle po){
    double a, lambda, k, h, ix, iy;
    reb_tools_particle_to_pal(G, po, primary, &a, &lambda, &k, &h, &ix, &iy);

    struct reb_particle np = {0.};
    double p=0.,q=0.;
    reb_tools_solve_kepler_pal(h, k, lambda, &p, &q);

    double slp = sin(lambda+p);
    double clp = cos(lambda+p);
    
    double l = 1.-sqrt(1.-h*h-k*k);
    double iz = sqrt(fabs(4.-ix*ix-iy*iy));
    np.x = 0.0;
    np.y = 0.0;
    np.z = 0.0;

    double dan_dmm = -0.25*sqrt(G/(a*(po.m+primary.m)*(po.m+primary.m)*(po.m+primary.m)));
    double ddxi_dmm  = dan_dmm/(1.-q)*(-slp+q/(2.-l)*h);
    double ddeta_dmm = dan_dmm/(1.-q)*(+clp-q/(2.-l)*k);

    double ddW_dmm = ddeta_dmm*ix-ddxi_dmm*iy;
    np.vx = ddxi_dmm+0.5*iy*ddW_dmm;
    np.vy = ddeta_dmm-0.5*ix*ddW_dmm;
    np.vz = 0.5*iz*ddW_dmm;

    return np;
}

struct reb_particle reb_derivatives_e(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double cosf = cos(o.f);
    double v0 = sqrt(G*(po.m+primary.m)/o.a/(1.-o.e*o.e)); 
    double dr = -o.a*(cosf*o.e*o.e+cosf+2.*o.e)/((cosf*o.e+1.)*(cosf*o.e+1.));
    double dv0 = sqrt(G*(po.m+primary.m)/o.a)*o.e/((1.-o.e*o.e)*sqrt(1.-o.e*o.e)); 

    double cO = cos(o.Omega);
    double sO = sin(o.Omega);
    double co = cos(o.omega);
    double so = sin(o.omega);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double ci = cos(o.inc);
    double si = sin(o.inc);
    
    p.x = dr*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci);
    p.y = dr*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci);
    p.z = dr*(so*cf+co*sf)*si;

    p.vx = dv0*((o.e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO));
    p.vy = dv0*((o.e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO));
    p.vz = dv0*((o.e+cf)*co*si - sf*si*so);
    
    p.vx += v0*(-ci*co*sO - cO*so);
    p.vy += v0*(ci*co*cO - sO*so);
    p.vz += v0*(co*si);

    return p;
}



struct reb_particle reb_derivatives_e_e(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double cosf = cos(o.f);
    double ddr = o.a*2.*(cosf*cosf-1.)/((cosf*o.e+1.)*(cosf*o.e+1.)*(cosf*o.e+1.));
    double v0 = sqrt(G*(po.m+primary.m)/o.a/(1.-o.e*o.e)); 
    double dv0 = o.e*v0/(1.-o.e*o.e); 
    double ddv0 = v0/((o.e*o.e-1.)*(o.e*o.e-1.)) * (2.*o.e*o.e+1.);

    double cO = cos(o.Omega);
    double sO = sin(o.Omega);
    double co = cos(o.omega);
    double so = sin(o.omega);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double ci = cos(o.inc);
    double si = sin(o.inc);
    
    p.x = ddr*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci);
    p.y = ddr*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci);
    p.z = ddr*(so*cf+co*sf)*si;

    p.vx = ddv0*((o.e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO));
    p.vy = ddv0*((o.e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO));
    p.vz = ddv0*((o.e+cf)*co*si - sf*si*so);
    
    p.vx += 2.*dv0*(-ci*co*sO - cO*so);
    p.vy += 2.*dv0*(ci*co*cO - sO*so);
    p.vz += 2.*dv0*(co*si);

    return p;
}

struct reb_particle reb_derivatives_inc(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double r = o.a*(1.-o.e*o.e)/(1. + o.e*cos(o.f));
    double v0 = sqrt(G*(po.m+primary.m)/o.a/(1.-o.e*o.e)); 

    double cO = cos(o.Omega);
    double sO = sin(o.Omega);
    double co = cos(o.omega);
    double so = sin(o.omega);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double dci = -sin(o.inc);
    double dsi = cos(o.inc);
    
    p.x = r*(- sO*(so*cf+co*sf)*dci);
    p.y = r*(+ cO*(so*cf+co*sf)*dci);
    p.z = r*(so*cf+co*sf)*dsi;

    p.vx = v0*((o.e+cf)*(-dci*co*sO) - sf*(- dci*so*sO));
    p.vy = v0*((o.e+cf)*(dci*co*cO)  - sf*(dci*so*cO));
    p.vz = v0*((o.e+cf)*co*dsi - sf*dsi*so);
    

    return p;
}

struct reb_particle reb_derivatives_inc_inc(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double r = o.a*(1.-o.e*o.e)/(1. + o.e*cos(o.f));
    double v0 = sqrt(G*(po.m+primary.m)/o.a/(1.-o.e*o.e)); 

    double cO = cos(o.Omega);
    double sO = sin(o.Omega);
    double co = cos(o.omega);
    double so = sin(o.omega);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double ddci = -cos(o.inc);
    double ddsi = -sin(o.inc);
    
    p.x = r*(- sO*(so*cf+co*sf)*ddci);
    p.y = r*(+ cO*(so*cf+co*sf)*ddci);
    p.z = r*(so*cf+co*sf)*ddsi;

    p.vx = v0*((o.e+cf)*(-ddci*co*sO) - sf*(- ddci*so*sO));
    p.vy = v0*((o.e+cf)*(ddci*co*cO)  - sf*(ddci*so*cO));
    p.vz = v0*((o.e+cf)*co*ddsi - sf*ddsi*so);

    return p;
}

struct reb_particle reb_derivatives_Omega(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double r = o.a*(1.-o.e*o.e)/(1. + o.e*cos(o.f));
    double v0 = sqrt(G*(po.m+primary.m)/o.a/(1.-o.e*o.e)); 

    double dcO = -sin(o.Omega);
    double dsO = cos(o.Omega);
    double co = cos(o.omega);
    double so = sin(o.omega);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double ci = cos(o.inc);
    
    p.x = r*(dcO*(co*cf-so*sf) - dsO*(so*cf+co*sf)*ci);
    p.y = r*(dsO*(co*cf-so*sf) + dcO*(so*cf+co*sf)*ci);
    p.z = 0.;

    p.vx = v0*((o.e+cf)*(-ci*co*dsO - dcO*so) - sf*(co*dcO - ci*so*dsO));
    p.vy = v0*((o.e+cf)*(ci*co*dcO - dsO*so)  - sf*(co*dsO + ci*so*dcO));
    p.vz = 0.;

    return p;
}

struct reb_particle reb_derivatives_Omega_Omega(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double r = o.a*(1.-o.e*o.e)/(1. + o.e*cos(o.f));
    double v0 = sqrt(G*(po.m+primary.m)/o.a/(1.-o.e*o.e)); 

    double ddcO = -cos(o.Omega);
    double ddsO = -sin(o.Omega);
    double co = cos(o.omega);
    double so = sin(o.omega);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double ci = cos(o.inc);
    
    p.x = r*(ddcO*(co*cf-so*sf) - ddsO*(so*cf+co*sf)*ci);
    p.y = r*(ddsO*(co*cf-so*sf) + ddcO*(so*cf+co*sf)*ci);
    p.z = 0.;

    p.vx = v0*((o.e+cf)*(-ci*co*ddsO - ddcO*so) - sf*(co*ddcO - ci*so*ddsO));
    p.vy = v0*((o.e+cf)*(ci*co*ddcO - ddsO*so)  - sf*(co*ddsO + ci*so*ddcO));
    p.vz = 0.;

    return p;
}

struct reb_particle reb_derivatives_omega(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double r = o.a*(1.-o.e*o.e)/(1. + o.e*cos(o.f));
    double v0 = sqrt(G*(po.m+primary.m)/o.a/(1.-o.e*o.e)); 

    double cO = cos(o.Omega);
    double sO = sin(o.Omega);
    double dco = -sin(o.omega);
    double dso = cos(o.omega);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double ci = cos(o.inc);
    double si = sin(o.inc);
    
    p.x = r*(cO*(dco*cf-dso*sf) - sO*(dso*cf+dco*sf)*ci);
    p.y = r*(sO*(dco*cf-dso*sf) + cO*(dso*cf+dco*sf)*ci);
    p.z = r*(dso*cf+dco*sf)*si;

    p.vx = v0*((o.e+cf)*(-ci*dco*sO - cO*dso) - sf*(dco*cO - ci*dso*sO));
    p.vy = v0*((o.e+cf)*(ci*dco*cO - sO*dso)  - sf*(dco*sO + ci*dso*cO));
    p.vz = v0*((o.e+cf)*dco*si - sf*si*dso);
    
    return p;
}

struct reb_particle reb_derivatives_omega_omega(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double r = o.a*(1.-o.e*o.e)/(1. + o.e*cos(o.f));
    double v0 = sqrt(G*(po.m+primary.m)/o.a/(1.-o.e*o.e)); 

    double cO = cos(o.Omega);
    double sO = sin(o.Omega);
    double ddco = -cos(o.omega);
    double ddso = -sin(o.omega);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double ci = cos(o.inc);
    double si = sin(o.inc);
    
    p.x = r*(cO*(ddco*cf-ddso*sf) - sO*(ddso*cf+ddco*sf)*ci);
    p.y = r*(sO*(ddco*cf-ddso*sf) + cO*(ddso*cf+ddco*sf)*ci);
    p.z = r*(ddso*cf+ddco*sf)*si;

    p.vx = v0*((o.e+cf)*(-ci*ddco*sO - cO*ddso) - sf*(ddco*cO - ci*ddso*sO));
    p.vy = v0*((o.e+cf)*(ci*ddco*cO - sO*ddso)  - sf*(ddco*sO + ci*ddso*cO));
    p.vz = v0*((o.e+cf)*ddco*si - sf*si*ddso);
    
    return p;
}

struct reb_particle reb_derivatives_f(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double r = o.a*(1.-o.e*o.e)/(1. + o.e*cos(o.f));
    double dr = o.a*(1.-o.e*o.e)/((1. + o.e*cos(o.f))*(1. + o.e*cos(o.f)))*o.e*sin(o.f);
    double v0 = sqrt(G*(po.m+primary.m)/o.a/(1.-o.e*o.e)); 

    double cO = cos(o.Omega);
    double sO = sin(o.Omega);
    double co = cos(o.omega);
    double so = sin(o.omega);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double dcf = -sin(o.f);
    double dsf = cos(o.f);
    double ci = cos(o.inc);
    double si = sin(o.inc);
    
    p.x = dr*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci);
    p.y = dr*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci);
    p.z = dr*(so*cf+co*sf)*si;
    
    p.x += r*(cO*(co*dcf-so*dsf) - sO*(so*dcf+co*dsf)*ci);
    p.y += r*(sO*(co*dcf-so*dsf) + cO*(so*dcf+co*dsf)*ci);
    p.z += r*(so*dcf+co*dsf)*si;

    p.vx = v0*(dcf*(-ci*co*sO - cO*so) - dsf*(co*cO - ci*so*sO));
    p.vy = v0*(dcf*(ci*co*cO - sO*so)  - dsf*(co*sO + ci*so*cO));
    p.vz = v0*(dcf*co*si - dsf*si*so);
    
    return p;
}

struct reb_particle reb_derivatives_f_f(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double r = o.a*(1.-o.e*o.e)/(1. + o.e*cos(o.f));
    double dr = o.a*(1.-o.e*o.e)/((1. + o.e*cos(o.f))*(1. + o.e*cos(o.f)))*o.e*sin(o.f);
    double ddr = 2.*o.a*(1.-o.e*o.e)/((1. + o.e*cos(o.f))*(1. + o.e*cos(o.f))*(1. + o.e*cos(o.f)))*o.e*o.e*sin(o.f)*sin(o.f) + o.a*(1.-o.e*o.e)*o.e*cos(o.f)/((1. + o.e*cos(o.f))*(1. + o.e*cos(o.f)));
    double v0 = sqrt(G*(po.m+primary.m)/o.a/(1.-o.e*o.e)); 

    double cO = cos(o.Omega);
    double sO = sin(o.Omega);
    double co = cos(o.omega);
    double so = sin(o.omega);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double dcf = -sin(o.f);
    double dsf = cos(o.f);
    double ddcf = -cos(o.f);
    double ddsf = -sin(o.f);
    double ci = cos(o.inc);
    double si = sin(o.inc);
    
    p.x = ddr*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci);
    p.y = ddr*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci);
    p.z = ddr*(so*cf+co*sf)*si;
    
    p.x += 2.*dr*(cO*(co*dcf-so*dsf) - sO*(so*dcf+co*dsf)*ci);
    p.y += 2.*dr*(sO*(co*dcf-so*dsf) + cO*(so*dcf+co*dsf)*ci);
    p.z += 2.*dr*(so*dcf+co*dsf)*si;
    
    p.x += r*(cO*(co*ddcf-so*ddsf) - sO*(so*ddcf+co*ddsf)*ci);
    p.y += r*(sO*(co*ddcf-so*ddsf) + cO*(so*ddcf+co*ddsf)*ci);
    p.z += r*(so*ddcf+co*ddsf)*si;

    p.vx = v0*(ddcf*(-ci*co*sO - cO*so) - ddsf*(co*cO - ci*so*sO));
    p.vy = v0*(ddcf*(ci*co*cO - sO*so)  - ddsf*(co*sO + ci*so*cO));
    p.vz = v0*(ddcf*co*si - ddsf*si*so);
    
    return p;
}


struct reb_particle reb_derivatives_a_e(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double cosf = cos(o.f);
    double ddr = -(cosf*o.e*o.e+cosf+2.*o.e)/((cosf*o.e+1.)*(cosf*o.e+1.));
    double dv0_da = -0.5/sqrt(G*(po.m+primary.m)/o.a/(1.-o.e*o.e))*G*(po.m+primary.m)/(o.a*o.a)/(1.-o.e*o.e); 
    
    double dv0_da_de = o.e*dv0_da/(1.-o.e*o.e); 

    double cO = cos(o.Omega);
    double sO = sin(o.Omega);
    double co = cos(o.omega);
    double so = sin(o.omega);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double ci = cos(o.inc);
    double si = sin(o.inc);
    
    p.x = ddr*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci);
    p.y = ddr*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci);
    p.z = ddr*(so*cf+co*sf)*si;

    p.vx = dv0_da_de*((o.e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO));
    p.vy = dv0_da_de*((o.e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO));
    p.vz = dv0_da_de*((o.e+cf)*co*si - sf*si*so);
    
    p.vx += dv0_da*(-ci*co*sO - cO*so);
    p.vy += dv0_da*(ci*co*cO - sO*so);
    p.vz += dv0_da*(co*si);

    return p;
}

struct reb_particle reb_derivatives_a_inc(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double dr = (1.-o.e*o.e)/(1. + o.e*cos(o.f));
    double dv0 = -0.5/sqrt(o.a*o.a*o.a)*sqrt(G*(po.m+primary.m)/(1.-o.e*o.e)); 

    double cO = cos(o.Omega);
    double sO = sin(o.Omega);
    double co = cos(o.omega);
    double so = sin(o.omega);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double dci = -sin(o.inc);
    double dsi = cos(o.inc);
    
    p.x = dr*(- sO*(so*cf+co*sf)*dci);
    p.y = dr*(+ cO*(so*cf+co*sf)*dci);
    p.z = dr*(so*cf+co*sf)*dsi;

    p.vx = dv0*((o.e+cf)*(-dci*co*sO) - sf*(- dci*so*sO));
    p.vy = dv0*((o.e+cf)*(dci*co*cO)  - sf*(dci*so*cO));
    p.vz = dv0*((o.e+cf)*co*dsi - sf*dsi*so);
    
    return p;
}

struct reb_particle reb_derivatives_a_Omega(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double dr = (1.-o.e*o.e)/(1. + o.e*cos(o.f));
    double dv0 = -0.5/sqrt(o.a*o.a*o.a)*sqrt(G*(po.m+primary.m)/(1.-o.e*o.e)); 

    double dcO = -sin(o.Omega);
    double dsO = cos(o.Omega);
    double co = cos(o.omega);
    double so = sin(o.omega);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double ci = cos(o.inc);
    
    p.x = dr*(dcO*(co*cf-so*sf) - dsO*(so*cf+co*sf)*ci);
    p.y = dr*(dsO*(co*cf-so*sf) + dcO*(so*cf+co*sf)*ci);

    p.vx = dv0*((o.e+cf)*(-ci*co*dsO - dcO*so) - sf*(co*dcO - ci*so*dsO));
    p.vy = dv0*((o.e+cf)*(ci*co*dcO - dsO*so)  - sf*(co*dsO + ci*so*dcO));
    
    return p;
}

struct reb_particle reb_derivatives_a_omega(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double dr = (1.-o.e*o.e)/(1. + o.e*cos(o.f));
    double dv0 = -0.5/sqrt(o.a*o.a*o.a)*sqrt(G*(po.m+primary.m)/(1.-o.e*o.e)); 

    double cO = cos(o.Omega);
    double sO = sin(o.Omega);
    double dco = -sin(o.omega);
    double dso = cos(o.omega);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double ci = cos(o.inc);
    double si = sin(o.inc);
    
    p.x = dr*(cO*(dco*cf-dso*sf) - sO*(dso*cf+dco*sf)*ci);
    p.y = dr*(sO*(dco*cf-dso*sf) + cO*(dso*cf+dco*sf)*ci);
    p.z = dr*(dso*cf+dco*sf)*si;

    p.vx = dv0*((o.e+cf)*(-ci*dco*sO - cO*dso) - sf*(dco*cO - ci*dso*sO));
    p.vy = dv0*((o.e+cf)*(ci*dco*cO - sO*dso)  - sf*(dco*sO + ci*dso*cO));
    p.vz = dv0*((o.e+cf)*dco*si - sf*si*dso);
    
    return p;
}

struct reb_particle reb_derivatives_a_f(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double dr = (1.-o.e*o.e)/(1. + o.e*cos(o.f));
    double ddr = o.e*sin(o.f)*(1.-o.e*o.e)/(1. + o.e*cos(o.f))/(1. + o.e*cos(o.f));
    double dv0 = -0.5/sqrt(o.a*o.a*o.a)*sqrt(G*(po.m+primary.m)/(1.-o.e*o.e));

    double cO = cos(o.Omega);
    double sO = sin(o.Omega);
    double co = cos(o.omega);
    double so = sin(o.omega);
    double dcf = -sin(o.f);
    double dsf = cos(o.f);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double ci = cos(o.inc);
    double si = sin(o.inc);
    
    p.x = dr*(cO*(co*dcf-so*dsf) - sO*(so*dcf+co*dsf)*ci);
    p.y = dr*(sO*(co*dcf-so*dsf) + cO*(so*dcf+co*dsf)*ci);
    p.z = dr*(so*dcf+co*dsf)*si;
    
    p.x += ddr*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci);
    p.y += ddr*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci);
    p.z += ddr*(so*cf+co*sf)*si;

    p.vx = dv0*(dcf*(-ci*co*sO - cO*so) - dsf*(co*cO - ci*so*sO));
    p.vy = dv0*(dcf*(ci*co*cO - sO*so)  - dsf*(co*sO + ci*so*cO));
    p.vz = dv0*(dcf*co*si - dsf*si*so);
    
    return p;
}

struct reb_particle reb_derivatives_e_inc(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double cosf = cos(o.f);
    double v0 = sqrt(G*(po.m+primary.m)/o.a/(1.-o.e*o.e)); 
    double dr = -o.a*(cosf*o.e*o.e+cosf+2.*o.e)/((cosf*o.e+1.)*(cosf*o.e+1.));
    double dv0 = sqrt(G*(po.m+primary.m)/o.a)*o.e/((1.-o.e*o.e)*sqrt(1.-o.e*o.e)); 

    double cO = cos(o.Omega);
    double sO = sin(o.Omega);
    double co = cos(o.omega);
    double so = sin(o.omega);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double dci = -sin(o.inc);
    double dsi = cos(o.inc);
    
    p.x = dr*(- sO*(so*cf+co*sf)*dci);
    p.y = dr*(+ cO*(so*cf+co*sf)*dci);
    p.z = dr*(so*cf+co*sf)*dsi;

    p.vx = dv0*((o.e+cf)*(-dci*co*sO) - sf*(- dci*so*sO));
    p.vy = dv0*((o.e+cf)*(dci*co*cO)  - sf*(+ dci*so*cO));
    p.vz = dv0*((o.e+cf)*co*dsi - sf*dsi*so);
    
    p.vx += v0*(-dci*co*sO);
    p.vy += v0*(dci*co*cO);
    p.vz += v0*(co*dsi);

    return p;
}


struct reb_particle reb_derivatives_e_Omega(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double cosf = cos(o.f);
    double v0 = sqrt(G*(po.m+primary.m)/o.a/(1.-o.e*o.e)); 
    double dr = -o.a*(cosf*o.e*o.e+cosf+2.*o.e)/((cosf*o.e+1.)*(cosf*o.e+1.));
    double dv0 = sqrt(G*(po.m+primary.m)/o.a)*o.e/((1.-o.e*o.e)*sqrt(1.-o.e*o.e)); 

    double dcO = -sin(o.Omega);
    double dsO = cos(o.Omega);
    double co = cos(o.omega);
    double so = sin(o.omega);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double ci = cos(o.inc);
    
    p.x = dr*(dcO*(co*cf-so*sf) - dsO*(so*cf+co*sf)*ci);
    p.y = dr*(dsO*(co*cf-so*sf) + dcO*(so*cf+co*sf)*ci);

    p.vx = dv0*((o.e+cf)*(-ci*co*dsO - dcO*so) - sf*(co*dcO - ci*so*dsO));
    p.vy = dv0*((o.e+cf)*(ci*co*dcO - dsO*so)  - sf*(co*dsO + ci*so*dcO));
    
    p.vx += v0*(-ci*co*dsO - dcO*so);
    p.vy += v0*(ci*co*dcO - dsO*so);

    return p;
}

struct reb_particle reb_derivatives_e_omega(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double cosf = cos(o.f);
    double v0 = sqrt(G*(po.m+primary.m)/o.a/(1.-o.e*o.e)); 
    double dr = -o.a*(cosf*o.e*o.e+cosf+2.*o.e)/((cosf*o.e+1.)*(cosf*o.e+1.));
    double dv0 = sqrt(G*(po.m+primary.m)/o.a)*o.e/((1.-o.e*o.e)*sqrt(1.-o.e*o.e)); 

    double cO = cos(o.Omega);
    double sO = sin(o.Omega);
    double dco = -sin(o.omega);
    double dso = cos(o.omega);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double ci = cos(o.inc);
    double si = sin(o.inc);
    
    p.x = dr*(cO*(dco*cf-dso*sf) - sO*(dso*cf+dco*sf)*ci);
    p.y = dr*(sO*(dco*cf-dso*sf) + cO*(dso*cf+dco*sf)*ci);
    p.z = dr*(dso*cf+dco*sf)*si;

    p.vx = dv0*((o.e+cf)*(-ci*dco*sO - cO*dso) - sf*(dco*cO - ci*dso*sO));
    p.vy = dv0*((o.e+cf)*(ci*dco*cO - sO*dso)  - sf*(dco*sO + ci*dso*cO));
    p.vz = dv0*((o.e+cf)*dco*si - sf*si*dso);
    
    p.vx += v0*(-ci*dco*sO - cO*dso);
    p.vy += v0*(ci*dco*cO - sO*dso);
    p.vz += v0*(dco*si);

    return p;
}
struct reb_particle reb_derivatives_e_f(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double cosf = cos(o.f);
    double dr = -o.a*(cosf*o.e*o.e+cosf+2.*o.e)/((cosf*o.e+1.)*(cosf*o.e+1.));
    double ddr = -o.a*(-sin(o.f)*o.e*o.e-sin(o.f))/((cosf*o.e+1.)*(cosf*o.e+1.))
                -2.*o.e*sin(o.f) * o.a*(cosf*o.e*o.e+cosf+2.*o.e)/((cosf*o.e+1.)*(cosf*o.e+1.)*(cosf*o.e+1.));
    double dv0 = sqrt(G*(po.m+primary.m)/o.a)*o.e/((1.-o.e*o.e)*sqrt(1.-o.e*o.e)); 

    double cO = cos(o.Omega);
    double sO = sin(o.Omega);
    double co = cos(o.omega);
    double so = sin(o.omega);
    double dcf = -sin(o.f);
    double dsf = cos(o.f);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double ci = cos(o.inc);
    double si = sin(o.inc);
    
    p.x = dr*(cO*(co*dcf-so*dsf) - sO*(so*dcf+co*dsf)*ci);
    p.y = dr*(sO*(co*dcf-so*dsf) + cO*(so*dcf+co*dsf)*ci);
    p.z = dr*(so*dcf+co*dsf)*si;
    
    p.x += ddr*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci);
    p.y += ddr*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci);
    p.z += ddr*(so*cf+co*sf)*si;
    
    p.vx = dv0*(dcf*(-ci*co*sO - cO*so) - dsf*(co*cO - ci*so*sO));
    p.vy = dv0*(dcf*(ci*co*cO - sO*so)  - dsf*(co*sO + ci*so*cO));
    p.vz = dv0*(dcf*co*si - dsf*si*so);
    
    return p;
}
struct reb_particle reb_derivatives_m_e(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double dv0m = 0.5*G/o.a/(1.-o.e*o.e)/sqrt(G*(po.m+primary.m)/o.a/(1.-o.e*o.e)); 
    double dv0ea = 0.5*G/o.a/sqrt(G*(po.m+primary.m)/o.a)*o.e/((1.-o.e*o.e)*sqrt(1.-o.e*o.e)); 

    double cO = cos(o.Omega);
    double sO = sin(o.Omega);
    double co = cos(o.omega);
    double so = sin(o.omega);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double ci = cos(o.inc);
    double si = sin(o.inc);
    
    p.vx = dv0ea*((o.e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO));
    p.vy = dv0ea*((o.e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO));
    p.vz = dv0ea*((o.e+cf)*co*si - sf*si*so);
    
    p.vx += dv0m*(-ci*co*sO - cO*so);
    p.vy += dv0m*(ci*co*cO - sO*so);
    p.vz += dv0m*(co*si);

    return p;
}

struct reb_particle reb_derivatives_inc_Omega(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double r = o.a*(1.-o.e*o.e)/(1. + o.e*cos(o.f));
    double v0 = sqrt(G*(po.m+primary.m)/o.a/(1.-o.e*o.e)); 

    double dcO = -sin(o.Omega);
    double dsO = cos(o.Omega);
    double co = cos(o.omega);
    double so = sin(o.omega);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double dci = -sin(o.inc);
    
    p.x = r*(- dsO*(so*cf+co*sf)*dci);
    p.y = r*(+ dcO*(so*cf+co*sf)*dci);

    p.vx = v0*((o.e+cf)*(-dci*co*dsO) - sf*(- dci*so*dsO));
    p.vy = v0*((o.e+cf)*(dci*co*dcO)  - sf*(dci*so*dcO));

    return p;
}

struct reb_particle reb_derivatives_inc_omega(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double r = o.a*(1.-o.e*o.e)/(1. + o.e*cos(o.f));
    double v0 = sqrt(G*(po.m+primary.m)/o.a/(1.-o.e*o.e)); 

    double cO = cos(o.Omega);
    double sO = sin(o.Omega);
    double dco = -sin(o.omega);
    double dso = cos(o.omega);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double dci = -sin(o.inc);
    double dsi = cos(o.inc);
    
    p.x = r*(- sO*(dso*cf+dco*sf)*dci);
    p.y = r*(+ cO*(dso*cf+dco*sf)*dci);
    p.z = r*(dso*cf+dco*sf)*dsi;

    p.vx = v0*((o.e+cf)*(-dci*dco*sO) - sf*(- dci*dso*sO));
    p.vy = v0*((o.e+cf)*(dci*dco*cO)  - sf*(dci*dso*cO));
    p.vz = v0*((o.e+cf)*dco*dsi - sf*dsi*dso);

    return p;
}

struct reb_particle reb_derivatives_inc_f(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double r = o.a*(1.-o.e*o.e)/(1. + o.e*cos(o.f));
    double dr = o.e*sin(o.f)*o.a*(1.-o.e*o.e)/(1. + o.e*cos(o.f))/(1. + o.e*cos(o.f));
    double v0 = sqrt(G*(po.m+primary.m)/o.a/(1.-o.e*o.e)); 

    double cO = cos(o.Omega);
    double sO = sin(o.Omega);
    double co = cos(o.omega);
    double so = sin(o.omega);
    double dcf = -sin(o.f);
    double dsf = cos(o.f);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double dci = -sin(o.inc);
    double dsi = cos(o.inc);
    
    p.x = r*(- sO*(so*dcf+co*dsf)*dci);
    p.y = r*(+ cO*(so*dcf+co*dsf)*dci);
    p.z = r*(so*dcf+co*dsf)*dsi;
    
    p.x += dr*(- sO*(so*cf+co*sf)*dci);
    p.y += dr*(+ cO*(so*cf+co*sf)*dci);
    p.z += dr*(so*cf+co*sf)*dsi;

    p.vx = v0*(dcf*(-dci*co*sO) - dsf*(- dci*so*sO));
    p.vy = v0*(dcf*(dci*co*cO)  - dsf*(dci*so*cO));
    p.vz = v0*(dcf*co*dsi - dsf*dsi*so);

    return p;
}

struct reb_particle reb_derivatives_m_inc(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double dv0 = 0.5/sqrt(po.m+primary.m)*sqrt(G/o.a/(1.-o.e*o.e)); 

    double cO = cos(o.Omega);
    double sO = sin(o.Omega);
    double co = cos(o.omega);
    double so = sin(o.omega);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double dci = -sin(o.inc);
    double dsi = cos(o.inc);
    
    p.vx = dv0*((o.e+cf)*(-dci*co*sO) - sf*(- dci*so*sO));
    p.vy = dv0*((o.e+cf)*(dci*co*cO)  - sf*(dci*so*cO));
    p.vz = dv0*((o.e+cf)*co*dsi - sf*dsi*so);
    
    return p;
}

struct reb_particle reb_derivatives_omega_Omega(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double r = o.a*(1.-o.e*o.e)/(1. + o.e*cos(o.f));
    double v0 = sqrt(G*(po.m+primary.m)/o.a/(1.-o.e*o.e)); 

    double dcO = -sin(o.Omega);
    double dsO = cos(o.Omega);
    double dco = -sin(o.omega);
    double dso = cos(o.omega);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double ci = cos(o.inc);
    
    p.x = r*(dcO*(dco*cf-dso*sf) - dsO*(dso*cf+dco*sf)*ci);
    p.y = r*(dsO*(dco*cf-dso*sf) + dcO*(dso*cf+dco*sf)*ci);

    p.vx = v0*((o.e+cf)*(-ci*dco*dsO - dcO*dso) - sf*(dco*dcO - ci*dso*dsO));
    p.vy = v0*((o.e+cf)*(ci*dco*dcO - dsO*dso)  - sf*(dco*dsO + ci*dso*dcO));

    return p;
}

struct reb_particle reb_derivatives_Omega_f(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double r = o.a*(1.-o.e*o.e)/(1. + o.e*cos(o.f));
    double dr = o.e*sin(o.f)*o.a*(1.-o.e*o.e)/(1. + o.e*cos(o.f))/(1. + o.e*cos(o.f));
    double v0 = sqrt(G*(po.m+primary.m)/o.a/(1.-o.e*o.e)); 

    double dcO = -sin(o.Omega);
    double dsO = cos(o.Omega);
    double co = cos(o.omega);
    double so = sin(o.omega);
    double dcf = -sin(o.f);
    double dsf = cos(o.f);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double ci = cos(o.inc);
    
    p.x = r*(dcO*(co*dcf-so*dsf) - dsO*(so*dcf+co*dsf)*ci);
    p.y = r*(dsO*(co*dcf-so*dsf) + dcO*(so*dcf+co*dsf)*ci);
    
    p.x += dr*(dcO*(co*cf-so*sf) - dsO*(so*cf+co*sf)*ci);
    p.y += dr*(dsO*(co*cf-so*sf) + dcO*(so*cf+co*sf)*ci);

    p.vx = v0*((dcf)*(-ci*co*dsO - dcO*so) - dsf*(co*dcO - ci*so*dsO));
    p.vy = v0*((dcf)*(ci*co*dcO - dsO*so)  - dsf*(co*dsO + ci*so*dcO));

    return p;
}

struct reb_particle reb_derivatives_m_Omega(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double dv0 = 0.5/sqrt(po.m+primary.m)*sqrt(G/o.a/(1.-o.e*o.e)); 

    double dcO = -sin(o.Omega);
    double dsO = cos(o.Omega);
    double co = cos(o.omega);
    double so = sin(o.omega);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double ci = cos(o.inc);
    
    p.vx = dv0*((o.e+cf)*(-ci*co*dsO - dcO*so) - sf*(co*dcO - ci*so*dsO));
    p.vy = dv0*((o.e+cf)*(ci*co*dcO - dsO*so)  - sf*(co*dsO + ci*so*dcO));

    return p;
}

struct reb_particle reb_derivatives_omega_f(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double r = o.a*(1.-o.e*o.e)/(1. + o.e*cos(o.f));
    double dr = o.e*sin(o.f)*o.a*(1.-o.e*o.e)/(1. + o.e*cos(o.f))/(1. + o.e*cos(o.f));
    double v0 = sqrt(G*(po.m+primary.m)/o.a/(1.-o.e*o.e)); 

    double cO = cos(o.Omega);
    double sO = sin(o.Omega);
    double dco = -sin(o.omega);
    double dso = cos(o.omega);
    double dcf = -sin(o.f);
    double dsf = cos(o.f);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double ci = cos(o.inc);
    double si = sin(o.inc);
    
    p.x = r*(cO*(dco*dcf-dso*dsf) - sO*(dso*dcf+dco*dsf)*ci);
    p.y = r*(sO*(dco*dcf-dso*dsf) + cO*(dso*dcf+dco*dsf)*ci);
    p.z = r*(dso*dcf+dco*dsf)*si;
    
    p.x += dr*(cO*(dco*cf-dso*sf) - sO*(dso*cf+dco*sf)*ci);
    p.y += dr*(sO*(dco*cf-dso*sf) + cO*(dso*cf+dco*sf)*ci);
    p.z += dr*(dso*cf+dco*sf)*si;

    p.vx = v0*((dcf)*(-ci*dco*sO - cO*dso) - dsf*(dco*cO - ci*dso*sO));
    p.vy = v0*((dcf)*(ci*dco*cO - sO*dso)  - dsf*(dco*sO + ci*dso*cO));
    p.vz = v0*((dcf)*dco*si - dsf*si*dso);
    
    return p;
}

struct reb_particle reb_derivatives_m_omega(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double dv0 = 0.5*sqrt(G/o.a/(1.-o.e*o.e))/sqrt(po.m+primary.m); 

    double cO = cos(o.Omega);
    double sO = sin(o.Omega);
    double dco = -sin(o.omega);
    double dso = cos(o.omega);
    double cf = cos(o.f);
    double sf = sin(o.f);
    double ci = cos(o.inc);
    double si = sin(o.inc);
    
    p.vx = dv0*((o.e+cf)*(-ci*dco*sO - cO*dso) - sf*(dco*cO - ci*dso*sO));
    p.vy = dv0*((o.e+cf)*(ci*dco*cO - sO*dso)  - sf*(dco*sO + ci*dso*cO));
    p.vz = dv0*((o.e+cf)*dco*si - sf*si*dso);
    
    return p;
}

struct reb_particle reb_derivatives_m_f(double G, struct reb_particle primary, struct reb_particle po){
    struct reb_orbit o = reb_tools_particle_to_orbit(G, po, primary);
    struct reb_particle p = {0};
    double dv0 = 0.5*sqrt(G/o.a/(1.-o.e*o.e))/sqrt(po.m+primary.m); 

    double cO = cos(o.Omega);
    double sO = sin(o.Omega);
    double co = cos(o.omega);
    double so = sin(o.omega);
    double dcf = -sin(o.f);
    double dsf = cos(o.f);
    double ci = cos(o.inc);
    double si = sin(o.inc);
    
    p.vx = dv0*(dcf*(-ci*co*sO - cO*so) - dsf*(co*cO - ci*so*sO));
    p.vy = dv0*(dcf*(ci*co*cO - sO*so)  - dsf*(co*sO + ci*so*cO));
    p.vz = dv0*(dcf*co*si - dsf*si*so);
    
    return p;
}
