/**
 * @file    integrator_eos.c
 * @brief   Embedded Operator Splitting (EOS) method
 * @author  Hanno Rein
 * 
 * @section LICENSE
 * Copyright (c) 2019 Hanno Rein
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
#include <string.h>
#include "rebound.h"
#include "integrator.h"
#include "gravity.h"
#include "integrator_eos.h"
#include "tools.h"
#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< Returns the minimum of a and b
#define MAX(a, b) ((a) > (b) ? (a) : (b))    ///< Returns the maximum of a and b

static const double lf4_a = 0.675603595979828817023843904485;

static const double lf6_a[5] = {0.1867, 0.5554970237124784, 0.1294669489134754, -0.843265623387734, 0.9432033015235604};

static const double lf8_a[9] = {0.128865979381443, 0.581514087105251, -0.410175371469850, 0.1851469357165877, -0.4095523434208514, 0.1444059410800120, 0.2783355003936797, 0.3149566839162949, -0.6269948254051343979}; 

static const double lf4_2_a = 0.211324865405187117745425609749;

static const double lf8_6_4_a[4] = {0.0711334264982231177779387300061549964174, 0.241153427956640098736487795326289649618, 0.521411761772814789212136078067994229991, -0.333698616227678005726562603400438876027};
static const double lf8_6_4_b[4] = {0.183083687472197221961703757166430291072, 0.310782859898574869507522291054262796375, -0.0265646185119588006972121379164987592663, 0.0653961422823734184559721793911134363710};

static const double pmlf6_a[2] = {-0.0682610383918630,0.568261038391863038121699}; 
static const double pmlf6_b[2] = {0.2621129352517028, 0.475774129496594366806050}; 
static const double pmlf6_c[2] = {0., 0.0164011128160783}; 
static const double pmlf6_z[6] = { 0.07943288242455420, 0.02974829169467665, -0.7057074964815896, 0.3190423451260838, -0.2869147334299646, 0.564398710666239478150885};
static const double pmlf6_y[6] = {1.3599424487455264, -0.6505973747535132, -0.033542814598338416, -0.040129915275115030, 0.044579729809902803, -0.680252073928462652752103};
static const double pmlf6_v[6] = {-0.034841228074994859, 0.031675672097525204, -0.005661054677711889, 0.004262222269023640, 0.005, -0.005};

static const double pmlf4_y[3] = {0.1859353996846055, 0.0731969797858114, -0.1576624269298081};
static const double pmlf4_z[3] = {0.8749306155955435, -0.237106680151022, -0.5363539829039128};

static const double plf7_6_4_a[2] = {0.5600879810924619,-0.060087981092461900000};
static const double plf7_6_4_b[2] = {1.5171479707207228, -2.0342959414414456000};
static const double plf7_6_4_z[6] = {-0.3346222298730800, 1.0975679907321640, -1.0380887460967830, 0.6234776317921379, -1.1027532063031910, -0.0141183222088869};
static const double plf7_6_4_y[6] = {-1.6218101180868010, 0.0061709468110142, 0.8348493592472594, -0.0511253369989315, 0.5633782670698199, -0.5};
                
static inline void reb_integrator_eos_interaction_shell0(struct reb_simulation* r, double y, double v){
    // Calculate gravity using standard gravity routine
    r->gravity_ignore_terms = 2;
    r->gravity = REB_GRAVITY_BASIC;
    reb_update_acceleration(r);
    if (v!=0.){
        reb_calculate_and_apply_jerk(r,v);
    }
    // Apply acceleration (jerk already applied)
    struct reb_particle* restrict const particles = r->particles;
    const int N = r->N;
    for (int i=0;i<N;i++){
        particles[i].vx += y*particles[i].ax;
        particles[i].vy += y*particles[i].ay;
        particles[i].vz += y*particles[i].az;
    }
}

static inline void reb_integrator_eos_interaction_shell1(struct reb_simulation* r, double y, double v){
    const int N = r->N;
	const int N_real   = N - r->N_var;
    const int N_active = r->N_active==-1?N_real:r->N_active;
    const int testparticle_type   = r->testparticle_type;
    struct reb_particle* restrict const particles = r->particles;

    const double G = r->G;
    
    if (v!=0.){ // is jerk even used?
        // Normal force calculation 
        particles[0].ax = 0;
        particles[0].ay = 0;
        particles[0].az = 0;
        // Interactions between central object and all other active particles
        for (int j=1; j<N_active; j++){
            const double dx = particles[0].x - particles[j].x;
            const double dy = particles[0].y - particles[j].y;
            const double dz = particles[0].z - particles[j].z;
            const double dr = sqrt(dx*dx + dy*dy + dz*dz);

            const double prefact = G/(dr*dr*dr);
            const double prefactj = -prefact*particles[j].m;
            particles[0].ax    += prefactj*dx;
            particles[0].ay    += prefactj*dy;
            particles[0].az    += prefactj*dz;
            const double prefacti = prefact*particles[0].m;
            particles[j].ax    = prefacti*dx;
            particles[j].ay    = prefacti*dy;
            particles[j].az    = prefacti*dz;
        }
        // Interactions between central object and all test particles
        for (int j=N_active; j<N_real; j++){
            const double dx = particles[0].x - particles[j].x;
            const double dy = particles[0].y - particles[j].y;
            const double dz = particles[0].z - particles[j].z;
            const double dr = sqrt(dx*dx + dy*dy + dz*dz);

            const double prefact = G/(dr*dr*dr);
            const double prefacti = prefact*particles[0].m;
            particles[j].ax    = prefacti*dx;
            particles[j].ay    = prefacti*dy;
            particles[j].az    = prefacti*dz;
            if (testparticle_type){
                const double prefactj = -prefact*particles[j].m;
                particles[0].ax    += prefactj*dx;
                particles[0].ay    += prefactj*dy;
                particles[0].az    += prefactj*dz;
            }
        }
        // Jerk calculation
        // Interactions between central object and all other active particles
        for (int i=1; i<N_active; i++){
            const double dx = particles[0].x - particles[i].x; 
            const double dy = particles[0].y - particles[i].y; 
            const double dz = particles[0].z - particles[i].z; 
            
            const double dax = particles[0].ax - particles[i].ax; 
            const double day = particles[0].ay - particles[i].ay; 
            const double daz = particles[0].az - particles[i].az; 

            const double dr = sqrt(dx*dx + dy*dy + dz*dz);
            const double alphasum = dax*dx+day*dy+daz*dz;
            const double prefact2 = 2.*v*G /(dr*dr*dr);
            const double prefact2i = prefact2*particles[i].m;
            const double prefact2j = prefact2*particles[0].m;
            const double prefact1 = alphasum*prefact2/dr *3./dr;
            const double prefact1i = prefact1*particles[i].m;
            const double prefact1j = prefact1*particles[0].m;
            particles[0].vx    += -dax*prefact2i + dx*prefact1i;
            particles[0].vy    += -day*prefact2i + dy*prefact1i;
            particles[0].vz    += -daz*prefact2i + dz*prefact1i;
            particles[i].vx    += y*particles[i].ax + dax*prefact2j - dx*prefact1j;
            particles[i].vy    += y*particles[i].ay + day*prefact2j - dy*prefact1j;
            particles[i].vz    += y*particles[i].az + daz*prefact2j - dz*prefact1j;
        }
        // Interactions between central object and all test particles
        for (int i=N_active; i<N_real; i++){
            const double dx = particles[0].x - particles[i].x; 
            const double dy = particles[0].y - particles[i].y; 
            const double dz = particles[0].z - particles[i].z; 
            
            const double dax = particles[0].ax - particles[i].ax; 
            const double day = particles[0].ay - particles[i].ay; 
            const double daz = particles[0].az - particles[i].az; 

            const double dr = sqrt(dx*dx + dy*dy + dz*dz);
            const double alphasum = dax*dx+day*dy+daz*dz;
            const double prefact2 = 2.*v*G /(dr*dr*dr);
            const double prefact2j = prefact2*particles[0].m;
            const double prefact1 = alphasum*prefact2/dr *3./dr;
            const double prefact1j = prefact1*particles[0].m;
            if (testparticle_type){
                const double prefact2i = prefact2*particles[i].m;
                const double prefact1i = prefact1*particles[i].m;
                particles[0].vx    += -dax*prefact2i + dx*prefact1i;
                particles[0].vy    += -day*prefact2i + dy*prefact1i;
                particles[0].vz    += -daz*prefact2i + dz*prefact1i;
            }
            particles[i].vx    += y*particles[i].ax + dax*prefact2j - dx*prefact1j;
            particles[i].vy    += y*particles[i].ay + day*prefact2j - dy*prefact1j;
            particles[i].vz    += y*particles[i].az + daz*prefact2j - dz*prefact1j;
        }
        particles[0].vx += y*particles[0].ax;
        particles[0].vy += y*particles[0].ay;
        particles[0].vz += y*particles[0].az;
    }else{
        // Normal force calculation 
        // Interactions between central object and all other active particles
        for (int j=1; j<N_active; j++){
            const double dx = particles[0].x - particles[j].x;
            const double dy = particles[0].y - particles[j].y;
            const double dz = particles[0].z - particles[j].z;
            const double dr = sqrt(dx*dx + dy*dy + dz*dz);

            const double prefact = y*G/(dr*dr*dr);
            const double prefactj = -prefact*particles[j].m;
            particles[0].vx    += prefactj*dx;
            particles[0].vy    += prefactj*dy;
            particles[0].vz    += prefactj*dz;
            const double prefacti = prefact*particles[0].m;
            particles[j].vx    += prefacti*dx;
            particles[j].vy    += prefacti*dy;
            particles[j].vz    += prefacti*dz;
        }
        // Interactions between central object and all test particles
        for (int j=N_active; j<N_real; j++){
            const double dx = particles[0].x - particles[j].x;
            const double dy = particles[0].y - particles[j].y;
            const double dz = particles[0].z - particles[j].z;
            const double dr = sqrt(dx*dx + dy*dy + dz*dz);

            const double prefact = y*G/(dr*dr*dr);
            const double prefacti = prefact*particles[0].m;
            particles[j].vx    += prefacti*dx;
            particles[j].vy    += prefacti*dy;
            particles[j].vz    += prefacti*dz;
            if (testparticle_type){
                const double prefactj = -prefact*particles[j].m;
                particles[0].vx    += prefactj*dx;
                particles[0].vy    += prefactj*dy;
                particles[0].vz    += prefactj*dz;
            }
        }
        for (int v=0;v<r->var_config_N;v++){
            struct reb_variational_configuration const vc = r->var_config[v];
            if (vc.order==1){
                //////////////////
                /// 1st order  ///
                //////////////////
                struct reb_particle* const particles_var1 = particles + vc.index;
                if (vc.testparticle<0){
                    for (int j=1; j<N_active; j++){
                        const double dx = particles[0].x - particles[j].x;
                        const double dy = particles[0].y - particles[j].y;
                        const double dz = particles[0].z - particles[j].z;
                        const double r2 = dx*dx + dy*dy + dz*dz;
                        const double _r  = sqrt(r2);
                        const double r3inv = 1./(r2*_r);
                        const double r5inv = 3.*r3inv/r2;
                        const double ddx = particles_var1[0].x - particles_var1[j].x;
                        const double ddy = particles_var1[0].y - particles_var1[j].y;
                        const double ddz = particles_var1[0].z - particles_var1[j].z;
                        const double Gmi = y*G * particles[0].m;
                        const double Gmj = y*G * particles[j].m;

                        // Variational equations
                        const double dxdx = dx*dx*r5inv - r3inv;
                        const double dydy = dy*dy*r5inv - r3inv;
                        const double dzdz = dz*dz*r5inv - r3inv;
                        const double dxdy = dx*dy*r5inv;
                        const double dxdz = dx*dz*r5inv;
                        const double dydz = dy*dz*r5inv;
                        const double dax =   ddx * dxdx + ddy * dxdy + ddz * dxdz;
                        const double day =   ddx * dxdy + ddy * dydy + ddz * dydz;
                        const double daz =   ddx * dxdz + ddy * dydz + ddz * dzdz;

                        // Variational mass contributions
                        const double dGmi = y*G*particles_var1[0].m;
                        const double dGmj = y*G*particles_var1[j].m;

                        particles_var1[0].vx += Gmj * dax - dGmj*r3inv*dx;
                        particles_var1[0].vy += Gmj * day - dGmj*r3inv*dy;
                        particles_var1[0].vz += Gmj * daz - dGmj*r3inv*dz;

                        particles_var1[j].vx -= Gmi * dax - dGmi*r3inv*dx;
                        particles_var1[j].vy -= Gmi * day - dGmi*r3inv*dy;
                        particles_var1[j].vz -= Gmi * daz - dGmi*r3inv*dz; 
                    }
                }else{ //testparticle
                    int i = vc.testparticle;
                    const double dx = particles[i].x - particles[0].x;
                    const double dy = particles[i].y - particles[0].y;
                    const double dz = particles[i].z - particles[0].z;
                    const double r2 = dx*dx + dy*dy + dz*dz;
                    const double _r  = sqrt(r2);
                    const double r3inv = 1./(r2*_r);
                    const double r5inv = 3.*r3inv/r2;
                    const double ddx = particles_var1[0].x;
                    const double ddy = particles_var1[0].y;
                    const double ddz = particles_var1[0].z;
                    const double Gmj = y*G*particles[0].m;

                    // Variational equations
                    const double dxdx = dx*dx*r5inv - r3inv;
                    const double dydy = dy*dy*r5inv - r3inv;
                    const double dzdz = dz*dz*r5inv - r3inv;
                    const double dxdy = dx*dy*r5inv;
                    const double dxdz = dx*dz*r5inv;
                    const double dydz = dy*dz*r5inv;
                    const double dax =   ddx * dxdx + ddy * dxdy + ddz * dxdz;
                    const double day =   ddx * dxdy + ddy * dydy + ddz * dydz;
                    const double daz =   ddx * dxdz + ddy * dydz + ddz * dzdz;

                    // No variational mass contributions for test particles!
                    particles_var1[0].vx += Gmj * dax;
                    particles_var1[0].vy += Gmj * day;
                    particles_var1[0].vz += Gmj * daz;
                }
            }
        }
    }

}
static inline void reb_integrator_eos_preprocessor(struct reb_simulation* const r, double dt, enum REB_EOS_TYPE type, void (*drift_step)(struct reb_simulation* const r, double a), void (*interaction_step)(struct reb_simulation* const r, double y, double v)){
    switch(type){
        case REB_EOS_PMLF6:
            for (int i=0;i<6;i++){
                drift_step(r, dt*pmlf6_z[i]);
                interaction_step(r, dt*pmlf6_y[i], dt*dt*dt*pmlf6_v[i]);
            }
            break;
        case REB_EOS_PMLF4:
            for (int i=0;i<3;i++){
                interaction_step(r, dt*pmlf4_y[i], 0.);
                drift_step(r, dt*pmlf4_z[i]);
            }
            break;
        case REB_EOS_PLF7_6_4:
            for (int i=0;i<6;i++){
                drift_step(r, dt*plf7_6_4_z[i]);
                interaction_step(r, dt*plf7_6_4_y[i], 0.);
            }
            break;
        default:
            break;
    }
}
static inline void reb_integrator_eos_postprocessor(struct reb_simulation* const r, double dt, enum REB_EOS_TYPE type, void (*drift_step)(struct reb_simulation* const r, double a), void (*interaction_step)(struct reb_simulation* const r, double y, double v)){
    switch(type){
        case REB_EOS_PMLF6:
            for (int i=5;i>=0;i--){
                interaction_step(r, -dt*pmlf6_y[i], -dt*dt*dt*pmlf6_v[i]); 
                drift_step(r, -dt*pmlf6_z[i]);
             }
            break;
        case REB_EOS_PMLF4:
            for (int i=2;i>=0;i--){
                drift_step(r, -dt*pmlf4_z[i]);
                interaction_step(r, -dt*pmlf4_y[i], 0.); 
             }
            break;
        case REB_EOS_PLF7_6_4:
            for (int i=5;i>=0;i--){
                interaction_step(r, -dt*plf7_6_4_y[i], 0.);
                drift_step(r, -dt*plf7_6_4_z[i]);
            }
            break;
        default:
            break;
    }
}
static void reb_integrator_eos_drift_shell1(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;
    unsigned int N = r->N;
    for (int i=0;i<N;i++){  
        particles[i].x += dt*particles[i].vx;
        particles[i].y += dt*particles[i].vy;
        particles[i].z += dt*particles[i].vz;
    } 
}

static void reb_integrator_eos_drift_shell0(struct reb_simulation* const r, double _dt){
    struct reb_simulation_integrator_eos* const reos = &(r->ri_eos);
    const int n = reos->n;
    const double dt = _dt/n;
    reb_integrator_eos_preprocessor(r, dt, reos->phi1, reb_integrator_eos_drift_shell1, reb_integrator_eos_interaction_shell1);
    switch(reos->phi1){
        case REB_EOS_LF:
            reb_integrator_eos_drift_shell1(r, dt*0.5); 
            for (int i=0;i<n;i++){
                reb_integrator_eos_interaction_shell1(r, dt, 0.);
                if (i<n-1){
                    reb_integrator_eos_drift_shell1(r, dt);
                }
            }
            reb_integrator_eos_drift_shell1(r, dt*0.5);
            break;
        case REB_EOS_LF4:
            reb_integrator_eos_drift_shell1(r, dt*lf4_a);
            for (int i=0;i<n;i++){
                reb_integrator_eos_interaction_shell1(r, dt*2.*lf4_a, 0.);
                reb_integrator_eos_drift_shell1(r, dt*(0.5-lf4_a));
                reb_integrator_eos_interaction_shell1(r, dt*(1.-4.*lf4_a), 0.);
                reb_integrator_eos_drift_shell1(r, dt*(0.5-lf4_a));
                reb_integrator_eos_interaction_shell1(r, dt*2.*lf4_a, 0.);
                if (i<n-1){
                    reb_integrator_eos_drift_shell1(r, dt*2.*lf4_a);
                }
            }
            reb_integrator_eos_drift_shell1(r, dt*lf4_a);
            break;
        case REB_EOS_LF6:
            reb_integrator_eos_drift_shell1(r, dt*lf6_a[0]*0.5);
            for (int i=0;i<n;i++){
                reb_integrator_eos_interaction_shell1(r, dt*lf6_a[0], 0.);
                reb_integrator_eos_drift_shell1(r, dt*(lf6_a[0]+lf6_a[1])*0.5);
                reb_integrator_eos_interaction_shell1(r, dt*lf6_a[1], 0.);
                reb_integrator_eos_drift_shell1(r, dt*(lf6_a[1]+lf6_a[2])*0.5);
                reb_integrator_eos_interaction_shell1(r, dt*lf6_a[2], 0.);
                reb_integrator_eos_drift_shell1(r, dt*(lf6_a[2]+lf6_a[3])*0.5);
                reb_integrator_eos_interaction_shell1(r, dt*lf6_a[3], 0.);
                reb_integrator_eos_drift_shell1(r, dt*(lf6_a[3]+lf6_a[4])*0.5);
                reb_integrator_eos_interaction_shell1(r, dt*lf6_a[4], 0.);
                reb_integrator_eos_drift_shell1(r, dt*(lf6_a[3]+lf6_a[4])*0.5);
                reb_integrator_eos_interaction_shell1(r, dt*lf6_a[3], 0.);
                reb_integrator_eos_drift_shell1(r, dt*(lf6_a[2]+lf6_a[3])*0.5);
                reb_integrator_eos_interaction_shell1(r, dt*lf6_a[2], 0.);
                reb_integrator_eos_drift_shell1(r, dt*(lf6_a[1]+lf6_a[2])*0.5);
                reb_integrator_eos_interaction_shell1(r, dt*lf6_a[1], 0.);
                reb_integrator_eos_drift_shell1(r, dt*(lf6_a[0]+lf6_a[1])*0.5);
                reb_integrator_eos_interaction_shell1(r, dt*lf6_a[0], 0.);
                if (i<n-1){
                    reb_integrator_eos_drift_shell1(r, dt*lf6_a[0]);
                }
            }
            reb_integrator_eos_drift_shell1(r, dt*lf6_a[0]*0.5);
            break; 
        case REB_EOS_LF8: 
            reb_integrator_eos_drift_shell1(r, dt*lf8_a[0]*0.5);
            for (int i=0;i<n;i++){
                reb_integrator_eos_interaction_shell1(r, dt*lf8_a[0], 0.);
                reb_integrator_eos_drift_shell1(r, dt*(lf8_a[0]+lf8_a[1])*0.5);
                reb_integrator_eos_interaction_shell1(r, dt*lf8_a[1], 0.);
                reb_integrator_eos_drift_shell1(r, dt*(lf8_a[1]+lf8_a[2])*0.5);
                reb_integrator_eos_interaction_shell1(r, dt*lf8_a[2], 0.);
                reb_integrator_eos_drift_shell1(r, dt*(lf8_a[2]+lf8_a[3])*0.5);
                reb_integrator_eos_interaction_shell1(r, dt*lf8_a[3], 0.);
                reb_integrator_eos_drift_shell1(r, dt*(lf8_a[3]+lf8_a[4])*0.5);
                reb_integrator_eos_interaction_shell1(r, dt*lf8_a[4], 0.);
                reb_integrator_eos_drift_shell1(r, dt*(lf8_a[4]+lf8_a[5])*0.5);
                reb_integrator_eos_interaction_shell1(r, dt*lf8_a[5], 0.);
                reb_integrator_eos_drift_shell1(r, dt*(lf8_a[5]+lf8_a[6])*0.5);
                reb_integrator_eos_interaction_shell1(r, dt*lf8_a[6], 0.);
                reb_integrator_eos_drift_shell1(r, dt*(lf8_a[6]+lf8_a[7])*0.5);
                reb_integrator_eos_interaction_shell1(r, dt*lf8_a[7], 0.);
                reb_integrator_eos_drift_shell1(r, dt*(lf8_a[7]+lf8_a[8])*0.5);
                reb_integrator_eos_interaction_shell1(r, dt*lf8_a[8], 0.);
                reb_integrator_eos_drift_shell1(r, dt*(lf8_a[7]+lf8_a[8])*0.5);
                reb_integrator_eos_interaction_shell1(r, dt*lf8_a[7], 0.);
                reb_integrator_eos_drift_shell1(r, dt*(lf8_a[6]+lf8_a[7])*0.5);
                reb_integrator_eos_interaction_shell1(r, dt*lf8_a[6], 0.);
                reb_integrator_eos_drift_shell1(r, dt*(lf8_a[5]+lf8_a[6])*0.5);
                reb_integrator_eos_interaction_shell1(r, dt*lf8_a[5], 0.);
                reb_integrator_eos_drift_shell1(r, dt*(lf8_a[4]+lf8_a[5])*0.5);
                reb_integrator_eos_interaction_shell1(r, dt*lf8_a[4], 0.);
                reb_integrator_eos_drift_shell1(r, dt*(lf8_a[3]+lf8_a[4])*0.5);
                reb_integrator_eos_interaction_shell1(r, dt*lf8_a[3], 0.);
                reb_integrator_eos_drift_shell1(r, dt*(lf8_a[2]+lf8_a[3])*0.5);
                reb_integrator_eos_interaction_shell1(r, dt*lf8_a[2], 0.);
                reb_integrator_eos_drift_shell1(r, dt*(lf8_a[1]+lf8_a[2])*0.5);
                reb_integrator_eos_interaction_shell1(r, dt*lf8_a[1], 0.);
                reb_integrator_eos_drift_shell1(r, dt*(lf8_a[0]+lf8_a[1])*0.5);
                reb_integrator_eos_interaction_shell1(r, dt*lf8_a[0], 0.);
                if (i<n-1){
                    reb_integrator_eos_drift_shell1(r, dt*lf8_a[0]);
                }
            }
            reb_integrator_eos_drift_shell1(r, dt*lf8_a[0]*0.5);
            break;
        case REB_EOS_LF4_2: 
            reb_integrator_eos_drift_shell1(r, dt*lf4_2_a); 
            for (int i=0;i<n;i++){
                reb_integrator_eos_interaction_shell1(r, dt*0.5, 0.); 
                reb_integrator_eos_drift_shell1(r, dt*(1.-2.*lf4_2_a));
                reb_integrator_eos_interaction_shell1(r, dt*0.5, 0.); 
                if (i<n-1){
                    reb_integrator_eos_drift_shell1(r, 2.*dt*lf4_2_a); 
                }
            }
            reb_integrator_eos_drift_shell1(r, dt*lf4_2_a); 
            break;
        case REB_EOS_LF8_6_4:
            reb_integrator_eos_drift_shell1(r, dt*lf8_6_4_a[0]);   
            for (int i=0;i<n;i++){
                reb_integrator_eos_interaction_shell1(r, lf8_6_4_b[0]*dt,0);
                reb_integrator_eos_drift_shell1(r, lf8_6_4_a[1]*dt);   
                reb_integrator_eos_interaction_shell1(r, lf8_6_4_b[1]*dt,0);
                reb_integrator_eos_drift_shell1(r, lf8_6_4_a[2]*dt);   
                reb_integrator_eos_interaction_shell1(r, lf8_6_4_b[2]*dt,0);
                reb_integrator_eos_drift_shell1(r, lf8_6_4_a[3]*dt);   
                reb_integrator_eos_interaction_shell1(r, lf8_6_4_b[3]*dt,0);
                reb_integrator_eos_drift_shell1(r, lf8_6_4_a[3]*dt);   
                reb_integrator_eos_interaction_shell1(r, lf8_6_4_b[2]*dt,0);
                reb_integrator_eos_drift_shell1(r, lf8_6_4_a[2]*dt);   
                reb_integrator_eos_interaction_shell1(r, lf8_6_4_b[1]*dt,0);
                reb_integrator_eos_drift_shell1(r, lf8_6_4_a[1]*dt);   
                reb_integrator_eos_interaction_shell1(r, lf8_6_4_b[0]*dt,0);
                if (i<n-1){
                    reb_integrator_eos_drift_shell1(r, 2.*dt*lf8_6_4_a[0]);   
                }
            }
            reb_integrator_eos_drift_shell1(r, dt*lf8_6_4_a[0]);   
            break;
        case REB_EOS_PMLF4:
            reb_integrator_eos_drift_shell1(r, dt*0.5); 
            for (int i=0;i<n;i++){
                reb_integrator_eos_interaction_shell1(r, dt, dt*dt*dt/24.); 
                if (i<n-1){
                    reb_integrator_eos_drift_shell1(r, dt);
                }
            }
            reb_integrator_eos_drift_shell1(r, dt*0.5);
            break;
        case REB_EOS_PMLF6:
            reb_integrator_eos_drift_shell1(r, dt*pmlf6_a[0]); 
            for (int i=0;i<n;i++){
                reb_integrator_eos_interaction_shell1(r, dt*pmlf6_b[0], dt*dt*dt*pmlf6_c[0]); 
                reb_integrator_eos_drift_shell1(r, dt*pmlf6_a[1]);
                reb_integrator_eos_interaction_shell1(r, dt*pmlf6_b[1], dt*dt*dt*pmlf6_c[1]); 
                reb_integrator_eos_drift_shell1(r, dt*pmlf6_a[1]);
                reb_integrator_eos_interaction_shell1(r, dt*pmlf6_b[0], dt*dt*dt*pmlf6_c[0]);
                if (i<n-1){
                    reb_integrator_eos_drift_shell1(r, 2.*dt*pmlf6_a[0]);
                }
            }
            reb_integrator_eos_drift_shell1(r, dt*pmlf6_a[0]);
            break;
        case REB_EOS_PLF7_6_4:
            reb_integrator_eos_drift_shell1(r, dt*plf7_6_4_a[0]);   
            for (int i=0;i<n;i++){
                reb_integrator_eos_interaction_shell1(r, plf7_6_4_b[0]*dt,0);
                reb_integrator_eos_drift_shell1(r, plf7_6_4_a[1]*dt);   
                reb_integrator_eos_interaction_shell1(r, plf7_6_4_b[1]*dt,0);
                reb_integrator_eos_drift_shell1(r, plf7_6_4_a[1]*dt);   
                reb_integrator_eos_interaction_shell1(r, plf7_6_4_b[0]*dt,0);
                if (i<n-1){
                    reb_integrator_eos_drift_shell1(r, 2.*dt*plf7_6_4_a[0]);   
                }
            }
            reb_integrator_eos_drift_shell1(r, dt*plf7_6_4_a[0]);   
            break;
    }
    reb_integrator_eos_postprocessor(r, dt, reos->phi1, reb_integrator_eos_drift_shell1, reb_integrator_eos_interaction_shell1);
}

void reb_integrator_eos_part1(struct reb_simulation* r){
    if (r->gravity != REB_GRAVITY_BASIC){
        reb_warning(r,"EOS only supports the BASIC gravity routine.");
    }
    // No force calculation needed between part1 and part2 of the integrator. 
    // eos_interaction() routine will set r->gravity to BASIC later. 
    r->gravity = REB_GRAVITY_NONE;

}

void reb_integrator_eos_part2(struct reb_simulation* const r){
    struct reb_simulation_integrator_eos* const reos = &(r->ri_eos);
    const double dt = r->dt;

    double dtfac = 1.;
    if (reos->is_synchronized){
        reb_integrator_eos_preprocessor(r, r->dt, reos->phi0, reb_integrator_eos_drift_shell0, reb_integrator_eos_interaction_shell0);
    }else{
        dtfac = 2.;
    }
    switch(reos->phi0){
        case REB_EOS_LF:
            reb_integrator_eos_drift_shell0(r, dt*0.5*dtfac);
            reb_integrator_eos_interaction_shell0(r, dt, 0.);
            break;
        case REB_EOS_LF4:
            reb_integrator_eos_drift_shell0(r, dt*lf4_a*dtfac);
            reb_integrator_eos_interaction_shell0(r, dt*2.*lf4_a, 0.);
            reb_integrator_eos_drift_shell0(r, dt*(0.5-lf4_a));
            reb_integrator_eos_interaction_shell0(r, dt*(1.-4.*lf4_a), 0.);
            reb_integrator_eos_drift_shell0(r, dt*(0.5-lf4_a));
            reb_integrator_eos_interaction_shell0(r, dt*2.*lf4_a, 0.);
            break;
        case REB_EOS_LF6:
            reb_integrator_eos_drift_shell0(r, dt*lf6_a[0]*0.5*dtfac);
            reb_integrator_eos_interaction_shell0(r, dt*lf6_a[0], 0.);
            reb_integrator_eos_drift_shell0(r, dt*(lf6_a[0]+lf6_a[1])*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*lf6_a[1], 0.);
            reb_integrator_eos_drift_shell0(r, dt*(lf6_a[1]+lf6_a[2])*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*lf6_a[2], 0.);
            reb_integrator_eos_drift_shell0(r, dt*(lf6_a[2]+lf6_a[3])*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*lf6_a[3], 0.);
            reb_integrator_eos_drift_shell0(r, dt*(lf6_a[3]+lf6_a[4])*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*lf6_a[4], 0.);
            reb_integrator_eos_drift_shell0(r, dt*(lf6_a[3]+lf6_a[4])*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*lf6_a[3], 0.);
            reb_integrator_eos_drift_shell0(r, dt*(lf6_a[2]+lf6_a[3])*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*lf6_a[2], 0.);
            reb_integrator_eos_drift_shell0(r, dt*(lf6_a[1]+lf6_a[2])*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*lf6_a[1], 0.);
            reb_integrator_eos_drift_shell0(r, dt*(lf6_a[0]+lf6_a[1])*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*lf6_a[0], 0.);
            break; 
        case REB_EOS_LF8: 
            reb_integrator_eos_drift_shell0(r, dt*lf8_a[0]*0.5*dtfac);
            reb_integrator_eos_interaction_shell0(r, dt*lf8_a[0], 0.);
            reb_integrator_eos_drift_shell0(r, dt*(lf8_a[0]+lf8_a[1])*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*lf8_a[1], 0.);
            reb_integrator_eos_drift_shell0(r, dt*(lf8_a[1]+lf8_a[2])*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*lf8_a[2], 0.);
            reb_integrator_eos_drift_shell0(r, dt*(lf8_a[2]+lf8_a[3])*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*lf8_a[3], 0.);
            reb_integrator_eos_drift_shell0(r, dt*(lf8_a[3]+lf8_a[4])*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*lf8_a[4], 0.);
            reb_integrator_eos_drift_shell0(r, dt*(lf8_a[4]+lf8_a[5])*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*lf8_a[5], 0.);
            reb_integrator_eos_drift_shell0(r, dt*(lf8_a[5]+lf8_a[6])*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*lf8_a[6], 0.);
            reb_integrator_eos_drift_shell0(r, dt*(lf8_a[6]+lf8_a[7])*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*lf8_a[7], 0.);
            reb_integrator_eos_drift_shell0(r, dt*(lf8_a[7]+lf8_a[8])*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*lf8_a[8], 0.);
            reb_integrator_eos_drift_shell0(r, dt*(lf8_a[7]+lf8_a[8])*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*lf8_a[7], 0.);
            reb_integrator_eos_drift_shell0(r, dt*(lf8_a[6]+lf8_a[7])*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*lf8_a[6], 0.);
            reb_integrator_eos_drift_shell0(r, dt*(lf8_a[5]+lf8_a[6])*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*lf8_a[5], 0.);
            reb_integrator_eos_drift_shell0(r, dt*(lf8_a[4]+lf8_a[5])*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*lf8_a[4], 0.);
            reb_integrator_eos_drift_shell0(r, dt*(lf8_a[3]+lf8_a[4])*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*lf8_a[3], 0.);
            reb_integrator_eos_drift_shell0(r, dt*(lf8_a[2]+lf8_a[3])*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*lf8_a[2], 0.);
            reb_integrator_eos_drift_shell0(r, dt*(lf8_a[1]+lf8_a[2])*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*lf8_a[1], 0.);
            reb_integrator_eos_drift_shell0(r, dt*(lf8_a[0]+lf8_a[1])*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*lf8_a[0], 0.);
            break;
        case REB_EOS_LF4_2: 
            reb_integrator_eos_drift_shell0(r, dt*lf4_2_a*dtfac); 
            reb_integrator_eos_interaction_shell0(r, dt*0.5, 0.); 
            reb_integrator_eos_drift_shell0(r, dt*(1.-2.*lf4_2_a));
            reb_integrator_eos_interaction_shell0(r, dt*0.5, 0.); 
            break;
        case REB_EOS_LF8_6_4:
            reb_integrator_eos_drift_shell0(r, dt*lf8_6_4_a[0]*dtfac);   
            reb_integrator_eos_interaction_shell0(r, lf8_6_4_b[0]*dt,0);
            reb_integrator_eos_drift_shell0(r, lf8_6_4_a[1]*dt);   
            reb_integrator_eos_interaction_shell0(r, lf8_6_4_b[1]*dt,0);
            reb_integrator_eos_drift_shell0(r, lf8_6_4_a[2]*dt);   
            reb_integrator_eos_interaction_shell0(r, lf8_6_4_b[2]*dt,0);
            reb_integrator_eos_drift_shell0(r, lf8_6_4_a[3]*dt);   
            reb_integrator_eos_interaction_shell0(r, lf8_6_4_b[3]*dt,0);
            reb_integrator_eos_drift_shell0(r, lf8_6_4_a[3]*dt);   
            reb_integrator_eos_interaction_shell0(r, lf8_6_4_b[2]*dt,0);
            reb_integrator_eos_drift_shell0(r, lf8_6_4_a[2]*dt);   
            reb_integrator_eos_interaction_shell0(r, lf8_6_4_b[1]*dt,0);
            reb_integrator_eos_drift_shell0(r, lf8_6_4_a[1]*dt);   
            reb_integrator_eos_interaction_shell0(r, lf8_6_4_b[0]*dt,0);
            break;
        case REB_EOS_PMLF4:
            reb_integrator_eos_drift_shell0(r, dt*0.5*dtfac); 
            reb_integrator_eos_interaction_shell0(r, dt, dt*dt*dt/24.); 
            break;
        case REB_EOS_PMLF6:
            reb_integrator_eos_drift_shell0(r, dt*pmlf6_a[0]*dtfac); 
            reb_integrator_eos_interaction_shell0(r, dt*pmlf6_b[0], dt*dt*dt*pmlf6_c[0]); 
            reb_integrator_eos_drift_shell0(r, dt*pmlf6_a[1]);
            reb_integrator_eos_interaction_shell0(r, dt*pmlf6_b[1], dt*dt*dt*pmlf6_c[1]); 
            reb_integrator_eos_drift_shell0(r, dt*pmlf6_a[1]);
            reb_integrator_eos_interaction_shell0(r, dt*pmlf6_b[0], dt*dt*dt*pmlf6_c[0]);
            break;
        case REB_EOS_PLF7_6_4:
            reb_integrator_eos_drift_shell0(r, dt*plf7_6_4_a[0]*dtfac);   
            reb_integrator_eos_interaction_shell0(r, plf7_6_4_b[0]*dt,0);
            reb_integrator_eos_drift_shell0(r, plf7_6_4_a[1]*dt);   
            reb_integrator_eos_interaction_shell0(r, plf7_6_4_b[1]*dt,0);
            reb_integrator_eos_drift_shell0(r, plf7_6_4_a[1]*dt);   
            reb_integrator_eos_interaction_shell0(r, plf7_6_4_b[0]*dt,0);
            break;
    }

    reos->is_synchronized = 0;
    if (reos->safe_mode){
        reb_integrator_eos_synchronize(r);
    }
    
    r->t+=r->dt;
    r->dt_last_done = r->dt;
    
    if (r->calculate_megno){
        r->gravity_ignore_terms = 0;
        reb_calculate_acceleration_var(r);
        double dY = r->dt * 2. * r->t * reb_tools_megno_deltad_delta(r);
        reb_tools_megno_update(r, dY);
    }

}

void reb_integrator_eos_synchronize(struct reb_simulation* r){
    struct reb_simulation_integrator_eos* const reos = &(r->ri_eos);
    const double dt = r->dt;
    if (reos->is_synchronized == 0){
        switch(reos->phi0){
            case REB_EOS_PMLF4:
            case REB_EOS_LF:
                reb_integrator_eos_drift_shell0(r, dt*0.5); 
                break;
            case REB_EOS_PMLF6:
                reb_integrator_eos_drift_shell0(r, dt*pmlf6_a[0]); 
                break;
            case REB_EOS_LF4:
                reb_integrator_eos_drift_shell0(r, dt*lf4_a);
                break;
            case REB_EOS_LF4_2:
                reb_integrator_eos_drift_shell0(r, dt*lf4_2_a); 
                break;
            case REB_EOS_PLF7_6_4:
                reb_integrator_eos_drift_shell0(r, dt*plf7_6_4_a[0]);   
                break;
            case REB_EOS_LF8_6_4:
                reb_integrator_eos_drift_shell0(r, dt*lf8_6_4_a[0]);   
                break;
            case REB_EOS_LF6:
                reb_integrator_eos_drift_shell0(r, dt*lf6_a[0]*0.5);
                break;
            case REB_EOS_LF8: 
                reb_integrator_eos_drift_shell0(r, dt*lf8_a[0]*0.5);
                break;
        }
        reb_integrator_eos_postprocessor(r, r->dt, reos->phi0, reb_integrator_eos_drift_shell0, reb_integrator_eos_interaction_shell0);
        reos->is_synchronized = 1;
    }
}

void reb_integrator_eos_reset(struct reb_simulation* r){
    r->ri_eos.n = 2;
    r->ri_eos.phi0 = REB_EOS_LF;
    r->ri_eos.phi1 = REB_EOS_LF;
    r->ri_eos.safe_mode = 1;
    r->ri_eos.is_synchronized = 1;
}

