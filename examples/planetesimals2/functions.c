//
//  functions.c
//  
//
//  Created by Ari Silburt on 2015-06-12.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include "functions.h"
#include "../../src/rebound.h"
#include "../../src/integrator_whfast.h"

void legend(char* planetdir, char* legenddir, char* removeddir, char* CEprint, struct reb_simulation* r, double tmax, double m_planetesimal, double total_planetesimal_mass, int N_planetesimals, double inner, double outer, double powerlaw, double Ms, double drh, double epsilon, int seed, int HYBRID_ON){
    
    int N_active = r->N_active, N = r->N;
    
    //system("rm -v output/orbit*.txt");
    
    char* us = "_";
    char* txt = ".txt";
    char* intgrtr;
    
    char str[110] = {0};
    if(r->integrator == REB_INTEGRATOR_WHFAST){
        if(HYBRID_ON == 1)intgrtr = "HYBRID"; else intgrtr = "WHFAST";
        strcat(str, intgrtr);
        char* teq = "_t";
        strcat(str, teq);
        char dtstr[15];
        sprintf(dtstr, "%.0f", tmax);
        strcat(str, dtstr); //planet directory
        char* Neq = "_Np";
        char Nstr[15];
        sprintf(Nstr, "%d", N_planetesimals);
        strcat(str, Neq);
        strcat(str, Nstr);
        char seedstr[15];
        char* Seq = "_sd";
        sprintf(seedstr, "%d", seed);
        strcat(str, Seq);
        strcat(str, seedstr);
        
        //Epsilon
        //char* epsln = "_Ep";
        //char eps[15];
        //sprintf(eps,"%.0e",epsilon);
        //strcat(str, epsln);
        //strcat(str, eps);
        
        //HSR and dRHill
        //strcat(str, us);
        //char str3[15];
        //sprintf(str3, "%.0f", r->ri_hybrid.switch_ratio);
        //strcat(str, str3);
        //strcat(str, us);
        //char str2[15];
        //sprintf(str2, "%.2f", drh);
        //strcat(str, str2);
        
        //char strtime[10];
        //sprintf(strtime, "%d", hybrid_rint);
        
    } else{ //pure IAS15
        intgrtr = "IAS15";
        char* teq = "_t";
        strcat(str, intgrtr);
        strcat(str, teq);
        char dtstr[15];
        sprintf(dtstr, "%.0f", tmax);
        strcat(str, dtstr); //planet directory
    }
    
    strcat(legenddir, str);
    strcat(planetdir, str);
    strcat(planetdir, txt);
    
    //particles that have been removed due to ejection or collision
    strcat(removeddir,str);
    char* err = "_removed";
    strcat(removeddir,err);
    strcat(removeddir,txt);
    
    //CE print
    strcat(CEprint,str);
    char* CEs = "_CEs";
    strcat(CEprint,CEs);
    strcat(CEprint,txt);
    
    char* file = "Properties.txt";
    strcat(legenddir, us);
    strcat(legenddir, file);
    FILE *ff;
    ff=fopen(legenddir, "w");
    fprintf(ff,"General:\ndt, tmax, Stellar mass, N_active, ri_hybrid.switch_ratio, dRHill, ias_epsilon, Seed, Integrator\n");
    fprintf(ff,"%f,%.1f,%f,%d,%f,%f,%.2e,%d,%s \n\n",r->dt,tmax,Ms,N_active,r->ri_hybrid.switch_ratio,drh,epsilon,seed,intgrtr);
    fprintf(ff,"Planetesimal:\nN_planetesimals, Mtot_planetsimal, m_planetesimal, planetesimal boundary conditions: inner/outer edge, powerlaw\n");
    fprintf(ff,"%d,%.10f,%.15f,%f,%f,%f\n\n",N-N_active,total_planetesimal_mass, m_planetesimal, inner, outer, powerlaw);
    fclose(ff);
    
    //remove any old planet files if there are
    char rmv[100] = {0};
    char* rm_v = "rm -v ";
    strcat(rmv, rm_v);
    strcat(rmv, planetdir);
    system(rmv);
    
    char rmv2[100] = {0};
    strcat(rmv2, rm_v);
    strcat(rmv2, removeddir);
    system(rmv2);
    
    char rmv3[100] = {0};
    strcat(rmv3, rm_v);
    strcat(rmv3, CEprint);
    system(rmv3);
    
}

void output_to_mercury_swifter(struct reb_simulation* r, double HSR, double tmax, int n_output, int movie_output_interval){
    struct reb_particle* restrict const particles = r->particles;
    struct reb_particle p0 = particles[0];
    int N = r->N;
    int N_active = r->N_active;
    
    //Need Hill radii for swifter too.
    FILE* swifter = fopen("swifter_mercury_output/swifter_pl.in","w");
    FILE* swifterparams = fopen("swifter_mercury_output/param.in","w");
    FILE* mercuryb = fopen("swifter_mercury_output/mercury_big.in","w");
    FILE* mercurys = fopen("swifter_mercury_output/mercury_small.in","w");
    FILE* mercuryparams = fopen("swifter_mercury_output/mercury_param.in","w");
    
    //conversion options - swifter
    int alt_units = 0;
    double mass_conv = 1, vel_conv = 1, time_conv = 1;
    if(alt_units == 1){
        mass_conv = 2.959139768995959e-04;  //solar masses to this unit
        vel_conv = 0.017202424;             //converts [v] = AU/(yr/2pi) -> AU/day
        time_conv = 58.09155423;            //converts [yr/2pi] -> days
    }
    
    //swifter initial - Nbodies and sun:
    fprintf(swifter," %d\n",N);
    fprintf(swifter," 1 %.16f\n",particles[0].m*mass_conv);
    fprintf(swifter," .0 .0 .0\n");
    fprintf(swifter," .0 .0 .0\n");
    
    //SWIFTER - heliocentric coords
    for(int i=1;i<N;i++){
        struct reb_particle p = particles[i];
        double m; if(i >= N_active) m = planetesimal_mass*mass_conv; else m = p.m*mass_conv;
        double r = sqrt((p.x-p0.x)*(p.x-p0.x) + (p.y-p0.y)*(p.y-p0.y) + (p.z-p0.z)*(p.z-p0.z));
        fprintf(swifter," %d %.16f %f\n",i+1,m,r*sqrt(Hill2[i]));
        fprintf(swifter," %f\n",p.r);
        fprintf(swifter," %.16f %.16f %.16f\n",p.x - p0.x, p.y - p0.y, p.z - p0.z);
        fprintf(swifter," %.16f %.16f %.16f\n",(p.vx - p0.vx)*vel_conv,(p.vy - p0.vy)*vel_conv,(p.vz - p0.vz)*vel_conv);
    }
    
    //SWIFTER - Other params (time, dt, etc.)
    fprintf(swifterparams,"! \n");
    fprintf(swifterparams,"! Parameter file for Swifter, with N=%d total bodies. \n",r->N);
    fprintf(swifterparams,"! \n! \n");
    fprintf(swifterparams,"T0             0.0E0 \n");
    fprintf(swifterparams,"TSTOP          %e        !In units where G=1\n",tmax);
    fprintf(swifterparams,"DT             %e        !In units where G=1\n",r->dt);
    fprintf(swifterparams,"PL_IN          swifter_pl.in\n");
    fprintf(swifterparams,"!TP_IN         tp.in     !Commented out for now, no test par\n");
    fprintf(swifterparams,"IN_TYPE        ASCII\n");
    fprintf(swifterparams,"ISTEP_OUT      %d        !# timesteps between outputs \n",movie_output_interval);
    fprintf(swifterparams,"BIN_OUT        out.dat\n");
    fprintf(swifterparams,"OUT_TYPE       REAL8\n");
    fprintf(swifterparams,"OUT_FORM       XV\n");
    fprintf(swifterparams,"OUT_STAT       NEW\n");
    fprintf(swifterparams,"ISTEP_DUMP     10000     !Dump parameters (incase of crash)\n");
    fprintf(swifterparams,"J2             0.0E0\n");
    fprintf(swifterparams,"J4             0.0E0\n");
    fprintf(swifterparams,"CHK_CLOSE      yes\n");
    fprintf(swifterparams,"CHK_RMIN       -1.0\n");
    fprintf(swifterparams,"CHK_RMAX       1000.0\n");
    fprintf(swifterparams,"CHK_EJECT      -1.0\n");
    fprintf(swifterparams,"CHK_QMIN       -1.0\n");
    fprintf(swifterparams,"!CHK_QMIN_COORD HELIO\n");
    fprintf(swifterparams,"!CHK_QMIN_RANGE 1.0 1000.0\n");
    fprintf(swifterparams,"ENC_OUT        enc.dat\n");
    fprintf(swifterparams,"EXTRA_FORCE    no\n");
    fprintf(swifterparams,"BIG_DISCARD    yes\n");
    fprintf(swifterparams,"RHILL_PRESENT  yes\n");
    
    //mercury initial:
    double day_zero = 2451179.5;
    fprintf(mercuryb,")O+_06 Big-body initial data  (WARNING: Do not delete this line!!)\n");
    fprintf(mercuryb,") Lines beginning with `)' are ignored.\n");
    fprintf(mercuryb,")---------------------------------------------------------------------\n");
    fprintf(mercuryb," style (Cartesian, Asteroidal, Cometary) = Cartesian\n");
    fprintf(mercuryb," epoch (in days) = %f\n",day_zero);
    fprintf(mercuryb,")---------------------------------------------------------------------\n");
    fprintf(mercurys,")O+_06 Big-body initial data  (WARNING: Do not delete this line!!)\n");
    fprintf(mercurys,") Lines beginning with `)' are ignored.\n");
    fprintf(mercurys,")---------------------------------------------------------------------\n");
    fprintf(mercurys," style (Cartesian, Asteroidal, Cometary) = Cartesian\n");
    fprintf(mercurys,")---------------------------------------------------------------------\n");
    
    //MERCURY - heliocentric coords
    //massive planets
    double AU_d = 0.01720242383; //converts [v] = AU/(yr/2pi) -> AU/day
    for(int i=1;i<N_active;i++){
        struct reb_particle p = particles[i];
        double r = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
        fprintf(mercuryb," BODY%d      m=%.16f r=%f\n",i,p.m,HSR*r*sqrt(Hill2[i]));
        fprintf(mercuryb," %.16f %.16f %.16f\n",p.x - p0.x,p.y - p0.y,p.z - p0.z);
        fprintf(mercuryb," %.16f %.16f %.16f\n",(p.vx - p0.vx)*AU_d,(p.vy - p0.vy)*AU_d,(p.vz - p0.vz)*AU_d);   //AU/day
        fprintf(mercuryb," 0. 0. 0.\n");
    }
    //mini bodies
    for(int i=N_active;i<N;i++){
        struct reb_particle p = particles[i];
        double r = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
        fprintf(mercurys," BODY%d      m=%.16f r=%f\n",i,planetesimal_mass,HSR*r*sqrt(Hill2[i]));
        fprintf(mercurys," %.16f %.16f %.16f\n",p.x - p0.x,p.y - p0.y,p.z - p0.z);     //AU, heliocentric
        fprintf(mercurys," %.16f %.16f %.16f\n",(p.vx - p0.vx)*AU_d,(p.vy - p0.vy)*AU_d,(p.vz - p0.vz)*AU_d);   //AU/day
        fprintf(mercurys," 0. 0. 0.\n");
    }
    
    //Mercury param file
    //int mercury_timestep = r->dt/AU_d;
    fprintf(mercuryparams,")O+_06 Big-body initial data  (WARNING: Do not delete this line!!)\n");
    fprintf(mercuryparams,") Lines beginning with `)' are ignored.\n");
    fprintf(mercuryparams,")---------------------------------------------------------------------\n");
    fprintf(mercuryparams,") Important integration parameters:\n");
    fprintf(mercuryparams,")---------------------------------------------------------------------\n");
    fprintf(mercuryparams," algorithm (MVS, BS, BS2, RADAU, HYBRID etc) = hyb\n");
    fprintf(mercuryparams," start time (days)= %f\n",day_zero);
    fprintf(mercuryparams," stop time (days) =%.1f\n",tmax/AU_d + day_zero);
    fprintf(mercuryparams," output interval (days) = %.2fd0\n",movie_output_interval*r->dt/AU_d);
    fprintf(mercuryparams," timestep (days) = %f\n",r->dt/AU_d);
    fprintf(mercuryparams," accuracy parameter=1.d-12\n");
    fprintf(mercuryparams,")---------------------------------------------------------------------\n");
    fprintf(mercuryparams,") Integration options:\n");
    fprintf(mercuryparams,")---------------------------------------------------------------------\n");
    fprintf(mercuryparams," stop integration after a close encounter = no\n");
    fprintf(mercuryparams," allow collisions to occur = yes\n");
    fprintf(mercuryparams," include collisional fragmentation = no\n");
    fprintf(mercuryparams," express time in days or years = years\n");
    fprintf(mercuryparams," express time relative to integration start time = no\n");
    fprintf(mercuryparams," output precision = medium\n");
    fprintf(mercuryparams," < not used at present >\n");
    fprintf(mercuryparams," include relativity in integration= no\n");
    fprintf(mercuryparams," include user-defined force = no\n");
    fprintf(mercuryparams,")---------------------------------------------------------------------\n");
    fprintf(mercuryparams,") These parameters do not need to be adjusted often:\n");
    fprintf(mercuryparams,")---------------------------------------------------------------------\n");
    fprintf(mercuryparams," ejection distance (AU)= 100\n");
    fprintf(mercuryparams," radius of central body (AU) = 0.005\n");
    fprintf(mercuryparams," central mass (solar) = 1.0\n");
    fprintf(mercuryparams," central J2 = 0\n");
    fprintf(mercuryparams," central J4 = 0\n");
    fprintf(mercuryparams," central J6 = 0\n");
    fprintf(mercuryparams," < not used at present >\n");
    fprintf(mercuryparams," < not used at present >\n");
    fprintf(mercuryparams," Hybrid integrator changeover (Hill radii) = 3.\n");
    fprintf(mercuryparams," number of timesteps between data dumps = 500\n");
    fprintf(mercuryparams," number of timesteps between periodic effects = 100\n");
    
    fclose(mercuryb);
    fclose(mercurys);
    fclose(swifter);
    fclose(swifterparams);
    fclose(mercuryparams);
}

double calc_dt(struct reb_simulation* r, double mp, double Ms, double a, double dRHill, double dt_prev){
    if(dRHill > r->ri_hybrid.switch_ratio){
        printf("\033[1mWarning!\033[0m dRhill !> N_RHill. Setting dRhill = N_Rhill/2 \n");
        dRHill = 0.5*r->ri_hybrid.switch_ratio;
    }
    double e_max = 0.3;  //max hypothesized eccentricity that the planet/esimals could have
    double Hill = a*(1 - e_max)*pow(mp/(3*Ms),1./3.);   //Hill radius of massive body
    double vmax = sqrt(r->G*(Ms + planetesimal_mass)*(1 + e_max)/(a*(1 - e_max)));   //perihelion speed of planetesimal
    double dt = dRHill*Hill/vmax;
    
    if(dt_prev < dt) dt = dt_prev;  //make sure I'm taking the smallest dt.
    
    return dt;
}

//Calc Hill Sphere (for speed in check_for_encounter)
void calc_Hill2(struct reb_simulation* r){
    struct reb_particle* restrict const particles = r->particles;
    struct reb_particle p0 = particles[0];
    for(int i=1;i<r->N;i++){
        struct reb_particle body = particles[i];
        double mp;
        if(i>=r->N_active) mp = planetesimal_mass; else mp = body.m;
        Hill2[i] = pow((mp/(p0.m*3.)), 2./3.);
    }
}

double calc_Etot(struct reb_simulation* a, double soft, double dE_collision){
    double m1,m2;
    const int N = a->N;
    const int N_active = a->N_active;
    const double G = a->G;
    double U = 0, K = 0;
    struct reb_particle* const particles = a->particles;
    for(int i=0;i<N;i++){
        struct reb_particle par = particles[i];
        if(i < N_active) m1 = par.m; else m1 = planetesimal_mass;
        const double dvx = par.vx;
        const double dvy = par.vy;
        const double dvz = par.vz;
        const double dx = par.x;
        const double dy = par.y;
        const double dz = par.z;
        
        //L_tot = m*(r x v)
        //const double hx = dy*dvz - dz*dvy;
        //const double hy = dz*dvx - dx*dvz;
        //const double hz = dx*dvy - dy*dvx;
        //L += m1*sqrt ( hx*hx + hy*hy + hz*hz );
        
        //E_tot
        K += 0.5*m1*(dvx*dvx + dvy*dvy + dvz*dvz);
        
        if(i<N_active){//ignore dE/dx = forces between planetesimals
            for(int j=i+1;j<N;j++){
                struct reb_particle par2 = particles[j];
                if(j < N_active) m2 = par2.m; else m2 = planetesimal_mass;
                double ddx = dx - par2.x;
                double ddy = dy - par2.y;
                double ddz = dz - par2.z;
                U -= G*m1*m2/sqrt(ddx*ddx + ddy*ddy + ddz*ddz + soft*soft);
            }
        }
    }
    
    return K + U + dE_collision;
}

//Calculates 'a' and 'e' of planet each output.
void calc_ae(double* a, double* e, double* d_out, struct reb_simulation* r, int i, double t){
    struct reb_particle* const particles = r->particles;
    struct reb_particle com = reb_get_com(r);
    struct reb_particle* par = &(particles[i]); //output planets only.
    const double G = r->G;
    const double m = par->m;
    const double mu = G*(com.m + m);
    const double dvx = par->vx-com.vx;
    const double dvy = par->vy-com.vy;
    const double dvz = par->vz-com.vz;
    const double dx = par->x-com.x;
    const double dy = par->y-com.y;
    const double dz = par->z-com.z;
    
    const double vv = dvx*dvx + dvy*dvy + dvz*dvz;
    const double d = sqrt( dx*dx + dy*dy + dz*dz );    //distance
    const double dinv = 1./d;
    const double muinv = 1./mu;
    const double vr = (dx*dvx + dy*dvy + dz*dvz)*dinv;
    const double term1 = vv-mu*dinv;
    const double term2 = d*vr;
    const double ex = muinv*( term1*dx - term2*dvx );
    const double ey = muinv*( term1*dy - term2*dvy );
    const double ez = muinv*( term1*dz - term2*dvz );
    *e = sqrt(ex*ex + ey*ey + ez*ez);   // eccentricity
    *a = -mu/( vv - 2.*mu*dinv );
    *d_out = d;
    
}

void planetesimal_forces_global(struct reb_simulation *r){
    const double G = r->G;
    const int N = r->N;
    const int N_active = r->N_active;
    struct reb_particle* const particles = r->particles;
    const double Gm1 = G*planetesimal_mass;
    
    for(int j=N_active;j<N;j++){//add planetesimal forces to massive bodies
        struct reb_particle p = particles[j];
        _Bool skip = 0;
        for(int k=0;k<N_encounters_previous && skip == 0;k++){
            if(p.id == previous_encounter_index[k]){
                skip++;
            }
        }
            
        if(skip == 0){//If CE don't calculate planetesimal forces in global
            for(int i=0;i<N_active;i++){
                struct reb_particle* body = &(particles[i]);
                const double dx = body->x - p.x;
                const double dy = body->y - p.y;
                const double dz = body->z - p.z;
                
                const double rijinv2 = 1.0/(dx*dx + dy*dy + dz*dz + soft*soft);
                const double ac = -Gm1*rijinv2*sqrt(rijinv2);
                
                body->ax += ac*dx;      //perturbation on planets due to planetesimals.
                body->ay += ac*dy;
                body->az += ac*dz;
            }
        }
    }
}

void planetesimal_forces_mini(struct reb_simulation *s){
    const double G = s->G;
    const int N = s->N;
    const int N_active = s->N_active;
    struct reb_particle* mini = s->particles;
    
    const double Gm1 = G*planetesimal_mass;
    for(int i=0;i<N_active;i++){
        struct reb_particle* body = &(mini[i]);
        for(int j=N_active;j<N;j++){//add planetesimal forces to massive bodies
            struct reb_particle p = mini[j];
            
            const double dx = body->x - p.x;
            const double dy = body->y - p.y;
            const double dz = body->z - p.z;
            
            const double rijinv2 = 1.0/(dx*dx + dy*dy + dz*dz + soft*soft);
            const double ac = -Gm1*rijinv2*sqrt(rijinv2);
            
            body->ax += ac*dx;      //perturbation on planets due to planetesimals.
            body->ay += ac*dy;
            body->az += ac*dz;
        }
    }
    
    //forces from global into mini
    struct reb_particle* const global = r->particles;
    const double timefac = (s->t - t_prev)/(r->t - t_prev);
    int rN_active = r->N_active;
    for(int i=rN_active;i<r->N;i++){    //planetesimals
        if(x_prev[i] != 0){             //find planetesimals with !=0 values, i.e. part of global but not mini
            const double ix = x_prev[i] + timefac*(global[i].x - x_prev[i]); //interpolated values
            const double iy = y_prev[i] + timefac*(global[i].y - y_prev[i]);
            const double iz = z_prev[i] + timefac*(global[i].z - z_prev[i]);
            for(int j=0;j<rN_active;j++){//massive bodies
                struct reb_particle* body = &(mini[j]);
                const double ddx = body->x - ix;
                const double ddy = body->y - iy;
                const double ddz = body->z - iz;
                
                const double rijinv2 = 1.0/(ddx*ddx + ddy*ddy + ddz*ddz + soft*soft);
                const double ac = -Gm1*rijinv2*sqrt(rijinv2);
                
                body->ax += ac*ddx;     //perturbation on planets due to planetesimals.
                body->ay += ac*ddy;
                body->az += ac*ddz;
            }
        }
    }
}

//initialize mini-simulation for close encounters
void ini_mini(struct reb_simulation* const r, struct reb_simulation* s, double ias_epsilon, int turn_planetesimal_forces_on, double ias_timestep, double soft){
    s->N_active = r->N_active;
    s->integrator = REB_INTEGRATOR_IAS15;
    if(turn_planetesimal_forces_on==1)s->additional_forces = planetesimal_forces_mini;
    s->exact_finish_time = 1;
    s->ri_ias15.epsilon = ias_epsilon;
    s->dt = ias_timestep;
    s->softening = soft;
    
    struct reb_particle* restrict const particles = r->particles;
    for(int k=0; k<s->N_active; k++){
        struct reb_particle p = particles[k];
        reb_add(s,p);
    }
}

//collect the id/array number of all planetesimals involved in a close encounter
void check_for_encounter(struct reb_simulation* r, struct reb_simulation* s, int* N_encounters, int N_encounters_previous, double* min_r, double* max_val, char* removeddir, int* output_it, double* dE_collision, double soft, double ejection_distance2){
    const int rN = r->N;
    const int rN_active = r->N_active;
    struct reb_particle* global = r->particles;
    struct reb_particle p0 = global[0];
    int num_encounters = 0;
    for (int i=0; i<rN_active; i++){
        struct reb_particle* body = &(global[i]);
        const double dxi = p0.x - body->x;
        const double dyi = p0.y - body->y;
        const double dzi = p0.z - body->z;
        const double r0i2 = dxi*dxi + dyi*dyi + dzi*dzi;
        const double rhi = r0i2*Hill2[i];
        
        for (int j=i+1; j<rN; j++){
            struct reb_particle pj = global[j];
            double HSR = r->ri_hybrid.switch_ratio;
            
            _Bool found_in_mini = 0;
            for(int k=0; k<N_encounters_previous && found_in_mini == 0; k++){
                if(global[j].id == previous_encounter_index[k]){
                    HSR *= 1.05;    //HSR(mini) is a bit bigger so no constant enter/leave
                    found_in_mini = 1;
                }
            }
            
            const double dxj = p0.x - pj.x;
            const double dyj = p0.y - pj.y;
            const double dzj = p0.z - pj.z;
            const double r0j2 = dxj*dxj + dyj*dyj + dzj*dzj;
            const double rhj = r0j2*Hill2[j];
            
            const double dx = body->x - pj.x;
            const double dy = body->y - pj.y;
            const double dz = body->z - pj.z;
            const double rij2 = dx*dx + dy*dy + dz*dz;
            const double ratio = rij2/(rhi+rhj);    //(p-p distance/Hill radii)^2
            
            if(ratio < HSR){
                double radius2 = body->r*body->r;
                if(rij2 < 10000*radius2){//Super close encounter, check for collision by predicting the future step
                    double dx1 = (pj.vx - body->vx)*r->dt; //xf - xi = distance travelled in dt relative to body
                    double dy1 = (pj.vy - body->vy)*r->dt;
                    double dz1 = (pj.vz - body->vz)*r->dt;
                    double dx2 = pj.x - body->x;
                    double dy2 = pj.y - body->y;
                    double dz2 = pj.z - body->z;
                    double x = dy1*dz2 - dz1*dy2;
                    double y = dz1*dx2 - dx1*dz2;
                    double z = dx1*dy2 - dy1*dx2;
                    double d2 = (x*x + y*y + z*z)/(dx1*dx1 + dy1*dy1 + dz1*dz1);
                    //if(pj.id == 159) printf("\nParticle 159: t=%f,d=%f,rij=%f,radius=%f\n",r->t,sqrt(d2),sqrt(rij2),body->r);
                    //if(pj.id == 159) printf("\nParticle 159: t=%f,predictxyz=%f,%f,%f, xyz=%f,%f,%f\n",r->t,dx1+pj.x,dy1+pj.y,dz1+pj.z,pj.x,pj.y,pj.z);
                    if(d2 < radius2){//Collision will happen next time step.
                         if(j < rN_active){
                         fprintf(stderr,"\n\033[1mCollision at t=%f between %d and %d, both are Massive bodies. Can't deal with this collisional physics right now. Exiting. \033[0m \n",r->t,pj.id,body->id);
                         exit(0);
                         }
                         double massive_mass = body->m;
                         double invmass = 1.0/(massive_mass + planetesimal_mass);
                         double E_i = calc_Etot(r, soft, 0);
                         
                         body->vx = (body->vx*massive_mass + pj.vx*planetesimal_mass)*invmass;
                         body->vy = (body->vy*massive_mass + pj.vy*planetesimal_mass)*invmass;
                         body->vz = (body->vz*massive_mass + pj.vz*planetesimal_mass)*invmass;
                         body->m += planetesimal_mass;
                         struct reb_particle* mini = s->particles;
                         mini[i] = *body;     //need to update mini accordingly
                         
                         fprintf(stderr,"\n\033[1mCollision at t=%.16f!\033[0m between Particle %d and Planet %d, r=%f, planet radius=%f.\n",r->t,pj.id,body->id,sqrt(rij2),sqrt(radius2));
                         FILE* ff;
                         ff = fopen(removeddir,"a");
                         fprintf(ff,"Collision at t=%f between Particle %d and Planet %d, r=%f.\n",r->t,pj.id,body->id,sqrt(rij2));
                         fclose(ff);
                         *output_it = 1;
                         
                         reb_remove(r,j,1);
                         
                         double E_f = calc_Etot(r, soft, 0);
                         *dE_collision += E_i - E_f;
                         
                         //Update Hill radii and xyz_prev arrays
                         Hill2[i] = pow((body->m/(p0.m*3.)), 2./3.);
                         for(int k=j;k<rN-1;k++){
                         x_prev[k] = x_prev[k+1];
                         y_prev[k] = y_prev[k+1];
                         z_prev[k] = z_prev[k+1];
                         Hill2[k] = Hill2[k+1];
                         }
                         Hill2 = realloc(Hill2,(rN-1)*sizeof(double));
                         x_prev = realloc(x_prev,(rN-1)*sizeof(double));
                         y_prev = realloc(y_prev,(rN-1)*sizeof(double));
                         z_prev = realloc(z_prev,(rN-1)*sizeof(double));
                    }
                } else {//add to CE array
                    num_encounters++;
                    if(num_encounters == 1) encounter_index[0] = pj.id;
                    else if(num_encounters > 1){//multiple close encounters
                        encounter_index = realloc(encounter_index,num_encounters*sizeof(int));
                        encounter_index[num_encounters - 1] = pj.id;
                    }
                }
            }
            
            if(i==0 && rij2 > ejection_distance2){//Ejection
                FILE* ff;
                ff = fopen(removeddir,"a");
                fprintf(ff,"Ejection at t=%f for particle %d, r=%f.\n",r->t,pj.id,sqrt(rij2));
                fclose(ff);
                fprintf(stderr,"\n\033[1mEjected Particle %d at t=%f!\033[0m Particle too far from sun r=%f.\n",pj.id,r->t,sqrt(rij2));
                *output_it = 1;
                
                double E_i = calc_Etot(r, soft, 0);
                reb_remove(r,j,1);
                double E_f = calc_Etot(r, soft, 0);
                *dE_collision += E_i - E_f;
            }
            
            //calculate dt*(vrel/rmin)
            double vx = body->vx - pj.vx;
            double vy = body->vy - pj.vy;
            double vz = body->vz - pj.vz;
            double vrel = sqrt(vx*vx + vy*vy + vz*vz);
            double rr = sqrt(rij2);
            double val = r->dt*vrel/rr;
            if(rr < *min_r) *min_r = rr;
            if(val > *max_val) *max_val = val;
        }
    }
    *N_encounters = num_encounters;
}
 
//Just after mini has been integrated up to r->t, update global.
void update_global(struct reb_simulation* const s, struct reb_simulation* r, int N_encounters_previous){
    int N_active = s->N_active;
    struct reb_particle* global = r->particles;
    struct reb_particle* const mini = s->particles;
    
    //update massive and planetesimal particles
    for(int i=0; i<N_active; i++) global[i] = mini[i];  //massive particles, always in same order
    for(int j=0; j<N_encounters_previous; j++){
        int PEI = previous_encounter_index[j];          //encounter index == global[EI].id
        _Bool found_mini = 0;
        _Bool particle_update = 0;
        int mini_index;
        for(int k=N_active;found_mini==0 && k<s->N;k++){
            if(mini[k].id == PEI){ mini_index = k; found_mini = 1; }
        }
        for(int k=N_active;particle_update==0 && k<r->N;k++){
            if(global[k].id == PEI && found_mini == 1){
                particle_update = 1;
                global[k] = mini[mini_index];
            }
        }
        
        if(particle_update == 0){
            fprintf(stderr,"\n\033[1mAlert!\033[0m Particle %d couldn't be found in update_global. Exiting. \n",PEI);
            printf("N_encounters_previous=%d,size=%lu,",N_encounters_previous,sizeof(previous_encounter_index)/sizeof(previous_encounter_index[0]));
            for(int i=0;i<N_encounters_previous;i++) printf("PEI(%d)=%d,",i,previous_encounter_index[i]);
            printf("\n");
            printf("Mini:");
            for(int i=0;i<N_encounters_previous;i++) printf("mini[%d].id=%d,",i,mini[N_active+i].id);
            printf("\n");
            exit(0);
        }
    }
    
}

void add_or_subtract_particles(struct reb_simulation* r, struct reb_simulation* s, int N_encounters, int N_encounters_previous, char* CEprint, double soft, double dE_collision, double E0){
    int N_active = s->N_active;
    struct reb_particle* mini = s->particles;
    struct reb_particle* global = r->particles;
    int dN = N_encounters - N_encounters_previous;
    
    //check to add particles
    for(int i=0;i<N_encounters;i++){
        _Bool index_found = 0;
        int EI = encounter_index[i];
        if(EI >= N_active){//don't want to add/remove massive bodies. Already in mini
            for(int j=0;index_found == 0 && j<N_encounters_previous;j++){
                if(EI == previous_encounter_index[j]) index_found = 1;
            }
            if(index_found == 0){//couldn't find index
                _Bool added_particle = 0;
                for(int k=r->N_active && added_particle == 0;k<r->N;k++){
                    if(global[k].id == EI){
                        struct reb_particle pt = global[k];
                        reb_add(s,pt);
                        N_encounters_tot++;
                        added_particle = 1;
                        
                        double E1 = calc_Etot(r, soft, dE_collision);
                        FILE *output;
                        output = fopen(CEprint, "a");
                        fprintf(output,"t=%f,%f particle %d added. dN == %d, N_close_encounters=%d, dE/E=%.16f\n",r->t,s->t,EI,dN,N_encounters,fabs((E1 - E0)/E0));
                        for(int i=0;i<N_encounters;i++)fprintf(output,"EI[%d]=%d,",i,encounter_index[i]);
                        fprintf(output,"\n");
                        for(int i=0;i<N_encounters_previous;i++)fprintf(output,"PEI[%d]=%d,",i,previous_encounter_index[i]);
                        fprintf(output,"\n");
                        fclose(output);
                    }
                }
            }
        }
    }
    
    //check to remove particles
    for(int i=0;i<N_encounters_previous;i++){
        _Bool index_found = 0;
        int PEI = previous_encounter_index[i];
        if(PEI >= N_active){//don't want to add/remove massive bodies. Already in mini
            for(int j=0;index_found == 0 && j<N_encounters;j++){
                if(PEI == encounter_index[j]) index_found = 1;
            }
            if(index_found == 0){//couldn't find index, remove particle
                int removed_particle = 0;
                for(int k=N_active;removed_particle==0 && k<s->N;k++){
                    if(mini[k].id == PEI){
                        removed_particle = reb_remove(s,k,1);    //remove particle
                        
                        double E1 = calc_Etot(r, soft, dE_collision);
                        FILE *output;
                        output = fopen(CEprint, "a");
                        fprintf(output,"t=%f,%f particle %d leaving. dN == %d, N_close_encounters=%d, dE/E=%.16f.\n",r->t,s->t,PEI,dN,N_encounters,fabs((E1 - E0)/E0));
                        for(int i=0;i<N_encounters;i++)fprintf(output,"EI[%d]=%d,",i,encounter_index[i]);
                        fprintf(output,"\n");
                        for(int i=0;i<N_encounters_previous;i++)fprintf(output,"PEI[%d]=%d,",i,previous_encounter_index[i]);
                        fprintf(output,"\n");
                        fclose(output);
                    }
                }
            }
        }
    }
}

void update_previous_global_positions(struct reb_simulation* r, int N_encounters){
    struct reb_particle* global = r->particles;
    t_prev = r->t;
    for(int i=r->N_active;i<r->N;i++){
        int ID = global[i].id;
        _Bool found_particle = 0;
        for(int j=0;j<N_encounters && found_particle == 0;j++){
            if(ID == encounter_index[j]){
                found_particle = 1;
                x_prev[i] = 0.;         //reset planetesimals involved in mini to 0
                y_prev[i] = 0.;
                z_prev[i] = 0.;
            }
        }
        if(found_particle == 0){        //planetesimal not involved in mini, update position for later interp.
            x_prev[i] = global[i].x;
            y_prev[i] = global[i].y;
            z_prev[i] = global[i].z;
        }
    }
}

//transfer values from encounter_index to previous_encounter_index
void update_encounter_indices(int* N_encounters, int* N_encounters_previous){
    int size;
    if(*N_encounters == 0) size = 1; else size = *N_encounters;
    previous_encounter_index = realloc(previous_encounter_index,size*sizeof(int));
    if(*N_encounters == 0) previous_encounter_index[0] = 0; else {
        for(int i=0;i<*N_encounters;i++) previous_encounter_index[i] = encounter_index[i];
    }
    
    //reset encounter index
    encounter_index = realloc(encounter_index,sizeof(int));   //reset to single element
    encounter_index[0] = 0;
    
    //reset encounter counters.
    *N_encounters_previous = *N_encounters;
    *N_encounters = 0;
}

void output_frame_per_body(struct reb_particle* particles, char* dir, int N, double t){
    struct reb_particle p0 = particles[0];
    for(int i=0;i<N;i++){
        char str[50] = {0};
        char temp[7];
        strcat(str, dir);
        sprintf(temp, "%d",particles[i].id);
        strcat(str,temp);
        strcat(str,".txt");
        FILE *output;
        output = fopen(str,"a");
        fprintf(output, "%.12f,%d,%.16f,%.16f,%.16f\n",t,particles[i].id,particles[i].x-p0.x,particles[i].y-p0.y,particles[i].z-p0.z);
        fclose(output);
    }
}

void output_frame_per_time(struct reb_particle* particles, char* name, int N, double t, int* movie_counter){
        char str[50] = {0};
        char temp[7];
        strcat(str, name);
        int mc = *movie_counter;
        sprintf(temp, "%d",mc);
        strcat(str,temp);
        strcat(str,".txt");
        FILE *output;
        output = fopen(str,"w");
        for(int i=0;i<N;i++) fprintf(output, "%f,%d,%.16f,%.16f,%.16f\n",t,particles[i].id,particles[i].x,particles[i].y,particles[i].z);
        fclose(output);
        mc++;
        *movie_counter = mc;
}

time_t clock_start(){
    char buf[64];
    time_t t_ini = time(NULL);
    struct tm *tmp = gmtime(&t_ini);
    strftime(buf, sizeof(buf), "%j:%H:%M:%S\n", tmp);
    printf("start time (GMT): %s\n",buf);
    
    return t_ini;
}

void clock_finish(clock_t t_ini, int N_encounters, int N, char* legenddir){
    char buf[64];
    time_t t_fini = time(NULL);
    struct tm *tmp = gmtime(&t_fini);
    strftime(buf, sizeof(buf), "%j:%H:%M:%S\n", tmp);
    printf("\nfinish time (GMT): %s\n",buf);
    
    double time = t_fini - t_ini;
    
    FILE *ff;
    ff=fopen(legenddir, "a");
    fprintf(ff,"Elapsed simulation time is %.2f s, with %d close encounters and %d planetesimals remaining.\n",time,N_encounters,N);
    printf("\nSimulation complete. Elapsed simulation time is %.2f s, with %d close encounters and %d planetesimals remaining.\n\n",time,N_encounters, N);
}

void global_free(){
    free(encounter_index);
    free(previous_encounter_index);
    free(Hill2);
    free(x_prev);
    free(y_prev);
    free(z_prev);
}
