#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rebound.h"


double min_frag_mass = 1.4e-8;
int tot_no_frags = 0;  //if restarting a simulation this needs to be changed to the last number of frags in the simulation, otherwise new fragments added will rewrite exisiting frags
int trace_close_encounters = 0;
int mercurius_close_encounters = 0;


#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< Returns the minimum of a and b
#define MAX(a, b) ((a) > (b) ? (a) : (b))    // Returns the maximum of a and b


struct collision_params
{
    int target;
    int projectile;
    double dx;
    double dy;
    double dz;
    double b;
    double Vix;
    double Viy;
    double Viz;
    double Vi;
    double l;
    double rho1;
    double cstar;
    double mu;
    double QR;
    double QpRD;
    double V_esc;
    double separation_distance;
    double Mlr;
    double Mslr;
    double Q;
    double Mlr_dag;
    double Q_star;
    double vrel;
    double xrel;
    int collision_type;
    int no_frags;
}; 


void make_vector(double x1, double y1, double z1, double x2, double y2, double z2, double *x, double*y, double*z){   //Galilean transform
    *x = x1-x2;
    *y = y1-y2;
    *z = z1-z2;
}

double get_dot(double x1, double y1, double z1, double x2, double y2, double z2){ 
    return (x1*x2)+(y1*y2)+(z1*z2);
}  //return dot product of two vectors

double get_mag(double x, double y, double z){
    return sqrt(pow(x,2)+pow(y,2)+pow(z,2));
}   //return magnitude of vector
double get_radii(double m, double rho){
    return pow((3*m)/(4*M_PI*rho),1./3.);
}

void add_fragments(struct reb_simulation* const r, struct reb_collision c, struct collision_params *params){
    struct reb_particle* target = &(r->particles[params->target]);
    struct reb_particle* projectile = &(r->particles[params->projectile]);
    struct reb_particle com = reb_particle_com_of_pair(*target, *projectile);
    double initial_mass = target -> m + projectile -> m;
    double remaining_mass = initial_mass - params->Mlr;
    double rho = target->m/(4./3*M_PI*pow(target ->r, 3));
    double rtot = target -> r + projectile -> r;

    int big_frags = 0;
    if (params->Mslr > 0){
        remaining_mass = remaining_mass -  params->Mslr;
        big_frags = 1;
    }

    int no_frags = remaining_mass/min_frag_mass;  //fragments are broken up into equal sizes
    double frag_mass = remaining_mass/no_frags;


    int new_bodies = no_frags + big_frags;
    params->no_frags = new_bodies;


    char hash[10];
    double mxsum[3] = {0,0,0};
    double mvsum[3] = {0,0,0};
    //target gets mass of Mlr and is assigned COM position and velocity;
    target -> last_collision = r->t;
    target -> m = params->Mlr;
    target -> r = get_radii(params->Mlr, rho);
    target->x = com.x;
    target->y = com.y;
    target->z = com.z;

    target->vx = com.vx;
    target->vy = com.vy;
    target->vz = com.vz;

    if (no_frags == 1 && params->Mlr <= frag_mass){
        target->m = frag_mass;
        frag_mass = params->Mlr;
    }

    mxsum[0] = mxsum[0] + target->m*target->x;
    mxsum[1] = mxsum[1] + target->m*target->y;  
    mxsum[2] = mxsum[2] + target->m*target->z;

    mvsum[0] = mvsum[0] + target->m*target->vx;
    mvsum[1] = mvsum[1] + target->m*target->vy; 
    mvsum[2] = mvsum[2] + target->m*target->vz;         

    double theta_inc = (2.*M_PI)/new_bodies;


    double unit_vix, unit_viy, unit_viz, zx, zy, zz, z, ox, oy, oz, o;

    unit_vix = params->Vix/params->vrel;  //unit vector parallel to target velocity
    unit_viy = params->Viy/params->vrel;
    unit_viz = params->Viz/params->vrel;

    zx = (params->Viy*params->dz - params->Viz*params->dy);                     // vector normal to the collision plane; vrel cross xrel
    zy = (params->Viz*params->dx - params->Vix*params->dz);
    zz = (params->Vix*params->dy - params->Viy*params->dx);

    z = get_mag(zx, zy, zz);

    zx = zx/z;          //unit vector
    zy = zy/z;
    zz = zz/z;


    ox = (zy*params->Viz - zz*params->Viy);                   // vector normal to target velocity in collision plane; z cross vrel
    oy = (zz*params->Vix - zx*params->Viz);
    oz = (zx*params->Viy - zy*params->Vix);

    o = get_mag(ox, oy, oz);

    ox = ox/o;      //unit vector
    oy = oy/o;
    oz = oz/o;

    double fragment_velocity =sqrt(1.1*pow(params->V_esc,2) - 2*r->G*initial_mass*(1./rtot - 1./params->separation_distance));

    if (big_frags == 1){  //assign radii, positions and velocities to second largest remnant, theta=0
        struct reb_particle Slr1 = {0};
        Slr1.m = params->Mslr;
        Slr1.x = com.x + params->separation_distance*unit_vix;
        Slr1.y = com.y + params->separation_distance*unit_viy;
        Slr1.z = com.z + params->separation_distance*unit_viz;

        Slr1.vx = com.vx + fragment_velocity*unit_vix;
        Slr1.vy = com.vy + fragment_velocity*unit_viy;
        Slr1.vz = com.vz + fragment_velocity*unit_viz;

        Slr1.r = get_radii(Slr1.m, rho);
        sprintf(hash,"FRAG%d", tot_no_frags+1);
        Slr1.hash = reb_hash(hash);
        //printf("%s hash, mass:      %u %e\n", hash, Slr1.hash, Slr1.m);
        mxsum[0] += Slr1.m*Slr1.x;
        mxsum[1] += Slr1.m*Slr1.y;    
        mxsum[2] += Slr1.m*Slr1.z;

        mvsum[0] += Slr1.m*Slr1.vx;
        mvsum[1] += Slr1.m*Slr1.vy;   
        mvsum[2] += Slr1.m*Slr1.vz;
        Slr1.last_collision = r->t;
        reb_simulation_add(r, Slr1);
    }



    int new_beginning_frag_index = tot_no_frags+big_frags+1;
    for (int i=(new_beginning_frag_index); i<(new_beginning_frag_index+no_frags); i++){          //add fragments
        struct reb_particle fragment = {0};
        int j = i - new_beginning_frag_index+1;
        fragment.m = frag_mass;                  
        fragment.x = com.x + params->separation_distance*(cos(theta_inc*j)*unit_vix + sin(theta_inc*j)*ox);
        fragment.y = com.y + params->separation_distance*(cos(theta_inc*j)*unit_viy + sin(theta_inc*j)*oy);
        fragment.z = com.z + params->separation_distance*(cos(theta_inc*j)*unit_viz + sin(theta_inc*j)*oz);
        fragment.vx = com.vx + fragment_velocity*(cos(theta_inc*j)*unit_vix + sin(theta_inc*j)*ox);
        fragment.vy = com.vy + fragment_velocity*(cos(theta_inc*j)*unit_viy + sin(theta_inc*j)*oy);
        fragment.vz = com.vz + fragment_velocity*(cos(theta_inc*j)*unit_viz + sin(theta_inc*j)*oz);

        fragment.r = get_radii(frag_mass, rho);
        fragment.last_collision = r->t;
        sprintf(hash, "FRAG%d", i);
        fragment.hash = reb_hash(hash);
        //printf("%s hash, mass:      %u %e\n", hash, fragment.hash, fragment.m);
        mxsum[0] +=fragment.m*fragment.x;
        mxsum[1] += fragment.m*fragment.y;    
        mxsum[2] += fragment.m*fragment.z;

        mvsum[0] += fragment.m*fragment.vx;
        mvsum[1] += fragment.m*fragment.vy;    
        mvsum[2] += fragment.m*fragment.vz;

        reb_simulation_add(r, fragment); 
    }
    tot_no_frags += big_frags+no_frags;


    //Ensure momentum is conserved



    double xoff[3] = {com.x - mxsum[0]/initial_mass, com.y - mxsum[1]/initial_mass, com.z - mxsum[2]/initial_mass};
    double voff[3] = {com.vx - mvsum[0]/initial_mass, com.vy - mvsum[1]/initial_mass, com.vz - mvsum[2]/initial_mass};


    target -> x +=  xoff[0]*target->m/initial_mass; 
    target -> y += xoff[1]*target->m/initial_mass; 
    target -> z += xoff[2]*target->m/initial_mass; 
    target -> vx += voff[0]*target->m/initial_mass; 
    target -> vy += voff[1]*target->m/initial_mass; 
    target -> vz += voff[2]*target->m/initial_mass; 

    for (int i=(tot_no_frags-new_bodies)+1; i<(tot_no_frags+1); i++){ 
        char frag[10];
        sprintf(frag, "FRAG%d", i);
        double mass_fraction = reb_simulation_particle_by_hash(r, reb_hash(frag))->m/initial_mass;
        reb_simulation_particle_by_hash(r, reb_hash(frag))->x += xoff[0]*mass_fraction;
        reb_simulation_particle_by_hash(r, reb_hash(frag))->y += xoff[1]*mass_fraction;
        reb_simulation_particle_by_hash(r, reb_hash(frag))->z += xoff[2]*mass_fraction;

        reb_simulation_particle_by_hash(r, reb_hash(frag))->vx += voff[0]*mass_fraction;
        reb_simulation_particle_by_hash(r, reb_hash(frag))->vy += voff[1]*mass_fraction;
        reb_simulation_particle_by_hash(r, reb_hash(frag))->vz += voff[2]*mass_fraction;
    }

    return;
}


void merge(struct reb_simulation* const r, struct reb_collision c, struct collision_params *params){
    struct reb_particle* pi = &(r->particles[params->target]);
    struct reb_particle* pj = &(r->particles[params->projectile]);

    double invmass = 1.0/(pi->m + pj->m);
    double targ_rho = pi->m/(4./3*M_PI*pow(pi->r,3));  //new body recieves density of the target
                                                       // Merge by conserving mass, volume and momentum
    pi->vx = (pi->vx*pi->m + pj->vx*pj->m)*invmass;
    pi->vy = (pi->vy*pi->m + pj->vy*pj->m)*invmass;
    pi->vz = (pi->vz*pi->m + pj->vz*pj->m)*invmass;
    pi->x  = (pi->x*pi->m + pj->x*pj->m)*invmass;
    pi->y  = (pi->y*pi->m + pj->y*pj->m)*invmass;
    pi->z  = (pi->z*pi->m + pj->z*pj->m)*invmass;
    pi->m  = pi->m + pj->m;
    pi->r  = pow((3*pi->m)/(4*M_PI*targ_rho),1./3.);
    pi->last_collision = r->t;


    return; // 
}

int hit_and_run(struct reb_simulation* const r, struct reb_collision c, struct collision_params *params){  //also includes partial accretion.  Mlr = M_target.  Projectile is erroded.
    struct reb_particle* target = &(r->particles[params->target]);
    struct reb_particle* projectile = &(r->particles[params->projectile]);


    int swap = 2;
    int i = c.p1;
    int j = c.p2;   //make sure projectile is the particle being removed
    struct reb_particle* pi = &(r->particles[i]);
    struct reb_particle* pj = &(r->particles[j]);
    if (pi->m < pj->m){
        swap = 1;
    }

    double phi = 2*acos((params->l-projectile->r)/projectile->r);
    double A_interact = pow(projectile->r, 2)*((M_PI-(phi-sin(phi))/2.));  //Leinhardt Eq. 46;
    double L_interact = 2.*pow(pow(target->r,2)-(pow(target->r-params->l/2.,2)), .5);   //Leinhardt Eq. 47
    double beta = (A_interact*L_interact)/target->m;  //Chambers Eq. 11
    double Rc1 = pow(3./(4.*M_PI*params->rho1)*(beta*target->m + projectile->m), 1./3.); 
    double Q0 = .8*params->cstar*M_PI*params->rho1*r->G*pow(Rc1, 2);
    double gamma = (beta*target->m)/projectile->m;
    double Q_star = (pow(1+gamma, 2)/4*gamma)* Q0;

    double mu = (beta*target->m*projectile->m)/(beta*target->m+projectile->m);  //Chambers Eq. 13
    double Q = .5*(mu*pow(params->Vi,2))/(beta*target->m+projectile->m); //Chambers Eq. 12

    double c1 = 2.43;
    double c2 = -0.0408;
    double c3 = 1.86;
    double c4 = 1.08;

    double targ_m = target->m;
    double imp_m = projectile->m;
    double zeta = pow((targ_m - imp_m)/(targ_m + imp_m),2);
    double fac = pow(1-params->b/(target->r + projectile->r),2.5);
    double v_crit = params->V_esc*(c1*zeta*fac + c2*zeta +c3*fac + c4);

    if (params->Vi <= v_crit){             //if impact velocity is low, the hit-and-run results in a merger.
                                           //printf("GRAZE AND MERGE\n");  
        params->collision_type = 1;          
        merge(r,c,params);
        return swap;
    } else{ //vi>v_crit
        if (params->Mlr<targ_m){ //Target is being eroded, projectile should also fragment
            if (targ_m+imp_m <= 2*min_frag_mass){ //not enough mass to produce new fragments
                                                  //printf("ELASTIC BOUNCE\n");
                params->collision_type=0;
                reb_collision_resolve_hardsphere(r,c);
                swap = 0;
            } else{
                params->Mlr = MAX(params->Mlr, min_frag_mass);
                //printf("GRAZING PARTIAL EROSION\n");
                params->collision_type = 3;
                add_fragments(r,c,params);
            }
        } else{ //Mlr > Mt, either a hit and run or an elastic bounce
            double Mlr_dag = (beta*target->m + projectile->m)/10 * pow(Q/(1.8*Q_star), -1.5);
            if (Q < 1.8*Q_star){
                Mlr_dag = (beta*targ_m + imp_m)*(1 - Q/ (2*Q_star));
            }

            double projectile_mass_accreted = params->Mlr - targ_m;
            double new_projectile_mass = projectile->m - projectile_mass_accreted;
            Mlr_dag = MAX(Mlr_dag, min_frag_mass);
            if (new_projectile_mass-Mlr_dag < min_frag_mass){
                //printf("ELASTIC BOUNCE\n");
                params->collision_type=0;
                reb_collision_resolve_hardsphere(r,c);
                swap = 0;
            }
            else{
                params->Mslr = Mlr_dag;
                //printf("HIT AND RUN\n");
                params->collision_type = 2;
                add_fragments(r,c,params);
            }
        }
        return swap;
    }
}

//*** PRODUCING AN EJECTION FILE ***

//The following snippet of code shows one way to keep track of ejected particles during a REBOUND simulation using the C heartbeat function. 
//Every 100 years of simulation time, the heartbeat checks to see if any objects have a semi-major axis that exceeds 100 AU. 
//If a particle meeets this condition, then its removed from the simulation and its hash and time of ejection are recorded.

void heartbeat(struct reb_simulation* sim){
    trace_close_encounters += sim->ri_trace.encounter_N;
    mercurius_close_encounters += sim->ri_mercurius.encounter_N;
    if (reb_simulation_output_check(sim, 100)){   //heartbeat is at every 100 years
        for (int i=0;i<sim->N;i++){
            struct reb_orbit o =  reb_orbit_from_particle(sim->G, sim->particles[i], sim->particles[0]);
            if (o.a > 100.){
                int keepSorted = 1;
                int removed_hash = sim->particles[i].hash;
                reb_simulation_remove_particle_by_hash(sim, removed_hash, keepSorted);
            }
        }
    }
    reb_simulation_move_to_com(sim);
}

struct collision_params* create_collision_params(){
    struct collision_params* params = calloc(1,sizeof(struct reb_simulation));
    return params;
}


int reb_collision_resolve_fragment(struct reb_simulation* const r, struct reb_collision c){
    if (r->particles[c.p1].last_collision==r->t || r->particles[c.p2].last_collision==r->t) return 0;
    int i = c.p1;
    int j = c.p2; 
    if (i<j) return 0;      //only return one collision callback

    int swap = 2;
    if (r->particles[i].m < r->particles[j].m){        //unless swap is redfined as 0, projectile is going to be removed.
        swap =1;
        i = c.p2;
        j = c.p1;
    }

    struct reb_particle* particles = r->particles;
    struct collision_params* params = create_collision_params();

    double imp_r = particles[j].r;
    double targ_r = particles[i].r;
    double R_tot = imp_r + targ_r;

    double imp_m = particles[j].m;
    double targ_m = particles[i].m;

    double M_tot = imp_m + targ_m;
    double G = r->G;
    double Mlr,dx,dy,dz,Vix,Viy,Viz;
    double x2rel, xrel, v2rel, v2imp, Vi;
    double hx,hy,hz,h2,b;
    make_vector(particles[i].x, particles[i].y, particles[i].z, particles[j].x, particles[j].y, particles[j].z, &dx,&dy,&dz);  //find relative coordinates dx, dy,dz
    x2rel = get_dot(dx,dy,dz,dx,dy,dz); 
    make_vector(particles[i].vx, particles[i].vy, particles[i].vz, particles[j].vx, particles[j].vy, particles[j].vz, &Vix,&Viy,&Viz);  //find relative velocity
    v2rel = get_dot(Vix,Viy,Viz,Vix,Viy,Viz);

    xrel = sqrt(x2rel);  //distance between the centers of the projectile and target


    hx = (dy*Viz - dz*Viy);                     //angular momentum vector xrel X Vrel
    hy = (dz*Vix - dx*Viz);
    hz = (dx*Viy - dy*Vix);

    h2 = get_dot(hx,hy,hz,hx,hy,hz);

    v2imp = v2rel + 2*G*M_tot*(1./R_tot - 1./xrel); //impact velocity with gravitational focusing at time of detected collision

    if (1./R_tot - 1./xrel < 0){v2imp = v2rel;}  //if collision is detected after physical contact

    Vi = sqrt(v2imp);  //magnitude of impact velocity vector
    b = sqrt(h2/v2imp);  //impact parameter, b=R_tot*sin(theta)
    if (b != b){
        //printf("NAN b \n");
        exit(0);}
    //Stewart & Leinhardt 2012 parameters
    double mu = (targ_m*imp_m)/M_tot;  //Chambers Eq. 2, reduced mass
    double l = R_tot-b;  //Leinhardt Eq. 7, the projected length of the projectile overlapping the target
    l = MIN(l, 2*imp_r);
    double alpha = (pow(l,2)*(3*imp_r-l))/(4*pow(imp_r, 3)); //Leinhardt Eq. 11, interacting mass fraction
    alpha = MIN(1., alpha);
    double Q = .5*v2imp*targ_m*imp_m/pow(M_tot,2);  //specific energy per unit mass
    double V_esc = pow(2.*G*M_tot/R_tot, .5); //mutal escape velocity as defined in Wallace et al 2018 and Chambers 2013
    double alphamu = (alpha*targ_m*imp_m)/(alpha*imp_m + targ_m);  //Leinhardt Eq. 12, reduced interacting mass for fraction alpha.
    double gamma = imp_m/targ_m;  //Chambers Eq. 6

    const double cstar = 1.8;      //may be a user defined variable, default taken from paper

    double rho1;         //constant density

    if (G==6.674e-8){rho1 =1;} //CGS
    if (G==6.674e-11){rho1 =1000;} //SI
    if (G==39.476926421373 || G==1){rho1 = 1.684e6;}  //Msun/AU^3
    double Rc1 = pow((M_tot*3)/(4.*M_PI*rho1), 1./3.);  //Chambers Eq. 4, combined radius of target and projectile with constant density
    double Q0 = .8*cstar*M_PI*rho1*G*pow(Rc1,2);  //Chambers Eq. 3, critical value of impact energy for head-on collisions
    double Q_star = pow(mu/alphamu, 1.5)*(pow(1+gamma, 2)/ (4*gamma))*Q0;  //Chambers Eq. 5, critical value for oblique or different mass collisons.  
    if (alpha == 0.0){Q_star = 6364136223846793005.0;}
    if (b == 0 && imp_m == targ_m){
        Q_star = Q0;
    }
    double qratio = Q/Q_star;
    if (qratio < 1.8){
        Mlr = M_tot*(1.0-.5*qratio);
    }
    else{
        Mlr = .1*M_tot*pow(qratio/1.8, -1.5);  //Chambers Eq.8
    }

    double separation_distance = 4 * R_tot;  //seperation distance of fragments.  Should be userdefined but this is what chambers uses
                                             ///POPULATE STRUCT OBJECTS
    params->target = i;
    params->projectile =j;
    params->dx = dx;
    params->dy = dy;
    params->dz = dz;
    params->b = b;
    params->Vix = Vix;
    params->Viy = Viy;
    params->Viz = Viz;
    params->Vi = Vi;
    params->l = l;
    params->rho1 = rho1;
    params->cstar = cstar;
    params->mu = mu;
    params->Q = Q;
    params->separation_distance = separation_distance;
    params->V_esc = V_esc;
    params->vrel = sqrt(v2rel);
    params->Mslr = 0;
    params->xrel = xrel;
    params->Mlr = Mlr;


    if (Vi <= V_esc){
        params->collision_type = 1;
        //printf("SIMPLY MERGED\n");
        merge(r,c, params);
    }
    else{  //Vi > V_esc
        if (b<targ_r){ //non-grazing regime
            if (M_tot - params->Mlr < min_frag_mass){
                params->collision_type = 1;
                //printf("EFFECTIVELY MERGED\n");
                merge(r,c,params);
            }
            else{ // M_tot - params->Mlr >= min_frag_mass; fragments will be produced unless it is a graze and merge or elastic bounce 
                if (params->Mlr < targ_m){
                    if (params->Mlr <= 0.1*targ_m){
                        //printf("SUPER-CATASTROPHIC\n");
                        params->collision_type = 4;
                        params->Mlr = MAX(Mlr, min_frag_mass);
                        add_fragments(r,c,params);
                    }
                    else{
                        //printf("PARTIAL EROSION\n");
                        params->collision_type = 3;
                        params->Mlr = MAX(Mlr, min_frag_mass);
                        add_fragments(r,c,params);
                    }
                }

                else{  //(params->Mlr >= targ_m)
                       //printf("PARTIAL ACCRETION\n");
                    params->collision_type = 2;
                    add_fragments(r,c,params);
                }
            }
        }
        else{ // b > b_crit, grazing regime
            swap = hit_and_run(r,c,params); //swap gets redefined here as it may be set to 0 in the case of a bounce
        }
    }

    return swap;
}


void heartbeat(struct reb_simulation* r);

int add_particles_from_file(struct reb_simulation* r, char* filename, double m){
    double rho = 5.05e6; //3 g/cm^3
    size_t len = 0;
    char *line = NULL;
    FILE* file = fopen(filename, "r");
    int added = 0;
    while (getline(&line, &len, file) != -1) {
        double a;
        sscanf(line, "%lf", &a);
        double inc = reb_random_uniform(r,0,0.0175);
        double ecc = reb_random_uniform(r,0,0.01);
        double omega = reb_random_uniform(r,0,2*M_PI);
        double Omega = reb_random_uniform(r,0,2*M_PI);
        double f = reb_random_uniform(r,0,2*M_PI);
        struct reb_particle p = reb_particle_from_orbit(r->G, r->particles[0], m, a, ecc, inc, Omega, omega, f);
        p.r = get_radii(m, rho);
        p.hash = r->N;
        reb_simulation_add(r, p); 
        added += 1;
    }
    fclose(file);
    return added;
}


int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    reb_simulation_start_server(r, 1234);
    r->G = 39.476926421373;
    r->dt = 6./365.;

    if (argc < 3){
        printf("Usage: ./rebound INTEGRATOR SEED\n");
        return -1;
    }
    r->rand_seed = atoi(argv[2]);
    if (strcmp(argv[1],"mercurius")==0){
        r->integrator = REB_INTEGRATOR_MERCURIUS;
    }
    if (strcmp(argv[1],"trace")==0){
        r->integrator = REB_INTEGRATOR_TRACE;
        r->ri_trace.r_crit_hill = 3.*1.21;
        r->ri_trace.S_peri = reb_integrator_trace_switch_peri_none;
    }

    char filename[1024];
    sprintf(filename, "%s_%2d.bin",argv[2],r->rand_seed);

    r->heartbeat = heartbeat;
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_fragment;



    struct reb_particle star = {0};
    star.m = 1.00;
    star.r = 0.1; 
    reb_simulation_add(r, star);

    // You can substitute your own input file of semi major axis here. 
    // Here, I am using the disk from Chambers et al. 2013
    add_particles_from_file(r, "semis_emb_chambers.txt", 2.8e-7);
    add_particles_from_file(r, "semis_pl_chambers.txt", 2.8e-8);


    //Add Jupiter and Saturn
    double rho = 5.05e6; //3 g/cm^3
    struct reb_particle Jup = reb_particle_from_orbit(r->G, r->particles[0], 9.543e-4, 5.20349, 0.048381, 0.365*(M_PI/180), 0.0, 68.3155*(M_PI/180), 227.0537*(M_PI/180));
    Jup.r = get_radii(Jup.m, rho); 
    Jup.hash = reb_hash("JUPITER");
    reb_simulation_add(r,Jup);

    struct reb_particle Sat = reb_particle_from_orbit(r->G, r->particles[0], 0.0002857, 9.54309, 0.052519, 0.8892*(M_PI/180), M_PI, 324.5263*(M_PI/180), 256.9188*(M_PI/180));
    Sat.r = get_radii(Sat.m, rho);
    Sat.hash = reb_hash("SATURN");
    reb_simulation_add(r,Sat);


    reb_simulation_move_to_com(r);  // This makes sure the planetary systems stays within the computational domain and doesn't drift.
    double run_time = 2e7;
    reb_simulation_save_to_file_interval(r,filename,1.e5);
    reb_simulation_integrate(r, run_time);

    //printf("Total Steps: %d Rejected Steps: %ld MERCRUIUS Encounters: %d TRACE encounters: %d\n", r->steps_done, r->ri_trace.step_rejections, mercurius_close_encounters, trace_close_encounters);
}
