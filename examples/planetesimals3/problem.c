/**
 * A.S. This is my planetesimal disk with special integrator for close encounters.
 *      Particle id's: 0 = star, 1 = massive body, 2 = planetesimal, 3 = CLOSE ENCOUNTER
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);
double E0, t_output, t_log_output;
char* output_name;

//swifter/mercury compare
void output_to_mercury_swifter(struct reb_simulation* r, double HSR, double tmax, int n_output);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();

    double planetesimal_mass = 3e-8;
    double amin = 0.45, amax = 0.75;        //for planetesimal disk
    double powerlaw = 0.5;
    double tmax = 1000;
    int N_planetesimals = 100;
    int seed = 30;
    output_name = "output/Energy.txt";
    //double tmax = atof(argv[1]);
    //int N_planetesimals = atoi(argv[2]);
    //int seed = atoi(argv[3]);
    //output_name = argv[4];
    
	//Simulation Setup
	r->integrator	= REB_INTEGRATOR_HYBARID;
	//r->integrator	= REB_INTEGRATOR_IAS15;
	//r->integrator	= REB_INTEGRATOR_WHFAST;
    r->ri_hybarid.switch_ratio = sqrt(3.);  //Hill radii
    r->ri_hybarid.CE_radius = 18.;          //X*radius
    r->ri_hybarid.ejection_distance = 10.;  //AU
    r->ri_hybarid.collisions = 1;
    r->testparticle_type = 1;
	r->heartbeat	= heartbeat;
    r->dt = 0.0015;
    
	// Initial conditions
	struct reb_particle star = {0};
	star.m 		= 1;
    star.r		= 0.005;        // Radius of particle is in AU!
	reb_add(r, star);
   
    srand(seed);
    int n_output = 10000;
    t_log_output = pow(tmax + 1, 1./(n_output - 1));
    t_output = r->dt;
    
    //planet 1
    {
        double a=0.5, m=5e-4, e=0, inc = reb_random_normal(0.00001);
        struct reb_particle p1 = {0};
        p1 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, 0);
        p1.r = 1.6e-4;              //radius of particle is in AU!
        p1.id = r->N;
        reb_add(r, p1);
    }
    
    //planet 2
    {
        double a=0.7, m=5e-4, e=0.01, inc=reb_random_normal(0.00001);
        struct reb_particle p2 = {0};
        p2 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, 0);
        p2.r = 1.6e-4;
        p2.id = r->N;
        reb_add(r, p2);
    }
    
    r->N_active = r->N;
    reb_move_to_com(r);
    
    //planetesimals
    while(r->N<N_planetesimals + r->N_active){
		struct reb_particle pt = {0};
		double a	= reb_random_powerlaw(amin,amax,powerlaw);
        double phi 	= reb_random_uniform(0,2.*M_PI);
        double inc = reb_random_normal(0.0001);
        double Omega = reb_random_uniform(0,2.*M_PI);
        double apsis = reb_random_uniform(0,2.*M_PI);
        pt = reb_tools_orbit_to_particle(r->G, star, planetesimal_mass, a, 0., inc, Omega, apsis,phi);
		pt.r 		= 4e-5;
        pt.id = r->N;
		reb_add(r, pt);
    }
   
    //energy
    E0 = reb_tools_energy(r);
    system("rm -v output/Energy.txt");
    
    //swifter/mercury compare
    output_to_mercury_swifter(r, r->ri_hybarid.switch_ratio, tmax, n_output);
    
    time_t t_ini = time(NULL);
    struct tm *tmp = gmtime(&t_ini);
    
    //Integrate!
    reb_integrate(r, tmax);
    
    time_t t_fini = time(NULL);
    struct tm *tmp2 = gmtime(&t_fini);
    double time = t_fini - t_ini;
    printf("\nSimulation complete. Elapsed simulation time is %.2f s, \n\n",time);
    
}

void heartbeat(struct reb_simulation* r){
    if(r->t > t_output){//log output
        t_output = r->t*t_log_output;
        
        double E = reb_tools_energy(r) + r->ri_hybarid.dE_offset;
        double dE = fabs((E-E0)/E0);
        reb_output_timing(r, 0);
        printf("    dE=%e",dE);
        
        FILE *append;
        append = fopen(output_name, "a");
        fprintf(append, "%.16f,%.16f\n",r->t,dE);
        fclose(append);
    }
}

void output_to_mercury_swifter(struct reb_simulation* r, double HSR, double tmax, int n_output){
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
        double m = p.m*mass_conv;
        double rr = sqrt((p.x-p0.x)*(p.x-p0.x) + (p.y-p0.y)*(p.y-p0.y) + (p.z-p0.z)*(p.z-p0.z));
        fprintf(swifter," %d %.16f %f\n",i+1,m,rr*pow(p.m/(particles[0].m*3.),2./3.));
        fprintf(swifter," %f\n",p.r);
        fprintf(swifter," %.16f %.16f %.16f\n",p.x - p0.x, p.y - p0.y, p.z - p0.z);
        fprintf(swifter," %.16f %.16f %.16f\n",(p.vx - p0.vx)*vel_conv,(p.vy - p0.vy)*vel_conv,(p.vz - p0.vz)*vel_conv);
    }
    
    //SWIFTER - Other params (time, dt, etc.)
    int output_rate = tmax/n_output;
    fprintf(swifterparams,"! \n");
    fprintf(swifterparams,"! Parameter file for Swifter, with N=%d total bodies. \n",r->N);
    fprintf(swifterparams,"! \n! \n");
    fprintf(swifterparams,"T0             0.0E0 \n");
    fprintf(swifterparams,"TSTOP          %e        !In units where G=1\n",tmax);
    fprintf(swifterparams,"DT             %e        !In units where G=1\n",r->dt);
    fprintf(swifterparams,"PL_IN          swifter_pl.in\n");
    fprintf(swifterparams,"!TP_IN         tp.in     !Commented out for now, no test par\n");
    fprintf(swifterparams,"IN_TYPE        ASCII\n");
    fprintf(swifterparams,"ISTEP_OUT      %d        !# timesteps between outputs \n",output_rate);
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
        double rr = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
        fprintf(mercuryb," BODY%d      m=%.16f r=%f\n",i,p.m,HSR*rr*pow(p.m/(particles[0].m*3.),2./3.));
        fprintf(mercuryb," %.16f %.16f %.16f\n",p.x - p0.x,p.y - p0.y,p.z - p0.z);
        fprintf(mercuryb," %.16f %.16f %.16f\n",(p.vx - p0.vx)*AU_d,(p.vy - p0.vy)*AU_d,(p.vz - p0.vz)*AU_d);   //AU/day
        fprintf(mercuryb," 0. 0. 0.\n");
    }
    //mini bodies
    for(int i=N_active;i<N;i++){
        struct reb_particle p = particles[i];
        double rr = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
        fprintf(mercurys," BODY%d      m=%.16f r=%f\n",i,p.m,HSR*rr*pow(p.m/(particles[0].m*3.),2./3.));
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
    fprintf(mercuryparams," output interval (days) = %.2fd0\n",(tmax/n_output)*365);
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
