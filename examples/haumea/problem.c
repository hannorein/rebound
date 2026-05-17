#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#define N_p 10000 // particle number

double coefficient_of_restitution_bridges(const struct reb_simulation* const r, double v);
void heartbeat(struct reb_simulation* const r);
void force_J2_C22(struct reb_simulation* r);
void force_J2_inertial(struct reb_simulation* r);

const double Mplanet = 3.95e21;// kg
const double Rplanet = 7.8e5; //m mean radius of Haumea
double OMEGA = 0.0001485;    // 1/s T = 3*3.92 hours (assuming 3:1 ratio between T_ring and T_haumea)
const double OMEGA_H = 0.0004455;
const double J2planet = 0.24; // dimensionless
const double C22planet = 0.05; // dimensionless
const int t_c = 400; // time during which impact velocities are measured (measure in periods)
const int tmax_wo_c = 1100; // time before measuring collisions

int main(int argc, char* argv[]) {
	struct reb_simulation* r = reb_simulation_create();
	
	// This allows you to connect to the simulation using
	// a web browser. Simply go to http://localhost:1234
	reb_simulation_start_server(r, 1234);
	// Setup constants
	r->integrator        = REB_INTEGRATOR_WHFAST;
    r->N_active = 1;
	r->G                 = 6.67428e-11;         // N / kg^2 m^2
    double a_min = 2.235e6;
    double a_max = 2.34e6;
    double a_avg = (a_min + a_max)/2;
    OMEGA = sqrt(r->G*Mplanet/pow(a_avg, 3));
	r->dt                = 1e-3*2*M_PI/OMEGA;  // s
	r->heartbeat         = heartbeat;           // function pointer for heartbeat
	double surfacedensity          = 640;     // kg/m^2
	double particle_density        = 400;     // kg/m^3
	double particle_radius_min     = 1;       // m
	double particle_radius_max     = 2;       // m
	double particle_radius_slope   = -3;    

	printf("Toomre wavelength: %f\n",(4.*M_PI*M_PI*surfacedensity*r->G)/(OMEGA*OMEGA));
	
    struct reb_particle planet = {0};
    planet.m  = Mplanet;
    //planet.r = Rplanet;
    reb_simulation_add(r, planet);
    float a, e, v, corr;
    float _omega_0[N_p];
    struct reb_orbit o_0;
    for (int i = 0; i < N_p; i++) {
        a = reb_random_uniform(r, a_min, a_max);
        corr = sqrt(1 + ((3.0/2.0) * J2planet + 3 * C22planet) * (Rplanet * Rplanet/a/a));
        struct reb_particle p = {0};          // test particle
        e = 0.0;//reb_random_uniform(r, 1e-5, 1e-4);
        float inc = reb_random_uniform(r, 0.0, 9e-3);
        float radius = reb_random_uniform(r, 0.1, 1.0);
        //float m = 4/3*M_PI*pow(radius, 3)*particle_density;
        reb_simulation_add_fmt(r, "a e inc primary", a, e, inc, r->particles[0]);
        r->particles[i+1].vx *= corr;
        r->particles[i+1].vy *= corr;
    }

    reb_simulation_move_to_com(r);
    r->additional_forces = force_J2_C22;
    float period = 2*M_PI/OMEGA;

    reb_simulation_integrate(r, tmax_wo_c*period);
    double _e[N_p], _a[N_p], _Omega[N_p], _omega[N_p], _inc[N_p], _x[N_p], _y[N_p];

    for (int i = 1; i < N_p + 1; i++) {
        
        struct reb_orbit o= reb_orbit_from_particle(r->G, r->particles[i], r->particles[0]);
        corr = sqrt(1 + ((3.0/2.0) * J2planet + 3 * C22planet) * (Rplanet * Rplanet/a/a));
        //r->particles[i].vx /= corr;
        //r->particles[i].vy /= corr;
        o= reb_orbit_from_particle(r->G, r->particles[i], r->particles[0]);
        _e[i-1] = o.e;
        _a[i-1] = o.a / 1e6;
        _Omega[i-1] = o.Omega;
        _omega[i-1] = o.omega;
        _inc[i-1] = o.inc;
        _x[i-1] = o.d*cos(o.omega + o.f)/1e3;
        _y[i-1] = o.d*sin(o.omega + o.f)/1e3;
    }

    FILE *fp;
    FILE *data;
    fp = fopen("commands.gplot", "w");

    //fprintf(fp, "set terminal qt persist\n");
    //fprintf(fp, "set lmargin 10\n");
    //fprintf(fp, "set rmargin 2\n");
    //fprintf(fp, "set tmargin 2\n");
    //fprintf(fp, "set bmargin 4\n");
    //fprintf(fp, "set multiplot layout 2,2\n");
    //fprintf(fp, "set xlabel font ', 10'\n");
    //fprintf(fp, "set ylabel font ', 10'\n");
    //fprintf(fp, "set xtics font ',12'\n");
    //fprintf(fp, "set ytics font ',12'\n");
    //fprintf(fp, "set key font ',12'\n");
    //fprintf(fp, "set xlabel 'a (10^3 km)'\n");
    //fprintf(fp, "set ylabel 'Omega'\n");
    //fprintf(fp, "plot [ %f : %f ] \'fn1.dat\' w points, 0 \n", 2.0, 2.6);     
    
    //fprintf(fp, "set xlabel 'a (10^3 km)'\n");
    //fprintf(fp, "set ylabel 'i (rad)'\n");
    //fprintf(fp, "plot [ %f : %f ] \'fn2.dat\' w points, 0 \n", 2.0, 2.6);

    fprintf(fp, "set terminal qt persist size 1400,1000\n");

    fprintf(fp, "set multiplot layout 2,1 rowsfirst spacing 0.08,0.08 margins 0.10,0.95,0.10,0.95\n");

    fprintf(fp, "set tics font ',8'\n");
    fprintf(fp, "set xlabel font ',12'\n");
    fprintf(fp, "set ylabel font ',12'\n");
    fprintf(fp, "set title font ',11'\n");
    fprintf(fp, "set xtics font ',12'\n");
    fprintf(fp, "set ytics font ',12'\n");
    fprintf(fp, "set key off\n");
    fprintf(fp, "set grid\n");
    /*
    fprintf(fp, "set xlabel 'a (10^3 km)'\n");
    fprintf(fp, "set ylabel 'e '\n");
    fprintf(fp, "plot [ %f : %f ] \'fn3.dat\' w points, 0 \n", 2.0, 2.6);
    
    fprintf(fp, "set xlabel 'a (10^3 km)'\n");
    fprintf(fp, "set ylabel 'omega'\n");
    fprintf(fp, "plot [ %f : %f ] \'fn4.dat\' w points, 0 \n", 2.0, 2.6);
    */

    //fprintf(fp, "set xlabel 'x (km)'\n");
    //fprintf(fp, "set ylabel 'y (km)'\n");
    //fprintf(fp, "set yrange [-3000.:3000.]\n");
    //fprintf(fp, "plot [ %f : %f ] \'fn5.dat\' w points, 0 \n", -3000., 3000.);

    fprintf(fp, "set xlabel 't (period)'\n");
    fprintf(fp, "set ylabel 'Impact velocity (m/s)'\n");
    fprintf(fp, "set yrange [*:*]\n");
    fprintf(fp, "plot [ %f : %f ] \'fn6.dat\' w points, 0 \n", (float)tmax_wo_c, (float)(tmax_wo_c + t_c));
    
    fprintf(fp, "set xlabel 't (period)'\n");
    fprintf(fp, "set ylabel 'Number of collisions'\n");
    fprintf(fp, "plot [ %f : %f ] \'fn7.dat\' w points, 0 \n", (float)tmax_wo_c, (float)(tmax_wo_c + t_c));

    fprintf(fp, "unset multiplot\n");    
    fclose(fp);
    /*
    data = fopen("fn1.dat", "w");
    for (int i = 0; i < N_p; i++) {
        fprintf(data, "%25.15f  %25.15f\n",(float)_a[i], (float)_Omega[i]);
    }
    fclose(data);
    
    data = fopen("fn2.dat", "w");
    for (int i = 0; i < N_p; i++) {
        fprintf(data, "%25.15f  %25.15f\n",(float)_a[i], (float)_inc[i]);
    }
    fclose(data);
    */
    /*
    data = fopen("fn3.dat", "w");
    for (int i = 0; i < N_p; i++) {
        fprintf(data, "%25.15f  %25.15f\n",(float)_a[i], (float)_e[i]);
    }
    fclose(data);
    
    data = fopen("fn4.dat", "w");
    for (int i = 0; i < N_p; i++) {
        fprintf(data, "%25.15f  %25.15f\n",(float)_a[i], (float)_omega[i]);
    }
    fclose(data);
    */
    /*
    data = fopen("fn5.dat", "w");
    
    for (int i = 0; i < N_p; i++) {
        fprintf(data, "%25.15f  %25.15f\n",(float)_x[i], (float)_y[i]);
    }
    fclose(data);
    */
    
    double rx, ry, rz, r2, vx, vy, vz, v2;
    double delta_r = r->dt*340./2.2; // two times inflated radius
    int _t_c = 5, n;
    for (int i = 0; i < t_c/_t_c; i++) {
        if (i == 0) {
            data = fopen("fn6.dat", "w");
        } else {
            data = fopen("fn6.dat", "a");
        }
        reb_simulation_integrate(r, (tmax_wo_c + _t_c*(i+1))*period);
        n = 0;
        for (int j = 0; j < N_p; j++) {
            for (int k = j + 1; k < N_p; k++) {
                rx = r->particles[j + 1].x - r->particles[k + 1].x;
                ry = r->particles[j + 1].y - r->particles[k + 1].y;
                rz = r->particles[j + 1].z - r->particles[k + 1].z;
                r2 = rx*rx + ry*ry + rz*rz;
                if (sqrt(r2) < delta_r) {
                    vx = r->particles[j + 1].vx - r->particles[k + 1].vx;
                    vy = r->particles[j + 1].vy - r->particles[k + 1].vy;
                    vz = r->particles[j + 1].vz - r->particles[k + 1].vz;
                    v2 = vx*vx + vy*vy + vz*vz;
                    printf("v: %f\n", sqrt(v2));
                    printf("r: %f\n", sqrt(r2));
                    fprintf(data, "%25.15f  %25.15f\n",(float)(r->t/period), (float)(sqrt(v2)));
                    n += 1;
                }
            }
        }
        fclose(data);
        if (i == 0) {
            data = fopen("fn7.dat", "w");
        } else {
            data = fopen("fn7.dat", "a");
        }
        fprintf(data, "%25.15f  %25.15f\n",(float)(r->t/period), (float)(n));
        fclose(data);
    }
    fclose(data);
    system("gnuplot commands.gplot");
    reb_simulation_free(r);
}

/// @brief Zonal and tesseral contributions are computed in the body-fixed, rotating frame and then converted back
/// into the inertial frame.
/// @param r 
void force_J2_C22(struct reb_simulation* r){
    if (J2planet==0 || r->particles[0].m == 0) return;

    const struct reb_particle planet = r->particles[0];     // cache
    const int N = r->N;
#pragma omp parallel for
    for (int i=1;i<N;i++){
        const struct reb_particle p = r->particles[i];      // cache
        // setup rotation around the z axis
        struct reb_vec3d axis = {.x = 0, .y = 0, .z = 1};
        struct reb_rotation R1 = reb_rotation_init_angle_axis(-OMEGA_H*r->t, axis);

        struct reb_vec3d _r = {p.x-planet.x, p.y-planet.y, p.z-planet.z}; // position in the inertial frame
        struct reb_vec3d pr = reb_vec3d_rotate(_r, R1); // position in the rotating frame
        
        const double pr2  = pr.x*pr.x + pr.y*pr.y + pr.z*pr.z;        // distance^2 relative to planet
        const double fac  = -3.*r->G*planet.m*Rplanet*Rplanet/pow(pr2,3.5);
        
        // J2 term
        double pax  = J2planet*fac*pr.x*(pr.x*pr.x + pr.y*pr.y - 4.*pr.z*pr.z)/2.;
        double pay  = J2planet*fac*pr.y*(pr.x*pr.x + pr.y*pr.y - 4.*pr.z*pr.z)/2.;
        double paz  = J2planet*fac*pr.z*(3.*(pr.x*pr.x + pr.y*pr.y) - 2.*pr.z*pr.z)/2.;
        // C22 term
        pax += C22planet*fac*pr.x*(2.*pr.z*pr.z + 7.*pr.y*pr.y - 3.*pr.x*pr.x);
        pay += C22planet*fac*pr.y*(3.*pr.y*pr.y - 2.*pr.z*pr.z - 7.*pr.x*pr.x);
        paz += C22planet*fac*pr.z*(5.*pr.y*pr.y - 5.*pr.x*pr.x);
        
        struct reb_rotation R2 = reb_rotation_init_angle_axis(OMEGA_H*r->t, axis);
        struct reb_vec3d pa = {pax,pay,paz}; // acceleration in the rotating frame
        struct reb_vec3d a = reb_vec3d_rotate(pa, R2); // acceleration in the inertial frame

        r->particles[i].ax += a.x;
        r->particles[i].ay += a.y;
        r->particles[i].az += a.z;

        const double mfac = r->particles[i].m/r->particles[0].m;

        r->particles[0].ax -= mfac*a.x;
        r->particles[0].ay -= mfac*a.y;
        r->particles[0].az -= mfac*a.z;
    }
}

void force_J2_inertial(struct reb_simulation* r){
    if (J2planet==0 || r->particles[0].m == 0) return;

    const struct reb_particle planet = r->particles[0];     // cache
    const int N = r->N;
#pragma omp parallel for
    for (int i=1;i<N;i++){
        const struct reb_particle p = r->particles[i];      // cache
        // Coordinates in inertial frame
        const double prx = p.x-planet.x;
        const double pry = p.y-planet.y;
        const double prz = p.z-planet.z;

        const double pr2  = prx*prx + pry*pry + prz*prz;        // distance^2 relative to planet
        const double fac  = -3.*r->G*planet.m*Rplanet*Rplanet/pow(pr2,3.5);
        
        // J2 term
        double pax  = J2planet*fac*prx*(prx*prx + pry*pry - 4.*prz*prz)/2.;
        double pay  = J2planet*fac*pry*(prx*prx + pry*pry - 4.*prz*prz)/2.;
        double paz  = J2planet*fac*prz*(3.*(prx*prx + pry*pry) - 2.*prz*prz)/2.;

        r->particles[i].ax += pax;
        r->particles[i].ay += pay;
        r->particles[i].az += paz;

        const double mfac = r->particles[i].m/r->particles[0].m;

        r->particles[0].ax -= mfac*pax;
        r->particles[0].ay -= mfac*pay;
        r->particles[0].az -= mfac*paz;
    }
}

// This example is using a custom velocity dependend coefficient of restitution
double coefficient_of_restitution_bridges(const struct reb_simulation* const r, double v){
	// assumes v in units of [m/s]
	double eps = 0.32*pow(fabs(v)*100.,-0.234);
	if (eps>1) eps=1;
	if (eps<0) eps=0;
	return eps;
}

void heartbeat(struct reb_simulation* const r){
	if (reb_simulation_output_check(r, 1e-1*2.*M_PI/OMEGA)){
		reb_simulation_output_timing(r, 0);
		//reb_output_append_velocity_dispersion("veldisp.txt");
	}
}