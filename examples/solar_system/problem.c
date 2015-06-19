/**
 * @file 	problem.c
 * @brief 	Example problem: solar system.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This example integrates all planets of the Solar
 * System. The data comes from the NASA HORIZONS system. 
 * 
 * @section 	LICENSE
 * Copyright (c) 2015 Hanno Rein, Shangfei Liu, Dave Spiegel
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
#include "main.h"
#include "output.h"
#include "tools.h"
#include "particle.h"
#include "boundaries.h"
#include "integrator.h"
#include "integrator_whfast.h"

double ss_pos[10][3] = 
{
	{3.256101656448802E-03  , -1.951205394420489E-04 , -1.478264728548705E-04},  
	{-1.927589645545195E-01 , 2.588788361485397E-01  , 3.900432597062033E-02 }, 
	{-5.976537074581466E-01 , 3.918678996109574E-01  , 3.990356741282203E-02 }, 
	{-7.986189029000561E-01 , -6.086873314992410E-01 , -1.250824315650566E-04}, 
	{7.897942807177173E-01  , 1.266671734964037E+00  , 7.092292179885432E-03 }, 
	{-4.314503046344270E+00 , 3.168094294126697E+00  , 8.331048545353310E-02 }, 
	{-4.882304833383455E+00 , -8.689263067189865E+00 , 3.453930436208210E-01 }, 
	{1.917757033372740E+01  , 5.671738750949031E+00  , -2.273858614425555E-01},  
	{2.767031517959636E+01  , -1.150331645280942E+01 , -4.008018419157927E-01},  
	{7.765250227278298E+00  , -3.190996242617413E+01 , 1.168394015703735E+00 }, 

};
double ss_vel[10][3] = 
{
	{3.039963463108432E-06 ,  6.030576499910942E-06 ,  -7.992931269075703E-08}, 
	{-2.811550184725887E-02,  -1.586532995282261E-02,  1.282829413699522E-03 }, 
	{-1.113090630745269E-02,  -1.703310700277280E-02,  4.089082927733997E-04 },
	{1.012305635253317E-02 ,  -1.376389620972473E-02,  3.482505080431706E-07 }, 
	{-1.135279609707971E-02,  8.579013475676980E-03 ,  4.582774369441005E-04 }, 
	{-4.555986691913995E-03,  -5.727124269621595E-03,  1.257262404884127E-04 }, 
	{4.559352462922572E-03 ,  -2.748632232963112E-03,  -1.337915989241807E-04}, 
	{-1.144087185031310E-03,  3.588282323722787E-03 ,  2.829006644043203E-05 }, 
	{1.183702780101068E-03 ,  2.917115980784960E-03 ,  -8.714411604869349E-05}, 
	{3.112825364672655E-03 ,  1.004673400082409E-04 ,  -9.111652976208292E-04},
};

double ss_mass[10] =
{
	1.988544e30,
	3.302e23,
	48.685e23,
	6.0477246e24,
	6.4185e23,
	1898.13e24,
	5.68319e26,
	86.8103e24,
	102.41e24,
	1.4639248e+22,
};

double energy();
double e_init;

void problem_init(int argc, char* argv[]){
	// Setup constants
	dt 		= 4;				// in days
	tmax		= 7.3e10;			// 200 Myr
	G		= 1.4880826e-34;		// in AU^3 / kg / day^2.
	init_boxwidth(200); 				// Init box with width 200 astronomical units
	integrator_whfast_safe_mode = 0;		// Turn off safe mode. Need to call integrator_synchronize() before outputs. 
	integrator_whfast_corrector = 11;		// 11th order symplectic corrector
	integrator_force_is_velocitydependent = 0;	// Force only depends on positions. 
	//integrator	= WH;
	integrator	= WHFAST;
	//integrator	= IAS15;

	// Initial conditions
	for (int i=0;i<10;i++){
		struct particle p;
		p.x  = ss_pos[i][0]; 		p.y  = ss_pos[i][1];	 	p.z  = ss_pos[i][2];
		p.vx = ss_vel[i][0]; 		p.vy = ss_vel[i][1];	 	p.vz = ss_vel[i][2];
		p.ax = 0; 			p.ay = 0; 			p.az = 0;
		p.m  = ss_mass[i];
		particles_add(p); 
	}
	if (integrator==WH){
		// Move to heliocentric frame (required by WH integrator)
		for (int i=1;i<N;i++){
			particles[i].x -= particles[0].x;	particles[i].y -= particles[0].y;	particles[i].z -= particles[0].z;
			particles[i].vx -= particles[0].vx;	particles[i].vy -= particles[0].vy;	particles[i].vz -= particles[0].vz;
		}
		particles[0].x = 0;	particles[0].y = 0;	particles[0].z = 0;
		particles[0].vx= 0;	particles[0].vy= 0;	particles[0].vz= 0;
	}else{
		tools_move_to_center_of_momentum();
	}
	//tools_megno_init(1e-16);
	e_init = energy();
	system("rm -f energy.txt");
}

double energy(){
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

void problem_output(){
	if (output_check(10000.)){
		output_timing();
		integrator_synchronize();
		FILE* f = fopen("energy.txt","a");
		double e = energy();
		fprintf(f,"%e %e %e\n",t, fabs((e-e_init)/e_init), tools_megno());
		fclose(f);
	}
}

void problem_finish(){
}
