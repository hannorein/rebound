/**
 * @file 	integrator.c
 * @brief 	IAS15 integrator.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This file implements the IAS15 integration scheme.  
 * IAS stands for Integrator with Adaptive Step-size control, 15th 
 * order. This scheme is a fifteenth order integrator well suited for 
 * high accuracy orbit integration with non-conservative forces.
 * For more details see Rein & Spiegel 2014. Also see Everhart, 1985,
 * ASSL Vol. 115, IAU Colloq. 83, Dynamics of Comets, Their Origin 
 * and Evolution, 185 for the original implementation by Everhart.
 * Part of this code is based a function from the ORSE package.
 * See orsa.sourceforge.net for more details on their implementation.
 *
 * 
 * @section 	LICENSE
 * Copyright (c) 2011-2012 Hanno Rein, Dave Spiegel.
 * Copyright (c) 2002-2004 Pasquale Tricarico.
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
// Uncomment the following line to generate numerical constants with extended precision.
//#define GENERATE_CONSTANTS
#ifdef GENERATE_CONSTANTS
#include <gmp.h>
void integrator_generate_constants(void);
#endif // GENERATE_CONSTANTS
#include "particle.h"
#include "main.h"
#include "gravity.h"
#include "problem.h"
#include "tools.h"
#include "integrator.h"
#include "integrator_ias15.h"

#ifdef TREE
#error IAS15 integrator not working with TREE module.
#endif
#ifdef MPI
#error IAS15 integrator not working with MPI.
#endif

unsigned int integrator_ias15_epsilon_global	= 1;	// if 1: estimate the fractional error by max(acceleration_error)/max(acceleration), where max is take over all particles.
double integrator_ias15_epsilon 		= 1e-9;
double integrator_ias15_min_dt 			= 0;
							// if 0: estimate the fractional error by max(acceleration_error/acceleration).
unsigned long integrator_iterations_max_exceeded= 0;	// Count how many times the iteration did not converge
const double safety_factor 			= 0.25;	// Maximum increase/deacrease of consecutve timesteps.

// Gauss Radau spacings
const double h[8]	= { 0.0, 0.0562625605369221464656521910, 0.1802406917368923649875799428, 0.3526247171131696373739077702, 0.5471536263305553830014485577, 0.7342101772154105410531523211, 0.8853209468390957680903597629, 0.9775206135612875018911745004}; 
// Other constants
const double r[28] = {0.0562625605369221464656522, 0.1802406917368923649875799, 0.1239781311999702185219278, 0.3526247171131696373739078, 0.2963621565762474909082556, 0.1723840253762772723863278, 0.5471536263305553830014486, 0.4908910657936332365357964, 0.3669129345936630180138686, 0.1945289092173857456275408, 0.7342101772154105410531523, 0.6779476166784883945875001, 0.5539694854785181760655724, 0.3815854601022409036792446, 0.1870565508848551580517038, 0.8853209468390957680903598, 0.8290583863021736216247076, 0.7050802551022034031027798, 0.5326962297259261307164520, 0.3381673205085403850889112, 0.1511107696236852270372074, 0.9775206135612875018911745, 0.9212580530243653554255223, 0.7972799218243951369035946, 0.6248958964481178645172667, 0.4303669872307321188897259, 0.2433104363458769608380222, 0.0921996667221917338008147};
const double c[21] = {-0.0562625605369221464656522, 0.0101408028300636299864818, -0.2365032522738145114532321, -0.0035758977292516175949345, 0.0935376952594620658957485, -0.5891279693869841488271399, 0.0019565654099472210769006, -0.0547553868890686864408084, 0.4158812000823068616886219, -1.1362815957175395318285885, -0.0014365302363708915610919, 0.0421585277212687082291130, -0.3600995965020568162530901, 1.2501507118406910366792415, -1.8704917729329500728817408, 0.0012717903090268677658020, -0.0387603579159067708505249, 0.3609622434528459872559689, -1.4668842084004269779203515, 2.9061362593084293206895457, -2.7558127197720458409721005};
const double d[21] = {0.0562625605369221464656522, 0.0031654757181708292499905, 0.2365032522738145114532321, 0.0001780977692217433881125, 0.0457929855060279188954539, 0.5891279693869841488271399, 0.0000100202365223291272096, 0.0084318571535257015445000, 0.2535340690545692665214616, 1.1362815957175395318285885, 0.0000005637641639318207610, 0.0015297840025004658189490, 0.0978342365324440053653648, 0.8752546646840910912297246, 1.8704917729329500728817408, 0.0000000317188154017613665, 0.0002762930909826476593130, 0.0360285539837364596003871, 0.5767330002770787313544596, 2.2485887607691598182153473, 2.7558127197720458409721005};

// The following values will be set dynamically.
double s[9];				// Summation coefficients 

int N3allocated	= 0; 			// Size of allocated arrays.

double* at   	= NULL;			// Temporary buffer for acceleration
double* x0  	= NULL;			// Temporary buffer for position (used for initial values at h=0) 
double* v0  	= NULL;			//                      velocity
double* a0  	= NULL;			//                      acceleration
double* csx  	= NULL;			//                      compensated summation
double* csv  	= NULL;			//                      compensated summation

double* g[7] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL} ;
double* b[7] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL} ;
double* e[7] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL} ;

// The following values are used for resetting the b and e coefficients if a timestep gets rejected
double* br[7] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL} ;
double* er[7] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL} ;
double dt_last_success = 0.;			// Last accepted timestep (corresponding to br and er)
// Helper functions for resetting the b and e coefficients
void copybuffers(double* _a[7], double* _b[7], int N3);
void predict_next_step(double ratio, int N3, double* _e[7], double* _b[7]);

// MEGNO helper routines
// Weights for integration of a first order differential equation (Note: interval length = 2) 
const double w[8] = {0.03125, 0.185358154802979278540728972807180754479812609, 0.304130620646785128975743291458180383736715043, 0.376517545389118556572129261157225608762708603, 0.391572167452493593082499533303669362149363727, 0.347014795634501068709955597003528601733139176, 0.249647901329864963257869294715235590174262844, 0.114508814744257199342353731044292225247093225};

// Do nothing here. This is only used in a leapfrog-like DKD integrator. IAS15 performs one complete timestep.
void integrator_ias15_part1(void){
}

int integrator_ias15_step(void); // Does the actual timestep.

void integrator_ias15_part2(void){
#ifdef GENERATE_CONSTANTS
	integrator_generate_constants();
#endif  // GENERATE_CONSTANTS
	// Try until a step was successful.
	while(!integrator_ias15_step());
}
 
int integrator_ias15_step(void) {
	const int N3 = 3*N;
	if (N3 > N3allocated) {
		for (int l=0;l<7;++l) {
			g[l] = realloc(g[l],sizeof(double)*N3);
			b[l] = realloc(b[l],sizeof(double)*N3);
			e[l] = realloc(e[l],sizeof(double)*N3);
			br[l] = realloc(br[l],sizeof(double)*N3);
			er[l] = realloc(er[l],sizeof(double)*N3);
			for (int k=0;k<N3;k++){
				b[l][k] = 0;
				e[l][k] = 0;
				br[l][k] = 0;
				er[l][k] = 0;
			}
		}
		at = realloc(at,sizeof(double)*N3);
		x0 = realloc(x0,sizeof(double)*N3);
		v0 = realloc(v0,sizeof(double)*N3);
		a0 = realloc(a0,sizeof(double)*N3);
		csx= realloc(csx,sizeof(double)*N3);
		csv= realloc(csv,sizeof(double)*N3);
		for (int i=0;i<N3;i++){
			// Kill compensated summation coefficients
			csx[i] = 0;
			csv[i] = 0;
		}
		N3allocated = N3;
	}
	
	// integrator_update_acceleration(); // Not needed. Forces are already calculated in main routine.

	for(int k=0;k<N;k++) {
		x0[3*k]   = particles[k].x;
		x0[3*k+1] = particles[k].y;
		x0[3*k+2] = particles[k].z;
		v0[3*k]   = particles[k].vx;
		v0[3*k+1] = particles[k].vy;
		v0[3*k+2] = particles[k].vz;
		a0[3*k]   = particles[k].ax;
		a0[3*k+1] = particles[k].ay;  
		a0[3*k+2] = particles[k].az;
	}

	for(int k=0;k<N3;k++) {
		g[0][k] = b[6][k]*d[15] + b[5][k]*d[10] + b[4][k]*d[6] + b[3][k]*d[3]  + b[2][k]*d[1]  + b[1][k]*d[0]  + b[0][k];
		g[1][k] = b[6][k]*d[16] + b[5][k]*d[11] + b[4][k]*d[7] + b[3][k]*d[4]  + b[2][k]*d[2]  + b[1][k];
		g[2][k] = b[6][k]*d[17] + b[5][k]*d[12] + b[4][k]*d[8] + b[3][k]*d[5]  + b[2][k];
		g[3][k] = b[6][k]*d[18] + b[5][k]*d[13] + b[4][k]*d[9] + b[3][k];
		g[4][k] = b[6][k]*d[19] + b[5][k]*d[14] + b[4][k];
		g[5][k] = b[6][k]*d[20] + b[5][k];
		g[6][k] = b[6][k];
	}

	double integrator_megno_thisdt = 0.;
	double integrator_megno_thisdt_init = 0.;
	if (N_megno){
		integrator_megno_thisdt_init = w[0]* t * tools_megno_deltad_delta();
	}

	double t_beginning = t;
	double predictor_corrector_error = 1e300;
	double predictor_corrector_error_last = 2;
	int iterations = 0;	
	// Predictor corrector loop
	// Stops if one of the following conditions is satisfied: 
	//   1) predictor_corrector_error better than 1e-16 
	//   2) predictor_corrector_error starts to oscillate
	//   3) more than 12 iterations
	while(1){
		if(predictor_corrector_error<1e-16){
			break;
		}
		if(iterations > 2 && predictor_corrector_error_last <= predictor_corrector_error){
			break;
		}
		if (iterations>=12){
			integrator_iterations_max_exceeded++;
			const int integrator_iterations_warning = 10;
			if (integrator_iterations_max_exceeded==integrator_iterations_warning ){
				fprintf(stderr,"\n\033[1mWarning!\033[0m At least %d predictor corrector loops in integrator_ias15.c did not converge. This is typically an indication of the timestep being too large.\n",integrator_iterations_warning);
			}
			break;								// Quit predictor corrector loop
		}
		predictor_corrector_error_last = predictor_corrector_error;
		predictor_corrector_error = 0;
		iterations++;

		integrator_megno_thisdt = integrator_megno_thisdt_init;

		for(int n=1;n<8;n++) {							// Loop over interval using Gauss-Radau spacings

			s[0] = dt * h[n];
			s[1] = s[0] * s[0] / 2.;
			s[2] = s[1] * h[n] / 3.;
			s[3] = s[2] * h[n] / 2.;
			s[4] = 3. * s[3] * h[n] / 5.;
			s[5] = 2. * s[4] * h[n] / 3.;
			s[6] = 5. * s[5] * h[n] / 7.;
			s[7] = 3. * s[6] * h[n] / 4.;
			s[8] = 7. * s[7] * h[n] / 9.;
			
			t = t_beginning + s[0];

			// Prepare particles arrays for force calculation
			for(int i=0;i<N;i++) {						// Predict positions at interval n using b values
				const int k0 = 3*i+0;
				const int k1 = 3*i+1;
				const int k2 = 3*i+2;

				double xk0  = csx[k0] + (s[8]*b[6][k0] + s[7]*b[5][k0] + s[6]*b[4][k0] + s[5]*b[3][k0] + s[4]*b[2][k0] + s[3]*b[1][k0] + s[2]*b[0][k0] + s[1]*a0[k0] + s[0]*v0[k0] );
				particles[i].x = xk0 + x0[k0];
				double xk1  = csx[k1] + (s[8]*b[6][k1] + s[7]*b[5][k1] + s[6]*b[4][k1] + s[5]*b[3][k1] + s[4]*b[2][k1] + s[3]*b[1][k1] + s[2]*b[0][k1] + s[1]*a0[k1] + s[0]*v0[k1] );
				particles[i].y = xk1 + x0[k1];
				double xk2  = csx[k2] + (s[8]*b[6][k2] + s[7]*b[5][k2] + s[6]*b[4][k2] + s[5]*b[3][k2] + s[4]*b[2][k2] + s[3]*b[1][k2] + s[2]*b[0][k2] + s[1]*a0[k2] + s[0]*v0[k2] );
				particles[i].z = xk2 + x0[k2];
			}
			if (N_megno || ((problem_additional_forces || problem_additional_forces_with_parameters)&& integrator_force_is_velocitydependent)){
				s[0] = dt * h[n];
				s[1] =      s[0] * h[n] / 2.;
				s[2] = 2. * s[1] * h[n] / 3.;
				s[3] = 3. * s[2] * h[n] / 4.;
				s[4] = 4. * s[3] * h[n] / 5.;
				s[5] = 5. * s[4] * h[n] / 6.;
				s[6] = 6. * s[5] * h[n] / 7.;
				s[7] = 7. * s[6] * h[n] / 8.;

				for(int i=0;i<N;i++) {					// Predict velocities at interval n using b values
					const int k0 = 3*i+0;
					const int k1 = 3*i+1;
					const int k2 = 3*i+2;

					double vk0 =  csv[k0] + s[7]*b[6][k0] + s[6]*b[5][k0] + s[5]*b[4][k0] + s[4]*b[3][k0] + s[3]*b[2][k0] + s[2]*b[1][k0] + s[1]*b[0][k0] + s[0]*a0[k0];
					particles[i].vx = vk0 + v0[k0];
					double vk1 =  csv[k1] + s[7]*b[6][k1] + s[6]*b[5][k1] + s[5]*b[4][k1] + s[4]*b[3][k1] + s[3]*b[2][k1] + s[2]*b[1][k1] + s[1]*b[0][k1] + s[0]*a0[k1];
					particles[i].vy = vk1 + v0[k1];
					double vk2 =  csv[k2] + s[7]*b[6][k2] + s[6]*b[5][k2] + s[5]*b[4][k2] + s[4]*b[3][k2] + s[3]*b[2][k2] + s[2]*b[1][k2] + s[1]*b[0][k2] + s[0]*a0[k2];
					particles[i].vz = vk2 + v0[k2];
				}
			}


			integrator_update_acceleration();				// Calculate forces at interval n
			if (N_megno){
				integrator_megno_thisdt += w[n] * t * tools_megno_deltad_delta();
			}

			for(int k=0;k<N;++k) {
				at[3*k]   = particles[k].ax;
				at[3*k+1] = particles[k].ay;  
				at[3*k+2] = particles[k].az;
			}
			switch (n) {							// Improve b and g values
				case 1: 
					for(int k=0;k<N3;++k) {
						double tmp = g[0][k];
						g[0][k]  = (at[k] - a0[k]) / r[0];
						b[0][k] += g[0][k] - tmp;
					} break;
				case 2: 
					for(int k=0;k<N3;++k) {
						double tmp = g[1][k];
						const double gk = at[k] - a0[k];
						g[1][k] = (gk/r[1] - g[0][k])/r[2];
						tmp = g[1][k] - tmp;
						b[0][k] += tmp * c[0];
						b[1][k] += tmp;
					} break;
				case 3: 
					for(int k=0;k<N3;++k) {
						double tmp = g[2][k];
						const double gk = at[k] - a0[k];
						g[2][k] = ((gk/r[3] - g[0][k])/r[4] - g[1][k])/r[5];
						tmp = g[2][k] - tmp;
						b[0][k] += tmp * c[1];
						b[1][k] += tmp * c[2];
						b[2][k] += tmp;
					} break;
				case 4:
					for(int k=0;k<N3;++k) {
						double tmp = g[3][k];
						const double gk = at[k] - a0[k];
						g[3][k] = (((gk/r[6] - g[0][k])/r[7] - g[1][k])/r[8] - g[2][k])/r[9];
						tmp = g[3][k] - tmp;
						b[0][k] += tmp * c[3];
						b[1][k] += tmp * c[4];
						b[2][k] += tmp * c[5];
						b[3][k] += tmp;
					} break;
				case 5:
					for(int k=0;k<N3;++k) {
						double tmp = g[4][k];
						const double gk = at[k] - a0[k];
						g[4][k] = ((((gk/r[10] - g[0][k])/r[11] - g[1][k])/r[12] - g[2][k])/r[13] - g[3][k])/r[14];
						tmp = g[4][k] - tmp;
						b[0][k] += tmp * c[6];
						b[1][k] += tmp * c[7];
						b[2][k] += tmp * c[8];
						b[3][k] += tmp * c[9];
						b[4][k] += tmp;
					} break;
				case 6:
					for(int k=0;k<N3;++k) {
						double tmp = g[5][k];
						const double gk = at[k] - a0[k];
						g[5][k] = (((((gk/r[15] - g[0][k])/r[16] - g[1][k])/r[17] - g[2][k])/r[18] - g[3][k])/r[19] - g[4][k])/r[20];
						tmp = g[5][k] - tmp;
						b[0][k] += tmp * c[10];
						b[1][k] += tmp * c[11];
						b[2][k] += tmp * c[12];
						b[3][k] += tmp * c[13];
						b[4][k] += tmp * c[14];
						b[5][k] += tmp;
					} break;
				case 7:
				{
					double maxak = 0.0;
					double maxb6ktmp = 0.0;
					for(int k=0;k<N3;++k) {
						double tmp = g[6][k];
						const double gk = at[k] - a0[k];
						g[6][k] = ((((((gk/r[21] - g[0][k])/r[22] - g[1][k])/r[23] - g[2][k])/r[24] - g[3][k])/r[25] - g[4][k])/r[26] - g[5][k])/r[27];
						tmp = g[6][k] - tmp;	
						b[0][k] += tmp * c[15];
						b[1][k] += tmp * c[16];
						b[2][k] += tmp * c[17];
						b[3][k] += tmp * c[18];
						b[4][k] += tmp * c[19];
						b[5][k] += tmp * c[20];
						b[6][k] += tmp;
						
						// Monitor change in b[6][k] relative to at[k]. The predictor corrector scheme is converged if it is close to 0.
						if (integrator_ias15_epsilon_global){
							const double ak  = fabs(at[k]);
							if (isnormal(ak) && ak>maxak){
								maxak = ak;
							}
							const double b6ktmp = fabs(tmp);  // change of b6ktmp coefficient
							if (isnormal(b6ktmp) && b6ktmp>maxb6ktmp){
								maxb6ktmp = b6ktmp;
							}
						}else{
							const double ak  = at[k];
							const double b6ktmp = tmp; 
							const double errork = fabs(b6ktmp/ak);
							if (isnormal(errork) && errork>predictor_corrector_error){
								predictor_corrector_error = errork;
							}
						}
					} 
					if (integrator_ias15_epsilon_global){
						predictor_corrector_error = maxb6ktmp/maxak;
					}
					
					break;
				}
			}
		}
	}
	// Set time back to initial value (will be updated below) 
	t = t_beginning;
	// Find new timestep
	const double dt_done = dt;
	
	if (integrator_ias15_epsilon>0){
		// Estimate error (given by last term in series expansion) 
		// There are two options:
		// integrator_ias15_epsilon_global==1  (default)
		//   First, we determine the maximum acceleration and the maximum of the last term in the series. 
		//   Then, the two are divided.
		// integrator_ias15_epsilon_global==0
		//   Here, the fractional error is calculated for each particle individually and we use the maximum of the fractional error.
		//   This might fail in cases where a particle does not experience any (physical) acceleration besides roundoff errors. 
		double integrator_error = 0.0;
		if (integrator_ias15_epsilon_global){
			double maxak = 0.0;
			double maxb6k = 0.0;
			for(int i=0;i<N;i++){ // Looping over all particles and all 3 components of the acceleration. 
				const double v2 = particles[i].vx*particles[i].vx+particles[i].vy*particles[i].vy+particles[i].vz*particles[i].vz;
				const double x2 = particles[i].x*particles[i].x+particles[i].y*particles[i].y+particles[i].z*particles[i].z;
				// Skip slowly varying accelerations
				if (fabs(v2*dt*dt/x2) < 1e-16) continue;
				for(int k=3*i;k<3*(i+1);k++) { 
					const double ak  = fabs(at[k]);
					if (isnormal(ak) && ak>maxak){
						maxak = ak;
					}
					const double b6k = fabs(b[6][k]); 
					if (isnormal(b6k) && b6k>maxb6k){
						maxb6k = b6k;
					}
				}
			}
			integrator_error = maxb6k/maxak;
		}else{
			for(int k=0;k<N3;k++) {
				const double ak  = at[k];
				const double b6k = b[6][k]; 
				const double errork = fabs(b6k/ak);
				if (isnormal(errork) && errork>integrator_error){
					integrator_error = errork;
				}
			}
		}

		double dt_new;
		if  (isnormal(integrator_error)){ 	
			// if error estimate is available increase by more educated guess
		 	dt_new = pow(integrator_ias15_epsilon/integrator_error,1./7.)*dt_done;
		}else{					// In the rare case that the error estimate doesn't give a finite number (e.g. when all forces accidentally cancel up to machine precission).
		 	dt_new = dt_done/safety_factor; // by default, increase timestep a little
		}
		
		if (fabs(dt_new)<integrator_ias15_min_dt) dt_new = copysign(integrator_ias15_min_dt,dt_new);
		
		if (fabs(dt_new/dt_done) < safety_factor) {	// New timestep is significantly smaller.
			// Reset particles
			for(int k=0;k<N;++k) {
				particles[k].x = x0[3*k+0];	// Set inital position
				particles[k].y = x0[3*k+1];
				particles[k].z = x0[3*k+2];

				particles[k].vx = v0[3*k+0];	// Set inital velocity
				particles[k].vy = v0[3*k+1];
				particles[k].vz = v0[3*k+2];
			}
			dt = dt_new;
			if (dt_last_success!=0.){		// Do not predict next e/b values if this is the first time step.
				double ratio = dt/dt_last_success;
				predict_next_step(ratio, N3, er, br);
			}
			
			return 0; // Step rejected. Do again. 
		}		
		if (fabs(dt_new/dt_done) > 1.0) {	// New timestep is larger.
			if (dt_new/dt_done > 1./safety_factor) dt_new = dt_done /safety_factor;	// Don't increase the timestep by too much compared to the last one.
		}
		dt = dt_new;
	}

	// Find new position and velocity values at end of the sequence
	const double dt_done2 = dt_done * dt_done;
	for(int k=0;k<N3;++k) {
		{
			double a = x0[k];
			csx[k]  +=  (b[6][k]/72. + b[5][k]/56. + b[4][k]/42. + b[3][k]/30. + b[2][k]/20. + b[1][k]/12. + b[0][k]/6. + a0[k]/2.) 
					* dt_done2 + v0[k] * dt_done;
			x0[k]    = a + csx[k];
			csx[k]  += a - x0[k]; 
		}
		{
			double a = v0[k]; 
			csv[k]  += (b[6][k]/8. + b[5][k]/7. + b[4][k]/6. + b[3][k]/5. + b[2][k]/4. + b[1][k]/3. + b[0][k]/2. + a0[k])
					* dt_done;
			v0[k]    = a + csv[k];
			csv[k]  += a - v0[k];
		}
	}

	t += dt_done;

	if (N_megno){
		double dY = dt_done*integrator_megno_thisdt;
		tools_megno_update(dY);
	}

	// Swap particle buffers
	for(int k=0;k<N;++k) {
		particles[k].x = x0[3*k+0];	// Set final position
		particles[k].y = x0[3*k+1];
		particles[k].z = x0[3*k+2];

		particles[k].vx = v0[3*k+0];	// Set final velocity
		particles[k].vy = v0[3*k+1];
		particles[k].vz = v0[3*k+2];
	}
	dt_last_success = dt_done;
	copybuffers(e,er,N3);		
	copybuffers(b,br,N3);		
	double ratio = dt/dt_done;
	predict_next_step(ratio, N3, e, b);
	return 1; // Success.
}

void predict_next_step(double ratio, int N3, double* _e[7], double* _b[7]){
	// Predict new B values to use at the start of the next sequence. The predicted
	// values from the last call are saved as E. The correction, BD, between the
	// actual and predicted values of B is applied in advance as a correction.
	//
	const double q1 = ratio;
	const double q2 = q1 * q1;
	const double q3 = q1 * q2;
	const double q4 = q2 * q2;
	const double q5 = q2 * q3;
	const double q6 = q3 * q3;
	const double q7 = q3 * q4;

	for(int k=0;k<N3;++k) {
		double be0 = _b[0][k] - _e[0][k];
		double be1 = _b[1][k] - _e[1][k];
		double be2 = _b[2][k] - _e[2][k];
		double be3 = _b[3][k] - _e[3][k];
		double be4 = _b[4][k] - _e[4][k];
		double be5 = _b[5][k] - _e[5][k];
		double be6 = _b[6][k] - _e[6][k];


		e[0][k] = q1*(_b[6][k]* 7.0 + _b[5][k]* 6.0 + _b[4][k]* 5.0 + _b[3][k]* 4.0 + _b[2][k]* 3.0 + _b[1][k]*2.0 + _b[0][k]);
		e[1][k] = q2*(_b[6][k]*21.0 + _b[5][k]*15.0 + _b[4][k]*10.0 + _b[3][k]* 6.0 + _b[2][k]* 3.0 + _b[1][k]);
		e[2][k] = q3*(_b[6][k]*35.0 + _b[5][k]*20.0 + _b[4][k]*10.0 + _b[3][k]* 4.0 + _b[2][k]);
		e[3][k] = q4*(_b[6][k]*35.0 + _b[5][k]*15.0 + _b[4][k]* 5.0 + _b[3][k]);
		e[4][k] = q5*(_b[6][k]*21.0 + _b[5][k]* 6.0 + _b[4][k]);
		e[5][k] = q6*(_b[6][k]* 7.0 + _b[5][k]);
		e[6][k] = q7* _b[6][k];
		

		b[0][k] = e[0][k] + be0;
		b[1][k] = e[1][k] + be1;
		b[2][k] = e[2][k] + be2;
		b[3][k] = e[3][k] + be3;
		b[4][k] = e[4][k] + be4;
		b[5][k] = e[5][k] + be5;
		b[6][k] = e[6][k] + be6;
	}
}

void copybuffers(double* _a[7], double* _b[7], int N3){
	for (int i=0;i<N3;i++){	
		_b[0][i] = _a[0][i];
		_b[1][i] = _a[1][i];
		_b[2][i] = _a[2][i];
		_b[3][i] = _a[3][i];
		_b[4][i] = _a[4][i];
		_b[5][i] = _a[5][i];
		_b[6][i] = _a[6][i];
	}
// The above code seems faster than the code below, probably due to some compiler optimizations. 
//	for (int i=0;i<7;i++){	
//		memcpy(_b[i],_a[i], sizeof(double)*N3);
//	}
}
void integrator_ias15_synchronize(void){
}

void integrator_ias15_reset(void){
	N3allocated 	= 0;
	dt_last_success = 0;
	for (int l=0;l<7;++l) {
		free(g[l]);
		g[l] = NULL;
		free(b[l]);
		b[l] = NULL;
		free(e[l]);
		e[l] = NULL;
		free(br[l]);
		br[l] = NULL;
		free(er[l]);
		er[l] = NULL;
	}
	free(at);
	at =  NULL;
	free(x0);
	x0 =  NULL;
	free(v0);
	v0 =  NULL;
	free(a0);
	a0 =  NULL;
	free(csx);
	csx=  NULL;
	free(csv);
	csv=  NULL;
}

#ifdef GENERATE_CONSTANTS
void integrator_generate_constants(void){
	printf("Generaring constants.\n\n");
	mpf_set_default_prec(512);
	mpf_t* _h = malloc(sizeof(mpf_t)*8);
	for (int i=0;i<8;i++){
		mpf_init(_h[i]);
	}
	mpf_t* _r = malloc(sizeof(mpf_t)*28);
	for (int i=0;i<28;i++){
		mpf_init(_r[i]);
	}
	mpf_t* _c = malloc(sizeof(mpf_t)*21);
	mpf_t* _d = malloc(sizeof(mpf_t)*21);
	for (int i=0;i<21;i++){
		mpf_init(_c[i]);
		mpf_init(_d[i]);
	}
	mpf_set_str(_h[0],"0.0",10);
	mpf_set_str(_h[1],"0.0562625605369221464656521910",10);
	mpf_set_str(_h[2],"0.1802406917368923649875799428",10);
	mpf_set_str(_h[3],"0.3526247171131696373739077702",10);
	mpf_set_str(_h[4],"0.5471536263305553830014485577",10);
	mpf_set_str(_h[5],"0.7342101772154105410531523211",10);
	mpf_set_str(_h[6],"0.8853209468390957680903597629",10);
	mpf_set_str(_h[7],"0.9775206135612875018911745004",10);

	int l=0;
	for (int j=1;j<8;++j) {
		for(int k=0;k<j;++k) {
			// r[l] = h[j] - h[k];
			mpf_sub(_r[l],_h[j],_h[k]);
			++l;
		}
	}
	//c[0] = -h[1];
	mpf_neg(_c[0],_h[1]);
	//d[0] =  h[1];
	mpf_set(_d[0],_h[1]);
	l=0;
	for (int j=2;j<7;++j) {
		++l;
		// c[l] = -h[j] * c[l-j+1];
		mpf_mul(_c[l], _h[j], _c[l-j+1]);
		mpf_neg(_c[l], _c[l]);
		//d[l] =  h[1] * d[l-j+1];
		mpf_mul(_d[l], _h[1], _d[l-j+1]);
		for(int k=2;k<j;++k) {
			++l;
			//c[l] = c[l-j] - h[j] * c[l-j+1];
			mpf_mul(_c[l], _h[j], _c[l-j+1]);
			mpf_sub(_c[l], _c[l-j], _c[l]);
			//d[l] = d[l-j] + h[k] * d[l-j+1];
			mpf_mul(_d[l], _h[k], _d[l-j+1]);
			mpf_add(_d[l], _d[l-j], _d[l]);
		}
		++l;
		//c[l] = c[l-j] - h[j];
		mpf_sub(_c[l], _c[l-j], _h[j]);
		//d[l] = d[l-j] + h[j]; 
		mpf_add(_d[l], _d[l-j], _h[j]);
	}

	// Output	
	printf("double r[28] = {");
	for (int i=0;i<28;i++){
	     gmp_printf ("%.*Ff", 25, _r[i]);
	     if (i!=27) printf(", ");
	}
	printf("};\n");
	printf("double c[21] = {");
	for (int i=0;i<21;i++){
	     gmp_printf ("%.*Ff", 25, _c[i]);
	     if (i!=20) printf(", ");
	}
	printf("};\n");
	printf("double d[21] = {");
	for (int i=0;i<21;i++){
	     gmp_printf ("%.*Ff", 25, _d[i]);
	     if (i!=20) printf(", ");
	}
	printf("};\n");
	exit(0);
}
#endif // GENERATE_CONSTANTS
