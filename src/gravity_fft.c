/**
 * @file 	gravity.c
 * @brief 	Gravity calculation on a grid using an FFT.
 * @author 	Hanno Rein <hanno@hanno-rein.de>, Geoffroy Lesur <geoffroy.lesur@obs.ujf-grenoble.fr>
 *
 * @details 	This is a 2D FFT poisson solver for periodic and shearing sheet boxes.
 * 		It has not been well tested yet, so use with caution.
 *		Furthermore, it is not parallelized yet.
 *		The number of grid points is set by root_nx and root_ny. 
 *
 * 
 * @section LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu, Geoffroy Lesur
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
#include "boundaries.h"
#include <fftw3.h>

#ifdef MPI
#error GRAVITY_FFT not compatible with MPI yet
#endif

int grid_NX_COMPLEX;
int grid_NY_COMPLEX;
int grid_NCOMPLEX;
double dx,dy;		/**< Grid spacing */
double* kx;		/**< Wave vector */
double* ky;		/**< Wave vector */
double* kxt;		/**< Time dependent wave vector (shearing sheet only) */
double* k;		/**< Magnitude of wave vector */
double* density;	/**< Complex density field */
double* density_r;	/**< Real density field  */
double* fx;		/**< Force in x direction */
double* fy;		/**< Force in y direction */
fftw_plan r2cfft;	/**< FFT plan real to complex */
fftw_plan c2rfft;	/**< FFT plan complex to real */
#ifdef INTEGRATOR_SEI
double* w1d;		/**< Temporary 1D arrary for remapping (shearing sheet only) */
fftw_plan for1dfft;	/**< FFT plan for remapping (1D, shearing sheet only) */
fftw_plan bac1dfft;	/**< FFT plan for remapping (1D, shearing sheet only) */
#endif 	// INTEGRATOR_SEI

int gravity_fft_init_done = 0;	/**< Flag if arrays and plans are initialized */
void gravity_fft_init();
void gravity_fft_grid2p(struct particle* p);
void gravity_fft_p2grid();
void gravity_fft_remap(double* wi, const double direction);
double shift_shear = 0;

void gravity_calculate_acceleration(){
	// Setting up the grid
	if (gravity_fft_init_done==0){
		gravity_fft_init();
		gravity_fft_init_done=1;
	}
#pragma omp parallel for schedule(guided)
	for (int i=0; i<N; i++){
		particles[i].ax = 0; 
		particles[i].ay = 0; 
		particles[i].az = 0; 
	}
#ifdef INTEGRATOR_SEI
	struct ghostbox gb = boundaries_get_ghostbox(1,0,0);
	shift_shear = gb.shifty;
#endif 	// INTEGRATOR_SEI
	gravity_fft_p2grid();
	
#ifdef INTEGRATOR_SEI
	// Remap in fourier space to deal with shearing sheet boundary conditions.
	gravity_fft_remap(density_r, 1);
#endif 	// INTEGRATOR_SEI
	
	fftw_execute_dft_r2c(r2cfft, density_r, (fftw_complex*)density);
	
	
	// Inverse Poisson equation
	
	for(int i = 0 ; i < grid_NCOMPLEX ; i++) {
#ifdef INTEGRATOR_SEI
		// Compute time-dependent wave-vectors
		kxt[i] = kx[i] + shift_shear/boxsize_y * ky[i];
		k[i]  = sqrt( kxt[i]*kxt[i] + ky[i] * ky[i]);
		// we will use 1/k, that prevents singularity 
		// (the k=0 is set to zero by renormalization...)
		if ( k[i] == 0.0 ) k[i] = 1.0; 
#endif 	// INTEGRATOR_SEI
		double q0 = - 2.0 * M_PI * density[2*i] / (k[i] * root_nx * root_ny);
		double q1 = - 2.0 * M_PI * density[2*i+1] / (k[i] * root_nx * root_ny);
		double sinkxt = sin(kxt[i] * dx);
		double sinky  = sin(ky[i] * dy);
		fx[2*i]		=   q1 * sinkxt / dx;		// Real part of Fx
		fx[2*i+1] 	= - q0 * sinkxt / dx;		// Imaginary part of Fx
		fy[2*i]		=   q1 * sinky  / dy;	
		fy[2*i+1] 	= - q0 * sinky  / dy;
	}
	
	// Transform back the force field
	fftw_execute_dft_c2r(c2rfft, (fftw_complex*)fx, fx);
	fftw_execute_dft_c2r(c2rfft, (fftw_complex*)fy, fy);
	
#ifdef INTEGRATOR_SEI
	// Remap in fourier space to deal with shearing sheet boundary conditions.
	gravity_fft_remap(fx, -1);
	gravity_fft_remap(fy, -1);
#endif	// INTEGRATOR_SEI

	for(int i=0;i<N;i++){
		gravity_fft_grid2p(&(particles[i]));
	}
}

void gravity_fft_init() {

	// dimension definition
	grid_NX_COMPLEX	= root_nx;		
	grid_NY_COMPLEX	= (root_ny / 2 + 1);
	grid_NCOMPLEX	= grid_NX_COMPLEX * grid_NY_COMPLEX;
	dx		= boxsize_x / root_nx;		
	dy		= boxsize_y / root_ny;	
	
	// Array allocation
	kx  = (double *) fftw_malloc( sizeof(double) * grid_NCOMPLEX);
	ky  = (double *) fftw_malloc( sizeof(double) * grid_NCOMPLEX);
#ifdef INTEGRATOR_SEI
	kxt = (double *) fftw_malloc( sizeof(double) * grid_NCOMPLEX);
	w1d = (double *) fftw_malloc( sizeof(double) * root_ny * 2 );
#else	// INTEGRATOR_SEI
	kxt = kx; 	// No time dependent wave vectors.
#endif 	// INTEGRATOR_SEI
	k   = (double *) fftw_malloc( sizeof(double) * grid_NCOMPLEX);
	density = (double *) fftw_malloc( sizeof(double) * grid_NCOMPLEX * 2);
	density_r = (double *) fftw_malloc( sizeof(double) * grid_NCOMPLEX * 2);
	fx  = (double *) fftw_malloc( sizeof(double) * grid_NCOMPLEX * 2);
	fy  = (double *) fftw_malloc( sizeof(double) * grid_NCOMPLEX * 2);
	
	// Init wavevectors
	
	for(int i = 0; i < grid_NX_COMPLEX; i++) {
		for(int j =0; j < grid_NY_COMPLEX; j++) {
			int IDX2D = i * grid_NY_COMPLEX + j;
			kx[IDX2D] = (2.0 * M_PI) / boxsize_x * (
					fmod( (double) (i + (grid_NX_COMPLEX/2.0 )), (double) grid_NX_COMPLEX)
					 - (double) grid_NX_COMPLEX / 2.0 );
			ky[IDX2D] = (2.0 * M_PI) / boxsize_y * ((double) j);
#ifndef INTEGRATOR_SEI
			k[IDX2D]  = pow( kx[IDX2D]*kx[IDX2D] + ky[IDX2D] * ky[IDX2D], 0.5);
			
			// we will use 1/k, that prevents singularity 
			// (the k=0 is set to zero by renormalization...)
			if ( k[IDX2D] == 0.0 ) k[IDX2D] = 1.0; 
#endif	// INTEGRATOR_SEI
		}
	}
	
	// Init ffts (use in place fourier transform for efficient memory usage)
	r2cfft = fftw_plan_dft_r2c_2d( root_nx, root_ny, density, (fftw_complex*)density, FFTW_MEASURE);
	c2rfft = fftw_plan_dft_c2r_2d( root_nx, root_ny, (fftw_complex*)density, density, FFTW_MEASURE);
#ifdef INTEGRATOR_SEI
	for1dfft = fftw_plan_dft_1d(root_ny, (fftw_complex*)w1d, (fftw_complex*)w1d, FFTW_FORWARD, FFTW_MEASURE);
	bac1dfft = fftw_plan_dft_1d(root_ny, (fftw_complex*)w1d, (fftw_complex*)w1d, FFTW_BACKWARD, FFTW_MEASURE);
#endif	// INTEGRATOR_SEI
}

// Assignement function (TSC Scheme) 
// See Hockney and Eastwood (1981), Computer Simulations Using Particles
double W(double x){
	if (fabs(x)<=0.5) return 0.75 - x*x;
	if (fabs(x)>=0.5 && fabs(x)<=3./2.) return 0.5*(3./2.-fabs(x))*(3./2.-fabs(x));
	return 0; 
}

#ifdef INTEGRATOR_SEI
void gravity_fft_remap(double* wi, const double direction) {
	double phase, rew, imw;
	
	for(int i = 0 ; i < root_nx ; i++) {
		for(int j = 0 ; j < root_ny ; j++) {
			w1d[ 2 * j ] = wi[j + (root_ny + 2) * i];		// w1d is supposed to be a complex array. 
			w1d[ 2 * j + 1 ] = 0.0;
		}
		
		// Transform w1d, which will be stored in w2d
		fftw_execute(for1dfft);
					
		for(int j = 0 ; j < root_ny ; j++) {
			// phase = ky * (-shift_shear)
			phase =  - direction * (2.0 * M_PI) / boxsize_y * ((j + (root_ny / 2)) % root_ny - root_ny / 2) * shift_shear * ((double) i) / ((double) root_nx);
			
			rew = w1d[2 * j];
			imw = w1d[2 * j + 1];
			
			// Real part
			w1d[2 * j    ] = rew * cos(phase) - imw * sin(phase);
			// Imaginary part
			w1d[2 * j + 1] = rew * sin(phase) + imw * cos(phase);
			
			// Throw the Nyquist Frequency (should be useless anyway)
			if(j==root_ny/2) {
				w1d[2 * j    ] =0.0;
				w1d[2 * j + 1] = 0.0;
			}
		}
		
		fftw_execute(bac1dfft);
		
		for(int j = 0 ; j < root_ny ; j++) {
			wi[j + (root_ny + 2) * i] = w1d[ 2 * j ] / root_ny;
		}
	}
}
#endif // INTEGRATOR_SEI

void gravity_fft_p2grid(){
		
	// clean the current density
	for(int i = 0 ; i < root_nx * (root_ny + 2) ; i++) {
		density_r[i] = 0.0;			// density is used to store the surface density
	}
	
	for (int i=0; i<N; i++){
		struct particle p = particles[i];
		// I'm sorry to say I have to keep these traps. Something's wrong if these traps are called.
		
		int x = (int) floor((p.x / boxsize_x + 0.5) * root_nx);
		int y = (int) floor((p.y / boxsize_y + 0.5) * root_ny);
		
		// Formally, pos.x is in the interval [-size/2 , size/2 [. Therefore, x and y should be in [0 , grid_NJ-1]
		
		// xp1, xm1... might be out of bound. They are however the relevant coordinates for the interpolation.

		int xp1 = x + 1;
		int xm1 = x - 1;
		int ym1 = y - 1;
		int yp1 = y + 1;

		
		// Target according to boundary conditions.
		// Although xTarget and yTarget are not relevant here, they will be relevant with shear
		// We have to use all these fancy variables since y depends on x because of the shearing path
		// Any nicer solution is welcome
		
		int xTarget = x;
		int xp1Target = xp1;
		int xm1Target = xm1;
		
		int ym1_xm1Target = (ym1 + root_ny) % root_ny;
		int ym1_xTarget   = ym1_xm1Target;
		int ym1_xp1Target = ym1_xm1Target;
		
		int y_xm1Target = y % root_ny;
		int y_xTarget   = y_xm1Target;
		int y_xp1Target = y_xm1Target;
		
		int yp1_xm1Target = yp1 % root_ny;
		int yp1_xTarget   = yp1_xm1Target;
		int yp1_xp1Target = yp1_xm1Target;
		
		double tx, ty;
		double q0 =  G* p.m /(dx*dy);
		
		// Shearing patch trick
		// This is only an **approximate** mapping
		// one should use an exact interpolation scheme here (Fourier like).
		
		if(xp1Target>=root_nx) {
			xp1Target -= root_nx;							// X periodicity
			y_xp1Target = y_xp1Target + round((shift_shear/boxsize_y) * root_ny);
			y_xp1Target = (y_xp1Target + root_ny) % root_ny;		// Y periodicity
			yp1_xp1Target = yp1_xp1Target + round((shift_shear/boxsize_y) * root_ny);
			yp1_xp1Target = (yp1_xp1Target + root_ny) % root_ny;
			ym1_xp1Target = ym1_xp1Target + round((shift_shear/boxsize_y) * root_ny);
			ym1_xp1Target = (ym1_xp1Target + root_ny) % root_ny;
		}
		
		if(xm1Target<0) {
			xm1Target += root_nx;
			y_xm1Target = y_xm1Target - round((shift_shear/boxsize_x) * root_ny);
			y_xm1Target = (y_xm1Target + root_ny) % root_ny;		// Y periodicity
			yp1_xm1Target = yp1_xm1Target - round((shift_shear/boxsize_x) * root_ny);
			yp1_xm1Target = (yp1_xm1Target + root_ny) % root_ny;
			ym1_xm1Target = ym1_xm1Target - round((shift_shear/boxsize_x) * root_ny);
			ym1_xm1Target = (ym1_xm1Target + root_ny) % root_ny;
		}

		// Distribute density to the 9 nearest cells
		
		tx = ((double)xm1 +0.5) * boxsize_x / root_nx -0.5*boxsize_x - p.x;
		ty = ((double)ym1 +0.5) * boxsize_y / root_ny -0.5*boxsize_y - p.y;
		density_r[(root_ny+2) * xm1Target + ym1_xm1Target] += q0 * W(tx/dx)*W(ty/dy); 

		tx = ((double)x   +0.5) * boxsize_x / root_nx -0.5*boxsize_x - p.x;
		ty = ((double)ym1 +0.5) * boxsize_y / root_ny -0.5*boxsize_y - p.y;
		density_r[(root_ny+2) * xTarget   + ym1_xTarget] += q0 * W(tx/dx)*W(ty/dy); 
		
		tx = ((double)xp1 +0.5) * boxsize_x / root_nx -0.5*boxsize_x - p.x;
		ty = ((double)ym1 +0.5) * boxsize_y / root_ny -0.5*boxsize_y - p.y;
		density_r[(root_ny+2) * xp1Target + ym1_xp1Target] += q0 * W(tx/dx)*W(ty/dy); 

		tx = ((double)xm1 +0.5) * boxsize_x / root_nx -0.5*boxsize_x - p.x;
		ty = ((double)y   +0.5) * boxsize_y / root_ny -0.5*boxsize_y - p.y;
		density_r[(root_ny+2) * xm1Target + y_xm1Target  ] += q0 * W(tx/dx)*W(ty/dy); 

		tx = ((double)x   +0.5) * boxsize_x / root_nx -0.5*boxsize_x - p.x;
		ty = ((double)y   +0.5) * boxsize_y / root_ny -0.5*boxsize_y - p.y;
		density_r[(root_ny+2) * xTarget   + y_xTarget  ] += q0 * W(tx/dx)*W(ty/dy); 

		tx = ((double)xp1 +0.5) * boxsize_x / root_nx -0.5*boxsize_x - p.x;
		ty = ((double)y   +0.5) * boxsize_y / root_ny -0.5*boxsize_y - p.y;
		density_r[(root_ny+2) * xp1Target + y_xp1Target  ] += q0 * W(tx/dx)*W(ty/dy); 

		tx = ((double)xm1 +0.5) * boxsize_x / root_nx -0.5*boxsize_x - p.x;
		ty = ((double)yp1 +0.5) * boxsize_y / root_ny -0.5*boxsize_y - p.y;
		density_r[(root_ny+2) * xm1Target + yp1_xm1Target] += q0 * W(tx/dx)*W(ty/dy); 

		tx = ((double)x   +0.5) * boxsize_x / root_nx -0.5*boxsize_x - p.x;
		ty = ((double)yp1 +0.5) * boxsize_y / root_ny -0.5*boxsize_y - p.y;
		density_r[(root_ny+2) * xTarget   + yp1_xTarget] += q0 * W(tx/dx)*W(ty/dy); 

		tx = ((double)xp1 +0.5) * boxsize_x / root_nx -0.5*boxsize_x - p.x;
		ty = ((double)yp1 +0.5) * boxsize_y / root_ny -0.5*boxsize_y - p.y;
		density_r[(root_ny+2) * xp1Target + yp1_xp1Target] += q0 * W(tx/dx)*W(ty/dy); 
	}
}


void gravity_fft_grid2p(struct particle* p){
	
	// I'm sorry to say I have to keep these traps. Something's wrong if these traps are called.
		
	int x = (int) floor((p->x / boxsize_x + 0.5) * root_nx);
	int y = (int) floor((p->y / boxsize_y + 0.5) * root_ny);
		
	// Formally, pos.x is in the interval [-size/2 , size/2 [. Therefore, x and y should be in [0 , grid_NJ-1]
		
	
	// xp1, xm1... might be out of bound. They are however the relevant coordinates for the interpolation.
	int xp1 = x + 1;
	int xm1 = x - 1;
	int ym1 = y - 1;
	int yp1 = y + 1;

		
	// Target according to boundary conditions.
	// Although xTarget and yTarget are not relevant here, they will be relevant with shear
	// We have to use all these fancy variables since y depends on x because of the shearing path
	// Any nicer solution is welcome
	
	int xTarget = x;
	int xp1Target = xp1;
	int xm1Target = xm1;
	
	int ym1_xm1Target = (ym1 + root_ny) % root_ny;
	int ym1_xTarget   = ym1_xm1Target;
	int ym1_xp1Target = ym1_xm1Target;
		
	int y_xm1Target = y % root_ny;
	int y_xTarget   = y_xm1Target;
	int y_xp1Target = y_xm1Target;
		
	int yp1_xm1Target = yp1 % root_ny;
	int yp1_xTarget   = yp1_xm1Target;
	int yp1_xp1Target = yp1_xm1Target;


	double tx, ty;
	
	// Shearing patch trick
	// This is only an **approximate** mapping
	// one should use an exact interpolation scheme here (Fourier like).
		
	if(xp1Target>=root_nx) {
		xp1Target -= root_nx;							// X periodicity
		y_xp1Target = y_xp1Target + round((shift_shear/boxsize_y) * root_ny);
		y_xp1Target = (y_xp1Target + root_ny) % root_ny;			// Y periodicity
		yp1_xp1Target = yp1_xp1Target + round((shift_shear/boxsize_y) * root_ny);
		yp1_xp1Target = (yp1_xp1Target + root_ny) % root_ny;
		ym1_xp1Target = ym1_xp1Target + round((shift_shear/boxsize_y) * root_ny);
		ym1_xp1Target = (ym1_xp1Target + root_ny) % root_ny;
	}
		
	if(xm1Target<0) {
		xm1Target += root_nx;
		y_xm1Target = y_xm1Target - round((shift_shear/boxsize_y) * root_ny);
		y_xm1Target = (y_xm1Target + root_ny) % root_ny;			// Y periodicity
		yp1_xm1Target = yp1_xm1Target - round((shift_shear/boxsize_y) * root_ny);
		yp1_xm1Target = (yp1_xm1Target + root_ny) % root_ny;
		ym1_xm1Target = ym1_xm1Target - round((shift_shear/boxsize_y) * root_ny);
		ym1_xm1Target = (ym1_xm1Target + root_ny) % root_ny;
	}


	tx = ((double)xm1 +0.5) * boxsize_x / root_nx -0.5*boxsize_x - p->x;
	ty = ((double)ym1 +0.5) * boxsize_y / root_ny -0.5*boxsize_y - p->y;

	p->ax += fx[(root_ny+2) * xm1Target + ym1_xm1Target] * W(-tx/dx)*W(-ty/dy);
	p->ay += fy[(root_ny+2) * xm1Target + ym1_xm1Target] * W(-tx/dx)*W(-ty/dy);

	tx = ((double)x   +0.5) * boxsize_x / root_nx -0.5*boxsize_x - p->x;
	ty = ((double)ym1 +0.5) * boxsize_y / root_ny -0.5*boxsize_y - p->y;
	
	p->ax += fx[(root_ny+2) * xTarget   + ym1_xTarget] * W(-tx/dx)*W(-ty/dy);
	p->ay += fy[(root_ny+2) * xTarget   + ym1_xTarget] * W(-tx/dx)*W(-ty/dy);
	
	tx = ((double)xp1 +0.5) * boxsize_x / root_nx -0.5*boxsize_x - p->x;
	ty = ((double)ym1 +0.5) * boxsize_y / root_ny -0.5*boxsize_y - p->y;
	
	p->ax += fx[(root_ny+2) * xp1Target + ym1_xp1Target] * W(-tx/dx)*W(-ty/dy);
	p->ay += fy[(root_ny+2) * xp1Target + ym1_xp1Target] * W(-tx/dx)*W(-ty/dy);

	tx = ((double)xm1 +0.5) * boxsize_x / root_nx -0.5*boxsize_x - p->x;
	ty = ((double)y   +0.5) * boxsize_y / root_ny -0.5*boxsize_y - p->y;
	
	p->ax += fx[(root_ny+2) * xm1Target + y_xm1Target  ] * W(-tx/dx)*W(-ty/dy);
	p->ay += fy[(root_ny+2) * xm1Target + y_xm1Target  ] * W(-tx/dx)*W(-ty/dy);

	tx = ((double)x   +0.5) * boxsize_x / root_nx -0.5*boxsize_x - p->x;
	ty = ((double)y   +0.5) * boxsize_y / root_ny -0.5*boxsize_y - p->y;

	p->ax += fx[(root_ny+2) * xTarget   + y_xTarget  ] * W(-tx/dx)*W(-ty/dy);
	p->ay += fy[(root_ny+2) * xTarget   + y_xTarget  ] * W(-tx/dx)*W(-ty/dy); 

	tx = ((double)xp1 +0.5) * boxsize_x / root_nx -0.5*boxsize_x - p->x;
	ty = ((double)y   +0.5) * boxsize_y / root_ny -0.5*boxsize_y - p->y;
	
	p->ax += fx[(root_ny+2) * xp1Target + y_xp1Target  ] * W(-tx/dx)*W(-ty/dy);
	p->ay += fy[(root_ny+2) * xp1Target + y_xp1Target  ] * W(-tx/dx)*W(-ty/dy); 

	tx = ((double)xm1 +0.5) * boxsize_x / root_nx -0.5*boxsize_x - p->x;
	ty = ((double)yp1 +0.5) * boxsize_y / root_ny -0.5*boxsize_y - p->y;

	p->ax += fx[(root_ny+2) * xm1Target + yp1_xm1Target] * W(-tx/dx)*W(-ty/dy);
	p->ay += fy[(root_ny+2) * xm1Target + yp1_xm1Target] * W(-tx/dx)*W(-ty/dy);  

	tx = ((double)x   +0.5) * boxsize_x / root_nx -0.5*boxsize_x - p->x;
	ty = ((double)yp1 +0.5) * boxsize_y / root_ny -0.5*boxsize_y - p->y;

	p->ax += fx[(root_ny+2) * xTarget   + yp1_xTarget] * W(-tx/dx)*W(-ty/dy);
	p->ay += fy[(root_ny+2) * xTarget   + yp1_xTarget] * W(-tx/dx)*W(-ty/dy);  

	tx = ((double)xp1 +0.5) * boxsize_x / root_nx -0.5*boxsize_x - p->x;
	ty = ((double)yp1 +0.5) * boxsize_y / root_ny -0.5*boxsize_y - p->y;

	p->ax += fx[(root_ny+2) * xp1Target + yp1_xp1Target] * W(-tx/dx)*W(-ty/dy);
	p->ay += fy[(root_ny+2) * xp1Target + yp1_xp1Target] * W(-tx/dx)*W(-ty/dy);  
	
}


