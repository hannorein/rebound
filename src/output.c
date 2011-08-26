/**
 * @file 	output.c
 * @brief 	Output routines.
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
#include <sys/time.h>
#include "particle.h"
#include "main.h"
#include "output.h"
#include "communication_mpi.h"
#ifdef OPENGL
#include "display.h"
#ifdef LIBPNG
#include <png.h>
#endif // LIBPNG
#ifdef _APPLE
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif  // _APPLE
#endif  // OPENGL
#ifdef MPI
#include "mpi.h"
#endif // MPI

#ifdef INTEGRATOR_SEI 	// Shearing sheet
extern double OMEGA;
#endif

// Check if output is needed

int output_check(double interval){
	return output_check_phase(interval,0);
}

int output_check_phase(double interval,double phase){
	double shift = t+interval*phase;
	if (floor(shift/interval)!=floor((shift-dt)/interval)){
		return 1;
	}
	// Output at beginning or end of simulation
	if (t==0||t==tmax){
		return 1;
	}
	return 0;
}


/**
 * 3D vector struct.
 */
struct vec3 {
	double x;
	double y;
	double z;
};

/**
 * Struct representing a Keplerian orbit.
 */
struct opv_orbit {
	double a;
	double e;
	double inc;
	double Omega; 	// longitude of ascending node
	double omega; 	// argument of perihelion
	double M;  	// mean anomaly
	double E;  	// eccentric anomaly
	double f; 	// true anomaly
};

void posvel2orbit(struct opv_orbit* o, struct particle* pv, double gmsum);


double output_timing_last = -1; 	/**< Time when output_timing() was called the last time. */
void output_timing(){
#ifdef MPI
	int N_tot = 0;
	MPI_Reduce(&N, &N_tot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
	if (mpi_id!=0) return;
#else
	int N_tot = N;
#endif
	struct timeval tim;
	gettimeofday(&tim, NULL);
	double temp = tim.tv_sec+(tim.tv_usec/1000000.0);
	if (output_timing_last==-1){
		output_timing_last = temp;
	}else{
		printf("\r");
	}
	if (tmax>0){
		printf("N_tot= %- 9d  t= %- 9f  cpu= %- 9f s  t/tmax= %5.2f%%",N_tot,t,temp-output_timing_last,t/tmax*100);
	}else{
		printf("N_tot= %- 9d  t= %- 9f  cpu= %- 9f s ",N_tot,t,temp-output_timing_last);
	}
	fflush(stdout);
	output_timing_last = temp;
}

void output_append_ascii(char* filename){
#ifdef MPI
	char filename_mpi[1024];
	sprintf(filename_mpi,"%s_%d",filename,mpi_id);
	FILE* of = fopen(filename_mpi,"a"); 
#else // MPI
	FILE* of = fopen(filename,"a"); 
#endif // MPI
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
		fprintf(of,"%e\t%e\t%e\t%e\t%e\t%e\n",p.x,p.y,p.z,p.vx,p.vy,p.vz);
	}
	fclose(of);
}

void output_ascii(char* filename){
#ifdef MPI
	char filename_mpi[1024];
	sprintf(filename_mpi,"%s_%d",filename,mpi_id);
	FILE* of = fopen(filename_mpi,"w"); 
#else // MPI
	FILE* of = fopen(filename,"w"); 
#endif // MPI
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
		fprintf(of,"%e\t%e\t%e\t%e\t%e\t%e\n",p.x,p.y,p.z,p.vx,p.vy,p.vz);
	}
	fclose(of);
}

void output_orbits(char* filename){
#ifdef MPI
	char filename_mpi[1024];
	sprintf(filename_mpi,"%s_%d",filename,mpi_id);
	FILE* of = fopen(filename_mpi,"w"); 
#else // MPI
	FILE* of = fopen(filename,"w"); 
#endif // MPI
	for (int i=1;i<N;i++){
		struct particle* p = &(particles[i]);
		struct opv_orbit o;
		posvel2orbit(&o,p,1.);
		fprintf(of,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",o.a,o.e,o.inc,o.Omega,o.omega,o.M,o.E,o.f);
	}
	fclose(of);
}


void output_binary(char* filename){
#ifdef MPI
	char filename_mpi[1024];
	sprintf(filename_mpi,"%s_%d",filename,mpi_id);
	FILE* of = fopen(filename_mpi,"wb"); 
#else // MPI
	FILE* of = fopen(filename,"wb"); 
#endif // MPI
	fwrite(&N,sizeof(int),1,of);
	fwrite(&t,sizeof(double),1,of);
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
		fwrite(&(p),sizeof(struct particle),1,of);
	}
	fclose(of);
}

void output_binary_positions(char* filename){
#ifdef MPI
	char filename_mpi[1024];
	sprintf(filename_mpi,"%s_%d",filename,mpi_id);
	FILE* of = fopen(filename_mpi,"wb"); 
#else // MPI
	FILE* of = fopen(filename,"wb"); 
#endif // MPI
	for (int i=0;i<N;i++){
		struct vec3 v;
		v.x = particles[i].x;
		v.y = particles[i].y;
		v.z = particles[i].z;
		fwrite(&(v),sizeof(struct vec3),1,of);
	}
	fclose(of);
}

void output_append_velocity_dispersion(char* filename){
	// Algorithm with reduced roundoff errors (see wikipedia)
	struct vec3 A = {.x=0, .y=0, .z=0};
	struct vec3 Q = {.x=0, .y=0, .z=0};
	for (int i=0;i<N;i++){
		struct vec3 Aim1 = A;
		struct particle p = particles[i];
		A.x = A.x + (p.vx-A.x)/(double)(i+1);
#ifdef INTEGRATOR_SEI 	// Shearing sheet
		A.y = A.y + (p.vy+1.5*OMEGA*p.x-A.y)/(double)(i+1);
#else
		A.y = A.y + (p.vy-A.y)/(double)(i+1);
#endif
		A.z = A.z + (p.vz-A.z)/(double)(i+1);
		Q.x = Q.x + (p.vx-Aim1.x)*(p.vx-A.x);
#ifdef INTEGRATOR_SEI 	// Shearing sheet
		Q.y = Q.y + (p.vy+1.5*OMEGA*p.x-Aim1.y)*(p.vy+1.5*OMEGA*p.x-A.y);
#else
		Q.y = Q.y + (p.vy-Aim1.y)*(p.vy-A.y);
#endif
		Q.z = Q.z + (p.vz-Aim1.z)*(p.vz-A.z);
	}
#ifdef MPI
	int N_tot = 0;
	struct vec3 A_tot = {.x=0, .y=0, .z=0};
	struct vec3 Q_tot = {.x=0, .y=0, .z=0};
	MPI_Reduce(&N, &N_tot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
	MPI_Reduce(&A, &A_tot, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
	MPI_Reduce(&Q, &Q_tot, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
	if (mpi_id!=0) return;
#else
	int N_tot = N;
	struct vec3 A_tot = A;
	struct vec3 Q_tot = Q;
#endif
	Q_tot.x = sqrt(Q_tot.x/(double)N_tot);
	Q_tot.y = sqrt(Q_tot.y/(double)N_tot);
	Q_tot.z = sqrt(Q_tot.z/(double)N_tot);
	FILE* of = fopen(filename,"a"); 
	fprintf(of,"%e\t%e\t%e\t%e\t%e\t%e\t%e\n",t,A_tot.x,A_tot.y,A_tot.z,Q_tot.x,Q_tot.y,Q_tot.z);
	fclose(of);
}

#ifdef OPENGL
#ifdef LIBPNG
unsigned char* 	imgdata = NULL;
int output_png_num = 0;
void output_png(char* dirname){
	char filename[1024];
	sprintf(filename,"%s%09d.png",dirname,output_png_num);
	output_png_num++;
	output_png_single(filename);
}

void output_png_single(char* filename){
	// Read Image
	if (display_init_done==0) return;
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	int width = viewport[2];
	int height = viewport[3];
	glReadBuffer(GL_BACK);
	//glReadBuffer(GL_FRONT);
	if (imgdata==NULL){
		imgdata = calloc(width*height*3,sizeof(unsigned char));
	}
	png_byte* row_pointers[height];
	for (int h = 0; h < height; h++) {
		row_pointers[height-h-1] = (png_bytep) &imgdata[width*3*h];
	}

	glReadPixels(0,0,width,height,GL_RGB,GL_UNSIGNED_BYTE,imgdata);

	/* open the file */
	FILE *fp;
	png_structp png_ptr;
	png_infop info_ptr;
	fp = fopen(filename, "wb");
	if (fp == NULL) return;

	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

	if (png_ptr == NULL) {
		fclose(fp);
		return;
	}

	/* Allocate/initialize the image information data.  REQUIRED */
	info_ptr = png_create_info_struct(png_ptr);
	if (info_ptr == NULL) {
		fclose(fp);
		png_destroy_write_struct(&png_ptr,  (png_infopp)NULL);
		return;
	}

	/* Set error handling.  REQUIRED if you aren't supplying your own
	* error hadnling functions in the png_create_write_struct() call.
	*/
	/*
	if (setjmp(png_ptr->jmpbuf))
	{
	fclose(fp);
	png_destroy_write_struct(&png_ptr,  (png_infopp)NULL);
	return;
	}
	*/

	/* I/O initialization functions is REQUIRED */
	/* set up the output control if you are using standard C streams */
	png_init_io(png_ptr, fp);

	/* Set the image information here.  Width and height are up to 2^31,
	* bit_depth is one of 1, 2, 4, 8, or 16, but valid values also depend on
	* the color_type selected. color_type is one of PNG_COLOR_TYPE_GRAY,
	* PNG_COLOR_TYPE_GRAY_ALPHA, PNG_COLOR_TYPE_PALETTE, PNG_COLOR_TYPE_RGB,
	* or PNG_COLOR_TYPE_RGB_ALPHA.  interlace is either PNG_INTERLACE_NONE or
	* PNG_INTERLACE_ADAM7, and the compression_type and filter_type MUST
	* currently be PNG_COMPRESSION_TYPE_BASE and PNG_FILTER_TYPE_BASE. REQUIRED
	*/
	png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB,
	PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

	/* Write the file  */
	png_write_info(png_ptr, info_ptr);
	png_write_image(png_ptr, row_pointers);
	png_write_end(png_ptr, info_ptr);

	/* if you allocated any text comments, free them here */

	/* clean up after the write, and free any memory allocated */
	png_destroy_write_struct(&png_ptr, (png_infopp)NULL);

	/* close the file */
	fclose(fp);
}
#endif
#endif

#define TINY 1.0e-12
void posvel2orbit(struct opv_orbit* o, struct particle* pv, double gmsum){
	// Compute the angular momentum H, and thereby the inclination INC.
	double Omega, a, e, M, E=0, f=0, omega, inc; 	// orbital paramaers

	double u; int ialpha;			// internals
	
	double hx = pv->y*pv->vz - pv->z*pv->vy;
	double hy = pv->z*pv->vx - pv->x*pv->vz;
	double hz = pv->x*pv->vy - pv->y*pv->vx;
	double h2 = hx*hx + hy*hy +hz*hz;
	double h  = sqrt(h2);
	inc = acos(hz/h);

	// Compute longitude of ascending node CAPOM and the argument of latitude u.
	double fac = sqrt(hx*hx + hy*hy)/h;
	
	if(fac < TINY ){
	  	Omega = 0.;
	  	u = atan2(pv->y,pv->x);
	  	if(fabs(inc - M_PI) < 10.*TINY){
			u = -u;
	  	}
	}else{
	  	Omega = atan2(hx,-hy); 
	  	u = atan2 ( pv->z/sin(inc) , pv->x*cos(Omega) + pv->y*sin(Omega));
	}

	while(Omega < 0.) Omega = Omega + 2.*M_PI;
	
	while(u < 0.) u = u + 2.*M_PI;
	

	//  Compute the radius R and velocity squared V2, and the dot product RDOTV, the energy per unit mass ENERGY .

	double r = sqrt(pv->x*pv->x + pv->y*pv->y + pv->z*pv->z);
	double v2 = pv->vx*pv->vx + pv->vy*pv->vy + pv->vz*pv->vz;
	double vdotr = pv->x*pv->vx + pv->y*pv->vy + pv->z*pv->vz;
	double energy = 0.5*v2 - gmsum/r;

	//  Determine type of conic section and label it via IALPHA
	if(fabs(energy*r/gmsum) < sqrt(TINY)){
		ialpha = 0;
	}else{
	   	if(energy < 0.) ialpha = -1;
	   	if(energy > 0.) ialpha = +1;
	}

	// Depending on the conic type, determine the remaining elements

	// ELLIPSE :
	if(ialpha == -1){
		a = -0.5*gmsum/energy;
	  	fac = 1. - h2/(gmsum*a);
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
	  	a = 0.5*gmsum/energy;
	  	fac = h2/(gmsum*a);
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
	    		double tmpf = 0.5*h2/gmsum;
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
		a =  0.5*h2/gmsum;
	  	e = 1.;
	  	w = acos(2.*a/r -1.);
	  	if ( vdotr < 0.) w = 2.*M_PI - w;
	  	double tmpf = tan(0.5 * w);
	  	M = tmpf* (1. + tmpf*tmpf/3.);
	  	omega = u - w;
	  	if(omega < 0.) omega = omega + 2.*M_PI;
	  	omega = omega - floor(omega/(2.*M_PI))*2.*M_PI; 	 
	}
	
	o->Omega 	= Omega;
	o->omega 	= omega;
	o->M		= M;
	o->f		= f;
	o->E		= E;
	o->inc		= inc;
	o->e		= e;
	o->a		= a;
}
