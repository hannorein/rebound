#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "main.h"
#include "output.h"
#ifdef OPENGL
#include "opengl.h"
#ifdef LIBPNG
#include <png.h>
#endif
#ifdef _APPLE
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif 
#endif 

#ifdef INTEGRATOR_SEI 	// Shearing sheet
extern double OMEGA;
#endif

// Check if output is needed

int output_check(double interval){
	if (floor(t/interval)!=floor((t-dt)/interval)){
		return 1;
	}
	// Output at beginning or end of simulation
	if (t==0||t==tmax){
		return 1;
	}
	return 0;
}


// Various output routines that are used often.

struct vec3 {
	double x;
	double y;
	double z;
} vec3;

void output_append_ascii(char* filename){
	FILE* of = fopen(filename,"a"); 
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
		fprintf(of,"%e\t%e\t%e\t%e\t%e\t%e\n",p.x,p.y,p.z,p.vx,p.vy,p.vz);
	}
	fclose(of);
}

void output_ascii(char* filename){
	FILE* of = fopen(filename,"w"); 
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
		fprintf(of,"%e\t%e\t%e\t%e\t%e\t%e\n",p.x,p.y,p.z,p.vx,p.vy,p.vz);
	}
	fclose(of);
}


void output_binary(char* filename){
	FILE* of = fopen(filename,"wb"); 
	fwrite(&N,sizeof(int),1,of);
	fwrite(&t,sizeof(double),1,of);
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
		fwrite(&(p),sizeof(struct particle),1,of);
	}
	fclose(of);
}

void output_binary_positions(char* filename){
	FILE* of = fopen(filename,"wb"); 
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
	Q.x = sqrt(Q.x/(double)N);
	Q.y = sqrt(Q.y/(double)N);
	Q.z = sqrt(Q.z/(double)N);
	FILE* of = fopen(filename,"a"); 
	fprintf(of,"%e\t%e\t%e\t%e\t%e\t%e\t%e\n",t,A.x,A.y,A.z,Q.x,Q.y,Q.z);
	fclose(of);
}

#ifdef OPENGL
#ifdef LIBPNG
unsigned char* 	imgdata = NULL;
void output_png(char* filename){
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
