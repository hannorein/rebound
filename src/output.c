#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "main.h"
#include "output.h"

#ifdef INTEGRATOR_SEI 	// Shearing sheet
extern double OMEGA;
#endif

// Check if output is needed

int output_check(double interval){
	if (floor(t/interval)!=floor((t+dt)/interval)){
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
