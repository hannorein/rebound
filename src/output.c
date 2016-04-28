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
#include <string.h>
#include <sys/time.h>
#include "particle.h"
#include "rebound.h"
#include "tools.h"
#include "output.h"
#include "integrator_sei.h"
#include "input.h"
#ifdef MPI
#include "communication_mpi.h"
#include "mpi.h"
#endif // MPI

/**
 * @brief Same as reb_output_check but with a phase argument
 */
int reb_output_check_phase(struct reb_simulation* r, double interval,double phase){
	double shift = r->t+interval*phase;
	if (floor(shift/interval)!=floor((shift-r->dt)/interval)){
		return 1;
	}
	// Output at beginning 
	if (r->t==0){
		return 1;
	}
	return 0;
}

int reb_output_check(struct reb_simulation* r, double interval){
	return reb_output_check_phase(r, interval,0);
}


#ifdef PROFILING
#warning PROFILING enabled. Rebound is NOT thread-safe.
double profiling_time_sum[PROFILING_CAT_NUM];
double profiling_time_initial 	= 0;
double profiling_timing_initial	= 0;
double profiling_time_final 	= 0;
void profiling_start(void){
	struct timeval tim;
	gettimeofday(&tim, NULL);
	profiling_time_initial = tim.tv_sec+(tim.tv_usec/1000000.0);
}
void profiling_stop(int cat){
	struct timeval tim;
	gettimeofday(&tim, NULL);
	profiling_time_final = tim.tv_sec+(tim.tv_usec/1000000.0);
	profiling_time_sum[cat] += profiling_time_final - profiling_time_initial;
}
#endif // PROFILING

void reb_output_timing(struct reb_simulation* r, const double tmax){
	const int N = r->N;
#ifdef MPI
	int N_tot = 0;
	MPI_Reduce(&N, &N_tot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
	if (r->mpi_id!=0) return;
#else
	int N_tot = N;
#endif
	struct timeval tim;
	gettimeofday(&tim, NULL);
	double temp = tim.tv_sec+(tim.tv_usec/1000000.0);
	if (r->output_timing_last==-1){
		r->output_timing_last = temp;
	}else{
		printf("\r");
#ifdef PROFILING
		fputs("\033[A\033[2K",stdout);
		for (int i=0;i<=PROFILING_CAT_NUM;i++){
			fputs("\033[A\033[2K",stdout);
		}
#endif // PROFILING
	}
	printf("N_tot= %- 9d  ",N_tot);
	if (r->integrator==REB_INTEGRATOR_SEI){
		printf("t= %- 9f [orb]  ",r->t*r->ri_sei.OMEGA/2./M_PI);
	}else{
		printf("t= %- 9f  ",r->t);
	}
	printf("dt= %- 9f  ",r->dt);
	if (r->integrator==REB_INTEGRATOR_HYBRID){
		printf("INT= %- 1d  ",r->ri_hybrid.mode);
	}
	printf("cpu= %- 9f [s]  ",temp-r->output_timing_last);
	if (tmax>0){
		printf("t/tmax= %5.2f%%",r->t/tmax*100.0);
	}
#ifdef PROFILING
	if (profiling_timing_initial==0){
		struct timeval tim;
		gettimeofday(&tim, NULL);
		profiling_timing_initial = tim.tv_sec+(tim.tv_usec/1000000.0);
	}
	printf("\nCATEGORY       TIME \n");
	double _sum = 0;
	for (int i=0;i<=PROFILING_CAT_NUM;i++){
		switch (i){
			case PROFILING_CAT_INTEGRATOR:
				printf("Integrator     ");
				break;
			case PROFILING_CAT_BOUNDARY:
				printf("Boundary check ");
				break;
			case PROFILING_CAT_GRAVITY:
				printf("Gravity/Forces ");
				break;
			case PROFILING_CAT_COLLISION:
				printf("Collisions     ");
				break;
#ifdef OPENGL
			case PROFILING_CAT_VISUALIZATION:
				printf("Visualization  ");
				break;
#endif // OPENGL
			case PROFILING_CAT_NUM:
				printf("Other          ");
				break;
		}
		if (i==PROFILING_CAT_NUM){
			printf("%5.2f%%",(1.-_sum/(profiling_time_final - profiling_timing_initial))*100.);
		}else{
			printf("%5.2f%%\n",profiling_time_sum[i]/(profiling_time_final - profiling_timing_initial)*100.);
			_sum += profiling_time_sum[i];
		}
	}
#endif // PROFILING
	fflush(stdout);
	r->output_timing_last = temp;
}


void reb_output_ascii(struct reb_simulation* r, char* filename){
	const int N = r->N;
#ifdef MPI
	char filename_mpi[1024];
	sprintf(filename_mpi,"%s_%d",filename,r->mpi_id);
	FILE* of = fopen(filename_mpi,"a"); 
#else // MPI
	FILE* of = fopen(filename,"a"); 
#endif // MPI
	if (of==NULL){
		reb_exit("Can not open file.");
	}
	for (int i=0;i<N;i++){
		struct reb_particle p = r->particles[i];
		fprintf(of,"%e\t%e\t%e\t%e\t%e\t%e\n",p.x,p.y,p.z,p.vx,p.vy,p.vz);
	}
	fclose(of);
}

void reb_output_orbits(struct reb_simulation* r, char* filename){
	const int N = r->N;
#ifdef MPI
	char filename_mpi[1024];
	sprintf(filename_mpi,"%s_%d",filename,r->mpi_id);
	FILE* of = fopen(filename_mpi,"a"); 
#else // MPI
	FILE* of = fopen(filename,"a"); 
#endif // MPI
	if (of==NULL){
		reb_exit("Can not open file.");
	}
	struct reb_particle com = r->particles[0];
	for (int i=1;i<N;i++){
		struct reb_orbit o = reb_tools_particle_to_orbit(r->G, r->particles[i],com);
		fprintf(of,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",r->t,o.a,o.e,o.inc,o.Omega,o.omega,o.l,o.P,o.f);
		com = reb_get_com_of_pair(com,r->particles[i]);
	}
	fclose(of);
}

void reb_output_binary(struct reb_simulation* r, char* filename){
#ifdef MPI
	char filename_mpi[1024];
	sprintf(filename_mpi,"%s_%d",filename,r->mpi_id);
	FILE* of = fopen(filename_mpi,"wb"); 
#else // MPI
	FILE* of = fopen(filename,"wb"); 
#endif // MPI
	if (of==NULL){
		reb_exit("Can not open file.");
	}
	fwrite(r,sizeof(struct reb_simulation),1,of);
	fwrite(r->particles,sizeof(struct reb_particle),r->N,of);
    if (r->var_config_N){
	    fwrite(r->var_config,sizeof(struct reb_variational_configuration),r->var_config_N,of);
    }
	fclose(of);
}

void reb_output_binary_positions(struct reb_simulation* r, char* filename){
	const int N = r->N;
#ifdef MPI
	char filename_mpi[1024];
	sprintf(filename_mpi,"%s_%d",filename,r->mpi_id);
	FILE* of = fopen(filename_mpi,"wb"); 
#else // MPI
	FILE* of = fopen(filename,"wb"); 
#endif // MPI
	if (of==NULL){
		reb_exit("Can not open file.");
	}
	for (int i=0;i<N;i++){
		struct reb_vec3d v;
		v.x = r->particles[i].x;
		v.y = r->particles[i].y;
		v.z = r->particles[i].z;
		fwrite(&(v),sizeof(struct reb_vec3d),1,of);
	}
	fclose(of);
}

void reb_output_velocity_dispersion(struct reb_simulation* r, char* filename){
	const int N = r->N;
	// Algorithm with reduced roundoff errors (see wikipedia)
	struct reb_vec3d A = {.x=0, .y=0, .z=0};
	struct reb_vec3d Q = {.x=0, .y=0, .z=0};
	for (int i=0;i<N;i++){
		struct reb_vec3d Aim1 = A;
		struct reb_particle p = r->particles[i];
		A.x = A.x + (p.vx-A.x)/(double)(i+1);
		if (r->integrator==REB_INTEGRATOR_SEI){
			A.y = A.y + (p.vy+1.5*r->ri_sei.OMEGA*p.x-A.y)/(double)(i+1);
		}else{
			A.y = A.y + (p.vy-A.y)/(double)(i+1);
		}
		A.z = A.z + (p.vz-A.z)/(double)(i+1);
		Q.x = Q.x + (p.vx-Aim1.x)*(p.vx-A.x);
		if (r->integrator==REB_INTEGRATOR_SEI){
			Q.y = Q.y + (p.vy+1.5*r->ri_sei.OMEGA*p.x-Aim1.y)*(p.vy+1.5*r->ri_sei.OMEGA*p.x-A.y);
		}else{
			Q.y = Q.y + (p.vy-Aim1.y)*(p.vy-A.y);
		}
		Q.z = Q.z + (p.vz-Aim1.z)*(p.vz-A.z);
	}
#ifdef MPI
	int N_tot = 0;
	struct reb_vec3d A_tot = {.x=0, .y=0, .z=0};
	struct reb_vec3d Q_tot = {.x=0, .y=0, .z=0};
	MPI_Reduce(&N, &N_tot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
	MPI_Reduce(&A, &A_tot, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
	MPI_Reduce(&Q, &Q_tot, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
	if (r->mpi_id!=0) return;
#else
	int N_tot = N;
	struct reb_vec3d A_tot = A;
	struct reb_vec3d Q_tot = Q;
#endif
	Q_tot.x = sqrt(Q_tot.x/(double)N_tot);
	Q_tot.y = sqrt(Q_tot.y/(double)N_tot);
	Q_tot.z = sqrt(Q_tot.z/(double)N_tot);
	FILE* of = fopen(filename,"a"); 
	if (of==NULL){
		reb_exit("Can not open file.");
	}
	fprintf(of,"%e\t%e\t%e\t%e\t%e\t%e\t%e\n",r->t,A_tot.x,A_tot.y,A_tot.z,Q_tot.x,Q_tot.y,Q_tot.z);
	fclose(of);
}

	
