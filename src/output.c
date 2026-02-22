/**
 * @file    output.c
 * @brief   ASCII output and profiling routines.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section     LICENSE
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
#include "rebound.h"
#include "rebound_internal.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "particle.h"
#include "tools.h"
#include "output.h"
#include "server.h"

#ifdef MPI
#include "communication_mpi.h"
#include "mpi.h"
#endif // MPI



// Same as reb_simulation_output_check but with a phase argument
int reb_simulation_output_check_phase(struct reb_simulation* r, double interval,double phase){
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

int reb_simulation_output_check(struct reb_simulation* r, double interval){
    return reb_simulation_output_check_phase(r, interval,0);
}

#ifdef PROFILING
#warning PROFILING enabled. Rebound is NOT thread-safe.
double profiling_time_sum[PROFILING_CAT_NUM];
double profiling_time_initial   = 0;
double profiling_timing_initial = 0;
double profiling_time_final     = 0;
void profiling_start(void){
    struct reb_timeval tim;
    gettimeofday(&tim, NULL);
    profiling_time_initial = tim.tv_sec+(tim.tv_usec/1000000.0);
}
void profiling_stop(int cat){
    struct reb_timeval tim;
    gettimeofday(&tim, NULL);
    profiling_time_final = tim.tv_sec+(tim.tv_usec/1000000.0);
    profiling_time_sum[cat] += profiling_time_final - profiling_time_initial;
}
#endif // PROFILING

#ifdef __EMSCRIPTEN__
// fflush does not work in emscripten. Workaround.
EM_JS(void, reb_remove_last_line, (), {
        var output = document.getElementById("output");
        if (output){
        const lastIndex1 = output.value.lastIndexOf("\n");
        const lastIndex2 = output.value.lastIndexOf("\n",lastIndex1-1);
        const lastIndexNtot = output.value.lastIndexOf("N_tot=");
        if(lastIndex1>0 && lastIndex2<lastIndexNtot){
        output.value = output.value.substring(0, lastIndex2+1);
        }
        }
        });
#endif

int reb_simulation_output_screenshot(struct reb_simulation* r, const char* filename){
#ifdef SERVER
    if (!r->server_data){
        reb_simulation_error(r, "To take a screenshot, call reb_simulation_start_server() and connect a web browser.");
        return 0;
    }

    r->server_data->status_before_screenshot = r->status;
    // Tell client to take screenshot
    r->status = REB_STATUS_SCREENSHOT;

    // Release mutex so client can pull simulation
    if (r->server_data->mutex_locked_by_integrate){
#ifdef _WIN32
        ReleaseMutex(r->server_data->mutex);
#else // _WIN32
        pthread_mutex_unlock(&(r->server_data->mutex));
#endif // _WIN32
    }

    // Wait until screenshot arrives
    while (!r->server_data->screenshot && r->status <0){
        usleep(100);
        if (reb_sigint > 1){
            r->status = REB_STATUS_SIGINT;
        }
    }

    // Lock mutex again before continuing
    if (r->server_data->mutex_locked_by_integrate){
#ifdef _WIN32
        WaitForSingleObject(r->server_data->mutex, INFINITE);
#else // _WIN32
        pthread_mutex_lock(&(r->server_data->mutex)); 
#endif // _WIN32
    }

    r->status = r->server_data->status_before_screenshot;

    if (r->server_data->screenshot){
        FILE* f = fopen(filename,"wb");
        if (!f){
            reb_simulation_error(r, "Error opening output file for screenshot.");
            free(r->server_data->screenshot);
            r->server_data->screenshot = 0;
            r->server_data->N_screenshot = 0;
            return 0;
        }else{
            fwrite(r->server_data->screenshot, r->server_data->N_screenshot, 1, f);
            fclose(f);
            free(r->server_data->screenshot);
            r->server_data->screenshot = 0;
            r->server_data->N_screenshot = 0;
            return 1;
        }
    }
#else //SERVER
    reb_simulation_error(r, "To take a screenshot compile with SERVER=1, call reb_simulation_start_server(), and connect with a web browser.");
#endif //SERVER
    return 0;
}


void reb_simulation_output_timing(struct reb_simulation* r, const double tmax){
    const int N = r->N;
#ifdef MPI
    int N_tot = 0;
    MPI_Reduce(&N, &N_tot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
    if (r->mpi_id!=0) return;
#else
    int N_tot = N;
#endif
    struct reb_timeval tim;
    gettimeofday(&tim, NULL);
    double temp = tim.tv_sec+(tim.tv_usec/1000000.0);
    if (r->output_timing_last==-1){
        r->output_timing_last = temp;
    }else{
#ifdef __EMSCRIPTEN__
        reb_remove_last_line();
#else
        printf("\r");
#endif
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
    printf("cpu= %- 9f [s]  ",temp-r->output_timing_last);
    if (tmax>0){
        printf("t/tmax= %5.2f%%",r->t/tmax*100.0);
    }
#ifdef PROFILING
    if (profiling_timing_initial==0){
        struct reb_timeval tim;
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
        }else if (i==PROFILING_CAT_INTEGRATOR){
            printf("%5.2f%%\n",(profiling_time_sum[PROFILING_CAT_INTEGRATOR]-profiling_time_sum[PROFILING_CAT_GRAVITY])/(profiling_time_final - profiling_timing_initial)*100.);
        }else{
            printf("%5.2f%%\n",profiling_time_sum[i]/(profiling_time_final - profiling_timing_initial)*100.);
        }
        if (i!=PROFILING_CAT_NUM && i!=PROFILING_CAT_GRAVITY){ // already in INTEGRATOR
            _sum += profiling_time_sum[i];
        }
    }
#endif // PROFILING
#ifdef __EMSCRIPTEN__
    printf("\n");
#else
    fflush(stdout);
#endif
    r->output_timing_last = temp;
}


void reb_simulation_output_ascii(struct reb_simulation* r, char* filename){
    const int N = r->N;
#ifdef MPI
    char filename_mpi[1024];
    sprintf(filename_mpi,"%s_%d",filename,r->mpi_id);
    FILE* of = fopen(filename_mpi,"ab"); 
#else // MPI
    FILE* of = fopen(filename,"ab"); 
#endif // MPI
    if (of==NULL){
        reb_simulation_error(r, "Can not open file.");
        return;
    }
    for (int i=0;i<N;i++){
        struct reb_particle p = r->particles[i];
        fprintf(of,"%e\t%e\t%e\t%e\t%e\t%e\n",p.x,p.y,p.z,p.vx,p.vy,p.vz);
    }
    fclose(of);
}

void reb_simulation_output_orbits(struct reb_simulation* r, char* filename){
    const int N = r->N;
#ifdef MPI
    char filename_mpi[1024];
    sprintf(filename_mpi,"%s_%d",filename,r->mpi_id);
    FILE* of = fopen(filename_mpi,"ab"); 
#else // MPI
    FILE* of = fopen(filename,"ab"); 
#endif // MPI
    if (of==NULL){
        reb_simulation_error(r, "Can not open file.");
        return;
    }
    struct reb_particle com = r->particles[0];
    for (int i=1;i<N;i++){
        struct reb_orbit o = reb_orbit_from_particle(r->G, r->particles[i],com);
        fprintf(of,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",r->t,o.a,o.e,o.inc,o.Omega,o.omega,o.l,o.P,o.f);
        com = reb_particle_com_of_pair(com,r->particles[i]);
    }
    fclose(of);
}


void reb_simulation_output_velocity_dispersion(struct reb_simulation* r, char* filename){
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
    FILE* of = fopen(filename,"ab"); 
    if (of==NULL){
        reb_simulation_error(r, "Can not open file.");
        return;
    }
    fprintf(of,"%e\t%e\t%e\t%e\t%e\t%e\t%e\n",r->t,A_tot.x,A_tot.y,A_tot.z,Q_tot.x,Q_tot.y,Q_tot.z);
    fclose(of);
}

