/**
 * @file    output.c
 * @brief   Output routines.
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
#include "integrator.h"
#include "integrator_sei.h"
#include "input.h"
#ifdef MPI
#include "communication_mpi.h"
#include "mpi.h"
#endif // MPI

/** 
 * @brief Replacement for open_memstream
 */
void reb_output_stream_write(char** bufp, size_t* allocatedsize, size_t* sizep, void* restrict data, size_t size){
    // Increase size
    int increased = 0;
    while (*allocatedsize==0 || (*sizep)+size>(*allocatedsize)){
        increased = 1;
	    *allocatedsize = (*allocatedsize) ? (*allocatedsize) * 2 : 32;
    }
    if (increased){
        *bufp = realloc(*bufp,*allocatedsize);
    }
    // Copy data to buffer
    memcpy((*bufp)+(*sizep),data,size);
    *sizep += size;
}

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
double profiling_time_initial   = 0;
double profiling_timing_initial = 0;
double profiling_time_final     = 0;
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

void static inline reb_save_dp7(struct reb_dp7* dp7, const int N3, char** bufp, size_t* sizep, size_t* allocatedsize){
    reb_output_stream_write(bufp, allocatedsize, sizep, dp7->p0,sizeof(double)*N3);
    reb_output_stream_write(bufp, allocatedsize, sizep, dp7->p1,sizeof(double)*N3);
    reb_output_stream_write(bufp, allocatedsize, sizep, dp7->p2,sizeof(double)*N3);
    reb_output_stream_write(bufp, allocatedsize, sizep, dp7->p3,sizeof(double)*N3);
    reb_output_stream_write(bufp, allocatedsize, sizep, dp7->p4,sizeof(double)*N3);
    reb_output_stream_write(bufp, allocatedsize, sizep, dp7->p4,sizeof(double)*N3);
    reb_output_stream_write(bufp, allocatedsize, sizep, dp7->p6,sizeof(double)*N3);
}

// Macro to write a single field to a binary file.
// Memset forces padding to be set to 0 (not necessary but
// helps when comparing binary files)
#define WRITE_FIELD(typename, value, length) {\
        struct reb_binary_field field;\
        memset(&field,0,sizeof(struct reb_binary_field));\
        field.type = REB_BINARY_FIELD_TYPE_##typename;\
        field.size = (length);\
        reb_output_stream_write(bufp, &allocatedsize, sizep, &field,sizeof(struct reb_binary_field));\
        reb_output_stream_write(bufp, &allocatedsize, sizep, value,field.size);\
    }


void reb_output_binary_to_stream(struct reb_simulation* r, char** bufp, size_t* sizep){
    size_t allocatedsize = 0;
    *bufp = NULL;
    *sizep = 0;
    // Init integrators. This helps with bit-by-bit reproducibility.
    reb_integrator_init(r);

    // Output header.
    char header[64] = "\0";
    int cwritten = sprintf(header,"REBOUND Binary File. Version: %s",reb_version_str);
    snprintf(header+cwritten+1,64-cwritten-1,"%s",reb_githash_str);
    reb_output_stream_write(bufp, &allocatedsize, sizep, header,sizeof(char)*64);
   
    WRITE_FIELD(T,                  &r->t,                              sizeof(double));
    WRITE_FIELD(G,                  &r->G,                              sizeof(double));
    WRITE_FIELD(SOFTENING,          &r->softening,                      sizeof(double));
    WRITE_FIELD(DT,                 &r->dt,                             sizeof(double));
    WRITE_FIELD(N,                  &r->N,                              sizeof(int));
    WRITE_FIELD(NVAR,               &r->N_var,                          sizeof(int));
    WRITE_FIELD(VARCONFIGN,         &r->var_config_N,                   sizeof(int));
    WRITE_FIELD(NACTIVE,            &r->N_active,                       sizeof(int));
    WRITE_FIELD(TESTPARTICLETYPE,   &r->testparticle_type,              sizeof(int));
    WRITE_FIELD(HASHCTR,            &r->hash_ctr,                       sizeof(int));
    WRITE_FIELD(OPENINGANGLE2,      &r->opening_angle2,                 sizeof(double));
    WRITE_FIELD(STATUS,             &r->status,                         sizeof(int));
    WRITE_FIELD(EXACTFINISHTIME,    &r->exact_finish_time,              sizeof(int));
    WRITE_FIELD(FORCEISVELOCITYDEP, &r->force_is_velocity_dependent,    sizeof(unsigned int));
    WRITE_FIELD(GRAVITYIGNORETERMS, &r->gravity_ignore_terms,           sizeof(unsigned int));
    WRITE_FIELD(OUTPUTTIMINGLAST,   &r->output_timing_last,             sizeof(double));
    WRITE_FIELD(SAVEMESSAGES,       &r->save_messages,                  sizeof(int));
    WRITE_FIELD(EXITMAXDISTANCE,    &r->exit_max_distance,              sizeof(double));
    WRITE_FIELD(EXITMINDISTANCE,    &r->exit_min_distance,              sizeof(double));
    WRITE_FIELD(USLEEP,             &r->usleep,                         sizeof(double));
    WRITE_FIELD(TRACKENERGYOFFSET,  &r->track_energy_offset,            sizeof(int));
    WRITE_FIELD(ENERGYOFFSET,       &r->energy_offset,                  sizeof(double));
    WRITE_FIELD(BOXSIZE,            &r->boxsize,                        sizeof(struct reb_vec3d));
    WRITE_FIELD(BOXSIZEMAX,         &r->boxsize_max,                    sizeof(double));
    WRITE_FIELD(ROOTSIZE,           &r->root_size,                      sizeof(double));
    WRITE_FIELD(ROOTN,              &r->root_n,                         sizeof(int));
    WRITE_FIELD(ROOTNX,             &r->root_nx,                        sizeof(int));
    WRITE_FIELD(ROOTNY,             &r->root_ny,                        sizeof(int));
    WRITE_FIELD(ROOTNZ,             &r->root_nz,                        sizeof(int));
    WRITE_FIELD(NGHOSTX,            &r->nghostx,                        sizeof(int));
    WRITE_FIELD(NGHOSTY,            &r->nghosty,                        sizeof(int));
    WRITE_FIELD(NGHOSTZ,            &r->nghostz,                        sizeof(int));
    WRITE_FIELD(COLLISIONRESOLVEKEEPSORTED, &r->collision_resolve_keep_sorted, sizeof(int));
    WRITE_FIELD(MINIMUMCOLLISIONVELOCITY, &r->minimum_collision_velocity, sizeof(double));
    WRITE_FIELD(COLLISIONSPLOG,     &r->collisions_plog,                sizeof(double));
    WRITE_FIELD(MAXRADIUS,          &r->max_radius,                     2*sizeof(double));
    WRITE_FIELD(COLLISIONSNLOG,     &r->collisions_Nlog,                sizeof(long));
    WRITE_FIELD(CALCULATEMEGNO,     &r->calculate_megno,                sizeof(int));
    WRITE_FIELD(MEGNOYS,            &r->megno_Ys,                       sizeof(double));
    WRITE_FIELD(MEGNOYSS,           &r->megno_Yss,                      sizeof(double));
    WRITE_FIELD(MEGNOCOVYT,         &r->megno_cov_Yt,                   sizeof(double));
    WRITE_FIELD(MEGNOVART,          &r->megno_var_t,                    sizeof(double));
    WRITE_FIELD(MEGNOMEANT,         &r->megno_mean_t,                   sizeof(double));
    WRITE_FIELD(MEGNOMEANY,         &r->megno_mean_Y,                   sizeof(double));
    WRITE_FIELD(MEGNON,             &r->megno_n,                        sizeof(long));
    WRITE_FIELD(SAVERSION,          &r->simulationarchive_version,      sizeof(int));
    WRITE_FIELD(SASIZESNAPSHOT,     &r->simulationarchive_size_snapshot,sizeof(long));
    WRITE_FIELD(SAAUTOINTERVAL,     &r->simulationarchive_auto_interval, sizeof(double));
    WRITE_FIELD(SAAUTOWALLTIME,     &r->simulationarchive_auto_walltime, sizeof(double));
    WRITE_FIELD(SANEXT,             &r->simulationarchive_next,         sizeof(double));
    WRITE_FIELD(WALLTIME,           &r->walltime,                       sizeof(double));
    WRITE_FIELD(COLLISION,          &r->collision,                      sizeof(int));
    WRITE_FIELD(VISUALIZATION,      &r->visualization,                  sizeof(int));
    WRITE_FIELD(INTEGRATOR,         &r->integrator,                     sizeof(int));
    WRITE_FIELD(BOUNDARY,           &r->boundary,                       sizeof(int));
    WRITE_FIELD(GRAVITY,            &r->gravity,                        sizeof(int));
    WRITE_FIELD(SEI_OMEGA,          &r->ri_sei.OMEGA,                   sizeof(double));
    WRITE_FIELD(SEI_OMEGAZ,         &r->ri_sei.OMEGAZ,                  sizeof(double));
    WRITE_FIELD(SEI_LASTDT,         &r->ri_sei.lastdt,                  sizeof(double));
    WRITE_FIELD(SEI_SINDT,          &r->ri_sei.sindt,                   sizeof(double));
    WRITE_FIELD(SEI_TANDT,          &r->ri_sei.tandt,                   sizeof(double));
    WRITE_FIELD(SEI_SINDTZ,         &r->ri_sei.sindtz,                  sizeof(double));
    WRITE_FIELD(SEI_TANDTZ,         &r->ri_sei.tandtz,                  sizeof(double));
    WRITE_FIELD(WHFAST_CORRECTOR,   &r->ri_whfast.corrector,            sizeof(unsigned int));
    WRITE_FIELD(WHFAST_RECALCJAC,   &r->ri_whfast.recalculate_coordinates_this_timestep, sizeof(unsigned int));
    WRITE_FIELD(WHFAST_SAFEMODE,    &r->ri_whfast.safe_mode,            sizeof(unsigned int));
    WRITE_FIELD(WHFAST_KEEPUNSYNC,  &r->ri_whfast.keep_unsynchronized,  sizeof(unsigned int));
    WRITE_FIELD(WHFAST_ISSYNCHRON,  &r->ri_whfast.is_synchronized,      sizeof(unsigned int));
    WRITE_FIELD(WHFAST_TIMESTEPWARN,&r->ri_whfast.timestep_warning,     sizeof(unsigned int));
    WRITE_FIELD(WHFAST_PJ,          r->ri_whfast.p_jh,                  sizeof(struct reb_particle)*r->ri_whfast.allocated_N);
    WRITE_FIELD(WHFAST_COORDINATES, &r->ri_whfast.coordinates,          sizeof(int));
    WRITE_FIELD(IAS15_EPSILON,      &r->ri_ias15.epsilon,               sizeof(double));
    WRITE_FIELD(IAS15_MINDT,        &r->ri_ias15.min_dt,                sizeof(double));
    WRITE_FIELD(IAS15_EPSILONGLOBAL,&r->ri_ias15.epsilon_global,        sizeof(unsigned int));
    WRITE_FIELD(IAS15_ITERATIONSMAX,&r->ri_ias15.iterations_max_exceeded,sizeof(unsigned long));
    WRITE_FIELD(IAS15_ALLOCATEDN,   &r->ri_ias15.allocatedN,            sizeof(int));
    WRITE_FIELD(JANUS_SCALEPOS,     &r->ri_janus.scale_pos,             sizeof(double));
    WRITE_FIELD(JANUS_SCALEVEL,     &r->ri_janus.scale_vel,             sizeof(double));
    WRITE_FIELD(JANUS_ORDER,        &r->ri_janus.order,                 sizeof(unsigned int));
    WRITE_FIELD(JANUS_ALLOCATEDN,   &r->ri_janus.allocated_N,           sizeof(unsigned int));
    WRITE_FIELD(JANUS_RECALC,       &r->ri_janus.recalculate_integer_coordinates_this_timestep, sizeof(unsigned int));
    WRITE_FIELD(JANUS_PINT,         r->ri_janus.p_int,                  sizeof(struct reb_particle_int)*r->ri_janus.allocated_N);
    WRITE_FIELD(MERCURIUS_HILLFAC,  &r->ri_mercurius.hillfac,           sizeof(double));
    WRITE_FIELD(MERCURIUS_SAFEMODE, &r->ri_mercurius.safe_mode,         sizeof(unsigned int));
    WRITE_FIELD(MERCURIUS_ISSYNCHRON, &r->ri_mercurius.is_synchronized, sizeof(unsigned int));
    WRITE_FIELD(MERCURIUS_DCRIT,    r->ri_mercurius.dcrit,              sizeof(double)*r->ri_mercurius.dcrit_allocatedN);
    WRITE_FIELD(MERCURIUS_COMPOS,   &(r->ri_mercurius.com_pos),         sizeof(struct reb_vec3d));
    WRITE_FIELD(MERCURIUS_COMVEL,   &(r->ri_mercurius.com_vel),         sizeof(struct reb_vec3d));
    WRITE_FIELD(PYTHON_UNIT_L,      &r->python_unit_l,                  sizeof(uint32_t));
    WRITE_FIELD(PYTHON_UNIT_M,      &r->python_unit_m,                  sizeof(uint32_t));
    WRITE_FIELD(PYTHON_UNIT_T,      &r->python_unit_t,                  sizeof(uint32_t));
    WRITE_FIELD(STEPSDONE,          &r->steps_done,                     sizeof(unsigned long long));
    WRITE_FIELD(SAAUTOSTEP,         &r->simulationarchive_auto_step,    sizeof(unsigned long long));
    WRITE_FIELD(SANEXTSTEP,         &r->simulationarchive_next_step,    sizeof(unsigned long long));
    int functionpointersused = 0;
    if (r->coefficient_of_restitution ||
        r->collision_resolve ||
        r->additional_forces ||
        r->heartbeat ||
        r->post_timestep_modifications ||
        r->free_particle_ap){
        functionpointersused = 1;
    }
    WRITE_FIELD(FUNCTIONPOINTERS,   &functionpointersused,              sizeof(int));
    {
        struct reb_binary_field field;
        memset(&field,0,sizeof(struct reb_binary_field));
        field.type = REB_BINARY_FIELD_TYPE_PARTICLES;
        field.size = sizeof(struct reb_particle)*r->N;
        reb_output_stream_write(bufp, &allocatedsize, sizep, &field,sizeof(struct reb_binary_field));
        // output one particle at a time to sanitize pointers.
        for (int l=0;l<r->N;l++){
            struct reb_particle op = r->particles[l];
            op.c = NULL;
            op.ap = NULL;
            op.sim = NULL;
            reb_output_stream_write(bufp, &allocatedsize, sizep, &op,sizeof(struct reb_particle));
        }
    } 
    if (r->var_config){
        WRITE_FIELD(VARCONFIG,      r->var_config,                      sizeof(struct reb_variational_configuration)*r->var_config_N);
    }
    if (r->ri_ias15.allocatedN){
        int N3 = r->ri_ias15.allocatedN;
        WRITE_FIELD(IAS15_AT,   r->ri_ias15.at,     sizeof(double)*N3);
        WRITE_FIELD(IAS15_X0,   r->ri_ias15.x0,     sizeof(double)*N3);
        WRITE_FIELD(IAS15_V0,   r->ri_ias15.v0,     sizeof(double)*N3);
        WRITE_FIELD(IAS15_A0,   r->ri_ias15.a0,     sizeof(double)*N3);
        WRITE_FIELD(IAS15_CSX,  r->ri_ias15.csx,    sizeof(double)*N3);
        WRITE_FIELD(IAS15_CSV,  r->ri_ias15.csv,    sizeof(double)*N3);
        WRITE_FIELD(IAS15_CSA0, r->ri_ias15.csa0,   sizeof(double)*N3);
        {
            struct reb_binary_field field = {.type = REB_BINARY_FIELD_TYPE_IAS15_G, .size = sizeof(double)*N3*7};
            reb_output_stream_write(bufp, &allocatedsize, sizep, &field,sizeof(struct reb_binary_field));
            reb_save_dp7(&(r->ri_ias15.g),N3,bufp,sizep,&allocatedsize);
        }
        {
            struct reb_binary_field field = {.type = REB_BINARY_FIELD_TYPE_IAS15_B, .size = sizeof(double)*N3*7};
            reb_output_stream_write(bufp, &allocatedsize, sizep, &field,sizeof(struct reb_binary_field));
            reb_save_dp7(&(r->ri_ias15.b),N3,bufp,sizep,&allocatedsize);
        }
        {
            struct reb_binary_field field = {.type = REB_BINARY_FIELD_TYPE_IAS15_CSB, .size = sizeof(double)*N3*7};
            reb_output_stream_write(bufp, &allocatedsize, sizep, &field,sizeof(struct reb_binary_field));
            reb_save_dp7(&(r->ri_ias15.csb),N3,bufp,sizep,&allocatedsize);
        }
        {
            struct reb_binary_field field = {.type = REB_BINARY_FIELD_TYPE_IAS15_E, .size = sizeof(double)*N3*7};
            reb_output_stream_write(bufp, &allocatedsize, sizep, &field,sizeof(struct reb_binary_field));
            reb_save_dp7(&(r->ri_ias15.e),N3,bufp,sizep,&allocatedsize);
        }
        {
            struct reb_binary_field field = {.type = REB_BINARY_FIELD_TYPE_IAS15_BR, .size = sizeof(double)*N3*7};
            reb_output_stream_write(bufp, &allocatedsize, sizep, &field,sizeof(struct reb_binary_field));
            reb_save_dp7(&(r->ri_ias15.br),N3,bufp,sizep,&allocatedsize);
        }
        {
            struct reb_binary_field field = {.type = REB_BINARY_FIELD_TYPE_IAS15_ER, .size = sizeof(double)*N3*7};
            reb_output_stream_write(bufp, &allocatedsize, sizep, &field,sizeof(struct reb_binary_field));
            reb_save_dp7(&(r->ri_ias15.er),N3,bufp,sizep,&allocatedsize);
        }
    }
    // To output size of binary file, need to calculate it first. 
    r->simulationarchive_size_first = (*sizep)+sizeof(struct reb_binary_field)*2+sizeof(long)+sizeof(struct reb_simulationarchive_blob);
    WRITE_FIELD(SASIZEFIRST,        &r->simulationarchive_size_first,   sizeof(long));
    int end_null = 0;
    WRITE_FIELD(END, &end_null, 0);
    struct reb_simulationarchive_blob blob = {0};
    reb_output_stream_write(bufp, &allocatedsize, sizep, &blob, sizeof(struct reb_simulationarchive_blob));
}

void reb_output_binary(struct reb_simulation* r, const char* filename){
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
    char* bufp;
    size_t sizep;
    reb_output_binary_to_stream(r, &bufp,&sizep);
    fwrite(bufp,sizep,1,of);
    free(bufp);
    fclose(of);
}

void reb_output_binary_positions(struct reb_simulation* r, const char* filename){
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

    
