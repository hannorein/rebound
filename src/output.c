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
#include <stddef.h>
#include "particle.h"
#include "rebound.h"
#include "tools.h"
#include "output.h"
#include "integrator.h"
#include "integrator_sei.h"
#include "integrator_tes.h"

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
        reb_error(r, "Can not open file.");
        return;
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
        reb_error(r, "Can not open file.");
        return;
    }
    struct reb_particle com = r->particles[0];
    for (int i=1;i<N;i++){
        struct reb_orbit o = reb_tools_particle_to_orbit(r->G, r->particles[i],com);
        fprintf(of,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",r->t,o.a,o.e,o.inc,o.Omega,o.omega,o.l,o.P,o.f);
        com = reb_get_com_of_pair(com,r->particles[i]);
    }
    fclose(of);
}

static inline void reb_save_controlVars(controlVars* dp7, char** bufp, size_t* sizep, size_t* allocatedsize){
    reb_output_stream_write(bufp, allocatedsize, sizep, &dp7->size, sizeof(uint32_t));
    reb_output_stream_write(bufp, allocatedsize, sizep, dp7->p0, dp7->size);
    reb_output_stream_write(bufp, allocatedsize, sizep, dp7->p1, dp7->size);
    reb_output_stream_write(bufp, allocatedsize, sizep, dp7->p2, dp7->size);
    reb_output_stream_write(bufp, allocatedsize, sizep, dp7->p3, dp7->size);
    reb_output_stream_write(bufp, allocatedsize, sizep, dp7->p4, dp7->size);
    reb_output_stream_write(bufp, allocatedsize, sizep, dp7->p5, dp7->size);
    reb_output_stream_write(bufp, allocatedsize, sizep, dp7->p6, dp7->size);
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


    // Compress data if possible
    // This does not affect future calculation, but might trigger a realloc.
    if (r->ri_ias15.allocatedN > 3*r->N){
        r->ri_ias15.allocatedN = 3*r->N;
    }

    /// Output all fields
    int i=0;
    while (reb_binary_field_descriptor_list[i].type!=9999){
        int dtype = reb_binary_field_descriptor_list[i].dtype;
        // Simple data types:
        if (dtype == REB_DOUBLE || dtype == REB_INT || dtype == REB_UINT || dtype == REB_UINT32 ||
                dtype == REB_LONG || dtype == REB_ULONG || dtype == REB_ULONGLONG ||
                dtype == REB_PARTICLE || dtype == REB_VEC3D ){
            struct reb_binary_field field;
            memset(&field,0,sizeof(struct reb_binary_field));
            field.type = reb_binary_field_descriptor_list[i].type;
            switch (dtype){
                case REB_DOUBLE: 
                    field.size = sizeof(double);
                    break;
                case REB_INT: 
                    field.size = sizeof(int);
                    break;
                case REB_UINT: 
                    field.size = sizeof(unsigned int);
                    break;
                case REB_UINT32: 
                    field.size = sizeof(uint32_t);
                    break;
                case REB_LONG:
                    field.size = sizeof(long);
                    break;
                case REB_ULONG:
                    field.size = sizeof(unsigned long);
                    break;
                case REB_ULONGLONG:
                    field.size = sizeof(unsigned long long);
                    break;
                case REB_VEC3D:
                    field.size = sizeof(struct reb_vec3d);
                    break;
                case REB_PARTICLE:
                    field.size = sizeof(struct reb_particle);
                    break;
            }
            reb_output_stream_write(bufp, &allocatedsize, sizep, &field, sizeof(struct reb_binary_field));
            char* pointer = (char*)r + reb_binary_field_descriptor_list[i].offset;
            reb_output_stream_write(bufp, &allocatedsize, sizep, pointer, field.size);
        }
        // Pointer data types
        if (dtype == REB_POINTER || dtype == REB_POINTER_ALIGNED ){
            struct reb_binary_field field;
            memset(&field,0,sizeof(struct reb_binary_field));
            field.type = reb_binary_field_descriptor_list[i].type;
            unsigned int* pointer_N = (unsigned int*)((char*)r + reb_binary_field_descriptor_list[i].offset_N);
            field.size = (*pointer_N) * reb_binary_field_descriptor_list[i].element_size;
                
            if (field.size){
                reb_output_stream_write(bufp, &allocatedsize, sizep, &field, sizeof(struct reb_binary_field));
                char* pointer = (char*)r + reb_binary_field_descriptor_list[i].offset;
                pointer = *(char**)pointer;
                reb_output_stream_write(bufp, &allocatedsize, sizep, pointer, field.size);
            }
        }
        // Special datatype for IAS15. Similar to POINTER
        if (dtype == REB_DP7 ){
            struct reb_binary_field field;
            memset(&field,0,sizeof(struct reb_binary_field));
            field.type = reb_binary_field_descriptor_list[i].type;
            unsigned int* pointer_N = (unsigned int*)((char*)r + reb_binary_field_descriptor_list[i].offset_N);
            field.size = (*pointer_N) * reb_binary_field_descriptor_list[i].element_size;
                
            if (field.size){
                reb_output_stream_write(bufp, &allocatedsize, sizep, &field, sizeof(struct reb_binary_field));
                char* pointer = (char*)r + reb_binary_field_descriptor_list[i].offset;
                struct reb_dp7* dp7 = (struct reb_dp7*)pointer;
                reb_output_stream_write(bufp, &allocatedsize, sizep, dp7->p0,field.size/7);
                reb_output_stream_write(bufp, &allocatedsize, sizep, dp7->p1,field.size/7);
                reb_output_stream_write(bufp, &allocatedsize, sizep, dp7->p2,field.size/7);
                reb_output_stream_write(bufp, &allocatedsize, sizep, dp7->p3,field.size/7);
                reb_output_stream_write(bufp, &allocatedsize, sizep, dp7->p4,field.size/7);
                reb_output_stream_write(bufp, &allocatedsize, sizep, dp7->p5,field.size/7);
                reb_output_stream_write(bufp, &allocatedsize, sizep, dp7->p6,field.size/7);
            }
        }
        i++;
    }

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


    // Output fields for TES integrator.
    WRITE_FIELD(TES_COM,                &r->ri_tes.COM,                   3*sizeof(double));
    WRITE_FIELD(TES_COM_DOT,            &r->ri_tes.COM_dot,               3*sizeof(double));
    WRITE_FIELD(TES_ALLOCATED_N,        &r->ri_tes.allocated_N,           sizeof(uint32_t));

    
    if(r->ri_tes.allocated_N)
    {
        uint32_t N = r->ri_tes.allocated_N;
        WRITE_FIELD(TES_PARTICLES_DH,       r->ri_tes.particles_dh,          N*sizeof(struct reb_particle));
        WRITE_FIELD(TES_MASS,               r->ri_tes.mass,                  N*sizeof(double));
        WRITE_FIELD(TES_X_DH,               r->ri_tes.X_dh,                  6*N*sizeof(double));

        // Kepler solver vars.
        WRITE_FIELD(TES_UVARS_SV_SIZE, &r->ri_tes.uVars->stateVectorSize, sizeof(uint32_t));
        WRITE_FIELD(TES_UVARS_T0, r->ri_tes.uVars->t0, r->ri_tes.uVars->controlVectorSize);
        WRITE_FIELD(TES_UVARS_TLAST, r->ri_tes.uVars->tLast, r->ri_tes.uVars->controlVectorSize);
        WRITE_FIELD(TES_UVARS_CSQ, r->ri_tes.uVars->uv_csq, r->ri_tes.uVars->stateVectorSize);
        WRITE_FIELD(TES_UVARS_CSP, r->ri_tes.uVars->uv_csp, r->ri_tes.uVars->stateVectorSize);
        WRITE_FIELD(TES_UVARS_CSV, r->ri_tes.uVars->uv_csv, r->ri_tes.uVars->stateVectorSize);
        WRITE_FIELD(TES_UVARS_Q0, r->ri_tes.uVars->Q0, r->ri_tes.uVars->stateVectorSize);
        WRITE_FIELD(TES_UVARS_V0, r->ri_tes.uVars->V0, r->ri_tes.uVars->stateVectorSize);
        WRITE_FIELD(TES_UVARS_P0, r->ri_tes.uVars->P0, r->ri_tes.uVars->stateVectorSize);
        WRITE_FIELD(TES_UVARS_Q1, r->ri_tes.uVars->Q1, r->ri_tes.uVars->stateVectorSize);
        WRITE_FIELD(TES_UVARS_V1, r->ri_tes.uVars->V1, r->ri_tes.uVars->stateVectorSize);
        WRITE_FIELD(TES_UVARS_P1, r->ri_tes.uVars->P1, r->ri_tes.uVars->stateVectorSize);
        WRITE_FIELD(TES_UVARS_X, r->ri_tes.uVars->X, r->ri_tes.uVars->controlVectorSize);
        WRITE_FIELD(TES_UVARS_Q0_NORM, r->ri_tes.uVars->Q0_norm, r->ri_tes.uVars->controlVectorSize);
        WRITE_FIELD(TES_UVARS_BETA, r->ri_tes.uVars->beta, r->ri_tes.uVars->controlVectorSize);
        WRITE_FIELD(TES_UVARS_ETA, r->ri_tes.uVars->eta, r->ri_tes.uVars->controlVectorSize);
        WRITE_FIELD(TES_UVARS_ZETA, r->ri_tes.uVars->zeta, r->ri_tes.uVars->controlVectorSize);
        WRITE_FIELD(TES_UVARS_PERIOD, r->ri_tes.uVars->period, r->ri_tes.uVars->controlVectorSize);
        WRITE_FIELD(TES_UVARS_XPERIOD, r->ri_tes.uVars->Xperiod, r->ri_tes.uVars->controlVectorSize);
        WRITE_FIELD(TES_UVARS_STUMPF_C0, r->ri_tes.uVars->C.c0, r->ri_tes.uVars->controlVectorSize);
        WRITE_FIELD(TES_UVARS_STUMPF_C1, r->ri_tes.uVars->C.c1, r->ri_tes.uVars->controlVectorSize);
        WRITE_FIELD(TES_UVARS_STUMPF_C2, r->ri_tes.uVars->C.c2, r->ri_tes.uVars->controlVectorSize);
        WRITE_FIELD(TES_UVARS_STUMPF_C3, r->ri_tes.uVars->C.c3, r->ri_tes.uVars->controlVectorSize);
        WRITE_FIELD(TES_UVARS_MU, &r->ri_tes.uVars->mu, sizeof(double));

        // Integrator vars
        WRITE_FIELD(TES_RADAU_DX, r->ri_tes.radau->dX, r->ri_tes.stateVectorSize);
        WRITE_FIELD(TES_RADAU_XOUT, r->ri_tes.radau->Xout, r->ri_tes.stateVectorSize);
        WRITE_FIELD(TES_RADAU_RECTI_ARRAY, r->ri_tes.radau->rectifiedArray, sizeof(uint32_t)*r->ri_tes.stateVectorLength);
        WRITE_FIELD(TES_RADAU_PREDICTORS, r->ri_tes.radau->predictors, r->ri_tes.stateVectorSize);
        WRITE_FIELD(TES_RADAU_DSTATE0, r->ri_tes.radau->dState0, r->ri_tes.stateVectorSize);
        WRITE_FIELD(TES_RADAU_DDSTATE0, r->ri_tes.radau->ddState0, r->ri_tes.stateVectorSize);
        WRITE_FIELD(TES_RADAU_DSTATE, r->ri_tes.radau->dState, r->ri_tes.stateVectorSize);
        WRITE_FIELD(TES_RADAU_DDSTATE, r->ri_tes.radau->ddState, r->ri_tes.stateVectorSize);
        WRITE_FIELD(TES_RADAU_CS_DSTATE0, r->ri_tes.radau->cs_dState0, r->ri_tes.stateVectorSize);
        WRITE_FIELD(TES_RADAU_CS_DDSTATE0, r->ri_tes.radau->cs_ddState0, r->ri_tes.stateVectorSize);
        WRITE_FIELD(TES_RADAU_CS_DSTATE, r->ri_tes.radau->cs_dState, r->ri_tes.stateVectorSize);
        WRITE_FIELD(TES_RADAU_CS_DDSTATE, r->ri_tes.radau->cs_ddState, r->ri_tes.stateVectorSize);
        WRITE_FIELD(TES_RADAU_CS_DX, r->ri_tes.radau->cs_dX, r->ri_tes.stateVectorSize);
        WRITE_FIELD(TES_RADAU_FCALLS, &r->ri_tes.radau->fCalls, sizeof(uint64_t));
        WRITE_FIELD(TES_RADAU_RECTIS, &r->ri_tes.radau->rectifications, sizeof(uint64_t));
        WRITE_FIELD(TES_RADAU_ITERS, &r->ri_tes.radau->convergenceIterations, sizeof(uint32_t));
        WRITE_FIELD(TES_RADAU_B6, r->ri_tes.radau->b6_store, r->ri_tes.stateVectorSize);

        {
            uint32_t array_size = r->ri_tes.radau->B.size;
            struct reb_binary_field field = {.type = REB_BINARY_FIELD_TYPE_TES_RADAU_B, .size = 7*array_size+sizeof(uint32_t)};
            reb_output_stream_write(bufp, &allocatedsize, sizep, &field,sizeof(struct reb_binary_field));
            reb_save_controlVars(&r->ri_tes.radau->B, bufp,sizep,&allocatedsize);
        }
        {
            uint32_t array_size = r->ri_tes.radau->Blast.size;
            struct reb_binary_field field = {.type = REB_BINARY_FIELD_TYPE_TES_RADAU_BLAST, .size = 7*array_size+sizeof(uint32_t)};
            reb_output_stream_write(bufp, &allocatedsize, sizep, &field,sizeof(struct reb_binary_field));
            reb_save_controlVars(&r->ri_tes.radau->Blast,bufp,sizep,&allocatedsize);
        }
        {
            uint32_t array_size = r->ri_tes.radau->B_1st.size;
            struct reb_binary_field field = {.type = REB_BINARY_FIELD_TYPE_TES_RADAU_B_1ST, .size = 7*array_size+sizeof(uint32_t)};
            reb_output_stream_write(bufp, &allocatedsize, sizep, &field,sizeof(struct reb_binary_field));
            reb_save_controlVars(&r->ri_tes.radau->B_1st,bufp,sizep,&allocatedsize);
        }
        {
            uint32_t array_size = r->ri_tes.radau->Blast_1st.size;
            struct reb_binary_field field = {.type = REB_BINARY_FIELD_TYPE_TES_RADAU_BLAST_1ST, .size = 7*array_size+sizeof(uint32_t)};
            reb_output_stream_write(bufp, &allocatedsize, sizep, &field,sizeof(struct reb_binary_field));
            reb_save_controlVars(&r->ri_tes.radau->Blast_1st,bufp,sizep,&allocatedsize);
        }
        {
            uint32_t array_size = r->ri_tes.radau->cs_B.size;
            struct reb_binary_field field = {.type = REB_BINARY_FIELD_TYPE_TES_RADAU_CS_B, .size = 7*array_size+sizeof(uint32_t)};
            reb_output_stream_write(bufp, &allocatedsize, sizep, &field,sizeof(struct reb_binary_field));
            reb_save_controlVars(&r->ri_tes.radau->cs_B,bufp,sizep,&allocatedsize);
        }
        {
            uint32_t array_size = r->ri_tes.radau->cs_B1st.size;
            struct reb_binary_field field = {.type = REB_BINARY_FIELD_TYPE_TES_RADAU_CS_B_1ST, .size = 7*array_size+sizeof(uint32_t)};
            reb_output_stream_write(bufp, &allocatedsize, sizep, &field,sizeof(struct reb_binary_field));
            reb_save_controlVars(&r->ri_tes.radau->cs_B1st,bufp,sizep,&allocatedsize);
        }          
        {
            uint32_t array_size = r->ri_tes.radau->G.size;
            struct reb_binary_field field = {.type = REB_BINARY_FIELD_TYPE_TES_RADAU_G, .size = 7*array_size+sizeof(uint32_t)};
            reb_output_stream_write(bufp, &allocatedsize, sizep, &field,sizeof(struct reb_binary_field));
            reb_save_controlVars(&r->ri_tes.radau->G,bufp,sizep,&allocatedsize);
        }
        {
            uint32_t array_size = r->ri_tes.radau->G_1st.size;
            struct reb_binary_field field = {.type = REB_BINARY_FIELD_TYPE_TES_RADAU_G_1ST, .size = 7*array_size+sizeof(uint32_t)};
            reb_output_stream_write(bufp, &allocatedsize, sizep, &field,sizeof(struct reb_binary_field));
            reb_save_controlVars(&r->ri_tes.radau->G_1st,bufp,sizep,&allocatedsize);
        }  
        // force model vars
        WRITE_FIELD(TES_DHEM_XOSC_STORE, r->ri_tes.rhs->XoscStore, 9*r->ri_tes.stateVectorSize);
        WRITE_FIELD(TES_DHEM_XOSC_PRED_STORE, r->ri_tes.rhs->XoscPredStore, 9*r->ri_tes.stateVectorSize);
        WRITE_FIELD(TES_DHEM_XOSC_CS_STORE, r->ri_tes.rhs->XoscStore_cs, 9*r->ri_tes.stateVectorSize);
        WRITE_FIELD(TES_DHEM_XOSC_DOT_STORE, r->ri_tes.rhs->Xosc_dotStore, 9*r->ri_tes.stateVectorSize);        
        WRITE_FIELD(TES_DHEM_X, r->ri_tes.rhs->X, r->ri_tes.stateVectorSize);
        WRITE_FIELD(TES_DHEM_M_INV, r->ri_tes.rhs->m_inv, r->ri_tes.controlVectorSize);
        WRITE_FIELD(TES_DHEM_M_TOTAL, &r->ri_tes.rhs->mTotal, sizeof(double));
        WRITE_FIELD(TES_DHEM_RECTI_TIME, r->ri_tes.rhs->rectifyTimeArray, r->ri_tes.controlVectorSize);
        WRITE_FIELD(TES_DHEM_RECTI_PERIOD, r->ri_tes.rhs->rectificationPeriod, r->ri_tes.controlVectorSize);
    
    }
   
    // To output size of binary file, need to calculate it first. 
    if (r->simulationarchive_version<3){ // to be removed in a future release
        r->simulationarchive_size_first = (*sizep)+sizeof(struct reb_binary_field)*2+sizeof(long)+sizeof(struct reb_simulationarchive_blob16);
    }else{
        r->simulationarchive_size_first = (*sizep)+sizeof(struct reb_binary_field)*2+sizeof(long)+sizeof(struct reb_simulationarchive_blob);
    }
    WRITE_FIELD(SASIZEFIRST,        &r->simulationarchive_size_first,   sizeof(long));
    int end_null = 0;
    WRITE_FIELD(END, &end_null, 0);
    if (r->simulationarchive_version<3){ // to be removed in a future release
        struct reb_simulationarchive_blob16 blob = {0};
        reb_output_stream_write(bufp, &allocatedsize, sizep, &blob, sizeof(struct reb_simulationarchive_blob16));
    }else{
        struct reb_simulationarchive_blob blob = {0};
        reb_output_stream_write(bufp, &allocatedsize, sizep, &blob, sizeof(struct reb_simulationarchive_blob));
    }
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
        reb_error(r, "Can not open file.");
        return;
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
        reb_error(r, "Can not open file.");
        return;
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
        reb_error(r, "Can not open file.");
        return;
    }
    fprintf(of,"%e\t%e\t%e\t%e\t%e\t%e\t%e\n",r->t,A_tot.x,A_tot.y,A_tot.z,Q_tot.x,Q_tot.y,Q_tot.z);
    fclose(of);
}

    
