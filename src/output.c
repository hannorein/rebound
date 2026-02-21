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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "particle.h"
#include "rebound.h"
#include "tools.h"
#include "tree.h"
#include "output.h"
#include "server.h"
#include "simulationarchive.h"
#include "integrator.h"
#include "integrator_sei.h"
#include "integrator_whfast512.h"

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


// Macro to read a single field from a binary file.
#define CASE(typename, value) case REB_BINARY_FIELD_TYPE_##typename: \
{\
    fread(value, field.size,1,inf);\
    goto next_field;\
}\
break;

#define CASE_CONTROL_VARS(typename, valueref) case REB_BINARY_FIELD_TYPE_##typename: \
{\
    fread(&valueref->size, sizeof(uint32_t),1,inf);\
    fread(valueref->p0, valueref->size,1,inf);\
    fread(valueref->p1, valueref->size,1,inf);\
    fread(valueref->p2, valueref->size,1,inf);\
    fread(valueref->p3, valueref->size,1,inf);\
    fread(valueref->p4, valueref->size,1,inf);\
    fread(valueref->p5, valueref->size,1,inf);\
    fread(valueref->p6, valueref->size,1,inf);\
    goto next_field;\
}\
break;        


void reb_input_fields(struct reb_simulation* r, FILE* inf, enum reb_simulation_binary_error_codes* warnings){
    struct reb_binary_field field;
    // A few fields need special treatment. Find their descriptors first.
    struct reb_binary_field_descriptor fd_header = reb_binary_field_descriptor_for_name("header");
    struct reb_binary_field_descriptor fd_end = reb_binary_field_descriptor_for_name("end");
    struct reb_binary_field_descriptor fd_functionpointers = reb_binary_field_descriptor_for_name("functionpointers");

next_field:
    // Loop over all fields
    while(1){

        int numread = (int)fread(&field,sizeof(struct reb_binary_field),1,inf);
        if (numread<1){
            goto finish_fields; // End of file
        }
        if (field.type==fd_end.type){
            goto finish_fields; // End of snapshot
        }
        int i=0;

        // Loop over field descriptor list. Simple datatypes and pointers will be read in this loop.
        while (reb_binary_field_descriptor_list[i].dtype!=REB_FIELD_END){
            struct reb_binary_field_descriptor fd = reb_binary_field_descriptor_list[i];
            if (fd.type==field.type){
                // Read simple data types
                if (fd.dtype == REB_DOUBLE || fd.dtype == REB_INT || fd.dtype == REB_UINT 
                        || fd.dtype == REB_UINT32 || fd.dtype == REB_INT64 
                        || fd.dtype == REB_UINT64 || fd.dtype == REB_PARTICLE 
                        || fd.dtype == REB_PARTICLE4 || fd.dtype == REB_VEC3D ){
                    char* pointer = (char*)r + reb_binary_field_descriptor_list[i].offset;
                    fread(pointer, field.size, 1, inf);
                    goto next_field;
                }
                // Read a pointer data type. 
                // 1) reallocate memory
                // 2) read data into memory
                // 3) set N_allocated variable
                if (fd.dtype == REB_POINTER || fd.dtype == REB_POINTER_ALIGNED){
                    if (field.size % reb_binary_field_descriptor_list[i].element_size){
                        reb_simulation_warning(r, "Inconsistent size encountered in binary field.");
                    }
                    char* pointer = (char*)r + reb_binary_field_descriptor_list[i].offset;
                    if (fd.dtype == REB_POINTER_ALIGNED){
                        if (*(char**)pointer) free(*(char**)pointer);
#if defined(_WIN32) || !defined(AVX512)
                        // WHFast512 not supported on Windows!
                        *(char**)pointer = malloc(sizeof(struct reb_particle_avx512));
#else 
                        *(char**)pointer = aligned_alloc(64,sizeof(struct reb_particle_avx512));
#endif // _WIN32
                    }else{ // normal malloc
                        *(char**)pointer = realloc(*(char**)pointer, field.size);
                    }
                    fread(*(char**)pointer, field.size,1,inf);

                    unsigned int* pointer_N = (unsigned int*)((char*)r + reb_binary_field_descriptor_list[i].offset_N);
                    *pointer_N = (unsigned int)field.size/reb_binary_field_descriptor_list[i].element_size;

                    goto next_field;
                }
                if (fd.dtype == REB_POINTER_FIXED_SIZE){
                    if (field.size != reb_binary_field_descriptor_list[i].element_size){
                        reb_simulation_warning(r, "Inconsistent size encountered in binary field (fixed pointer size).");
                    }
                    char* pointer = (char*)r + reb_binary_field_descriptor_list[i].offset;
                    *(char**)pointer = realloc(*(char**)pointer, field.size);
                    fread(*(char**)pointer, field.size,1,inf);

                    goto next_field;
                }
                if (fd.dtype == REB_CHARP_LIST){
                    size_t serialized_size = field.size;
                    char* serialized_strings = malloc(serialized_size);
                    fread(serialized_strings, serialized_size,1,inf);
                    // Process strings back into a list
                    char*** pointer = (char***)((char*)r + reb_binary_field_descriptor_list[i].offset);
                    unsigned int* pointer_N = (unsigned int*)((char*)r + reb_binary_field_descriptor_list[i].offset_N);
                    size_t current_pos = 0;
                    while (current_pos < serialized_size){
                        char* current_string = serialized_strings + current_pos;
                        // character count + NULL character + original pointer address
                        size_t current_string_len = strlen(current_string)+1+sizeof(char*);
                        current_pos += current_string_len;
                        // Add current_string to list
                        (*pointer_N)++;
                        *pointer = realloc(*pointer,sizeof(char*)*(*pointer_N));
                        (*pointer)[(*pointer_N)-1] = malloc(sizeof(char)*(current_string_len));
                        memcpy((*pointer)[(*pointer_N)-1], current_string, current_string_len);
                    }
                    free(serialized_strings);
                    goto next_field;
                }
                // Special datatype for ias15. Similar to REB_POINTER. 
                if (fd.dtype == REB_DP7){
                    if (field.size % reb_binary_field_descriptor_list[i].element_size){
                        reb_simulation_warning(r, "Inconsistent size encountered in binary field.");
                    }
                    char* pointer = (char*)r + reb_binary_field_descriptor_list[i].offset;
                    struct reb_dp7* dp7 = (struct reb_dp7*)pointer;

                    dp7->p0 = realloc(dp7->p0,field.size/7);
                    dp7->p1 = realloc(dp7->p1,field.size/7);
                    dp7->p2 = realloc(dp7->p2,field.size/7);
                    dp7->p3 = realloc(dp7->p3,field.size/7);
                    dp7->p4 = realloc(dp7->p4,field.size/7);
                    dp7->p5 = realloc(dp7->p5,field.size/7);
                    dp7->p6 = realloc(dp7->p6,field.size/7);
                    fread(dp7->p0, field.size/7, 1, inf);
                    fread(dp7->p1, field.size/7, 1, inf);
                    fread(dp7->p2, field.size/7, 1, inf);
                    fread(dp7->p3, field.size/7, 1, inf);
                    fread(dp7->p4, field.size/7, 1, inf);
                    fread(dp7->p5, field.size/7, 1, inf);
                    fread(dp7->p6, field.size/7, 1, inf);

                    unsigned int* pointer_N = (unsigned int*)((char*)r + reb_binary_field_descriptor_list[i].offset_N);
                    *pointer_N = (unsigned int)field.size/reb_binary_field_descriptor_list[i].element_size;

                    goto next_field;
                }
                // If we're here then it was not a simple or pointer datatype. 
                // Can skip the iteration trough the descriptor list.
                break;
            }
            i++;
        }

        // Fields with types that require special handling
        if (field.type == fd_functionpointers.type){
            // Warning for when function pointers were used. 
            // No effect on simulation.
            int fpwarn;
            fread(&fpwarn, field.size,1,inf);
            if (fpwarn && warnings){
                *warnings |= REB_SIMULATION_BINARY_WARNING_POINTERS;
            }
            goto next_field;
        }
        if (field.type == fd_header.type){
            // Check header.
            int64_t objects = 0;
            const size_t bufsize = 64 - sizeof(struct reb_binary_field);
            char readbuf[64], curvbuf[64];
            const char* header = "REBOUND Binary File. Version: ";
            sprintf(curvbuf,"%s%s",header+sizeof(struct reb_binary_field), reb_version_str);

            objects += fread(readbuf,sizeof(char),bufsize,inf);
            if (objects < 1){
                *warnings |= REB_SIMULATION_BINARY_WARNING_CORRUPTFILE;
            }else{
                // Note: following compares version, but ignores githash.
                if(strncmp(readbuf,curvbuf,bufsize)!=0){
                    *warnings |= REB_SIMULATION_BINARY_WARNING_VERSION;
                }
            }
            goto next_field;
        }

        // We should never get here. If so, it's an unknown field type.
        *warnings |= REB_SIMULATION_BINARY_WARNING_FIELD_UNKNOWN;
        int err = fseek(inf, field.size, SEEK_CUR);
        if (err){
            // Even worse, can't seek to end of field.
            *warnings |= REB_SIMULATION_BINARY_WARNING_CORRUPTFILE;
        }
    } 

finish_fields:
    // Some final initialization
    for (unsigned int l=0;l<r->N_var_config;l++){
        r->var_config[l].sim = r;
    }
    r->N_allocated = r->N; // This used to be different. Now only saving N.
    for (unsigned int l=0;l<r->N_allocated;l++){
        r->particles[l].c = NULL;
        r->particles[l].ap = NULL;
        r->particles[l].sim = r;
#ifndef MPI
        // Restore names
        if (r->particles[l].name){
            int name_found = 0;
            for (int n=0;n<r->N_name_list;n++){
                char* original_pointer = *((char**)(r->name_list[n]+strlen(r->name_list[n])+1));
                if (r->particles[l].name==original_pointer){
                    name_found = 1;
                    r->particles[l].name = r->name_list[n];
                }
            }
            if (!name_found){
                reb_simulation_warning(r,"A name for a particle was not stored in the Simulationarchive.");
                r->particles[l].name = NULL;
            }
        }
#endif // MPI
    }
    reb_tree_delete(r);
    if (r->gravity==REB_GRAVITY_TREE || r->collision==REB_COLLISION_TREE || r->collision==REB_COLLISION_LINETREE){
        for (unsigned int l=0;l<r->N_allocated;l++){
            reb_tree_add_particle_to_tree(r, l);
        }
    }
    // Commented out on Nov 26 2024. Not sure why this was added. Might be for an older SA version.
    // if (r->ri_ias15.at){ 
    //     // Assume that all arrays were saved whenever ri_ias15.at was saved.
    //     // Only 3*N entries got saved. 
    //     r->ri_ias15.N_allocated = 3*r->N;
    // }
    r->ri_whfast512.recalculate_constants = 1;
}

struct reb_simulation* reb_input_process_warnings(struct reb_simulation* r, enum reb_simulation_binary_error_codes warnings){
    if (warnings & REB_SIMULATION_BINARY_ERROR_NOFILE){
        reb_simulation_error(r,"Cannot read binary file. Check filename and file contents.");
        if (r) free(r);
        return NULL;
    }
    if (warnings & REB_SIMULATION_BINARY_WARNING_VERSION){
        reb_simulation_warning(r,"Binary file was saved with a different version of REBOUND. Binary format might have changed.");
    }
    if (warnings & REB_SIMULATION_BINARY_WARNING_POINTERS){
        reb_simulation_warning(r,"You have to reset function pointers after creating a reb_simulation struct with a binary file.");
    }
    if (warnings & REB_SIMULATION_BINARY_WARNING_PARTICLES){
        reb_simulation_warning(r,"Binary file might be corrupted. Number of particles found does not match expected number.");
    }
    if (warnings & REB_SIMULATION_BINARY_ERROR_FILENOTOPEN){
        reb_simulation_error(r,"Error while reading binary file (file was not open).");
        if (r) free(r);
        return NULL;
    }
    if (warnings & REB_SIMULATION_BINARY_ERROR_OUTOFRANGE){
        reb_simulation_error(r,"Index out of range.");
        if (r) free(r);
        return NULL;
    }
    if (warnings & REB_SIMULATION_BINARY_ERROR_SEEK){
        reb_simulation_error(r,"Error while trying to seek file.");
        if (r) free(r);
        return NULL;
    }
    if (warnings & REB_SIMULATION_BINARY_WARNING_FIELD_UNKNOWN){
        reb_simulation_warning(r,"Unknown field found in binary file.");
    }
    if (warnings & REB_SIMULATION_BINARY_ERROR_NOFILE){
        reb_simulation_error(r,"Cannot read binary file. Check filename and file contents.");
        if (r) free(r);
        return NULL;
    }
    if (warnings & REB_SIMULATION_BINARY_ERROR_OLD){
        reb_simulation_error(r,"Reading old Simulationarchives (version < 2) is no longer supported. If you need to read such an archive, use a REBOUND version <= 3.26.3");
        if (r) free(r);
        return NULL;
    }
    if (warnings & REB_SIMULATION_BINARY_WARNING_CORRUPTFILE){
        reb_simulation_warning(r,"The binary file seems to be corrupted. An attempt has been made to read the uncorrupted parts of it.");
    }
    return r;
}

