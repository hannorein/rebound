/**
 * @file    input.c
 * @brief   Parse command line options and read retart files.
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
#include <getopt.h>
#include <string.h>
#include "particle.h"
#include "rebound.h"
#include "collision.h"
#include "input.h"
#include "tree.h"
#include "simulationarchive.h"
#include "integrator_tes.h"

#ifdef MPI
#include "communication_mpi.h"
#endif

static size_t reb_fread(void *restrict ptr, size_t size, size_t nitems, FILE *restrict stream, char **restrict mem_stream){
    if (mem_stream!=NULL){
        // read from memory
        memcpy(ptr,*mem_stream,size*nitems);
        *mem_stream = (char*)(*mem_stream)+ size*nitems;
        return size*nitems;
    }else if(stream!=NULL){
        // read from file
        return fread(ptr,size,nitems,stream);
    }
    return 0; 
}

// The reb_fseek function is currently not used, but provides functionality along the lines of reb_fread
// static int reb_fseek(FILE *stream, long offset, int whence, char **restrict mem_stream){
//     if (mem_stream!=NULL){
//         // read from memory
//         if (whence==SEEK_CUR){
//             *mem_stream = (char*)(*mem_stream)+offset;
//             return 0;
//         }
//         return -1;
//     }else if(stream!=NULL){
//         // read from file
//         return fseek(stream,offset,whence);
//     }
//     return -1;
// }


// Macro to read a single field from a binary file.
#define CASE(typename, value) case REB_BINARY_FIELD_TYPE_##typename: \
    {\
        reb_fread(value, field.size,1,inf,mem_stream);\
    }\
    break;

#define CASE_CONTROL_VARS(typename, valueref) case REB_BINARY_FIELD_TYPE_##typename: \
    {\
        reb_fread(&valueref->size, sizeof(uint32_t),1,inf,mem_stream);\
        reb_fread(valueref->p0, valueref->size,1,inf,mem_stream);\
        reb_fread(valueref->p1, valueref->size,1,inf,mem_stream);\
        reb_fread(valueref->p2, valueref->size,1,inf,mem_stream);\
        reb_fread(valueref->p3, valueref->size,1,inf,mem_stream);\
        reb_fread(valueref->p4, valueref->size,1,inf,mem_stream);\
        reb_fread(valueref->p5, valueref->size,1,inf,mem_stream);\
        reb_fread(valueref->p6, valueref->size,1,inf,mem_stream);\
    }\
    break;        

void reb_input_field_finish(struct reb_simulation* r, enum reb_input_binary_messages* warnings){
    for (int l=0;l<r->var_config_N;l++){
        r->var_config[l].sim = r;
    }
    r->allocatedN = r->N; // This used to be different. Now only saving N.
    for (unsigned int l=0;l<r->allocatedN;l++){
        r->particles[l].c = NULL;
        r->particles[l].ap = NULL;
        r->particles[l].sim = r;
    }
    if (r->gravity==REB_GRAVITY_TREE || r->collision==REB_COLLISION_TREE || r->collision==REB_COLLISION_LINETREE){
        for (unsigned int l=0;l<r->allocatedN;l++){
            reb_tree_add_particle_to_tree(r, l);
        }
    }
    if (r->ri_ias15.at){ 
        // Assume that all arrays were saved whenever ri_ias15.at was saved.
        // Only 3*N entries got saved. 
        r->ri_ias15.allocatedN = 3*r->N;
    }
    r->ri_whfast512.recalculate_constants = 1;
}


int reb_input_field(struct reb_simulation* r, FILE* inf, enum reb_input_binary_messages* warnings, char **restrict mem_stream){
    struct reb_binary_field field;
    int numread = reb_fread(&field,sizeof(struct reb_binary_field),1,inf,mem_stream);
    if (numread<1){
        return 0; // End of file
    }
    if (field.type==9999){
        return 0; // End of snapshot
    }
        // only here for testding. delete TODO 
        // only here for testding. delete TODO 
    if (field.type==REB_BINARY_FIELD_TYPE_TES_ALLOCATED_N){
        reb_fread(&r->ri_tes.allocated_N, field.size, 1, inf, mem_stream);
        // Allocate all memory for loading the simulation archive.
        if (r->ri_tes.allocated_N) {
            reb_integrator_tes_allocate_memory(r);
        }
        printf("allocate new memory %d\n", r->ri_tes.allocated_N);
        return 1;
    }
    int i=0;
    while (reb_binary_field_descriptor_list[i].type!=9999){
        int type = reb_binary_field_descriptor_list[i].type;
        int dtype = reb_binary_field_descriptor_list[i].dtype;
        if (type==field.type){
            // Read simple data types
            if (dtype == REB_DOUBLE || dtype == REB_INT || dtype == REB_UINT || dtype == REB_UINT32 ||
                    dtype == REB_LONG || dtype == REB_ULONG || dtype == REB_ULONGLONG || 
                    dtype == REB_PARTICLE || dtype == REB_VEC3D ){
                char* pointer = (char*)r + reb_binary_field_descriptor_list[i].offset;
                reb_fread(pointer, field.size, 1, inf ,mem_stream);
                return 1;
            }
            // Read a pointer data type. 
            // 1) reallocate memory
            // 2) read data into memory
            // 3) set allocated_N variable
            if (dtype == REB_POINTER || dtype == REB_POINTER_ALIGNED){
                if (field.size % reb_binary_field_descriptor_list[i].element_size){
                    reb_warning(r, "Inconsistent size encountered in binary field.");
                }
                char* pointer = (char*)r + reb_binary_field_descriptor_list[i].offset;
                if (dtype == REB_POINTER_ALIGNED){
                    if (*(char**)pointer) free(*(char**)pointer);
                    *(char**)pointer = aligned_alloc(64,sizeof(struct reb_particle_avx512));
                }else{ // normal malloc
                    *(char**)pointer = realloc(*(char**)pointer, field.size);
                }
                reb_fread(*(char**)pointer, field.size,1,inf,mem_stream);
                
                unsigned int* pointer_N = (unsigned int*)((char*)r + reb_binary_field_descriptor_list[i].offset_N);
                *pointer_N = field.size/reb_binary_field_descriptor_list[i].element_size;

                return 1;
            }
            // Special datatype for ias15. Similar to REB_POINTER. 
            if (dtype == REB_DP7){
                if (field.size % reb_binary_field_descriptor_list[i].element_size){
                    reb_warning(r, "Inconsistent size encountered in binary field.");
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
                reb_fread(dp7->p0, field.size/7, 1, inf, mem_stream);
                reb_fread(dp7->p1, field.size/7, 1, inf, mem_stream);
                reb_fread(dp7->p2, field.size/7, 1, inf, mem_stream);
                reb_fread(dp7->p3, field.size/7, 1, inf, mem_stream);
                reb_fread(dp7->p4, field.size/7, 1, inf, mem_stream);
                reb_fread(dp7->p5, field.size/7, 1, inf, mem_stream);
                reb_fread(dp7->p6, field.size/7, 1, inf, mem_stream);
            
                unsigned int* pointer_N = (unsigned int*)((char*)r + reb_binary_field_descriptor_list[i].offset_N);
                *pointer_N = field.size/reb_binary_field_descriptor_list[i].element_size;

                return 1;
            }
        }
        i++;
    }

    switch (field.type){
        case 35:
            { // Only kept for backwards compatability. Can be removed in future version.
                double max_radius[2];
                reb_fread(&max_radius, field.size,1,inf,mem_stream);
                r->max_radius0 = max_radius[0];
                r->max_radius1 = max_radius[1];
            }
            break;
        case 45: // simulationarchive_size_first was manually written. reading it manually here.
            reb_fread(&r->simulationarchive_size_first, field.size,1,inf,mem_stream);
            break;
        case REB_BINARY_FIELD_TYPE_FUNCTIONPOINTERS:
            {
                int fpwarn;
                reb_fread(&fpwarn, field.size,1,inf,mem_stream);
                if (fpwarn && warnings){
                    *warnings |= REB_INPUT_BINARY_WARNING_POINTERS;
                }
            }
            break;
        case REB_BINARY_FIELD_TYPE_HEADER:
            {
                long objects = 0;
                // Input header.
                const long bufsize = 64 - sizeof(struct reb_binary_field);
                char readbuf[bufsize], curvbuf[bufsize];
                const char* header = "REBOUND Binary File. Version: ";
                sprintf(curvbuf,"%s%s",header+sizeof(struct reb_binary_field), reb_version_str);
                
                objects += reb_fread(readbuf,sizeof(char),bufsize,inf,mem_stream);
                if (objects < 1){
                    *warnings |= REB_INPUT_BINARY_WARNING_CORRUPTFILE;
                }else{
                    // Note: following compares version, but ignores githash.
                    if(strncmp(readbuf,curvbuf,bufsize)!=0){
                        *warnings |= REB_INPUT_BINARY_WARNING_VERSION;
                    }
                }
            }
            break;

        // TES integrator variables
        CASE(TES_COM,                &r->ri_tes.COM);
        CASE(TES_COM_DOT,            &r->ri_tes.COM_dot);   

        

        CASE(TES_PARTICLES_DH, r->ri_tes.particles_dh);
        //CASE(TES_MASS, r->ri_tes.mass);
        CASE(TES_X_DH, r->ri_tes.X_dh); 
        
        // TES Kepler vars
        CASE(TES_UVARS_SV_SIZE, &r->ri_tes.uVars->stateVectorSize);
        CASE(TES_UVARS_T0, r->ri_tes.uVars->t0);
        CASE(TES_UVARS_TLAST, r->ri_tes.uVars->tLast);
        CASE(TES_UVARS_CSQ, r->ri_tes.uVars->uv_csq);
        CASE(TES_UVARS_CSP, r->ri_tes.uVars->uv_csp);
        CASE(TES_UVARS_CSV, r->ri_tes.uVars->uv_csv);
        CASE(TES_UVARS_Q0, r->ri_tes.uVars->Q0);
        CASE(TES_UVARS_V0, r->ri_tes.uVars->V0);
        CASE(TES_UVARS_P0, r->ri_tes.uVars->P0);
        CASE(TES_UVARS_Q1, r->ri_tes.uVars->Q1);
        CASE(TES_UVARS_V1, r->ri_tes.uVars->V1);
        CASE(TES_UVARS_P1, r->ri_tes.uVars->P1);
        CASE(TES_UVARS_X, r->ri_tes.uVars->X);
        CASE(TES_UVARS_Q0_NORM, r->ri_tes.uVars->Q0_norm);
        CASE(TES_UVARS_BETA, r->ri_tes.uVars->beta);
        CASE(TES_UVARS_ETA, r->ri_tes.uVars->eta);
        CASE(TES_UVARS_ZETA, r->ri_tes.uVars->zeta);
        CASE(TES_UVARS_PERIOD, r->ri_tes.uVars->period);
        CASE(TES_UVARS_XPERIOD, r->ri_tes.uVars->Xperiod);
        CASE(TES_UVARS_STUMPF_C0, r->ri_tes.uVars->C.c0);
        CASE(TES_UVARS_STUMPF_C1, r->ri_tes.uVars->C.c1);
        CASE(TES_UVARS_STUMPF_C2, r->ri_tes.uVars->C.c2);
        CASE(TES_UVARS_STUMPF_C3, r->ri_tes.uVars->C.c3);
        CASE(TES_UVARS_MU, &r->ri_tes.uVars->mu);

        // TES Radau vars
        CASE(TES_RADAU_DX, r->ri_tes.radau->dX);
        CASE(TES_RADAU_XOUT, r->ri_tes.radau->Xout);
        CASE(TES_RADAU_RECTI_ARRAY, r->ri_tes.radau->rectifiedArray);
        CASE(TES_RADAU_PREDICTORS, r->ri_tes.radau->predictors);
        CASE(TES_RADAU_DSTATE0, r->ri_tes.radau->dState0);
        CASE(TES_RADAU_DDSTATE0, r->ri_tes.radau->ddState0);
        CASE(TES_RADAU_DSTATE, r->ri_tes.radau->dState);
        CASE(TES_RADAU_DDSTATE, r->ri_tes.radau->ddState);
        CASE(TES_RADAU_CS_DSTATE0, r->ri_tes.radau->cs_dState0);
        CASE(TES_RADAU_CS_DDSTATE0, r->ri_tes.radau->cs_ddState0);
        CASE(TES_RADAU_CS_DSTATE, r->ri_tes.radau->cs_dState);
        CASE(TES_RADAU_CS_DDSTATE, r->ri_tes.radau->cs_ddState);
        CASE(TES_RADAU_CS_DX, r->ri_tes.radau->cs_dX);
        CASE(TES_RADAU_FCALLS, &r->ri_tes.radau->fCalls);
        CASE(TES_RADAU_RECTIS, &r->ri_tes.radau->rectifications);
        CASE(TES_RADAU_ITERS, &r->ri_tes.radau->convergenceIterations);
        CASE(TES_RADAU_B6, r->ri_tes.radau->b6_store);
        CASE_CONTROL_VARS(TES_RADAU_B, (&(r->ri_tes.radau->B)));
        CASE_CONTROL_VARS(TES_RADAU_BLAST, (&(r->ri_tes.radau->Blast)));
        CASE_CONTROL_VARS(TES_RADAU_B_1ST, (&(r->ri_tes.radau->B_1st)));
        CASE_CONTROL_VARS(TES_RADAU_BLAST_1ST, (&(r->ri_tes.radau->Blast_1st)));
        CASE_CONTROL_VARS(TES_RADAU_CS_B, (&(r->ri_tes.radau->cs_B)));
        CASE_CONTROL_VARS(TES_RADAU_CS_B_1ST, (&(r->ri_tes.radau->cs_B1st)));
        CASE_CONTROL_VARS(TES_RADAU_G, (&(r->ri_tes.radau->G)));
        CASE_CONTROL_VARS(TES_RADAU_G_1ST, (&(r->ri_tes.radau->G_1st)));

        // TES force model vars
        CASE(TES_DHEM_XOSC_STORE, r->ri_tes.rhs->XoscStore);
        CASE(TES_DHEM_XOSC_PRED_STORE, r->ri_tes.rhs->XoscPredStore);
        CASE(TES_DHEM_XOSC_CS_STORE, r->ri_tes.rhs->XoscStore_cs);
        CASE(TES_DHEM_XOSC_DOT_STORE, r->ri_tes.rhs->Xosc_dotStore);
        CASE(TES_DHEM_X, r->ri_tes.rhs->X);
        CASE(TES_DHEM_M_INV, r->ri_tes.rhs->m_inv);
        CASE(TES_DHEM_M_TOTAL, &r->ri_tes.rhs->mTotal);
        CASE(TES_DHEM_RECTI_TIME, r->ri_tes.rhs->rectifyTimeArray);
        CASE(TES_DHEM_RECTI_PERIOD, r->ri_tes.rhs->rectificationPeriod);
    
    }
    return 1;
} 

struct reb_simulation* reb_input_process_warnings(struct reb_simulation* r, enum reb_input_binary_messages warnings){
    if (warnings & REB_INPUT_BINARY_ERROR_NOFILE){
        reb_error(r,"Cannot read binary file. Check filename and file contents.");
        if (r) free(r);
        return NULL;
    }
    if (warnings & REB_INPUT_BINARY_WARNING_VERSION){
        reb_warning(r,"Binary file was saved with a different version of REBOUND. Binary format might have changed.");
    }
    if (warnings & REB_INPUT_BINARY_WARNING_POINTERS){
        reb_warning(r,"You have to reset function pointers after creating a reb_simulation struct with a binary file.");
    }
    if (warnings & REB_INPUT_BINARY_WARNING_PARTICLES){
        reb_warning(r,"Binary file might be corrupted. Number of particles found does not match expected number.");
    }
    if (warnings & REB_INPUT_BINARY_ERROR_FILENOTOPEN){
        reb_error(r,"Error while reading binary file (file was not open).");
        if (r) free(r);
        return NULL;
    }
    if (warnings & REB_INPUT_BINARY_ERROR_OUTOFRANGE){
        reb_error(r,"Index out of range.");
        if (r) free(r);
        return NULL;
    }
    if (warnings & REB_INPUT_BINARY_ERROR_SEEK){
        reb_error(r,"Error while trying to seek file.");
        if (r) free(r);
        return NULL;
    }
    if (warnings & REB_INPUT_BINARY_WARNING_FIELD_UNKOWN){
        reb_warning(r,"Unknown field found in binary file.");
    }
    if (warnings & REB_INPUT_BINARY_ERROR_NOFILE){
        reb_error(r,"Cannot read binary file. Check filename and file contents.");
        if (r) free(r);
        return NULL;
    }
    if (warnings & REB_INPUT_BINARY_ERROR_OLD){
        reb_error(r,"Reading old SimulationArchives (version < 2) is no longer supported. If you need to read such an archive, use a REBOUND version <= 3.26.3");
        if (r) free(r);
        return NULL;
    }
    if (warnings & REB_INPUT_BINARY_WARNING_CORRUPTFILE){
        reb_warning(r,"The binary file seems to be corrupted. An attempt has been made to read the uncorrupted parts of it.");
    }
    return r;
}

struct reb_simulation* reb_create_simulation_from_binary(char* filename){
    enum reb_input_binary_messages warnings = REB_INPUT_BINARY_WARNING_NONE;
    struct reb_simulation* r = reb_create_simulation();
    
    struct reb_simulationarchive* sa = malloc(sizeof(struct reb_simulationarchive)); 
    reb_read_simulationarchive_with_messages(sa, filename, NULL, &warnings);
    if (warnings & REB_INPUT_BINARY_ERROR_NOFILE){
        // Don't output an error if file does not exist, just return NULL.
        free(sa);
        return NULL;
    }else{
        reb_input_process_warnings(NULL, warnings);
    }
    reb_create_simulation_from_simulationarchive_with_messages(r, sa, -1, &warnings);
    reb_close_simulationarchive(sa);
    r = reb_input_process_warnings(r, warnings);
    return r;
}

