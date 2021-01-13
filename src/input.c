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
#ifdef MPI
#include "communication_mpi.h"
#endif

double reb_read_double(int argc, char** argv, const char* argument, double _default){
    char* value = reb_read_char(argc,argv,argument);
    if (value){
        return atof(value);
    }
    return _default;
}

int reb_read_int(int argc, char** argv, const char* argument, int _default){
    char* value = reb_read_char(argc,argv,argument);
    if (value){
        return atoi(value);
    }
    return _default;
}


char* reb_read_char(int argc, char** argv, const char* argument){
    opterr = 0;
    optind = 1;
    while (1) {
        struct option long_options[] = {
            {NULL, required_argument, 0, 'a'},
            {0,0,0,0}
        };

        long_options[0].name = argument;

        /* getopt_long stores the option index here.   */
        int option_index = 0;
        //              short options. format abc:d::
        int c = getopt_long (argc, argv, "", long_options, &option_index);

        /* Detect the end of the options.   */
        if (c == -1) break;

        switch (c){
            case 'a':
                return optarg;
                break;
            default:
                break;
        }
    }
    return NULL;
}

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

static int reb_fseek(FILE *stream, long offset, int whence, char **restrict mem_stream){
    if (mem_stream!=NULL){
        // read from memory
        if (whence==SEEK_CUR){
            *mem_stream = (char*)(*mem_stream)+offset;
            return 0;
        }
        return -1;
    }else if(stream!=NULL){
        // read from file
        return fseek(stream,offset,whence);
    }
    return -1;
}


void reb_read_dp7(struct reb_dp7* dp7, const int N3, FILE* inf, char **restrict mem_stream){
    reb_fread(dp7->p0,sizeof(double),N3,inf,mem_stream);
    reb_fread(dp7->p1,sizeof(double),N3,inf,mem_stream);
    reb_fread(dp7->p2,sizeof(double),N3,inf,mem_stream);
    reb_fread(dp7->p3,sizeof(double),N3,inf,mem_stream);
    reb_fread(dp7->p4,sizeof(double),N3,inf,mem_stream);
    reb_fread(dp7->p5,sizeof(double),N3,inf,mem_stream);
    reb_fread(dp7->p6,sizeof(double),N3,inf,mem_stream);
}

// Macro to read a single field from a binary file.
#define CASE(typename, value) case REB_BINARY_FIELD_TYPE_##typename: \
    {\
        reb_fread(value, field.size,1,inf,mem_stream);\
    }\
    break;

#define CASE_MALLOC(typename, valueref) case REB_BINARY_FIELD_TYPE_##typename: \
    {\
        valueref = malloc(field.size);\
        reb_fread(valueref, field.size,1,inf,mem_stream);\
    }\
    break;

#define CASE_MALLOC_DP7(typename, valueref) case REB_BINARY_FIELD_TYPE_##typename: \
    {\
        valueref.p0 = malloc(field.size/7);\
        valueref.p1 = malloc(field.size/7);\
        valueref.p2 = malloc(field.size/7);\
        valueref.p3 = malloc(field.size/7);\
        valueref.p4 = malloc(field.size/7);\
        valueref.p5 = malloc(field.size/7);\
        valueref.p6 = malloc(field.size/7);\
        reb_fread(valueref.p0, field.size/7,1,inf,mem_stream);\
        reb_fread(valueref.p1, field.size/7,1,inf,mem_stream);\
        reb_fread(valueref.p2, field.size/7,1,inf,mem_stream);\
        reb_fread(valueref.p3, field.size/7,1,inf,mem_stream);\
        reb_fread(valueref.p4, field.size/7,1,inf,mem_stream);\
        reb_fread(valueref.p5, field.size/7,1,inf,mem_stream);\
        reb_fread(valueref.p6, field.size/7,1,inf,mem_stream);\
    }\
    break;
    
int reb_input_field(struct reb_simulation* r, FILE* inf, enum reb_input_binary_messages* warnings, char **restrict mem_stream){
    struct reb_binary_field field;
    int numread = reb_fread(&field,sizeof(struct reb_binary_field),1,inf,mem_stream);
    if (numread<1){
        return 0; // End of file
    }
    switch (field.type){
        CASE(T,                  &r->t);
        CASE(G,                  &r->G);
        CASE(SOFTENING,          &r->softening);
        CASE(DT,                 &r->dt);
        CASE(DTLASTDONE,         &r->dt_last_done);
        CASE(N,                  &r->N);
        CASE(NVAR,               &r->N_var);
        CASE(VARCONFIGN,         &r->var_config_N);
        CASE(NACTIVE,            &r->N_active);
        CASE(TESTPARTICLETYPE,   &r->testparticle_type);
        CASE(HASHCTR,            &r->hash_ctr);
        CASE(OPENINGANGLE2,      &r->opening_angle2);
        CASE(STATUS,             &r->status);
        CASE(EXACTFINISHTIME,    &r->exact_finish_time);
        CASE(FORCEISVELOCITYDEP, &r->force_is_velocity_dependent);
        CASE(GRAVITYIGNORETERMS, &r->gravity_ignore_terms);
        CASE(OUTPUTTIMINGLAST,   &r->output_timing_last);
        CASE(SAVEMESSAGES,       &r->save_messages);
        CASE(EXITMAXDISTANCE,    &r->exit_max_distance);
        CASE(EXITMINDISTANCE,    &r->exit_min_distance);
        CASE(USLEEP,             &r->usleep);
        CASE(TRACKENERGYOFFSET,  &r->track_energy_offset);
        CASE(ENERGYOFFSET,       &r->energy_offset);
        CASE(BOXSIZE,            &r->boxsize);
        CASE(BOXSIZEMAX,         &r->boxsize_max);
        CASE(ROOTSIZE,           &r->root_size);
        CASE(ROOTN,              &r->root_n);
        CASE(ROOTNX,             &r->root_nx);
        CASE(ROOTNY,             &r->root_ny);
        CASE(ROOTNZ,             &r->root_nz);
        CASE(NGHOSTX,            &r->nghostx);
        CASE(NGHOSTY,            &r->nghosty);
        CASE(NGHOSTZ,            &r->nghostz);
        CASE(COLLISIONRESOLVEKEEPSORTED, &r->collision_resolve_keep_sorted);
        CASE(MINIMUMCOLLISIONVELOCITY, &r->minimum_collision_velocity);
        CASE(COLLISIONSPLOG,     &r->collisions_plog);
        CASE(MAXRADIUS,          &r->max_radius);
        CASE(COLLISIONSNLOG,     &r->collisions_Nlog);
        CASE(CALCULATEMEGNO,     &r->calculate_megno);
        CASE(MEGNOYS,            &r->megno_Ys);
        CASE(MEGNOYSS,           &r->megno_Yss);
        CASE(MEGNOCOVYT,         &r->megno_cov_Yt);
        CASE(MEGNOVART,          &r->megno_var_t);
        CASE(MEGNOMEANT,         &r->megno_mean_t);
        CASE(MEGNOMEANY,         &r->megno_mean_Y);
        CASE(MEGNON,             &r->megno_n);
        CASE(SAVERSION,          &r->simulationarchive_version);
        CASE(SASIZEFIRST,        &r->simulationarchive_size_first);
        CASE(SASIZESNAPSHOT,     &r->simulationarchive_size_snapshot);
        CASE(SAAUTOINTERVAL,     &r->simulationarchive_auto_interval);
        CASE(SAAUTOWALLTIME,     &r->simulationarchive_auto_walltime);
        CASE(SANEXT,             &r->simulationarchive_next);
        CASE(WALLTIME,           &r->walltime);
        CASE(COLLISION,          &r->collision);
        CASE(VISUALIZATION,      &r->visualization);
        case REB_BINARY_FIELD_TYPE_INTEGRATOR: 
        {
            reb_fread(&r->integrator, field.size,1,inf,mem_stream);\
            // This is for backwards compatibility. To be removed in the future.
            if (r->integrator == REB_INTEGRATOR_IAS15){
                r->ri_ias15.neworder = 0;
            }
        }
        break;
        CASE(BOUNDARY,           &r->boundary);
        CASE(GRAVITY,            &r->gravity);
        CASE(SEI_OMEGA,          &r->ri_sei.OMEGA);
        CASE(SEI_OMEGAZ,         &r->ri_sei.OMEGAZ);
        CASE(SEI_LASTDT,         &r->ri_sei.lastdt);
        CASE(SEI_SINDT,          &r->ri_sei.sindt);
        CASE(SEI_TANDT,          &r->ri_sei.tandt);
        CASE(SEI_SINDTZ,         &r->ri_sei.sindtz);
        CASE(SEI_TANDTZ,         &r->ri_sei.tandtz);
        CASE(WHFAST_CORRECTOR,   &r->ri_whfast.corrector);
        CASE(WHFAST_RECALCJAC,   &r->ri_whfast.recalculate_coordinates_this_timestep);
        CASE(WHFAST_SAFEMODE,    &r->ri_whfast.safe_mode);
        CASE(WHFAST_KEEPUNSYNC,  &r->ri_whfast.keep_unsynchronized);
        CASE(WHFAST_ISSYNCHRON,  &r->ri_whfast.is_synchronized);
        CASE(WHFAST_TIMESTEPWARN,&r->ri_whfast.timestep_warning);
        CASE(WHFAST_COORDINATES, &r->ri_whfast.coordinates);
        CASE(IAS15_EPSILON,      &r->ri_ias15.epsilon);
        CASE(IAS15_MINDT,        &r->ri_ias15.min_dt);
        CASE(IAS15_EPSILONGLOBAL,&r->ri_ias15.epsilon_global);
        CASE(IAS15_ITERATIONSMAX,&r->ri_ias15.iterations_max_exceeded);
        CASE(IAS15_ALLOCATEDN,   &r->ri_ias15.allocatedN);
        CASE(IAS15_NEWORDER,     &r->ri_ias15.neworder);
        CASE(JANUS_SCALEPOS,     &r->ri_janus.scale_pos);
        CASE(JANUS_SCALEVEL,     &r->ri_janus.scale_vel);
        CASE(JANUS_ORDER,        &r->ri_janus.order);
        CASE(JANUS_ALLOCATEDN,   &r->ri_janus.allocated_N);
        CASE(JANUS_RECALC,       &r->ri_janus.recalculate_integer_coordinates_this_timestep);
        CASE(MERCURIUS_HILLFAC,  &r->ri_mercurius.hillfac);
        CASE(MERCURIUS_SAFEMODE, &r->ri_mercurius.safe_mode);
        CASE(MERCURIUS_ISSYNCHRON, &r->ri_mercurius.is_synchronized);
        CASE(MERCURIUS_COMPOS,   &r->ri_mercurius.com_pos);
        CASE(MERCURIUS_COMVEL,   &r->ri_mercurius.com_vel);
        CASE(PYTHON_UNIT_L,      &r->python_unit_l);
        CASE(PYTHON_UNIT_M,      &r->python_unit_m);
        CASE(PYTHON_UNIT_T,      &r->python_unit_t);
        CASE(STEPSDONE,          &r->steps_done);
        CASE(SAAUTOSTEP,         &r->simulationarchive_auto_step);
        CASE(SANEXTSTEP,         &r->simulationarchive_next_step);
        CASE(SABA_TYPE,          &r->ri_saba.type);
        CASE(SABA_KEEPUNSYNC,    &r->ri_saba.keep_unsynchronized);
        CASE(EOS_PHI0,           &r->ri_eos.phi0);
        CASE(EOS_PHI1,           &r->ri_eos.phi1);
        CASE(EOS_N,              &r->ri_eos.n);
        CASE(EOS_SAFEMODE,       &r->ri_eos.safe_mode);
        CASE(EOS_ISSYNCHRON,     &r->ri_eos.is_synchronized);
        CASE(RAND_SEED,          &r->rand_seed);
        // temporary solution for depreciated SABA k and corrector variables.
        // can be removed in future versions
        case 138: 
            {
            unsigned int k = 0;
            reb_fread(&k, field.size,1,inf,mem_stream);
            r->ri_saba.type/=0x100;
            r->ri_saba.type += k-1;
            }
            break;
        case 139: 
            {
            unsigned int corrector = 0;
            reb_fread(&corrector, field.size,1,inf,mem_stream);
            r->ri_saba.type%=0x100;
            r->ri_saba.type += 0x100*corrector;
            }
            break;

        CASE(SABA_SAFEMODE,      &r->ri_saba.safe_mode);
        CASE(SABA_ISSYNCHRON,    &r->ri_saba.is_synchronized);
        CASE(WHFAST_CORRECTOR2,  &r->ri_whfast.corrector2);
        CASE(WHFAST_KERNEL,      &r->ri_whfast.kernel);
        case REB_BINARY_FIELD_TYPE_PARTICLES:
            if(r->particles){
                free(r->particles);
            }
            r->allocatedN = (int)(field.size/sizeof(struct reb_particle));
            if (field.size){
                r->particles = malloc(field.size);
                reb_fread(r->particles, field.size,1,inf,mem_stream);
            }
            if (r->allocatedN<r->N && warnings){
                *warnings |= REB_INPUT_BINARY_WARNING_PARTICLES;
            }
            for (int l=0;l<r->allocatedN;l++){
                r->particles[l].c = NULL;
                r->particles[l].ap = NULL;
                r->particles[l].sim = r;
            }
            if (r->gravity==REB_GRAVITY_TREE || r->collision==REB_COLLISION_TREE || r->collision==REB_COLLISION_LINETREE){
                for (int l=0;l<r->allocatedN;l++){
                    reb_tree_add_particle_to_tree(r, l);
                }
            }
            break;
        case REB_BINARY_FIELD_TYPE_WHFAST_PJ:
            if(r->ri_whfast.p_jh){
                free(r->ri_whfast.p_jh);
            }
            r->ri_whfast.allocated_N = (int)(field.size/sizeof(struct reb_particle));
            if (field.size){
                r->ri_whfast.p_jh = malloc(field.size);
                reb_fread(r->ri_whfast.p_jh, field.size,1,inf,mem_stream);
            }
            break;
        case REB_BINARY_FIELD_TYPE_JANUS_PINT:
            if(r->ri_janus.p_int){
                free(r->ri_janus.p_int);
            }
            r->ri_janus.allocated_N = (int)(field.size/sizeof(struct reb_particle_int));
            if (field.size){
                r->ri_janus.p_int = malloc(field.size);
                reb_fread(r->ri_janus.p_int, field.size,1,inf,mem_stream);
            }
            break;
        case REB_BINARY_FIELD_TYPE_VARCONFIG:
            if (r->var_config){
                free(r->var_config);
            }
            if (r->var_config_N>0){
                r->var_config = malloc(field.size);
                reb_fread(r->var_config, field.size,1,inf,mem_stream);
                for (int l=0;l<r->var_config_N;l++){
                    r->var_config[l].sim = r;
                }
            }
            break;
        case REB_BINARY_FIELD_TYPE_MERCURIUS_DCRIT:
            if(r->ri_mercurius.dcrit){
                free(r->ri_mercurius.dcrit);
            }
            r->ri_mercurius.dcrit_allocatedN = (int)(field.size/sizeof(double));
            if (field.size){
                r->ri_mercurius.dcrit = malloc(field.size);
                reb_fread(r->ri_mercurius.dcrit, field.size,1,inf,mem_stream);
            }
            break;
        CASE_MALLOC(IAS15_AT,     r->ri_ias15.at);
        CASE_MALLOC(IAS15_X0,     r->ri_ias15.x0);
        CASE_MALLOC(IAS15_V0,     r->ri_ias15.v0);
        CASE_MALLOC(IAS15_A0,     r->ri_ias15.a0);
        CASE_MALLOC(IAS15_CSX,    r->ri_ias15.csx);
        CASE_MALLOC(IAS15_CSV,    r->ri_ias15.csv);
        CASE_MALLOC(IAS15_CSA0,   r->ri_ias15.csa0);
        CASE_MALLOC_DP7(IAS15_G,  r->ri_ias15.g);
        CASE_MALLOC_DP7(IAS15_B,  r->ri_ias15.b);
        CASE_MALLOC_DP7(IAS15_CSB,r->ri_ias15.csb);
        CASE_MALLOC_DP7(IAS15_E,  r->ri_ias15.e);
        CASE_MALLOC_DP7(IAS15_BR, r->ri_ias15.br);
        CASE_MALLOC_DP7(IAS15_ER, r->ri_ias15.er);
        case REB_BINARY_FIELD_TYPE_END:
            return 0;
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
                // Note: following compares version, but ignores githash.
                if(strncmp(readbuf,curvbuf,bufsize)!=0){
                    *warnings |= REB_INPUT_BINARY_WARNING_VERSION;
                }
            }
            break;
        default:
            if (warnings){
                *warnings |= REB_INPUT_BINARY_WARNING_FIELD_UNKOWN;
            }
            reb_fseek(inf,field.size,SEEK_CUR,mem_stream);
            break;
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

