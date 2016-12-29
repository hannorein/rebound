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

void reb_read_dp7(struct reb_dp7* dp7, const int N3, FILE* inf){
    dp7->p0 = malloc(sizeof(double)*N3);
    dp7->p1 = malloc(sizeof(double)*N3);
    dp7->p2 = malloc(sizeof(double)*N3);
    dp7->p3 = malloc(sizeof(double)*N3);
    dp7->p4 = malloc(sizeof(double)*N3);
    dp7->p5 = malloc(sizeof(double)*N3);
    dp7->p6 = malloc(sizeof(double)*N3);
    fread(dp7->p0,sizeof(double),N3,inf);
    fread(dp7->p1,sizeof(double),N3,inf);
    fread(dp7->p2,sizeof(double),N3,inf);
    fread(dp7->p3,sizeof(double),N3,inf);
    fread(dp7->p4,sizeof(double),N3,inf);
    fread(dp7->p5,sizeof(double),N3,inf);
    fread(dp7->p6,sizeof(double),N3,inf);
}

// Macro to read a single field from a binary file.
#define CASE(typename, value) case REB_BINARY_FIELD_TYPE_##typename: \
    {\
        fread(value, field.size,1,inf);\
    }\
    break;

#define CASE_MALLOC(typename, valueref) case REB_BINARY_FIELD_TYPE_##typename: \
    {\
        valueref = malloc(field.size);\
        fread(valueref, field.size,1,inf);\
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
        fread(valueref.p0, field.size/7,1,inf);\
        fread(valueref.p1, field.size/7,1,inf);\
        fread(valueref.p2, field.size/7,1,inf);\
        fread(valueref.p3, field.size/7,1,inf);\
        fread(valueref.p4, field.size/7,1,inf);\
        fread(valueref.p5, field.size/7,1,inf);\
        fread(valueref.p6, field.size/7,1,inf);\
    }\
    break;
    

void reb_create_simulation_from_binary_with_messages(struct reb_simulation* r, char* filename, enum reb_input_binary_messages* warnings){
    FILE* inf = fopen(filename,"rb"); 
    
    if (!inf){
        *warnings |= REB_INPUT_BINARY_ERROR_NOFILE;
        return;
    }

    long objects = 0;
    // Input header.
    const char str[] = "REBOUND Binary File. Version: ";
    const char zero = '\0';
    char readbuf[65], curvbuf[65];
    sprintf(curvbuf,"%s%s",str,reb_version_str);
    memcpy(curvbuf+strlen(curvbuf)+1,reb_githash_str,sizeof(char)*(62-strlen(curvbuf)));
    curvbuf[63] = zero;
    
    objects += fread(readbuf,sizeof(char),64,inf);
    // Note: following compares version, but ignores githash.
    if(strcmp(readbuf,curvbuf)!=0){
        *warnings |= REB_INPUT_BINARY_WARNING_VERSION;
    }
    
    reb_reset_temporary_pointers(r);
    reb_reset_function_pointers(r);
    r->simulationarchive_filename = NULL;

    int reading_fields = 1;
    while(reading_fields){
        struct reb_binary_field field;
        fread(&field,sizeof(struct reb_binary_field),1,inf);
        switch (field.type){
            CASE(T,                  &r->t);
            CASE(G,                  &r->G);
            CASE(SOFTENING,          &r->softening);
            CASE(DT,                 &r->dt);
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
            CASE(SASIZEFIRST,        &r->simulationarchive_size_first);
            CASE(SASIZESNAPSHOT,     &r->simulationarchive_size_snapshot);
            CASE(SAINTERVAL,         &r->simulationarchive_interval);
            CASE(SAINTERVALWALLTIME, &r->simulationarchive_interval_walltime);
            CASE(SANEXT,             &r->simulationarchive_next);
            CASE(SAWALLTIME,         &r->simulationarchive_walltime);
            CASE(COLLISION,          &r->collision);
            CASE(INTEGRATOR,         &r->integrator);
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
            CASE(WHFAST_RECALCJAC,   &r->ri_whfast.recalculate_jacobi_this_timestep);
            CASE(WHFAST_SAFEMODE,    &r->ri_whfast.safe_mode);
            CASE(WHFAST_KEEPUNSYNC,  &r->ri_whfast.keep_unsynchronized);
            CASE(WHFAST_ISSYNCHRON,  &r->ri_whfast.is_synchronized);
            CASE(WHFAST_TIMESTEPWARN,&r->ri_whfast.timestep_warning);
            CASE(IAS15_EPSILON,      &r->ri_ias15.epsilon);
            CASE(IAS15_MINDT,        &r->ri_ias15.min_dt);
            CASE(IAS15_EPSILONGLOBAL,&r->ri_ias15.epsilon_global);
            CASE(IAS15_ITERATIONSMAX,&r->ri_ias15.iterations_max_exceeded);
            CASE(IAS15_ALLOCATEDN,   &r->ri_ias15.allocatedN);
            CASE(HERMES_HSF,         &r->ri_hermes.hill_switch_factor);
            CASE(HERMES_SSF,         &r->ri_hermes.solar_switch_factor);
            CASE(HERMES_ADAPTIVE,    &r->ri_hermes.adaptive_hill_switch_factor);
            CASE(HERMES_TIMESTEPWARN,&r->ri_hermes.timestep_too_large_warning);
            CASE(HERMES_STEPS,       &r->ri_hermes.steps);
            CASE(HERMES_STEPS_MA,    &r->ri_hermes.steps_miniactive);
            CASE(HERMES_STEPS_MN,    &r->ri_hermes.steps_miniN);
            CASE(WHFASTH_CORRECTOR,  &r->ri_whfasthelio.corrector);
            CASE(WHFASTH_RECALCHELIO,&r->ri_whfasthelio.recalculate_heliocentric_this_timestep);
            CASE(WHFASTH_SAFEMODE,   &r->ri_whfasthelio.safe_mode);
            CASE(WHFASTH_ISSYNCHRON, &r->ri_whfasthelio.is_synchronized);
            CASE(WHFASTH_KEEPUNSYNC, &r->ri_whfasthelio.keep_unsynchronized);
            case REB_BINARY_FIELD_TYPE_PARTICLES:
                if(r->particles){
                    free(r->particles);
                }
                r->particles = malloc(field.size);
                r->allocatedN = (int)(field.size/sizeof(struct reb_particle));
                if (r->allocatedN<r->N){
                    *warnings |= REB_INPUT_BINARY_WARNING_PARTICLES;
                }
                fread(r->particles, field.size,1,inf);
                for (int l=0;l<r->allocatedN;l++){
                    r->particles[l].c = NULL;
                    r->particles[l].ap = NULL;
                    r->particles[l].sim = r;
                }
                break;
            case REB_BINARY_FIELD_TYPE_WHFAST_PJ:
                if(r->ri_whfast.p_j){
                    free(r->ri_whfast.p_j);
                }
                r->ri_whfast.p_j = malloc(field.size);
                r->ri_whfast.allocated_N = (int)(field.size/sizeof(struct reb_particle));
                fread(r->ri_whfast.p_j, field.size,1,inf);
                break;
            case REB_BINARY_FIELD_TYPE_WHFAST_ETA:
                if(r->ri_whfast.eta){
                    free(r->ri_whfast.eta);
                }
                r->ri_whfast.eta = malloc(field.size);
                r->ri_whfast.allocated_N = (int)(field.size/sizeof(double));
                fread(r->ri_whfast.eta, field.size,1,inf);
                break;
            case REB_BINARY_FIELD_TYPE_WHFASTH_PH:
                if(r->ri_whfasthelio.p_h){
                    free(r->ri_whfasthelio.p_h);
                }
                r->ri_whfasthelio.p_h = malloc(field.size);
                r->ri_whfasthelio.allocated_N = (int)(field.size/sizeof(struct reb_particle));
                fread(r->ri_whfasthelio.p_h, field.size,1,inf);
                break;
            case REB_BINARY_FIELD_TYPE_VARCONFIG:
                if (r->var_config){
                    free(r->var_config);
                }
                r->var_config = malloc(field.size);
                fread(r->var_config, field.size,1,inf);
                if (r->var_config_N>0){
                    for (int l=0;l<r->var_config_N;l++){
                        r->var_config[l].sim = r;
                    }
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
                reading_fields = 0;
                break;
            case REB_BINARY_FIELD_TYPE_FUNCTIONPOINTERS:
                {
                    int fpwarn;
                    fread(&fpwarn, field.size,1,inf);
                    if (fpwarn){
                        *warnings |= REB_INPUT_BINARY_WARNING_POINTERS;
                    }
                }
                break;
            default:
                *warnings |= REB_INPUT_BINARY_WARNING_FIELD_UNKOWN;
                fseek(inf,field.size,SEEK_CUR);
                break;
        }
    }
    fclose(inf);
}

struct reb_simulation* reb_create_simulation_from_binary(char* filename){
    enum reb_input_binary_messages warnings = REB_INPUT_BINARY_WARNING_NONE;
    struct reb_simulation* r = reb_create_simulation();
    reb_create_simulation_from_binary_with_messages(r,filename,&warnings);
    if (warnings & REB_INPUT_BINARY_WARNING_VERSION){
        reb_warning(r,"Binary file was saved with a different version of REBOUND. Binary format might have changed.");
    }
    if (warnings & REB_INPUT_BINARY_WARNING_FIELD_UNKOWN){
        reb_warning(r,"Unknown field found in binary file.");
    }
    if (warnings & REB_INPUT_BINARY_WARNING_PARTICLES){
        reb_warning(r,"Binary file might be corrupted. Number of particles found does not match expected number.");
    }
    if (warnings & REB_INPUT_BINARY_WARNING_POINTERS){
        reb_warning(r,"You have to reset function pointers after creating a reb_simulation struct with a binary file.");
    }
    if (warnings & REB_INPUT_BINARY_ERROR_NOFILE){
        reb_error(r,"Cannot read binary file. Check filename and file contents.");
        free(r);
        r = NULL;
    }
    return r;
}

