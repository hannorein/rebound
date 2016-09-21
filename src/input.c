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
    

void reb_create_simulation_from_binary_with_messages(struct reb_simulation* r, char* filename, enum reb_input_binary_messages* warnings){
    FILE* inf = fopen(filename,"rb"); 
    
    if (!inf){
        *warnings |= REB_INPUT_BINARY_ERROR_NOFILE;
        return;
    }

    long objects = 0;
    // Input header.
    const char str[] = "REBOUND Binary File. Version: ";
    char readbuf[65], curvbuf[65];
    sprintf(curvbuf,"%s%s",str,reb_version_str);
    for(size_t j=strlen(curvbuf);j<63;j++){
        curvbuf[j] = ' ';
    }
    curvbuf[63] = '\0';
    objects += fread(readbuf,sizeof(char),64,inf);
    if(strcmp(readbuf,curvbuf)!=0){
        *warnings |= REB_INPUT_BINARY_WARNING_VERSION;
    }
    
    reb_reset_temporary_pointers(r);
    r->tree_root = NULL;

    // Read two longs to get size of entire restart file and location of first 
    // particle data.
    
    
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
            CASE(SASEEKFIRST,        &r->simulationarchive_seek_first);
            CASE(SASEEKBLOB,         &r->simulationarchive_seek_blob);
            CASE(SAINTERVAL,         &r->simulationarchive_interval);
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
            case REB_BINARY_FIELD_TYPE_PARTICLES:
                r->particles = malloc(field.size);
                r->allocatedN = field.size/sizeof(struct reb_particle);
                fread(r->particles, field.size,1,inf);
                for (int l=0;l<r->N;l++){
                    r->particles[l].c = NULL;
                    r->particles[l].ap = NULL;
                    r->particles[l].sim = r;
                }
                break;
            case REB_BINARY_FIELD_TYPE_END:
                reading_fields = 0;
                break;
            default:
                *warnings |= REB_INPUT_BINARY_WARNING_FIELD_UNKOWN;
                fseek(inf,field.size,SEEK_CUR);
                break;
        }
    }

    int ri_ias15_allocatedN = r->ri_ias15.allocatedN;
    if(reb_reset_function_pointers(r)){
        *warnings |= REB_INPUT_BINARY_WARNING_POINTERS;
    }

    // Read variational config structures
    if (r->var_config_N>0){
        r->var_config = malloc(sizeof(struct reb_variational_configuration)*r->var_config_N);
        if (r->var_config){
            objects = fread(r->var_config,sizeof(struct reb_variational_configuration),r->var_config_N,inf);
            if (objects==r->var_config_N){
                for (int l=0;l<r->var_config_N;l++){
                    r->var_config[l].sim = r;
                }
            }else{
                *warnings |= REB_INPUT_BINARY_WARNING_VARCONFIG;
            }
        }else{
            *warnings |= REB_INPUT_BINARY_WARNING_VARCONFIG;
        }
    }

    // Read temporary arrays for IAS15 (needed for bit-by-bit reproducability)
    if (ri_ias15_allocatedN && !(*warnings & REB_INPUT_BINARY_WARNING_PARTICLES)){
        int N3 = ri_ias15_allocatedN;
        r->ri_ias15.allocatedN = N3;
        r->ri_ias15.at = malloc(sizeof(double)*N3);
        fread(r->ri_ias15.at,sizeof(double),N3,inf);
        r->ri_ias15.x0 = malloc(sizeof(double)*N3);
        fread(r->ri_ias15.x0,sizeof(double),N3,inf);
        r->ri_ias15.v0 = malloc(sizeof(double)*N3);
        fread(r->ri_ias15.v0,sizeof(double),N3,inf);
        r->ri_ias15.a0 = malloc(sizeof(double)*N3);
        fread(r->ri_ias15.a0,sizeof(double),N3,inf);
        r->ri_ias15.csx = malloc(sizeof(double)*N3);
        fread(r->ri_ias15.csx,sizeof(double),N3,inf);
        r->ri_ias15.csv = malloc(sizeof(double)*N3);
        fread(r->ri_ias15.csv,sizeof(double),N3,inf);
        r->ri_ias15.csa0 = malloc(sizeof(double)*N3);
        fread(r->ri_ias15.csa0,sizeof(double),N3,inf);
        reb_read_dp7(&(r->ri_ias15.g)  ,N3,inf);
        reb_read_dp7(&(r->ri_ias15.b)  ,N3,inf);
        reb_read_dp7(&(r->ri_ias15.csb),N3,inf);
        reb_read_dp7(&(r->ri_ias15.e)  ,N3,inf);
        reb_read_dp7(&(r->ri_ias15.br) ,N3,inf);
        reb_read_dp7(&(r->ri_ias15.er) ,N3,inf);
    }
    
    r->simulationarchive_filename = NULL;
    fclose(inf);
}

struct reb_simulation* reb_create_simulation_from_binary(char* filename){
    enum reb_input_binary_messages warnings = REB_INPUT_BINARY_WARNING_NONE;
    struct reb_simulation* r = malloc(sizeof(struct reb_simulation));
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
    if (warnings & REB_INPUT_BINARY_WARNING_VARCONFIG){
        reb_warning(r,"Binary file might be corrupted. Number of variational config structs found does not match number of variational config structs expected.");
    }
    if (warnings & REB_INPUT_BINARY_WARNING_POINTERS){
        reb_warning(r,"You have to reset function pointers after creating a reb_simulation struct with a binary file.");
    }
    if (warnings & REB_INPUT_BINARY_ERROR_NOFILE){
        reb_error(r,"Cannot read binary file. Check filename and file contents.");
        free(r);
        r = NULL;
    }
    if (warnings & REB_INPUT_BINARY_ERROR_SEEK){
        reb_error(r,"Error occured during recovery attempt.");
        free(r);
        r = NULL;
    }
    return r;
}

