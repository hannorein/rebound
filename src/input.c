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

static void reb_read_dp7(struct reb_dp7* dp7, const int N3, FILE* inf){
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

struct reb_simulation* reb_create_simulation_from_binary_with_messages(char* filename, enum reb_input_binary_messages* warnings){
    FILE* inf = fopen(filename,"rb"); 
    if (inf){
        struct reb_simulation* r = malloc(sizeof(struct reb_simulation));
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

        // Read main simulation oject.
        objects += fread(r,sizeof(struct reb_simulation),1,inf);
        int ri_ias15_allocatedN = r->ri_ias15.allocatedN;
        if(reb_reset_function_pointers(r)){
            *warnings |= REB_INPUT_BINARY_WARNING_POINTERS;
        }
        reb_reset_temporary_pointers(r);
        r->allocatedN = r->N;
        r->tree_root = NULL;

        // Read particles
        if (r->N>0){
            r->particles = malloc(sizeof(struct reb_particle)*r->N);
            if (r->particles){
                objects = fread(r->particles,sizeof(struct reb_particle),r->N,inf);
                if (objects==r->N){
                    for (int l=0;l<r->N;l++){
                        r->particles[l].c = NULL;
                        r->particles[l].ap = NULL;
                        r->particles[l].sim = r;
                    }
                }else{
                    *warnings |= REB_INPUT_BINARY_WARNING_PARTICLES;
                }
            }else{
                *warnings |= REB_INPUT_BINARY_WARNING_PARTICLES;
            }
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
        fclose(inf);
        
        return r;
    }
    *warnings |= REB_INPUT_BINARY_ERROR_NOFILE;
    return NULL;
}

struct reb_simulation* reb_create_simulation_from_binary(char* filename){
    enum reb_input_binary_messages warnings = REB_INPUT_BINARY_WARNING_NONE;
    struct reb_simulation* r = reb_create_simulation_from_binary_with_messages(filename,&warnings);
    if (warnings & REB_INPUT_BINARY_WARNING_VERSION){
        reb_warning(r,"Binary file was saved with a different version of REBOUND. Binary format might have changed.");
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
    }
    return r;
}

