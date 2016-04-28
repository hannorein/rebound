/**
 * @file 	input.c
 * @brief 	Parse command line options and read retart files.
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
		//				short options. format abc:d::
      		int c = getopt_long (argc, argv, "", long_options, &option_index);

      		/* Detect the end of the options.   */
      		if (c == -1) break;

      		switch (c)
		{
			case 'a':
				return optarg;
				break;
			default:
				break;
		}
  	}
	return NULL;
}

struct reb_simulation* reb_create_simulation_from_binary(char* filename){
	reb_warning("You have to reset function pointers after creating a reb_simulation struct with a binary file.");
	struct reb_simulation* r = malloc(sizeof(struct reb_simulation));
#ifdef MPI
	char filename_mpi[1024];
#warning following code not working yet. mpi_id will be random number.
	sprintf(filename_mpi,"%s_%d",filename,r->mpi_id);
	FILE* inf = fopen(filename_mpi,"rb"); 
#else // MPI
	FILE* inf = fopen(filename,"rb"); 
#endif // MPI
	if (inf){
		long objects = 0;
		objects += fread(r,sizeof(struct reb_simulation),1,inf);
		reb_reset_temporary_pointers(r);
		reb_reset_function_pointers(r);
		r->allocatedN = r->N;
		r->tree_root = NULL;
		r->particles = malloc(sizeof(struct reb_particle)*r->N);
		objects += fread(r->particles,sizeof(struct reb_particle),r->N,inf);
        for (int l=0;l<r->N;l++){
            r->particles[l].sim = r;
        }
#ifdef MPI
		printf("Found %d particles in file '%s'. \n",r->N,filename_mpi);
#else // MPI
		printf("Found %d particles in file '%s'. \n",r->N,filename);
#endif // MPI
		if (r->var_config_N){
            r->var_config = malloc(sizeof(struct reb_variational_configuration)*r->var_config_N);
	        objects +=fread(r->var_config,sizeof(struct reb_variational_configuration),r->var_config_N,inf);
            for (int l=0;l<r->var_config_N;l++){
                r->var_config[l].sim = r;
            }
        }
		fclose(inf);
	}else{
		printf("Can not open file '%s'\n.",filename);
        free(r);
		return NULL;
	}
    for(int i=0; i<r->N; i++){
		r->particles[i].sim = r;
	}
	return r;
}

