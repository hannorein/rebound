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
#include "main.h"
#include "input.h"
#include "communication_mpi.h"

char input_arguments[4096]; // This is a bit of an arbitrary number. Should be dynamic.

void input_append_input_arguments_with_int(const char* argument, int value){
	if (strlen(input_arguments)){
		strcat(input_arguments,"__");
	}
	char addition[2048];
	sprintf(addition,"%s_%d",argument,value);
	strcat(input_arguments,addition);
}

void input_append_input_arguments_with_double(const char* argument, double value){
	if (strlen(input_arguments)){
		strcat(input_arguments,"__");
	}
	char addition[2048];
	sprintf(addition,"%s_%.3e",argument,value);
	strcat(input_arguments,addition);
}

int input_check_restart(int argc, char** argv){
	char filename[1024];
	int restart = 0;
	opterr = 0;
	optind = 1;
  	while (1) {
      		static struct option long_options[] = {
	  		{"restart", required_argument, 0, 'r'},
			{0,0,0,0}
		};

      		/* getopt_long stores the option index here.   */
      		int option_index = 0;
		//				short options. format abc:d::
      		int c = getopt_long (argc, argv, "", long_options, &option_index);

      		/* Detect the end of the options.   */
      		if (c == -1) break;

      		switch (c)
		{
			case 'r':
				restart = 1;
				strcpy(filename, optarg);
				break;
			default:
				break;
		}
  	}
	if (restart==1){
		input_binary(filename);
	}
	return restart;
}

double input_get_double(int argc, char** argv, const char* argument, double _default){
	char* value = input_get_argument(argc,argv,argument);
	if (value){
		input_append_input_arguments_with_double(argument,atof(value));
		return atof(value);
	}
	return _default;
}

int input_get_int(int argc, char** argv, const char* argument, int _default){
	char* value = input_get_argument(argc,argv,argument);
	if (value){
		input_append_input_arguments_with_int(argument,atoi(value));
		return atoi(value);
	}
	return _default;
}


char* input_get_argument(int argc, char** argv, const char* argument){
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

void input_binary(char* filename){
#ifdef MPI
	char filename_mpi[1024];
	sprintf(filename_mpi,"%s_%d",filename,mpi_id);
	FILE* inf = fopen(filename_mpi,"rb"); 
#else // MPI
	FILE* inf = fopen(filename,"rb"); 
#endif // MPI
	long objects = 0;
	int _N;
	objects += fread(&_N,sizeof(int),1,inf);
	objects += fread(&t,sizeof(double),1,inf);
#ifdef MPI
	printf("Found %d particles in file '%s'. ",_N,filename_mpi);
#else // MPI
	printf("Found %d particles in file '%s'. ",_N,filename);
#endif // MPI
	for (int i=0;i<_N;i++){
		struct particle p;
		objects += fread(&p,sizeof(struct particle),1,inf);
		particles_add(p);
	}
	fclose(inf);
	printf("%ld objects read. Restarting at time t=%f\n",objects,t);
}

