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

int input_check_restart(int argc, char** argv){
	char filename[1024];
	int restart = 0;
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

void input_binary(char* filename){
#ifdef MPI
	char filename_mpi[1024];
	sprintf("%s_%d",filename,mpi_id);
	FILE* inf = fopen(filename_mpi,"rb"); 
#else // MPI
	FILE* inf = fopen(filename,"rb"); 
#endif // MPI
	long bytes = 0;
	int _N;
	bytes += fread(&_N,sizeof(int),1,inf);
	bytes += fread(&t,sizeof(double),1,inf);
	printf("Found %d particles in file '%s'. ",_N,filename);
	for (int i=0;i<_N;i++){
		struct particle p;
		bytes += fread(&p,sizeof(struct particle),1,inf);
		particles_add(p);
	}
	fclose(inf);
	printf("%ld bytes read. Restarting at time t=%f\n",bytes,t);
}

