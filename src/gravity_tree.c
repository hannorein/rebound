#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "main.h"
#include "tree.h"
#include "boundaries.h"

void calculate_forces(){
	if (root==NULL){
		init_tree();
		root = update_tree(root);
	} else {
		check_boundaries();
		root = update_tree(root);
	}
	printf("Number of particles is %d\n",-root->pt);
	for (int i=0; i<N; i++){
		particles[i].ax = 0; 
		particles[i].ay = 0; 
		particles[i].az = 0; 
	}
	// Summing over all Ghost Boxes
	for (int gbx=-nghostx; gbx<=nghostx; gbx++){
	for (int gby=-nghosty; gby<=nghosty; gby++){
	for (int gbz=-nghostz; gbz<=nghostz; gbz++){
		struct ghostbox gb = get_ghostbox(gbx,gby,gbz);
		// Summing over all particle pairs
		for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){
			if (i==j) continue;
			double dx = (gb.shiftx+particles[i].x) - particles[j].x;
			double dy = (gb.shifty+particles[i].y) - particles[j].y;
			double dz = (gb.shiftz+particles[i].z) - particles[j].z;
			double r = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
			double prefact = -G/(r*r*r)*particles[j].m;
			particles[i].ax += prefact*dx; 
			particles[i].ay += prefact*dy; 
			particles[i].az += prefact*dz; 
		}
		}
	}
	}
	}
}
