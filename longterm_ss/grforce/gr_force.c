#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <sys/stat.h>
#include "rebound.h"


void gr_force(struct reb_simulation* const r){
    // From REBOUNDx
    const struct reb_particle source = r->particles[0];
    const double C2 = 1.0130251e+08;
    const double prefac1 = 6.*(r->G*source.m)*(r->G*source.m)/C2;
    
	for (int i=1;i<r->N;i++){
        const struct reb_particle p = r->particles[i];
        const double dx = p.x - source.x;
        const double dy = p.y - source.y;
        const double dz = p.z - source.z;
        const double r2 = dx*dx + dy*dy + dz*dz;
        const double prefac = prefac1/(r2*r2);
        r->particles[i].ax -= prefac*dx;
        r->particles[i].ay -= prefac*dy;
        r->particles[i].az -= prefac*dz;
        r->particles[0].ax += p.m/source.m*prefac*dx;
        r->particles[0].ay += p.m/source.m*prefac*dy;
        r->particles[0].az += p.m/source.m*prefac*dz;
    }
}


