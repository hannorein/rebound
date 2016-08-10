#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>
#include <sys/stat.h>
#include "rebound.h"


const long output_N_log = 10000;
double output_log;
double output_next_log = 1e5*2.*M_PI;

void heartbeat(struct reb_simulation* const r);
double energy0;

void append_binary(struct reb_simulation* const r, const char* const filename){
    reb_integrator_synchronize(r);
    FILE* of = fopen(filename,"a");
    fwrite(&(r->t),sizeof(double),1, of);
    double energy = fabs((energy0-reb_tools_energy(r))/energy0);
    fwrite(&energy,sizeof(double),1, of);
    fclose(of);
}

int proc_i;

int main(int argc, char* argv[]) {
    int procs_N = 1;
    if (argc==2){
        procs_N = atoi(argv[1]);    
    }
    for(proc_i=0; proc_i<procs_N-1; proc_i++) {
        int pid = fork();
        if (pid == 0) {
            break;
        }
    }

    // Read initial conditions
    struct reb_simulation* r = reb_create_simulation_from_binary("../ss2-2016-06-18.bin");
    r->particles[2].x += 6.6845871e-12*(double)proc_i; // Shift by 1 m * proc_i
    reb_move_to_com(r);

    // Setup simulation
    // 10 kyrs takes ~1s
    double tmax = 2.*M_PI*1e9*100.; //100 Gyr  
    r->dt = 8./365.25*2.*M_PI;      // 8 days initial timestep
    r->ri_whfast.safe_mode = 0;     // Turn of safe mode. Need to call integrator_synchronize() before outputs.
    r->ri_whfast.corrector = 5;    // Turn on symplectic correctors 
    r->heartbeat = heartbeat;
    r->integrator = REB_INTEGRATOR_WHFAST;
    //r->gravity = REB_GRAVITY_COMPENSATED;

    energy0 = reb_tools_energy(r);
    output_log = exp((log(tmax)-log(output_next_log))/(double)output_N_log);
    reb_integrate(r, tmax);
}

void heartbeat(struct reb_simulation* const r){
    if (output_next_log <= r->t){
        output_next_log *= output_log;
        char filename[1024];
        sprintf(filename, "output_log_%04d.bin",proc_i);
        append_binary(r, filename);
    }
}
