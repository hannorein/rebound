/**
 * @file 	simulationarchive.c
 * @brief 	Tools for creating and reading Simulation Archive binary files.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2016 Hanno Rein
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
#include <string.h>
#include <sys/time.h>
#include <stdint.h>
#include "particle.h"
#include "rebound.h"
#include "tools.h"



int reb_simulationarchive_load_blob(struct reb_simulation* r, char* filename, long blob){
    if (access(filename, F_OK) == -1) return -1;
    if (!r) return -2;
    
    FILE* fd = fopen(filename,"r");
    int fseekret = 0;
    if (blob<0){
        fseekret = fseek(fd,-r->simulationarchive_seek_blob,SEEK_END);
    }else if(blob==0){
        // load original binary file
        enum reb_input_binary_messages warnings = 0;
        reb_create_simulation_from_binary_with_messages(r,filename,&warnings);
        if (warnings & REB_INPUT_BINARY_ERROR_NOFILE){
            reb_error(r,"Cannot read binary file. Check filename and file contents.");
        }
        return 0;
    }else{
        fseekret = fseek(fd,r->simulationarchive_seek_first + (blob-1)*r->simulationarchive_seek_blob,SEEK_SET);
    }
    if (fseekret){
        fclose(fd);
        return -3;
    }
    fread(&(r->t),sizeof(double),1,fd);
    fread(&(r->simulationarchive_walltime),sizeof(double),1,fd);
    gettimeofday(&r->simulationarchive_time,NULL);
    while (r->simulationarchive_next<=r->t){
        r->simulationarchive_next += r->simulationarchive_interval;
    }
    switch (r->gravity){
        case REB_GRAVITY_BASIC:
        case REB_GRAVITY_NONE:
            break;
        default:
            reb_error(r,"Restart not implemented for this gravity module.");
            break;
    }
    switch (r->integrator){
        case REB_INTEGRATOR_WHFAST:
            {
                // Recreate Jacobi arrrays
                struct reb_particle* ps = r->particles;
                if (r->ri_whfast.safe_mode==0){
                    // If same mode is off, store unsynchronized Jacobi coordinates
                    if (r->ri_whfast.allocated_N<r->N){
                        if (r->ri_whfast.p_j){
                            free(r->ri_whfast.p_j);
                        }
                        if (r->ri_whfast.eta){
                            free(r->ri_whfast.eta);
                        }
                        r->ri_whfast.p_j= malloc(sizeof(struct reb_particle)*r->N);
                        r->ri_whfast.eta= malloc(sizeof(double)*r->N);
                        r->ri_whfast.allocated_N = r->N;
                    }
                    ps = r->ri_whfast.p_j;
                }
                for(int i=0;i<r->N;i++){
                    fread(&(r->particles[i].m),sizeof(double),1,fd);
                    fread(&(ps[i].x),sizeof(double),1,fd);
                    fread(&(ps[i].y),sizeof(double),1,fd);
                    fread(&(ps[i].z),sizeof(double),1,fd);
                    fread(&(ps[i].vx),sizeof(double),1,fd);
                    fread(&(ps[i].vy),sizeof(double),1,fd);
                    fread(&(ps[i].vz),sizeof(double),1,fd);
                }
                if (r->ri_whfast.safe_mode==0){
                    // Assume we are not synchronized
                    r->ri_whfast.is_synchronized=0.;
                    // Recalculate Jacobi masses
                    r->ri_whfast.eta[0] = r->particles[0].m;
                    r->ri_whfast.p_j[0].m = r->particles[0].m;
                    for (unsigned int i=1;i<r->N;i++){
                        r->ri_whfast.eta[i] = r->ri_whfast.eta[i-1] + r->particles[i].m;
                        r->ri_whfast.p_j[i].m = r->particles[i].m;
                    }
                }
            }
            break;
        default:
            reb_error(r,"Restart not implemented for this integrator.");
            break;
    }
    fclose(fd);

    return 0;
}

struct reb_simulation* reb_simulationarchive_restart(char* filename){
    if (access(filename, F_OK) == -1) return NULL;
    struct reb_simulation* r = reb_create_simulation_from_binary(filename);
    if (r){
        reb_simulationarchive_load_blob(r, filename, -1);
    }
    return r;
}

void reb_simulationarchive_heartbeat(struct reb_simulation* const r){
    if (r->t==0){
        // First output
        r->simulationarchive_seek_blob = sizeof(double)*2+sizeof(double)*7*r->N;
        reb_output_binary(r,r->simulationarchive_filename);
        r->simulationarchive_next += r->simulationarchive_interval;
        r->simulationarchive_walltime = 0.;
        gettimeofday(&r->simulationarchive_time,NULL);
    }else{
        // Appending outputs
        if (r->simulationarchive_next <= r->t){
            r->simulationarchive_next += r->simulationarchive_interval;
            
            struct timeval time_now;
            gettimeofday(&time_now,NULL);
            r->simulationarchive_walltime += time_now.tv_sec-r->simulationarchive_time.tv_sec+(time_now.tv_usec-r->simulationarchive_time.tv_usec)/1e6;
            r->simulationarchive_time = time_now;

            FILE* of = fopen(r->simulationarchive_filename,"a");
            fwrite(&(r->t),sizeof(double),1, of);
            fwrite(&(r->simulationarchive_walltime),sizeof(double),1, of);
            switch (r->integrator){
                case REB_INTEGRATOR_WHFAST:
                    {
                        struct reb_particle* ps = r->particles;
                        if (r->ri_whfast.safe_mode==0){
                            ps = r->ri_whfast.p_j;
                        }
                        for(int i=0;i<r->N;i++){
                            fwrite(&(r->particles[i].m),sizeof(double),1,of);
                            fwrite(&(ps[i].x),sizeof(double),1,of);
                            fwrite(&(ps[i].y),sizeof(double),1,of);
                            fwrite(&(ps[i].z),sizeof(double),1,of);
                            fwrite(&(ps[i].vx),sizeof(double),1,of);
                            fwrite(&(ps[i].vy),sizeof(double),1,of);
                            fwrite(&(ps[i].vz),sizeof(double),1,of);
                        }
                    }
                    break;
                default:
                    reb_error(r,"Restart not implemented for this integrator.");
                    break;
            }
            fclose(of);
        }
    }
}
