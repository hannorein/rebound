/**
 * @file 	simulationarchive.c
 * @brief 	Tools for creating and reading Simulationarchive binary files.
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
#include <string.h>
#include <sys/stat.h>
#include <stdint.h>
#include "particle.h"
#include "rebound.h"
#include "fmemopen.h"
#include "binarydiff.h"
#include "output.h"
#include "tools.h"
#include "input.h"
#include "output.h"
#include "integrator_ias15.h"
#ifdef MPI
#include "communication_mpi.h"
#endif


void reb_simulation_create_from_simulationarchive_with_messages(struct reb_simulation* r, struct reb_simulationarchive* sa, int64_t snapshot, enum reb_simulation_binary_error_codes* warnings){
    FILE* inf = sa->inf;
    if (inf == NULL){
        *warnings |= REB_SIMULATION_BINARY_ERROR_FILENOTOPEN;
        return;
    }
    if (snapshot<0) snapshot += sa->nblobs;
    if (snapshot>=sa->nblobs || snapshot<0){
        *warnings |= REB_SIMULATION_BINARY_ERROR_OUTOFRANGE;
        return;
    }
    
    // load original binary file
    reb_simulation_free_pointers(r);
    memset(r,0,sizeof(struct reb_simulation));
    reb_simulation_init(r);
#ifdef MPI
    reb_communication_mpi_init(r, 0, NULL);
#endif //MPI
    r->simulationarchive_filename = NULL;
    // reb_simulation_create sets simulationarchive_version to 3 by default.
    // This will break reading in old version.
    // Set to old version by default. Will be overwritten if new version was used.
    r->simulationarchive_version = 0;

    fseek(inf, 0, SEEK_SET);
    reb_input_fields(r, inf, warnings);

    // Done?
    if (snapshot==0) return;

    // Read SA snapshot
    if(fseek(inf, sa->offset[snapshot], SEEK_SET)){
        *warnings |= REB_SIMULATION_BINARY_ERROR_SEEK;
        //reb_simulation_free(r);
        return;
    }
    if (r->simulationarchive_version<2){ 
        *warnings |= REB_SIMULATION_BINARY_ERROR_OLD;
        //reb_simulation_free(r);
        return;
    }else{
        // Version 2 or higher
        reb_input_fields(r, inf, warnings);
    }
    return;
}

struct reb_simulation* reb_simulation_create_from_simulationarchive(struct reb_simulationarchive* sa, int64_t snapshot){
    if (sa==NULL) return NULL;
    enum reb_simulation_binary_error_codes warnings = REB_SIMULATION_BINARY_WARNING_NONE;
    struct reb_simulation* r = reb_simulation_create();
    reb_simulation_create_from_simulationarchive_with_messages(r, sa, snapshot, &warnings);
    r = reb_input_process_warnings(r, warnings);
    return r; // might be null if error occured
}

// Old 16 bit offsets. Used only to read old files.
struct reb_simulationarchive_blob16 {
    int32_t index;
    int16_t offset_prev;
    int16_t offset_next;
};

void reb_read_simulationarchive_from_stream_with_messages(struct reb_simulationarchive* sa, struct reb_simulationarchive* sa_index, enum reb_simulation_binary_error_codes* warnings){
    // Assumes sa->inf is set to an open stream
    const int debug = 0;
    if (sa->inf==NULL){
        *warnings |= REB_SIMULATION_BINARY_ERROR_NOFILE;
        return;
    }
    
    // Get version
    fseek(sa->inf, 0, SEEK_SET);  
    struct reb_binary_field field = {0};
    sa->version = 0;
    double t0 = 0;
    sa->reb_version_major = 0;
    sa->reb_version_minor = 0;
    sa->reb_version_patch = 0;
    int uses32bitoffsets = 1; 
    // Cache descriptors
    struct reb_binary_field_descriptor fd_header = reb_binary_field_descriptor_for_name("header");
    struct reb_binary_field_descriptor fd_t = reb_binary_field_descriptor_for_name("t");
    struct reb_binary_field_descriptor fd_sa_version = reb_binary_field_descriptor_for_name("simulationarchive_version");
    struct reb_binary_field_descriptor fd_sa_auto_walltime = reb_binary_field_descriptor_for_name("simulationarchive_auto_walltime");
    struct reb_binary_field_descriptor fd_sa_auto_interval = reb_binary_field_descriptor_for_name("simulationarchive_auto_interval");
    struct reb_binary_field_descriptor fd_sa_auto_step = reb_binary_field_descriptor_for_name("simulationarchive_auto_step");
    struct reb_binary_field_descriptor fd_end = reb_binary_field_descriptor_for_name("end");


    do{
        int didReadField = (int)fread(&field,sizeof(struct reb_binary_field),1,sa->inf);
        if (!didReadField){
            *warnings |= REB_SIMULATION_BINARY_WARNING_CORRUPTFILE;
            break;
        }
        if (field.type == fd_header.type){
            int64_t objects = 0;
            // Input header.
            const int64_t bufsize = 64 - sizeof(struct reb_binary_field);
            char readbuf[64], curvbuf[64];
            const char* header = "REBOUND Binary File. Version: ";
            sprintf(curvbuf,"%s%s",header+sizeof(struct reb_binary_field), reb_version_str);

            objects += fread(readbuf,sizeof(char),bufsize,sa->inf);
            // Finding version_major/version_minor version
            int c1=0, c2=0, c3=0; 
            for (int c=0; c<bufsize; c++){
                if (c2 != 0 && c3 == 0 && readbuf[c] == '.'){
                    c3 = c;
                }
                if (c1 != 0 && c2 == 0 && readbuf[c] == '.'){
                    c2 = c;
                }
                if (c1 == 0 && readbuf[c] == ':'){
                    c1 = c;
                }
            }
            if (c1==0 || c2==0 || c3==0){
                if (debug) printf("SA Error. Cannot determine version.\n");
                *warnings |= REB_SIMULATION_BINARY_WARNING_CORRUPTFILE;
            }else{
                char cpatch[64];
                char cminor[64];
                char cmajor[64];
                strncpy(cpatch, readbuf+c3+1, 3);
                cminor[4] = '\0';
                strncpy(cminor, readbuf+c2+1, c3-c2-1);
                cminor[c3-c2-1] = '\0';
                strncpy(cmajor, readbuf+c1+1, c2-c1-1);
                cmajor[c2-c1-1] = '\0';
                sa->reb_version_patch = atoi(cpatch);
                sa->reb_version_minor = atoi(cminor);
                sa->reb_version_major = atoi(cmajor);
                if (sa->reb_version_major <= 3 && sa->reb_version_minor < 18){
                    uses32bitoffsets = 0; // fallback to 16 bit 
                }
            }
            if (objects < 1){
                *warnings |= REB_SIMULATION_BINARY_WARNING_CORRUPTFILE;
            }else{
                // Note: following compares version, but ignores githash.
                if(strncmp(readbuf,curvbuf,bufsize)!=0){
                    *warnings |= REB_SIMULATION_BINARY_WARNING_VERSION;
                }
            }

        }else if (field.type == fd_t.type){
            fread(&t0, field.size, 1, sa->inf);
        }else if (field.type == fd_sa_version.type){
            fread(&(sa->version), field.size, 1, sa->inf);
        }else if (field.type == fd_sa_auto_walltime.type){
            fread(&(sa->auto_walltime), field.size, 1, sa->inf);
        }else if (field.type == fd_sa_auto_interval.type){
            fread(&(sa->auto_interval), field.size, 1, sa->inf);
        }else if (field.type == fd_sa_auto_step.type){
            fread(&(sa->auto_step), field.size, 1, sa->inf);
        }else{
            fseek(sa->inf,field.size,SEEK_CUR);
        }
    }while(field.type!=fd_end.type);

    // Make index
    if (sa->version<2){
        // Version 1 no longer supported
        free(sa->filename);
        sa->filename = NULL;
        fclose(sa->inf);
        sa->inf = NULL;
        *warnings |= REB_SIMULATION_BINARY_ERROR_OLD;
        return;
    }else{
        // New version
        if (debug) printf("=============\n");
        if (debug) printf("SA Version: 2\n");
        if (sa_index == NULL){ // Need to construct offset index from file.
            int64_t nblobsmax = 1024;
            sa->t = calloc(nblobsmax,sizeof(double));
            sa->offset = calloc(nblobsmax,sizeof(uint64_t));
            fseek(sa->inf, 0, SEEK_SET);  
            sa->nblobs = 0;
            int read_error = 0;
            struct reb_binary_field_descriptor fd_header = reb_binary_field_descriptor_for_name("header");
            struct reb_binary_field_descriptor fd_t = reb_binary_field_descriptor_for_name("t");
            struct reb_binary_field_descriptor fd_end = reb_binary_field_descriptor_for_name("end");
            for(int64_t i=0;i<nblobsmax;i++){
                struct reb_binary_field field = {0};
                sa->offset[i] = ftell(sa->inf);
                int blob_finished = 0;
                do{
                    size_t r1 = fread(&field,sizeof(struct reb_binary_field),1,sa->inf);
                    if (r1==1){
                        if (field.type == fd_header.type){
                            if (debug) printf("SA Field. type=HEADER\n");
                            int s1 = fseek(sa->inf,64 - sizeof(struct reb_binary_field),SEEK_CUR);
                            if (s1){
                                read_error = 1;
                            }
                        }else if (field.type == fd_t.type){
                            size_t r2 = fread(&(sa->t[i]), field.size,1,sa->inf);
                            if (debug) printf("SA Field. type=TIME      value=%.10f\n",sa->t[1]);
                            if (r2!=1){
                                read_error = 1;
                            }
                        }else if (field.type == fd_end.type){
                            if (debug) printf("SA Field. type=END   =====\n");
                            blob_finished = 1;
                        }else{
                            int s2 = fseek(sa->inf,field.size,SEEK_CUR);
                            if (debug) printf("SA Field. type=%-6d    size=%" PRIu64 "\n",field.type,(uint64_t)field.size);
                            if (s2){
                                read_error = 1;
                            }
                        }
                    }else{
                        read_error = 1;
                    }
                }while(blob_finished==0 && read_error==0);
                if (read_error){
                    if (debug) printf("SA Error. Error while reading current blob.\n");
                    // Error during reading. Current snapshot is corrupt.
                    break;
                }
                // Everything looks normal so far. Attempt to read next blob
                struct reb_simulationarchive_blob blob = {0};
                size_t r3;
                if (uses32bitoffsets){
                    r3 = fread(&blob, sizeof(struct reb_simulationarchive_blob), 1, sa->inf);
                }else{
                    // Workaround for versions < 3.18
                    struct reb_simulationarchive_blob16 blob16 = {0};
                    r3 = fread(&blob16, sizeof(struct reb_simulationarchive_blob16), 1, sa->inf);
                    blob.index = blob16.index;
                    blob.offset_prev = blob16.offset_prev;
                    blob.offset_next = blob16.offset_next;
                }
                int next_blob_is_corrupted = 0;
                if (r3!=1){ // Next snapshot is definitly corrupted. 
                            // Assume we have reached the end of the file. 
                            // Won't be able to do checksum. 
                    if (debug) printf("SA Error. Error while reading next blob.\n");
                    next_blob_is_corrupted = 1;
                    *warnings |= REB_SIMULATION_BINARY_WARNING_CORRUPTFILE;
                }
                if (i>0){
                    size_t blobsize;
                    if (uses32bitoffsets){
                        blobsize = sizeof(struct reb_simulationarchive_blob);
                    }else{
                        blobsize = sizeof(struct reb_simulationarchive_blob16);
                    }
                    // Checking the offsets. Acts like a checksum.
                    if (((int64_t)blob.offset_prev )+ ((int64_t)blobsize) != ftell(sa->inf) - ((int64_t)sa->offset[i]) ){
                        // Offsets don't work. Next snapshot is definitly corrupted. Assume current one as well.
                        if (debug) printf("SA Error. Offset mismatch: %lu != %" PRIu64 ".\n",blob.offset_prev + blobsize, (uint64_t)(ftell(sa->inf) - sa->offset[i]) );
                        read_error = 1;
                        break;
                    }
                }
                // All tests passed. Accept current snapshot. Increase blob count.
                sa->nblobs = i+1;
                if (blob.offset_next==0 || next_blob_is_corrupted){
                    // Last blob. 
                    if (debug) printf("SA Reached final blob.\n");
                    break;
                }
                if (i==nblobsmax-1){ // Increase 
                    nblobsmax += 1024;
                    sa->t = realloc(sa->t,sizeof(double)*nblobsmax);
                    sa->offset = realloc(sa->offset,sizeof(uint64_t)*nblobsmax);
                }
            }
            if (read_error){
                if (sa->nblobs>0){
                    *warnings |= REB_SIMULATION_BINARY_WARNING_CORRUPTFILE;
                }else{
                    fclose(sa->inf);
                    sa->inf = NULL;
                    free(sa->filename);
                    sa->filename = NULL;
                    free(sa->t);
                    sa->t = NULL;
                    free(sa->offset);
                    sa->offset = NULL;
                    free(sa);
                    *warnings |= REB_SIMULATION_BINARY_ERROR_SEEK;
                    return;
                }
            }

        }else{ // reuse index from other SA
            // This is an optimzation for loading many large SAs.
            // It assumes the structure of this SA is *exactly* the same as in sa_index.
            // Unexpected behaviour if the shape is not the same.
            sa->nblobs = sa_index->nblobs;
            sa->t = malloc(sizeof(double)*sa->nblobs);
            sa->offset = malloc(sizeof(uint64_t)*sa->nblobs);
            fseek(sa->inf, 0, SEEK_SET);
            // No need to read the large file, just copying the index.
            memcpy(sa->offset, sa_index->offset, sizeof(uint64_t)*sa->nblobs);
            memcpy(sa->t, sa_index->t, sizeof(double)*sa->nblobs);
        }
    }
}

void reb_simulationarchive_create_from_file_with_messages(struct reb_simulationarchive* sa, const char* filename,  struct reb_simulationarchive* sa_index, enum reb_simulation_binary_error_codes* warnings){
    // Somewhat complicated calls for backwards compatability.
#ifdef MPI
    int initialized;
    MPI_Initialized(&initialized);
    if (!initialized){
        int argc = 0;
        char** argv = NULL;
        MPI_Init(&argc, &argv);
    }
    int mpi_id=0;
    MPI_Comm_rank(MPI_COMM_WORLD,&mpi_id);
    char filename_mpi[1024];
    sprintf(filename_mpi,"%s_%d",filename,mpi_id);
    sa->inf = fopen(filename_mpi,"rb");
#else // MPI
    sa->inf = fopen(filename,"rb");
#endif // MPI
    sa->filename = malloc(strlen(filename)+1);
    strcpy(sa->filename,filename);
    reb_read_simulationarchive_from_stream_with_messages(sa, sa_index, warnings);
}

void reb_simulationarchive_init_from_buffer_with_messages(struct reb_simulationarchive* sa, char* buf, size_t size, struct reb_simulationarchive* sa_index, enum reb_simulation_binary_error_codes* warnings){
    // Somewhat complicated calls for backwards compatability.
    sa->inf = reb_fmemopen(buf,size,"rb");
    sa->filename = NULL;
    reb_read_simulationarchive_from_stream_with_messages(sa, sa_index, warnings);
}

struct reb_simulationarchive* reb_simulationarchive_create_from_file(const char* filename){
    struct reb_simulationarchive* sa = malloc(sizeof(struct reb_simulationarchive));
    enum reb_simulation_binary_error_codes warnings = REB_SIMULATION_BINARY_WARNING_NONE;
    reb_simulationarchive_create_from_file_with_messages(sa, filename, NULL, &warnings);
    if (warnings & REB_SIMULATION_BINARY_ERROR_NOFILE){
        // Don't output an error if file does not exist, just return NULL.
        free(sa);
        sa = NULL;
    }else{
        reb_input_process_warnings(NULL, warnings);
    }
    return sa;
}

void reb_simulationarchive_free(struct reb_simulationarchive* sa){
    reb_simulationarchive_free_pointers(sa);
    free(sa);
}

void reb_simulationarchive_free_pointers(struct reb_simulationarchive* sa){
    if (sa==NULL) return;
    if (sa->inf){
        fclose(sa->inf);
    }
    free(sa->filename);
    free(sa->t);
    free(sa->offset);
}
    
void reb_simulationarchive_heartbeat(struct reb_simulation* const r){
    if (r->simulationarchive_filename!=NULL){
        int modes = 0;
        if (r->simulationarchive_auto_interval!=0) modes++;
        if (r->simulationarchive_auto_walltime!=0.) modes++;
        if (r->simulationarchive_auto_step!=0) modes++;
        if (modes>1){
            reb_simulation_error(r,"Only use one of simulationarchive_auto_interval, simulationarchive_auto_walltime, or simulationarchive_auto_step");
        }
        if (r->simulationarchive_auto_interval!=0.){
            const double sign = r->dt>0.?1.:-1;
            if (sign*r->simulationarchive_next <= sign*r->t){
                r->simulationarchive_next += sign*r->simulationarchive_auto_interval;
                //Snap
                reb_simulation_save_to_file(r, NULL);
            }
        }
        if (r->simulationarchive_auto_step!=0.){
            if (r->simulationarchive_next_step <= r->steps_done){
                r->simulationarchive_next_step += r->simulationarchive_auto_step;
                //Snap
                reb_simulation_save_to_file(r, NULL);
            }
        }
        if (r->simulationarchive_auto_walltime!=0.){
            if (r->simulationarchive_next <= r->walltime){
                r->simulationarchive_next += r->simulationarchive_auto_walltime;
                //Snap
                reb_simulation_save_to_file(r, NULL);
            }
        } 
    }
}

void reb_simulation_save_to_file(struct reb_simulation* const r, const char* filename){
    if (r->simulationarchive_version<3){
        reb_simulation_error(r, "Writing Simulationarchives with a version < 3 is no longer supported.\n");
        return;
    }
    if (filename==NULL) filename = r->simulationarchive_filename;
    struct stat buffer;
#ifdef MPI
    char filename_mpi[1024];
    sprintf(filename_mpi,"%s_%d",filename,r->mpi_id);
    if (stat(filename_mpi, &buffer) < 0){
#else // MPI
    if (stat(filename, &buffer) < 0){
#endif // MPI
        // File does not exist. Output binary.
#ifdef MPI
        FILE* of = fopen(filename_mpi,"wb"); 
#else // MPI
        FILE* of = fopen(filename,"wb"); 
#endif // MPI
        if (of==NULL){
            reb_simulation_error(r, "Can not open file.");
            return;
        }
        char* bufp;
        size_t sizep;
        reb_simulation_save_to_stream(r, &bufp,&sizep);
        fwrite(bufp,sizep,1,of);
        free(bufp);
        fclose(of);
    }else{
        // File exists, append snapshot.
        struct reb_binary_field_descriptor fd_end = reb_binary_field_descriptor_for_name("end");
        // Create buffer containing original binary file
#ifdef MPI
        char filename_mpi[1024];
        sprintf(filename_mpi,"%s_%d",filename,r->mpi_id);
        FILE* of = fopen(filename_mpi,"r+b");
#else // MPI
        FILE* of = fopen(filename,"r+b");
#endif // MPI
        fseek(of, 64, SEEK_SET); // Header
        struct reb_binary_field field = {0};
        struct reb_simulationarchive_blob blob = {0};
        int bytesread;
        do{
            bytesread = (int)fread(&field,sizeof(struct reb_binary_field),1,of);
            fseek(of, field.size, SEEK_CUR);
        }while(field.type!=fd_end.type && bytesread);
        int64_t size_old = ftell(of);
        if (bytesread!=1){
            reb_simulation_warning(r, "Simulationarchive appears to be corrupted. A recovery attempt has failed. No snapshot has been saved.\n");
            return;
        }

        bytesread = (int)fread(&blob,sizeof(struct reb_simulationarchive_blob),1,of);
        if (bytesread!=1){
            reb_simulation_warning(r, "Simulationarchive appears to be corrupted. A recovery attempt has failed. No snapshot has been saved.\n");
            return;
        }
        int archive_contains_more_than_one_blob = 0;
        if (blob.offset_next>0){
            archive_contains_more_than_one_blob = 1;
        }


        char* buf_old = malloc(size_old);
        fseek(of, 0, SEEK_SET);  
        fread(buf_old, size_old,1,of);

        // Create buffer containing current binary file
        char* buf_new;
        size_t size_new;
        reb_simulation_save_to_stream(r, &buf_new, &size_new);

        // Create buffer containing diff
        char* buf_diff;
        size_t size_diff;
        reb_binary_diff(buf_old, size_old, buf_new, size_new, &buf_diff, &size_diff, 0);

        int file_corrupt = 0;
        int seek_ok = fseek(of, -sizeof(struct reb_simulationarchive_blob), SEEK_END);
        int blobs_read = (int)fread(&blob, sizeof(struct reb_simulationarchive_blob), 1, of);
        if (seek_ok !=0 || blobs_read != 1){ // cannot read blob
            file_corrupt = 1;
        }
        if ( (archive_contains_more_than_one_blob && blob.offset_prev <=0) || blob.offset_next != 0){ // blob contains unexpected data. Note: First blob is all zeros.
            file_corrupt = 1;
        }
        if (file_corrupt==0 && archive_contains_more_than_one_blob ){
            // Check if last two blobs are consistent.
            seek_ok = fseek(of, - sizeof(struct reb_simulationarchive_blob) - sizeof(struct reb_binary_field), SEEK_CUR);  
            bytesread = (int)fread(&field, sizeof(struct reb_binary_field), 1, of);
            if (seek_ok!=0 || bytesread!=1){
                file_corrupt = 1;
            }
            if (field.type != fd_end.type || field.size !=0){
                // expected an END field
                file_corrupt = 1;
            }
            seek_ok = fseek(of, -blob.offset_prev - sizeof(struct reb_simulationarchive_blob), SEEK_CUR);  
            struct reb_simulationarchive_blob blob2 = {0};
            blobs_read = (int)fread(&blob2, sizeof(struct reb_simulationarchive_blob), 1, of);
            if (seek_ok!=0 || blobs_read!=1 || blob2.offset_next != blob.offset_prev){
                file_corrupt = 1;
            }
        }

        if (file_corrupt){
            // Find last valid snapshot to allow for restarting and appending to archives where last snapshot was cut off
            reb_simulation_warning(r, "Simulationarchive appears to be corrupted. REBOUND will attempt to fix it before appending more snapshots.\n");
            int seek_ok;
            seek_ok = fseek(of, size_old, SEEK_SET);
            int64_t last_blob = size_old + sizeof(struct reb_simulationarchive_blob);
            do
            {
                seek_ok = fseek(of, -sizeof(struct reb_binary_field), SEEK_CUR);
                if (seek_ok != 0){
                    break;
                }
                bytesread = (int)fread(&field, sizeof(struct reb_binary_field), 1, of);
                if (bytesread != 1 || field.type != fd_end.type){ // could be EOF or corrupt snapshot
                    break;
                }
                bytesread = (int)fread(&blob, sizeof(struct reb_simulationarchive_blob), 1, of);
                if (bytesread != 1){
                    break;
                }
                last_blob = ftell(of);
                if (blob.offset_next>0){
                    seek_ok = fseek(of, blob.offset_next, SEEK_CUR);
                }else{
                    break;
                }
                if (seek_ok != 0){
                    break;
                }
            } while(1);

            // To append diff, seek to last valid location (=EOF if all snapshots valid)
            fseek(of, last_blob, SEEK_SET);
        }else{
            // File is not corrupt. Start at end to save time.
            fseek(of, 0, SEEK_END);  
        }

        // Update blob info and Write diff to binary file
        fseek(of, -sizeof(struct reb_simulationarchive_blob), SEEK_CUR);  
        fread(&blob, sizeof(struct reb_simulationarchive_blob), 1, of);
        blob.offset_next = (int32_t)size_diff+sizeof(struct reb_binary_field);
        fseek(of, -sizeof(struct reb_simulationarchive_blob), SEEK_CUR);  
        fwrite(&blob, sizeof(struct reb_simulationarchive_blob), 1, of);
        fwrite(buf_diff, size_diff, 1, of); 
        field.type = fd_end.type;
        field.size = 0;
        fwrite(&field,sizeof(struct reb_binary_field), 1, of);
        blob.index++;
        blob.offset_prev = blob.offset_next;
        blob.offset_next = 0;
        fwrite(&blob, sizeof(struct reb_simulationarchive_blob), 1, of);

        fclose(of);
        free(buf_new);
        free(buf_old);
        free(buf_diff);
    }
}

static int _reb_simulationarchive_automate_set_filename(struct reb_simulation* const r, const char* filename){
    if (r==NULL) return -1;
    if (filename==NULL){
        reb_simulation_error(r, "Filename missing.");
        return -1;
    }
    struct stat buffer;
#ifdef MPI
    char filename_mpi[1024];
    sprintf(filename_mpi,"%s_%d",filename,r->mpi_id);
    if (stat(filename_mpi, &buffer) == 0){
#else // MPI
    if (stat(filename, &buffer) == 0){
#endif // MPI
        reb_simulation_warning(r, "File in use for Simulationarchive already exists. Snapshots will be appended.");
    }
    free(r->simulationarchive_filename);
    r->simulationarchive_filename = malloc((strlen(filename)+1)*sizeof(char));
    strcpy(r->simulationarchive_filename, filename);
    return 0;
}

void reb_simulation_save_to_file_interval(struct reb_simulation* const r, const char* filename, double interval){
    if(_reb_simulationarchive_automate_set_filename(r,filename)<0) return;
    if(r->simulationarchive_auto_interval != interval){
        // Only update simulationarchive_next if interval changed. 
        // This ensures that interrupted simulations will continue
        // after being restarted from a simulationarchive
        r->simulationarchive_auto_interval = interval;
        r->simulationarchive_next = r->t;
    }
}

void reb_simulation_save_to_file_walltime(struct reb_simulation* const r, const char* filename, double walltime){
    if(_reb_simulationarchive_automate_set_filename(r,filename)<0) return;
    // Note that this will create two snapshots if restarted.
    r->simulationarchive_auto_walltime = walltime;
    r->simulationarchive_next = r->walltime;
}

void reb_simulation_save_to_file_step(struct reb_simulation* const r, const char* filename, uint64_t step){
    if(_reb_simulationarchive_automate_set_filename(r,filename)<0) return;
    if(r->simulationarchive_auto_step != step){
        // Only update simulationarchive_next if interval changed. 
        // This ensures that interrupted simulations will continue
        // after being restarted from a simulationarchive
        r->simulationarchive_auto_step = step;
        r->simulationarchive_next_step = r->steps_done;
    }
}
