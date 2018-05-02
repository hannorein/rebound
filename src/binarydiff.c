/**
 * @file    binarydiff.c
 * @brief   Binary diff allows to compare binary snapshots.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section     LICENSE
 * Copyright (c) 2018 Hanno Rein
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
#include "particle.h"
#include "rebound.h"
#include "tools.h"
#include "binarydiff.h"

void reb_binary_diff(FILE* f1, FILE* f2, char** bufp, size_t* sizep){
    if (!f1 || !f2){
        printf("Cannot read binary file.\n");
        return;
    }
    
    long f10 = ftell(f1);
    long f20 = ftell(f2);

    // Create buffer which is large enough (note that files could have different fields)
    FILE* diff = open_memstream(bufp, sizep);//fmemopen(NULL,f1length+f2length,"w+b");
    if (diff==0){
        printf("fmemopen failed\n");
        return;
    }
    

    long objects1 = 0;
    long objects2 = 0;
    // Input header.
    int bufN = 128;
    char* readbuf1 = malloc(sizeof(char)*bufN);
    char* readbuf2 = malloc(sizeof(char)*bufN);
    objects1 += fread(readbuf1,sizeof(char),64,f1);
    objects2 += fread(readbuf2,sizeof(char),64,f2);
    if(strcmp(readbuf1,readbuf2)!=0){
        printf("Header in binary files are different %s.\n", readbuf2);
    }
    
    while(1){
        int bytesread = 0;
        struct reb_binary_field field1;
        bytesread = fread(&field1,sizeof(struct reb_binary_field),1,f1);
        if (bytesread==0){
            break;
        }
        if (field1.type==REB_BINARY_FIELD_TYPE_END){
            break;
        }
        struct reb_binary_field field2;
        bytesread = fread(&field2,sizeof(struct reb_binary_field),1,f2);
        
        // Fields might not be in the same order.
        if (field1.type!=field2.type){
            // Will search for element in f2, starting at beginning just past header
            // Note that we ignore all ADDITIONAL fields in f2 that were not present in f1 
            fseek(f2, f20+64, SEEK_SET);
            int notfound = 0; 
            while(1) {
                bytesread = fread(&field2,sizeof(struct reb_binary_field),1,f2);
                if (bytesread == 0){
                    notfound = 1;
                    break;
                }
                if(field2.type==REB_BINARY_FIELD_TYPE_END){
                    notfound = 1;
                    break;
                }
                if (field2.type==field1.type){
                    break; // found!!
                }else{
                    fseek(f2,field2.size,SEEK_CUR); //skip
                }
            };
            if (notfound == 1){
                // printf("No match found for field %d in f2.\n",field1.type);
                // Output field with size 0
                fseek(f1, field1.size, SEEK_CUR); // For next search
                fseek(f2, f20+64, SEEK_SET); // For next search
                field1.size = 0;
                fwrite(&field1,sizeof(struct reb_binary_field),1,diff);
                continue;
            }
        }
        // Can assume field1.type == field2.type from here on
        if (field1.size>bufN || field2.size>bufN){
            while (field1.size>bufN || field2.size>bufN){
                bufN *= 2;
            }
            readbuf1 = realloc(readbuf1, sizeof(char)*bufN);
            readbuf2 = realloc(readbuf2, sizeof(char)*bufN);
        }
        bytesread = fread(readbuf1, field1.size,1,f1);
        if (bytesread==0 && field1.size!=0) printf("Corrupt binary file f1.\n");
        bytesread = fread(readbuf2, field2.size,1,f2);
        if (bytesread==0 && field2.size!=0) printf("Corrupt binary file f2.\n");

        int fields_differ = 0;
        if (field1.size==field2.size){
            if (memcmp(readbuf1,readbuf2,field1.size)!=0){
                fields_differ = 1;
            }
        }else{
            fields_differ = 1;
        }
        if(fields_differ){
            // printf("Field %d %d differs.\n",field1.type,field2.type);
            fwrite(&field2,sizeof(struct reb_binary_field),1,diff);
            fwrite(readbuf2,field2.size,1,diff);
        }
    }
    // Search for fields which are present in f2 but not in f1
    fseek(f1, f10+64, SEEK_SET); 
    fseek(f2, f20+64, SEEK_SET);
    while(1){
        int bytesread = 0;
        struct reb_binary_field field2;
        bytesread = fread(&field2,sizeof(struct reb_binary_field),1,f2);
        if (bytesread==0){
            break;
        }
        if (field2.type==REB_BINARY_FIELD_TYPE_END){
            break;
        }
        struct reb_binary_field field1;
        bytesread = fread(&field1,sizeof(struct reb_binary_field),1,f1);
        
        if (field1.type==field2.type){
            // Not a new field. Skip.
            fseek(f1, field1.size, SEEK_CUR); // For next search
            fseek(f2, field2.size, SEEK_CUR); // For next search
            continue;
        }
        // Fields might not be in the same order.
        // Will search for element in f1, starting at beginning just past header
        fseek(f1, f10+64, SEEK_SET);
        int notfound = 0; 
        while(1) {
            bytesread = fread(&field1,sizeof(struct reb_binary_field),1,f1);
            if (bytesread == 0){
                notfound = 1;
                break;
            }
            if(field1.type==REB_BINARY_FIELD_TYPE_END){
                notfound = 1;
                break;
            }
            if (field2.type==field1.type){
                notfound = 0;
                break; // found!!
            }else{
                notfound = 1;
                fseek(f1,field1.size,SEEK_CUR); //skip
            }
        };
        if (notfound == 0){
            // Not a new field. Skip.
            fseek(f1, field1.size, SEEK_CUR); // For next search
            fseek(f2, field2.size, SEEK_CUR); // For next search
            continue;
        }
        if (field2.size>bufN){
            while (field2.size>bufN){
                bufN *= 2;
            }
            readbuf1 = realloc(readbuf1, sizeof(char)*bufN);
            readbuf2 = realloc(readbuf2, sizeof(char)*bufN);
        }
        bytesread = fread(readbuf2, field2.size,1,f2);
        if (bytesread==0 && field2.size!=0) printf("Corrupt binary file f2.\n");

        //printf("Field %d new in f2.\n",field2.type);
        fwrite(&field2,sizeof(struct reb_binary_field),1,diff);
        fwrite(readbuf2,field2.size,1,diff);
    }
    free(readbuf1);
    free(readbuf2);
    fclose(diff);
    return;
}
