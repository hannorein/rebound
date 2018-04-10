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

FILE* reb_binary_diff(FILE* f1, FILE* f2){
    if (!f1 || !f2){
        printf("Cannot read binary file.\n");
        return NULL;
    }
    
    long f10 = ftell(f1);
    fseek(f1, 0 , SEEK_END);
    long f1length = ftell(f1)-f10;
    fseek(f1, f10 , SEEK_SET);
    long f20 = ftell(f2);
    fseek(f2, 0 , SEEK_END);
    long f2length = ftell(f2)-f20;
    fseek(f2, f20 , SEEK_SET);

    FILE* diff = fmemopen(NULL,512,"w+");
    if (diff==0){
        printf("fmemopen failed\n");
        return NULL;
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
        printf("Header in binary files are different.\n");
    }
    
    int reading_fields1 = 1;
    int reading_fields2 = 1;
    while(reading_fields1==1 && reading_fields2==1){
        struct reb_binary_field field1;
        fread(&field1,sizeof(struct reb_binary_field),1,f1);
        if (field1.type==REB_BINARY_FIELD_TYPE_END){
            reading_fields1 = 0;
        }
        struct reb_binary_field field2;
        fread(&field2,sizeof(struct reb_binary_field),1,f2);
        if (field2.type==REB_BINARY_FIELD_TYPE_END){
            reading_fields2 = 0;
        }
        // Fields are not in the same order.
        if (field1.type!=field2.type){
            // Will search for element in f2, starting at beginning just past header
            // Note that we ignore all ADDITIONAL fields in f2 that were not present in f1 
            fseek(f2, 64, SEEK_SET);
            int notfound = 1; 
            int bytesread = 0;
            do {
                bytesread = fread(&field2,sizeof(struct reb_binary_field),1,f2);
                if (bytesread == 0){
                    notfound = 0;
                    break;
                }
                if(field2.type==REB_BINARY_FIELD_TYPE_END){
                    notfound = 0;
                    break;
                }
            }while(field2.type!=field1.type);
            if (notfound == 1){
                printf("No matching field for %d found in f2.\n",field1.type);
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
        fread(readbuf1, field1.size,1,f1);
        fread(readbuf2, field2.size,1,f2);

        int fields_differ = 0;
        if (field1.size==field2.size){
            if (memcmp(readbuf1,readbuf2,field1.size)!=0){
                fields_differ = 1;
            }
        }else{
            fields_differ = 1;
        }
        if(fields_differ){
            printf("Field %d differs.\n",field1.type);
            fwrite(&field2,sizeof(struct reb_binary_field),1,diff);
            fwrite(readbuf2,field2.size,1,diff);
        }
    }
    free(readbuf1);
    free(readbuf2);
    return diff;
}
