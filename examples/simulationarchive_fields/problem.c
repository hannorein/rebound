/**
 * Simulationarchive Fields
 *
 * This program outputs all fields of a Simulationarchive
 * in human readable form. Note that the output can be extensive
 * for large files. This is only intended for debugging and to
 * illustrate the nature of the binary file system.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rebound.h"

#define ifprintf(...) if (blob_index==blob_requested || blob_requested == -1) {printf(__VA_ARGS__);}
#define CASE(DT) {\
    case DT:\
            ifprintf("\tdtype:  " #DT "\n");\
    break;\
}

// Old 16 bit offsets. Used only to read old files.
struct reb_simulationarchive_blob16 {
    int32_t index;
    int16_t offset_prev;
    int16_t offset_next;
};

void print_particle(struct reb_particle v, char* padding){
    printf("%sx=%.16g y=%.16g z=%.16g\n", padding, v.x, v.y, v.z);
    printf("%svx=%.16g vy=%.16g vz=%.16g\n", padding, v.vx, v.vy, v.vz);
    printf("%sax=%.16g ay=%.16g az=%.16g\n", padding, v.ax, v.ay, v.az);
    printf("%sm=%.16g r=%.16g\n", padding, v.m, v.r);
    printf("%slast_collision=%.16g hash=%d\n", padding, v.last_collision, v.hash);
}

int main(int argc, char* argv[]) {
    if (argc<2 || argc > 3){
        printf("Usage: rebound filename [snapshot]\n");
        printf("If snapshot is not given, then all snapshots are dumped.\n");
        return 1;
    }
    int32_t blob_index = 0;
    int32_t blob_requested = -1;
    if (argc==3){
        blob_requested = atoi(argv[2]);
    }

    FILE* sa = fopen(argv[1], "rb");
    if (!sa){
        printf("Error opening file \"%s\"\n", argv[1]);
        return 1;
    }


    int uses32bitoffsets = 1; 
    struct reb_binary_field field = {0};
    struct reb_simulationarchive_blob blob = {0};
    struct reb_binary_field_descriptor fd_header = reb_binary_field_descriptor_for_name("header");
    struct reb_binary_field_descriptor fd_particles = reb_binary_field_descriptor_for_name("particles");
    struct reb_binary_field_descriptor fd_ri_whfast_p_jh = reb_binary_field_descriptor_for_name("ri_whfast.p_jh");
    struct reb_binary_field_descriptor fd_end = reb_binary_field_descriptor_for_name("end");


    printf("==== START OF FILE ====\n");
    do{
        ifprintf("==== START OF BLOB[%d] ====\n", blob_index);
        do{
            int didReadField = (int)fread(&field,sizeof(struct reb_binary_field),1,sa);
            if (!didReadField){
                printf("ERROR. Unable to read from file.\n");
                return 1;
            }
            struct reb_binary_field_descriptor fd = reb_binary_field_descriptor_for_type(field.type);
            if (field.type == fd_end.type){
                ifprintf("FIELD\n");
                ifprintf("\tname:   %s\n", fd.name);
                ifprintf("\ttype:   %d\n", fd.type);
                ifprintf("\tsize:   %llu bytes\n", field.size);
            }else if (field.type == fd_header.type){
                // Input header.
                const int64_t bufsize = 64 - sizeof(struct reb_binary_field);
                char readbuf[64];
                fread(readbuf,sizeof(char),bufsize,sa);
                printf("HEADER\n");
                // Finding version_major/version_minor version
                int c1=0, c2=0, c3=0; 
                for (int c=0; c<bufsize; c++){
                    if (c2 != 0 && c3 == 0 && readbuf[c] == '.'){ c3 = c; }
                    if (c1 != 0 && c2 == 0 && readbuf[c] == '.'){ c2 = c; }
                    if (c1 == 0 && readbuf[c] == ':'){ c1 = c; }
                }
                if (c1==0 || c2==0 || c3==0){
                    printf("ERROR. Cannot determine version.\n");
                    return 1;
                }
                char cpatch[64];
                char cminor[64];
                char cmajor[64];
                strncpy(cpatch, readbuf+c3+1, 3);
                cminor[4] = '\0';
                strncpy(cminor, readbuf+c2+1, c3-c2-1);
                cminor[c3-c2-1] = '\0';
                strncpy(cmajor, readbuf+c1+1, c2-c1-1);
                cmajor[c2-c1-1] = '\0';
                int reb_version_patch = atoi(cpatch);
                int reb_version_minor = atoi(cminor);
                int reb_version_major = atoi(cmajor);
                printf("\tversion: %d.%d.%d\n", reb_version_major, reb_version_minor, reb_version_patch);
                if (reb_version_major <= 3 && reb_version_minor < 18){
                    uses32bitoffsets = 0; // fallback to 16 bit 
                }
                if (uses32bitoffsets){
                    printf("\toffsets: 32 bit\n");
                }else{
                    printf("\toffsets: 16 bit (legacy)\n");
                }
            }else{
                ifprintf("FIELD\n");
                ifprintf("\tname:   %s\n", fd.name);
                ifprintf("\ttype:   %d\n", fd.type);
                ifprintf("\tsize:   %llu bytes\n", field.size);
                switch (fd.dtype){
                    CASE(REB_DOUBLE);
                    CASE(REB_INT);
                    CASE(REB_UINT);
                    CASE(REB_UINT32);
                    CASE(REB_INT64);
                    CASE(REB_UINT64);
                    CASE(REB_VEC3D);
                    CASE(REB_PARTICLE);
                    CASE(REB_POINTER);
                    CASE(REB_POINTER_ALIGNED);
                    CASE(REB_DP7);
                    CASE(REB_OTHER);
                    CASE(REB_PARTICLE4);
                    CASE(REB_POINTER_FIXED_SIZE);
                default:
                    ifprintf("\tdtype:  <other>\n");
                    break;
                }

                switch (fd.dtype){
                    case REB_DOUBLE:
                        {
                            double v;
                            fread(&v,field.size,1,sa);
                            ifprintf("\tvalue:  %.16g\n",v);
                        }
                        break;
                    case REB_INT:
                        {
                            int v;
                            fread(&v,field.size,1,sa);
                            ifprintf("\tvalue:  %d\n",v);
                        }
                        break;
                    case REB_UINT:
                        {
                            unsigned int v;
                            fread(&v,field.size,1,sa);
                            ifprintf("\tvalue:  %d\n",v);
                        }
                        break;
                    case REB_UINT32:
                        {
                            uint32_t v;
                            fread(&v,field.size,1,sa);
                            ifprintf("\tvalue:  %d\n",v);
                        }
                        break;
                    case REB_INT64:
                        {
                            int64_t v;
                            fread(&v,field.size,1,sa);
                            ifprintf("\tvalue:  %lld\n",v);
                        }
                        break;
                    case REB_UINT64:
                        {
                            uint64_t v;
                            fread(&v,field.size,1,sa);
                            ifprintf("\tvalue:  %llu\n",v);
                        }
                        break;
                    case REB_VEC3D:
                        {
                            struct reb_vec3d v;
                            fread(&v,field.size,1,sa);
                            ifprintf("\tvalue:  x=%.16g y=%.16g z=%.16g\n",v.x, v.y, v.z);
                        }
                        break;
                    case REB_PARTICLE:
                        {
                            struct reb_particle v;
                            fread(&v,field.size,1,sa);
                            print_particle(v, "value:\t");
                            print_particle(v, "\t\t");
                        }
                        break;
                    case REB_POINTER:
                        if (field.type == fd_particles.type || field.type == fd_ri_whfast_p_jh.type){
                            ifprintf("\tvalue:  \n");
                            int N = field.size/sizeof(struct reb_particle);
                            struct reb_particle* vp = malloc(field.size);
                            if (blob_index==blob_requested || blob_requested == -1) {
                                fread(vp,field.size,1,sa);
                                for (int i=0; i<N; i++){
                                    printf("\t\tPARTICLE[%d]\n", i);
                                    print_particle(vp[i],"\t\t\t");
                                }
                                free(vp);
                            }else{
                                fseek(sa,field.size,SEEK_CUR);
                            }
                            break;
                        }else{
                            ifprintf("\tvalue:  <not shown>\n");
                            fseek(sa,field.size,SEEK_CUR);
                        }
                        break;
                    default:
                        ifprintf("\tvalue:  <not shown>\n");
                        fseek(sa,field.size,SEEK_CUR);
                        break;
                }

            }
        }while(field.type!=fd_end.type);
        int r3=0;
        if (uses32bitoffsets){
            r3 = fread(&blob, sizeof(struct reb_simulationarchive_blob), 1, sa);
        }else{
            // Workaround for versions < 3.18
            struct reb_simulationarchive_blob16 blob16 = {0};
            r3 = fread(&blob16, sizeof(struct reb_simulationarchive_blob16), 1, sa);
            blob.index = blob16.index;
            blob.offset_prev = blob16.offset_prev;
            blob.offset_next = blob16.offset_next;
        }
        if (!r3){
            printf("ERROR. Unable to read next blob from file.\n");
            return 1;
        }
        ifprintf("BLOB[%d]\n", blob_index);
        ifprintf("\tindex:        %d\n", blob.index);
        ifprintf("\toffset_prev:  %d bytes\n", blob.offset_prev);
        ifprintf("\toffset_next:  %d bytes\n", blob.offset_next);
        blob_index++;
    }while(blob.offset_next!=0);
    printf("==== END OF FILE ====\n");

    fclose(sa);
    return 0; 
}
