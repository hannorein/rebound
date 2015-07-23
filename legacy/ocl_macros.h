#if !defined OCL_MACROS_H
#define OCL_MACROS_H

#define LOG_OCL_ERROR(x, STRING )  if(x!=CL_SUCCESS) {printf( "\nLine No: %d ", __LINE__ ); printf(STRING); printf("\n    Error= %d\n",x); exit(-1); }

#define LOG_OCL_COMPILER_ERROR(PROGRAM, DEVICE)                                          \
        {                                                                                \
            cl_int logStatus;                                                            \
            char * buildLog = NULL;                                                      \
            size_t buildLogSize = 0;                                                     \
            logStatus = clGetProgramBuildInfo(PROGRAM,                                   \
                                              DEVICE,                                    \
                                              CL_PROGRAM_BUILD_LOG,                      \
                                              buildLogSize,                              \
                                              buildLog,                                  \
                                              &buildLogSize);                            \
            if(logStatus != CL_SUCCESS)                                                  \
            {                                                                            \
                printf( "Error # %d logStatus", logStatus );                             \
                printf( ":: clGetProgramBuildInfo<CL_PROGRAM_BUILD_LOG> failed.");       \
                exit(1);                                                                 \
            }                                                                            \
                                                                                         \
            buildLog = (char*)malloc(buildLogSize);                                      \
            if(buildLog == NULL)                                                         \
            {                                                                            \
                printf("Failed to allocate host memory. (buildLog)\n");                  \
                exit(1);                                                                 \
            }                                                                            \
            memset(buildLog, 0, buildLogSize);                                           \
                                                                                         \
            logStatus = clGetProgramBuildInfo(PROGRAM,                                   \
                                              DEVICE,                                    \
                                              CL_PROGRAM_BUILD_LOG,                      \
                                              buildLogSize,                              \
                                              buildLog,                                  \
                                              NULL);                                     \
            if(logStatus != CL_SUCCESS)                                                  \
            {                                                                            \
                printf( "Error # %d logStatus ", logStatus);                             \
                printf( ":: clGetProgramBuildInfo<CL_PROGRAM_BUILD_LOG> failed.");       \
                exit(1);                                                                 \
            }                                                                            \
                                                                                         \
            printf(" \n\t\t\tBUILD LOG\n");                                              \
            printf(" ************************************************\n");               \
            printf("%s",buildLog);                                                            \
            printf(" ************************************************\n");               \
            free(buildLog);                                                              \
            exit(1);                                                              \
        } 

/* Get platform information and set up the Platform for the defined vendor*/                                                            
#define OCL_CREATE_PLATFORMS( PLATFORM )                                      \
    cl_uint     num_platforms;                                                        \
    if ((clGetPlatformIDs(0, NULL, &num_platforms)) == CL_SUCCESS)                    \
    {                                                                                 \
        PLATFORM = (cl_platform_id *)malloc(sizeof(cl_platform_id)*num_platforms);    \
        if(clGetPlatformIDs(num_platforms, PLATFORM, NULL) != CL_SUCCESS)             \
        {                                                                             \
            free(PLATFORM);                                                           \
            exit(-1);                                                                 \
        }                                                                             \
    }                                                                                 

/*Release the Allocated Platforms*/
#define OCL_RELEASE_PLATFORMS( PLATFORM )                                             \
    free(PLATFORM);

#define OCL_CREATE_DEVICE( PLATFORM, DEVICE_TYPE, DEVICES )                                 \
    cl_uint     num_devices;                                                                \
    if (clGetDeviceIDs( PLATFORM, DEVICE_TYPE, 0,                                           \
            NULL, &num_devices) == CL_SUCCESS)                                              \
    {                                                                                       \
        DEVICES = (cl_device_id *)malloc(sizeof(cl_device_id)*num_devices);                 \
        if (clGetDeviceIDs( PLATFORM, DEVICE_TYPE, num_devices,                             \
            DEVICES, NULL) != CL_SUCCESS)                                                   \
        {                                                                                   \
            free(DEVICES);                                                                  \
            exit(-1);                                                                       \
        }                                                                                   \
    }
    
/*Release the Allocated Device*/
#define OCL_RELEASE_DEVICES( DEVICES )                                                 \
    free(DEVICES);
    
#endif
