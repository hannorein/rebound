/**
 * @file 	gravity.c
 * @brief 	Direct gravity calculation, O(N^2).
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @details 	This is the crudest implementation of an N-body code
 * which sums up every pair of particles. It is only useful very small 
 * particle numbers (N<~100) as it scales as O(N^2). Note that the MPI
 * implementation is not well tested and only works for very specific
 * problems. This should be resolved in the future. 
 *
 * 
 * @section LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
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
#include "particle.h"
#include "main.h"
#include "boundaries.h"
#include "communication_mpi.h"
#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif
#include "ocl_macros.h"

#ifdef MPI
#warning GRAVITY_OPENCL might not work with MPI for your problem. 
#endif


#define CONFIG_USE_DOUBLE 1
#if CONFIG_USE_DOUBLE

#if defined(cl_khr_fp64)  // Khronos extension available?
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#define DOUBLE_SUPPORT_AVAILABLE
#elif defined(cl_amd_fp64)  // AMD extension available?
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#define DOUBLE_SUPPORT_AVAILABLE
#endif

#endif // CONFIG_USE_DOUBLE

#if defined(DOUBLE_SUPPORT_AVAILABLE)

// double
typedef double real_t;
#define PI 3.14159265358979323846

#else

#warning Using single precission.
// float
typedef float real_t;
#define PI 3.14159265359f

#endif





const char *src_kernel = 
"__kernel    \n"
"void gravity_opencl_kernel(				\n"
"				const float G, 		\n"
"				const float softening,	\n"
"				__global float* r,	\n"
"				__global float* a,	\n"
"				__global float* m)	\n"
"{ 							\n"
"	int i	= get_global_id(0);			\n"
"	int N  	= get_num_groups(0)*get_local_size(0);	\n"
"       a[i*3+0] = 0;					\n"
"       a[i*3+1] = 0;					\n"
"       a[i*3+2] = 0;					\n"
"	for (int j=0;j<N;j++){				\n"
"		if (i!=j){				\n"
"			float dx = r[i*3+0]-r[j*3+0];	\n"
"			float dy = r[i*3+1]-r[j*3+1];	\n"
"			float dz = r[i*3+2]-r[j*3+2];	\n"
"			float r  = sqrt(		\n"
"			  	  dx*dx			\n"
"				+ dy*dy			\n"
"				+ dz*dz			\n"
"				+ softening*softening);	\n"
"			float prefact = -G/(r*r*r) 	\n"
"				* m[j];			\n"
"							\n"
"			a[i*3+0] += prefact*dx;		\n"
"			a[i*3+1] += prefact*dy;		\n"
"			a[i*3+2] += prefact*dz;		\n"
"							\n"
"		}					\n"
"	}						\n"
"							\n"
"							\n"
"							\n"
"}							\n"
"							\n"
"    							\n";


cl_int 			clStatus;
cl_platform_id*		platforms 	= NULL;
cl_uint			num_devices	= 0;
cl_device_id*		device_list	= NULL;
cl_context		context;
cl_command_queue	command_queue;
cl_program		program;
cl_kernel		kernel;
size_t 			wgs;
size_t 			preferred_wgs_multiple;
cl_mem			R_clmem;
cl_mem			A_clmem;
cl_mem			M_clmem;
cl_event 		gravity_event;
real_t* 			r_host;
real_t* 			a_host;
real_t* 			m_host;

void gravity_calculate_acceleration(){
	if (platforms==NULL){
		r_host = malloc(3*N*sizeof(real_t));
		a_host = malloc(3*N*sizeof(real_t));
		m_host = malloc(N*sizeof(real_t));
		OCL_CREATE_PLATFORMS(platforms);
		OCL_CREATE_DEVICE(platforms[0],CL_DEVICE_TYPE_GPU,device_list);
		clGetDeviceIDs(platforms[0],CL_DEVICE_TYPE_DEFAULT,0,NULL,&num_devices);
		device_list = malloc(sizeof(cl_device_id)*num_devices);
		clGetDeviceIDs(platforms[0],CL_DEVICE_TYPE_DEFAULT,num_devices,device_list,NULL);
		
	
		printf("\n\nOpenCL Setup.\n----------\n");
		for (int i=0;i<num_devices;i++){	
			char buffer[10240];
			clGetDeviceInfo(device_list[0], CL_DEVICE_NAME, sizeof(buffer), buffer, NULL);
			printf("(device %d)\t Name:                  \t%s\n",i, buffer);
			clGetDeviceInfo(device_list[0], CL_DEVICE_VENDOR, sizeof(buffer), buffer, NULL);
			printf("(device %d)\t Vendor:                \t%s\n",i, buffer);
			clGetDeviceInfo(device_list[0], CL_DEVICE_VERSION, sizeof(buffer), buffer, NULL);
			printf("(device %d)\t Version:               \t%s\n",i, buffer);
			clGetDeviceInfo(device_list[0], CL_DEVICE_EXTENSIONS, sizeof(buffer), buffer, NULL);
			printf("(device %d)\t Extensions:            \t%s\n",i, buffer);
		}
		cl_context_properties props[3] = {
			CL_CONTEXT_PLATFORM,
			(cl_context_properties)platforms[0],
			0
		};
		context = clCreateContext(NULL, num_devices, device_list, NULL, NULL, &clStatus);
		LOG_OCL_ERROR(clStatus, "clCreateContext failed.");
		command_queue = clCreateCommandQueue(context, device_list[0],0,&clStatus);
		LOG_OCL_ERROR(clStatus, "clCreateCommandQueue failed.");
		


		R_clmem = clCreateBuffer(context, CL_MEM_READ_ONLY, 3*N*sizeof(real_t), NULL, &clStatus);
		M_clmem = clCreateBuffer(context, CL_MEM_READ_ONLY, N*sizeof(real_t), NULL, &clStatus);
		A_clmem = clCreateBuffer(context, CL_MEM_WRITE_ONLY, 3*N*sizeof(real_t), NULL, &clStatus);
		LOG_OCL_ERROR(clStatus, "clEnqueueWriteBuffer failed.");
		program = clCreateProgramWithSource(context, 1, (const char**)&src_kernel, NULL, &clStatus);
		LOG_OCL_ERROR(clStatus, "clCreateProgramWithSource failed.");
		clStatus = clBuildProgram(program, 1, device_list, NULL, NULL, NULL);
		if (clStatus != CL_SUCCESS){
			LOG_OCL_COMPILER_ERROR(program, device_list[0]);
		}
		kernel = clCreateKernel(program,"gravity_opencl_kernel",&clStatus);
		real_t f_G = G;
		real_t f_softening = softening;
		clStatus  = clSetKernelArg(kernel, 0, sizeof(real_t), (void*)&f_G);
		clStatus |= clSetKernelArg(kernel, 1, sizeof(real_t), (void*)&f_softening);
		clStatus |= clSetKernelArg(kernel, 2, sizeof(cl_mem), (void*)&R_clmem);
		clStatus |= clSetKernelArg(kernel, 3, sizeof(cl_mem), (void*)&A_clmem);
		clStatus |= clSetKernelArg(kernel, 4, sizeof(cl_mem), (void*)&M_clmem);
		LOG_OCL_ERROR(clStatus, "clSetKernelArg failed.");


		clGetKernelWorkGroupInfo(kernel, device_list[0], CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(wgs), &wgs, NULL);
		clGetKernelWorkGroupInfo(kernel, device_list[0], CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(preferred_wgs_multiple), &preferred_wgs_multiple, NULL);
		
		
		printf("(kernel)\t Work group size:                    \t%d\n", wgs);
		printf("(kernel)\t Preferred Work group size multiple: \t%d\n", preferred_wgs_multiple);

		printf("----------\n\n");
	}

	// Setting up host memory
	
	for (int i=0;i<N;i++){
		r_host[i*3+0] 	= particles[i].x;
		r_host[i*3+1] 	= particles[i].y;
		r_host[i*3+2] 	= particles[i].z;
		m_host[i] 	= particles[i].m;
	}
	clStatus = clEnqueueWriteBuffer(command_queue, R_clmem, CL_TRUE, 0, 3*N*sizeof(real_t), r_host, 0,NULL,NULL); 
	clStatus = clEnqueueWriteBuffer(command_queue, M_clmem, CL_TRUE, 0, N*sizeof(real_t), m_host, 0,NULL,NULL); 

	size_t global_size 	= N;
	size_t local_size 	= preferred_wgs_multiple;
	
	clStatus = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global_size, &local_size, 0, NULL, &gravity_event);
	LOG_OCL_ERROR(clStatus, "clEnqueueNDRangeKernel failed.");
	
	clStatus = clEnqueueReadBuffer(command_queue, A_clmem, CL_TRUE, 0, 3*N*sizeof(real_t), a_host, 1, &gravity_event, NULL);
	clStatus = clFinish(command_queue);

	for (int i=0; i<N; i++){
		particles[i].ax = a_host[i*3+0]; 
		particles[i].ay = a_host[i*3+1]; 
		particles[i].az = a_host[i*3+2]; 
	}
}
