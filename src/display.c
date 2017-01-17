/**
 * @file    display.c
 * @brief   Realtime OpenGL visualization.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * @details These functions provide real time visualizations
 * using OpenGL. 
 * 
 * @section LICENSE
 * Copyright (c) 2015 Hanno Rein, Shangfei Liu
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
#ifdef OPENGL
#define DEG2RAD (M_PI/180.)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <pthread.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#define GLFW_INCLUDE_NONE
#include "glad.h"
#include <GLFW/glfw3.h>
#include "rebound.h"
#include "tools.h"
#include "particle.h"
#include "boundary.h"
#include "display.h"
#include "output.h"
#include "integrator.h"
#include "simplefont.h"

static void reb_display(GLFWwindow* window);
                
static const char* onscreenhelp[] = { 
                "REBOUND OPENGL mouse and keyboard commands",
                "----------------------------------------------------",
                " To rotate the view, simply drag the simulation",
                " with the mouse. To zoom in, press the shift key ",
                " and then drag the simulation with the mouse.",
                "----------------------------------------------------",
                " h       | Show/hide this page",
                " q       | Quit simulation",
                " (space) | Pause simulation",
                " d       | Pause real-time visualization", 
                "         | (the simulation continues)",
                " s       | Toggle three dimensional spheres ",
                "         | (looks better)/points (draws faster)",
                " g       | Toggle ghost boxes",
                " r       | Reset view. Press multiple times to ",
                "         | change orientation",
                " x/X     | Move to a coordinate system centred ",
                "         | on a particle (note: does not work if", 
                "         | particle array is resorted)",
                " c       | Toggle clear screen after each time-step",
                " m       | Toggle multisampling",
                " w       | Draw orbits as wires",
                " t       | Show/hide logo, time, timestep and number ",
                "         | of particles.",
                "----------------------------------------------------"
};

struct reb_particle_opengl {
    double x,y,z;
    double vx,vy,vz;
    double r;
};

struct reb_orbit_opengl {
    double x,y,z;
    double a, e, f;
    double omega, Omega, inc;
};

static struct reb_quaternion normalize(struct reb_quaternion quat) {
    float L = sqrtf(quat.x * quat.x + quat.y * quat.y + quat.z * quat.z + quat.w * quat.w);
    quat.x /= L; quat.y /= L; quat.z /= L; quat.w /= L;
    return quat;
}
static struct reb_quaternion conjugate(struct reb_quaternion quat) {
    quat.x = -quat.x; quat.y = -quat.y; quat.z = -quat.z;
    return quat;
}
static struct reb_quaternion mult(struct reb_quaternion A, struct reb_quaternion B) {
    struct reb_quaternion C;
    C.x = A.w*B.x + A.x*B.w + A.y*B.z - A.z*B.y;
    C.y = A.w*B.y - A.x*B.z + A.y*B.w + A.z*B.x;
    C.z = A.w*B.z + A.x*B.y - A.y*B.x + A.z*B.w;
    C.w = A.w*B.w - A.x*B.x - A.y*B.y - A.z*B.z;
    return C;
}
static void matscale(float mat[16], float s){
    mat[0] = s; mat[1] = 0.; mat[2] = 0.; mat[3] = 0.; 
    mat[4] = 0.; mat[5] = s; mat[6] = 0.; mat[7] = 0.; 
    mat[8] = 0.; mat[9] = 0.; mat[10] = s; mat[11] = 0.;
    mat[12] = 0.; mat[13] = 0.; mat[14] = 0.; mat[15] = 1.;
}

static void mattranslate(float mat[16], float x, float y, float z){
    mat[0] = 1.; mat[1] = 0.; mat[2] = 0.; mat[3] = x; 
    mat[4] = 0.; mat[5] = 1.; mat[6] = 0.; mat[7] = y; 
    mat[8] = 0.; mat[9] = 0.; mat[10] = 1.; mat[11] = z;
    mat[12] = 0.; mat[13] = 0.; mat[14] = 0.; mat[15] = 1.;
}

static void matortho(float mat[16], float l, float r, float b, float t, float n, float f) {
    mat[0] = 2.f/(r-l); mat[1] = 0.; mat[2] = 0.; mat[3] = -(r+l)/(r-l);
    mat[4] = 0.; mat[5] = 2.f/(t-b); mat[6] = 0.; mat[7] = -(t+b)/(t-b);
    mat[8] = 0.; mat[9] = 0.; mat[10] = -2.f/(f-n); mat[11] = -(f+n)/(f-n);
    mat[12] = 0.; mat[13] = 0.; mat[14] = 0.; mat[15] = 1.f;
}
static void quat2mat(struct reb_quaternion A, float mat[16]){
    float xx = A.x * A.x; float xy = A.x * A.y; float xz = A.x * A.z;
    float xw = A.x * A.w; float yy = A.y * A.y; float yz = A.y * A.z;
    float yw = A.y * A.w; float zz = A.z * A.z; float zw = A.z * A.w;
    mat[0] = 1.-2.*(yy+zz);
    mat[1] =    2.*(xy-zw);
    mat[2] =    2.*(xz+yw);
    mat[4] =    2.*(xy+zw);
    mat[5] = 1.-2.*(xx+zz);
    mat[6] =    2.*(yz-xw);
    mat[8] =    2.*(xz-yw);
    mat[9] =    2.*(yz+xw);
    mat[10]= 1.-2.*(xx+yy);
    mat[3] = mat[7] = mat[11] = mat[12] = mat[13] = mat[14] = 0; mat[15]= 1;
}
static void multvec(struct reb_quaternion A, float B[3], float vecr[3]) {
    float mat[16];
    quat2mat(A,mat);
    vecr[0] = mat[0]*B[0] + mat[1]*B[1] + mat[2]*B[2];
    vecr[1] = mat[4]*B[0] + mat[5]*B[1] + mat[6]*B[2];
    vecr[2] = mat[8]*B[0] + mat[9]*B[1] + mat[10]*B[2];
}

static void matmult(float A[16], float B[16], float C[16]) {
    for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
        C[i+4*j] = 0.;
    for(int k=0;k<4;k++){
        C[i+4*j] += A[k+4*j]*B[i+4*k];
    }}}
}

static int convertLine(const char* in, float* out){
    int j = 0;
    while(in[j]!=0&&in[j]!=10){  // end on new line or \0
        out[j*2+0] = (float)((((int)in[j]))%16);
        out[j*2+1] = (float)((((int)in[j]))/16);
        j++;
    }
    return j;
}

static unsigned int compileShader(int* shader, int type, const char* source){
    GLint status;
    
    *shader = glCreateShader(type);
    glShaderSource(*shader, 1, &source, NULL);
    glCompileShader(*shader);
    
    GLint logLength;
    glGetShaderiv(*shader, GL_INFO_LOG_LENGTH, &logLength);
    if (logLength > 0) {
        GLchar *log = (GLchar *)malloc(logLength);
        glGetShaderInfoLog(*shader, logLength, &logLength, log);
        printf("\n\n%s\n\n",log);
        free(log);
    }
    
    glGetShaderiv(*shader, GL_COMPILE_STATUS, &status);
    if (status == 0) {
        glDeleteShader(*shader);
        return 0;
    }
    
    return 1;
}

static unsigned int linkProgram(int prog){
    GLint status;
    glLinkProgram(prog);
    
    GLint logLength;
    glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &logLength);
    if (logLength > 0) {
        GLchar *log = (GLchar *)malloc(logLength);
        glGetProgramInfoLog(prog, logLength, &logLength, log);
        printf("\n\n%s\n\n",log);
        free(log);
    }
    
    glGetProgramiv(prog, GL_LINK_STATUS, &status);
    if (status == 0) {
        return 0;
    }
    
    return 1;
}

static int loadShader(const char* vert_source, const char* frag_source){
    GLint vertShader, fragShader;
    
    GLuint _program = glCreateProgram();
    
    if(!compileShader(&vertShader, GL_VERTEX_SHADER, vert_source)) {
        printf("Failed to compile vertex shader.\n");
        return -1;
    }
    
    if (!compileShader(&fragShader, GL_FRAGMENT_SHADER, frag_source)){
        printf("Failed to compile fragment shader.\n");
        return -1;
    }
    
    glAttachShader(_program, vertShader);
    glAttachShader(_program, fragShader);
    
    if (!linkProgram(_program)) {
        printf("Failed to link shader.\n");
        return -1;
    }

    if (vertShader) {
        glDetachShader(_program, vertShader);
        glDeleteShader(vertShader);
    }
    if (fragShader) {
        glDetachShader(_program, fragShader);
        glDeleteShader(fragShader);
    }
    return _program;
}

static void reb_glfw_error_callback(int error, const char* description){
        fprintf(stderr, "GLFW Error: %s\n", description);
}

static void reb_display_mouse_button(GLFWwindow* window, int button, int action, int mods){
    struct reb_display_data* data = glfwGetWindowUserPointer(window);
    data->mouse_action = action;
}

static void reb_display_resize(GLFWwindow* window, int x, int y){
    reb_display(window);
}

static void reb_display_cursor(GLFWwindow* window, double x, double y){
    struct reb_display_data* data = glfwGetWindowUserPointer(window);
    int width, height;
    glfwGetWindowSize(window, &width, &height);
    if (data->mouse_action==GLFW_RELEASE){
        // No button pressed
        data->mouse_x = INFINITY;
        data->mouse_y = INFINITY;
        return;
    }
    if (isinf(data->mouse_x)){
        // New drag event
        data->mouse_x = x;
        data->mouse_y = y;
        return;
    }
    if (data->mouse_action==GLFW_PRESS){
        if ((data->key_mods&GLFW_MOD_SHIFT)==0){
            // Drag 
            float dx = 3.*(x-data->mouse_x)/width;
            float dy = 3.*(y-data->mouse_y)/height;
            data->mouse_x = x;
            data->mouse_y = y;

            struct reb_quaternion inv = conjugate(data->view);
            float up[3] = {0.,1.,0.};
            float right[3] = {1.,0.,0.};
            float inv_right[3];
            float inv_up[3];
            multvec(inv,right,inv_right);   
            multvec(inv,up,inv_up); 

            float sin_dy = sin(dy);
            struct reb_quaternion rot_dy;
            rot_dy.x    = inv_right[0]*sin_dy;
            rot_dy.y    = inv_right[1]*sin_dy;
            rot_dy.z    = inv_right[2]*sin_dy;
            rot_dy.w    = cos(dy);
            rot_dy = normalize( rot_dy );
            data->view = mult(data->view,rot_dy);
            
            float sin_dx = sin(dx);
            struct reb_quaternion rot_dx;
            rot_dx.x    = inv_up[0]*sin_dx;
            rot_dx.y    = inv_up[1]*sin_dx;
            rot_dx.z    = inv_up[2]*sin_dx;
            rot_dx.w    = cos(dx);
            rot_dx = normalize(rot_dx);
            data->view = mult(data->view,rot_dx);
        }else{
            // Zoom
            float dx = 3.*(x-data->mouse_x)/width;
            float dy = 3.*(y-data->mouse_y)/height;
            data->mouse_x = x;
            data->mouse_y = y;
            data->scale *= (1.+dx+dy);

        }
    }
}

static void reb_display_keyboard(GLFWwindow* window, int key, int scancode, int action, int mods){
    struct reb_display_data* data = glfwGetWindowUserPointer(window);
    data->key_mods = mods;
    if (action==GLFW_PRESS){
        switch(key){
            case 'H':
                data->onscreenhelp = !data->onscreenhelp;
                break;
            case 'Q':
                data->r->status = REB_EXIT_USER;
                break;
            case ' ':
                if (data->r->status == REB_RUNNING_PAUSED){
                    printf("Resume.\n");
                    data->r->status = REB_RUNNING;
                }else{
                    printf("Pause.\n");
                    data->r->status = REB_RUNNING_PAUSED;
                }
                break;
            case 'S':
                data->spheres = (data->spheres+1)%3;
                break;
            case 'G':
                data->ghostboxes = !data->ghostboxes;
                break;
            case 'M':
                data->multisample = !data->multisample;
                if (data->multisample){
                    glEnable(GL_MULTISAMPLE); 
                }else{
                    glDisable(GL_MULTISAMPLE); 
                }
                break;
            case 'R':
                data->view.x = 0.;
                data->view.y = 0.;
                data->view.z = 0.;
                data->view.w = 1.;
                data->reference     = -1;
                data->scale = data->r->boxsize_max/2.;
                break;
            case 'D':
                data->pause = !data->pause;
                break;
            case 'W':
                data->wire = !data->wire;
                break;
            case 'T':
                data->onscreentext = !data->onscreentext;
                break;
            case 'C':
                data->clear = !data->clear;
                break;
            case 'X': 
                if (mods!=GLFW_MOD_SHIFT){
                    data->reference++;
                    if (data->reference>data->r->N) data->reference = -1;
                    printf("Reference particle: %d.\n",data->reference);
                }else{
                    data->reference--;
                    if (data->reference<-1) data->reference = data->r->N-1;
                    printf("Reference particle: %d.\n",data->reference);
                }
                break;
        }
        reb_display(window);
    }
}

static void reb_display(GLFWwindow* window){
    struct reb_display_data* data = glfwGetWindowUserPointer(window);
    if (!data){
        // No user pointer available
        return;
    }
    if (data->pause){
        return;
    }
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    double ratio = (double)width/(double)height;
    glViewport(0,0,width,height);
    if (data->clear){
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT );
    }
    glPointSize(15.*data->retina);
   
    // Precalculate matricies 
    float projection[16];
    matortho( projection,
            -1.6*ratio*data->scale, 1.6*ratio*data->scale, 
            -1.6*data->scale,1.6*data->scale,
            -2.5*data->scale,2.5*data->scale);
    float view[16];
    float tmp1[16];
    float tmp2[16];
    float tmp3[16];
    if (data->reference>=0){
        struct reb_particle p = data->particles_copy[data->reference];
        mattranslate(tmp2,-p.x,-p.y,-p.z);
        quat2mat(data->view,tmp1);
        matmult(tmp1,tmp2,view);
    }else{
        quat2mat(data->view,view);
    }
    
    for (int i=-data->ghostboxes*data->r->nghostx;i<=data->ghostboxes*data->r->nghostx;i++){
    for (int j=-data->ghostboxes*data->r->nghosty;j<=data->ghostboxes*data->r->nghosty;j++){
    for (int k=-data->ghostboxes*data->r->nghostz;k<=data->ghostboxes*data->r->nghostz;k++){
        struct reb_ghostbox gb = reb_boundary_get_ghostbox(data->r_copy, i,j,k);
        { // Particles
            mattranslate(tmp2,gb.shiftx,gb.shifty,gb.shiftz);
            matmult(view,tmp2,tmp1);
            matmult(projection,tmp1,tmp2);
            if(data->spheres>0){
                // Solid Spheres
                glEnable(GL_DEPTH_TEST);
                glUseProgram(data->sphere_shader_program);
                glBindVertexArray(data->sphere_shader_particle_vao);
                glUniformMatrix4fv(data->sphere_shader_mvp_location, 1, GL_TRUE, (GLfloat*) tmp2);
                glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, data->sphere_shader_vertex_count, data->r_copy->N);
                glBindVertexArray(0);
                glDisable(GL_DEPTH_TEST);
            }

            if(data->spheres%2==0){
                // Points
                glUseProgram(data->point_shader_program);
                glBindVertexArray(data->point_shader_particle_vao);
                glUniform4f(data->point_shader_color_location, 1.,1.,0.,0.8);
                glUniformMatrix4fv(data->point_shader_mvp_location, 1, GL_TRUE, (GLfloat*) tmp2);
                glDrawArrays(GL_POINTS, 0, data->r_copy->N);
                glBindVertexArray(0);
            }
            if (data->wire){
                // Orbits
                glUseProgram(data->orbit_shader_program);
                glBindVertexArray(data->orbit_shader_particle_vao);
                glUniformMatrix4fv(data->orbit_shader_mvp_location, 1, GL_TRUE, (GLfloat*) tmp2);
                glDrawArraysInstanced(GL_LINE_STRIP, 0, data->orbit_shader_vertex_count, data->r->N-1);
                glBindVertexArray(0);
            }
        }
        { // Box
            glUseProgram(data->box_shader_program);
            if (data->r->boundary == REB_BOUNDARY_NONE){
                glBindVertexArray(data->box_shader_cross_vao);
            }else{
                glBindVertexArray(data->box_shader_box_vao);
            }
            glUniform4f(data->box_shader_color_location, 1.,0.,0.,1.);
            matscale(tmp1,data->r->boxsize_max/2.);
            mattranslate(tmp2,gb.shiftx,gb.shifty,gb.shiftz);
            matmult(tmp2,tmp1,tmp3);
            matmult(view,tmp3,tmp1);
            matmult(projection,tmp1,tmp2);
            glUniformMatrix4fv(data->box_shader_mvp_location, 1, GL_TRUE, (GLfloat*) tmp2);
            glDrawArrays(GL_LINES, 0, 24);
            glBindVertexArray(0);
        }
    }}}
    if (data->onscreentext){ // On screen text
        glUseProgram(data->simplefont_shader_program);
        glBindVertexArray(data->simplefont_shader_vao);
        glBindTexture(GL_TEXTURE_2D,data->simplefont_tex);
        glUniform1i(glGetUniformLocation(data->simplefont_shader_program, "tex"), 0);
        glUniform2f(data->simplefont_shader_pos_location, -0.96,-0.72);
        glUniform1f(data->simplefont_shader_aspect_location, 1.9);
        glUniform1f(data->simplefont_shader_scale_location, 0.01);
        glBindBuffer(GL_ARRAY_BUFFER, data->simplefont_shader_charval_buffer);
        float val[200] = {0.};
        for (int i=0;i<sizeof(reb_logo)/sizeof(reb_logo[0]);i++){
            int j = convertLine(reb_logo[i],val);
            glUniform1f(data->simplefont_shader_ypos_location, (float)i);
            glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(val), val);
            glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, j);
        }
        
        char str[256];
        int ypos = 0;
        glUniform2f(data->simplefont_shader_pos_location, -0.70,-269./350.);
        glUniform1f(data->simplefont_shader_aspect_location, 1.4545);
        glUniform1f(data->simplefont_shader_scale_location, 16./350.);
        
        glUniform1f(data->simplefont_shader_ypos_location, ypos++);
        sprintf(str,"REBOUND v%s",reb_version_str);
        int j = convertLine(str,val);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(val), val);
        glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, j);
        
        if (data->r_copy->status == REB_RUNNING){
            sprintf(str, "Simulation is running  ");
        }else if (data->r_copy->status == REB_RUNNING_PAUSED){
            sprintf(str, "Simulation is paused   ");
        }
        glUniform1f(data->simplefont_shader_ypos_location, ypos++);
        j = convertLine(str,val);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(val), val);
        glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, j);
        
        sprintf(str, "Press h for help ");
        glUniform1f(data->simplefont_shader_ypos_location, ypos++);
        j = convertLine(str,val);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(val), val);
        glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, j);
        
        sprintf(str, "N = %d ",data->r_copy->N);
        glUniform1f(data->simplefont_shader_ypos_location, ypos++);
        j = convertLine(str,val);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(val), val);
        glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, j);
        
        glUniform1f(data->simplefont_shader_ypos_location, ypos++);
        if (data->r_copy->integrator==REB_INTEGRATOR_SEI){
            sprintf(str, "t = %f [orb]  ", data->r_copy->t*data->r_copy->ri_sei.OMEGA/2./M_PI);
        }else{
            sprintf(str, "t = %f  ", data->r_copy->t);
        }
        j = convertLine(str,val);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(val), val);
        glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, j);

        glBindVertexArray(0);
        glBindTexture(GL_TEXTURE_2D,0);
    }
    if (data->onscreenhelp){ // On screen help
        glUseProgram(data->simplefont_shader_program);
        glBindVertexArray(data->simplefont_shader_vao);
        glBindTexture(GL_TEXTURE_2D,data->simplefont_tex);
        glUniform1i(glGetUniformLocation(data->simplefont_shader_program, "tex"), 0);
        glUniform2f(data->simplefont_shader_pos_location, -0.67,0.7);
        glUniform1f(data->simplefont_shader_aspect_location, 1.4);
        glUniform1f(data->simplefont_shader_scale_location, 0.035);
        glBindBuffer(GL_ARRAY_BUFFER, data->simplefont_shader_charval_buffer);
        float val[200] = {0.};
        for (int i=0;i<sizeof(onscreenhelp)/sizeof(onscreenhelp[0]);i++){
            int j = convertLine(onscreenhelp[i],val);
            glUniform1f(data->simplefont_shader_ypos_location, (float)i);
            glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(val), val);
            glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, j);
        }
    }

    glfwSwapBuffers(window);
    return;
}



void reb_display_init(struct reb_display_data *data){
    struct reb_simulation* r = data->r;
    if (!glfwInit()){
        reb_error(r, "GLFW initialization failed.");
        return;
    }
    glfwSetErrorCallback(reb_glfw_error_callback);
    glfwWindowHint(GLFW_SAMPLES, 4);
#ifdef _APPLE
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, 1);
#endif // _APPLE

    GLFWwindow*  window = glfwCreateWindow(700, 700, "rebound", NULL, NULL);
    if (!window){
        reb_error(r,"GLFW window creation failed.");
        return;
    }

    glfwMakeContextCurrent(window);
    gladLoadGLLoader((GLADloadproc) glfwGetProcAddress);
    glfwSwapInterval(1);
    
    glfwSetWindowUserPointer(window,data); 

    // Default parameters
    { // Check if we have a retina display
        int wwidth, wheight, fwidth, fheight;
        glfwGetWindowSize(window, &wwidth, &wheight);
        glfwGetFramebufferSize(window, &fwidth, &fheight);
        data->retina = (double)fwidth/(double)wwidth;
    }
    data->spheres       = 0; 
    data->pause         = 0; 
    data->multisample = 1; 
    if (data->r->integrator==REB_INTEGRATOR_WHFAST || data->r->integrator==REB_INTEGRATOR_WHFASTHELIO){
        data->wire          = 1; 
    }else{
        data->wire          = 0; 
    }
    data->onscreentext  = 1; 
    data->clear         = 1; 
    data->ghostboxes    = 0; 
    data->reference     = -1;
    data->view.w        = 1.;
    data->p_h_copy      = NULL;
    data->p_j_copy      = NULL;
    data->eta_copy      = NULL;
    data->allocated_N_whfast = 0;
    data->allocated_N_whfasthelio = 0;
    data->scale         = r->boxsize_max/2.;

    glfwSetKeyCallback(window,reb_display_keyboard);
    glfwGetInputMode(window, GLFW_STICKY_MOUSE_BUTTONS);
    glfwSetMouseButtonCallback(window, reb_display_mouse_button);
    glfwSetCursorPosCallback(window, reb_display_cursor);
    glfwSetWindowSizeCallback(window, reb_display_resize);
    glDepthMask(GL_TRUE);
    glEnable(GL_MULTISAMPLE); 
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBlendEquation(GL_FUNC_ADD);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_FRONT); 
    {
        // Load simplefont
        unsigned char image[65536];
        for(int i=0;i<8192;i++){
            unsigned char byte = simplefont[i];
            for(int b=0;b<8;b++){
                image[i*8+b] = ((byte >> b) & 0x01)? 0xFF : 0x00;
            }
        }
        glGenTextures(1, &data->simplefont_tex);
        glBindTexture(GL_TEXTURE_2D,data->simplefont_tex);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
        float color[] = {1.,0.,0.,1.}; // This is black, given we use a red texture
        glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, color);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, 256,256, 0, GL_RED, GL_UNSIGNED_BYTE, image);
        glBindTexture(GL_TEXTURE_2D,0);
    }

    {
        const char* vertex_shader =
            "#version 330\n"
            "in vec2 vp;\n"
            "in float charpos;\n"
            "uniform float scale;\n"
            "uniform float aspect;\n"
            "uniform float ypos;\n"
            "uniform vec2 pos;\n"
            "in vec2 charval;\n"
            "in vec2 texcoord;\n"
            "out vec2 Texcoord;\n"
            "void main() {\n"
            "  gl_Position = vec4(pos,0.,0.)+vec4(scale*(vp.x+charpos/aspect+charpos/16.),scale*(vp.y-ypos),0., 1.);\n"
            "  Texcoord = vec2((charval.s+texcoord.s)/16.,(charval.t+texcoord.t)/16.00);\n"
            "}\n";
        const char* fragment_shader =
            "#version 330\n"
            "out vec4 outcolor;\n"
            "uniform sampler2D tex;\n"
            "in vec2 Texcoord;\n"
            "void main() {\n"
            "  outcolor =   vec4(0.5,0.5,0.5,(1.-texture(tex, Texcoord).r)); \n"
            "}\n";

        data->simplefont_shader_program = loadShader(vertex_shader, fragment_shader);
        data->simplefont_shader_ypos_location = glGetUniformLocation(data->simplefont_shader_program, "ypos");
        data->simplefont_shader_pos_location = glGetUniformLocation(data->simplefont_shader_program, "pos");
        data->simplefont_shader_scale_location = glGetUniformLocation(data->simplefont_shader_program, "scale");
        data->simplefont_shader_aspect_location = glGetUniformLocation(data->simplefont_shader_program, "aspect");
    }
    
    {
        const char* vertex_shader =
            "#version 330\n"
            "in vec3 vp;\n"
            "uniform mat4 mvp;\n"
            "uniform vec4 vc;\n"
            "out vec4 color;\n"
            "void main() {\n"
            "  gl_Position = mvp*vec4(vp, 1.0);\n"
            "  color = vc;\n"
            "}\n";
        const char* fragment_shader =
            "#version 330\n"
            "in vec4 color;\n"
            "out vec4 outcolor;\n"
            "void main() {\n"
            "  vec2 rel = gl_PointCoord.st;\n"
            "  rel.s -=0.5f;\n"
            "  rel.t -=0.5f;\n"
            "  if (length(rel)>0.25f){\n"
            "     outcolor = vec4(0.f,0.f,0.f,0.f); \n"
            "  }else{\n"
            "     vec4 cmod = color;\n"
            "     cmod.a*= min(1.,1.-4.*(length(rel)/0.25-0.75));\n"
            "     outcolor = cmod;\n"
            "  }\n"
            "}\n";

        data->point_shader_program = loadShader(vertex_shader, fragment_shader);
        data->point_shader_mvp_location = glGetUniformLocation(data->point_shader_program, "mvp");
        data->point_shader_color_location = glGetUniformLocation(data->point_shader_program, "vc");
    }

    {
        const char* vertex_shader =
            "#version 330\n"
            "in vec3 vp;\n"
            "uniform mat4 mvp;\n"
            "uniform vec4 vc;\n"
            "out vec4 color;\n"
            "void main() {\n"
            "  gl_Position = mvp*vec4(vp, 1.0);\n"
            "  color = vc;\n"
            "}\n";
        const char* fragment_shader =
            "#version 330\n"
            "in vec4 color;\n"
            "out vec4 outcolor;\n"
            "void main() {\n"
            "  outcolor = color;\n"
            "}\n";

        data->box_shader_program = loadShader(vertex_shader, fragment_shader);
        data->box_shader_mvp_location = glGetUniformLocation(data->box_shader_program, "mvp");
        data->box_shader_color_location = glGetUniformLocation(data->box_shader_program, "vc");
    }
    
    {
        const char* vertex_shader =
            "#version 330\n"
            "in vec3 vp;\n"
            "in float sr;\n"
            "in vec3 sp;\n"
            "out vec3 normal;\n"
            "uniform mat4 mvp;\n"
            "void main() {\n"
            "  gl_Position = mvp*(vec4(sr*vp, 1.0)+vec4(sp,0.));\n"
            "  normal = vp;\n"
            "}\n";
        const char* fragment_shader =
            "#version 330\n"
            "out vec4 outcolor;\n"
            "in vec3 normal;\n"
            "void main() {\n"
            "  vec3 lightdir = vec3(1.,1.,1.);\n"
            "  float intensity = 0.5+max(0.,0.5*dot(normalize(lightdir),normalize(normal)));\n"
            "  outcolor = vec4(intensity,intensity,intensity,1.);\n"
            "}\n";

        data->sphere_shader_program = loadShader(vertex_shader, fragment_shader);
        data->sphere_shader_mvp_location = glGetUniformLocation(data->sphere_shader_program, "mvp");
    }
    
    {
        const char* vertex_shader =
            "#version 330\n"
            "in vec3 focus;\n"
            "in vec3 aef;\n"
            "in vec3 omegaOmegainc;\n"
            "in float lintwopi;\n"
            "out float lin;\n"
            "uniform mat4 mvp;\n"
            "const float M_PI = 3.14159265359;\n"
            "void main() {\n"
            "   float a = aef.x;\n"
            "   float e = aef.y;\n"
            "   float f = aef.z+lintwopi;\n"
            "   lin = lintwopi/(M_PI*2.);\n"
            "   if (e>1.){\n"
            "       float theta_max = acos(-1./e);\n"
            "       f = 0.0001-theta_max+1.9998*lin*theta_max;\n"
            "       lin = sqrt(min(0.5,lin));\n"
            "   }\n"
            "   float omega = omegaOmegainc.x;\n"
            "   float Omega = omegaOmegainc.y;\n"
            "   float inc = omegaOmegainc.z;\n"
            "   float r = a*(1.-e*e)/(1. + e*cos(f));\n"
            "   float cO = cos(Omega);\n"
            "   float sO = sin(Omega);\n"
            "   float co = cos(omega);\n"
            "   float so = sin(omega);\n"
            "   float cf = cos(f);\n"
            "   float sf = sin(f);\n"
            "   float ci = cos(inc);\n"
            "   float si = sin(inc);\n"
            "   vec3 pos = vec3(r*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci),r*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci),+ r*(so*cf+co*sf)*si);\n"
            "    gl_Position = mvp*(vec4(focus+pos, 1.0));\n"
            "}\n";
        const char* fragment_shader =
            "#version 330\n"
            "out vec4 outcolor;\n"
            "in float lin;\n"
            "void main() {\n"
            "  outcolor = vec4(1.,1.,1.,sqrt(lin));\n"
            "}\n";

        data->orbit_shader_program = loadShader(vertex_shader, fragment_shader);
        data->orbit_shader_mvp_location = glGetUniformLocation(data->orbit_shader_program, "mvp");
    }
    
    // Create simplefont mesh
    glUseProgram(data->simplefont_shader_program);
    glGenVertexArrays(1, &data->simplefont_shader_vao);
    glBindVertexArray(data->simplefont_shader_vao);
    GLuint sfvp = glGetAttribLocation(data->simplefont_shader_program,"vp");
    glEnableVertexAttribArray(sfvp);
    GLuint sftexcoordp = glGetAttribLocation(data->simplefont_shader_program,"texcoord");
    glEnableVertexAttribArray(sftexcoordp);
    GLuint simplefont_shader_charval_location = glGetAttribLocation(data->simplefont_shader_program,"charval");
    glEnableVertexAttribArray(simplefont_shader_charval_location);
    GLuint simplefont_shader_charpos_location = glGetAttribLocation(data->simplefont_shader_program,"charpos");
    glEnableVertexAttribArray(simplefont_shader_charpos_location);
    
    float simplefont_data[] = {
        0., 0., 0., 1., 
        0., 1., 0., 0.0,
        1., 0., 1., 1.,
        1., 1., 1., 0.0
    };
    GLuint simplefont_buffer;
    glGenBuffers(1, &simplefont_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, simplefont_buffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(simplefont_data), simplefont_data, GL_STATIC_DRAW);
    
    glVertexAttribPointer(sfvp, 2, GL_FLOAT, GL_FALSE, sizeof(float)*4, NULL);
    glVertexAttribPointer(sftexcoordp, 2, GL_FLOAT, GL_FALSE, sizeof(float)*4, (void *)(sizeof(float)*2));

    GLuint simplefont_shader_charpos_buffer;
    glGenBuffers(1, &simplefont_shader_charpos_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, simplefont_shader_charpos_buffer);
    float charpos[100];
    for(int i=0;i<100;i++){
        charpos[i] = (float)i;
    }
    glBufferData(GL_ARRAY_BUFFER, sizeof(charpos), charpos, GL_STATIC_DRAW);
    glVertexAttribPointer(simplefont_shader_charpos_location, 1, GL_FLOAT, GL_FALSE, sizeof(float)*1, NULL);

    glGenBuffers(1, &data->simplefont_shader_charval_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, data->simplefont_shader_charval_buffer);
    glBufferData(GL_ARRAY_BUFFER, 200*sizeof(float), NULL,GL_DYNAMIC_DRAW);
    glVertexAttribPointer(simplefont_shader_charval_location, 2, GL_FLOAT, GL_FALSE, sizeof(float)*2, NULL);

    glVertexAttribDivisor(sfvp, 0); 
    glVertexAttribDivisor(sftexcoordp, 0); 
    glVertexAttribDivisor(simplefont_shader_charpos_location,1);
    glVertexAttribDivisor(simplefont_shader_charval_location,1);

    glBindVertexArray(0);


    // Particle data is dynamic
    glUseProgram(data->point_shader_program);
    glGenVertexArrays(1, &data->point_shader_particle_vao);
    glBindVertexArray(data->point_shader_particle_vao);
    GLuint pvp = glGetAttribLocation(data->point_shader_program,"vp");
    glEnableVertexAttribArray(pvp);

    int particle_data_allocated_N = 1;
    struct reb_particle_opengl* particle_data = malloc(sizeof(struct reb_particle_opengl)*particle_data_allocated_N);
    GLuint particle_buffer;
    glGenBuffers(1, &particle_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, particle_buffer);
    glBufferData(GL_ARRAY_BUFFER, particle_data_allocated_N*sizeof(struct reb_particle_opengl), NULL, GL_STATIC_DRAW);
    
    glVertexAttribPointer(pvp, 3, GL_DOUBLE, GL_FALSE, sizeof(double)*7, NULL);
    glBindVertexArray(0);
    
    // Create cross mesh
    glUseProgram(data->box_shader_program);
    glGenVertexArrays(1, &data->box_shader_cross_vao);
    glBindVertexArray(data->box_shader_cross_vao);
    GLuint cvp = glGetAttribLocation(data->box_shader_program,"vp");
    glEnableVertexAttribArray(cvp);
    
    float cross_data[] = {
        -0.04,0.,0., +0.04,0.,0., 0.,-0.04,0., 0.,+0.04,0., 0.,0.,-0.04, 0.,0.,+0.04,
    };
    GLuint cross_buffer;
    glGenBuffers(1, &cross_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, cross_buffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cross_data), cross_data, GL_STATIC_DRAW);
    
    glVertexAttribPointer(cvp, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    glBindVertexArray(0);
    
    // Create box mesh
    glGenVertexArrays(1, &data->box_shader_box_vao);
    glBindVertexArray(data->box_shader_box_vao);
    GLuint bvp = glGetAttribLocation(data->box_shader_program,"vp");
    glEnableVertexAttribArray(bvp);

    float box_data[] = {
        -1,-1,-1, 1,-1,-1, 1,-1,-1, 1, 1,-1, 1, 1,-1, -1, 1,-1, -1, 1,-1, -1,-1,-1,
        -1,-1, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1, 1, -1, 1, 1, -1, 1, 1, -1,-1, 1,
        -1,-1,-1, -1,-1, 1, -1, 1,-1, -1, 1, 1, 1, 1,-1, 1, 1, 1, 1,-1,-1, 1,-1, 1,
    };
    GLuint box_buffer;
    glGenBuffers(1, &box_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, box_buffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(box_data), box_data, GL_STATIC_DRAW);
    
    glVertexAttribPointer(bvp, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    glBindVertexArray(0);
    
                
    // Sphere data
    glUseProgram(data->sphere_shader_program);
    glGenVertexArrays(1, &data->sphere_shader_particle_vao);
    glBindVertexArray(data->sphere_shader_particle_vao);
    GLuint svp = glGetAttribLocation(data->sphere_shader_program,"vp");
    glEnableVertexAttribArray(svp);
    GLuint ssp = glGetAttribLocation(data->sphere_shader_program,"sp");
    glEnableVertexAttribArray(ssp);
    GLuint ssr = glGetAttribLocation(data->sphere_shader_program,"sr");
    glEnableVertexAttribArray(ssr);
   
    float* sphere_data = malloc(sizeof(float)*5000);
    int count = 0;
    int ni = 20;
    int nj = 20;
    for(int i=0;i<ni;i++){
        for(int j=0;j<nj;j++){
            sphere_data[count*3+0] = sinf((float)i/(float)ni*M_PI)*sinf((float)j/(float)nj*2.*M_PI);
            sphere_data[count*3+1] = sinf((float)i/(float)ni*M_PI)*cosf((float)j/(float)nj*2.*M_PI);
            sphere_data[count*3+2] = cosf((float)i/(float)ni*M_PI);
            count++;
            sphere_data[count*3+0] = sinf((float)(i+1)/(float)ni*M_PI)*sinf((float)j/(float)nj*2.*M_PI);
            sphere_data[count*3+1] = sinf((float)(i+1)/(float)ni*M_PI)*cosf((float)j/(float)nj*2.*M_PI);
            sphere_data[count*3+2] = cosf((float)(i+1)/(float)ni*M_PI);
            count++;
        }
    }
    data->sphere_shader_vertex_count = count;
    GLuint sphere_vertex_buffer;
    glGenBuffers(1, &sphere_vertex_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, sphere_vertex_buffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*3*count, sphere_data, GL_STATIC_DRAW);
    free(sphere_data);
    glVertexAttribPointer(svp, 3, GL_FLOAT, GL_FALSE, 0, NULL);

    glBindBuffer(GL_ARRAY_BUFFER, particle_buffer);
    glVertexAttribPointer(ssp, 3, GL_DOUBLE, GL_FALSE, sizeof(double)*7, NULL);
    glVertexAttribPointer(ssr, 1, GL_DOUBLE, GL_FALSE, sizeof(double)*7, (void*)(sizeof(double)*6));

    glVertexAttribDivisor(svp, 0); 
    glVertexAttribDivisor(data->sphere_shader_mvp_location, 0); 
    glVertexAttribDivisor(ssp, 1);
    glVertexAttribDivisor(ssr, 1);
    
    glBindVertexArray(0);

    
    // Orbit data
    glUseProgram(data->orbit_shader_program);
    glGenVertexArrays(1, &data->orbit_shader_particle_vao);
    glBindVertexArray(data->orbit_shader_particle_vao);
    GLuint olintwopip = glGetAttribLocation(data->orbit_shader_program,"lintwopi");
    glEnableVertexAttribArray(olintwopip);
    GLuint ofocusp = glGetAttribLocation(data->orbit_shader_program,"focus");
    glEnableVertexAttribArray(ofocusp);
    GLuint oaefp = glGetAttribLocation(data->orbit_shader_program,"aef");
    glEnableVertexAttribArray(oaefp);
    GLuint oomegaOmegaincp = glGetAttribLocation(data->orbit_shader_program,"omegaOmegainc");
    glEnableVertexAttribArray(oomegaOmegaincp);
   
    data->orbit_shader_vertex_count = 500;
    float* lin_data = malloc(sizeof(float)*data->orbit_shader_vertex_count);
    for(int i=0;i<data->orbit_shader_vertex_count;i++){
        lin_data[i] = (float)i/(float)(data->orbit_shader_vertex_count-1)*2.*M_PI;
    }
    GLuint orbit_vertex_buffer;
    glGenBuffers(1, &orbit_vertex_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, orbit_vertex_buffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*data->orbit_shader_vertex_count, lin_data, GL_STATIC_DRAW);
    free(lin_data);
    glVertexAttribPointer(olintwopip, 1, GL_FLOAT, GL_FALSE, 0, NULL);


    struct reb_orbit_opengl* orbit_data = malloc(sizeof(struct reb_orbit_opengl)*particle_data_allocated_N);
    GLuint orbit_buffer;
    glGenBuffers(1, &orbit_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, orbit_buffer);
    glBufferData(GL_ARRAY_BUFFER, particle_data_allocated_N*sizeof(struct reb_orbit_opengl), NULL, GL_STATIC_DRAW);
    glVertexAttribPointer(ofocusp, 3, GL_DOUBLE, GL_FALSE, sizeof(double)*9, NULL);
    glVertexAttribPointer(oaefp, 3, GL_DOUBLE, GL_FALSE, sizeof(double)*9, (void*)(sizeof(double)*3));
    glVertexAttribPointer(oomegaOmegaincp, 3, GL_DOUBLE, GL_FALSE, sizeof(double)*9, (void*)(sizeof(double)*6));

    glVertexAttribDivisor(olintwopip, 0); 
    glVertexAttribDivisor(data->orbit_shader_mvp_location, 0); 
    glVertexAttribDivisor(ofocusp, 1);
    glVertexAttribDivisor(oaefp, 1);
    glVertexAttribDivisor(oomegaOmegaincp, 1);
    
    glBindVertexArray(0);


    // Main display loop
    while(!glfwWindowShouldClose(window) && data->return_status<0){
        { // lock mutex for update
            pthread_mutex_lock(data->mutex);    
            if (r->N>data->allocated_N){
                data->allocated_N = r->N;
                data->r_copy = realloc(data->r_copy,sizeof(struct reb_simulation));
                data->particles_copy = realloc(data->particles_copy,r->N*sizeof(struct reb_particle));
            }
            memcpy(data->r_copy, r, sizeof(struct reb_simulation));
            memcpy(data->particles_copy, r->particles, sizeof(struct reb_particle)*r->N);
            data->r_copy->particles = data->particles_copy;
            if (r->integrator==REB_INTEGRATOR_WHFAST && r->ri_whfast.is_synchronized==0){
                if (r->ri_whfast.allocated_N > data->allocated_N_whfast){
                    data->allocated_N_whfast = r->ri_whfast.allocated_N;
                    data->p_j_copy = realloc(data->p_j_copy,data->allocated_N_whfast*sizeof(struct reb_particle));
                    data->eta_copy = realloc(data->eta_copy,data->allocated_N_whfast*sizeof(double));
                }
                memcpy(data->p_j_copy, r->ri_whfast.p_j, data->allocated_N_whfast*sizeof(struct reb_particle));
                memcpy(data->eta_copy, r->ri_whfast.eta, data->allocated_N_whfast*sizeof(double));
            }
            data->r_copy->ri_whfast.p_j= data->p_j_copy;
            data->r_copy->ri_whfast.eta= data->eta_copy;
            if (r->integrator==REB_INTEGRATOR_WHFASTHELIO && r->ri_whfasthelio.is_synchronized==0){
                if (r->ri_whfasthelio.allocated_N > data->allocated_N_whfast){
                    data->allocated_N_whfasthelio = r->ri_whfasthelio.allocated_N;
                    data->p_h_copy = realloc(data->p_h_copy,data->allocated_N_whfasthelio*sizeof(struct reb_particle));
                }
                memcpy(data->p_h_copy, r->ri_whfasthelio.p_h, data->allocated_N_whfasthelio*sizeof(struct reb_particle));
            }
            data->r_copy->ri_whfasthelio.p_h= data->p_h_copy;
            pthread_mutex_unlock(data->mutex);  
        }

        // this only does something for WHFAST + WHFASTHELIO
        reb_integrator_synchronize(data->r_copy);
           
        // Update data on GPU 
        if (data->r_copy->N > particle_data_allocated_N){
            particle_data_allocated_N = data->r_copy->N;
            particle_data = realloc(particle_data, particle_data_allocated_N*sizeof(struct reb_particle_opengl));
            glBindBuffer(GL_ARRAY_BUFFER, particle_buffer);
            glBufferData(GL_ARRAY_BUFFER, particle_data_allocated_N*sizeof(struct reb_particle_opengl), NULL, GL_STATIC_DRAW);
            orbit_data = realloc(orbit_data, particle_data_allocated_N*sizeof(struct reb_orbit_opengl));
            glBindBuffer(GL_ARRAY_BUFFER, orbit_buffer);
            glBufferData(GL_ARRAY_BUFFER, particle_data_allocated_N*sizeof(struct reb_orbit_opengl), NULL, GL_STATIC_DRAW);
        }
        for (int i=0;i<data->r_copy->N;i++){
            struct reb_particle p = data->r_copy->particles[i];
            particle_data[i].x  = p.x;
            particle_data[i].y  = p.y;
            particle_data[i].z  = p.z;
            particle_data[i].vx = p.vx;
            particle_data[i].vy = p.vy;
            particle_data[i].vz = p.vz;
            particle_data[i].r  = p.r;
        }
        glBindBuffer(GL_ARRAY_BUFFER, particle_buffer);
        glBufferSubData(GL_ARRAY_BUFFER, 0, data->r_copy->N*sizeof(struct reb_particle_opengl), particle_data);

        if (data->wire){
            struct reb_particle com = data->r_copy->particles[0];
            for (int i=1;i<r->N;i++){
                struct reb_particle p = data->r_copy->particles[i];
                orbit_data[i-1].x  = com.x;
                orbit_data[i-1].y  = com.y;
                orbit_data[i-1].z  = com.z;
                struct reb_orbit o = reb_tools_particle_to_orbit(data->r_copy->G, p,com);
                orbit_data[i-1].a = o.a;
                orbit_data[i-1].e = o.e;
                orbit_data[i-1].f = o.f;
                orbit_data[i-1].omega = o.omega;
                orbit_data[i-1].Omega = o.Omega;
                orbit_data[i-1].inc = o.inc;
                com = reb_get_com_of_pair(p,com);
            }
            glBindBuffer(GL_ARRAY_BUFFER, orbit_buffer);
            glBufferSubData(GL_ARRAY_BUFFER, 0, (data->r_copy->N-1)*sizeof(struct reb_orbit_opengl), orbit_data);
        }

        // Do actual drawing
        reb_display(window);
        glfwPollEvents();
    }
    glfwDestroyWindow(window);
    glfwTerminate();

    data->r->status = REB_EXIT_USER;
    data->return_status = REB_EXIT_USER;
    free(data->r_copy);
    free(data->particles_copy);
    free(data->eta_copy);
    free(data->p_j_copy);
    free(data->p_h_copy);
    free(particle_data);
    free(orbit_data);

}

#endif // OPENGL
