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
#define DEG2RAD (M_PI/180.)
#include <stdio.h>
#include <stdlib.h>
#ifndef _WIN32
#include <pthread.h>
#endif // _WIN32
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include "rebound.h"
#include "display.h"
#include "tools.h"
#include "particle.h"
#include "boundary.h"
#include "display.h"
#include "output.h"
#include "integrator.h"
#define MAX(a, b) ((a) < (b) ? (b) : (a))       ///< Returns the maximum of a and b

#ifdef OPENGL
#include "simplefont.h"
#define GLFW_INCLUDE_NONE
#include "glad.h"
#include <GLFW/glfw3.h>

#ifdef __EMSCRIPTEN__
#include <emscripten/fetch.h>
// Need to use emscripten_ version of these functions because types are wrong otherwise
void emscripten_glVertexAttribDivisor(GLuint index, GLuint divisor);
void emscripten_glDrawArraysInstanced(GLenum mode, GLint first, GLsizei count, GLsizei instancecount);
#define reb_glVertexAttribDivisor emscripten_glVertexAttribDivisor
#define reb_glDrawArraysInstanced emscripten_glDrawArraysInstanced

// Getting size from canvas element
EM_JS(int, canvas_get_width, (), {
  return document.getElementById("canvas").scrollWidth;
});
EM_JS(int, canvas_get_height, (), {
  return document.getElementById("canvas").scrollHeight;
});

EM_JS(void, reb_overlay_update, (const char* text, int status), {
    var overlaytext = document.getElementById("overlaytext");
    if (overlaytext){
        overlaytext.innerHTML = UTF8ToString(text);
    }
    var overlay = document.getElementById("overlay");
    if (overlay){
        if (status==-3){ // Pause
            overlay.style.backgroundColor = "rgba(100.0, 100.0, 0.0, 0.5)";
        }else if (status==0 || status==5){ // Finished.
            overlay.style.backgroundColor = "rgba(0.0, 255.0, 0.0, 0.5)";
        }else if (status==10){ // Connection error.
            overlay.style.backgroundColor = "rgba(255.0, 0.0, 0.0, 0.5)";
        }else{
            overlay.style.backgroundColor = "rgba(0, 0, 0, 0.5)";
        }
    }
});
EM_JS(void, reb_overlay_help_set_text, (const char* text), {
    var overlaytext = document.getElementById("overlaytext-help");
    if (overlaytext){
        overlaytext.innerHTML = UTF8ToString(text);
    }
});
EM_JS(int, reb_overlay_help_show, (int show), {
    var overlaytoggle = document.getElementById("overlay-toggle");
    if (overlaytoggle){
        if (overlaytoggle.innerHTML == "1"){
            overlaytoggle.innerHTML = "";
            show = !show;
        }
    }
    var overlayhelp = document.getElementById("overlay-help");
    if (show){
        overlayhelp.style.display = "block";
    }else{
        overlayhelp.style.display = "none";
    }
    return show;
});
EM_JS(void, reb_overlay_hide, (int hide), {
    var overlay = document.getElementById("overlay");
    if (hide){
        overlay.style.display = "none";
    }else{
        overlay.style.display = "block";
    }

});
#else
#define reb_glVertexAttribDivisor glVertexAttribDivisor
#define reb_glDrawArraysInstanced glDrawArraysInstanced
#endif

void reb_render_frame(void* p);
static void reb_display_set_default_scale(struct reb_simulation* const r);
                
static const char* onscreenhelp[] = { 
                "REBOUND mouse and keyboard commands",
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
                " x/X     | Move to a coordinate system centerd ",
                "         | on a particle (note: does not work if", 
                "         | particle array is resorted)",
                " c       | Toggle clear screen after each time-step",
                " m       | Toggle multisampling",
                " w       | Draw orbits as wires",
                " t       | Show/hide logo, time, timestep and number",
                "         | of particles.",
                "----------------------------------------------------"
};


static void matscale(float mat[16], float s){
    mat[0] = s; mat[1] = 0.; mat[2] = 0.; mat[3] = 0.;
    mat[4] = 0.; mat[5] = s; mat[6] = 0.; mat[7] = 0.;
    mat[8] = 0.; mat[9] = 0.; mat[10] = s; mat[11] = 0.;
    mat[12] = 0.; mat[13] = 0.; mat[14] = 0.; mat[15] = 1.;
}

static void matscale3(float mat[16], float s[3]){
    mat[0] = s[0]; mat[1] = 0.; mat[2] = 0.; mat[3] = 0.;
    mat[4] = 0.; mat[5] = s[1]; mat[6] = 0.; mat[7] = 0.;
    mat[8] = 0.; mat[9] = 0.; mat[10] = s[2]; mat[11] = 0.;
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
static void rotation2mat(struct reb_rotation A, float mat[16]){
    float xx = A.ix * A.ix; float xy = A.ix * A.iy; float xz = A.ix * A.iz;
    float xw = A.ix * A.r; float yy = A.iy * A.iy; float yz = A.iy * A.iz;
    float yw = A.iy * A.r; float zz = A.iz * A.iz; float zw = A.iz * A.r;
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
    if (logLength > 1) { // 0 should work - bug in emscripten?
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
    if (!data){
        printf("Error accessing data in reb_display_mouse_button\n");
        return;
    }
    data->mouse_action = action;
}

static void reb_display_resize(GLFWwindow* window, int x, int y){
    struct reb_display_data* data = glfwGetWindowUserPointer(window);
    if (!data){
        printf("Error accessing data in reb_display_resize\n");
        return;
    }
}

static void reb_display_cursor(GLFWwindow* window, double x, double y){
    struct reb_display_data* data = glfwGetWindowUserPointer(window);
    if (!data){
        printf("Error accessing data in reb_display_cursor\n");
        return;
    }
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
        data->mouse_scale = data->scale;
        return;
    }
    if (data->mouse_action==GLFW_PRESS){
        if ((data->key_mods&GLFW_MOD_SHIFT)==0){
            // Drag 
            float dx = 3.*(x-data->mouse_x)/width;
            float dy = 3.*(y-data->mouse_y)/height;
            data->mouse_x = x;
            data->mouse_y = y;

            struct reb_rotation inv = reb_rotation_conjugate(data->view);
            struct reb_vec3d up = {.x=0.,.y=1.,.z=0.};
            struct reb_vec3d right = {.x=1.,.y=0.,.z=0.};
            struct reb_vec3d inv_right = reb_vec3d_rotate(right, inv);
            struct reb_vec3d inv_up = reb_vec3d_rotate(up, inv);

            float sin_dy = sin(dy);
            struct reb_rotation rot_dy;
            rot_dy.ix    = inv_right.x*sin_dy;
            rot_dy.iy    = inv_right.y*sin_dy;
            rot_dy.iz    = inv_right.z*sin_dy;
            rot_dy.r    = cos(dy);
            rot_dy = reb_rotation_normalize( rot_dy );
            data->view = reb_rotation_mul(data->view,rot_dy);
            
            float sin_dx = sin(dx);
            struct reb_rotation rot_dx;
            rot_dx.ix    = inv_up.x*sin_dx;
            rot_dx.iy    = inv_up.y*sin_dx;
            rot_dx.iz    = inv_up.z*sin_dx;
            rot_dx.r    = cos(dx);
            rot_dx = reb_rotation_normalize(rot_dx);
            data->view = reb_rotation_mul(data->view,rot_dx);
        }else{
            // Zoom
            float ix = data->mouse_x/width-0.5;
            float iy = data->mouse_y/height-0.5;
            float ir = sqrt(ix*ix + iy*iy);
            float nx = x/width-0.5;
            float ny = y/height-0.5;
            float nr = sqrt(nx*nx + ny*ny);
            data->scale = data->mouse_scale*ir/nr;

        }
    }
}

void reb_display_keyboard(GLFWwindow* window, int key, int scancode, int action, int mods){
    struct reb_display_data* data = glfwGetWindowUserPointer(window);
    if (!data){
        printf("Error accessing data in reb_display_keyboard\n");
        return;
    }
    if (!data->r){
        printf("Error accessing data->r in reb_display_keyboard\n");
        return;
    }
    // User defined keys:
    int skip_default_keys = 0;
    if (data->r->key_callback){
        skip_default_keys = data->r->key_callback(data->r, key);
    } 
    if (skip_default_keys){
        return;
    }
    // Default keys:
    data->key_mods = mods;
    if (action==GLFW_PRESS){
        switch(key){
            case 'H':
                data->onscreenhelp = !data->onscreenhelp;
                break;
            case 'Q':
                data->r->status = REB_STATUS_USER;
                break;
            case ' ':
                if (data->r->status == REB_STATUS_PAUSED){
                    printf("Resume.\n");
                    data->r->status = REB_STATUS_RUNNING;
                }else if (data->r->status == REB_STATUS_RUNNING || data->r->status == REB_STATUS_LAST_STEP){
                    printf("Pause.\n");
                    data->r->status = REB_STATUS_PAUSED;
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
                if (data->view.r ==1.){
                    data->view.ix = 1./sqrt(2.);
                    data->view.iy = 0.;
                    data->view.iz = 0.;
                    data->view.r = 1./sqrt(2.);
                }else if (data->view.ix == 1./sqrt(2.)){
                    data->view.ix = 0.;
                    data->view.iy = -1./sqrt(2.);
                    data->view.iz = 0.;
                    data->view.r = 1./sqrt(2.);
                }else{
                    data->view.ix = 0.;
                    data->view.iy = 0.;
                    data->view.iz = 0.;
                    data->view.r = 1.;
                }
                data->reference     = -1;
                reb_display_set_default_scale(data->r);
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
    }
}

// Actual rendering
// Makes a copy of the simulation first.
void reb_render_frame(void* p){
    struct reb_display_data* data = (struct reb_display_data*)p;
    struct reb_simulation* r = data->r;
    if (!data){
        printf("reb_display_data undefinded in reb_render_frame().\n");
        return;
    }
    int width, height;
#ifdef __EMSCRIPTEN__
    // Need to query canvas size using JS, set window size, then read framebuffer size.
    width = canvas_get_width();
    height = canvas_get_height();
    int cwidth, cheight;
    glfwGetWindowSize(data->window, &cwidth, &cheight);
    if (cwidth!=width || cheight !=height){
        glfwSetWindowSize(data->window, width, height);
    }
#endif
    glfwGetFramebufferSize(data->window, &width, &height);

    struct reb_simulation* r_copy = r->display_data->r_copy;
    if (!r_copy){
        data->r_copy = reb_simulation_create();
        r_copy = data->r_copy;
    }
    
    // lock mutex for update
    data->need_copy = 1;
    int wait_count = 0;
    const int wait_count_max = 10;
    int ret_try = EBUSY;
    while (wait_count<wait_count_max && ret_try){// wait at most one frame.
        ret_try = pthread_mutex_trylock(&data->mutex);
        if (ret_try){ // not locked
            usleep(1./120.*1e6/wait_count_max); 
            wait_count++;
        }
    }

    if (!ret_try){
        // Copy if lock obtained. Otherwise use old data.
        enum reb_simulation_binary_error_codes warnings = REB_SIMULATION_BINARY_WARNING_NONE;
        reb_simulation_copy_with_messages(data->r_copy,r,&warnings);
        data->need_copy = 0;
        pthread_mutex_unlock(&(data->mutex));  
    }

    // prepare data (incl orbit calculation)
    const int N_real = r_copy->N - r_copy->N_var;
    if (N_real > data->N_allocated){
        data->N_allocated = N_real;
        data->particle_data = realloc(data->particle_data, data->N_allocated*sizeof(struct reb_particle_opengl));
        data->orbit_data = realloc(data->orbit_data, data->N_allocated*sizeof(struct reb_orbit_opengl));
        // Resize memory if needed
        glBindBuffer(GL_ARRAY_BUFFER, data->particle_buffer);
        glBufferData(GL_ARRAY_BUFFER, data->N_allocated*sizeof(struct reb_particle_opengl), NULL, GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, data->orbit_buffer);
        glBufferData(GL_ARRAY_BUFFER, data->N_allocated*sizeof(struct reb_orbit_opengl), NULL, GL_STATIC_DRAW);
    }

    // this only does something for WHFAST
    reb_simulation_synchronize(r_copy);
       
    // Update data on GPU 
    for (unsigned int i=0;i<N_real;i++){
        struct reb_particle p = r_copy->particles[i];
        data->particle_data[i].x  = (float)p.x;
        data->particle_data[i].y  = (float)p.y;
        data->particle_data[i].z  = (float)p.z;
        data->particle_data[i].vx = (float)p.vx;
        data->particle_data[i].vy = (float)p.vy;
        data->particle_data[i].vz = (float)p.vz;
        data->particle_data[i].r  = (float)p.r;
    }
    if (data->wire && N_real>1){
        struct reb_particle com = r_copy->particles[0];
        for (unsigned int i=1;i<N_real;i++){
            struct reb_particle p = r_copy->particles[i];
            data->orbit_data[i-1].x  = (float)com.x;
            data->orbit_data[i-1].y  = (float)com.y;
            data->orbit_data[i-1].z  = (float)com.z;
            struct reb_orbit o = reb_orbit_from_particle(r_copy->G, p,com);
            data->orbit_data[i-1].a = (float)o.a;
            data->orbit_data[i-1].e = (float)o.e;
            data->orbit_data[i-1].f = (float)o.f;
            data->orbit_data[i-1].omega = (float)o.omega;
            data->orbit_data[i-1].Omega = (float)o.Omega;
            data->orbit_data[i-1].inc = (float)o.inc;
            com = reb_particle_com_of_pair(p,com);
        }
    }
    if (N_real>0){
        // Fill memory (but not resize)
        glBindBuffer(GL_ARRAY_BUFFER, data->particle_buffer);
        glBufferSubData(GL_ARRAY_BUFFER, 0, N_real*sizeof(struct reb_particle_opengl), data->particle_data);
        glBindBuffer(GL_ARRAY_BUFFER, data->orbit_buffer);
        glBufferSubData(GL_ARRAY_BUFFER, 0, (N_real-1)*sizeof(struct reb_orbit_opengl), data->orbit_data);
    }

    // Do actual drawing
    double ratio = (double)width/(double)height;
    glViewport(0,0,width,height);
    if (data->clear){
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT );
    }

#ifndef __EMSCRIPTEN__
    glPointSize(15.*data->retina);
#endif
   
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
        struct reb_particle p = data->r_copy->particles[data->reference];
        mattranslate(tmp2,-p.x,-p.y,-p.z);
        rotation2mat(data->view,tmp1);
        matmult(tmp1,tmp2,view);
    }else{
        rotation2mat(data->view,view);
    }
    
    for (int i=-data->ghostboxes*data->r_copy->N_ghost_x;i<=data->ghostboxes*data->r_copy->N_ghost_x;i++){
    for (int j=-data->ghostboxes*data->r_copy->N_ghost_y;j<=data->ghostboxes*data->r_copy->N_ghost_y;j++){
    for (int k=-data->ghostboxes*data->r_copy->N_ghost_z;k<=data->ghostboxes*data->r_copy->N_ghost_z;k++){
        struct reb_vec6d gb = reb_boundary_get_ghostbox(data->r_copy, i,j,k);
        { // Particles
            mattranslate(tmp2,gb.x,gb.y,gb.z);
            matmult(view,tmp2,tmp1);
            matmult(projection,tmp1,tmp2);
            if(data->spheres>0){
                // Solid Spheres
                glEnable(GL_DEPTH_TEST);
                glUseProgram(data->sphere_shader_program);
                glBindVertexArray(data->sphere_shader_particle_vao);
                glUniformMatrix4fv(data->sphere_shader_mvp_location, 1, GL_TRUE, (GLfloat*) tmp2);
                reb_glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, data->sphere_shader_vertex_count, N_real);
                glBindVertexArray(0);
                glDisable(GL_DEPTH_TEST);
            }

            if(data->spheres%2==0){
                // Points
                glUseProgram(data->point_shader_program);
                glBindVertexArray(data->point_shader_particle_vao);
                glUniform4f(data->point_shader_color_location, 1.,1.,0.,0.8);
                glUniformMatrix4fv(data->point_shader_mvp_location, 1, GL_TRUE, (GLfloat*) tmp2);
                glDrawArrays(GL_POINTS, 0, N_real);
                glBindVertexArray(0);
            }
            if (data->wire){
                // Orbits
                glUseProgram(data->orbit_shader_program);
                glBindVertexArray(data->orbit_shader_particle_vao);
                glUniformMatrix4fv(data->orbit_shader_mvp_location, 1, GL_TRUE, (GLfloat*) tmp2);
                reb_glDrawArraysInstanced(GL_LINE_STRIP, 0, data->orbit_shader_vertex_count, N_real-1);
                glBindVertexArray(0);
            }
        }
        { // Box
            glUseProgram(data->box_shader_program);
            if (data->r_copy->boundary == REB_BOUNDARY_NONE){
                glBindVertexArray(data->box_shader_cross_vao);
                matscale(tmp1,data->scale);
            }else{
                glBindVertexArray(data->box_shader_box_vao);
                float boxsize[3] = {data->r_copy->boxsize.x/2., data->r_copy->boxsize.y/2., data->r_copy->boxsize.z/2.};
                matscale3(tmp1,boxsize);
            }
            glUniform4f(data->box_shader_color_location, 1.,0.,0.,1.);
            mattranslate(tmp2,gb.x,gb.y,gb.z);
            matmult(tmp2,tmp1,tmp3);
            matmult(view,tmp3,tmp1);
            matmult(projection,tmp1,tmp2);
            glUniformMatrix4fv(data->box_shader_mvp_location, 1, GL_TRUE, (GLfloat*) tmp2);
            if (data->r_copy->boundary == REB_BOUNDARY_NONE){
                glDrawArrays(GL_LINES, 0, 6);
            }else{
                glDrawArrays(GL_LINES, 0, 24);
            }
            glBindVertexArray(0);
        }
    }}}
#ifndef __EMSCRIPTEN__
    if (data->onscreentext){ // On screen text
        glUseProgram(data->simplefont_shader_program);
        glBindVertexArray(data->simplefont_shader_vao);
        glBindTexture(GL_TEXTURE_2D,data->simplefont_tex);
        glUniform1i(glGetUniformLocation(data->simplefont_shader_program, "tex"), 0);
        glUniform2f(data->simplefont_shader_pos_location, -0.96,-0.72);
        glUniform1f(data->simplefont_shader_aspect_location, 1.9);
        glUniform1f(data->simplefont_shader_screen_aspect_location, 1./ratio);
        glUniform1f(data->simplefont_shader_scale_location, 0.01);
        glBindBuffer(GL_ARRAY_BUFFER, data->simplefont_shader_charval_buffer);
        float val[200] = {0.};
        for (int i=0;i<sizeof(reb_logo)/sizeof(reb_logo[0]);i++){
            int j = convertLine(reb_logo[i],val);
            glUniform1f(data->simplefont_shader_ypos_location, (float)i);
            glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(val), val);
            reb_glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, j);
        }
        
        char str[256];
        int ypos = 0;
        glUniform2f(data->simplefont_shader_pos_location, -0.70,-269./350.);
        glUniform1f(data->simplefont_shader_aspect_location, 1.4);
        glUniform1f(data->simplefont_shader_scale_location, 16./350.);
        
        glUniform1f(data->simplefont_shader_ypos_location, ypos++);
        sprintf(str,"REBOUND v%s",reb_version_str);
        int j = convertLine(str,val);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(val), val);
        reb_glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, j);
        
        if (data->r_copy->status == REB_STATUS_RUNNING){
            sprintf(str, "Simulation is running  ");
        }else if (data->r_copy->status == REB_STATUS_PAUSED){
            sprintf(str, "Simulation is paused   ");
        }
        glUniform1f(data->simplefont_shader_ypos_location, ypos++);
        j = convertLine(str,val);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(val), val);
        reb_glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, j);
        
        sprintf(str, "Press h for help ");
        glUniform1f(data->simplefont_shader_ypos_location, ypos++);
        j = convertLine(str,val);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(val), val);
        reb_glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, j);
        
        sprintf(str, "N = %d ",data->r_copy->N);
        glUniform1f(data->simplefont_shader_ypos_location, ypos++);
        j = convertLine(str,val);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(val), val);
        reb_glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, j);
        
        glUniform1f(data->simplefont_shader_ypos_location, ypos++);
        if (data->r_copy->integrator==REB_INTEGRATOR_SEI){
            sprintf(str, "t = %f [orb]  ", data->r_copy->t*data->r_copy->ri_sei.OMEGA/2./M_PI);
        }else{
            sprintf(str, "t = %f  ", data->r_copy->t);
        }
        j = convertLine(str,val);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(val), val);
        reb_glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, j);

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
        glUniform1f(data->simplefont_shader_screen_aspect_location, 1./ratio);
        glUniform1f(data->simplefont_shader_scale_location, 0.035);
        glBindBuffer(GL_ARRAY_BUFFER, data->simplefont_shader_charval_buffer);
        float val[200] = {0.};
        for (int i=0;i<sizeof(onscreenhelp)/sizeof(onscreenhelp[0]);i++){
            int j = convertLine(onscreenhelp[i],val);
            glUniform1f(data->simplefont_shader_ypos_location, (float)i);
            glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(val), val);
            reb_glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, j);
        }
    }
#endif // __EMSCRIPTEN__

    glfwSwapBuffers(data->window);
    glfwPollEvents();
}

#ifdef __EMSCRIPTEN__
EM_BOOL reb_render_frame_emscripten(double time, void* p){
    struct reb_simulation* r = (struct reb_simulation*)p;
    struct reb_display_data* data = r->display_data;
    if (!data){
        return EM_TRUE;
    }
    if (!data->pause){
        reb_render_frame(data);
        reb_overlay_hide(!data->onscreentext);
        if (data->onscreentext){ 
            char str[10240] = "\0";
            char line[1024];
            sprintf(line,"<div class=\"reboundlogo\"></div>REBOUND v%s<br />",reb_version_str);
            strlcat(str, line, 10240);
            if (data->connection_status>=0){
                if (data->r_copy->status == REB_STATUS_RUNNING){
                    sprintf(line, "Simulation is running<br />");
                }else if (data->r_copy->status == REB_STATUS_PAUSED){
                    sprintf(line, "Simulation is paused<br />");
                }else if (data->r_copy->status == REB_STATUS_SUCCESS){
                    sprintf(line, "Simulation ready<br />");
                }else if (data->r_copy->status == REB_STATUS_USER){
                    sprintf(line, "Simulation canceled<br />");
                }else if (data->r_copy->status > 0){
                    sprintf(line, "Simulation error occured<br />");
                }
                strlcat(str, line, 10240);
                sprintf(line, "N = %d<br />",data->r_copy->N);
                strlcat(str, line, 10240);
                sprintf(line, "t = %g<br />",data->r_copy->t);
                strlcat(str, line, 10240);
                sprintf(line, "steps/s = %g<br />",1./data->r_copy->walltime_last_steps);
                strlcat(str, line, 10240);
                strlcat(str, "Press h or click for help<br />", 10240);
                reb_overlay_update(str, data->r_copy->status);
            }else{
                sprintf(line, "Unable to connect. Server might have shut down.");
                strlcat(str, line, 10240);
                reb_overlay_update(str, 10);
            }
        }
        data->onscreenhelp = reb_overlay_help_show(data->onscreenhelp);
        if (data->onscreenhelp){ 
            char str[10240] = "\0";
            for (int i=0;i<sizeof(onscreenhelp)/sizeof(onscreenhelp[0]);i++){
                strlcat(str, onscreenhelp[i], 10240);
                strlcat(str, "<br />", 10240);
                reb_overlay_help_set_text(str);
            }
        }
    }
    return EM_TRUE;
}
#endif

void reb_display_init(struct reb_simulation * const r){
    struct reb_display_data* data = r->display_data;
    if (!glfwInit()){
        reb_simulation_error(r, "GLFW initialization failed.");
        return;
    }

    glfwSetErrorCallback(reb_glfw_error_callback);
#ifdef __EMSCRIPTEN__
    glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_ES_API);
#else
    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#endif

    GLFWwindow*  window = glfwCreateWindow(700, 700, "rebound", NULL, NULL);
    if (!window){
        reb_simulation_error(r,"GLFW window creation failed.");
        return;
    }

    glfwMakeContextCurrent(window);
    gladLoadGLLoader((GLADloadproc) glfwGetProcAddress);
    
    glfwSetWindowUserPointer(window,data); 

    // Default parameters
    reb_display_set_default_scale(r);
    { // Check if we have a retina display
        int wwidth, wheight, fwidth, fheight;
        glfwGetWindowSize(window, &wwidth, &wheight);
        glfwGetFramebufferSize(window, &fwidth, &fheight);
        data->retina = (double)fwidth/(double)wwidth;
    }
    if (data->r->max_radius0 > 0.0){
        data->spheres       = 1; 
    }else{
        data->spheres       = 0; 
    }
    data->pause         = 0; 
    data->multisample = 1; 
    if (data->r->integrator==REB_INTEGRATOR_WHFAST){
        data->wire          = 1; 
    }else{
        data->wire          = 0; 
    }
    data->onscreentext  = 1; 
    data->clear         = 1; 
    data->ghostboxes    = 0; 
    data->reference     = -1;
    data->view.r        = 1.;
    data->window        = window;

    glfwSetKeyCallback(window,reb_display_keyboard);
    glfwGetInputMode(window, GLFW_STICKY_MOUSE_BUTTONS);
    glfwSetMouseButtonCallback(window, reb_display_mouse_button);
    glfwSetCursorPosCallback(window, reb_display_cursor);
    glfwSetWindowSizeCallback(window, reb_display_resize);
    glDepthMask(GL_TRUE);
#ifndef __EMSCRIPTEN__
    glEnable(GL_MULTISAMPLE); 
#endif
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
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_R8, 256, 256, 0, GL_RED, GL_UNSIGNED_BYTE, image);
        glBindTexture(GL_TEXTURE_2D,0);
    }

    {
        const char* vertex_shader =
#ifdef __EMSCRIPTEN__
            "#version 300 es\n"
#else
            "#version 330\n"
#endif
            "in vec2 vp;\n"
            "in float charpos;\n"
            "uniform float scale;\n"
            "uniform float aspect;\n"
            "uniform float screen_aspect;\n"
            "uniform float ypos;\n"
            "uniform vec2 pos;\n"
            "in vec2 charval;\n"
            "in vec2 texcoord;\n"
            "out vec2 Texcoord;\n"
            "void main() {\n"
            "  gl_Position = vec4(pos.x*screen_aspect,pos.y,0.,0.)+vec4(screen_aspect*scale*(vp.x+charpos/aspect+charpos/16.),scale*(vp.y-ypos),0., 1.);\n"
            "  Texcoord = vec2((charval.s+texcoord.s)/16.,(charval.t+texcoord.t)/16.00);\n"
            "}\n";
        const char* fragment_shader =
#ifdef __EMSCRIPTEN__
            "#version 300 es\n"
#else
            "#version 330\n"
#endif
            "precision highp float;"
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
        data->simplefont_shader_screen_aspect_location = glGetUniformLocation(data->simplefont_shader_program, "screen_aspect");
    }
    
    {
        const char* vertex_shader =
#ifdef __EMSCRIPTEN__
            "#version 300 es\n"
#else
            "#version 330\n"
#endif
            "in vec3 vp;\n"
            "uniform mat4 mvp;\n"
            "uniform vec4 vc;\n"
            "out vec4 color;\n"
            "void main() {\n"
            "  gl_Position = mvp*vec4(vp, 1.0);\n"
            " gl_PointSize = 15.0f;\n"
            "  color = vc;\n"
            "}\n";
        const char* fragment_shader =
#ifdef __EMSCRIPTEN__
            "#version 300 es\n"
#else
            "#version 330\n"
#endif
            "precision highp float;"
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
#ifdef __EMSCRIPTEN__
            "#version 300 es\n"
#else
            "#version 330\n"
#endif
            "in vec3 vp;\n"
            "uniform mat4 mvp;\n"
            "uniform vec4 vc;\n"
            "out vec4 color;\n"
            "void main() {\n"
            "  gl_Position = mvp*vec4(vp, 1.0);\n"
            "  color = vc;\n"
            "}\n";
        const char* fragment_shader =
#ifdef __EMSCRIPTEN__
            "#version 300 es\n"
#else
            "#version 330\n"
#endif
            "precision highp float;"
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
#ifdef __EMSCRIPTEN__
            "#version 300 es\n"
#else
            "#version 330\n"
#endif
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
#ifdef __EMSCRIPTEN__
            "#version 300 es\n"
#else
            "#version 330\n"
#endif
            "precision highp float;"
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
#ifdef __EMSCRIPTEN__
            "#version 300 es\n"
#else
            "#version 330\n"
#endif
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
#ifdef __EMSCRIPTEN__
            "#version 300 es\n"
#else
            "#version 330\n"
#endif
            "precision highp float;"
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

    reb_glVertexAttribDivisor(sfvp, 0); 
    reb_glVertexAttribDivisor(sftexcoordp, 0); 
    reb_glVertexAttribDivisor(simplefont_shader_charpos_location,1);
    reb_glVertexAttribDivisor(simplefont_shader_charval_location,1);

    glBindVertexArray(0);


    // Particle data is dynamic
    glUseProgram(data->point_shader_program);
    glGenVertexArrays(1, &data->point_shader_particle_vao);
    glBindVertexArray(data->point_shader_particle_vao);
    GLuint pvp = glGetAttribLocation(data->point_shader_program,"vp");
    glEnableVertexAttribArray(pvp);

    GLuint particle_buffer;
    glGenBuffers(1, &particle_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, particle_buffer);
    data->particle_buffer = particle_buffer;
    
    glVertexAttribPointer(pvp, 3, GL_FLOAT, GL_FALSE, sizeof(float)*7, NULL);
    glBindVertexArray(0);
    
    // Create cross mesh
    glUseProgram(data->box_shader_program);
    glGenVertexArrays(1, &data->box_shader_cross_vao);
    glBindVertexArray(data->box_shader_cross_vao);
    GLuint cvp = glGetAttribLocation(data->box_shader_program,"vp");
    glEnableVertexAttribArray(cvp);
    
    float cross_data[18] = {
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
    glVertexAttribPointer(ssp, 3, GL_FLOAT, GL_FALSE, sizeof(float)*7, NULL);
    glVertexAttribPointer(ssr, 1, GL_FLOAT, GL_FALSE, sizeof(float)*7, (void*)(sizeof(float)*6));

    reb_glVertexAttribDivisor(svp, 0); 
    reb_glVertexAttribDivisor(data->sphere_shader_mvp_location, 0); 
    reb_glVertexAttribDivisor(ssp, 1);
    reb_glVertexAttribDivisor(ssr, 1);
    
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


    GLuint orbit_buffer;
    glGenBuffers(1, &orbit_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, orbit_buffer);
    data->orbit_buffer = orbit_buffer;
    glVertexAttribPointer(ofocusp, 3, GL_FLOAT, GL_FALSE, sizeof(float)*9, NULL);
    glVertexAttribPointer(oaefp, 3, GL_FLOAT, GL_FALSE, sizeof(float)*9, (void*)(sizeof(float)*3));
    glVertexAttribPointer(oomegaOmegaincp, 3, GL_FLOAT, GL_FALSE, sizeof(float)*9, (void*)(sizeof(float)*6));

    GLuint divisor = 0;
    reb_glVertexAttribDivisor(olintwopip, divisor); 
    reb_glVertexAttribDivisor(data->orbit_shader_mvp_location, 0); 
    reb_glVertexAttribDivisor(ofocusp, 1);
    reb_glVertexAttribDivisor(oaefp, 1);
    reb_glVertexAttribDivisor(oomegaOmegaincp, 1);
    
    glBindVertexArray(0);



    // Main display loop
#ifdef __EMSCRIPTEN__
    // Will return 
    emscripten_request_animation_frame_loop(reb_render_frame_emscripten, r);
#else
    glfwSwapInterval(1);

    while(!glfwWindowShouldClose(window) && r->status<0){
        double t0 = glfwGetTime();
        if (!data->pause){
            reb_render_frame(data);
        }
        while (glfwGetTime()-t0 < 1.0/120.) { // Maxframerate 120Hz
            usleep(10);
            glfwPollEvents();
        }
    }
    glfwDestroyWindow(window);
    data->window = NULL;
    // Destroy particle buffers
    if (data->N_allocated){
        data->N_allocated = 0;
        free(data->particle_data);
        free(data->orbit_data);
        data->particle_data = NULL;
        data->orbit_data = NULL;
    }
    glfwTerminate();
#endif
}

static void reb_display_set_default_scale(struct reb_simulation* const r){
    // Need a scale for visualization
    if (r->root_size==-1){  
        r->display_data->scale = 0.;
        const struct reb_particle* p = r->particles;
        for (unsigned int i=0;i<r->N-r->N_var;i++){
            const double _r = sqrt(p[i].x*p[i].x+p[i].y*p[i].y+p[i].z*p[i].z);
            r->display_data->scale = MAX(r->display_data->scale, _r);
        }
        if(r->display_data->scale==0.){
            r->display_data->scale = 1.;
        }
        r->display_data->scale *= 1.1;
    }else{
        r->display_data->scale = r->boxsize_max/2.;
    }
}

#endif // OPENGL


