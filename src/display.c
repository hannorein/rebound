/**
 * @file 	display.c
 * @brief 	Realtime OpenGL visualization.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @details 	These functions provide real time visualizations
 * using OpenGL. 
 * 
 * @section 	LICENSE
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
#include <sys/stat.h>
#include <sys/time.h>
#include <semaphore.h>
#ifdef _APPLE
#include <GLUT/glut.h>
#else // _APPLE
#include <GL/glut.h>
#endif // _APPLE
#include "rebound.h"
#include "tools.h"
#include "zpr.h"
#include "particle.h"
#include "boundary.h"
#include "display.h"
#include "output.h"
#include "integrator.h"
#define WINWIDTH 700
#define WINHEIGHT 700


struct reb_display_config {
	int spheres;	/**< Switches between point sprite and real spheres. */
	int pause;	/**< Pauses visualization, but keep simulation running */
	int wire;	/**< Shows/hides orbit wires. */
	int onscreentext;	/**< Shows/hides onscreen text. */
	int clear;	/**< Toggles clearing the display on each draw. */
	int ghostboxes;	/**< Shows/hides ghost boxes. */
	int reference;	/**< reb_particle used as a reference for centering. */
	struct reb_simulation* r;	/**< Simulation to render */
	sem_t* mutex;			/**< Mutex to guarantee non-flickering */
#ifdef _APPLE
	GLuint dlist_sphere;		/**< Precalculated display list of a sphere. */
#endif // APPLE
};

struct reb_display_config reb_dc;

void reb_display_exit(struct reb_display_config* dc){
	dc->r = NULL;
#ifdef _APPLE
	glDeleteLists(dc->dlist_sphere,1);
#endif // _APPLE
	printf("Exiting vizualization.\n");
	exit(0);
}

void reb_display_timer(int value){
	if (reb_dc.r->status>=0){
		reb_display_exit(&reb_dc);
	}else{
		glutPostRedisplay();
	}
	glutTimerFunc(20,reb_display_timer,0); // 50 Hz refresh rate.
}

void reb_display(void){
	if (reb_dc.pause){
		return;
	}
	sem_wait(reb_dc.mutex);	
	const struct reb_particle* particles = reb_dc.r->particles;
	
	if (reb_dc.clear){
	        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	}
	glEnable(GL_POINT_SMOOTH);
	glVertexPointer(3, GL_DOUBLE, sizeof(struct reb_particle), particles);
	//int _N_active = ((N_active==-1)?N:N_active);
	if (reb_dc.reference>=0){
		glTranslatef(-particles[reb_dc.reference].x,-particles[reb_dc.reference].y,-particles[reb_dc.reference].z);
	}
	for (int i=-reb_dc.ghostboxes*reb_dc.r->nghostx;i<=reb_dc.ghostboxes*reb_dc.r->nghostx;i++){
	for (int j=-reb_dc.ghostboxes*reb_dc.r->nghosty;j<=reb_dc.ghostboxes*reb_dc.r->nghosty;j++){
	for (int k=-reb_dc.ghostboxes*reb_dc.r->nghostz;k<=reb_dc.ghostboxes*reb_dc.r->nghostz;k++){
		struct reb_ghostbox gb = reb_boundary_get_ghostbox(reb_dc.r, i,j,k);
		glTranslatef(gb.shiftx,gb.shifty,gb.shiftz);
		if (!(!reb_dc.clear&&reb_dc.wire)){
			if (reb_dc.spheres==0 || reb_dc.spheres==2){
				// Drawing Points
				glEnableClientState(GL_VERTEX_ARRAY);
				glPointSize(3.);
				glColor4f(1.0,1.0,1.0,0.5);
				//glDrawArrays(GL_POINTS, _N_active, N-_N_active);
				glColor4f(1.0,1.0,0.0,0.9);
				glPointSize(5.);
				glDrawArrays(GL_POINTS, 0, reb_dc.r->N-reb_dc.r->N_var);
				glDisableClientState(GL_VERTEX_ARRAY);
			}
			if (reb_dc.spheres){
				glDisable(GL_BLEND);                    
				glEnable(GL_DEPTH_TEST);
				glEnable(GL_LIGHTING);
				glEnable(GL_LIGHT0);
				GLfloat lightpos[] = {0, reb_dc.r->boxsize_max, reb_dc.r->boxsize_max, 0.f};
				glLightfv(GL_LIGHT0, GL_POSITION, lightpos);
				// Drawing Spheres
				glColor4f(1.0,1.0,1.0,1.0);
				for (int i=0;i<reb_dc.r->N-reb_dc.r->N_var;i++){
					struct reb_particle p = particles[i];
					if (p.r>0){
						glTranslatef(p.x,p.y,p.z);
						glScalef(p.r,p.r,p.r);
#ifdef _APPLE
						glCallList(reb_dc.dlist_sphere);
#else //_APPLE
						glutSolidSphere(1,40,10);
#endif //_APPLE
						glScalef(1./p.r,1./p.r,1./p.r);
						glTranslatef(-p.x,-p.y,-p.z);
					}
				}
				glEnable(GL_BLEND);                    
				glDisable(GL_DEPTH_TEST);
				glDisable(GL_LIGHTING);
				glDisable(GL_LIGHT0);
			}
		}
		// Drawing wires
		if (reb_dc.wire){
			if(reb_dc.r->integrator!=REB_INTEGRATOR_SEI){
				double radius = 0;
				struct reb_particle com = particles[0];
				for (int i=1;i<reb_dc.r->N-reb_dc.r->N_var;i++){
					struct reb_particle p = particles[i];
					if (reb_dc.r->N_active>0){
						// Different colors for active/test particles
						if (i>=reb_dc.r->N_active){
							glColor4f(0.9,1.0,0.9,0.9);
						}else{
							glColor4f(1.0,0.9,0.0,0.9);
						}
					}else{
						// Alternating colors
						if (i%2 == 1){
							glColor4f(0.0,1.0,0.0,0.9);
						}else{
							glColor4f(0.0,0.0,1.0,0.9);
						}
					}
					//if (reb_dc.r->integrator==REB_INTEGRATOR_WHFAST && reb_dc.r->ri_whfast.is_synchronized==0){
					//	double m = p.m;
					//	p = reb_dc.r->ri_whfast.p_j[i];
					//	p.m = m;
					//}
					struct reb_orbit o = reb_tools_particle_to_orbit(reb_dc.r->G, p,com);
					glPushMatrix();
					
					glTranslatef(com.x,com.y,com.z);
					glRotatef(o.Omega/DEG2RAD,0,0,1);
					glRotatef(o.inc/DEG2RAD,1,0,0);
					glRotatef(o.omega/DEG2RAD,0,0,1);
					
					glBegin(GL_LINE_LOOP);
					for (double trueAnom=0; trueAnom < 2.*M_PI; trueAnom+=M_PI/100.) {
						//convert degrees into radians
						radius = o.a * (1. - o.e*o.e) / (1. + o.e*cos(trueAnom));
						glVertex3f(radius*cos(trueAnom),radius*sin(trueAnom),0);
					}
					glEnd();
					glPopMatrix();
					com = reb_get_com_of_pair(p,com);
				}
			}else{
				for (int i=1;i<reb_dc.r->N;i++){
					struct reb_particle p = particles[i];
					glBegin(GL_LINE_LOOP);
					for (double _t=-100.*reb_dc.r->dt;_t<=100.*reb_dc.r->dt;_t+=20.*reb_dc.r->dt){
						double frac = 1.-fabs(_t/(120.*reb_dc.r->dt));
						glColor4f(1.0,(_t+100.*reb_dc.r->dt)/(200.*reb_dc.r->dt),0.0,frac);
						glVertex3f(p.x+p.vx*_t, p.y+p.vy*_t, p.z+p.vz*_t);
					}
					glEnd();
				}
			}
		}
		glTranslatef(-gb.shiftx,-gb.shifty,-gb.shiftz);
	}
	}
	}
	glColor4f(1.0,0.0,0.0,0.4);
	glScalef(reb_dc.r->boxsize.x,reb_dc.r->boxsize.y,reb_dc.r->boxsize.z);
	if (reb_dc.r->boundary == REB_BOUNDARY_NONE){
		glBegin(GL_LINES);
		glVertex3f(0,0,0.04);
		glVertex3f(0,0,-0.04);
		glVertex3f(0,0.04,0);
		glVertex3f(0,-0.04,0);
		glVertex3f(0.04,0,0);
		glVertex3f(-0.04,0,0);
		glEnd();

	}else{
		glutWireCube(1);
	}
	glScalef(1./reb_dc.r->boxsize.x,1./reb_dc.r->boxsize.y,1./reb_dc.r->boxsize.z);
	if (reb_dc.reference>=0){
		glTranslatef(particles[reb_dc.reference].x,particles[reb_dc.reference].y,particles[reb_dc.reference].z);
	}
	
    if (reb_dc.onscreentext){
        glMatrixMode(GL_PROJECTION);
        glPushMatrix(); // save
        glLoadIdentity();// and clear
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();
        glColor4f(1.0,1.0,1.0,0.5);
        glRasterPos2f(-0.98,-0.98);
        
        char str[4096] = "\0";
        sprintf(str, "REBOUND");
        if (reb_dc.r->status == REB_RUNNING){
            sprintf(str, "%s (running)  ", str);
        }else if (reb_dc.r->status == REB_RUNNING_PAUSED){
            sprintf(str, "%s (paused)   ", str);
        }
        sprintf(str, "%s  N_tot= %d  ",str, reb_dc.r->N);
        if (reb_dc.r->integrator==REB_INTEGRATOR_SEI){
            sprintf(str, "%st= %f [orb]  ",str, reb_dc.r->t*reb_dc.r->ri_sei.OMEGA/2./M_PI);
        }else{
            sprintf(str, "%st= %f  ",str, reb_dc.r->t);
        }
        sprintf(str,"%sdt= %f  ",str,reb_dc.r->dt);
        if (reb_dc.r->integrator==REB_INTEGRATOR_HYBRID){
            sprintf(str, "%s INT= %- 1d  ", str, reb_dc.r->ri_hybrid.mode);
        }
        
        const char* p = str;
        do glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *p); while(*(++p));
        
        glMatrixMode(GL_PROJECTION);
        glPopMatrix(); // revert back to the matrix
        glMatrixMode(GL_MODELVIEW );
        glPopMatrix();
    }

	glFlush();
	sem_post(reb_dc.mutex);	
}

void reb_display_keyboard(unsigned char key, int x, int y){
	switch(key){
		case 'q': case 'Q':
			reb_dc.r->status = REB_EXIT_USER;
			reb_display_exit(&reb_dc);
			break;
		case ' ':
			if (reb_dc.r->status == REB_RUNNING_PAUSED){
				printf("Resume.\n");
				reb_dc.r->status = REB_RUNNING;
			}else{
				printf("Pause.\n");
				reb_dc.r->status = REB_RUNNING_PAUSED;
			}
			break;
		case 's': case 'S':
			reb_dc.spheres = (reb_dc.spheres+1)%3;
			break;
		case 'g': case 'G':
			reb_dc.ghostboxes = !reb_dc.ghostboxes;
			break;
		case 'r': case 'R':
			zprReset(reb_dc.r->boxsize_max);
			break;
		case 'd': case 'D':
			reb_dc.pause = !reb_dc.pause;
			break;
		case 'w': case 'W':
			reb_dc.wire = !reb_dc.wire;
			break;
		case 't': case 'T':
			reb_dc.onscreentext = !reb_dc.onscreentext;
			break;
		case 'c': case 'C':
			reb_dc.clear = !reb_dc.clear;
			break;
		case 'x': 
			reb_dc.reference++;
			if (reb_dc.reference>reb_dc.r->N) reb_dc.reference = -1;
			printf("Reference particle: %d.\n",reb_dc.reference);
			break;
		case 'X': 
			reb_dc.reference--;
			if (reb_dc.reference<-1) reb_dc.reference = reb_dc.r->N-1;
			printf("Reference particle: %d.\n",reb_dc.reference);
			break;
	}
	reb_display();
}


void reb_display_init(int argc, char* argv[], struct reb_simulation* r, sem_t* mutex){
	reb_dc.r 			= r;
	reb_dc.mutex 		= mutex;
	// Default parameters
	reb_dc.spheres 		= 2; 
	reb_dc.pause 		= 0; 
	reb_dc.wire 		= 0; 
	reb_dc.onscreentext = 1; 
	reb_dc.clear 		= 1; 
	reb_dc.ghostboxes 	= 0; 
	reb_dc.reference 	= -1;

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH );
	glutInitWindowSize(WINWIDTH,WINHEIGHT);
	glutCreateWindow("rebound");
	zprInit(reb_dc.r->boxsize_max);
	glutDisplayFunc(reb_display);
	glutKeyboardFunc(reb_display_keyboard);
	glDepthMask(GL_TRUE);
	glEnable(GL_BLEND);                    
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);  
	
	// Sphere
#ifdef _APPLE
	reb_dc.dlist_sphere = glGenLists(1);
	GLUquadricObj *sphere;
	glNewList(reb_dc.dlist_sphere, GL_COMPILE);
	sphere = gluNewQuadric();
	gluSphere(sphere, 1.f, 20, 20);
	gluDeleteQuadric(sphere);
	glEndList();
#endif // _APPLE
  	
	// Setup lights

	glCullFace(GL_BACK);
	glShadeModel ( GL_SMOOTH );
	glEnable( GL_NORMALIZE );
	glEnable(GL_COLOR_MATERIAL);
	static GLfloat light[] = {0.7f, 0.7f, 0.7f, 1.f};
	static GLfloat lightspec[] = {0.2f, 0.2f, 0.2f, 1.f};
	static GLfloat lmodel_ambient[] = { 0.15, 0.14, 0.13, 1.0 };

	glLightfv(GL_LIGHT0, GL_DIFFUSE, light );
	glLightfv(GL_LIGHT0, GL_SPECULAR, lightspec );
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);

	static GLfloat sphere_mat[] = {0.8f, 0.8f, 0.8f, 1.f};
	static GLfloat sphere_spec[] = {1.0f, 1.0f, 1.0f, 1.f};
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, sphere_mat);
	glMaterialfv(GL_FRONT, GL_SPECULAR, sphere_spec);
	glMaterialf(GL_FRONT, GL_SHININESS, 80);

	glutTimerFunc(100,reb_display_timer,0);
	// Enter glut run loop and never come back.
	glutMainLoop();
}

#endif // OPENGL
