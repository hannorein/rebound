#ifdef OPENGL
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#ifdef _APPLE
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif 
#include "zpr.h"
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "tree.h"
#include "opengl.h"

GLuint DListSPHERE;
#ifndef COLLISIONS_NONE
int display_spheres = 1;
#else
int display_spheres = 0;
#endif

void displayKey(unsigned char key, int x, int y){
	switch(key){
		case 'q':
			exit(0);
			break;
		case 'b':
			printf("key caught\n");
			break;
		case 's':
			display_spheres = !display_spheres;
			break;
	}
}

void displayTree(struct cell *node){
	if (node == NULL) return;
	glTranslatef(node->x,node->y,node->z);
	glutWireCube(node->w);
	glTranslatef(-node->x,-node->y,-node->z);
	for (int i=0;i<8;i++) {
		displayTree(node->oct[i]);
	}
}

void display(){
	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	if (display_spheres){
		glDisable(GL_BLEND);                    
		glDepthMask(GL_TRUE);
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		GLfloat lightpos[] = {0, boxsize, boxsize, 0.f};
		glLightfv(GL_LIGHT0, GL_POSITION, lightpos);
	}else{
		glEnable(GL_BLEND);                    
		glDepthMask(GL_FALSE);
		glDisable(GL_DEPTH_TEST);
		glDisable(GL_LIGHTING);
		glDisable(GL_LIGHT0);
	}
	glTranslatef(0,0,-boxsize);
	glPointSize(5.);
	glEnable(GL_POINT_SMOOTH);
	glVertexPointer(3, GL_DOUBLE, sizeof(struct particle), particles);
	for (int i=-nghostx;i<=nghostx;i++){
	for (int j=-nghosty;j<=nghosty;j++){
	for (int k=-nghostz;k<=nghostz;k++){
		struct ghostbox gb = get_ghostbox(i,j,k);
		glTranslatef(gb.shiftx,gb.shifty,gb.shiftz);
		if (display_spheres){
			// Drawing Spheres
			glColor4f(1.0,1.0,1.0,1.0);
			for (int i=0;i<N;i++){
				struct particle p = particles[i];
				glTranslatef(p.x,p.y,p.z);
				glScalef(p.r,p.r,p.r);
				glCallList(DListSPHERE);
				glScalef(1./p.r,1./p.r,1./p.r);
				glTranslatef(-p.x,-p.y,-p.z);
			}
		}else{
			// Drawing Points
			glColor4f(1.0,1.0,0.0,0.6);
			glEnableClientState(GL_VERTEX_ARRAY);
			glDrawArrays(GL_POINTS, 0, N);
			glDisableClientState(GL_VERTEX_ARRAY);
		}
		// Drawing Tree
		glColor4f(1.0,0.0,0.0,0.4);
		displayTree(root);
		glTranslatef(-gb.shiftx,-gb.shifty,-gb.shiftz);
	}
	}
	}
	glColor4f(1.0,0.0,0.0,0.4);
	glutWireCube(boxsize);
	glutSwapBuffers();
	glTranslatef(0,0,boxsize);
}

void init_display(int argc, char* argv[]){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );
	glutInitWindowSize(700,700);
	glutCreateWindow("nbody");
	zprInit(0.7/boxsize);
	glutDisplayFunc(display);
	glutIdleFunc(iterate);
	glutKeyboardFunc(displayKey);
	glEnable(GL_BLEND);                    
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);  
	
	// Sphere
	DListSPHERE = glGenLists(1);
	GLUquadricObj *sphere;
	glNewList(DListSPHERE, GL_COMPILE);
	sphere = gluNewQuadric();
	gluSphere(sphere, 1.f, 20, 20);
	gluDeleteQuadric(sphere);
	glEndList();
  	
	// Light

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
	
	glutMainLoop();
}

#endif
