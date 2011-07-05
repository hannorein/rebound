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

int display_boundaries = 1;

void displayKey(unsigned char key, int x, int y){
	switch(key){
		case 'q':
			exit(0);
			break;
		case 'b':
			display_boundaries = !display_boundaries;
			break;
	}
}

void display(){
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0,0,-boxsize);
	glPointSize(5.);
	glEnable(GL_POINT_SMOOTH);
	glColor4f(1.0,1.0,0.0,0.6);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_DOUBLE, 10*sizeof(double), particles);
	if (display_boundaries==1){
		for (int i=-nghostx;i<=nghostx;i++){
		for (int j=-nghosty;j<=nghosty;j++){
		for (int k=-nghostz;k<=nghostz;k++){
			struct ghostbox gb = get_ghostbox(i,j,k);
			glTranslatef(gb.shiftx,gb.shifty,gb.shiftz);
			glDrawArrays(GL_POINTS, 0, N);
			glTranslatef(-gb.shiftx,-gb.shifty,-gb.shiftz);
		}
		}
		}
	}else{
		glDrawArrays(GL_POINTS, 0, N);
	}
	glDisableClientState(GL_VERTEX_ARRAY);
	glColor4f(1.0,0.0,0.0,0.4);
	glutWireCube(boxsize);
	glutSwapBuffers();
	glMatrixMode(GL_PROJECTION);
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
	glutMainLoop();
}

#endif
