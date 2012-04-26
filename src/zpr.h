#ifndef ZPR_H
/**
 * @file 	zpr.h
 * @brief 	Zoom-Pan-Rotate mouse manipulation module for GLUT.
 * @author 	Nigel Stewart
 * 		School of Computer Science and Information Technology
 * 		RMIT University
 * 		nigels@cs.rmit.edu.au
 * 
 * @section 	LICENSE
 * Copyright (c) 2003,2012 Nigel Stewart, Hanno Rein
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
 *
 */
#define ZPR_H


#ifdef OPENGL
#ifdef _APPLE
#include <GLUT/glut.h>
#else // _APPLE
#include <GL/glut.h>
#endif // _APPLE

/* Mouse Manipulation API */

void zprInit();
extern GLfloat zprReferencePoint[4];

/**
 * Reset the viewing angle, zooming scale.
 */
void zprReset();


/* Picking API (Optional) */

extern void zprSelectionFunc(void (*f)(void));      /* Selection-mode draw function */
extern void zprPickFunc(void (*f)(GLint name));     /* Pick event handling function */

#endif // OPENGL
#endif
