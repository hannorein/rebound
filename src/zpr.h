#ifndef ZPR_H
#define ZPR_H

/*
 * Zoom-pan-rotate mouse manipulation module for GLUT
 * Version 0.4, October 2003
 *
 * Nigel Stewart
 * School of Computer Science and Information Technology
 * RMIT University
 * nigels@cs.rmit.edu.au
 *
 * Instructions
 * ------------
 *
 * Call zprInit() immediately after your call to glutCreateWindow()
 *
 * The ZPR module handles glutReshapeFunc(), glutMouseFunc() and glutMotionFunc()
 * Applications should not bypass the ZPR handlers for reshape or mouse events.
 *
 * Mouse manipulation of the GLUT window via the modelview matrix:
 *
 * Left   button -> rotate
 * Middle button -> zoom
 * Right  button -> pan
 *
 * Picking is also provided via two configurable callbacks:
 *
 * void zprSelectionFunc(void (*f)(void))
 *
 *   The draw function to be called in OpenGL selection
 *   mode in response to a mouse-down button event.
 *
 *
 * Limitations
 * -----------
 *
 * Works best with zprReferencePoint appropriately configured.
 * Works best with ortho projection.
 * You may need to use glEnable(GL_NORMALIZATION) for correct lighting.
 * Near and far clip planes are hard-coded.
 * Zooming and rotation is centered on the origin.
 * Only one window can use the callbacks at one time.
 *
 */


/* Mouse Manipulation API */

void zprReset();
void zprInit();

#endif

