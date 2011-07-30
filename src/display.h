/**
 * @file 	display.h
 * @brief 	Realtime OpenGL visualization.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
 *
 * This file is part of nbody.
 *
 * nbody is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * nbody is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with nbody.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#ifndef _DISPLAY_H
#define _DISPLAY_H

/**
 * This routine is called by the glut run loop or manually, whenever the screen has to be redrawn 
 */
void display();

/**
 * This function initializes OpenGL and starts the run loop. It will never return.
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 */
void display_init(int argc, char* argv[]);

extern int display_init_done;	/**< Is set to one when the display is initialized and can be drawn. This prevents errors when output_png() is called, but display not initialized yet. */
#endif
