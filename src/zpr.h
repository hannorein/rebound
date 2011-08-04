/**
 * @file 	zpr.h
 * @brief 	Zoom-Pan-Rotate mouse manipulation module for GLUT.
 * @author 	Nigel Stewart
 * 		School of Computer Science and Information Technology
 * 		RMIT University
 * 		nigels@cs.rmit.edu.au
 * 
 * @section 	LICENSE
 * Copyright (c) 2003,2011 Nigel Stewart, Hanno Rein
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

#ifndef _ZPR_H
#define _ZPR_H

/**
 * Initialize the zoom/pan/rotating module.
 */
void zprInit();

/**
 * Reset the viewing angle, zooming scale.
 */
void zprReset();

#endif // _ZPR_H
