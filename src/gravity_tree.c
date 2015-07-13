/**
 * @file 	gravity.c
 * @brief 	Gravity calculation using an oct-tree, O(N log(N)).
 * @author 	Hanno Rein <hanno@hanno-rein.de>, Shangfei Liu <liushangfei@pku.edu.cn>
 *
 * @details 	The routines in this file implement a gravity calculation 
 * using the oct-tree defined in file tree.h. It can be run with
 * MPI using a distributed tree. In that case the locally essential tree 
 * is shared with every other node. The method scales as O(N log(N)) for 
 * large particles. For small particle numbers, a direct summation
 * might be faster, as it avoids having the overhead of a * complicated 
 * data structure. 
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
#include "rebound.h"
#include "tree.h"
#include "boundary.h"


