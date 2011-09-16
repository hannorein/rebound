/**
 * @file 	tree.h
 * @brief 	Tree header file. All the tree implementations needed by other routines.
 * @author 	Shangfei Liu <liushangfei@pku.edu.cn>, Hanno Rein <hanno@hanno-rein.de>
 *
 * @section 	LICENSE
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

#ifndef _TREE_H
#define _TREE_H
#if defined(GRAVITY_TREE) || defined(COLLISIONS_TREE)
#define TREE

struct cell; /**< The data structure of one node of a tree */

struct cell {
	double x; /**< The x position of the center of a cell */
	double y; /**< The y position of the center of a cell */
	double z; /**< The z position of the center of a cell */
	double w; /**< The width of a cell */
#ifdef GRAVITY_TREE
	double m; /**< The total mass of a cell */
	double mx; /**< The x position of the center of mass of a cell */
	double my; /**< The y position of the center of mass of a cell */
	double mz; /**< The z position of the center of mass of a cell */
#ifdef QUADRUPOLE
	double mxx; /**< The xx component of the quadrupole tensor of mass of a cell */
	double mxy; /**< The xy component of the quadrupole tensor of mass of a cell */
	double mxz; /**< The xz component of the quadrupole tensor of mass of a cell */
	double myy; /**< The yy component of the quadrupole tensor of mass of a cell */
	double myz; /**< The yz component of the quadrupole tensor of mass of a cell */
	double mzz; /**< The zz component of the quadrupole tensor of mass of a cell */
#endif // QUADRUPOLE
#endif // GRAVITY_TREE
	struct cell *oct[8]; /**< The pointer array to the octants of a cell */
	int pt;	/**< It has double usages: in a leaf node, it stores the index 
			  * of a particle; in a non-leaf node, it equals to (-1)*Total 
			  * Number of particles within that cell. */ 
};

extern struct cell** tree_root; /**< A public pointer to the roots of the trees. */

/**
  * The wrap function corresponds to initializing the trees when they don't exist and updating the structures of the trees by calling tree_update_cell. 
  */
void tree_update();

#ifdef GRAVITY_TREE
/**
  * The wrap function calls tree_update_gravity_data_in_cell() to for each tree.
  */
void tree_update_gravity_data();
#endif // GRAVITY_TREE

/**
  * The wrap function calls tree_add_particle_to_cell() to add the particle into one of the trees. If the tree_root doesn't exist, then it initializes the trees. 
  *
  * @param pt is the index of a particle.
  */
void tree_add_particle_to_tree(int pt);

#ifdef MPI
/**
  * Needs more comments!
  *
  * @param node is a pointer to a node cell.
  */
void tree_add_essential_node(struct cell* node);
#ifdef GRAVITY_TREE
/**
  * Needs more comments!
  */
void tree_prepare_essential_tree_for_gravity();
#endif //GRAVITY_TREE
#ifdef COLLISIONS_TREE
/**
  * Needs more comments!
  */
void tree_prepare_essential_tree_for_collisions();
#endif //COLLISIONS_TREE
#endif // MPI

/**
 * Particle between 0 and N_tree_fixed will not be shuffled around during tree-reconstruction.
 */
extern int N_tree_fixed;

#endif // TREE
#endif // _TREE_H
