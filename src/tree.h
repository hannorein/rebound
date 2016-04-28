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

struct reb_treecell; 

/**
 * @brief The data structure of one node of a tree 
 */
struct reb_treecell {
	double x; /**< The x position of the center of a cell */
	double y; /**< The y position of the center of a cell */
	double z; /**< The z position of the center of a cell */
	double w; /**< The width of a cell */
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
	struct reb_treecell *oct[8]; /**< The pointer array to the octants of a cell */
	int pt;		/**< It has double usages: in a leaf node, it stores the index 
			  * of a particle; in a non-leaf node, it equals to (-1)*Total 
			  * Number of particles within that cell. */ 
};

/**
  * @brief This function updates the tree.
  * @details The tree needs to be updated when particles move, this function does that.
  * @param r Rebound simulation to operate on
  */
void reb_tree_update(struct reb_simulation* const r);

/**
  * @brief The wrap function calls reb_tree_update_gravity_data_in_cell() for each tree.
  * @param r Rebound simulation to operate on
  */
void reb_tree_update_gravity_data(struct reb_simulation* const r);

/**
  * @brief The wrap function calls reb_tree_add_particle_to_cell() to add the particle into one of the trees. If the tree_root doesn't exist, then it initializes the tree. 
  * @param r Rebound simulation to operate on
  * @param pt Index of a particle.
  */
void reb_tree_add_particle_to_tree(struct reb_simulation* const r, int pt);

/**
 * @brief Free up all space occupied by the tree structure.
 * This will not modify particles.
  * @param r Rebound simulation to operate on
 */
void reb_tree_delete(struct reb_simulation* const r);

#ifdef MPI
/**
  * @brief MPI related function used to calculate gravity from nearby nodes
  * @param node is a pointer to a node cell.
  */
void reb_tree_add_essential_node(struct reb_simulation* const r, struct reb_treecell* node);
/**
  * @brief MPI related function used to calculate gravity from nearby nodes
  */
void reb_tree_prepare_essential_tree_for_gravity(struct reb_simulation* const r);
/**
  * @brief MPI related function used to calculate gravity from nearby nodes
  */
void reb_tree_prepare_essential_tree_for_collisions(struct reb_simulation* const r);
#endif // MPI

#endif // _TREE_H
