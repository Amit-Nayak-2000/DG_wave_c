#include "dg_hrefinement.h"
#include "dg_unit.h"
#include "dg_local_storage.h"
#include <cstdlib>	// random number
#include "dg_cantor_pairing.h"	// for the key
#include <cassert>
#include <unordered_map>
#include "dg_status_table.h"
#include "dg_param.h"
#include <array>
#include <cassert>
#include "dg_elem_length.h"
#include "dg_neighbour_list.h"
#include <algorithm>
#include <iostream>	// test

/// global variable
double xcoord_new[3]{};
double ycoord_new[3]{};

/// forward declaration ----------------------------------------------------------------------------
void Get_coordinates(int ith, double* xcoord, double* ycoord);
void Gen_index(int ith, int* index, int* index_new);
void Two_siblings(int new_key, int position);
void Four_children(int ith, int i, int j, int level, std::array<int, 4>& children);
void Non_sibling_interfaces(Unit* last, int parent);
void Form_one_direction(int key1, int key2, int parent, int facen);
void Match_neighbours(int parent, int local_key, int facen, std::vector<int>& neighbours);
void Flag_elem();
// --------------------------------------------------------------------------------------------------

/// @brief
/// loop through all the elements and flag then for refining or coarsening. 
void Flag_elem(){

	Unit* temp = local::head;

	for(int k = 0; k < local::local_elem_num; ++k){

		// generate random number
		int rand_num = rand() % 10 + 1;	// random number between [1, 10]

		bool check_h = ((temp -> index[2]) < grid::hlevel_max ) ? true : false;
		bool check_c = ((temp -> index[2]) > 0)	 ? true : false;


		if(check_h && rand_num <= 3){	// h-refinement

			temp -> hrefine = true;
		}
		else if(check_c && rand_num > 7 ){	// coarse
			
			temp -> coarsen = true;
		}
		temp = temp -> next;
	}
}

/// @brief
/// Random h-refinement scheme. Each element has 30% chance to split.
void h_refinement(){

	Unit* temp = local::head;
	Unit* temp2 = temp;

	int increment{};
	
	for(int k = 0; k < local::local_elem_num; ++k){

		if(temp -> hrefine){	// refine
			
			increment += 3; 
			
			// new coordinates
			xcoord_new[0] = temp -> xcoords[0]; xcoord_new[2] = temp -> xcoords[1];
			ycoord_new[0] = temp -> ycoords[0]; ycoord_new[2] = temp -> ycoords[1];
			xcoord_new[1] = ((temp -> xcoords[0]) + (temp -> xcoords[1])) / 2.0;
			ycoord_new[1] = ((temp -> ycoords[0]) + (temp -> ycoords[1])) / 2.0;
			
			// current key
			int old_key = Get_key_fun(temp->index[0], temp->index[1], temp->index[2]);
		
			// previous key between 4 children	
			int pre_key{};

			// create 4 unit in hash table
			for(int i = 0; i < 4; ++i){
				

				int position = Sibling_position(temp->status, i);
				
				int index_new[3]{};
				Gen_index(position, temp->index, index_new);
				
				int new_key = Get_key_fun(index_new[0], index_new[1], index_new[2]);

				// create unit
				local::Hash_elem[new_key] = new Unit();

				if(i == 0 && k == 0){	// first elem, update head
					local::head = local::Hash_elem[new_key];
					
				}
				else if(i == 0 && k != 0){
					temp2 -> next = local::Hash_elem[new_key];

				}
				
				// unit index
				local::Hash_elem[new_key] -> index[0] = index_new[0];
				local::Hash_elem[new_key] -> index[1] = index_new[1];
				local::Hash_elem[new_key] -> index[2] = index_new[2];

				// status
				local::Hash_elem[new_key] -> status = Status_table(temp->status, i );

				// coordinates
				Get_coordinates(position, local::Hash_elem[new_key] -> xcoords,
						 local::Hash_elem[new_key] -> ycoords);


				// ith child
				local::Hash_elem[new_key] -> child_position = position;

				// poly order
				local::Hash_elem[new_key] -> n = temp -> n;
				local::Hash_elem[new_key] -> m = temp -> m;
				
				// faces(here we construct between siblings)
				Two_siblings(new_key, position);
					

				// interpolate solutions
				local::Hash_elem[new_key] -> var = temp -> var;
				

				// form link
				if(i > 0){
					local::Hash_elem[pre_key] -> next = local::Hash_elem[new_key];
				}

				pre_key = new_key;
				
				temp2 = local::Hash_elem[new_key];
				
			}	
			// form then external interfaces between 4 children and updates their neighbour's faces			
			Non_sibling_interfaces(temp2, old_key);
			temp2 -> next = temp -> next;

			// erase the parent
			local::Hash_elem.erase(old_key);
			temp = temp2;
		}
		temp2 = temp;
		temp = temp -> next;
	}

	local::local_elem_num += increment;

}

/// @brief
/// Form the non-sibling interfaces of 4 children. 
/// @param last pointer to the 3th child.
/// @param parent parent key.
void Non_sibling_interfaces(Unit* last, int parent){

	// get four children's keys
	std::array<int, 4> children; 
	Four_children(last -> child_position, last -> index[0], last -> index[1], last -> index[2], children);

	// south
	Form_one_direction(children[0], children[3], parent, 0);
	
	// north
	Form_one_direction(children[1], children[2], parent, 1);
		
	// west
	Form_one_direction(children[0], children[1], parent, 2);
	
	// east
	Form_one_direction(children[3], children[2], parent, 3);
}

/// @brief
/// Loop through parent's face and sign out the nieghbours to the two children on one face.
/// Plus update the neighbours face info.
/// @param parent parent's key.
/// @param local_key child's key.
/// @parem facen face number.
/// @param neighbours vector of the possible neighbours.
void Match_neighbours(int parent, int local_key, int facen, std::vector<int>& neighbours){

	int l_tol{};
	int local_length = Elem_length(local::Hash_elem[local_key] -> index[2]);

	for(auto it = local::Hash_elem[parent] -> facen[facen].begin(); 
		it != local::Hash_elem[parent] -> facen[facen].end(); ++it){

		auto it_n = std::find(neighbours.begin(), neighbours.end(), it -> key);

		if(it_n != neighbours.end()){	// if find, stores
			
			Unit::Face obj = {it -> face_type, it -> hlevel, it -> porderx, it -> pordery, 
						it -> key, it -> rank};
				
			local::Hash_elem[local_key] -> facen[facen].emplace_back(obj);

			// update neighbours face info
			if(it -> face_type == 'L'){	// only local element
				
				int n_dir = Opposite_dir(facen);	// neighbour face direction
				int n_key = it -> key;	// neighbour's key

				// erase parent info, if exist
				for(auto it_nf = local::Hash_elem[n_key] -> facen[n_dir].begin();
					it_nf != local::Hash_elem[n_key] -> facen[n_dir].end(); ++it_nf){

					if(it_nf -> key == parent){	// find old parent info
				
						it_nf = local::Hash_elem[n_key] -> facen[n_dir].erase(it_nf);
						break;
					}

				}
				// child info
				Unit::Face obj2 = {'L', local::Hash_elem[local_key] -> index[2], 
							local::Hash_elem[local_key] -> n, 
							local::Hash_elem[local_key] -> m, 
							local_key, mpi::rank};
				local::Hash_elem[n_key] -> facen[n_dir].emplace_back(obj2);
			}
			

			l_tol += Elem_length(it -> hlevel);

			while(l_tol >= local_length){ break; }
		}
	
	}
	
}

/// @brief 
/// Form one direction's non-siblnig interfaces of two children.
/// @param key1 1st child's key
/// @param key2 2nd child's key (two children in index i/j ascending mode).
/// @param parent parent's key.
/// @param facen face direction. 
void Form_one_direction(int key1, int key2, int parent, int facen){
	
	// if on the physical boundary, only inherit. 
	if(local::Hash_elem[parent] -> facen[facen].front().face_type == 'B'){

		local::Hash_elem[key1] -> facen[facen].push_back(Unit::Face());
		local::Hash_elem[key1] -> facen[facen].front().face_type = 'B';
	
		local::Hash_elem[key2] -> facen[facen].push_back(Unit::Face());
		local::Hash_elem[key2] -> facen[facen].front().face_type = 'B';

		return; 
	}

	std::array<int, 2> two = {key1, key2};
	
	std::unordered_map<int, std::vector<int>> neighbours;

	// not on the physical boundary, form possible nieghbour list and search
	if(facen == 0){	// south

		for(auto& key : two){
	
			int i = local::Hash_elem[key] -> index[0];
			int j = local::Hash_elem[key] -> index[1];
			int k = local::Hash_elem[key] -> index[2];

			// form neighbour list
			Neighbours_array_x(i - 1, j, k, key, facen, neighbours);
			
			// form child's neighbours + updates neighbours' interfaces			
			Match_neighbours(parent, key, facen, neighbours[key]);
		}
	}
	else if(facen == 1){	// north 

		for(auto& key : two){
	
			int i = local::Hash_elem[key] -> index[0];
			int j = local::Hash_elem[key] -> index[1];
			int k = local::Hash_elem[key] -> index[2];

			// form neighbour list
			Neighbours_array_x(i + 1, j, k, key, facen, neighbours);
			
			// form child's neighbours + updates neighbours' interfaces			
			Match_neighbours(parent, key, facen, neighbours[key]);
		}

	}
	else if(facen == 2){	// west

		for(auto& key : two){
	
			int i = local::Hash_elem[key] -> index[0];
			int j = local::Hash_elem[key] -> index[1];
			int k = local::Hash_elem[key] -> index[2];

			// form neighbour list
			Neighbours_array_y(i, j - 1, k, key, facen, neighbours);
			
			// form child's neighbours + updates neighbours' interfaces			
			Match_neighbours(parent, key, facen, neighbours[key]);
		}

	}
	else{	// east

		for(auto& key : two){
	
			int i = local::Hash_elem[key] -> index[0];
			int j = local::Hash_elem[key] -> index[1];
			int k = local::Hash_elem[key] -> index[2];

			// form neighbour list
			Neighbours_array_y(i, j + 1, k, key, facen, neighbours);
			
			// form child's neighbours + updates neighbours' interfaces			
			Match_neighbours(parent, key, facen, neighbours[key]);
		}

	}

}


/// @brief
/// Input the last child's element index, then returns the array 
/// of 4 children's key in ascending child_position sequence.
/// @param ith ith children (only support 1 and 3).
/// @param i integer coordinate in x direction.
/// @param j integer coordinate in y direction.
/// @param level element hlevel.
/// @param children children key array. 
void Four_children(int ith, int i, int j, int level, std::array<int, 4>& children){

	assert((ith == 1 || ith == 3) && "Child_position must be 1 or 3");

	if(ith == 1){	// 1st child
	
		children[0] = Get_key_fun(i - 1, j, level);
		children[1] = Get_key_fun(i, j, level);
		children[2] = Get_key_fun(i, j + 1, level);
		children[3] = Get_key_fun(i - 1, j + 1, level);

	}
	else{
		children[0] = Get_key_fun(i , j - 1, level);
		children[1] = Get_key_fun(i + 1, j - 1, level);
		children[2] = Get_key_fun(i + 1, j, level);
		children[3] = Get_key_fun(i, j, level);
		
	}


}


/// @brief
/// Build faces between siblngs. Each element has two faces adjecant to siblings. 
/// @param new_key current child element's key
/// @param position ith child. 
void Two_siblings(int new_key, int position){
	
	// four faces
	for(int i = 0; i < 4; ++i){

		bool t = Sibling_table(position, i);

		if(t){
			local::Hash_elem[new_key] -> facen[i].push_back(Unit::Face());
			local::Hash_elem[new_key] -> facen[i][0].face_type = 'L';
			local::Hash_elem[new_key] -> facen[i][0].hlevel = local::Hash_elem[new_key] -> index[2];
			local::Hash_elem[new_key] -> facen[i][0].porderx = local::Hash_elem[new_key] -> n;
			local::Hash_elem[new_key] -> facen[i][0].pordery = local::Hash_elem[new_key] -> m;
			local::Hash_elem[new_key] -> facen[i][0].rank = mpi::rank;
			
			// key
			if(i == 0){	// south
				
				local::Hash_elem[new_key] -> facen[i][0].key = Get_key_fun(
										local::Hash_elem[new_key] -> index[0] - 1, 
											local::Hash_elem[new_key] -> index[1], 
											local::Hash_elem[new_key] -> index[2]);
		

			}
			else if(i == 1){	// north
				
				local::Hash_elem[new_key] -> facen[i][0].key = Get_key_fun(local::Hash_elem[new_key] -> index[0] + 1, 
											local::Hash_elem[new_key] -> index[1], 
											local::Hash_elem[new_key] -> index[2]);
		

			}
			else if(i == 2){	// west
				
				local::Hash_elem[new_key] -> facen[i][0].key = Get_key_fun(local::Hash_elem[new_key] -> index[0], 
											local::Hash_elem[new_key] -> index[1] - 1, 
											local::Hash_elem[new_key] -> index[2]);
		

			}
			else{	// east
				
				local::Hash_elem[new_key] -> facen[i][0].key = Get_key_fun(local::Hash_elem[new_key] -> index[0], 
											local::Hash_elem[new_key] -> index[1] + 1, 
											local::Hash_elem[new_key] -> index[2]);
		

			}
		}

	}

}



/// @brief
/// Generate the new unit index based on the child position and its old (parent) unit index
/// @param ith ith child
/// @param index parent's unit index
/// @param index_new child's index (size 3)
void Gen_index(int ith, int* index, int* index_new){

	assert(ith >=0 && ith<= 3 && "Error: The ith child's position is out of range.");
	
	if(ith == 0){

		index_new[0] = 2 * index[0];
		index_new[1] = 2 * index[1];
		index_new[2] = index[2] + 1;
	}
	else if(ith == 1){

		index_new[0] = 2 * index[0] + 1;
		index_new[1] = 2 * index[1];
		index_new[2] = index[2] + 1;

	}
	else if(ith == 2){

		index_new[0] = 2 * index[0] + 1;
		index_new[1] = 2 * index[1] + 1;
		index_new[2] = index[2] + 1;

	}
	else{


		index_new[0] = 2 * index[0];
		index_new[1] = 2 * index[1] + 1;
		index_new[2] = index[2] + 1;

	}

}

/// @brief 
/// Generate the new coordinates base on the child position inside the quard
/// @param ith ith child
/// @param xcoord new unit's xcoord
/// @param ycoord new unit's ycoord
void Get_coordinates(int ith, double* xcoord, double* ycoord){

	assert(ith >=0 && ith<= 3 && "Error: The ith child's position is out of range.");

	if(ith == 0){
		xcoord[0] = xcoord_new[0]; xcoord[1] = xcoord_new[1];
		ycoord[0] = ycoord_new[0]; ycoord[1] = ycoord_new[1];
	}
	else if(ith == 1){

		xcoord[0] = xcoord_new[1]; xcoord[1] = xcoord_new[2];
		ycoord[0] = ycoord_new[0]; ycoord[1] = ycoord_new[1];

	}
	else if(ith == 2){

		xcoord[0] = xcoord_new[1]; xcoord[1] = xcoord_new[2];
		ycoord[0] = ycoord_new[1]; ycoord[1] = ycoord_new[2];
	}
	else{

		xcoord[0] = xcoord_new[0]; xcoord[1] = xcoord_new[1];
		ycoord[0] = ycoord_new[1]; ycoord[1] = ycoord_new[2];
	}

}
