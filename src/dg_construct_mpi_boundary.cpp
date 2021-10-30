#include "dg_construct_mpi_boundary.h"
#include "dg_local_storage.h"
#include "dg_unit.h"
#include "dg_cantor_pairing.h"
#include "dg_nodal_2d_storage.h"
#include <unordered_map>
#include "dg_param.h"
#include "dg_search_rank.h"
#include "dg_put_into_mpi_table.h"
#include <iostream>	//test

/// @brief
/// Construct MPI boundaries and physical boundaries. 
/// Each element has 4 faces
/// \verbatim
///   ---f1----
///   |       |
/// f2|       |f3
///   |       |
///   ---f0----
/// \endverbatim
/// Each face is a vector of struct "mpi_table" that contains all the info required. 
void MPI_boundary_construct(){
	
	// start by the first element
	Unit* temp = local::head;
	
	// traverse the linked-list
	for(int k = 0; k < local::local_elem_num; ++k){

		// on the south physical boundary?
		if(temp -> index[0] == 0){

			temp -> facen[0].push_back(Unit::Face());

			temp -> facen[0][0].face_type = 'B'; // Yes
		}
		else{	// No, search south neighbour
			int ni = temp -> index[0] - 1;
			int nj = temp -> index[1];
			long long int nkey = Get_key_fun(ni, nj, 0); 	// before adapt

			temp -> facen[0].push_back(Unit::Face());

			std::unordered_map<long long int, Unit*>::const_iterator got = local::Hash_elem.find(nkey);
			// not found, so on the MPI boundary
			if(got == local::Hash_elem.end()){
			
				temp -> facen[0][0].face_type = 'M';
				temp -> facen[0][0].hlevel = 0;
				temp -> facen[0][0].porderx = grid::nmin;	// for uniform mesh, we can record
				temp -> facen[0][0].pordery = grid::nmin;	// for uniform mesh, we can record
				temp -> facen[0][0].key = nkey;		

				int target_rank = Target_rank(ni, nj);

				temp -> facen[0][0].rank = target_rank;

				// put this element into the mpi table
				auto it_face = temp -> facen[0].begin();
				Put_in_mpi_table(temp, it_face, hrefinement::south);
			}
			else{	// if found, record info

				temp -> facen[0][0].face_type = 'L';
				temp -> facen[0][0].hlevel = 0;
				temp -> facen[0][0].porderx = grid::nmin;	
				temp -> facen[0][0].pordery = grid::nmin;	
				temp -> facen[0][0].key = nkey;	
				temp -> facen[0][0].rank = mpi::rank;	
			

			}
		}

		// on the north physical boundary?
		temp -> facen[1].push_back(Unit::Face());
		if(temp -> index[0] == (SortMesh::num_of_element_x - 1)){
			
			temp -> facen[1][0].face_type = 'B';
		}
		else{	// no
			int ni = temp -> index[0] + 1;
			int nj = temp -> index[1];
			long long int nkey = Get_key_fun(ni, nj, 0); 	// before adapt
			
			std::unordered_map<long long int, Unit*>::const_iterator got = local::Hash_elem.find(nkey);
			// not found, so on the MPI boundary
			if(got == local::Hash_elem.end()){

				temp -> facen[1][0].face_type = 'M';
				temp -> facen[1][0].hlevel = 0;
				temp -> facen[1][0].porderx = grid::nmin;	
				temp -> facen[1][0].pordery = grid::nmin;	
				temp -> facen[1][0].key = nkey;	
				int target_rank = Target_rank(ni, nj);

				temp -> facen[1][0].rank = target_rank;
				
				// record in MPI boundary table
				auto it_face = temp -> facen[1].begin();
				Put_in_mpi_table(temp, it_face, hrefinement::north);
			}
			else{

				temp -> facen[1][0].face_type = 'L';
				temp -> facen[1][0].hlevel = 0;
				temp -> facen[1][0].porderx = grid::nmin;	
				temp -> facen[1][0].pordery = grid::nmin;	
				temp -> facen[1][0].key = nkey;	
				temp -> facen[1][0].rank = mpi::rank;	

			}
			
		}

		// on the west physical boundary?
		temp -> facen[2].push_back(Unit::Face());
		if(temp -> index[1] == 0){

			temp -> facen[2][0].face_type = 'B';

		}
		else{	// no
			int ni = temp -> index[0];
			int nj = temp -> index[1] - 1;
			long long int nkey = Get_key_fun(ni, nj, 0); 	// before adapt
			
			std::unordered_map<long long int, Unit*>::const_iterator got = local::Hash_elem.find(nkey);
			// not found, so on the MPI boundary
			if(got == local::Hash_elem.end()){
				
				temp -> facen[2][0].face_type = 'M';
				temp -> facen[2][0].hlevel = 0;
				temp -> facen[2][0].porderx = grid::nmin;	
				temp -> facen[2][0].pordery = grid::nmin;	
				temp -> facen[2][0].key = nkey;	
				int target_rank = Target_rank(ni, nj);

				temp -> facen[2][0].rank = target_rank;
				auto it_face = temp -> facen[2].begin();
				Put_in_mpi_table(temp, it_face, hrefinement::west);
			}
			else{	// locally
				
				temp -> facen[2][0].face_type = 'L';
				temp -> facen[2][0].hlevel = 0;
				temp -> facen[2][0].porderx = grid::nmin;	
				temp -> facen[2][0].pordery = grid::nmin;	
				temp -> facen[2][0].key = nkey;	
				temp -> facen[2][0].rank = mpi::rank;	

			}
			
		}

		// on the east physical boundary?
		temp -> facen[3].push_back(Unit::Face());
		if(temp -> index[1] == SortMesh::num_of_element_y - 1){

			temp -> facen[3][0].face_type = 'B';
		}
		else{	// no
			int ni = temp -> index[0];
			int nj = temp -> index[1] + 1;
			long long int nkey = Get_key_fun(ni, nj, 0); 	// before adapt
			
			std::unordered_map<long long int, Unit*>::const_iterator got = local::Hash_elem.find(nkey);
			// not found, so on the MPI boundary
			if(got == local::Hash_elem.end()){
				
				temp -> facen[3][0].face_type = 'M';
				temp -> facen[3][0].hlevel = 0;
				temp -> facen[3][0].porderx = grid::nmin;	
				temp -> facen[3][0].pordery = grid::nmin;	
				temp -> facen[3][0].key = nkey;	

				int target_rank = Target_rank(ni, nj);

				temp -> facen[3][0].rank = target_rank;
				auto it_face = temp -> facen[3].begin();
				Put_in_mpi_table(temp, it_face, hrefinement::east);
			}
			else{

				temp -> facen[3][0].face_type = 'L';
				temp -> facen[3][0].hlevel = 0;
				temp -> facen[3][0].porderx = grid::nmin;	
				temp -> facen[3][0].pordery = grid::nmin;	
				temp -> facen[3][0].key = nkey;	
				temp -> facen[3][0].rank = mpi::rank;	
				
			
			}
			
		}
		
		// pointer move to next element
		temp = temp -> next;
	}

}

