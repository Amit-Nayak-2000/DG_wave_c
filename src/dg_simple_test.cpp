#include "dg_simple_test.h"
#include "dg_unit.h"
#include <vector>
#include "dg_boundary_table.h"
#include "dg_local_storage.h"
#include "dg_boundary_table.h"
#include "dg_cantor_pairing.h"
#include <algorithm>	// std::sort
#include <mpi.h>
#include  <iostream>	// test
#include "dg_param.h"	//test

// forward declaration
void Construct_mpi_table(std::vector<table_elem>& north, std::vector<table_elem>& south);
void Update_mpi_boundaries(std::vector<table_elem>& north, std::vector<table_elem>& south);
void Accum_table(std::vector<table_elem>& south, std::vector<accum_elem>& south_accum);


void Simple_test(){

	// init, var = 0 at the beginning 
	

	// tarverse the hash table
	Unit* temp = local::head;

	// construct interface
	std::vector<int> n_interface(local::local_elem_num);
	std::vector<int> s_interface(local::local_elem_num);
	
	for(int k = 0; k < local::local_elem_num; ++k){

		n_interface[k] = temp -> var;

		temp = temp -> next;
	}


	// exchange info on mpi boundaries
	// first construct table
	std::vector<table_elem> north;
	std::vector<table_elem> south;
	
	Construct_mpi_table(north, south);


	// start to send info

//	if(! north.empty()){
//		for(auto& v : north ){
//			std::cout<< mpi::rank << " " << v.target_rank << "\n";
//		}
//	}
}

// only in x direction
void Construct_mpi_table(std::vector<table_elem>& north, std::vector<table_elem>& south){

	Unit* temp = local::head;

	for(int k = 0; k < local::local_elem_num; ++k){


		// south
		// interate through face 0
		for(auto& face_s : temp -> facen[0]){

			if(face_s.face_type == 'M'){	// if mpi boundary, record

				south.push_back(table_elem());
				south.back().local_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);
				south.back().target_rank = face_s.rank;
				south.back().coord = temp -> index[1];		// y coord
				south.back().hlevel = temp -> index[2]; 	// hlevel	
			}

		
		
		}
		
		// north
		// iterate through face 1
		for(auto& face_n : temp -> facen[1]){

			if(face_n.face_type == 'M'){	// if mpi boundary, record
				
				north.push_back(table_elem());
				north.back().local_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);
				north.back().target_rank = face_n.rank;
				north.back().coord = temp -> index[1];
				north.back().hlevel = temp -> index[2];	
			}

		
		
		}

		temp = temp -> next;

	}	

		// sort north and south table in the end
		std::sort(south.begin(), south.end(), compare_coord);
		std::sort(north.begin(), north.end(), compare_coord);
}



void Update_mpi_boundaries(std::vector<table_elem>& north, std::vector<table_elem>& south){

	// accumulates boundaries info
	std::vector<accum_elem> south_accum;
	std::vector<accum_elem> north_accum;

	Accum_table(south, south_accum);
	Accum_table(north, north_accum);

	int s = south_accum.size();
	int n = north_accum.size();
	
	// south send -----------------------------------------------------------------------------------------------------------------
	if(s > 0){	// there is thing to send
		MPI_Request s_request1[s], s_request2[s];	// for mpi_waitall
		MPI_Status s_status1[s], s_status2[s];		// mpi_waitall

		// first south send north receive
		int i{};
		int j{};
		for(auto& v : south_accum){	// south send
			MPI_Isend(&v.sum, 1, MPI_INT, v.rank, mpi::rank, MPI_COMM_WORLD, &s_request1[i]);	// tag = local rank
			
			std::vector<int> send_info(v.sum * 2);
			
			// serialization the struct
			for(int k = 0; k < v.sum; ++k){
				
				send_info[2 * k] = south[j].local_key;	// key
				send_info[2 * k + 1] = south[j].hlevel;	// hlevel
				++j;
	
			}
	
			MPI_Isend(&send_info, v.sum * 2, MPI_INT, v.rank, v.rank, MPI_COMM_WORLD, &s_request2[i]);
	
			++i;
			
		}
	
		MPI_Waitall(s, s_request1, s_status1);	// ensure all info recved
		MPI_Waitall(s, s_request2, s_status2);

	}
	//-----------------------------------------------------------------------------------------------------------------------
	

	// north recv ------------------------------------------------------------------------------------------------------------
	if(n > 0){
		MPI_Status status;		// dummy
		std::vector<int> recv_info;	// recv: key, hlevel
		std::vector<table_elem>::iterator it;	// declare an iterator
		it = north.begin();	// put the iterator at the begin of the north table
	
		for(auto& v : north_accum){
			int num;	// number of elem on the other side
			MPI_Recv(&num, 1, MPI_INT, v.rank, v.rank, MPI_COMM_WORLD, &status);
	
			recv_info = std::vector<int>(num * 2);
	
			MPI_Recv(&recv_info, num * 2, MPI_INT, v.rank, mpi::rank, MPI_COMM_WORLD, &status);
			
		}
	}
	//-------------------------------------------------------------------------------------------------------------------------
//
//MPI_Waitall(int count, MPI_Request array_of_requests[],
//    MPI_Status *array_of_statuses)

//MPI_Recv(void *buf, int count, MPI_Datatype datatype,
//    int source, int tag, MPI_Comm comm, MPI_Status *status)

//int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest,
//    int tag, MPI_Comm comm, MPI_Request *request)
}

/// @brief
/// Update facen in hash table.
/// @param recv_info recieved information vector.
/// @param table MPI direction table.
/// @param facei element ith face to be updates
/// @param num recieved element number.
/// @param it MPIdirection table iterator.
void Update_hash(std::vector<int>& recv_info, std::vector<table_elem>& table, 
			int facei, int num, std::vector<table_elem>& it){
	
	for(int k = 0; k < num; ++k){	// not table but number of recv elem
			
		// now update the facei neighbour
		for(auto& here : local::Hash_elem[v.local_key].facen[facei]){
			// erase old neightbours
			if(here.face_type == 'M' && here.rank == it -> target_rank)
			local::Hash_elem[it -> local_key].facen[1].erase();	


		}


	}	

}



void Accum_table(std::vector<table_elem>& south, std::vector<accum_elem>& south_accum){


	if(! south.empty()){	// if not empty
		
		south_accum.push_back(accum_elem());

		int rank1 = south.front().target_rank;
		south_accum.back().rank = rank1;

		for(auto& v : south){

			int rank2 = v.target_rank;

			if(rank2 == rank1){

				south_accum.back().sum += 1;
				
			}	
			else{
				south_accum.push_back(accum_elem());
				south_accum.back().rank = rank2;
				south_accum.back().sum += 1;
				rank1 = rank2;
		
			}		
			
		}
	}

}
