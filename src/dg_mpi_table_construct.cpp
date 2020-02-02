#include "dg_mpi_table_construct.h"
#include "dg_local_storage.h"
#include "dg_unit.h"
#include "dg_boundary_table.h"
#include "dg_cantor_pairing.h"
#include <algorithm>
#include <vector>
#include "dg_elem_length.h"
#include <mpi.h>
#include "dg_param.h"
#include <iostream>	// test

// forward declaration ------------------------------------------------------------------
void Accum_table(std::vector<table_elem>& south, std::vector<accum_elem>& south_accum);

void Erase_old_face(std::vector<Unit::Face>::iterator& it_face, std::vector<table_elem>::iterator& it, int facei);

void Sender_recver(int s, int n, std::vector<accum_elem>& south_accum, std::vector<accum_elem>& north_accum, 
			std::vector<table_elem>& south, std::vector<table_elem>& north, int update_dir);

void Update_hash(std::vector<int>& recv_info, std::vector<table_elem>& table, 
			int facei, int num1, std::vector<table_elem>::iterator& it);

void Clear_tables();
//---------------------------------------------------------------------------------------

/// @brief 
/// Construct MPI boundary tables.
/// @param north North direction table.
/// @param south South direction table.
// only in x direction
void Construct_mpi_table(std::vector<table_elem>& north, std::vector<table_elem>& south){

	Unit* temp = local::head;

	for(int k = 0; k < local::local_elem_num; ++k){

		int pre_rank = -1;

		// south
		// interate through face 0
		for(auto& face_s : temp -> facen[0]){


			if(face_s.face_type == 'M' && face_s.rank != pre_rank){	// if mpi boundary and rank changes, record

				south.push_back(table_elem());
				south.back().local_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);
				south.back().target_rank = face_s.rank;
				south.back().coord = temp -> ycoords[1];		// y coord
				south.back().hlevel = temp -> index[2]; 	// hlevel	
			}

		
			pre_rank = face_s.rank;
		
		}
	
		pre_rank = -1;
	
		// north
		// iterate through face 1
		for(auto& face_n : temp -> facen[1]){


			if(face_n.face_type == 'M' && face_n.rank != pre_rank){	// if mpi boundary, record
		
				north.push_back(table_elem());
				north.back().local_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);
				north.back().target_rank = face_n.rank;
				north.back().coord = temp -> ycoords[1];
				north.back().hlevel = temp -> index[2];	

			}

		
			pre_rank = face_n.rank;
		
		}

		temp = temp -> next;

	}	

		// sort north and south table in the end
		std::sort(south.begin(), south.end(), compare_coord);
		std::sort(north.begin(), north.end(), compare_coord);


}

/// @brief
/// Updates MPI boundaries
/// @param north north MPI table
/// @param south south MPI table
void Update_mpi_boundaries(std::vector<table_elem>& north, std::vector<table_elem>& south){

	Accum_table(south, hrefinement::south_accum);
	Accum_table(north, hrefinement::north_accum);

	int s = hrefinement::south_accum.size();
	int n = hrefinement::north_accum.size();
//if(mpi::rank == 0){
//	
//	std::cout<< "----------------------------- \n";
//	for(auto& v : north){
//		std::cout<< "local key: " << v.local_key << " target rank: " << v.target_rank << "\n";
//
//	}
//	std::cout<< "----------------------------- \n";
////
////	std::cout<< "----------------------------- \n";
////	for(auto& v : hrefinement::north_accum){
////
////		std::cout<< "rank: " << v.rank << "\n";
////		std::cout<< "sum: " << v.sum << "\n";
////	}
////
////	std::cout<< "----------------------------- \n";
//	//std::cout << s << " " << n << "\n";
//}
	// south send, north recv
	Sender_recver(s, n, hrefinement::south_accum, hrefinement::north_accum, south, north, 1);

	// north send, south recv. 
	Sender_recver(n, s, hrefinement::north_accum, hrefinement::south_accum, north, south, 0);

}


/// @brief
/// Send and recv info to update MPI boundaries. 
/// @param s number of element to send.
/// @param n number of element to recv.
/// @param south_accum sender's accum table.
/// @param north_accum recver's accum table.
/// @param south sender's MPI table.
/// @param north recver's MPI table.
/// @param update_dir update_direaction. Should be the recver's direction.
void Sender_recver(int s, int n, std::vector<accum_elem>& south_accum, std::vector<accum_elem>& north_accum, 
			std::vector<table_elem>& south, std::vector<table_elem>& north, int update_dir){


	if(s > 0){	// there is thing to send
		MPI_Request s_request[s];	// for mpi_waitall
		MPI_Status  s_status[s];		// mpi_waitall

		int i{};
		int j{};
		for(auto& v : south_accum){	// south send
		//	MPI_Isend(&v.sum, 1, MPI_INT, v.rank, mpi::rank, MPI_COMM_WORLD, &s_request1[i]); // tag = local rank
			

//if(mpi::rank == 3){
//
//	std::cout<< "rank 3 send : " << v.sum << "\n";
//
//
//}

			std::vector<int> send_info(v.sum * 4);	// key, hlevel, porderx, pordery
			
			// serialization the struct
			for(int k = 0; k < v.sum; ++k){
				
				send_info[4 * k] = south[j].local_key;	// key
				send_info[4 * k + 1] = south[j].hlevel;	// hlevel
				send_info[4 * k + 2] = local::Hash_elem[south[j].local_key] -> n;	// porderx
				send_info[4 * k + 3] = local::Hash_elem[south[j].local_key] -> m;	// pordery
				++j;
			}
	
			MPI_Isend(&send_info[0], v.sum * 4, MPI_INT, v.rank, mpi::rank, MPI_COMM_WORLD, &s_request[i]); 
//if(mpi::rank == 2){
//
//	
//
//}	
			++i;
			
		}
	
//		MPI_Waitall(s, s_request1, s_status1);	// ensure all info recved
		MPI_Waitall(s, s_request, s_status);

	}
	//-----------------------------------------------------------------------------------------------------------------------
	

	// north recv ------------------------------------------------------------------------------------------------------------
	if(n > 0){
	
		for(auto& v : north_accum){

			MPI_Status status1, status2;		// dummy

			int num;	// number of elem on the other side

			MPI_Probe(v.rank, v.rank, MPI_COMM_WORLD, &status1);

			MPI_Get_count(&status1, MPI_INT, &num);


		//	MPI_Recv(&num, 1, MPI_INT, v.rank, v.rank, MPI_COMM_WORLD, &status1);


//if(mpi::rank == 0){
//	
//	std::cout << "rank " << mpi::rank << " " << num << "\n";
//}
			std::vector<int> recv_info(num);	
			std::vector<table_elem>::iterator it;	// declare an iterator
			it = north.begin();	// put the iterator at the begin of the north table
	
			MPI_Recv(&recv_info[0], num, MPI_INT, v.rank, v.rank, MPI_COMM_WORLD, &status2);
//if(mpi::rank == 0){
//	for(int i = 0; i < num / 4; ++i){
//		std::cout<< "key: "<< recv_info[i * 4] << "hlevel: " << recv_info[i * 4 + 1] << "\n";
//
//	}	
//
//}			
			Update_hash(recv_info, north, update_dir, num, it);	// north recv
		}
	}
	//-----------------------------------------------------------------------------------------------------------------------

}

/// @brief
/// Update facen in hash table.
/// @param recv_info recieved information vector.
/// @param table MPI direction table.
/// @param facei element ith face to be updates
/// @param num recieved element number * 4.
/// @param it MPI direction table iterator.
void Update_hash(std::vector<int>& recv_info, std::vector<table_elem>& table, 
			int facei, int num1, std::vector<table_elem>::iterator& it){
	
	int l_tot{};

	int num = num1 / 4;
//if(mpi::rank == 0){
//
//
//	std::cout<< "num = " << num << "\n";
//
//}
	for(int k = 0; k < num; ){	// not table but number of recv elem


		int l_local = Elem_length(it -> hlevel);	// local elem length
		int l_n = Elem_length(recv_info[4 * k + 1]);		// neighbour elem length


		if(l_local == l_n){	// if same size

			// delete old face info			
			auto it_face = local::Hash_elem[it -> local_key] -> facen[facei].begin();
			Erase_old_face(it_face, it, facei);

			// type, hlevel, porder, key, rank	
			Unit::Face obj = {'M', recv_info[4 * k + 1], recv_info[4 * k + 2], recv_info[4 * k + 3], 
						recv_info[4 * k], it -> target_rank};	


			local::Hash_elem[it -> local_key] -> facen[facei].emplace(it_face, obj); 

			++k;
			++it;

		}
		else if(l_local < l_n){	// local is smaller
			l_tot = 0;
			// recv_info stall (k), it loop
			// traverse mpi table until face matches
			while(l_tot < l_n && it != table.end()){
				
				// erase old face info
				auto it_face = local::Hash_elem[it -> local_key] -> facen[facei].begin();
				Erase_old_face(it_face, it, facei);

				Unit::Face obj = {'M', recv_info[4 * k + 1], recv_info[4 * k + 2], 
							recv_info[4 * k + 3], recv_info[4 * k], it -> target_rank};

				it_face = local::Hash_elem[it -> local_key] -> facen[facei].emplace(it_face, obj);

				l_tot += l_local;
				
				++it;	// go to next local elem
			}

			++k;
			
		}
		else{	// local is larger
			l_tot = 0;
			// erase old face info
			auto it_face = local::Hash_elem[it -> local_key] -> facen[facei].begin();
			Erase_old_face(it_face, it, facei);
			
			// it stall, recv_info loop
			while(k < num && l_tot < l_local){
				Unit::Face obj = {'M', recv_info[4 * k + 1], recv_info[4 * k + 2], 
						recv_info[4 * k + 3], recv_info[4 * k], it -> target_rank};
			
				it_face = local::Hash_elem[it -> local_key] -> facen[facei].emplace(it_face, obj);
				l_tot += l_n;

				++it_face;
				++k;
			}

			++it;
		}
		
	}	

}

/// @brief
/// Erase the old mpi boundary info on facei
/// @param it_face iterator on facen.
/// @param it iterator of mpi table.
/// @param facei face direction. 
void Erase_old_face(std::vector<Unit::Face>::iterator& it_face, std::vector<table_elem>::iterator& it, int facei){

	int target_rank = it -> target_rank;
	for(; it_face != local::Hash_elem[it -> local_key] -> facen[facei].end(); ){
		// erase the old info	
		if(it_face -> face_type == 'M' && (it_face -> rank == target_rank)){
			
			it_face = local::Hash_elem[it -> local_key] -> facen[facei].erase(it_face); // it_face move to the next
	
			if(it_face -> rank != target_rank){break;}
	
		}
		else{
			++it_face;
		}
	}

}

/// @brief
/// Form rank Accumulate table.
/// @param south MPI south table.
/// @param north MPI north table. 
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

/// @brief
/// Remove all the elements inside the tables (MPI tables and accumlation table)
void Clear_tables(){

	hrefinement::south_accum.clear();
	hrefinement::north_accum.clear();
	hrefinement::south.clear();
	hrefinement::north.clear();


}
