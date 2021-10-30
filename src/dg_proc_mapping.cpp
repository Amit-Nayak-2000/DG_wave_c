#include "dg_proc_mapping.h"
#include "dg_unit.h"
#include "dg_local_storage.h"
#include "dg_param.h"
#include <cmath>
#include <mpi.h>
#include "dg_cantor_pairing.h"
#include <algorithm>
#include "dg_status_table.h"
#include "dg_elem_length.h"
#include "dg_load_struct.h"
#include <unordered_map>
#include "dg_boundary_table.h"
#include "dg_derived_datatype.h"
#include "dg_element_load.h"
#include <cassert>
#include <algorithm>
#include <iostream>	//test
//#include "dg_write_mpi_table.h"	//test
//#include "dg_write_send_list.h"	//test

// forward declaration -----------------------------------------
//double Elem_load(int porder);

void Build_mapping_table();
void Build_mapping_table_quality();

void Update_neighbours();

void Neighbour_change(int facei, long long int n_key, long long int my_key, int rank);

void Ownership_one_dir(std::unordered_map<int, std::vector<mpi_table>>& mtable);

void Send_recv_ownership(std::unordered_map<int, std::vector<mpi_table>>& sendo, 
			std::unordered_map<int, std::vector<mpi_table>>& recvo, int facei);

void Update_mpib(std::vector<owner_struct>& recv_info, std::unordered_map<int, std::vector<mpi_table>>& otable, 
		int facei, int num, int target_rank);

void Update_mpi_boundary();

void Change_face(int num, std::vector<owner_struct>& recv_info, std::vector<mpi_table>::iterator& ito, 
			std::vector<Unit::Face>::iterator& it_face);
// ------------------------------------------------------------

void Build_mapping_table_quality(){
	
	Unit* temp = local::head;

	std::vector<double> lprefix_load(local::local_elem_num);
	int pmapping;

	// local prefix sum of load
	lprefix_load[0] = Elem_load(temp -> n);
	temp = temp -> next;
	for(int k = 1; k < local::local_elem_num; ++k){
		lprefix_load[k] = Elem_load(temp -> n) + lprefix_load[k - 1];
		
		temp = temp -> next;

	}	
	
	double local_load_sum = lprefix_load.back();	// local computational load sum
	double exscan_sum{};	// the load of former processor


	// Global prefix sum of load
	MPI_Exscan(&local_load_sum, &exscan_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	// calculate for the average load
	if(mpi::rank == (mpi::num_proc - 1)){ // last proc does the job

		double load_tol = exscan_sum + local_load_sum;

		LB::load_average = load_tol / mpi::num_proc;

	}
	
	// broadcast average load
	MPI_Bcast(&LB::load_average, 1, MPI_DOUBLE, mpi::num_proc - 1, MPI_COMM_WORLD);
	
	// Global element number
	MPI_Exscan(&local::local_elem_num, &LB::elem_accum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	// form processor mapping table
	temp = local::head;
	LB::my_rank_first = local::head;
	int proc_pre = - 1;

	std::vector<double> table_prefix_sum; 	// global prefix sum of the elements in the table

	for(int k = 0; k < local::local_elem_num; ++k){
		
		lprefix_load[k] += exscan_sum;

		pmapping = std::floor((lprefix_load[k] - 0.01) / LB::load_average);

		assert(pmapping >= 0 && "processor mapping is smaller than 0.");	// check

		// form partial mapping table
		if(pmapping != proc_pre){

			LB::proc_mapping_table.push_back({pmapping, k + LB::elem_accum});
		
			table_prefix_sum.push_back(lprefix_load[k]);	// global prefix sum 
	
			proc_pre = pmapping;
		}	
		
		// form sending_envelope
		if(pmapping < mpi::rank){	// need to be moved to the former proc
			long long int key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);

			LB::Send.pre.push_back(key);

			LB::my_rank_first = temp -> next;
		}
		else if(pmapping > mpi::rank){

			long long int key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);

			LB::Send.next.push_back(key);
		}
		else{

			LB::my_rank_last = temp;
		}


		temp = temp -> next;
	}

	// send the last element's mapping number to the next rank
	// rank0 ~ rank_max-1 send
	MPI_Request request;
	if(mpi::rank != (mpi::num_proc - 1)){
	
		int last_rank = LB::proc_mapping_table.back().irank;
		MPI_Send(&last_rank, 1, MPI_INT, mpi::rank + 1, mpi::rank + 1, MPI_COMM_WORLD);// tag == recver's rank

	}
	
	// rank1 ~ rank_max recv
	if(mpi::rank != 0){

		int pre_rank;
		MPI_Status status1;

		MPI_Recv(&pre_rank, 1, MPI_INT, mpi::rank - 1, mpi::rank, MPI_COMM_WORLD, &status1);

		int first_rank = LB::proc_mapping_table.front().irank;

		if(first_rank == pre_rank){	// if equal than erase the first column

			LB::proc_mapping_table.erase(LB::proc_mapping_table.begin());
			table_prefix_sum.erase(table_prefix_sum.begin());	 
		}

	}
	
	// call allgather to gather the size of the mapping table on each proc
	int sizet = LB::proc_mapping_table.size();
	std::vector<int> sizea(mpi::num_proc);	// vector to store the sizet
	MPI_Allgather(&sizet, 1, MPI_INT, &sizea[0], 1, MPI_INT, MPI_COMM_WORLD);

	// prepare to build the whole mapping table
	std::vector<int> recvcounts(mpi::num_proc);
	std::vector<int> recvcounts_prefix(mpi::num_proc);
	std::vector<int> displs(mpi::num_proc);
	std::vector<int> displs_prefix(mpi::num_proc);
	for(int i = 0; i < mpi::num_proc; ++i){
		recvcounts[i] = sizea[i] * 2;
		recvcounts_prefix[i] = sizea[i];
		if(i > 0){
	
			displs[i] = displs[i - 1] + recvcounts[i - 1];
			displs_prefix[i] = displs_prefix[i - 1] + recvcounts_prefix[i - 1];
		}

	}
	std::vector<int> senda(sizet * 2);	// send buffer
	for(int i = 0; i < sizet; ++i){
		senda[2 * i] = LB::proc_mapping_table[i].irank;
		senda[2 * i + 1] = LB::proc_mapping_table[i].gnum;

	}

	// form the complete mapping table
	std::vector<int> recv_buff(mpi::num_proc * 2);
	MPI_Allgatherv(&senda[0], sizet * 2, MPI_INT, &recv_buff[0], &recvcounts[0], &displs[0], MPI_INT, MPI_COMM_WORLD);

	// ========================================================================================================
	// recv the global prefix sum of the first element in the new partition
	std::vector<double> recv_prefix_sum(mpi::num_proc);
	MPI_Allgatherv(&table_prefix_sum[0], sizet, MPI_DOUBLE, &recv_prefix_sum[0], &recvcounts_prefix[0],
	    		&displs_prefix[0], MPI_DOUBLE, MPI_COMM_WORLD);

	// get the optimal bottle neck
	assert(mpi::num_proc > 1 && "The processor number should be at least 2. \n");
	for(int i = 1; i < mpi::num_proc; ++i){

		double neck = recv_prefix_sum[i] - recv_prefix_sum[i - 1];

		LB::opt_bottleneck = std::max(neck, LB::opt_bottleneck);
		   
	}
	// ========================================================================================================


	// rebuild the mapping table
	LB::proc_mapping_table.clear();	// clear all the elements
	for(int i = 0; i < mpi::num_proc; ++i){
		
		LB::proc_mapping_table.push_back({recv_buff[2 * i], recv_buff[2 * i + 1]});

	}

	// Updates local facen based on the Sending list
	Update_neighbours();
	
	// get the opt bottleneck now turn off the switch to use this function.
	LB::first = false;
}


/// @brief Calculate the sum of the local computational load.
void Build_mapping_table(){
	
	Unit* temp = local::head;

	std::vector<double> lprefix_load(local::local_elem_num);
	int pmapping;

	// local prefix sum of load
	lprefix_load[0] = Elem_load(temp -> n);
	temp = temp -> next;
	for(int k = 1; k < local::local_elem_num; ++k){
		lprefix_load[k] = Elem_load(temp -> n) + lprefix_load[k - 1];
		
		temp = temp -> next;

	}	
	
	double local_load_sum = lprefix_load.back();	// local computational load sum
	double exscan_sum{};	// the load of former processor
	
	// Global prefix sum of load
	MPI_Exscan(&local_load_sum, &exscan_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	// Global element number
	MPI_Exscan(&local::local_elem_num, &LB::elem_accum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	// form processor mapping table
	temp = local::head;
	LB::my_rank_first = local::head;
	int proc_pre = - 1;
	for(int k = 0; k < local::local_elem_num; ++k){
		
		lprefix_load[k] += exscan_sum;

		pmapping = std::floor((lprefix_load[k] - 0.01) / LB::load_average);

		assert(pmapping >= 0 && "processor mapping is smaller than 0.");	// check

		// form partial mapping table
		if(pmapping != proc_pre){

			LB::proc_mapping_table.push_back({pmapping, k + LB::elem_accum});
			
			proc_pre = pmapping;
		}	
		
		// form sending_envelope
		if(pmapping < mpi::rank){	// need to be moved to the former proc
			long long int key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);

			LB::Send.pre.push_back(key);

			LB::my_rank_first = temp -> next;
		}
		else if(pmapping > mpi::rank){

			long long int key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);

			LB::Send.next.push_back(key);
		}
		else{

			LB::my_rank_last = temp;
		}


		temp = temp -> next;
	}

	// send the last element's mapping number to the next rank
	// rank0 ~ rank_max-1 send
	MPI_Request request;
	if(mpi::rank != (mpi::num_proc - 1)){
	
		int last_rank = LB::proc_mapping_table.back().irank;
		MPI_Send(&last_rank, 1, MPI_INT, mpi::rank + 1, mpi::rank + 1, MPI_COMM_WORLD);// tag == recver's rank

	}
	
	// rank1 ~ rank_max recv
	if(mpi::rank != 0){

		int pre_rank;
		MPI_Status status1;

		MPI_Recv(&pre_rank, 1, MPI_INT, mpi::rank - 1, mpi::rank, MPI_COMM_WORLD, &status1);

		int first_rank = LB::proc_mapping_table.front().irank;

		if(first_rank == pre_rank){	// if equal than erase the first column

			LB::proc_mapping_table.erase(LB::proc_mapping_table.begin());
		}

	}
	
	// call allgather to gather the size of the mapping table on each proc
	int sizet = LB::proc_mapping_table.size();
	std::vector<int> sizea(mpi::num_proc);	// vector to store the sizet
	MPI_Allgather(&sizet, 1, MPI_INT, &sizea[0], 1, MPI_INT, MPI_COMM_WORLD);

	// prepare to build the whole mapping table
	std::vector<int> recvcounts(mpi::num_proc);
	std::vector<int> displs(mpi::num_proc);
	for(int i = 0; i < mpi::num_proc; ++i){
		recvcounts[i] = sizea[i] * 2;
		if(i > 0){
	
			displs[i] = displs[i - 1] + recvcounts[i - 1];
		}

	}
	std::vector<int> senda(sizet * 2);	// send buffer
	for(int i = 0; i < sizet; ++i){
		senda[2 * i] = LB::proc_mapping_table[i].irank;
		senda[2 * i + 1] = LB::proc_mapping_table[i].gnum;

	}

	// form the complete mapping table
	std::vector<int> recv_buff(mpi::num_proc * 2);
	MPI_Allgatherv(&senda[0], sizet * 2, MPI_INT, &recv_buff[0], &recvcounts[0], &displs[0], MPI_INT, MPI_COMM_WORLD);

	// rebuild the mapping table
	LB::proc_mapping_table.clear();	// clear all the elements
	for(int i = 0; i < mpi::num_proc; ++i){
		
		LB::proc_mapping_table.push_back({recv_buff[2 * i], recv_buff[2 * i + 1]});

	}

	// Updates local facen based on the Sending list
	Update_neighbours();

	//----------------------
//	Write_send_list_all();
	//-----------------------
}

/// @brief
/// Updates the element neighbours based on the Send struct.
void Update_neighbours(){
	
	// first updates the Send pre list ----------------------------------------------------------------------------------
	for(auto& key : LB::Send.pre){

		for(int facei = 0; facei < 4; ++facei){

			auto it = local::Hash_elem[key] -> facen[facei].begin();

			// traverse the face
			for(; it != local::Hash_elem[key] -> facen[facei].end(); ++it){

				// skip 'M' and 'B'
				if(it -> face_type == 'L'){

					long long int n_key = it -> key;	// neighbour's key
					
					// find n_key in pre list
					if(std::find(LB::Send.pre.begin(), LB::Send.pre.end(), n_key) == LB::Send.pre.end()){	// if not find

						// find n_key in the next list
						if(std::find(LB::Send.next.begin(), LB::Send.next.end(), n_key) == LB::Send.next.end()){	// if not find
				
							// This element stays locally, updates
							it -> face_type = 'M';
				//			it -> rank = mpi::rank; //current rank	
							// updates neighbour	
							Neighbour_change(facei, n_key, key, mpi::rank - 1);
							
						}
						else{	// if find
							
							// This element will be sent to next proc
							it -> face_type = 'M';
							it -> rank = mpi::rank + 1;
						
							// updates neighbour
							Neighbour_change(facei, n_key, key, mpi::rank - 1);
						}

					}
					else{	// if find in pre list
						it -> rank = mpi::rank - 1;	// update the rank

					}
						
				}

			}		
		
		}
	} //-------------------------------------------------------------------------------------------------------------


	// updates the next list-----------------------------------------------------------------------------------------
	for(auto& key : LB::Send.next){

		for(int facei = 0; facei < 4; ++facei){
			
			auto it = local::Hash_elem[key] -> facen[facei].begin();

			// traverse the face
			for(; it != local::Hash_elem[key] -> facen[facei].end(); ++it){

				// skip 'M' and 'B'
				if(it -> face_type == 'L'){

					long long int n_key = it -> key;	// neighbour's key

					// find n_key in the next list
					if(std::find(LB::Send.next.begin(), LB::Send.next.end(), n_key) == LB::Send.next.end()){	// if not find
			
						// This element stays locally, updates
						it -> face_type = 'M';
						it -> rank = mpi::rank; //current rank	
						// updates neighbour	
						Neighbour_change(facei, n_key, key, mpi::rank + 1);
						
					}
					else{
						it -> rank = mpi::rank + 1;
					}

				}
			}		

		}

	}
	// --------------------------------------------------------------------------------------------------------------

}

/// @brief
/// Change neighbour's corresponding face. 
/// @param facei current face direction. 
/// @param n_key neighbour's key. 
/// @param my_key current element's key. 
/// @param rank Rank number that the current element will reside on. 
void Neighbour_change(int facei, long long int n_key, long long int my_key, int rank){

	// neighbour's face direction.
	int oface = Opposite_dir(facei);

	for(auto it = local::Hash_elem[n_key] -> facen[oface].begin(); it != local::Hash_elem[n_key] -> facen[oface].end(); ++it){

		if(it -> key == my_key){

			it -> face_type = 'M';
			it -> rank = rank;

			break;
		}

	}

}


/// @brief
/// Fill in ownership in MPI table in one direction.
/// @param mtable MPI boundary table of the corresponding direction. 
void Ownership_one_dir(std::unordered_map<int, std::vector<mpi_table>>& mtable){


	// keys are inherited from MPI boundary tables
	for(auto& v1 : mtable){

		for(auto& v2 : v1.second){
			if(std::find(LB::Send.pre.begin(), LB::Send.pre.end(), v2.local_key) != LB::Send.pre.end()){ // if find in pre list
				v2.owners_rank = mpi::rank - 1;
			}
			else if(std::find(LB::Send.next.begin(), LB::Send.next.end(), v2.local_key) != LB::Send.next.end()){	// if find in next list
	
				v2.owners_rank = mpi::rank + 1;
	
			}
			else{ // not inside the sending list, record directly
				v2.owners_rank = mpi::rank;
			}
		}

	}

}

/// @brief
/// Updates the MPI boundaries before repartitioning. 
void Update_mpi_boundary(){

	// form the element future ownership----------------------------------------------
	Ownership_one_dir(hrefinement::north);
	Ownership_one_dir(hrefinement::south);

//	Write_table_all(hrefinement::south, hrefinement::north);
	
	Ownership_one_dir(hrefinement::west);
	Ownership_one_dir(hrefinement::east);
	//-------------------------------------------------------------------------

	// x direction-----------------------------------------------------------------------------------
	
	// north send and south recv
	Send_recv_ownership(hrefinement::north, hrefinement::south, 0);
	// south send and north recv
	Send_recv_ownership(hrefinement::south, hrefinement::north, 1);
	//-----------------------------------------------------------------------------------------------

	// y direction-----------------------------------------------------------------------------------
	
	// west send and east recv
	Send_recv_ownership(hrefinement::west, hrefinement::east, 3);

	// east send and west recv
	Send_recv_ownership(hrefinement::east, hrefinement::west, 2);
	//-----------------------------------------------------------------------------------------------
}

/// @brief
/// Send and recv MPI boundary table to updates the MPI boundaries. 
/// @param sendo Sender's MPI boundary table.
/// @param recvo Recver's MPI boundary table.
/// @param send_accum Sender's accumulation table. 
/// @param recv_accum Receiver's accumulation table. 
/// @param facei the face direction to be updated. 
void Send_recv_ownership(std::unordered_map<int, std::vector<mpi_table>>& sendo, 
			std::unordered_map<int, std::vector<mpi_table>>& recvo, int facei){
	
	int size_s = sendo.size();
	int size_r = recvo.size();

	std::unordered_map<int, std::vector<owner_struct>> send_info;

	MPI_Request send_request[size_s];
	int isend{};
	
	// sender
	if(size_s > 0){	// there is something to send

		for(auto& v : sendo){

			int target_rank = v.first;

			int num_elem = v.second.size();	

			send_info[target_rank] = std::vector<owner_struct>(num_elem);

			auto it = v.second.begin();
	
			// serialization the struct
			for(int k = 0; k < num_elem; ++k){
				send_info[target_rank][k].local_key = it -> local_key;	
				send_info[target_rank][k].owners_rank = it -> owners_rank;
				++it;
			}

			MPI_Isend(&send_info[target_rank][0], num_elem, Hash::Owner_type, target_rank,
				    mpi::rank, MPI_COMM_WORLD, &send_request[isend]);

			++isend;

		}

	}


	// recver
	if(size_r > 0){

		for(auto& v : recvo){
			
			MPI_Status status1, status2;

			int num{};

			MPI_Probe(v.first, v.first, MPI_COMM_WORLD, &status1);

			MPI_Get_count(&status1, Hash::Owner_type, &num);

			std::vector<owner_struct> recv_info(num);
			MPI_Recv(&recv_info[0], num, Hash::Owner_type, v.first, v.first, MPI_COMM_WORLD, &status2);

			Update_mpib(recv_info, recvo, facei, num, v.first);
		}

	}

	// wait for isend
	if(size_s > 0){
		MPI_Status status1[size_s];

		MPI_Waitall(size_s, send_request, status1);

	}


}

/// @brief
/// Update MPI boundaries based on the ownership table.
/// @param recv_info Received infomation.
/// @param otable MPI boundary table.
/// @param ito iterator of the MPI boundary table.
/// @param num1 number of element received * 2.
void Update_mpib(std::vector<owner_struct>& recv_info, std::unordered_map<int, std::vector<mpi_table>>& otable, 
		int facei, int num, int target_rank){


	for(auto ito = otable[target_rank].begin(); ito != otable[target_rank].end(); ++ito){

		for(auto it_face = local::Hash_elem[ito -> local_key] -> facen[facei].begin();
			it_face != local::Hash_elem[ito -> local_key] -> facen[facei].end(); ++it_face){

			if(it_face -> face_type == 'M' && it_face -> rank == target_rank){

				Change_face(num, recv_info, ito, it_face);
			}

		}

	}

}

/// @brief
/// Compare the ownership of neighbour element and the stored rank to decide whether to update the facen.
/// @param k k-th element inside the received info.
/// @param recv_info Received information.
/// @param ito interator of the updating MPI boundary table.
/// @param it_face iterator of facen at the corresponding position. 
void Change_face(int num, std::vector<owner_struct>& recv_info, std::vector<mpi_table>::iterator& ito, 
			std::vector<Unit::Face>::iterator& it_face){

	for(int k = 0; k < num; ++k){	// loop to find the neighbour

		if(recv_info[k].local_key == it_face -> key){	// find it


			if(recv_info[k].owners_rank == (ito -> owners_rank)){	// if will be in the same rank
				
				// change 'M' to 'L'
				it_face -> face_type = 'L';
				it_face -> rank = ito -> owners_rank;
			}
			else{	// not in the same rank
				int rank_old = it_face -> rank;
		
				if(rank_old != recv_info[k].owners_rank){	// if this element will to assign to a new rank
					// change to target rank
					it_face -> rank = recv_info[k].owners_rank;
		
				}	// else no change
			}

		}
	}
}



