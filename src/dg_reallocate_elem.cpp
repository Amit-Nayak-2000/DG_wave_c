#include <mpi.h>
#include "dg_reallocate_elem.h"
#include "dg_local_storage.h"
#include "dg_derived_datatype.h"
#include "dg_param.h"
#include "dg_cantor_pairing.h"
#include <cassert>
#include <iostream> // test
#include <fstream>	// test
#include <cassert>	//test

// forward declaration-------------------------------------------------------------------
void Reallocate_elem(int kt);

void Send_pack(std::vector<info_pack>& send_info, std::vector<long long int>::iterator& it, int& solu_num);

void Recv_elem(int source, int tag, std::vector<info_pack>& recv_info, int& count);

void Enlarge_hash(std::vector<info_pack>& recv_info, char dir, int num_recv, std::vector<double>& solu);

void Face_pack(std::vector<face_pack>& face_info, std::vector<long long int>& send, int& num);

void Recv_face(int source, int tag, std::vector<face_pack>& recv_face);

void Fill_facen(std::vector<face_pack>& face_info);

void Erase_elem_old(std::vector<long long int>& send, char dir, int num);

void Solution_pack(std::vector<long long int>& send_list, int solu_num, std::vector<double>& solu_packed);

void Recv_solu(int source, int tag, std::vector<double>& solu);

void Write_send(int kt, std::vector<info_pack>& send_elem, int num_n, int target_rank); 	// test
void Write_recv(int kt, std::vector<info_pack>& recv_elem, int num_n, int target_rank);	//test
void Write_recv_face(int kt, std::vector<face_pack>& recv_face, int target_rank);	// test
void Write_send_face(int kt, std::vector<face_pack>& recv_face, int target_rank);
// --------------------------------------------------------------------------------------

/// @brief
/// After built the complete mapping table, now we decide how to reallocate the elements
void Reallocate_elem(int kt){

	int start = LB::elem_accum; 	// first elem global number
	int last = start + local::local_elem_num - 1;	// last elem global number

	int num_pre = LB::Send.pre.size();
	int num_next = LB::Send.next.size();

	// send to pre -----------------------------------------------------------
	std::vector<info_pack> send_elem_p(num_pre);	// element info

	std::vector<double> solu_packed_p;	// element solution

	std::vector<face_pack> face_info_p;	// face info
	// -----------------------------------------------------------------------

	// send to next ----------------------------------------------------------
	std::vector<info_pack> send_elem_n(num_next);	// element info

	std::vector<double> solu_packed_n;	// element solution

	std::vector<face_pack> face_info_n;	// face info
	//------------------------------------------------------------------------

	MPI_Request request_elem_pre, request_elem_next;
	MPI_Request request_solu_pre, request_solu_next;
	MPI_Request request_face_pre, request_face_next;


	if(num_pre > 0){	// something to send

		auto it = LB::Send.pre.begin();

		// pack info to send
		int solu_num{};
		Send_pack(send_elem_p, it, solu_num);	// element info

		Solution_pack(LB::Send.pre, solu_num, solu_packed_p);	// solutions

		int num_n{};
		Face_pack(face_info_p, LB::Send.pre, num_n);	// face info
		
		// test--------------------------------------------------------------
//		Write_send(kt, send_elem, num_pre, mpi::rank - 1); 	// test
//		Write_send_face(kt, face_info, mpi::rank - 1);
		//----------------------------------------------------------------------

		// ready to send 

		MPI_Isend(&send_elem_p[0], num_pre, Hash::Elem_type, mpi::rank - 1, 
				mpi::rank, MPI_COMM_WORLD, &request_elem_pre);	// tag = rank

		MPI_Isend(&solu_packed_p[0], solu_num, MPI_DOUBLE, mpi::rank - 1, 
				mpi::rank + 3, MPI_COMM_WORLD, &request_solu_pre);

		MPI_Isend(&face_info_p[0], num_n, Hash::Face_type, mpi::rank - 1, 
				mpi::rank + 1, MPI_COMM_WORLD, &request_face_pre);

		Erase_elem_old(LB::Send.pre, 'p', num_pre);
	}
	if(num_next > 0){	


		auto it = LB::Send.next.begin();

		int solu_num{};
		Send_pack(send_elem_n, it, solu_num);

		Solution_pack(LB::Send.next, solu_num, solu_packed_n);	// solutions
		int num_n{};
		Face_pack(face_info_n, LB::Send.next, num_n);

		// test		
//		Write_send(kt, send_elem, num_next, mpi::rank + 1); 	// test
//		Write_send_face(kt, face_info, mpi::rank + 1);

		// non-blocking ----------------------------------------------------------------------------------
		MPI_Isend(&send_elem_n[0], num_next, Hash::Elem_type, mpi::rank + 1, 
				mpi::rank, MPI_COMM_WORLD, &request_elem_next);

		MPI_Isend(&solu_packed_n[0], solu_num, MPI_DOUBLE, mpi::rank + 1, 
				mpi::rank + 4, MPI_COMM_WORLD, &request_solu_next);

		MPI_Isend(&face_info_n[0], num_n, Hash::Face_type, mpi::rank + 1, 
				mpi::rank + 1, MPI_COMM_WORLD, &request_face_next);
		// -----------------------------------------------------------------------------------------------

		Erase_elem_old(LB::Send.next, 'n', num_next);
	}
	
	// recv
	if(mpi::rank == 0){	// first proc
		
		if(LB::proc_mapping_table[1].gnum - 1 > last){	// recv from next

			std::vector<info_pack> recv_info;
			std::vector<face_pack> recv_face;
			std::vector<double> solu;

			int recv_num{};	
			Recv_elem(mpi::rank + 1, mpi::rank + 1, recv_info, recv_num);

			Recv_solu(mpi::rank + 1, mpi::rank + 4, solu);

			Recv_face(mpi::rank + 1, mpi::rank + 2, recv_face);

		//	Write_recv_face(kt, recv_face, mpi::rank + 1);	// test

			Enlarge_hash(recv_info, 'n', recv_num, solu);
			Fill_facen(recv_face);
		}

	}
	else if(mpi::rank == mpi::num_proc - 1){

		if(LB::proc_mapping_table[mpi::rank].gnum < start){	// recv from pre

			std::vector<info_pack> recv_info;
			std::vector<face_pack> recv_face;
			std::vector<double> solu;

			int recv_num{};	
			Recv_elem(mpi::rank - 1, mpi::rank - 1, recv_info, recv_num);
			
			Recv_solu(mpi::rank - 1, mpi::rank + 3, solu);

			Recv_face(mpi::rank - 1, mpi::rank, recv_face);

//			Write_recv_face(kt, recv_face, mpi::rank - 1);	// test

			Enlarge_hash(recv_info, 'p', recv_num, solu);
			Fill_facen(recv_face);

		}

	}
	else{	// proc in between

		if(LB::proc_mapping_table[mpi::rank].gnum < start){	// recv from pre
			std::vector<info_pack> recv_info;
			std::vector<face_pack> recv_face;
			std::vector<double> solu;
			int recv_num{};	
		
			Recv_elem(mpi::rank - 1, mpi::rank - 1, recv_info, recv_num);

			Recv_solu(mpi::rank - 1, mpi::rank + 3, solu);


//			Write_recv(kt, recv_info, recv_num, mpi::rank - 1);	//test

			Recv_face(mpi::rank - 1, mpi::rank, recv_face);

//			Write_recv_face(kt, recv_face, mpi::rank - 1);	// test

			Enlarge_hash(recv_info, 'p', recv_num, solu);
			Fill_facen(recv_face);
		}
		
		if(LB::proc_mapping_table[mpi::rank + 1].gnum - 1 > last){	// recv from next

			std::vector<info_pack> recv_info;
			std::vector<face_pack> recv_face;
			std::vector<double> solu;

			int recv_num{};	
			Recv_elem(mpi::rank + 1, mpi::rank + 1, recv_info, recv_num);
			
			Recv_solu(mpi::rank + 1, mpi::rank + 4, solu);

			Recv_face(mpi::rank + 1, mpi::rank + 2, recv_face);

//			Write_recv_face(kt, recv_face, mpi::rank + 1);	// test

			Enlarge_hash(recv_info, 'n', recv_num, solu);
			Fill_facen(recv_face);
		}
	}
	
	// wait
	if(num_pre > 0){

		MPI_Status s_elem_pre, s_solu_pre, s_face_pre;

		MPI_Wait(&request_elem_pre, &s_elem_pre);
		MPI_Wait(&request_solu_pre, &s_solu_pre);
		MPI_Wait(&request_face_pre, &s_face_pre);
	}

	if(num_next > 0){

		MPI_Status s_elem_next, s_solu_next, s_face_next;

		MPI_Wait(&request_elem_next, &s_elem_next);
		MPI_Wait(&request_solu_next, &s_solu_next);
		MPI_Wait(&request_face_next, &s_face_next);
	}


}

void Write_send_face(int kt, std::vector<face_pack>& recv_face, int target_rank){


	// generate the file name
	std::string my_rank = std::to_string(mpi::rank);
	std::string filename = "../send_face_info/rank" + my_rank + ".dat";
	std::ofstream myfile; 	// stream class to write on files	


	myfile.open(filename, std::ios::out | std::ios::app);

	if(!myfile){	// if this file does not exist

		myfile.open(filename, std::ios::out | std::ios::trunc);
	}

	// record
	myfile << "============================================== \n";
	myfile << "time " << kt << "\n";
	int num_face = recv_face.size();
//	myfile << "target_rank " << target_rank << " num_face " << num_face << "\n";
	myfile << "num_face " << num_face << "\n";

	for(auto& v : recv_face){

//		myfile << "owners_key" << v.owners_key << " facei "<<v.facei << " face_type "<< v.face_type << 
//		" neighbour "<< v.key << " n_rank "<< v.rank<< "\n";
		myfile << v.owners_key << " "<<v.facei << " "<< v.face_type << 
		" "<< v.key << " "<< v.rank<< "\n";
	}
	

	myfile.close();

}
void Write_recv_face(int kt, std::vector<face_pack>& recv_face, int target_rank){


	// generate the file name
	std::string my_rank = std::to_string(mpi::rank);
	std::string filename = "../recv_face_info/rank" + my_rank + ".dat";
	std::ofstream myfile; 	// stream class to write on files	


	myfile.open(filename, std::ios::out | std::ios::app);

	if(!myfile){	// if this file does not exist

		myfile.open(filename, std::ios::out | std::ios::trunc);
//		myfile << "my rank " << mpi::rank << "\n";
	}

	// record
	myfile << "============================================== \n";
	myfile << "time " << kt << "\n";

	int num_face = recv_face.size();
	myfile << "num_face " << num_face << "\n";
//	myfile << "target_rank " << target_rank << " num_face " << num_face << "\n";

	for(auto& v : recv_face){

//		myfile << "owners_key" << v.owners_key << " facei "<<v.facei << " face_type "<< v.face_type << 
//		" neighbour "<< v.key << " n_rank "<< v.rank<< "\n";

		myfile << v.owners_key << " "<<v.facei << " "<< v.face_type << 
		" "<< v.key << " "<< v.rank << " ref_x " << v.ref_x[0] << " " << v.ref_x[1] 
		<< " ref_y " << v.ref_y[0] << " " << v.ref_y[1] 
		<< "\n";
	}
	

	myfile.close();

}

void Write_send(int kt, std::vector<info_pack>& send_elem, int num_n, int target_rank){

	// generate the file name
	std::string my_rank = std::to_string(mpi::rank);
	std::string filename = "../send_info/rank" + my_rank + "send.dat";
	std::ofstream myfile; 	// stream class to write on files	


	myfile.open(filename, std::ios::out | std::ios::app);

	if(!myfile){	// if this file does not exist

		myfile.open(filename, std::ios::out | std::ios::trunc);
		myfile << "my rank " << mpi::rank << "\n";

	}

	// record
	myfile << "============================================== \n";
	myfile << "time " << kt << "\n";
	myfile << "target_rank " << target_rank << " number elem "<< num_n << "\n";

//	int key_first = Get_key_fun(LB::my_rank_first -> index[0], LB::my_rank_first -> index[1], LB::my_rank_first -> index[2]);
	if((local::local_elem_num - num_n) > 0){
		myfile << "new first " << LB::my_rank_first -> index[0] << " " << 
					LB::my_rank_first -> index[1] << " " <<
					LB::my_rank_first -> index[2] << 
					" new_last "<< LB::my_rank_last -> index[0] 
					<< " "<<LB::my_rank_last -> index[1] << " "<< LB::my_rank_last -> index[2]<< "\n";
	}
//	int num{};
	for(auto& v : send_elem){

		long long int key = Get_key_fun(v.index[0], v.index[1], v.index[2]);

		myfile << key << " ";

		
	}
	myfile << "\n";
	myfile << "============================================== \n";

	myfile.close();

}


void Write_recv(int kt, std::vector<info_pack>& recv_elem, int num_n, int target_rank){

	// generate the file name
	std::string my_rank = std::to_string(mpi::rank);
	std::string filename = "../recv_info/rank" + my_rank + ".dat";
	std::ofstream myfile; 	// stream class to write on files	


	myfile.open(filename, std::ios::out | std::ios::app);

	if(!myfile){	// if this file does not exist

		myfile.open(filename, std::ios::out | std::ios::trunc);
		myfile << "my rank " << mpi::rank << "\n";

	}

	// record
	myfile << "============================================== \n";
	myfile << "time " << kt << "\n";
	myfile << "target_rank " << target_rank << " number elem "<< num_n << "\n";

	if(local::local_elem_num > 0){
		myfile << "new first " << LB::my_rank_first -> index[0] << " " << 
					LB::my_rank_first -> index[1] << " " <<
					LB::my_rank_first -> index[2] << 
					" new_last "<< LB::my_rank_last -> index[0] 
					<< " "<<LB::my_rank_last -> index[1] << " "<< LB::my_rank_last -> index[2]<< "\n";
	}
//	int num{};
	for(auto& v : recv_elem){

		long long int key = Get_key_fun(v.index[0], v.index[1], v.index[2]);

		myfile << key << " status " << v.status<< " c position "<< v.child_position<< " x p "<< v.xcoords[0]
			<< " " << v.xcoords[1] << " y "<< v.ycoords[0]<< v.ycoords[1] 
			<< " ref_x " << v.ref_x[0] << " " << v.ref_x[1] << " ref_y " << v.ref_y[0] 
			<< " " << v.ref_y[1] 
			<< "\n";

		
	}
	myfile << "\n";
	myfile << "============================================== \n";

	myfile.close();

}

/// @brief
/// After receive element information we create units in Hash table and store them.
/// @parma recv_info The received element info (without facen).
/// @param dir The info coming direction (pre: p, next: 'n').
/// @pram num_recv received element number. 
/// @param solu solution vector. 
void Enlarge_hash(std::vector<info_pack>& recv_info, char dir, int num_recv, std::vector<double>& solu){

	assert((dir == 'p' || dir == 'n') && "The sendign direction can only be 'p' or 'n'.");

	Unit* temp_head = nullptr;	// head pointer to the recved linked list 
	long long int pre_key;
	int nodei{};

	for(auto it = recv_info.begin(); it != recv_info.end(); ++it){

		long long int key = Get_key_fun((*it).index[0], (*it).index[1], (*it).index[2]);

		local::Hash_elem[key] = new Unit();

		if(it == recv_info.begin()){
			temp_head = local::Hash_elem[key];	// temp_head points to the first element. 
		}

		local::Hash_elem[key] -> n = (*it).n;		// for now n == m
		local::Hash_elem[key] -> m = (*it).n;		
		
		for(int i = 0; i < 3; ++i){
			local::Hash_elem[key] -> index[i] = (*it).index[i];
		}

		local::Hash_elem[key] -> status = (*it).status;

		local::Hash_elem[key] -> child_position = (*it).child_position; 

		for(int i = 0; i < 2; ++i){	
			local::Hash_elem[key] -> xcoords[i] = (*it).xcoords[i]; 
			local::Hash_elem[key] -> ycoords[i] = (*it).ycoords[i]; 
			local::Hash_elem[key] -> ref_x[i] = (*it).ref_x[i]; 
			local::Hash_elem[key] -> ref_y[i] = (*it).ref_y[i]; 
		}

		// solution	
		int num_solu = ((*it).n + 1) * ((*it).n + 1);
		int num_bound = 4*((*it).n + 1);
		for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

			local::Hash_elem[key] -> solution[equ] = std::vector<double> (num_solu);	// allocate
	
			for(int i = 0; i < num_solu; ++i){
				
				local::Hash_elem[key] -> solution[equ][i] = solu[nodei];
				++nodei;
			}
		}

		//mapped geometry class
		local::Hash_elem[key] -> holdmetrics.updateOrder((*it).n); //allocate

		for(int k = 0; k < 12; k++){
			switch(k) {
				case 0:
					//x_node
					for(int i = 0; i < num_solu; ++i){
						local::Hash_elem[key] -> holdmetrics.x_node[i] = solu[nodei]; //originally map[nodej]
						++nodei; //originally nodej
					}
					break;
				case 1:
					//y_node
					for(int i = 0; i < num_solu; ++i){
						local::Hash_elem[key] -> holdmetrics.y_node[i] = solu[nodei];
						++nodei;
					}
					break;
				case 2:
					//x_boundary
					for(int i = 0; i < num_bound; ++i){
						local::Hash_elem[key] -> holdmetrics.x_boundary[i] = solu[nodei];
						++nodei;
					}
					break;
				case 3:
					//y_boundary
					for(int i = 0; i < num_bound; ++i){
						local::Hash_elem[key] -> holdmetrics.y_boundary[i] = solu[nodei];
						++nodei;
					}
					break;
				case 4:
					//delx_delxi
					for(int i = 0; i < num_solu; ++i){
						local::Hash_elem[key] -> holdmetrics.delx_delxi[i] = solu[nodei];
						++nodei;
					}
					break;
				case 5:
					//delx_deleta
					for(int i = 0; i < num_solu; ++i){
						local::Hash_elem[key] -> holdmetrics.delx_deleta[i] = solu[nodei];
						++nodei;
					}
					break;
				case 6:
					//dely_delxi
					for(int i = 0; i < num_solu; ++i){
						local::Hash_elem[key] -> holdmetrics.dely_delxi[i] = solu[nodei];
						++nodei;
					}
					break;
				case 7:
					//dely_deleta
					for(int i = 0; i < num_solu; ++i){
						local::Hash_elem[key] -> holdmetrics.dely_deleta[i] = solu[nodei];
						++nodei;
					}
					break;
				case 8:
					//jacobian
					for(int i = 0; i < num_solu; ++i){
						local::Hash_elem[key] -> holdmetrics.jacobian[i] = solu[nodei];
						++nodei;
					}
					break;
				case 9:
					//normalx
					for(int i = 0; i < num_bound; ++i){
						local::Hash_elem[key] -> holdmetrics.boundary_normal[i] = solu[nodei];
						++nodei;
					}
					break;
				case 10:
					//normaly
					for(int i = 0; i < num_bound; ++i){
						local::Hash_elem[key] -> holdmetrics.boundary_normal_y[i] = solu[nodei];
						++nodei;
					}
					break;
				case 11:
					//scaling factor
					for(int i = 0; i < num_bound; ++i){
						local::Hash_elem[key] -> holdmetrics.scaling_factor[i] = solu[nodei];
						++nodei;
					}
					break;
				default:
					break;
			}
		}

	
		if(it != recv_info.begin()){
			local::Hash_elem[pre_key] -> next = local::Hash_elem[key];

		}

		pre_key = key;
	}

	if(local::local_elem_num == 0){	// if no local_elem

		local::head = temp_head;

	}
	else{

		if(dir == 'p'){		// put the elements at the beginning of the linked list
			
			local::Hash_elem[pre_key] -> next = LB::my_rank_first;
	
			local::head = temp_head;
	
		}	
		else{	// put the element at the end of the linked list
			LB::my_rank_last -> next = temp_head;
		}

	}

	// update local element number
	local::local_elem_num += num_recv;

}

/// @brief
/// Put the neighbours information inside the hash table. 
/// @param face_info Received face information. 
void Fill_facen(std::vector<face_pack>& face_info){

	for(auto& v : face_info){

		local::Hash_elem[v.owners_key] -> facen[v.facei].push_back(Unit::Face());

		local::Hash_elem[v.owners_key] -> facen[v.facei].back().face_type = v.face_type;

		local::Hash_elem[v.owners_key] -> facen[v.facei].back().hlevel = v.hlevel;

		local::Hash_elem[v.owners_key] -> facen[v.facei].back().porderx = v.porderx;

		local::Hash_elem[v.owners_key] -> facen[v.facei].back().pordery = v.pordery;

		local::Hash_elem[v.owners_key] -> facen[v.facei].back().key = v.key;

		local::Hash_elem[v.owners_key] -> facen[v.facei].back().rank = v.rank;

		for(int i = 0; i < 2; ++i){

			local::Hash_elem[v.owners_key] -> facen[v.facei].back().ref_x[i] = v.ref_x[i];
			local::Hash_elem[v.owners_key] -> facen[v.facei].back().ref_y[i] = v.ref_y[i];

		}


	}

}


/// @brief
/// Receive element solution from sender. 
/// @param source Sender's rank.
/// @param tag Tag of the message. 
/// @param solu Buffer for the message. 
void Recv_solu(int source, int tag, std::vector<double>& solu){

	int count{};

	MPI_Status status1, status2;

	MPI_Probe(source, tag, MPI_COMM_WORLD, &status1);

	MPI_Get_count(&status1, MPI_DOUBLE, &count);

	solu = std::vector<double>(count);

	MPI_Recv(&solu[0], count, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status2);
}

/// @brief
/// Receive element information from sender. 
/// @param source Sender's rank.
/// @param tag Tag of the message. 
/// @param recv_info Buffer for the message. 
/// @param count number of element received. 
void Recv_elem(int source, int tag, std::vector<info_pack>& recv_info, int& count){

	MPI_Status status1, status2;

	MPI_Probe(source, tag, MPI_COMM_WORLD, &status1);

	MPI_Get_count(&status1, Hash::Elem_type, &count);

	recv_info = std::vector<info_pack>(count);

	MPI_Recv(&recv_info[0], count, Hash::Elem_type, source, tag, MPI_COMM_WORLD, &status2);

}

/// @brief
/// Receive face info.
/// @param source source rank.
/// @param tag message tag.
/// @param recv_face vector to store the message. 
void Recv_face(int source, int tag, std::vector<face_pack>& recv_face){

	MPI_Status status1, status2;
	MPI_Probe(source, tag, MPI_COMM_WORLD, &status1);

	int count{};
	MPI_Get_count(&status1, Hash::Face_type, &count);
	
	recv_face = std::vector<face_pack>(count);

	MPI_Recv(&recv_face[0], count, Hash::Face_type, source, tag, MPI_COMM_WORLD, &status2);
}

/// @brief
/// Pack up the element solutions. No need to pre-allocate solu_packed. 
/// @param send_list the list of element to be send.
/// @param solu_num solution number, obtained by the Send_pack() function. 
/// @param solu_packed Pack up the solutions of all the sending elements together.  
void Solution_pack(std::vector<long long int>& send_list, int solu_num, std::vector<double>& solu_packed){

	// allocate solution package. 
	solu_packed = std::vector<double>(solu_num);

	int i{};

	//iterate thru list of elements
	for(auto& key : send_list){
		//pack solutions
		for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

			for(auto& v : local::Hash_elem[key] -> solution[equ]){

				solu_packed[i] = v;
				++i;
			}	
		}

		for(int k = 0; k < 12; k++){
			switch(k) {
				case 0:
					//x_node
					for(auto& v : local::Hash_elem[key] -> holdmetrics.x_node){
							solu_packed[i] = v;
							++i;
					}
					break;
				case 1:
					//y_node
					for(auto& v : local::Hash_elem[key] -> holdmetrics.y_node){
							solu_packed[i] = v;
							++i;
					}
					break;
				case 2:
					//x_boundary
					for(auto& v : local::Hash_elem[key] -> holdmetrics.x_boundary){
							solu_packed[i] = v;
							++i;
					}
					break;
				case 3:
					//y_boundary
					for(auto& v : local::Hash_elem[key] -> holdmetrics.y_boundary){
							solu_packed[i] = v;
							++i;
					}
					break;
				case 4:
					//delx_delxi
					for(auto& v : local::Hash_elem[key] -> holdmetrics.delx_delxi){
							solu_packed[i] = v;
							++i;
					}
					break;
				case 5:
					//delx_deleta
					for(auto& v : local::Hash_elem[key] -> holdmetrics.delx_deleta){
							solu_packed[i] = v;
							++i;
					}
					break;
				case 6:
					//dely_delxi
					for(auto& v : local::Hash_elem[key] -> holdmetrics.dely_delxi){
							solu_packed[i] = v;
							++i;
					}
					break;
				case 7:
					//dely_deleta
					for(auto& v : local::Hash_elem[key] -> holdmetrics.dely_deleta){
							solu_packed[i] = v;
							++i;
					}
					break;
				case 8:
					//jacobian
					for(auto& v : local::Hash_elem[key] -> holdmetrics.jacobian){
							solu_packed[i] = v;
							++i;
					}
					break;
				case 9:
					//normalx
					for(auto& v : local::Hash_elem[key] -> holdmetrics.boundary_normal){
							solu_packed[i] = v;
							++i;
					}
					break;
				case 10:
					//normaly
					for(auto& v : local::Hash_elem[key] -> holdmetrics.boundary_normal_y){
							solu_packed[i] = v;
							++i;
					}
					break;
				case 11:
					//scaling factor
					for(auto& v : local::Hash_elem[key] -> holdmetrics.scaling_factor){
							solu_packed[i] = v;
							++i;
					}
					break;
				default:
					break;
			}
		}

	}

}

/// @brief
/// Pack sending information. 
/// @param send_info The sending vector.
/// @param it iterator at the beginning of Sending list (pre or next). 
/// @param solu_num solution number. For next step to transfer the solutions. 
void Send_pack(std::vector<info_pack>& send_info, std::vector<long long int>::iterator& it, int& solu_num){

	int boundary_num{};

	for(auto& v : send_info){
		
		v.n = local::Hash_elem[*it] -> n;

		solu_num += (v.n + 1) * (v.n + 1); 	// assume same poly order in x and y direction. 

		boundary_num += 4*(v.n + 1);

		for(int i = 0; i < 3; ++i){
			v.index[i] = local::Hash_elem[*it] -> index[i];
		}

		v.status = local::Hash_elem[*it] -> status;

		v.child_position = local::Hash_elem[*it] -> child_position;

		for(int i = 0; i < 2; ++i){
		
			v.xcoords[i] = local::Hash_elem[*it] -> xcoords[i];
			v.ycoords[i] = local::Hash_elem[*it] -> ycoords[i];
			v.ref_x[i] = local::Hash_elem[*it] -> ref_x[i];
			v.ref_y[i] = local::Hash_elem[*it] -> ref_y[i];
		}

		++it;
	}

	solu_num *= dg_fun::num_of_equation;

	//after solu_num has been calculated for solutions, tack on the size of mapped geometry data
	solu_num += (7*(solu_num / dg_fun::num_of_equation)) + (5*boundary_num);

}

/// @brief
/// Pack face info. 
/// @param face_info Store all the face neighbours in a vector.
/// @param send Sending list. 
/// @param num Number of neighbour. 
void Face_pack(std::vector<face_pack>& face_info, std::vector<long long int>& send, int& num){

	for(auto& v : send){	// traverse the sending list

		for(int i = 0; i < 4; ++i){	// four faces

			for(auto it = local::Hash_elem[v] -> facen[i].begin(); 
				it != local::Hash_elem[v] -> facen[i].end(); ++it){

				face_info.push_back(face_pack());
				++num;

				face_info.back().owners_key = v; 
				face_info.back().facei = i;
				face_info.back().face_type = it -> face_type;
				face_info.back().hlevel = it -> hlevel;
				face_info.back().porderx = it -> porderx;
				face_info.back().pordery = it -> pordery;
				face_info.back().key = it -> key;
				face_info.back().rank = it -> rank;
				face_info.back().Copy_ref(it -> ref_x, it -> ref_y);

			}

		}

	}


}


/// @brief
/// Erase the elements in the sending list from the Hash table
/// @param send Sending list.
/// @param dir Sending direction (pre: 'p', next: 'n').
/// @parma num number of element to be send. 
void Erase_elem_old(std::vector<long long int>& send, char dir, int num){

	assert((dir == 'p' || dir == 'n') && "The sendign direction can only be 'p' or 'n'.");

	local::local_elem_num -= num;	// remaining elements.

	// take care of linked list
	if(local::local_elem_num == 0){	// nothing left

		local::head = nullptr;

	}
	else{

		if(dir == 'p'){

			long long int last = send.back();
			local::head = local::Hash_elem[last] -> next;
		}
		else{	// 'n'
			LB::my_rank_last -> next = nullptr;
		}		

	}
	
	// erase elements
	for(auto& v : send){

		local::Hash_elem.erase(v);

	}

}
