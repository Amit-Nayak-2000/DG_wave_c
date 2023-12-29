#include "dg_io.h"
#include <mpi.h>
#include "dg_param.h"
#include "dg_local_storage.h"
#include <iomanip>	// std::setfill
#include <string>
#include <iostream>
#include <sstream>	// std::stringstream
#include <fstream>	// read and write to file
#include "dg_unit.h"
#include "dg_cantor_pairing.h"
#include <vector>
#include <unordered_map>
#include "dg_interface_construct.h"
#include "dg_nodal_2d_storage.h"
#include "dg_basis.h"
#include "dg_single_index.h"
#include "dg_affine_map.h"
#include "dg_interpolate_to_new_points.h"

// forward declaration-----------------------------------------------------------------------------------
void Write_mesh(double t, int pre_elem);

void Interpolate_to_four_corner(Unit* temp, std::unordered_map<int, std::vector<double>>& four);
//-------------------------------------------------------------------------------------------------------

/// @brief
/// Output data in serial order.
/// @param t current time.
void Serial_io(double t){

	int pre_elem{};

	MPI_Exscan(&local::local_elem_num, &pre_elem, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	for(int k = 0; k < mpi::num_proc; ++k){

		if(mpi::rank == k){
			Write_mesh(t, pre_elem);

		}
		else{
			// std::cout << "We are here" << std::endl;
			MPI_Barrier(MPI_COMM_WORLD);		
		
		}
		

	}
	

}

// global file number
int file_num = 1;

/// @brief
/// Processor 1 create the file. All the processor write the data in order.
/// In the end processor 1 close the file.
/// @param t current time step
void Write_mesh(double t, int pre_elem){

	// generate the file name
	std::stringstream ss;
	ss << fileinfo::output_place <<"output" << std::setfill('0') << std::setw(5) << file_num << ".dat";
	// ss << fileinfo::output_place <<"output" << std::setfill('0') << std::setw(5) << ".csv." << file_num;

	std::string filename = 	ss.str();
	std::ofstream myfile; 	// stream class to write on files	

	// traverse the linked-list
	Unit* temp = local::head;
	int elem = pre_elem;

	// processor open file
	if(mpi::rank == 0){
		
		myfile.open(filename, std::ios::trunc);	// truncate the old file
		// myfile << "X,Y,P,U,V,L,PO" << "\n";
		myfile << "TITLE = \"DG Plot\"" << "\n";
		myfile << "VARIABLES = \"X\", \"Y\", \"P\", \"U\", \"V\", \"POLYNOMIAL ORDER\"" << "\n";
		
	}
	else{

		myfile.open(filename, std::ios::app);	// All output operations are performed at the end of the file
	}

	// write solutions
	for(int iel = 0; iel < local::local_elem_num; ++iel){

		myfile << std::fixed;
		myfile << std::setprecision(5);

		double x, y;

		// double del_x = temp->xcoords[1] - temp->xcoords[0];
		// double del_y = temp->ycoords[1] - temp->ycoords[0];

		int linnumpoints = 4*(temp->n+1);

		myfile << "ZONE T=\"" << iel + 1 << "\", I=" << linnumpoints << ", J=" << linnumpoints << ", DATAPACKING=POINT" << "\n";
		// myfile << "ZONE T=\"" << iel + 1 << "\", I=" << temp->n + 1 << ", J=" << temp->n + 1 << ", DATAPACKING=POINT" << "\n";

		++elem;

		//works
		std::vector<double> linx = linarray(-1, 1, linnumpoints);

		std::vector<double> T;	// only need one interpolation matrix (porderx = pordery)
		Polynomial_interpolate_matrix(nodal::gl_points[temp -> n], linx, T);

		std::unordered_map<int, std::vector<double>> fnew;

		std::unordered_map<int, std::vector<double>> newcoords;

		std::unordered_map<int, std::vector<double>> oldcoords;

		oldcoords[0] = std::vector<double>(temp->holdmetrics.x_node.size());
		oldcoords[1] = std::vector<double>(temp->holdmetrics.y_node.size());
		oldcoords[2] = std::vector<double>(temp->holdmetrics.y_node.size());

		for(int j = 0; j <= temp->n; j++){
			for(int i = 0; i <= temp->n; i++){
				int index = Get_single_index(i,  j, temp->n + 1);
				oldcoords[0][index] = temp->holdmetrics.x_node[index];
				oldcoords[1][index] = temp->holdmetrics.x_node[index];
				oldcoords[2][index] = temp->holdmetrics.y_node[index];


			}

			
		}

		

		//seems to work
		CoursetoFineInterp(temp->n , linnumpoints-1 , nodal::gl_points[temp->n], nodal::gl_points[temp->n], temp->solution, linx, linx, fnew);
		//testing dis
		CoursetoFineInterp(temp->n , linnumpoints-1, nodal::gl_points[temp->n], nodal::gl_points[temp->n], oldcoords, linx, linx, newcoords);

		// for(int i = 0; i < temp->holdmetrics.x_node.size(); i++){

		// }

		// std::cout << temp->holdmetrics.x_node.size()  << " | " << newxnode.size() << std::endl;
		// for(int j = 0; j <= temp->n; ++j){
		for(int j = 0; j <= linnumpoints - 1; ++j){
			// double y = Affine_mapping(nodal::gl_points[temp->n][j], temp -> ycoords[0], del_y);
			// for(int i = 0; i <= temp->n; ++i){
			for(int i = 0; i <= linnumpoints - 1; ++i){
				// double x = Affine_mapping(nodal::gl_points[temp->n][i], temp -> xcoords[0], del_x);
				// int index = Get_single_index(i,  j, temp->n + 1);
				int index = Get_single_index(i,  j, linnumpoints);
				// int MyIndex = Get_single_index(i, j, temp->n + 1);
				// int revindex = Get_single_index(temp->n - i,  temp->n - j, temp->n + 1);

				// x = temp->holdmetrics.x_node[index];
				// y = temp->holdmetrics.y_node[index];

				// x = newxnode[index];
				// y = newynode[index];

				// x = oldcoords[0][index];
				// y = oldcoords[1][index];

				// //output X, Y, P, U, V
				// myfile << x << " " << y << " " << temp -> solution[0][index] << " " << temp -> solution[1][index] << " " << temp -> solution[2][index] << " " << temp->n << "\n";

				//output X, Y, P, U, V
				// myfile << newcoords[0][index] << "," << newcoords[2][index] << "," << fnew[0][index] << "," <<  fnew[1][index] << "," << fnew[2][index] << "," << iel + 1 << "," << temp->n << "\n";
				myfile << newcoords[0][index] << " " << newcoords[2][index] << " " << fnew[0][index] << " " <<  fnew[1][index] << " " << fnew[2][index] << " " << temp->n <<  "\n";

			}
		}
		
		temp = temp -> next;
	}

	


	// close the file
	myfile.close();
	++file_num;

}

/// @brief
/// Interpolates the interior solutions to four corner points
/// @param temp pointer to the current element.
/// @param four solutions on the four corner points. <equ, four solutions>
// four solution sequnece
//     2 --------- 3
//     |           |
//     |           |
//     |           |
//     |           |
//     0 --------- 1
void Interpolate_to_four_corner(Unit* temp, std::unordered_map<int, std::vector<double>>& four){

	// construct interface in y direction.
	Construct_interface_y(temp);

	int a{};
	
	// then interpolates to the four corners
	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

		// copy the vector
		std::vector<double> s_left(temp -> n + 1);
		std::vector<double> s_right(temp -> n + 1);

		for(int i = 0; i <= (temp -> n); ++i){

			s_left[i] = temp -> solution_int_l[a];
			s_right[i] = temp -> solution_int_r[a];
			++a;	

		}

		// allocate space for corner hash
		four[equ] = std::vector<double> (4);

		// lower left
		four[equ][0] = Interpolate_to_boundary(temp -> n, s_left, nodal::lagrange_l[temp -> n]);		

		// lower right
		four[equ][1] = Interpolate_to_boundary(temp -> n, s_right, nodal::lagrange_l[temp -> n]);		

		// upper left
		four[equ][2] = Interpolate_to_boundary(temp -> n, s_left, nodal::lagrange_r[temp -> n]);		

		// upper right
		four[equ][3] = Interpolate_to_boundary(temp -> n, s_right, nodal::lagrange_r[temp -> n]);		

	}


	// clear solution_int_l/r
	(temp -> solution_int_l).clear();
	(temp -> solution_int_r).clear();
}


