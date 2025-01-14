#include "dg_interface_construct.h"
#include "dg_unit.h"
#include <vector>
#include <unordered_map>
#include "dg_param.h"
#include "dg_basis.h"
#include "dg_nodal_2d_storage.h"
#include "dg_single_index.h"
#include <iostream>

/// @brief
/// Use Lagrange interpolants to obtain the solution on the element left and right boundaries. X direction.
/// @param temp pointer points to the current element. 
void Construct_interface_x(Unit* temp){

	int m = temp -> m;
	int n = temp -> n;

	int size = dg_fun::num_of_equation * (m + 1);
	
	temp -> solution_int_l = std::vector<double>(size);
	temp -> solution_int_r = std::vector<double>(size);

	int now{};

	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

		for(int j = 0; j <= m; ++j ){

			std::vector<double> s_array(n + 1);

			for(int i = 0; i <= n; ++i){

				int nodei = Get_single_index(i, j, m + 1);

				s_array[i] = temp -> solution[equ][nodei];
			}

			temp -> solution_int_l[now] = Interpolate_to_boundary(n, s_array, nodal::lagrange_l[n]);
			temp -> solution_int_r[now] = Interpolate_to_boundary(n, s_array, nodal::lagrange_r[n]);


			++now;	
		}

	}

//if(mpi::rank == 1){
//
//		std::cout << temp -> index[0] << temp -> index[1]<< "\n";
//		std::cout.precision(17);
//
//		for(auto& solu : temp -> solution_int_l){
//
//			std::cout<< std::fixed<< solu << "\n";
////			std::cout<<solu << " ";
//
//		}
//		std::cout<< "\n";
//
//}

}

/// @brief
/// Use Lagrange interpolants to obtain the solution on the element left and right boundaries. Y direction.
/// @param temp pointer points to the current element. 
void Construct_interface_y(Unit* temp){

	int m = temp -> m;
	int n = temp -> n;

	int size = dg_fun::num_of_equation * (n + 1);
	temp -> solution_int_l = std::vector<double>(size);
	temp -> solution_int_r = std::vector<double>(size);

	int num{};
	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

		for(int i = 0; i <= n; ++i ){

			std::vector<double> s_array(m + 1);

			for(int j = 0; j <= m; ++j){

				int nodei = Get_single_index(i, j, m + 1);

				s_array[j] = temp -> solution[equ][nodei];
			}

			temp -> solution_int_l[num] = Interpolate_to_boundary(m, s_array, nodal::lagrange_l[m]);
			temp -> solution_int_r[num] = Interpolate_to_boundary(m, s_array, nodal::lagrange_r[m]);
			++num;	
		}

//if(mpi::rank == 0){
////std::cout << "equ " << equ << "\n";
//	if(equ == 2){
////		std::cout<< "i = "<< i << " nodei " << num << " " << temp -> solution_int_l[num] << " "
////				<< temp -> solution_int_r[num] << "\n";
//
//		std::cout << "i " << temp -> index[0] << " j "<< temp -> index[1]<< "\n";
//		for(auto& v : temp -> solution_int_r){
//
//			std::cout<< v << " ";
//
//		}
//		std::cout << "\n";
//	}
//
//}
	}


}
