#include "dg_interpolate_to_new_points.h"
#include "dg_local_storage.h"
#include "dg_unit.h"
#include "dg_affine_map.h"
#include "dg_nodal_2d_storage.h"
#include "dg_single_index.h"
#include "dg_single_index.h"
#include <vector>
#include "dg_param.h"
#include "dg_basis.h"
#include <iostream>	// test

// forward declaration------------------------------------------------------------------------------------------
void Polynomial_interpolate_matrix(std::vector<double>& x, std::vector<double>& xi, std::vector<double>& T);

void Solutions_to_child(int key, int p_key);

void Interpolate_to_new_points(int m, int n, std::vector<double>& T, 
				std::vector<double>& f, std::vector<double>& new_f, int start, int interval);
//--------------------------------------------------------------------------------------------------------------


/// @brief
///
/// @param key child's key.
/// @parma p_key parent's key. 
void Solutions_to_child(int key, int p_key){

	Unit* temp = local::Hash_elem[key];	// pointer to this elem

	double y_start{};

	if((temp -> child_position) == 0 || (temp -> child_position) == 1){

		y_start = -1.0;

	}
	else{

		y_start = 0.0;
	}

	for(int i = 0; i <= (temp -> n); ++i){

		std::vector<double> y(temp -> m + 1);	// location on parent

		// generate new sets of points
		for(int j = 0; j <= (temp -> m); ++j){
			// delta_y == 1.0 
			y[j] = Affine_mapping(nodal::gl_points[temp -> m][j], y_start, 1.0);
//std::cout<< "c_key "<< key << " j "<< j << " points " << y[j] << "\n";
		}

		// interpolate to new sets of points
		std::vector<double> T; 	// interpolation matrix
		// note: children porder == parent porder, so old points use child's poly order
		Polynomial_interpolate_matrix(nodal::gl_points[temp -> m], y, T);

		int start = Get_single_index(i, 0, temp -> m + 1);

		for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){
			// allocate space
			temp -> solution[equ] = std::vector<double> ((temp -> n + 1) * (temp -> m + 1));

			// interval == 1 since we are in teh y direction
			// restriction: children inderit parent's polynomial order. 
			Interpolate_to_new_points(temp -> m + 1, temp -> m + 1, T,
					local::Hash_elem[p_key] -> solution[equ], temp -> solution[equ], start, 1);
		}
	}

}

/// @brief
/// Matrix for inerpolation between two sets of points
/// @parma x Old set of points. Size [0, n].
/// @param xi New set of points. Size [0, m]. 
/// @param T The interpolation matrix. Size [0, m] * [0, n]. Do not need to pre-allocate space. 
void Polynomial_interpolate_matrix(std::vector<double>& x, std::vector<double>& xi, std::vector<double>& T){

	int m = xi.size();
	int n = x.size();

	T = std::vector<double>(m * n);	// allocate space

	std::vector<double> w(n);	// Barycentric weights for old sets of points. Size [0, n].

	BARW(n - 1, x, w);

//for(auto& v : w){
//	std::cout<< "porder " << n - 1 << " bary "<< v << "\n";
//}
//std::cout << " ------------------ \n";


	for(int k = 0; k < m; ++k){

		bool row_match = false;

		for(int j = 0; j < n; ++j){

			bool same_point = Almost_equal(x[j], xi[k]);

			if(same_point){	// new point is one of the collocation points

				row_match = true;

				int index = Get_single_index(k, j, n);
				T[index] = 1.0;
			}

		}

		if(!row_match){	// if not the same point

			double s{};

			for(int j = 0; j < n; ++j){

				double t = w[j] / (xi[k] - x[j]);

				int index = Get_single_index(k, j, n);

				T[index] = t;

				s += t;
//std::cout << "k" << k << " j "<< j<< " weight "<< w[j]<< " x[j] " << x[j] << " t = "<< t << " s = "<< s << "\n";
			}

			for(int j = 0; j < n; ++j){

				int index = Get_single_index(k, j, n);

				T[index] /= s;
				

			}
			

		}
	}

//for(int i = 0; i < m; ++i){
//
//
//	for(int j = 0; j < n; ++j){
//
//		int index = Get_single_index(i, j, n);
//		std::cout<< i << " " << j << " "<< T[index] << "\n";
//
//	}
//	
//}
//std::cout << " ========================= \n";

}

/// @brief
/// Interpolation between two sets of points by matrix multiplicaiton
/// @param m The number of new set of points.
/// @param n The number of old set of points.	
/// @param T Interpolating matrix. Size m * n.
/// @param f Solutions on the old set of points.
/// @param new_f Solusion on the new set of points. Size m. 
/// @param start The node number of the first value in new_f.
/// @param interval node index interval.
void Interpolate_to_new_points(int m, int n, std::vector<double>& T, 
				std::vector<double>& f, std::vector<double>& new_f, int start, int interval){

	int ni = start;	
	for(int i = 0; i < m; ++i){

		double t{};

		int fi = start;

		for(int j = 0; j < n; ++j){

			int index = Get_single_index(i, j, n);

			t += T[index] * f[fi];

			fi += interval;
		}

		new_f[ni] = t;
		ni += interval;

	}
}