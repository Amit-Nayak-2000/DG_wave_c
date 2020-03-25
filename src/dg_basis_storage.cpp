#include <vector>
#include "dg_nodal_2d_storage.h"
#include "dg_param.h"
#include "dg_basis.h"
#include "dg_basis_storage.h"
#include <iostream>	// test


/// @brief
/// Processor's local storages.
/// Generate the basis information (GL PIOINTS, WEIGHTS, DERIVATIVE 
/// MATRIX, etc.) for all polynomial level at the beginning,
/// and store it in hash tables. 
/// Therefore, after refinement start we do not need to 
/// re-generate the basis information, or store duplicated data. 
void Construct_basis_storage(){


	// generate the 2d storages
	for(int n = grid::nmin; n <= grid::nmax; n+=2){	// k = poly order
		
		// allocate space	
		nodal::gl_points[n] = std::vector<double>(n + 1);
		nodal::gl_weights[n] = std::vector<double>(n + 1);
	
		int der_size = (n + 1) * (n + 1);
		nodal::first_der[n] = std::vector<double>(der_size);
	
		// generate current gl_p, gl_w
		GL(n, nodal::gl_points[n], nodal::gl_weights[n]);
	
		std::vector<double> bary(n + 1);
	
		BARW(n, nodal::gl_points[n], bary);
	
		// first order derivative matrix
		Mth_order_polynomial_derivative_matrix(n, 1, nodal::gl_points[n], nodal::first_der[n], bary);

		// Lagrange interpolates on the boundaries. 
		nodal::lagrange_l[n] = std::vector<double>(n + 1);
		nodal::lagrange_r[n] = std::vector<double>(n + 1);
		Lagrange_interpolating_polynomial(n, -1.0, nodal::gl_points[n], bary, nodal::lagrange_l[n]);
		Lagrange_interpolating_polynomial(n,  1.0, nodal::gl_points[n], bary, nodal::lagrange_r[n]);

	}
	
}

