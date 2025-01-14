#include "dg_prefinement.h"
#include <vector>
#include "dg_unit.h"
#include "dg_param.h"
#include "dg_interpolate_to_new_points.h"
#include "dg_nodal_2d_storage.h"
#include <unordered_map>
#include "dg_local_storage.h"
#include "dg_status_table.h"
#include "dg_cantor_pairing.h"	
#include <cassert>
#include "dg_single_index.h"
#include <iostream>	// test

// forward declaration----------------------------------------
void p_refinement_apply(Unit* temp);

void p_refinement(int kt);

void p_coarsening_interpolate(Unit* temp);

void p_coarsening_L2(Unit* temp);

void Update_neighbours_facen(Unit* temp);
//------------------------------------------------------------

/// @brief
/// Apply p-refinement to the elements
void p_refinement(int kt){

	Unit* temp = local::head;

	for(int k = 0; k < local::local_elem_num; ++k){

		// test==============================================================================
//		long long int key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);
//
//		if(key == 16 && kt == 0){
//
//			p_refinement_apply(temp);
//
//		}
		//====================================================================================

		if(temp -> prefine){	// if flag refine and not exceed the maximum poly order

			p_refinement_apply(temp);

			temp -> prefine = false;	// turn off the switch after refine the element

		}

		temp = temp -> next;
	}

}

/// @brief
/// Increase polynomial order with 2 if does not exceed the maximum polynomial order. 
/// Here we increase polynomial order in x and y direction together.
/// @param temp Pointer to the current element.
void p_refinement_apply(Unit* temp){

	int new_order = temp -> n + 2;

	// impose p-refinement if new_order is inside the range
	if(new_order <= grid::nmax){

		// interpolate solution to new GL points
		// form interpolation matrix in x and y direction (assume refinement in x and y driection)-------------
		std::vector<double> T;	// only need one interpolation matrix (porderx = pordery)
		Polynomial_interpolate_matrix(nodal::gl_points[temp -> n], nodal::gl_points[new_order], T);
		// ----------------------------------------------------------------------------------------------------

		std::unordered_map<int, std::vector<double>> middle;	// intermidiate matrix

		// y direction
		for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){
			
			int start_old{};
			int start_new{};

			middle[equ] = std::vector<double> ((new_order + 1) * (temp -> n + 1));
			
			for(int i = 0; i <= (temp -> n); ++i){
	
				// interval == 1 since we are in teh y direction
				// restriction: children inderit parent's polynomial order. 
				Interpolate_to_new_points(new_order + 1, temp -> m + 1, T,
						temp -> solution[equ], middle[equ], start_old, start_new, 1);

				start_old += (temp -> m + 1);
				start_new += (new_order + 1);
			}

		}


		// clear the old solutions
		temp -> solution.clear();

		// x direction
		for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){
	
			int start_old{};
			int start_new{};
			
			// allocate space
			temp -> solution[equ] = std::vector<double> ((new_order + 1) * (new_order + 1));

			for(int j = 0; j <= new_order; ++j){
	
				Interpolate_to_new_points(new_order + 1,  temp -> n + 1, T,
						middle[equ], temp -> solution[equ], start_old, start_new, new_order + 1);
				
				++start_old;
				++start_new;
			}

		}

		//need to interpolate nodes, metrics and jacobian using the above method. 
		//need to interpolate boundary normals and scaling factors using 1d interpolation

		
		//Nodes
		//works!	
		coarsetofine(temp->holdmetrics.x_node, new_order, temp->n, T);
		coarsetofine(temp->holdmetrics.y_node, new_order, temp->n, T);


		//Metrics
		coarsetofine(temp->holdmetrics.delx_delxi, new_order, temp->n, T);
		coarsetofine(temp->holdmetrics.delx_deleta, new_order, temp->n, T);
		coarsetofine(temp->holdmetrics.dely_delxi, new_order, temp->n, T);
		coarsetofine(temp->holdmetrics.dely_deleta, new_order, temp->n, T);

		//Jacobian
		coarsetofine(temp->holdmetrics.jacobian, new_order, temp->n, T);


		coarsetofine1D(temp->holdmetrics.x_boundary, new_order, temp->n, T);
		// std::cout << "new x_boundary size: " << temp->holdmetrics.x_boundary.size() << std::endl;
		coarsetofine1D(temp->holdmetrics.y_boundary, new_order, temp->n, T);

		//Boundary Normals
		coarsetofine1D(temp->holdmetrics.boundary_normal, new_order, temp->n, T);
		coarsetofine1D(temp->holdmetrics.boundary_normal_y, new_order, temp->n, T);
		//Scaling Factors
		coarsetofine1D(temp->holdmetrics.scaling_factor, new_order, temp->n, T);

		// update the polynomial order
		temp -> n = new_order;
		temp -> m = new_order;

		// updates local neighbours face info
		Update_neighbours_facen(temp);

	}	

	temp -> prefine = false;	// turn off the switch

}

/// @brief
/// Decrease polynomial order by 2. Use L2 projection. 
/// @param temp pointer to the current element.
/// @note Assuming the polynomial orders are identical in x and y direction. 
void p_coarsening_L2(Unit* temp){

	std::cout << "P Course" << std::endl;

	assert((temp -> n - 2 >= grid::nmin) && 
			"The functionally coarsening cannot be less than the predefined minimum plonomial order ");

	assert((temp -> n == temp -> m) && "Only works for poly order equals in x and y direction. ");

	int n = temp -> n;
	int new_order = temp -> n - 2;	// assume the order in x and y direction are equal
	
	// interpolate solution to new GL points
	// form interpolation matrix in x and y direction (assume refinement in x and y driection)-------------
	std::vector<double> T;	// only need one interpolation matrix (porderx = pordery)
	Polynomial_interpolate_matrix(nodal::gl_points[new_order], nodal::gl_points[n], T); // use the transpose late
	// ----------------------------------------------------------------------------------------------------

	std::unordered_map<int, std::vector<double>> middle;	// intermidiate matrix

	// y direciton
	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

		middle[equ] = std::vector<double>((new_order + 1) * (n + 1));

		for(int xi = 0; xi <= n; ++xi ){	// loop in x direction

			for(int i = 0; i <= new_order; ++i ){

				int nodep = Get_single_index(xi, i, new_order + 1);

				for(int j = 0; j <= n; ++j){

					int nodei = Get_single_index(j, i, new_order + 1);

					int nodec = Get_single_index(xi, j, n + 1);

					middle[equ][nodep] +=  T[nodei] * 
								(nodal::gl_weights[n][j] / nodal::gl_weights[new_order][i])
								* (temp -> solution[equ][nodec]);
				}
			}

		}
	}

	// clear old solution 
	temp -> solution.clear();



	// x dir
	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

		temp -> solution[equ] = std::vector<double>((new_order + 1) * (new_order + 1));

		for(int yi = 0; yi <= new_order; ++yi){

			for(int i = 0; i <= new_order; ++i){

				int nodep = Get_single_index(i, yi, new_order + 1);

				for(int j = 0; j <= n; ++j){
				
					int nodei = Get_single_index(j, i, new_order + 1);

					int nodec = Get_single_index(j, yi, new_order + 1);

					temp -> solution[equ][nodep] += T[nodei] * 
								(nodal::gl_weights[n][j] / nodal::gl_weights[new_order][i])
								* (middle[equ][nodec]);

				}
			}

		}
		

	}

	//x_node
	// y dir
	std::vector<double> interarr = std::vector<double>((new_order + 1) * (n + 1));

	for(int xi = 0; xi <= n; ++xi ){	// loop in x direction

		for(int i = 0; i <= new_order; ++i ){

			int nodep = Get_single_index(xi, i, new_order + 1);

			for(int j = 0; j <= n; ++j){

				int nodei = Get_single_index(j, i, new_order + 1);

				int nodec = Get_single_index(xi, j, n + 1);

				interarr[nodep] +=  T[nodei] * 
							(nodal::gl_weights[n][j] / nodal::gl_weights[new_order][i])
							* (temp -> holdmetrics.x_node[nodec]);
			}
		}

	}

	temp->holdmetrics.x_node.clear();

	temp ->holdmetrics.x_node = std::vector<double>((new_order + 1) * (new_order + 1));

	//x dir
	for(int yi = 0; yi <= new_order; ++yi){

		for(int i = 0; i <= new_order; ++i){

			int nodep = Get_single_index(i, yi, new_order + 1);

			for(int j = 0; j <= n; ++j){
			
				int nodei = Get_single_index(j, i, new_order + 1);

				int nodec = Get_single_index(j, yi, new_order + 1);

				temp ->holdmetrics.x_node[nodep] += T[nodei] * 
							(nodal::gl_weights[n][j] / nodal::gl_weights[new_order][i])
							* (interarr[nodec]);

			}
		}

	}
	
	//y_node
	// y dir
	interarr.clear();
	interarr = std::vector<double>((new_order + 1) * (n + 1));

	for(int xi = 0; xi <= n; ++xi ){	// loop in x direction

		for(int i = 0; i <= new_order; ++i ){

			int nodep = Get_single_index(xi, i, new_order + 1);

			for(int j = 0; j <= n; ++j){

				int nodei = Get_single_index(j, i, new_order + 1);

				int nodec = Get_single_index(xi, j, n + 1);

				interarr[nodep] +=  T[nodei] * 
							(nodal::gl_weights[n][j] / nodal::gl_weights[new_order][i])
							* (temp -> holdmetrics.y_node[nodec]);
			}
		}

	}

	temp->holdmetrics.y_node.clear();

	temp ->holdmetrics.y_node = std::vector<double>((new_order + 1) * (new_order + 1));

	//x dir
	for(int yi = 0; yi <= new_order; ++yi){

		for(int i = 0; i <= new_order; ++i){

			int nodep = Get_single_index(i, yi, new_order + 1);

			for(int j = 0; j <= n; ++j){
			
				int nodei = Get_single_index(j, i, new_order + 1);

				int nodec = Get_single_index(j, yi, new_order + 1);

				temp ->holdmetrics.y_node[nodep] += T[nodei] * 
							(nodal::gl_weights[n][j] / nodal::gl_weights[new_order][i])
							* (interarr[nodec]);

			}
		}

	}

	//delx_delxi
	// y dir
	interarr.clear();
	interarr = std::vector<double>((new_order + 1) * (n + 1));

	for(int xi = 0; xi <= n; ++xi ){	// loop in x direction

		for(int i = 0; i <= new_order; ++i ){

			int nodep = Get_single_index(xi, i, new_order + 1);

			for(int j = 0; j <= n; ++j){

				int nodei = Get_single_index(j, i, new_order + 1);

				int nodec = Get_single_index(xi, j, n + 1);

				interarr[nodep] +=  T[nodei] * 
							(nodal::gl_weights[n][j] / nodal::gl_weights[new_order][i])
							* (temp -> holdmetrics.delx_delxi[nodec]);
			}
		}

	}

	temp->holdmetrics.delx_delxi.clear();

	temp ->holdmetrics.delx_delxi = std::vector<double>((new_order + 1) * (new_order + 1));

	//x dir
	for(int yi = 0; yi <= new_order; ++yi){

		for(int i = 0; i <= new_order; ++i){

			int nodep = Get_single_index(i, yi, new_order + 1);

			for(int j = 0; j <= n; ++j){
			
				int nodei = Get_single_index(j, i, new_order + 1);

				int nodec = Get_single_index(j, yi, new_order + 1);

				temp ->holdmetrics.delx_delxi[nodep] += T[nodei] * 
							(nodal::gl_weights[n][j] / nodal::gl_weights[new_order][i])
							* (interarr[nodec]);

			}
		}

	}

	//delx_deleta
	// y dir
	interarr.clear();
	interarr = std::vector<double>((new_order + 1) * (n + 1));

	for(int xi = 0; xi <= n; ++xi ){	// loop in x direction

		for(int i = 0; i <= new_order; ++i ){

			int nodep = Get_single_index(xi, i, new_order + 1);

			for(int j = 0; j <= n; ++j){

				int nodei = Get_single_index(j, i, new_order + 1);

				int nodec = Get_single_index(xi, j, n + 1);

				interarr[nodep] +=  T[nodei] * 
							(nodal::gl_weights[n][j] / nodal::gl_weights[new_order][i])
							* (temp -> holdmetrics.delx_deleta[nodec]);
			}
		}

	}

	temp->holdmetrics.delx_deleta.clear();

	temp ->holdmetrics.delx_deleta = std::vector<double>((new_order + 1) * (new_order + 1));

	//x dir
	for(int yi = 0; yi <= new_order; ++yi){

		for(int i = 0; i <= new_order; ++i){

			int nodep = Get_single_index(i, yi, new_order + 1);

			for(int j = 0; j <= n; ++j){
			
				int nodei = Get_single_index(j, i, new_order + 1);

				int nodec = Get_single_index(j, yi, new_order + 1);

				temp ->holdmetrics.delx_deleta[nodep] += T[nodei] * 
							(nodal::gl_weights[n][j] / nodal::gl_weights[new_order][i])
							* (interarr[nodec]);

			}
		}

	}

	//dely_delxi
	// y dir
	interarr.clear();
	interarr = std::vector<double>((new_order + 1) * (n + 1));

	for(int xi = 0; xi <= n; ++xi ){	// loop in x direction

		for(int i = 0; i <= new_order; ++i ){

			int nodep = Get_single_index(xi, i, new_order + 1);

			for(int j = 0; j <= n; ++j){

				int nodei = Get_single_index(j, i, new_order + 1);

				int nodec = Get_single_index(xi, j, n + 1);

				interarr[nodep] +=  T[nodei] * 
							(nodal::gl_weights[n][j] / nodal::gl_weights[new_order][i])
							* (temp -> holdmetrics.dely_delxi[nodec]);
			}
		}

	}

	temp->holdmetrics.dely_delxi.clear();

	temp ->holdmetrics.dely_delxi = std::vector<double>((new_order + 1) * (new_order + 1));

	//x dir
	for(int yi = 0; yi <= new_order; ++yi){

		for(int i = 0; i <= new_order; ++i){

			int nodep = Get_single_index(i, yi, new_order + 1);

			for(int j = 0; j <= n; ++j){
			
				int nodei = Get_single_index(j, i, new_order + 1);

				int nodec = Get_single_index(j, yi, new_order + 1);

				temp ->holdmetrics.dely_delxi[nodep] += T[nodei] * 
							(nodal::gl_weights[n][j] / nodal::gl_weights[new_order][i])
							* (interarr[nodec]);

			}
		}

	}

	//dely_deleta
	// y dir
	interarr.clear();
	interarr = std::vector<double>((new_order + 1) * (n + 1));

	for(int xi = 0; xi <= n; ++xi ){	// loop in x direction

		for(int i = 0; i <= new_order; ++i ){

			int nodep = Get_single_index(xi, i, new_order + 1);

			for(int j = 0; j <= n; ++j){

				int nodei = Get_single_index(j, i, new_order + 1);

				int nodec = Get_single_index(xi, j, n + 1);

				interarr[nodep] +=  T[nodei] * 
							(nodal::gl_weights[n][j] / nodal::gl_weights[new_order][i])
							* (temp -> holdmetrics.dely_deleta[nodec]);
			}
		}

	}

	temp->holdmetrics.dely_deleta.clear();

	temp ->holdmetrics.dely_deleta = std::vector<double>((new_order + 1) * (new_order + 1));

	//x dir
	for(int yi = 0; yi <= new_order; ++yi){

		for(int i = 0; i <= new_order; ++i){

			int nodep = Get_single_index(i, yi, new_order + 1);

			for(int j = 0; j <= n; ++j){
			
				int nodei = Get_single_index(j, i, new_order + 1);

				int nodec = Get_single_index(j, yi, new_order + 1);

				temp ->holdmetrics.dely_deleta[nodep] += T[nodei] * 
							(nodal::gl_weights[n][j] / nodal::gl_weights[new_order][i])
							* (interarr[nodec]);

			}
		}

	}

	//jacobian
	// y dir
	interarr.clear();
	interarr = std::vector<double>((new_order + 1) * (n + 1));

	for(int xi = 0; xi <= n; ++xi ){	// loop in x direction

		for(int i = 0; i <= new_order; ++i ){

			int nodep = Get_single_index(xi, i, new_order + 1);

			for(int j = 0; j <= n; ++j){

				int nodei = Get_single_index(j, i, new_order + 1);

				int nodec = Get_single_index(xi, j, n + 1);

				interarr[nodep] +=  T[nodei] * 
							(nodal::gl_weights[n][j] / nodal::gl_weights[new_order][i])
							* (temp -> holdmetrics.jacobian[nodec]);
			}
		}

	}

	temp->holdmetrics.jacobian.clear();

	temp ->holdmetrics.jacobian = std::vector<double>((new_order + 1) * (new_order + 1));

	//x dir
	for(int yi = 0; yi <= new_order; ++yi){

		for(int i = 0; i <= new_order; ++i){

			int nodep = Get_single_index(i, yi, new_order + 1);

			for(int j = 0; j <= n; ++j){
			
				int nodei = Get_single_index(j, i, new_order + 1);

				int nodec = Get_single_index(j, yi, new_order + 1);

				temp ->holdmetrics.jacobian[nodep] += T[nodei] * 
							(nodal::gl_weights[n][j] / nodal::gl_weights[new_order][i])
							* (interarr[nodec]);

			}
		}

	}

	//Boundary Points
	coarsetofine1D(temp->holdmetrics.x_boundary, new_order, temp->n, T);
	coarsetofine1D(temp->holdmetrics.y_boundary, new_order, temp->n, T);

	//Boundary Normals
	coarsetofine1D(temp->holdmetrics.boundary_normal, new_order, temp->n, T);
	coarsetofine1D(temp->holdmetrics.boundary_normal_y, new_order, temp->n, T);
	
	//Scaling Factors
	coarsetofine1D(temp->holdmetrics.scaling_factor, new_order, temp->n, T);


	// update the polynomial order
	temp -> n = new_order;
	temp -> m = new_order;

	// updates local neighbours face info
	Update_neighbours_facen(temp);

	temp -> coarsen = false;
}


/// @brief
/// Decrease polynomial order by 2. Use lagrange interpolation. 
/// @param temp pointer to the current element.
/// @note Assuming the polynomial orders are identical in x and y direction. 
void p_coarsening_interpolate(Unit* temp){

	std::cout << "P COARSE" << std::endl;

	assert((temp -> n - 2 >= grid::nmin) && 
			"The functionally coarsening cannot be less than the predefined minimum plonomial order ");

	assert((temp -> n == temp -> m) && "Only works for poly order equals in x and y direction. ");

	int new_order = temp -> n - 2;	// assume the order in x and y direction are equal

	// interpolate solution to new GL points
	// form interpolation matrix in x and y direction (assume refinement in x and y driection)-------------
	std::vector<double> T;	// only need one interpolation matrix (porderx = pordery)
	Polynomial_interpolate_matrix(nodal::gl_points[temp -> n], nodal::gl_points[new_order], T);
	// ----------------------------------------------------------------------------------------------------

	std::unordered_map<int, std::vector<double>> middle;	// intermidiate matrix

	// y direction
	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){
		
		int start_old{};
		int start_new{};

		middle[equ] = std::vector<double> ((new_order + 1) * (temp -> n + 1));
		
		for(int i = 0; i <= (temp -> n); ++i){

			// interval == 1 since we are in teh y direction
			// restriction: children inderit parent's polynomial order. 
			Interpolate_to_new_points(new_order + 1, temp -> m + 1, T,
					temp -> solution[equ], middle[equ], start_old, start_new, 1);

			start_old += (temp -> m + 1);
			start_new += (new_order + 1);
		}

	}


	// clear the old solutions
	temp -> solution.clear();

	// x direction
	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

		int start_old{};
		int start_new{};
		
		// allocate space
		temp -> solution[equ] = std::vector<double> ((new_order + 1) * (new_order + 1));

		for(int j = 0; j <= new_order; ++j){

			Interpolate_to_new_points(new_order + 1,  temp -> n + 1, T,
					middle[equ], temp -> solution[equ], start_old, start_new, new_order + 1);
			
			++start_old;
			++start_new;
		}

	}

	//Nodes
	coarsetofine(temp->holdmetrics.x_node, new_order, temp->n, T);
	coarsetofine(temp->holdmetrics.y_node, new_order, temp->n, T);


	//Metrics
	coarsetofine(temp->holdmetrics.delx_delxi, new_order, temp->n, T);
	coarsetofine(temp->holdmetrics.delx_deleta, new_order, temp->n, T);
	coarsetofine(temp->holdmetrics.dely_delxi, new_order, temp->n, T);
	coarsetofine(temp->holdmetrics.dely_deleta, new_order, temp->n, T);

	//Jacobian
	coarsetofine(temp->holdmetrics.jacobian, new_order, temp->n, T);

	//Boundaries
	coarsetofine1D(temp->holdmetrics.x_boundary, new_order, temp->n, T);
	coarsetofine1D(temp->holdmetrics.y_boundary, new_order, temp->n, T);

	//Boundary Normals
	coarsetofine1D(temp->holdmetrics.boundary_normal, new_order, temp->n, T);
	coarsetofine1D(temp->holdmetrics.boundary_normal_y, new_order, temp->n, T);
	
	//Scaling Factors
	coarsetofine1D(temp->holdmetrics.scaling_factor, new_order, temp->n, T);

	// update the polynomial order
	temp -> n = new_order;
	temp -> m = new_order;

	// updates local neighbours face info
	Update_neighbours_facen(temp);

	temp -> coarsen = false;	// turn off the coarsening switch

}


/// @brief
/// Since the current element's polynomial order changes, we need to update it's neighbours records.
/// @param temp Pointer to the current element.
void Update_neighbours_facen(Unit* temp){

	long long int my_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);

	// loop through four faces
	for(int i = 0; i < 4; ++i){

		auto it_face = temp -> facen[i].begin();
		
		for(; it_face != temp -> facen[i].end(); ++it_face){
	
			int opp_dir = Opposite_dir(i);

			// only change local neighbours
			if(it_face -> face_type == 'L'){

				// go to neighbour and updates the orders
				for(auto it_n = local::Hash_elem[it_face -> key] -> facen[opp_dir].begin();
					it_n != local::Hash_elem[it_face -> key] -> facen[opp_dir].end(); ++it_n ){

					if(it_n -> key == my_key){	// find it

						it_n -> porderx = temp -> n;
						it_n -> pordery = temp -> m;
					}
				}
			}
			
		}

	}

}
