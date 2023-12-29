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
#include <array>
#include <unordered_map>
#include <iostream>	// test

// forward declaration------------------------------------------------------------------------------------------
void Polynomial_interpolate_matrix(std::vector<double>& x, std::vector<double>& xi, std::vector<double>& T);

void Solution_back_to_parent(std::array<long long int, 4>& keys, long long int p_key);

void Solutions_to_children(std::array<long long int, 4>& keys, long long int p_key);

void Interpolate_to_new_points(int m, int n, std::vector<double>& T, 
				std::vector<double>& f, std::vector<double>& new_f, int start_old, int start_new, int interval);

void Form_new_set_of_points(int m, int start, std::vector<double>& y);

void Two_dir_inter(long long int p_key, Unit* c, std::vector<double>& T_x, std::vector<double>& T_y, int n, int m);

void Two_direc_interp(long long int p_key, Unit* c, std::vector<double>& T_x, std::vector<double>& T_y, int n, int m);

void Mortar_inter_back(Unit* c, Unit* p, std::vector<double>& Ty, std::vector<double>& Tx, double b);

void Mortar_inter_back_map(Unit* c, Unit* p, std::vector<double>& Ty, std::vector<double>& Tx, double b);

void Lagrange_inter_back(Unit* c, Unit* p, std::vector<double>& Ty, std::vector<double>& Tx, double b);

void Lagrange_inter_back_transpose(Unit* c, Unit* p, std::vector<double>& Ty, std::vector<double>& Tx, double b);

void coarsetofine(std::vector<double> &target, int new_order, int currentorder, std::vector<double> &T);

void coarsetofine1D(std::vector<double> &target, int new_order, int currentorder, const std::vector<double> &T);

void CoursetoFineInterp(int nold, int nnew, std::vector<double> &x, std::vector<double> &y, std::unordered_map<int, std::vector<double>> &f, 
std::vector<double> &xnew,  std::vector<double> &ynew, std::unordered_map<int, std::vector<double>> &fnew);

std::vector<double> linarray(double min, double max, int N);

void boundaryinterp1D(Unit* c, int n, int m);

void CoursetoFineInterpCoords(int nold, int nnew, std::vector<double> &x, std::vector<double> &y, std::unordered_map<int, std::vector<double>> &f, 
std::vector<double> &xnew,  std::vector<double> &ynew, std::unordered_map<int, std::vector<double>> &fnew);
//--------------------------------------------------------------------------------------------------------------


/// @brief
/// Interpolate the solution from parent to four children.
/// @param keys four children key. Sequence: SW -- NW -- NE -- SE
/// @parma p_key parent's key. 
/// @node this function can take care of different polynomial order in x and y direction. 
void Solutions_to_children(std::array<long long int, 4>& keys, long long int p_key){


	// pointers to the four children
	Unit* c0 = local::Hash_elem[keys[0]];
	Unit* c1 = local::Hash_elem[keys[1]];
	Unit* c2 = local::Hash_elem[keys[2]];
	Unit* c3 = local::Hash_elem[keys[3]];

	// poly order (four children are the same)
	int n = c0 -> n;
	int m = c0 -> m;

	// allocate solution space
	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){
		// allocate space
		c0 -> solution[equ] = std::vector<double> ((n + 1) * ( m + 1));
		c1 -> solution[equ] = std::vector<double> ((n + 1) * ( m + 1));
		c2 -> solution[equ] = std::vector<double> ((n + 1) * ( m + 1));
		c3 -> solution[equ] = std::vector<double> ((n + 1) * ( m + 1));

	}

	// four new sets of points---------------------------------------------------
	std::vector<double> xl(n + 1);	// location on parent
	std::vector<double> xr(n + 1);	// location on parent
	std::vector<double> yl(m + 1);	// location on parent
	std::vector<double> yr(m + 1);	// location on parent

	Form_new_set_of_points(n, -1.0, xl);	// left boundary is -1.0, right boundary is 0. 
	Form_new_set_of_points(n,  0.0, xr);
	Form_new_set_of_points(m, -1.0, yl);
	Form_new_set_of_points(m,  0.0, yr);
	//---------------------------------------------------------------------------

	// form interpolation matrix------------------------------------------------------
	std::vector<double> T_xl; 	// interpolation matrix
	std::vector<double> T_xr; 	// interpolation matrix
	std::vector<double> T_yl; 	// interpolation matrix
	std::vector<double> T_yr; 	// interpolation matrix


	// note: children porder == parent porder, so old points use child's poly order
	Polynomial_interpolate_matrix(nodal::gl_points[n], xl, T_xl);
	Polynomial_interpolate_matrix(nodal::gl_points[n], xr, T_xr);
	Polynomial_interpolate_matrix(nodal::gl_points[m], yl, T_yl);
	Polynomial_interpolate_matrix(nodal::gl_points[m], yr, T_yr);
	//--------------------------------------------------------------------------------

	// c0
	Two_dir_inter(p_key, c0, T_xl, T_yl, n, m);
	Two_direc_interp(p_key, c0, T_xl, T_yl, n, m);

	// c1
	Two_dir_inter(p_key, c1, T_xr, T_yl, n, m);
	Two_direc_interp(p_key, c1, T_xr, T_yl, n, m);

	// c2
	Two_dir_inter(p_key, c2, T_xr, T_yr, n, m);
	Two_direc_interp(p_key, c2, T_xr, T_yr, n, m);

	// c3
	Two_dir_inter(p_key, c3, T_xl, T_yr, n, m);
	Two_direc_interp(p_key, c3, T_xl, T_yr, n, m);

	// std::cout << "1st child" << std::endl;
	// c0->holdmetrics.outputMetrics();

	// std::cout << "2nd child" << std::endl;
	// c1->holdmetrics.outputMetrics();

	// std::cout << "3rd child" << std::endl;
	// c2->holdmetrics.outputMetrics();

	// std::cout << "4th child" << std::endl;
	// c3->holdmetrics.outputMetrics();

	
}


/// @brief
/// Use Lagrange interpolation to interpolate solution back to parent. 
/// @param c pointer to the child.
/// @param p pointer to the parent. 
/// @param Ty interpolation matrix in y direction. 
/// @param Tx interpolation matrix in x direction. 
/// @param b scaling factor. b = 1 / number if children.  
void Lagrange_inter_back(Unit* c, Unit* p, std::vector<double>& Ty, std::vector<double>& Tx, double b){

	int n = p -> n;
	int m = n;

	std::unordered_map<int, std::vector<double>> middle;	// intermidiate matrix. <function_num, solutions>
	std::unordered_map<int, std::vector<double>> middle2;	// intermidiate matrix. 

	// y direction
	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){
		
		middle[equ] = std::vector<double> ((n + 1) * (m + 1));
		for(int i = 0; i <= n; ++i){

			int start = Get_single_index(i, 0, m + 1);

			// interval == 1 since we are in the y direction
			// restriction: children inderit parent's polynomial order. 
			Interpolate_to_new_points(m + 1,  m + 1, Ty,
					c -> solution[equ], middle[equ], start, start, 1);
		}
	}

	// x direction
	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

		middle2[equ] = std::vector<double>((n + 1) * (m + 1));

		for(int j = 0; j <= m; ++j){

			int start = Get_single_index(0, j, m + 1);

			Interpolate_to_new_points(n + 1,  n + 1, Tx,
					middle[equ], middle2[equ], start, start, m + 1);
			
		}
	}

	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

		int nodei{};

		for(auto& v : middle2[equ]){

			p -> solution[equ][nodei] += b * v;

			++nodei;
		}

	}


}


/// @brief
/// Interpolate solution form four children back to parent. 
/// Here we use L2 projection.
/// @param keys four children key. Sequence: SW -- NW -- NE -- SE
/// @parma p_key parent's key. 
/// @note polynomial order should be identical in x and y direction. 
void Solution_back_to_parent(std::array<long long int, 4>& keys, long long int p_key){

	// pointers to the four children
	Unit* c0 = local::Hash_elem[keys[0]];
	Unit* c1 = local::Hash_elem[keys[1]];
	Unit* c2 = local::Hash_elem[keys[2]];
	Unit* c3 = local::Hash_elem[keys[3]];

	// pointer to parent
	Unit* temp = local::Hash_elem[p_key] ;

	// poly order (four children are the same)
	int n = temp -> n;
	int m = temp -> m;

	// allocate space for parent 
	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

		temp -> solution[equ] = std::vector<double>((n + 1) * (m + 1));

	}

	// four new sets of points---------------------------------------------------
	std::vector<double> pl(n + 1);	// location on parent
	std::vector<double> pr(n + 1);	// location on parent

	Form_new_set_of_points(n, -1.0, pl);	// left boundary is -1.0, right boundary is 0. 
	Form_new_set_of_points(n,  0.0, pr);
	//---------------------------------------------------------------------------

	// form interpolation matrix------------------------------------------------------
	std::vector<double> Tl; 	// interpolation matrix
	std::vector<double> Tr; 	// interpolation matrix

	//--------------------------------------------------------------------------------

	// project back to parent (L2 projection) =======================================
	Polynomial_interpolate_matrix(nodal::gl_points[n], pl, Tl);
	Polynomial_interpolate_matrix(nodal::gl_points[n], pr, Tr);

	Mortar_inter_back(c0, temp, Tl, Tl, 0.25);
	Mortar_inter_back(c1, temp, Tl, Tr, 0.25);
	Mortar_inter_back(c2, temp, Tr, Tr, 0.25);
	Mortar_inter_back(c3, temp, Tr, Tl, 0.25);

	//Project Nodes Back to parent
	Mortar_inter_back_map(c0, temp, Tl, Tl, 0.25);
	Mortar_inter_back_map(c1, temp, Tl, Tr, 0.25);
	Mortar_inter_back_map(c2, temp, Tr, Tr, 0.25);
	Mortar_inter_back_map(c3, temp, Tr, Tl, 0.25);

	//Boundary Coordinate of Parent
	boundaryinterp1D(temp, n, m);

	
	std::vector<double> one_x(n+1);
	std::vector<double> one_y(n+1);
	std::vector<double> two_x(n+1);
	std::vector<double> two_y(n+1);
	std::vector<double> three_x(n+1);
	std::vector<double> three_y(n+1);
	std::vector<double> four_x(n+1);
	std::vector<double> four_y(n+1);


	for(int i = 0; i <= n; i++){
		one_x[i] = temp->holdmetrics.x_boundary[Get_single_index(i, 0, 4)];
		one_y[i] = temp->holdmetrics.y_boundary[Get_single_index(i, 0, 4)]; 

		two_x[i] = temp->holdmetrics.x_boundary[Get_single_index(i, 1, 4)];
		two_y[i] = temp->holdmetrics.y_boundary[Get_single_index(i, 1, 4)];

		three_x[i] = temp->holdmetrics.x_boundary[Get_single_index(i, 2, 4)]; 
		three_y[i] = temp->holdmetrics.y_boundary[Get_single_index(i, 2, 4)];

		four_x[i] = temp->holdmetrics.x_boundary[Get_single_index(i, 3, 4)]; 
		four_y[i] = temp->holdmetrics.y_boundary[Get_single_index(i, 3, 4)];
	}

	CurveInterpolant bot(n, nodal::gl_points[n], one_x, one_y);
	CurveInterpolant right(n, nodal::gl_points[n], two_x, two_y);
	CurveInterpolant top(n, nodal::gl_points[n], three_x, three_y);
	CurveInterpolant left(n, nodal::gl_points[n], four_x, four_y);

	CurveInterpolant coarse_B[4] = {bot, right, top, left};

	//metrics, jacobians, normals, scaling factors of Parent
	temp->holdmetrics.recalcscalenorms(coarse_B, nodal::gl_points[n]);
	// ===============================================================================

	// use lagrange interpolation =====================================================
//	Polynomial_interpolate_matrix(pl, nodal::gl_points[n], Tl);
//	Polynomial_interpolate_matrix(pr, nodal::gl_points[n], Tr);
//	Lagrange_inter_back(c0, temp, Tl, Tl, 0.25);
//	Lagrange_inter_back(c1, temp, Tl, Tr, 0.25);
//	Lagrange_inter_back(c2, temp, Tr, Tr, 0.25);
//	Lagrange_inter_back(c3, temp, Tr, Tl, 0.25);


//	Polynomial_interpolate_matrix(nodal::gl_points[n], pl, Tl);
//	Polynomial_interpolate_matrix(nodal::gl_points[n], pr, Tr);
//	Lagrange_inter_back_transpose(c0, temp, Tl, Tl, 0.25);
//	Lagrange_inter_back_transpose(c1, temp, Tl, Tr, 0.25);
//	Lagrange_inter_back_transpose(c2, temp, Tr, Tr, 0.25);
//	Lagrange_inter_back_transpose(c3, temp, Tr, Tl, 0.25);

	// ===============================================================================


}

/// @brief
/// Use Lagrange interplation to interpolate solution back to parent. Use the transpose of the 
/// forward interpolation matrix. 
/// @note does not work. 
void Lagrange_inter_back_transpose(Unit* c, Unit* p, std::vector<double>& Ty, std::vector<double>& Tx, double b){

	
	int n = p -> n; // x, y dir poly order should be the same. 

	std::unordered_map<int, std::vector<double>> middle;

	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

		middle[equ] = std::vector<double>((n + 1) * (n + 1));

		// y direction
		for(int xi = 0; xi <= n; ++xi ){	// loop in x direction

			for(int i = 0; i <= n; ++i ){

				int nodep = Get_single_index(xi, i, n + 1);

				for(int j = 0; j <= n; ++j){

					int nodei = Get_single_index(j, i, n + 1);

					int nodec = Get_single_index(xi, j, n + 1);

					middle[equ][nodep] += Ty[nodei] * (c -> solution[equ][nodec]);
				}
			}

		}

		// x dir
		for(int yi = 0; yi <= n; ++yi){

			for(int i = 0; i <= n; ++i){

				int nodep = Get_single_index(i, yi, n + 1);

				for(int j = 0; j <= n; ++j){
				
					int nodei = Get_single_index(j, i, n + 1);

					int nodec = Get_single_index(j, yi, n + 1);

					p -> solution[equ][nodep] +=  b * Tx[nodei] * (middle[equ][nodec]);

				}
			}

		}
		

	}

}


/// @param c pointer to child.
/// @param p pointer to parent.
/// @param Ty y direction interpolatio matrix. 
/// @param Tx x direction interpolatio matrix. 
void Mortar_inter_back(Unit* c, Unit* p, std::vector<double>& Ty, std::vector<double>& Tx, double b){

	int n = p -> n; // x, y dir poly order should be the same. 

	std::unordered_map<int, std::vector<double>> middle;

	// y direction
	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

		middle[equ] = std::vector<double>((n + 1) * (n + 1));

		for(int xi = 0; xi <= n; ++xi ){	// loop in x direction

			for(int i = 0; i <= n; ++i ){

				int nodep = Get_single_index(xi, i, n + 1);

				for(int j = 0; j <= n; ++j){

					int nodei = Get_single_index(j, i, n + 1);

					int nodec = Get_single_index(xi, j, n + 1);

					middle[equ][nodep] +=  Ty[nodei] * 
								(nodal::gl_weights[n][j] / nodal::gl_weights[n][i])
								* (c -> solution[equ][nodec]);
				}
			}

		}

		// x dir
		for(int yi = 0; yi <= n; ++yi){

			for(int i = 0; i <= n; ++i){

				int nodep = Get_single_index(i, yi, n + 1);

				for(int j = 0; j <= n; ++j){
				
					int nodei = Get_single_index(j, i, n + 1);

					int nodec = Get_single_index(j, yi, n + 1);

					p -> solution[equ][nodep] +=  b * Tx[nodei] * 
								(nodal::gl_weights[n][j] / nodal::gl_weights[n][i])
								* (middle[equ][nodec]);

				}
			}

		}
		

	}


}



/// @param c pointer to child.
/// @param p pointer to parent.
/// @param Ty y direction interpolatio matrix. 
/// @param Tx x direction interpolatio matrix. 
void Mortar_inter_back_map(Unit* c, Unit* p, std::vector<double>& Ty, std::vector<double>& Tx, double b){

	int n = p -> n; // x, y dir poly order should be the same. 

	std::vector<double> middle;

	//x_node
	middle = std::vector<double>((n + 1) * (n + 1));

	for(int xi = 0; xi <= n; ++xi ){	// loop in x direction

		for(int i = 0; i <= n; ++i ){

			int nodep = Get_single_index(xi, i, n + 1);

			for(int j = 0; j <= n; ++j){

				int nodei = Get_single_index(j, i, n + 1);

				int nodec = Get_single_index(xi, j, n + 1);

				middle[nodep] +=  Ty[nodei] * 
							(nodal::gl_weights[n][j] / nodal::gl_weights[n][i])
							* (c -> holdmetrics.x_node[nodec]);
			}
		}

	}

	// x dir
	for(int yi = 0; yi <= n; ++yi){

		for(int i = 0; i <= n; ++i){

			int nodep = Get_single_index(i, yi, n + 1);

			for(int j = 0; j <= n; ++j){
			
				int nodei = Get_single_index(j, i, n + 1);

				int nodec = Get_single_index(j, yi, n + 1);

				p ->holdmetrics.x_node[nodep] +=  b * Tx[nodei] * 
							(nodal::gl_weights[n][j] / nodal::gl_weights[n][i])
							* (middle[nodec]);

			}
		}

	}

	middle.clear();


	//y_node
	middle = std::vector<double>((n + 1) * (n + 1));

	for(int xi = 0; xi <= n; ++xi ){	// loop in x direction

		for(int i = 0; i <= n; ++i ){

			int nodep = Get_single_index(xi, i, n + 1);

			for(int j = 0; j <= n; ++j){

				int nodei = Get_single_index(j, i, n + 1);

				int nodec = Get_single_index(xi, j, n + 1);

				middle[nodep] +=  Ty[nodei] * 
							(nodal::gl_weights[n][j] / nodal::gl_weights[n][i])
							* (c -> holdmetrics.y_node[nodec]);
			}
		}

	}

	// x dir
	for(int yi = 0; yi <= n; ++yi){

		for(int i = 0; i <= n; ++i){

			int nodep = Get_single_index(i, yi, n + 1);

			for(int j = 0; j <= n; ++j){
			
				int nodei = Get_single_index(j, i, n + 1);

				int nodec = Get_single_index(j, yi, n + 1);

				p ->holdmetrics.y_node[nodep] +=  b * Tx[nodei] * 
							(nodal::gl_weights[n][j] / nodal::gl_weights[n][i])
							* (middle[nodec]);

			}
		}

	}

	middle.clear();




	
	// //delx_delxi
	// middle = std::vector<double>((n + 1) * (n + 1));

	// for(int xi = 0; xi <= n; ++xi ){	// loop in x direction

	// 	for(int i = 0; i <= n; ++i ){

	// 		int nodep = Get_single_index(xi, i, n + 1);

	// 		for(int j = 0; j <= n; ++j){

	// 			int nodei = Get_single_index(j, i, n + 1);

	// 			int nodec = Get_single_index(xi, j, n + 1);

	// 			middle[nodep] +=  Ty[nodei] * 
	// 						(nodal::gl_weights[n][j] / nodal::gl_weights[n][i])
	// 						* (c -> holdmetrics.delx_delxi[nodec]);
	// 		}
	// 	}

	// }

	// // x dir
	// for(int yi = 0; yi <= n; ++yi){

	// 	for(int i = 0; i <= n; ++i){

	// 		int nodep = Get_single_index(i, yi, n + 1);

	// 		for(int j = 0; j <= n; ++j){
			
	// 			int nodei = Get_single_index(j, i, n + 1);

	// 			int nodec = Get_single_index(j, yi, n + 1);

	// 			p ->holdmetrics.delx_delxi[nodep] +=  b * Tx[nodei] * 
	// 						(nodal::gl_weights[n][j] / nodal::gl_weights[n][i])
	// 						* (middle[nodec]);

	// 		}
	// 	}

	// }

	// middle.clear();


	// //delx_deleta
	// middle = std::vector<double>((n + 1) * (n + 1));

	// for(int xi = 0; xi <= n; ++xi ){	// loop in x direction

	// 	for(int i = 0; i <= n; ++i ){

	// 		int nodep = Get_single_index(xi, i, n + 1);

	// 		for(int j = 0; j <= n; ++j){

	// 			int nodei = Get_single_index(j, i, n + 1);

	// 			int nodec = Get_single_index(xi, j, n + 1);

	// 			middle[nodep] +=  Ty[nodei] * 
	// 						(nodal::gl_weights[n][j] / nodal::gl_weights[n][i])
	// 						* (c -> holdmetrics.delx_deleta[nodec]);
	// 		}
	// 	}

	// }

	// // x dir
	// for(int yi = 0; yi <= n; ++yi){

	// 	for(int i = 0; i <= n; ++i){

	// 		int nodep = Get_single_index(i, yi, n + 1);

	// 		for(int j = 0; j <= n; ++j){
			
	// 			int nodei = Get_single_index(j, i, n + 1);

	// 			int nodec = Get_single_index(j, yi, n + 1);

	// 			p ->holdmetrics.delx_deleta[nodep] +=  b * Tx[nodei] * 
	// 						(nodal::gl_weights[n][j] / nodal::gl_weights[n][i])
	// 						* (middle[nodec]);

	// 		}
	// 	}

	// }

	// middle.clear();

	// //dely_delxi
	// middle = std::vector<double>((n + 1) * (n + 1));

	// for(int xi = 0; xi <= n; ++xi ){	// loop in x direction

	// 	for(int i = 0; i <= n; ++i ){

	// 		int nodep = Get_single_index(xi, i, n + 1);

	// 		for(int j = 0; j <= n; ++j){

	// 			int nodei = Get_single_index(j, i, n + 1);

	// 			int nodec = Get_single_index(xi, j, n + 1);

	// 			middle[nodep] +=  Ty[nodei] * 
	// 						(nodal::gl_weights[n][j] / nodal::gl_weights[n][i])
	// 						* (c -> holdmetrics.dely_delxi[nodec]);
	// 		}
	// 	}

	// }

	// // x dir
	// for(int yi = 0; yi <= n; ++yi){

	// 	for(int i = 0; i <= n; ++i){

	// 		int nodep = Get_single_index(i, yi, n + 1);

	// 		for(int j = 0; j <= n; ++j){
			
	// 			int nodei = Get_single_index(j, i, n + 1);

	// 			int nodec = Get_single_index(j, yi, n + 1);

	// 			p ->holdmetrics.dely_delxi[nodep] +=  b * Tx[nodei] * 
	// 						(nodal::gl_weights[n][j] / nodal::gl_weights[n][i])
	// 						* (middle[nodec]);

	// 		}
	// 	}

	// }

	// middle.clear();

	// //dely_deleta
	// middle = std::vector<double>((n + 1) * (n + 1));

	// for(int xi = 0; xi <= n; ++xi ){	// loop in x direction

	// 	for(int i = 0; i <= n; ++i ){

	// 		int nodep = Get_single_index(xi, i, n + 1);

	// 		for(int j = 0; j <= n; ++j){

	// 			int nodei = Get_single_index(j, i, n + 1);

	// 			int nodec = Get_single_index(xi, j, n + 1);

	// 			middle[nodep] +=  Ty[nodei] * 
	// 						(nodal::gl_weights[n][j] / nodal::gl_weights[n][i])
	// 						* (c -> holdmetrics.dely_deleta[nodec]);
	// 		}
	// 	}

	// }

	// // x dir
	// for(int yi = 0; yi <= n; ++yi){

	// 	for(int i = 0; i <= n; ++i){

	// 		int nodep = Get_single_index(i, yi, n + 1);

	// 		for(int j = 0; j <= n; ++j){
			
	// 			int nodei = Get_single_index(j, i, n + 1);

	// 			int nodec = Get_single_index(j, yi, n + 1);

	// 			p ->holdmetrics.dely_deleta[nodep] +=  b * Tx[nodei] * 
	// 						(nodal::gl_weights[n][j] / nodal::gl_weights[n][i])
	// 						* (middle[nodec]);

	// 		}
	// 	}

	// }

	// middle.clear();

	// //jacobian
	// middle = std::vector<double>((n + 1) * (n + 1));

	// for(int xi = 0; xi <= n; ++xi ){	// loop in x direction

	// 	for(int i = 0; i <= n; ++i ){

	// 		int nodep = Get_single_index(xi, i, n + 1);

	// 		for(int j = 0; j <= n; ++j){

	// 			int nodei = Get_single_index(j, i, n + 1);

	// 			int nodec = Get_single_index(xi, j, n + 1);

	// 			middle[nodep] +=  Ty[nodei] * 
	// 						(nodal::gl_weights[n][j] / nodal::gl_weights[n][i])
	// 						* (c -> holdmetrics.jacobian[nodec]);
	// 		}
	// 	}

	// }

	// // x dir
	// for(int yi = 0; yi <= n; ++yi){

	// 	for(int i = 0; i <= n; ++i){

	// 		int nodep = Get_single_index(i, yi, n + 1);

	// 		for(int j = 0; j <= n; ++j){
			
	// 			int nodei = Get_single_index(j, i, n + 1);

	// 			int nodec = Get_single_index(j, yi, n + 1);

	// 			p ->holdmetrics.jacobian[nodep] +=  b * Tx[nodei] * 
	// 						(nodal::gl_weights[n][j] / nodal::gl_weights[n][i])
	// 						* (middle[nodec]);

	// 		}
	// 	}

	// }

	// middle.clear();

}




/// @brief
/// Use the interpolation matrix in x and y direction and obtain the interpolated solutions.
/// @param p_key parents key.
/// @param c pointer to the current child.
/// @param T_x interpolation matrix in x direction.
/// @param T_y interpolation matrix in y direction.
/// @param n polynomial order in x direction.
/// @param m polynomial order in y direction.
void Two_dir_inter(long long int p_key, Unit* c, std::vector<double>& T_x, std::vector<double>& T_y, int n, int m){

	std::unordered_map<int, std::vector<double>> middle;	// intermidiate matrix. 

	// y direction
	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){
		
		middle[equ] = std::vector<double> ((n + 1) * (m + 1));
		
		for(int i = 0; i <= n; ++i){

			int start = Get_single_index(i, 0, m + 1);

			// interval == 1 since we are in teh y direction
			// restriction: children inderit parent's polynomial order. 
			Interpolate_to_new_points(m + 1,  m + 1, T_y,
					local::Hash_elem[p_key] -> solution[equ], middle[equ], start, start, 1);
		}
	}

	// x direction
	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

		for(int j = 0; j <= m; ++j){

			int start = Get_single_index(0, j, m + 1);

			Interpolate_to_new_points(n + 1,  n + 1, T_x,
					middle[equ], c -> solution[equ], start, start, m + 1);
			
		}
	}


}

//Interpolates nodes, metric terms and jacobians.
void Two_direc_interp(long long int p_key, Unit* c, std::vector<double>& T_x, std::vector<double>& T_y, int n, int m){

	//Interpolate Nodes
	std::vector<double> middle;	// intermidiate matrix. 

	//x_node:
	// y direction
	middle = std::vector<double> ((n + 1) * (m + 1));
	
	for(int i = 0; i <= n; ++i){

		int start = Get_single_index(i, 0, m + 1);

		// interval == 1 since we are in teh y direction
		// restriction: children inderit parent's polynomial order. 
		Interpolate_to_new_points(m + 1,  m + 1, T_y,
				local::Hash_elem[p_key] -> holdmetrics.x_node, middle, start, start, 1);
	}
	

	// x direction
	for(int j = 0; j <= m; ++j){

		int start = Get_single_index(0, j, m + 1);

		Interpolate_to_new_points(n + 1,  n + 1, T_x,
				middle, c -> holdmetrics.x_node, start, start, m + 1);
		
	}

	middle.clear();

	//y_node:
	// y direction
	middle = std::vector<double> ((n + 1) * (m + 1));
	
	for(int i = 0; i <= n; ++i){

		int start = Get_single_index(i, 0, m + 1);

		// interval == 1 since we are in teh y direction
		// restriction: children inderit parent's polynomial order. 
		Interpolate_to_new_points(m + 1,  m + 1, T_y,
				local::Hash_elem[p_key] -> holdmetrics.y_node, middle, start, start, 1);
	}
	

	// x direction
	for(int j = 0; j <= m; ++j){

		int start = Get_single_index(0, j, m + 1);

		Interpolate_to_new_points(n + 1,  n + 1, T_x,
				middle, c -> holdmetrics.y_node, start, start, m + 1);
		
	}

	middle.clear();

	//Calculate Boundaries
	boundaryinterp1D(c, n, m);

	std::vector<double> one_x(n+1);
	std::vector<double> one_y(n+1);
	std::vector<double> two_x(n+1);
	std::vector<double> two_y(n+1);
	std::vector<double> three_x(n+1);
	std::vector<double> three_y(n+1);
	std::vector<double> four_x(n+1);
	std::vector<double> four_y(n+1);


	for(int i = 0; i <= n; i++){
		one_x[i] = c->holdmetrics.x_boundary[Get_single_index(i, 0, 4)];
		one_y[i] = c->holdmetrics.y_boundary[Get_single_index(i, 0, 4)]; 

		two_x[i] = c->holdmetrics.x_boundary[Get_single_index(i, 1, 4)];
		two_y[i] = c->holdmetrics.y_boundary[Get_single_index(i, 1, 4)];

		three_x[i] = c->holdmetrics.x_boundary[Get_single_index(i, 2, 4)]; 
		three_y[i] = c->holdmetrics.y_boundary[Get_single_index(i, 2, 4)];

		four_x[i] = c->holdmetrics.x_boundary[Get_single_index(i, 3, 4)]; 
		four_y[i] = c->holdmetrics.y_boundary[Get_single_index(i, 3, 4)];
	}

	CurveInterpolant bot(n, nodal::gl_points[n], one_x, one_y);
	CurveInterpolant right(n, nodal::gl_points[n], two_x, two_y);
	CurveInterpolant top(n, nodal::gl_points[n], three_x, three_y);
	CurveInterpolant left(n, nodal::gl_points[n], four_x, four_y);

	CurveInterpolant refine_B[4] = {bot, right, top, left};

	//metrics, jacobians, normals, scaling factors
	c->holdmetrics.recalcscalenorms(refine_B, nodal::gl_points[n]);





	// //delx_delxi:
	// // y direction
	// middle = std::vector<double> ((n + 1) * (m + 1));
	
	// for(int i = 0; i <= n; ++i){

	// 	int start = Get_single_index(i, 0, m + 1);

	// 	// interval == 1 since we are in teh y direction
	// 	// restriction: children inderit parent's polynomial order. 
	// 	Interpolate_to_new_points(m + 1,  m + 1, T_y,
	// 			local::Hash_elem[p_key] -> holdmetrics.delx_delxi, middle, start, start, 1);
	// }
	

	// // x direction
	// for(int j = 0; j <= m; ++j){

	// 	int start = Get_single_index(0, j, m + 1);

	// 	Interpolate_to_new_points(n + 1,  n + 1, T_x,
	// 			middle, c -> holdmetrics.delx_delxi, start, start, m + 1);
		
	// }

	// middle.clear();

	// //delx_deleta:
	// // y direction
	// middle = std::vector<double> ((n + 1) * (m + 1));
	
	// for(int i = 0; i <= n; ++i){

	// 	int start = Get_single_index(i, 0, m + 1);

	// 	// interval == 1 since we are in teh y direction
	// 	// restriction: children inderit parent's polynomial order. 
	// 	Interpolate_to_new_points(m + 1,  m + 1, T_y,
	// 			local::Hash_elem[p_key] -> holdmetrics.delx_deleta, middle, start, start, 1);
	// }
	

	// // x direction
	// for(int j = 0; j <= m; ++j){

	// 	int start = Get_single_index(0, j, m + 1);

	// 	Interpolate_to_new_points(n + 1,  n + 1, T_x,
	// 			middle, c -> holdmetrics.delx_deleta, start, start, m + 1);
		
	// }

	// middle.clear();


	// //dely_delxi:
	// // y direction
	// middle = std::vector<double> ((n + 1) * (m + 1));
	
	// for(int i = 0; i <= n; ++i){

	// 	int start = Get_single_index(i, 0, m + 1);

	// 	// interval == 1 since we are in teh y direction
	// 	// restriction: children inderit parent's polynomial order. 
	// 	Interpolate_to_new_points(m + 1,  m + 1, T_y,
	// 			local::Hash_elem[p_key] -> holdmetrics.dely_delxi, middle, start, start, 1);
	// }
	

	// // x direction
	// for(int j = 0; j <= m; ++j){

	// 	int start = Get_single_index(0, j, m + 1);

	// 	Interpolate_to_new_points(n + 1,  n + 1, T_x,
	// 			middle, c -> holdmetrics.dely_delxi, start, start, m + 1);
		
	// }

	// middle.clear();


	// //dely_deleta:
	// // y direction
	// middle = std::vector<double> ((n + 1) * (m + 1));
	
	// for(int i = 0; i <= n; ++i){

	// 	int start = Get_single_index(i, 0, m + 1);

	// 	// interval == 1 since we are in teh y direction
	// 	// restriction: children inderit parent's polynomial order. 
	// 	Interpolate_to_new_points(m + 1,  m + 1, T_y,
	// 			local::Hash_elem[p_key] -> holdmetrics.dely_deleta, middle, start, start, 1);
	// }
	

	// // x direction
	// for(int j = 0; j <= m; ++j){

	// 	int start = Get_single_index(0, j, m + 1);

	// 	Interpolate_to_new_points(n + 1,  n + 1, T_x,
	// 			middle, c -> holdmetrics.dely_deleta, start, start, m + 1);
		
	// }

	// middle.clear();


	// //jacobian:
	// // y direction
	// middle = std::vector<double> ((n + 1) * (m + 1));
	
	// for(int i = 0; i <= n; ++i){

	// 	int start = Get_single_index(i, 0, m + 1);

	// 	// interval == 1 since we are in teh y direction
	// 	// restriction: children inderit parent's polynomial order. 
	// 	Interpolate_to_new_points(m + 1,  m + 1, T_y,
	// 			local::Hash_elem[p_key] -> holdmetrics.jacobian, middle, start, start, 1);
	// }
	

	// // x direction
	// for(int j = 0; j <= m; ++j){

	// 	int start = Get_single_index(0, j, m + 1);

	// 	Interpolate_to_new_points(n + 1,  n + 1, T_x,
	// 			middle, c -> holdmetrics.jacobian, start, start, m + 1);
		
	// }

	// middle.clear();

	// for(int j = 0; j <= m; j++){
	// 	for(int i = 0; i <= n; i++){
	// 		int num_p = Get_single_index(i, j, n + 1);
	// 		c->holdmetrics.delx_delxi[num_p] *= 0.5;
	// 		c->holdmetrics.dely_delxi[num_p] *= 0.5;
	// 		c->holdmetrics.delx_deleta[num_p] *= 0.5;
	// 		c->holdmetrics.dely_deleta[num_p] *= 0.5;
	// 		c->holdmetrics.jacobian[num_p] *= 0.25;
	// 	}
	// }

}


/// @brief
/// Map the new set of points from the reference space [-1, 1] to the space after h-refinement. 
/// @param m polynomial order.
/// @param start element left boundary coordinate.
/// @param y new set of points' coordinate. 
void Form_new_set_of_points(int m, int start, std::vector<double>& y){

	for(int j = 0; j <= m; ++j){
		
		y[j] = Affine_mapping(nodal::gl_points[m][j], start, 1.0);	// element size = 1.0
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
			}

			for(int j = 0; j < n; ++j){

				int index = Get_single_index(k, j, n);

				T[index] /= s;
				

			}
			

		}
	}


}

/// @brief
/// Interpolation between two sets of points by matrix multiplicaiton
/// @param m The number of new set of points.
/// @param n The number of old set of points.	
/// @param T Interpolating matrix. Size m * n.
/// @param f Solutions on the old set of points.
/// @param new_f Solusion on the new set of points. Size m. 
/// @param start_old The node number of the first value in f.
/// @param start_new The node number of the first value in new_f.
/// @param interval node index interval.
void Interpolate_to_new_points(int m, int n, std::vector<double>& T, 
				std::vector<double>& f, std::vector<double>& new_f, int start_old, int start_new, int interval){

	int ni = start_new;	
	for(int i = 0; i < m; ++i){

		double t{};

		int fi = start_old;

		for(int j = 0; j < n; ++j){

			int index = Get_single_index(i, j, n);

			t += T[index] * f[fi];
			fi += interval;
		}
		
		new_f[ni] = t;

		ni += interval;

	}

}

void coarsetofine(std::vector<double> &target, int new_order, int currentorder, std::vector<double> &T){
	int start_old = 0;
	int start_new = 0;

	std::vector<double> intermatrix =  std::vector<double>((new_order + 1) * (currentorder + 1));
	
	for(int i = 0; i <= (currentorder); ++i){
		Interpolate_to_new_points(new_order + 1, currentorder + 1, T,
				target, intermatrix, start_old, start_new, 1);

		start_old += (currentorder + 1);
		start_new += (new_order + 1);
	}

	target.clear();

	start_old = 0;
	start_new = 0;
	
	// allocate space
	target = std::vector<double> ((new_order + 1) * (new_order + 1));

	for(int j = 0; j <= new_order; ++j){
		Interpolate_to_new_points(new_order + 1,  currentorder + 1, T,
				intermatrix, target, start_old, start_new, new_order + 1);
		
		++start_old;
		++start_new;
	}
}

//m is new order, n is current order
//works
void coarsetofine1D(std::vector<double> &target, int new_order, int currentorder, const std::vector<double> &T){

	//intermediate vector
	std::vector<double> intermediate = std::vector<double> (target.size());

	//copy over orignal data to intermediate
	for(int i = 0; i < target.size(); i++){
		intermediate[i] = target[i];
	}

	//reallocate space
	target.clear();
	target = std::vector<double> ((new_order + 1)*4);

	for(int k = 0; k < 4; k++){

		for(int i = 0; i <= new_order; i++){
			double t = 0;
			int ndx = Get_single_index(i, k, 4);

			for(int j = 0; j <= currentorder; j++){
				int indx = Get_single_index(i, j, currentorder + 1);

				int odx = Get_single_index(j, k, 4);
				
				t += T[indx]*intermediate[odx];
				// std::cout << intermediate[odx];
			}

			target[ndx] = t;
		}

	}
	
	
}


//dont allocate for fnew
//poly order of nodes, number of new poitns-1, old x, old y, old solution, new x, new y, new soln
void CoursetoFineInterp(int nold, int nnew, std::vector<double> &x, std::vector<double> &y, std::unordered_map<int, std::vector<double>> &f, 
 std::vector<double> &xnew,  std::vector<double> &ynew, std::unordered_map<int, std::vector<double>> &fnew){

	std::vector<double> T1;

	Polynomial_interpolate_matrix(x, xnew, T1);

	std::unordered_map<int, std::vector<double>> middle;	// intermidiate matrix

	// y direction
	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){
		
		int start_old{};
		int start_new{};

		middle[equ] = std::vector<double> ((nnew + 1) * (nold + 1));
		
		for(int i = 0; i <= nold; ++i){

			// interval == 1 since we are in teh y direction
			// restriction: children inderit parent's polynomial order. 
			Interpolate_to_new_points(nnew+ 1, nold + 1, T1,
					f[equ], middle[equ], start_old, start_new, 1);

			start_old += (nold + 1);
			start_new += (nnew + 1);
			// ++start_old;
			// ++start_new;
		}

	}

	std::vector<double> T2;

	Polynomial_interpolate_matrix(y, ynew, T2);

	// x direction
	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

		int start_old{};
		int start_new{};
		
		// allocate space
		fnew[equ] = std::vector<double> ((nnew + 1) * (nnew + 1));

		for(int j = 0; j <= nnew; ++j){

			Interpolate_to_new_points(nnew + 1, nold + 1, T2,
					middle[equ], fnew[equ], start_old, start_new, nnew + 1);
			
			++start_old;
			++start_new;
			// start_old += (nold + 1);
			// start_new += (nnew + 1);
		}

	}
	
		
}


//generate a linearly spaced array (vector)
std::vector<double> linarray(double min, double max, int N){
    std::vector<double> result;

    double spacing = (max-min)/double(N-1);

    for(int i = 0; i < N; i++) {
        result.push_back(min + i * spacing);
    }

    return result;
}


//m is new order, n is current order
//works
void boundaryinterp1D(Unit* c, int n, int m){

	//X Boundaries
	for(int j = 0; j <= m; ++j ){

			std::vector<double> x_array(n + 1);
			std::vector<double> y_array(n + 1);

			for(int i = 0; i <= n; ++i){

				int nodei = Get_single_index(i, j, m + 1);

				x_array[i] = c->holdmetrics.x_node[nodei];
				y_array[i] = c->holdmetrics.y_node[nodei];
			}

			int boundaryindex1 = Get_single_index(j, 1, 4);
			int boundaryindex3 = Get_single_index(j, 3, 4);

			double xl = Interpolate_to_boundary(n, x_array, nodal::lagrange_l[n]);
			double xr = Interpolate_to_boundary(n, x_array, nodal::lagrange_r[n]);

			c->holdmetrics.x_boundary[boundaryindex3] = xl;
			c->holdmetrics.x_boundary[boundaryindex1] = xr;

			

			double yl = Interpolate_to_boundary(n, y_array, nodal::lagrange_l[n]);
			double yr = Interpolate_to_boundary(n, y_array, nodal::lagrange_r[n]);

			c->holdmetrics.y_boundary[boundaryindex3] = yl;
			c->holdmetrics.y_boundary[boundaryindex1] = yr;



	}
	

	//Y Boundaires
	for(int i = 0; i <= n; ++i ){

			std::vector<double> x_array(m + 1);
			std::vector<double> y_array(m + 1);

			for(int j = 0; j <= m; ++j){

				int nodei = Get_single_index(i, j, m + 1);

				x_array[j] = c->holdmetrics.x_node[nodei];
				y_array[j] = c->holdmetrics.y_node[nodei];

			}

			int boundaryindex0 = Get_single_index(i, 0, 4);
			int boundaryindex2 = Get_single_index(i, 2, 4);

			//check interpolate to boundary outputs
			c->holdmetrics.x_boundary[boundaryindex0] = Interpolate_to_boundary(n, x_array, nodal::lagrange_l[n]);
			c->holdmetrics.x_boundary[boundaryindex2] = Interpolate_to_boundary(n, x_array, nodal::lagrange_r[n]);

			c->holdmetrics.y_boundary[boundaryindex0] = Interpolate_to_boundary(n, y_array, nodal::lagrange_l[n]);
			c->holdmetrics.y_boundary[boundaryindex2] = Interpolate_to_boundary(n, y_array, nodal::lagrange_r[n]);	
	}

}

//dont allocate for fnew
//poly order of nodes, number of new poitns-1, old x, old y, old solution, new x, new y, new soln
void CoursetoFineInterpCoords(int nold, int nnew, std::vector<double> &x, std::vector<double> &y, std::unordered_map<int, std::vector<double>> &f, 
 std::vector<double> &xnew,  std::vector<double> &ynew, std::unordered_map<int, std::vector<double>> &fnew){

	std::vector<double> T1;

	Polynomial_interpolate_matrix(x, xnew, T1);

	std::unordered_map<int, std::vector<double>> middle;	// intermidiate matrix

	// y direction
	for(int equ = 0; equ < 2; ++equ){
		
		int start_old{};
		int start_new{};

		middle[equ] = std::vector<double> ((nnew + 1) * (nold + 1));
		
		for(int i = 0; i <= nold; ++i){

			// interval == 1 since we are in teh y direction
			// restriction: children inderit parent's polynomial order. 
			Interpolate_to_new_points(nnew+ 1, nold + 1, T1,
					f[equ], middle[equ], start_old, start_new, 1);

			start_old += (nold + 1);
			start_new += (nnew + 1);
			// ++start_old;
			// ++start_new;
		}

	}

	std::vector<double> T2;

	Polynomial_interpolate_matrix(y, ynew, T2);

	// x direction
	for(int equ = 0; equ < 2; ++equ){

		int start_old{};
		int start_new{};
		
		// allocate space
		fnew[equ] = std::vector<double> ((nnew + 1) * (nnew + 1));

		for(int j = 0; j <= nnew; ++j){

			Interpolate_to_new_points(nnew + 1, nold + 1, T2,
					middle[equ], fnew[equ], start_old, start_new, nnew + 1);
			
			++start_old;
			++start_new;
			// start_old += (nold + 1);
			// start_new += (nnew + 1);
		}

	}
	
		
}

