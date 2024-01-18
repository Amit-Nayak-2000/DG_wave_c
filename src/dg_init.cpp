#include "dg_local_storage.h"
#include "dg_unit.h"
#include "dg_param.h"
#include "dg_affine_map.h"
#include "dg_nodal_2d_storage.h"
#include <cmath>	// exp
#include "dg_single_index.h"
#include "dg_user_defined.h"
#include <iostream>	// test
#include "dg_transfinite_quad_map.h"

// forward declaration--------------------------------
void DG_init();

void DG_init_new();
//----------------------------------------------------

CurveInterpolant B[4];

/// @brief
/// Initialization all local elements based on the initial conditions.
void DG_init(){
	
	Unit* temp = local::head;

	//Building global bounds
	double PI = 4.0 * atan(1.0); 
	double rinf = 5.0;
	double ro = 0.5; //originally 0.5
	double w = 0.125; // originally 0.125
	double xs = 1.5; //originally 1.5
	double ys = 0.0; //originally 0.0

	std::vector<double> one_x(grid::nmin+1);
	std::vector<double> one_y(grid::nmin+1);
	std::vector<double> two_x(grid::nmin+1);
	std::vector<double> two_y(grid::nmin+1);
	std::vector<double> three_x(grid::nmin+1);
	std::vector<double> three_y(grid::nmin+1);
	std::vector<double> four_x(grid::nmin+1);
	std::vector<double> four_y(grid::nmin+1);

	for(int p = 0; p <= grid::nmin; p++){
		// // reference square 
		// one_x[p] = 1*nodal::gl_points[grid::nmin][p];
		// one_y[p] = -1; 

		// two_x[p] = 1; 
		// two_y[p] = 1*nodal::gl_points[grid::nmin][p];

		// three_x[p] = 1*nodal::gl_points[grid::nmin][p];
		// three_y[p] = 1; 

		// four_x[p] = -1; 
		// four_y[p] = 1*nodal::gl_points[grid::nmin][p];


		//bumpy lower boundary
		// one_x[p] = 1*nodal::gl_points[grid::nmin][p];
		// one_y[p] = (-0.25 * (std::sin(36*PI*nodal::gl_points[grid::nmin][p]) * std::sin(36*PI*nodal::gl_points[grid::nmin][p]))) - 1;

		// two_x[p] = 1; 
		// two_y[p] = 1*nodal::gl_points[grid::nmin][p];

		// three_x[p] = 1*nodal::gl_points[grid::nmin][p];
		// three_y[p] = 1; 

		// four_x[p] = -1; 
		// four_y[p] = 1*nodal::gl_points[grid::nmin][p];


		// half annulus
		// one_x[p] = ro + (rinf - ro)*(nodal::gl_points[grid::nmin][p]+1)/2;
		// one_y[p] = 0; 

		// two_x[p] = rinf*std::cos(PI*(nodal::gl_points[grid::nmin][p]+1)/2); 
		// two_y[p] = rinf*std::sin(PI*(nodal::gl_points[grid::nmin][p]+1)/2); 

		// three_x[p] = -ro - (rinf - ro)*(nodal::gl_points[grid::nmin][p]+1)/2; 
		// three_y[p] = 0;

		// four_x[p] = ro*std::cos(PI*(nodal::gl_points[grid::nmin][p]+1)/2); 
		// four_y[p] = ro*std::sin(PI*(nodal::gl_points[grid::nmin][p]+1)/2);

		

		// Curved Channel mine
		// one_x[p] = 3*nodal::gl_points[grid::nmin][p]/2;
		// one_y[p] = -(0.3 + 0.35*(std::tanh(2*nodal::gl_points[grid::nmin][p]) + 1 )); 

		// two_x[p] = 3/2; 
		// two_y[p] = (0.289 + 0.35*(std::tanh(3) + 1 ))*nodal::gl_points[grid::nmin][p]; 

		// three_x[p] = 3*nodal::gl_points[grid::nmin][p]/2;
		// three_y[p] = (0.3 + 0.35*(std::tanh(2*nodal::gl_points[grid::nmin][p]) + 1 ));

		// four_x[p] = -3/2; 
		// four_y[p] = (0.155*(std::tanh(3) + 1 ))*nodal::gl_points[grid::nmin][p];

		// Rectangle
		one_x[p] = 5*nodal::gl_points[grid::nmin][p];
		one_y[p] = 0; 

		two_x[p] = 5; 
		two_y[p] = 5*(nodal::gl_points[grid::nmin][p] + 1)/2;

		three_x[p] = 5*nodal::gl_points[grid::nmin][p];
		three_y[p] = 5; 

		four_x[p] = -5; 
		four_y[p] = 5*(nodal::gl_points[grid::nmin][p] + 1)/2;


	}


	CurveInterpolant bot(grid::nmin, nodal::gl_points[grid::nmin], one_x, one_y);
	CurveInterpolant right(grid::nmin, nodal::gl_points[grid::nmin], two_x, two_y);
	CurveInterpolant top(grid::nmin, nodal::gl_points[grid::nmin], three_x, three_y);
	CurveInterpolant left(grid::nmin, nodal::gl_points[grid::nmin], four_x, four_y);

	// CurveInterpolant B[4] = {bot, right, top, left};

	B[0] = bot;
	B[1] = right;
	B[2] = top;
	B[3] = left;

	//Finished building global bounds

	// traverse the linked list
	for(int k = 0; k < local::local_elem_num; ++k){

		// elemement size
		double del_x = (temp -> xcoords[1]) - (temp -> xcoords[0]);
		double del_y = (temp -> ycoords[1]) - (temp -> ycoords[0]);


		//setting local curvilinear boundaries
		std::vector<double> local_1x, local_1y, local_2x, local_2y, local_3x, local_3y, local_4x, local_4y;
		local_1x = std::vector<double>(grid::nmin + 1);
		local_1y = std::vector<double>(grid::nmin + 1);
		local_2x = std::vector<double>(grid::nmin + 1);
		local_2y = std::vector<double>(grid::nmin + 1);
		local_3x = std::vector<double>(grid::nmin + 1);
		local_3y = std::vector<double>(grid::nmin + 1);
		local_4x = std::vector<double>(grid::nmin + 1);
		local_4y = std::vector<double>(grid::nmin + 1);

		double xmin = Affine_mapping(nodal::gl_points[grid::nmin][0], temp -> xcoords[0], del_x);
		double xmax = Affine_mapping(nodal::gl_points[grid::nmin][grid::nmin], temp -> xcoords[0], del_x);
		double ymin = Affine_mapping(nodal::gl_points[grid::nmin][0], temp -> ycoords[0], del_y);
		double ymax = Affine_mapping(nodal::gl_points[grid::nmin][grid::nmin], temp -> ycoords[0], del_y);

		
		xmin = temp->xcoords[0];
		ymin = temp->ycoords[0];

		xmax = temp->xcoords[1];
		ymax = temp->ycoords[1];


		xmin *= 2;
		xmin -= 1;

		ymin *= 2;
		ymin -= 1;

		xmax *= 2;
		xmax -= 1;

		ymax *= 2;
		ymax-= 1;


		for(int p = 0; p <= grid::nmin; p++){
			double xlocal = Affine_mapping(nodal::gl_points[grid::nmin][p], temp -> xcoords[0], del_x);
			double ylocal = Affine_mapping(nodal::gl_points[grid::nmin][p], temp -> ycoords[0], del_y);

			//bring back to reference square
			xlocal *= 2;
			xlocal -= 1;

			ylocal *= 2;
			ylocal -= 1;

			// build local boundaries
			TransfiniteQuadMap(xlocal, ymin, B, local_1x[p], local_1y[p]);
			TransfiniteQuadMap(xmax, ylocal, B, local_2x[p], local_2y[p]);
			TransfiniteQuadMap(xlocal, ymax, B, local_3x[p], local_3y[p]);
			TransfiniteQuadMap(xmin, ylocal, B, local_4x[p], local_4y[p]);

		}


		CurveInterpolant local_bot(grid::nmin, nodal::gl_points[grid::nmin], local_1x, local_1y);
		CurveInterpolant local_right(grid::nmin, nodal::gl_points[grid::nmin], local_2x, local_2y);
		CurveInterpolant local_top(grid::nmin, nodal::gl_points[grid::nmin], local_3x, local_3y);
		CurveInterpolant local_left(grid::nmin, nodal::gl_points[grid::nmin], local_4x, local_4y);

		CurveInterpolant local_B[4] = {local_bot, local_right, local_top, local_left};
		

		//initialize metrics for the unit
		temp -> holdmetrics.initMetrics(local_B, nodal::gl_points[grid::nmin]);

		double x, y;

		//setting initial conditions
		for(int j = 0; j <= grid::nmin; ++j){
			
			for(int i = 0; i <= grid::nmin; ++i){
			
				int num_p = Get_single_index(i, j, grid::nmin + 1);

				x = temp->holdmetrics.x_node[num_p];
				y = temp->holdmetrics.y_node[num_p];


				// wave --------------------------------------------------------------------------------	
				// double inter = exp( - std::pow((user::kx * (x - user::xx0) + 
				// 			user::ky * (y - user::yy0)), 2) / std::pow(user::D, 2));
				
				// temp -> solution[0][num_p] = inter;
				// temp -> solution[1][num_p] = (user::kx / dg_fun::C * inter);
				// temp -> solution[2][num_p] = (user::ky / dg_fun::C * inter);

				// double inter2 = -1 * exp( - std::pow((-user::kx * (x + user::xx0) + 
				// 			-user::ky * (y + user::yy0)), 2) / std::pow(user::D, 2));
      			
				// temp -> solution[0][num_p] = inter + inter2;
				// temp -> solution[1][num_p] = (user::kx / dg_fun::C * inter) + (-user::kx / dg_fun::C * inter2);
				// temp -> solution[2][num_p] = (user::ky / dg_fun::C * inter) + (-user::ky / dg_fun::C * inter2);
				// // -------------------------------------------------------------------------------------
				

				// test case acoustic scattering --------------------------------------------------
				// p = exp( -ln(2)*( ((x-xs)^2 + y^2)/w^2 )
				temp -> solution[0][num_p] = exp( -log(2) * ( (std::pow((temp->holdmetrics.x_node[num_p] - xs), 2) + std::pow((temp->holdmetrics.y_node[num_p] - ys), 2) ) / std::pow(w, 2) ) );
				// std::cout << temp ->solution[0][num_p] << std::endl;
				//u = 0
				temp -> solution[1][num_p] = 0;
				//v = 0
				temp -> solution[2][num_p] = 0;
				//----------------------------------------------------------------------------------------

				// test case planar acoustic scattering --------------------------------------------------
				// p = exp( -ln(2)*( ((x-xs)^2 + y^2)/w^2 )
				// temp -> solution[0][num_p] = exp( -log(2) * ( (std::pow((temp->holdmetrics.x_node[num_p] - xs), 2) ) / std::pow(w, 2) ) );
				// // u = p;
				// temp -> solution[1][num_p] = exp( -log(2) * ( (std::pow((temp->holdmetrics.x_node[num_p] - xs), 2) ) / std::pow(w, 2) ) );
				// // v = 0;
				// temp -> solution[2][num_p] = 0;
				//----------------------------------------------------------------------------------------

			}
		}
		
		// std::cout << "Original" << std::endl;
		// temp->holdmetrics.outputMetrics();
	
		temp = temp -> next;		

	}


}

void DG_init2(){
	
	Unit* temp = local::head;

	// traverse the linked list
	for(int k = 0; k < local::local_elem_num; ++k){

		// elemement size
		double del_x = (temp -> xcoords[1]) - (temp -> xcoords[0]);
		double del_y = (temp -> ycoords[1]) - (temp -> ycoords[0]);

		for(int j = 0; j <= grid::nmin; ++j){
			
			// map reference location to physical localtion
			double gl_p_y = nodal::gl_points[grid::nmin][j];
			double y = Affine_mapping(gl_p_y, temp -> ycoords[0], del_y);
			
			for(int i = 0; i <= grid::nmin; ++i){
	
				double gl_p_x = nodal::gl_points[grid::nmin][i];
				double x = Affine_mapping(gl_p_x, temp -> xcoords[0], del_x);
			
				int num_p = Get_single_index(i, j, grid::nmin + 1);

				// wave --------------------------------------------------------------------------------	
				double kx = 0;	// sin(0) = 0
				double ky = 1.0;	// cos(0) = 1
				double inter = exp( - std::pow((kx * (x - user::xx0) + 
							ky * (y - user::yy0)), 2) / std::pow(user::D, 2));
      			
				temp -> solution[0][num_p] = inter;
				temp -> solution[1][num_p] = kx / dg_fun::C * inter;
				temp -> solution[2][num_p] = ky / dg_fun::C * inter;
				// -------------------------------------------------------------------------------------
				

			}
		}
	
		temp = temp -> next;		

	}


}

/// @brief
/// Initialize the domain. Not using the exact solutions. Uniformly 0. 
void DG_init_new(){
	
	Unit* temp = local::head;

	// traverse the linked list
	for(int k = 0; k < local::local_elem_num; ++k){

		for(int j = 0; j <= grid::nmin; ++j){
			
			
			for(int i = 0; i <= grid::nmin; ++i){
			
				int num_p = Get_single_index(i, j, grid::nmin + 1);

				temp -> solution[0][num_p] = 0.1;
				temp -> solution[1][num_p] = 0.1;
				temp -> solution[2][num_p] = 0.1;

			}
		}
	
		temp = temp -> next;		

	}


}
