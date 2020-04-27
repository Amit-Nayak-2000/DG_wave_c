#include "dg_error_estimator.h"
#include "dg_unit.h"
#include <vector>
#include <cassert>
#include "dg_param.h"
#include "dg_single_index.h"
#include <unordered_map>
#include "dg_basis.h"
#include "dg_nodal_2d_storage.h"
#include <cmath>
#include "dg_local_storage.h"

// forward declaration---------------------------------------------------------
double Decay_rate(std::vector<int>& porder, std::vector<double>& ap);

void Error_indicator(Unit* temp, std::vector<double>& sigma, std::vector<bool>& flag);

double Solution_l2_norm(int equ, Unit* temp);

void Refinement_flag();
//-----------------------------------------------------------------------------

void Refinement_flag(){

	Unit* temp = local::head;

	for(int k = 0; k < local::local_elem_num; ++k){

		std::vector<double> sigma(dg_fun::num_of_equation);
		std::vector<bool> flag(dg_fun::num_of_equation);

		Error_indicator(temp, sigma, flag);
		
		// now refine depends on pressure--------------
		if(flag.front()){
			if(sigma.front() < 1){	// h-refinemnt
				
				temp -> hrefine = true;
	
			}
			else{
	
				temp -> prefine = true;
			}
		}
		//-------------------------------------------------
		temp = temp -> next;
	}

}

/// @brief
/// Calculate the error indicator sigma, if sigma < 1 apply h-refinemnt, 
/// if sigma >= 1 apply prefinement. 
/// @param temp pointer to the estimated element.
/// @note The polynomial order in x and y direction should be identical.
void Error_indicator(Unit* temp, std::vector<double>& sigma, std::vector<bool>& flag){

	assert((temp -> n > 6) && "Polynomial order is too low to estimate error.");

	std::unordered_map<int, std::vector<double>> ap;	// discrete spectrum ap
	for(int i = 0; i < dg_fun::num_of_equation; ++i){	// allocate space

		ap[i] = std::vector<double> (dg_refine::fit_point_num);
	}


	int N = temp -> n;
	int M = temp -> m;
	int p = N;	// assume N == M

	std::vector<int> porder(dg_refine::fit_point_num);	// record the poilynomial order

	for(int node = dg_refine::fit_point_num - 1; node >= 0; --node){

		porder[node] = p;

		// x direction-----------------------------------------------------
		for(int i = 0; i <= p; ++i ){	// compute sum of coefficients: sum(a(i, p))
			
			// a(n, m) = a(i, p), apply Gauss quadrature
			for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

				for(int k = 0; k <= N; ++k){	// x dir
	
					// x
					// compute Legendre polynomial of order i at GL point k
					double legendre_x{};
					double legendre_dir{};
					Legendre_polynomial_and_derivative(i, nodal::gl_points[N][k], 
										legendre_x, legendre_dir);
					for(int l = 0; l <= M; ++l){	// y dir
	
						int index = Get_single_index(k, l, M + 1);

						double legendre_y{};
						
						// y
						// compute Legendre polynomial of order l at GL point l
						Legendre_polynomial_and_derivative(p, nodal::gl_points[M][l], 
										legendre_y, legendre_dir);

						ap[equ][node] += (temp -> solution[equ][index]) * 
								legendre_x * legendre_y * 
								nodal::gl_weights[N][k] * nodal::gl_weights[M][l];
								
						
					}
					
				}

				ap[equ][node] *= (2.0 * (double)i + 1.0) * (2.0 * (double)p + 1.0) / 4.0;
			}

		}
		//-----------------------------------------------------------------

		// y direction-----------------------------------------------------
		for(int j = 0; j < p; ++j ){	// here does not need to include the last point (j < p)
			
					
			for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

				for(int k = 0; k <= N; ++k){	// x dir
	
					// x
					// compute Legendre polynomial of order i at GL point k
					double legendre_x{};
					double legendre_dir{};
					Legendre_polynomial_and_derivative(p, nodal::gl_points[N][k], 
										legendre_x, legendre_dir);
					for(int l = 0; l <= M; ++l){	// y dir
	
						int index = Get_single_index(k, l, M + 1);

						double legendre_y{};
						
						// y
						// compute Legendre polynomial of order l at GL point l
						Legendre_polynomial_and_derivative(j, nodal::gl_points[M][l], 
										legendre_y, legendre_dir);

						ap[equ][node] += (temp -> solution[equ][index]) * 
								legendre_x * legendre_y * 
								nodal::gl_weights[N][k] * nodal::gl_weights[M][l];
								
						
					}
					
				}

				ap[equ][node] *= (2.0 * (double)j + 1.0) * (2.0 * (double)p + 1.0) / 4.0;

			}

		}
		//-----------------------------------------------------------------
		p -= 2;
	}
		
	// refinment criteria: if exceed the acceptable level then flag as needed refinement. 
	for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){

		// L2 norm of the current element 
		double u_norm = Solution_l2_norm(equ, temp);

		double sum = ap[equ].back();	// sum of the last spectrum

		if(sum > u_norm){	// need refine

			flag[equ] = true;

			// get decay indicator
			sigma[equ] = Decay_rate(porder, ap[equ]);
		}
		else{

			flag[equ] = false;
		}


	}
	

}

/// @brief
/// Calculate the L2 norm of solution of the current element of the demanding equaiton.
/// @param equ the equation number
/// @param temp pointer to the Current element.
double Solution_l2_norm(int equ, Unit* temp){

	double norm{};

	int nodei{};

	for(int i = 0; i <= temp -> n; ++i){

		for(int j = 0; j <= temp -> m; ++j){

			norm += std::pow(temp -> solution[equ][nodei], 2);
			++nodei;
		}

	}


	return (std::sqrt(norm));

}



/// @brief
/// Assume the decaying spectrum ap is an exponential function ( ap  = const * exp ^ (-sigma * n)).
/// Apply log arithmetic to the two side of the function and fit the ap by a linear regression line.
/// Then we get the decay factor (the absolute value of the line). 
/// @param porder the corresponding order of ap. (x)
/// @param ap the decaying spectrums. (y)
double Decay_rate(std::vector<int>& porder, std::vector<double>& ap){

	double x_avg{}; double y_avg{};

	for(int i = 0; i < dg_refine::fit_point_num; ++i){

		x_avg += (double)porder[i];
		y_avg += ap[i];

	}
	x_avg /= (double)dg_refine::fit_point_num;
	y_avg /= (double)dg_refine::fit_point_num;

	
	double sigma{};
	double numer{};
	double denumer{};

	for(int i = 0; i < dg_refine::fit_point_num; ++i){
		
		numer += ((double)porder[i] - x_avg) * (ap[i] - y_avg);
		denumer += std::pow(((double)porder[i] - x_avg), 2);
	
	}
		
	sigma = std::abs(numer / denumer);


	return sigma;

}