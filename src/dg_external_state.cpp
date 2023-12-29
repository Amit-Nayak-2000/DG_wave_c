#include "dg_external_state.h"
#include <vector>
#include "dg_user_defined.h"
#include <cmath>
#include "dg_param.h"


// forward declaration -------------------------------------------------------------------------------------------------
void External_state_Gaussian_exact(double t, double x, double y, std::vector<double>& q_ext, std::vector<int>& index);

void External_state_reflect_x(std::vector<double>& q_int, std::vector<double>& q_ext, 
				std::vector<int>& index);

void External_state_reflect_y(std::vector<double>& q_int, std::vector<double>& q_ext, 
				std::vector<int>& index);

void External_mirror_right_boundary(double t, double x, double y, std::vector<double>& q_ext, std::vector<int>& index);

void External_state_sin_exact(double t, double x, double y, std::vector<double>& q_ext, std::vector<int>& index);
//---------------------------------------------------------------------------------------------------------------------

/// @brief
/// Boundary conditions. User defined.
/// @param t current time.
/// @param x x coordinate (physical). 
/// @param y y coordinate (physical). 
/// @param q_ext exrernal solutions (outputs). 
/// @param index the index of three variables. 
void External_state_Gaussian_exact(double t, double x, double y, std::vector<double>& q_ext, std::vector<int>& index){


	double inter = exp( - pow((user::kx * ( x - user::xx0) + 
			user::ky * (y - user::yy0) - dg_fun::C * t), 2) / (user::D * user::D) );

	q_ext[index[0]] = inter;

	q_ext[index[1]] = user::kx / dg_fun::C * inter;

	q_ext[index[2]] = user::ky / dg_fun::C * inter;
}

// void External_state_Gaussian_exact(double t, double x, double y, std::vector<double>& q_ext, std::vector<int>& index){


// 	double inter = exp( - pow((user::kx * ( x - user::xx0) + 
// 			user::ky * (y - user::yy0) - dg_fun::C * t), 2) / (user::D * user::D) );

// 	double inter2 = -1 * exp( - pow((-user::kx * ( x + user::xx0) + 
// 			-user::ky * (y + user::yy0) - dg_fun::C * t), 2) / (user::D * user::D) );

// 	q_ext[index[0]] = inter + inter2;

// 	q_ext[index[1]] = (user::kx / dg_fun::C * inter) + (-user::kx / dg_fun::C * inter2);

// 	q_ext[index[2]] = (user::ky / dg_fun::C * inter) + (-user::ky / dg_fun::C * inter2);
// }

/// @brief
/// External state applied as a mirror image wave starts from the right boundary. 
/// Note that the domain should be [0, 1] * [0, 1]. 
/// The mirror image comes from [1, 2] * [0, 1]
/// @param t current time.
/// @param x x coordinate (physical). 
/// @param y y coordinate (physical). 
/// @param q_ext exrernal solutions (outputs). 
/// @param index the index of three variables. 
void External_mirror_right_boundary(double t, double x, double y, std::vector<double>& q_ext, std::vector<int>& index){


	double inter = exp( - pow((mirror::kx * (x - mirror::x0) + 
			mirror::ky * (y - mirror::y0) - dg_fun::C * t), 2) / (mirror::D * mirror::D) );

	q_ext[index[0]] += inter;

	q_ext[index[1]] += mirror::kx / dg_fun::C * inter;

	q_ext[index[2]] += mirror::ky / dg_fun::C * inter;

}



/// @brief
/// Reflect boundary conditions. 
void External_state_reflect_x(std::vector<double>& q_int, std::vector<double>& q_ext, 
				std::vector<int>& index){

	q_ext[index[0]] = q_int[index[0]];
	
	q_ext[index[1]] = (2.0 * q_int[index[0]] + dg_fun::C * q_int[index[1]]) / dg_fun::C;

	q_ext[index[2]] = q_int[index[2]];

}

void External_state_reflect_y(std::vector<double>& q_int, std::vector<double>& q_ext, 
				std::vector<int>& index){

	q_ext[index[0]] = q_int[index[0]];

	q_ext[index[1]] = 0.0;

	q_ext[index[2]] = (2.0 * q_int[index[0]] + dg_fun::C * q_int[index[2]]) / dg_fun::C;


}


// test
void External_state_sin_exact(double t, double x, double y, std::vector<double>& q_ext, std::vector<int>& index){

	double inter1 = std::sin(user::pi * x + t);
	double inter2 = std::sin(user::pi * y + t);

	q_ext[index[0]] = -1.0 / user::pi * inter1 - 1.0 / user::pi * inter2;

	q_ext[index[1]] = inter1;

	q_ext[index[2]] = inter2;
}

//Amit's BCs

//alpha is normalx, beta is normaly
void External_state_reflect(std::vector<double>& q_int, std::vector<double>& q_ext, 
				std::vector<int>& index, double alpha, double beta){

	//Pressure ext = Pressure int.
	q_ext[index[0]] = q_int[index[0]];
	
	//u ext = (beta^2 - alpha^2)*uint - 2*alpha*beta*vint
	q_ext[index[1]] = (beta*beta - alpha*alpha)*q_int[index[1]] - 2*alpha*beta*q_int[index[2]];

	//v ext = -2*alpha*beta*uint + (alpha^2 - beta^2)*vint
	q_ext[index[2]] = -2*alpha*beta*q_int[index[1]] + (alpha*alpha - beta*beta)*q_int[index[2]];

}

void External_state_radiation(std::vector<double>& q_int, std::vector<double>& q_ext, 
				std::vector<int>& index){

	q_ext[index[0]] = 0;
	
	q_ext[index[1]] = 0;

	q_ext[index[2]] = 0;

}
