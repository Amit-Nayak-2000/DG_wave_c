#ifndef DG_USER_DEFINED_H
#define DG_USER_DEFINED_H

#include <vector>
#include <unordered_map>

// user defined variable-------------------------------
namespace user{
	extern  double kx; 
	extern  double ky; 
	extern  double D;
	extern  double xx0;
	extern  double yy0; 
	extern int BC;
	
	extern const double pi;	// testing
};
//------------------------------------------------------

void Exact_solution_Gaussian(int n, int m, double x_l, double y_d,
				double del_x, double del_y, std::unordered_map<int, std::vector<double>>& e, double t);

void Exact_solution_Gaussian2(int n, int m, double x_l, double y_d,
				double del_x, double del_y, std::unordered_map<int, std::vector<double>>& e, double t);

void Exact_solution_sin(int n, int m, double x_l, double y_d,
				double del_x, double del_y, std::unordered_map<int, std::vector<double>>& e, double t);
#endif
