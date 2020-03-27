#include "dg_riemann_solver.h"
#include <vector>
#include <cassert>
#include "dg_param.h"

/// @brief
/// Riemann solver in x direction. Using upwind flux.
/// @param
void Riemann_solver_x(std::vector<double>& q_l, std::vector<double>& q_r, 
			std::vector<double>& n_flux, int normal, std::vector<int>& index){

	assert(normal == 1 || normal == -1 && "The normal unit vector should be 1 or -1");

	double p_l, p_r;
	double u_l, u_r;
	double v_l. v_r;
	double w_l, w_r;

	p_l = q_l[index[0]]; u_l = q_l[index[1]]; v_l = q_l[index[2]];
	p_r = q_r[index[0]]; u_r = q_r[index[1]]; v_r = q_r[index[2]];

	w_l = (p_l + dg_fun::C * u_l) / 2.0;
	w_r = (p_r - dg_fun::C * u_r) / 2.0;

	n_flux[index[0]] = dg_fun::C * (w_l - w_r) * normal;
	n_flux[index[1]] = (w_l + w_r) * normal;
	n_flux[index[2]] = 0.0;
	

}