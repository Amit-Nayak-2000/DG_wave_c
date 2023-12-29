#include "dg_flux_vector.h"
#include <vector>
#include "dg_param.h"
#include <iostream>

/// @brief
/// Compute the horizontal fluxes. X direction. 
/// @param q solutions.
/// @param xf horizontal numerical fluxes.
/// @parma solution indeices (3 equaitons).
void xflux(std::unordered_map<int, std::vector<double>>& q, 
		std::unordered_map<int, std::vector<double>>& xf, int index){

	xf[0].push_back(dg_fun::C * dg_fun::C * q[1][index]);

	xf[1].push_back(q[0][index]);

	xf[2].push_back(0.0);

	// std::cout << "Xflux" << std::endl;
	// std::cout << (dg_fun::C * dg_fun::C * q[1][index]) << ", " << q[0][index] << ", " << 0.0 << std::endl;

}

/// @brief
/// Compute the horizontal fluxes. Y direction. 
/// @param q solutions.
/// @param xf horizontal numerical fluxes.
/// @parma solution indeices (3 equaitons).
void yflux(std::unordered_map<int, std::vector<double>>& q, 
		std::unordered_map<int, std::vector<double>>& yf, int index){

	yf[0].push_back(dg_fun::C * dg_fun::C * q[2][index]);

	yf[1].push_back(0.0);

	yf[2].push_back(q[0][index]);

	// std::cout << "Yflux" << std::endl;
	// std::cout << (dg_fun::C * dg_fun::C * q[2][index]) << ", " << 0.0 << ", " << q[0][index] << std::endl;

}

void xflux_updated(std::unordered_map<int, std::vector<double>>& q, 
		std::unordered_map<int, std::vector<double>>& xf, int index, double y_eta, double x_eta){

	//y_eta * xflux - x_eta * yflux

	xf[0].push_back( (y_eta*(dg_fun::C * dg_fun::C * q[1][index])) - (x_eta*(dg_fun::C * dg_fun::C * q[2][index])) );

	xf[1].push_back( (y_eta*(q[0][index])) - (x_eta*(0.0)) );

	xf[2].push_back( (y_eta*(0.0)) - (x_eta*(q[0][index])) );

	// std::cout << "Xflux updated" << std::endl;
	// std::cout << ((y_eta*(dg_fun::C * dg_fun::C * q[1][index])) - (x_eta*(dg_fun::C * dg_fun::C * q[2][index]))) << ", " << ( (y_eta*(q[0][index])) - (x_eta*(0.0)) ) << ", " << ( (y_eta*(0.0)) - (x_eta*(q[0][index])) ) << std::endl;

}

void yflux_updated(std::unordered_map<int, std::vector<double>>& q, 
		std::unordered_map<int, std::vector<double>>& yf, int index, double y_xi, double x_xi){

	//-y_xi * xflux + x_xi * yflux

	yf[0].push_back( (-y_xi*(dg_fun::C * dg_fun::C * q[1][index])) + (x_xi*(dg_fun::C * dg_fun::C * q[2][index])) );

	yf[1].push_back( (-y_xi*(q[0][index])) + (x_xi*(0.0)) );

	yf[2].push_back( (-y_xi*(0.0)) + (x_xi*(q[0][index])) );

	// std::cout << "Yflux updated" << std::endl;
	// std::cout << ( (-y_xi*(dg_fun::C * dg_fun::C * q[1][index])) + (x_xi*(dg_fun::C * dg_fun::C * q[2][index])) ) << ", " << ( (-y_xi*(q[0][index])) + (x_xi*(0.0)) ) << ", " << ( (-y_xi*(0.0)) + (x_xi*(q[0][index])) ) << std::endl;

}
