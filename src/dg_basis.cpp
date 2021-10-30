#include <iostream>
#include <vector>
#include <limits>	// epsilon
#include <cmath>	// sqrt
#include "dg_single_index.h"
#include <cassert>


const double pi = 4.0 * atan(1.0); 

// forward declaration----------------------------
bool Almost_equal(double a, double b);

void Mth_order_polynomial_derivative_matrix(int n, int mth_der, std::vector<double>& x, std::vector<double>& der, 
						std::vector<double>& bary);

void GL(int n, std::vector<double>& gl_p, std::vector<double>& gl_w);

void BARW(int n, std::vector<double>& x, std::vector<double>& bary);

double Lagrange_interpolation(int n, double x, std::vector<double>& xi, 
				std::vector<double>& f, std::vector<double>& w, std::vector<int>& index);

void Lagrange_interpolating_polynomial(int n, double target_p, std::vector<double>& x, std::vector<double>& bary,
					 std::vector<double>& lag );

void Legendre_polynomial_and_derivative(int n, double x, double& q, double& dq);

double Interpolate_to_boundary(int n, std::vector<double>& q, std::vector<double>& lag);

void Matrix_vector_multiplication(int n, std::vector<double>& d, std::vector<double>& f, std::vector<double>& out);
//------------------------------------------------

/// @brief 
/// compute matrix vector multiplication.
/// d * f = der. Algorithm 19.
/// @param n polynomial order
/// @param d coefficient martrix (usually derivative matrix), d[n * n].
/// @param f vector, f[n].
/// @param out output. (usually the derivative of interpolation).
void Matrix_vector_multiplication(int n, std::vector<double>& d, std::vector<double>& f, std::vector<double>& out){

	for(int i = 0; i <= n; ++i){
	
		double t{};	// intermediate variable

		for(int j = 0; j <= n; ++j){

			int m = Get_single_index(i, j, n + 1);

			t += d[m] * f[j];
		}
		
		out[i] = t;
	}

}



/// @brief 
/// compute matrix vector multiplication.
/// d * f = der. Use BLAS library. 
/// @param n polynomial order
/// @param d coefficient martrix (usually derivative matrix), d[n * n].
/// @param f vector, f[n].
/// @param der output. (usually the derivative of interpolation), der[n].
void Matrix_vector_multiplication_blas(int n, double* d, double* f, double* der){


}


/// @brief
/// Legendre polynomial of degree k and its derivative using the three
/// term recursive. 
/// Algorithm 22
/// @param n polynomial order
/// @param x GL points
/// @param q Legendre polynomial of degree k
/// @param dq Derivative of Legendre polynomial  
void Legendre_polynomial_and_derivative(int n, double x, double& q, double& dq){

	if(n == 0){
		q = 1.0;
		dq = 0.0;
	}
	else if(n == 1){
		q = x;
		dq = 1.0;

	}
	else{
		double q_m2 = 1.0;
		double q_m1 = x;
		double dq_m2 = 0.0;
		double dq_m1 = 1.0;
		
		for(int k = 2; k <= n; ++k ){
			q = (double)(2 * k - 1) * x * q_m1 / (double)(k) - (double)(k -1) * q_m2 / (double)(k);
			dq = dq_m2 + (double)(2 * k - 1) * q_m1;
			q_m2 = q_m1;
			q_m1 = q;
			dq_m2 = dq_m1;
			dq_m1 = dq;
		}

	}
	
}


/// @brief
/// Compute the Gauss Legendre nodes and weights. Algprithm 23
/// @param n polynomial order
/// @param gl_p GL points
/// @param gl_w GL weigths
void GL(int n, std::vector<double>& gl_p, std::vector<double>& gl_w){

	double delta;
	double q, dq, tol;

	tol = 4.0 * std::numeric_limits<double>::epsilon();
	q = 0.0;
	dq = 0.0;
	
	if(n == 0){
		gl_p[0] = 0.0;
		gl_w[0] = 2.0;
	}
	else if(n == 1){

		gl_p[0] = - sqrt(1.0 / 3.0);
		gl_w[0] = 1.0;
		gl_p[1] = sqrt(1.0 / 3.0);
		gl_w[1] = 1.0;

	}
	else{
		for(int j = 0; j <= ((n+1)/2)-1; ++j ){
			// initial guess
			gl_p[j] = - cos(pi * (double)(2 * j + 1)/(double)(2 * n + 2));

			// iterative method
			delta = 10000000.0;
			while(true){
				
				Legendre_polynomial_and_derivative(n+1, gl_p[j], q, dq);

				delta = - q / dq;
				gl_p[j] = gl_p[j] + delta;
				
				if(std::abs(delta) <= (tol*std::abs(gl_p[j])))
					break;
			}

			Legendre_polynomial_and_derivative(n+1, gl_p[j], q, dq);
				
			gl_p[n - j] = - gl_p[j];
			gl_w[j] = 2.0 / ((1.0 - pow(gl_p[j], 2.0)) * pow(dq, 2.0));
			gl_w[n - j] = gl_w[j];

		}

	}


	if(n % 2 == 0){
		double gl_p_now = 0.0;
		Legendre_polynomial_and_derivative(n+1, gl_p_now, q, dq);

		gl_p[n/2] = 0.0;
		gl_w[n/2] = 2.0 / pow(dq, 2.0);

	}


}

/// @brief 
/// Barycentric weights for Lagrange Iterpolation. Algorithm 30.
/// @param n polynomial order
/// @param x spectral points
/// @param bary barycentric weights
void BARW(int n, std::vector<double>& x, std::vector<double>& bary){
	
	for(int i = 0; i <= n; ++i){
		bary[i] = 1.0;
	}

	for(int j = 1; j <= n; ++j){
		for(int k = 0; k <= j-1; ++k){
			bary[k] = bary[k] * (x[k] - x[j]);
			bary[j] = bary[j] * (x[j] - x[k]);

		}

	}

	for(int j = 0; j <= n; ++j){

		bary[j] = 1.0 / bary[j];

	}


}

/// @brief
/// Lagrange interpolant from Barycentric form. Algorithm 31.  
/// @param n polynomial order. 
/// @param x solution at the target point x (on reference space).
/// @param xi GL points. 
/// @param f solution set.
/// @param w barycentric weights. 
/// @param index solution indices. 
double Lagrange_interpolation(int n, double x, std::vector<double>& xi, 
				std::vector<double>& f, std::vector<double>& w, std::vector<int>& index){
		
	double numerator{};
	double denominator{};

	for(int j = 0; j <= n; ++j){

		if(Almost_equal(x, xi[j])){	// on GL point
			return f[index[j]];
		}

		double t = w[j] / (x - xi[j]);

		numerator += t * f[index[j]];

		denominator += t;
	}

	return numerator / denominator;
}


/// @brief
/// Lagrange interpolating polynoimial value at target point. Algorithm 34.
/// Return the interpolation value on point x. 
/// @param n polynomial order
/// @param target_p target point 
/// @param x GL points
/// @param bary barycentric points
/// @param lag lagrange interpolating values at target point. 
void Lagrange_interpolating_polynomial(int n, double target_p, std::vector<double>& x, std::vector<double>& bary,
					 std::vector<double>& lag ){

	double s = 0.0;
	bool match = false;

	for(int j = 0; j <= n; ++j){
		lag[j] = 0.0;
			
		bool flag = Almost_equal(target_p, x[j]);
		if(flag){
			lag[j] = 1.0;
			match = true;
		}
	}	
	

	if(match){
		return;	
	}


	for(int j = 0; j <= n; ++j){
		double t = bary[j] / (target_p - x[j]);
		lag[j] = t;
		s += t;
	}

	// maybe use blas dscal
	for(int j = 0; j <= n; ++j){
		lag[j] = lag[j] / s;

	}

	

}

/// @brief
/// M-th order derivative matrix. Algorithm 36-37
/// @param n polynomial order
/// @param mth_der m-th order polynomial derivative
/// @param x spectral points
/// @param der m-th order derivative matrix
/// @param bary Barycentric weights. 
void Mth_order_polynomial_derivative_matrix(int n, int mth_der, std::vector<double>& x, std::vector<double>& der, 
						std::vector<double>& bary){
	
	assert(mth_der == 1 && "Now the Mth_order_poynomial_derivative_maxtri() function can only compute the 1st"  
					" order derivative. ");

	// mth-order == 1
	for(int i = 0; i <= n; ++i){
		int inode = Get_single_index(i, i, n+1);
		der[inode] = 0.0;
		
		for(int j = 0; j <= n; ++j){
			if(j != i){
				int node1 = Get_single_index(i, j, n+1);
				int node2 = Get_single_index(i, i, n+1);

				der[node1] = bary[j] / bary[i] / (x[i] - x[j]);
				der[node2] = der[node2] - der[node1];

			}

		}

	}	
	

	if(mth_der == 1){
		return;
	}

	// mth_order > 1
	//------------
//	aux = der
//	//------------
//	for(int k = 2; k <= mth_der; ++k){
//		for(int i = 0; i <= n;, ++i){
//			
//			int inode = Get_single_index(i, i, n+1);
//			der[inode] = 0.0;
//
//			for(int j = 0; j <= n; ++j){
//				if(j != i){
//					
//					int node1 = Get_single_index(i, j, n+1);
//					int node2 = Get_single_index(i, i, n+1);
//
//					der[node1] = (double)(k) / (x[i] - x[j]) * 
//							(bary[j] / bary[i] * aux[node2] - aux[node1]);
//
//					der[node2] = der[node2] - der[node1];
//				}
//
//			}
//
//
//		}
//
//	}
	

}

/// @brief
/// Interpolate the interior solution to the element boundary/interface.
/// @param n polynomial order.
/// @param q solution array. 
/// @param lag Lagrange interpolation polynomial.
double Interpolate_to_boundary(int n, std::vector<double>& q, std::vector<double>& lag){

	double inter{};

	for(int j = 0; j < n + 1; ++j){

		inter += lag[j] * q[j];

	}

	return inter;

}


/// @brief
/// Testing equality of two floating point numbers. Algorithm 139
/// @param a number 1.
/// @param b number 2.
bool Almost_equal(double a, double b){
	
	if((a == 0.0) || (b == 0.0)){
 		double tol  = 2.0 * std::numeric_limits<double>::epsilon();
		if(std::abs(a - b) <= tol ){
			return true;

		}
		else{
			return false;
		}
	}
	else{
		double tol1 = std::abs(a) * std::numeric_limits<double>::epsilon();
		double tol2 = std::abs(b) * std::numeric_limits<double>::epsilon();
		if(((std::abs(a - b)) <= tol1) && ((std::abs(a - b)) <= tol2)){
			return true;

		}
		else{

			return false;
		}

	}

}
