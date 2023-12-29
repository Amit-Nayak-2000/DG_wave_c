#ifndef DG_INTERPOLATE_TO_NEW_POINT_H
#define DG_INTERPOLATE_TO_NEW_POINT_H

#include <array>
#include <vector>
#include <unordered_map>

void Polynomial_interpolate_matrix(std::vector<double>& x, std::vector<double>& xi, std::vector<double>& T);

void Solutions_to_children(std::array<long long int, 4>& keys, long long int p_key);

void Interpolate_to_new_points(int m, int n, std::vector<double>& T, 
				std::vector<double>& f, std::vector<double>& new_f, int start_old, int start_new, int interval);

void Solution_back_to_parent(std::array<long long int, 4>& keys, long long int p_key);

void coarsetofine(std::vector<double> &target, int new_order, int currentorder, std::vector<double> &T);

void coarsetofine1D(std::vector<double> &target, int new_order, int currentorder, const std::vector<double> &T);

void CoursetoFineInterp(int nold, int nnew, std::vector<double> &x, std::vector<double> &y, std::unordered_map<int, std::vector<double>> &f, 
std::vector<double> &xnew,  std::vector<double> &ynew, std::unordered_map<int, std::vector<double>> &fnew);

std::vector<double> linarray(double min, double max, int N);

// void CoursetoFineInterp(int nold, int nnew, std::vector<double> &x, std::vector<double> &y, std::vector<double> &f, 
// std::vector<double> &xnew,  std::vector<double> &ynew,  std::vector<double> &fnew);

void CoursetoFineInterpCoords(int nold, int nnew, std::vector<double> &x, std::vector<double> &y, std::unordered_map<int, std::vector<double>> &f, 
std::vector<double> &xnew,  std::vector<double> &ynew, std::unordered_map<int, std::vector<double>> &fnew);

#endif
