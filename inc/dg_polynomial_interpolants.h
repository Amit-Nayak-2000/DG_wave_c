#ifndef DG_POLYNOMIAL_INTERPOLANTS_H
#define DG_POLYNOMIAL_INTERPOLANTS_H

#include <vector>

//algo 30
void BaryWeights(std::vector<double> x, std::vector<double> &w);

//algo 31
double LagrangeInterp(double x, std::vector<double> xj, std::vector<double> fj, std::vector<double> wj);

//algo 36
double lagrangeInterpDeriv(double x, std::vector<double> xj, std::vector<double> fj, std::vector<double> wj);

#endif