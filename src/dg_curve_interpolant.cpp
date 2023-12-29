#include "dg_curve_interpolant.h"
#include "dg_basis.h"
#include <iostream>

//works
//nodes are GL points
//x, y are points along the actual curves
//w is barycentric weights
CurveInterpolant::CurveInterpolant(int N, std::vector<double> nodesparam, std::vector<double> xparam, std::vector<double> yparam){
    this->N = N;
    nodes = std::vector<double>(N+1);
    x = std::vector<double>(N+1);
    y = std::vector<double>(N+1);
    w = std::vector<double>(N+1);

    for(int i = 0; i <= N; i++){
        nodes[i] = nodesparam[i];
        x[i] = xparam[i];
        y[i] = yparam[i];
    }

    BARW(N, nodes, w);
}

void CurveInterpolant::EvalAt(double s, double &x, double &y){
    x = Lagrange_interpolation(this->N, s, this->nodes, this->x, this->w);
    y = Lagrange_interpolation(this->N, s, this->nodes, this->y, this->w);
}

void CurveInterpolant::DerivativeAt(double s, double &xprime, double &yprime){
    xprime = lagrangeInterpDeriv(this->N, s, this->nodes, this->x, this->w);
    yprime = lagrangeInterpDeriv(this->N, s, this->nodes, this->y, this->w);
}



CurveInterpolant::CurveInterpolant(){
    
}

