#include "dg_curve_interpolant.h"
#include "dg_polynomial_interpolants.h"


CurveInterpolant::CurveInterpolant(int N, std::vector<double> nodesparam, std::vector<double> xparam, std::vector<double> yparam){
    this->N = N;
    for(int i = 0; i <= N; i++){
        nodes.push_back(nodesparam[i]);
        x.push_back(xparam[i]);
        y.push_back(yparam[i]);
    }

    BaryWeights(nodes, w);
}

void CurveInterpolant::EvalAt(double s, double &x, double &y){
    x = LagrangeInterp(s, this->nodes, this->x, this->w);
    y = LagrangeInterp(s, this->nodes, this->y, this->w);
}

void CurveInterpolant::DerivativeAt(double s, double &xprime, double &yprime){
    xprime = lagrangeInterpDeriv(s, this->nodes, this->x, this->w);
    yprime = lagrangeInterpDeriv(s, this->nodes, this->y, this->w);
}





