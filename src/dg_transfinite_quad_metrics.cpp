#include "dg_transfinite_quad_metrics.h"
#include <cmath>

void TransfiniteQuadMetrics(double xi, double eta, CurveInterpolant boundaries[4], double &x_xi, double &x_eta, double &y_xi, double &y_eta){
    double x_1, y_1, x_2, y_2, x_3, y_3, x_4, y_4;
    double X1, Y1, X2, Y2, X3, Y3, X4, Y4;
    double X1_prime, Y1_prime, X2_prime, Y2_prime, X3_prime, Y3_prime, X4_prime, Y4_prime;


    boundaries[0].EvalAt(-1, x_1, y_1);
    boundaries[0].EvalAt(1, x_2, y_2);
    boundaries[2].EvalAt(1, x_3, y_3);
    boundaries[2].EvalAt(-1, x_4, y_4);

    boundaries[0].EvalAt(xi, X1, Y1);
    boundaries[1].EvalAt(eta, X2, Y2);
    boundaries[2].EvalAt(xi, X3, Y3);
    boundaries[3].EvalAt(eta, X4, Y4);

    boundaries[0].DerivativeAt(xi, X1_prime, Y1_prime);
    boundaries[1].DerivativeAt(eta, X2_prime, Y2_prime);
    boundaries[2].DerivativeAt(xi, X3_prime, Y3_prime);
    boundaries[3].DerivativeAt(eta, X4_prime, Y4_prime);

    x_xi = 0.5*(X2 - X4 + (1-eta)*X1_prime + (1+eta)*X3_prime) - 0.25*((1-eta)*(x_2 - x_1) + (1+eta)*(x_3 - x_4));
    y_xi = 0.5*(Y2 - Y4 + (1-eta)*Y1_prime + (1+eta)*Y3_prime) - 0.25*((1-eta)*(y_2 - y_1) + (1+eta)*(y_3 - y_4));

    x_eta = 0.5*((1-xi)*X4_prime + (1+xi)*X2_prime + X3 - X1) - 0.25*( (1-xi)*(x_4 - x_1) + (1+xi)*(x_3 - x_2) );
    y_eta = 0.5*((1-xi)*Y4_prime + (1+xi)*Y2_prime + Y3 - Y1) - 0.25*( (1-xi)*(y_4 - y_1) + (1+xi)*(y_3 - y_2) );

    //just added this
    if(std::abs(x_xi) <= 1e-10){
        x_xi = 0;
    }
    if(std::abs(x_eta) <= 1e-10){
        x_eta = 0;
    }
    if(std::abs(y_xi) <= 1e-10){
        y_xi = 0;
    }
    if(std::abs(y_eta) <= 1e-10){
        y_eta = 0;
    }
    
} 