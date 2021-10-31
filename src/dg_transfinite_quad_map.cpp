#include "dg_transfinite_quad_map.h"

//double check this
void TransfiniteQuadMap(double xi, double eta, CurveInterpolant boundaries[4], double &x, double &y){
    double x_1, y_1, x_2, y_2, x_3, y_3, x_4, y_4;
    double X1, Y1, X2, Y2, X3, Y3, X4, Y4;

    boundaries[0].EvalAt(-1, x_1, y_1);
    boundaries[1].EvalAt(1, x_2, y_2);
    boundaries[2].EvalAt(1, x_3, y_3);
    boundaries[3].EvalAt(-1, x_4, y_4);

    boundaries[0].EvalAt(xi, X1, Y1);
    boundaries[1].EvalAt(eta, X2, Y2);
    boundaries[2].EvalAt(xi, X3, Y3);
    boundaries[3].EvalAt(eta, X4, Y4);

    x = 0.5*((1-xi)*X4 + (1+xi)*X2 + (1-eta)*X1 + (1+eta)*X3) - 0.25*( (1-xi) *((1-eta)*x_1 + (1+eta*x_4)) + (1+xi) * ( (1-eta)*x_2 + (1+eta)*x_3 ) );
    y = 0.5*((1-xi)*Y4 + (1+xi)*Y2 + (1-eta)*Y1 + (1+eta)*Y3) - 0.25*( (1-xi) *((1-eta)*y_1 + (1+eta*y_4)) + (1+xi) * ( (1-eta)*y_2 + (1+eta)*y_3 ) );

}