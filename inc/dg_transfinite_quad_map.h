#ifndef DG_TRANSFINITE_QUAD_MAP_H
#define DG_TRANSFINITE_QUAD_MAP_H

#include "dg_curve_interpolant.h"

extern CurveInterpolant B[4];

//Algo 98
void TransfiniteQuadMap(double xi, double eta, CurveInterpolant boundaries[4], double &x, double &y);

#endif
