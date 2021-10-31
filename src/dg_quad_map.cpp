#include "dg_quad_map.h"

void QuadMap(double xj[], double yj[], double xi, double eta, double &x, double &y){
    x = 0.25 * (xj[0]*(1-xi)*(1-eta) + xj[1]*(xi+1)*(1-eta) + xj[2]*(xi+1)*(eta+1) + xj[3]*(1-xi)*(eta+1));
    y = 0.25 * (yj[0]*(1-xi)*(1-eta) + yj[1]*(xi+1)*(1-eta) + yj[2]*(xi+1)*(eta+1) + yj[3]*(1-xi)*(eta+1));
}