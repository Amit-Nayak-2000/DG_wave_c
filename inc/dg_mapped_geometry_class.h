#ifndef MAPPED_GEOMETRY_CLASS
#define MAPPED_GEOMETRY_CLASS

#include "dg_curve_interpolant.h"
#include <vector>

//Algo 101
class MappedGeometry {
    public:

   //polynomial order
    int N;

    //0 to N + 1 by 0 to N + 1
    std::vector<double> x_node, y_node;

    //0 to N + 1 by 0 to 3
    std::vector<double> x_boundary, y_boundary;

    //0 to N + 1 by 0 to N + 1
    std::vector<double> delx_delxi, delx_deleta, dely_delxi, dely_deleta;

    //0 to N + 1 by 0 to N + 1
    std::vector<double> jacobian;
    
    //0 to N + 1 by 0 to 3
    std::vector<double> boundary_normal;
    std::vector<double> boundary_normal_y;

    //0 to N + 1 by 0 to 3
    std::vector<double> scaling_factor;

    //constructors
    MappedGeometry();

    MappedGeometry(CurveInterpolant boundaries[4], std::vector<double> gl, int porder);

    void updateOrder(int N);
    
    void initMetrics(CurveInterpolant boundaries[4], std::vector<double> gl);

    void recalcscalenorms(CurveInterpolant boundaries[4], std::vector<double> gl);

    void outputMetrics();


};

#endif