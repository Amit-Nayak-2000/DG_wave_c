#ifndef DG_CURVE_INTERPOLANT_H
#define DG_CURVE_INTERPOLANT_H

#include <vector>

//page 227
//Algorithm 96 and 97
class CurveInterpolant{
    public:
        int N;
        std::vector<double> nodes;
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> w;

        CurveInterpolant(int N, std::vector<double> nodesparam, std::vector<double> xparam, std::vector<double> yparam);

        CurveInterpolant();

        void EvalAt(double s, double &x, double &y);

        void DerivativeAt(double s, double &xprime, double &yprime);

};

#endif