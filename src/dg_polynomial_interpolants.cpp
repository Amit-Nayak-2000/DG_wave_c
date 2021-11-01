#include "dg_polynomial_interpolants.h"
#include "dg_almost_equal.h"

void BaryWeights(std::vector<double> x, std::vector<double> &w){
    int N = x.size();
    for(int j = 0; j <= N; j++){
        w[j] = 1;
    }

    for(int j = 1; j <= N; j++){
        for(int k = 0; k <= j - 1; k++){
            w[k] = w[k]*(x[k] - x[j]);
            w[j] = w[j]*(x[j] - x[k]);
        }
    }

    for(int j = 0; j <= N; j++){
        w[j] = 1/w[j];
    }
}

double LagrangeInterp(double x, std::vector<double> xj, std::vector<double> fj, std::vector<double> wj){

    double numerator = 0;
    double denominator = 0;
    double t;

    for(int j = 0; j <= xj.size(); j++){
        if(AlmostEqual(x, xj[j])){
            return fj[j];
        }

        t = wj[j]/(x - xj[j]);

        numerator = numerator + t*fj[j];
        denominator = denominator + t;

    }

    return (numerator / denominator);

}

double lagrangeInterpDeriv(double x, std::vector<double> xj, std::vector<double> fj, std::vector<double> wj){
    bool atNode = false;
    double numerator = 0;
    double denominator;
    double p;
    int i = -1;

    for(int j = 0; j <= xj.size(); j++){
        if(AlmostEqual(x, xj[j])){
            atNode = true;
            p = fj[j];
            denominator = -wj[j];
            i = j;
        }

    }

    if(atNode){
        for(int j = 0; j <= xj.size(); j++){
            if(j != i){
                numerator = numerator + wj[j]*(p-fj[j])/(x-xj[j]);
            }
        }
    }
    else{
        denominator = 0;
        p = LagrangeInterp(x, xj, fj, wj);

        for(int j = 0; j <= xj.size(); j++){
            double t = wj[j] / (x - xj[j]);

            numerator = numerator + t*(p-fj[j])/(x-xj[j]);

            denominator = denominator + t;
        }
    }

    return (numerator / denominator);

}