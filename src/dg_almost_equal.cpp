#include "dg_almost_equal.h"
#include <cfloat>
#include <cmath>

bool AlmostEqual(double a, double b){
    bool result;
    if(a == 0 || b == 0){
        if(abs(a - b) <= 2*DBL_EPSILON){
            result = true;
        }
        else{
            result = false;
        }
    }
    else{
        if(abs(a - b) <= DBL_EPSILON*abs(a) && abs(a - b) <= DBL_EPSILON*abs(b)){
            result = true;
        }
        else{
            result = false;
        }
    }

    return result;
}