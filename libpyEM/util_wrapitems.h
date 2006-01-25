#include "util.h"
#include <boost/python.hpp>

#ifndef UTIL_WRAPPER_H
#define UTIL_WRAPPER_H
using EMAN::EMData;
boost::python::tuple 
util_Crosrng_e(EMData* circ1, EMData* circ2, 
               vector<int> numr, int neg){
    double qn = 0.;
    float tot = 0.f;
    boost::tie(qn, tot, neg) = EMAN::Util::Crosrng_e(circ1, circ2, numr, neg);
    return boost::python::make_tuple(qn, tot, neg);
}
#endif // UTIL_WRAPPER_H
