#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <cmath>
#include "emdata.h"
#include "util.h"

/* #include "fundamentals.h" */
/* #include "lapackblas.h" */

namespace {


/* def  ang_n(peakp, mode, maxrin): */
/* 	"""Calculate angle based on the position of the peak */
/* 	""" */
/* 	from math import fmod,cos,sin */
/* 	if (mode == 'f' or mode == 'F'): return fmod(((peakp-1.0) / maxrin+1.0)*360.0,360.0) */
/* 	else:                            return fmod(((peakp-1.0) / maxrin+1.0)*180.0,180.0) */

float ang_n(float peakp, string mode, int maxrin)
{
    if (mode == "f" || mode == "F") 
        return fmodf(((peakp-1.0) / maxrin+1.0)*360.0,360.0);
    else
        return fmodf(((peakp-1.0) / maxrin+1.0)*180.0,180.0);
}


EMAN::Dict apmq(EMAN::EMData* image, boost::python::list const &crefim_list,
                int xrng, int yrng, int step, string mode,
                vector< int >numr, int cnx, int cny) {

    // Determine shift and rotation between image and many reference
    // images (crefim, weights have to be applied) quadratic
    // interpolation  
    
    
    // Manually extract.
    vector< EMAN::EMData* > crefim;
    std::size_t crefim_len = PyObject_Length(crefim_list.ptr());
    crefim.reserve(crefim_len);

    for(std::size_t i=0;i<crefim_len;i++) {
        boost::python::extract<EMAN::EMData*> proxy(crefim_list[i]);
        crefim.push_back(proxy());
    }

    float peak = -1.0E23;
    int ky = int(2*yrng/step+0.5)/2; 
    int kx = int(2*xrng/step+0.5)/2;
    //for i in xrange(-ky, ky+1):
    int iref, nref=0, mirror=0, iy, ix, sx=0, sy=0;
    float ang=0.0f;
    for (int i = -ky; i <= ky; i++) {
        iy = i * step ;
        // for  j in xrange(-kx, kx+1):
        for (int j = -kx; j <= kx; j++) {
            ix = j*step ; 
            EMAN::EMData* cimage =
                EMAN::Util::Polar2Dm(image, cnx+ix, cny+iy, numr, mode);
            EMAN::Util::Frngs(cimage, numr);
            //  compare with all reference images
            // for iref in xrange(len(crefim)): 
            for ( iref = 0; iref < (int)crefim_len; iref++) {
                EMAN::Dict retvals =
                    EMAN::Util::Crosrng_ms(crefim[iref], cimage, numr);  
                double qn = retvals["qn"];
                double qm = retvals["qm"];
                if(qn >= peak || qm >= peak) {
                    sx = -ix;
                    sy = -iy;
                    nref = iref;
                    if (qn >= qm) {
                        ang = ang_n(retvals["tot"], mode, numr[numr.size()-1]);
                        peak = qn;
                        mirror = 0;
                        } 
                    else {
                        ang = ang_n(retvals["tmt"], mode, numr[numr.size()-1]);
                        peak = qm; 
                        mirror = 1;
                    }
                }
            }
        }
    }
    float co, so, sxs, sys;
    co =  cos(ang*pi/180.0);
    so = -sin(ang*pi/180.0);
    sxs = sx*co - sy*so;
    sys = sx*so + sy*co;

    EMAN::Dict retvals;
    retvals["ang"] = ang;
    retvals["sxs"] = sxs;
    retvals["sys"] = sys;
    retvals["mirror"] = mirror;
    retvals["peak"] = peak;
    retvals["nref"] = nref;
    return retvals;

}
}

namespace EMAN { 
namespace boost_python {

void wrap_apmq()
{
    using namespace boost::python;
    def("sparx_apmq", apmq);
}
}}
