#include "cmp_template.h"
#include "transform.h"


using namespace EMAN;

float XYZCmp::cmp(EMData *image, Transform *transform) const
{
    if (!image) {
	return 0;
    }
	    
#if 0
    EMData *with = params["with"];
    int param1 = params["param1"];
    float param2 = params["param2"];
#endif
    
    return 0;
}
