/**
 * $Id$
 */
#include "geometry.h"

using namespace EMAN;


bool Region::inside_region() const
{
    if (size.xsize >= 0 && size.ysize >= 0 && size.zsize >= 0) {
	return true;
    }

    return false;
}

bool Region::inside_region(const Point<float> & p) const
{
    if (p.get_ndim() == 2) {
	return inside_region(p.x, p.y);
    }

    return inside_region(p.x, p.y, p.z);
}

bool Region::inside_region(float x, float y) const
{
    if (size.xsize >= 0 && size.ysize >= 0 && 
	origin.x <= x && origin.y <= y &&
	(origin.x + size.xsize) >= x &&
	(origin.y + size.ysize) >= y) {
	return true;
    }
    return false;
}

bool Region::inside_region(float x, float y, float z) const
{
    if (size.xsize >= 0 && size.ysize >= 0 && size.zsize >= 0 && 
	origin.x <= x && origin.y <= y && origin.z <= z &&
	(origin.x + size.xsize) >= x &&
	(origin.y + size.ysize) >= y &&
	(origin.z + size.zsize) >= z) {
	return true;
    }
    return false;
}

string Region::get_string() const
{
    char str[1028];
    int ndim = origin.get_ndim();

    if (ndim == 2) {
	sprintf(str, "(%2.1f, %2.1f; %d, %d)", origin.x, origin.y, size.xsize, size.ysize);
    }
    else if (ndim == 3) {
	sprintf(str, "(%2.1f, %2.1f, %2.1f; %d, %d, %d)",
		origin.x, origin.y, origin.z, size.xsize, size.ysize, size.zsize);
    }

    return string(str);
}
