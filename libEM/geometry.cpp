/**
 * $Id$
 */
#include "geometry.h"

using namespace EMAN;


bool EMAN::operator<(const Pixel& p1, const Pixel& p2)
{
	if (p1.value < p2.value) {
		return true;
	}
	return false;
}

bool EMAN::operator==(const Pixel& p1, const Pixel& p2)
{
	if (p1.x == p2.x && p1.y == p2.y && p1.z == p2.z && p1.value == p2.value) {
		return true;
	}
	return false;
}

bool EMAN::operator!=(const Pixel& p1, const Pixel& p2)
{
	return !(p1 == p2);
}


bool Region::inside_region() const
{
	if (size[0] >= 0 && size[1] >= 0 && size[2] >= 0) {
		return true;
	}

	return false;
}

bool Region::inside_region(const FloatPoint & p) const
{
	if (p.get_ndim() == 1) {
		return inside_region(p[0]);
	}
  

	if (p.get_ndim() == 2) {
		return inside_region(p[0], p[1]);
	}

	return inside_region(p[0], p[1], p[2]);
}

bool Region::inside_region(float x) const
{
	if (size[0] >= 0 && origin[0] <= x && 
		(origin[0] + size[0]) >= x ) {
		return true;
	}
	return false;
}


bool Region::inside_region(float x, float y) const
{
	if (size[0] >= 0 && size[1] >= 0 &&
		origin[0] <= x && origin[1] <= y &&
		(origin[0] + size[0]) >= x && (origin[1] + size[1]) >= y) {
		return true;
	}
	return false;
}

bool Region::inside_region(float x, float y, float z) const
{
	if (size[0] >= 0 && size[1] >= 0 && size[2] >= 0 &&
		origin[0] <= x && origin[1] <= y && origin[2] <= z &&
		(origin[0] + size[0]) >= x &&
		(origin[1] + size[1]) >= y && (origin[2] + size[2]) >= z) {
		return true;
	}
	return false;
}


bool Region::is_region_in_box(const FloatSize & box) const
{
	if (size[0] >= 0 && size[1] >= 0 && size[2] >= 0 &&
		origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0 &&
		(origin[0] + size[0]) <= box[0] &&
		(origin[1] + size[1]) <= box[1] &&
		(origin[2] + size[2]) <= box[2]) {
		return true;
	}
	
	return false;	
}



string Region::get_string() const
{
	char str[1028];
	int ndim = origin.get_ndim();

	if (ndim == 2) {
		sprintf(str, "(%2.1f, %2.1f; %2.1f, %2.1f)",
				origin[0], origin[1], size[0], size[1]);
	}
	else if (ndim == 3) {
		sprintf(str, "(%2.1f, %2.1f, %2.1f; %2.1f, %2.1f, %2.1f)",
				origin[0], origin[1], origin[2], size[0], size[1], size[2]);
	}

	return string(str);
}
