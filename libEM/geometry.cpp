/**
 * $Id$
 */
#include "geometry.h"

using namespace EMAN;


bool EMAN::operator<(const Pixel& p1, const Pixel& p2)
{
	if (p1.value < p2.value) {
		return false;
	}
	return true;
}

bool Region::inside_region() const
{
	if (size.x >= 0 && size.y >= 0 && size.z >= 0) {
		return true;
	}

	return false;
}

bool Region::inside_region(const FloatPoint & p) const
{
	if (p.get_ndim() == 2) {
		return inside_region(p.x, p.y);
	}

	return inside_region(p.x, p.y, p.z);
}

bool Region::inside_region(float x, float y) const
{
	if (size.x >= 0 && size.y >= 0 &&
		origin.x <= x && origin.y <= y &&
		(origin.x + size.x) >= x && (origin.y + size.y) >= y) {
		return true;
	}
	return false;
}

bool Region::inside_region(float x, float y, float z) const
{
	if (size.x >= 0 && size.y >= 0 && size.z >= 0 &&
		origin.x <= x && origin.y <= y && origin.z <= z &&
		(origin.x + size.x) >= x &&
		(origin.y + size.y) >= y && (origin.z + size.z) >= z) {
		return true;
	}
	return false;
}

string Region::get_string() const
{
	char str[1028];
	int ndim = origin.get_ndim();

	if (ndim == 2) {
		sprintf(str, "(%2.1f, %2.1f; %d, %d)", origin.x, origin.y, size.x, size.y);
	}
	else if (ndim == 3) {
		sprintf(str, "(%2.1f, %2.1f, %2.1f; %d, %d, %d)",
				origin.x, origin.y, origin.z, size.x, size.y, size.z);
	}

	return string(str);
}
