#include "testutil.h"
#include "emdata.h"
#include "emobject.h"

#include <vector>
#include <string>
#include <map>

using std::vector;
using std::string;
using std::map;

using namespace EMAN;

const char* TestUtil::get_debug_image(const char* imagename)
{
	if (!imagename) {
		return "";
	}
	
	string fullpath = "";
	char * imgpath = getenv("DEBUG_IMAGE_PATH");
	if (imgpath) {
		fullpath = string(imgpath);
	}
	else {
		fullpath = string(getenv("HOME")) + "/images";
	}
	fullpath = fullpath + "/" + string(imagename);
	return fullpath.c_str();
}

void TestUtil::to_emobject(const Dict& d)
{
	if (d.has_key("farray")) {
		vector<float> array = d["farray"].get_farray();
		for (size_t i = 0; i < array.size(); i++) {
			printf("farray[%d] = %f\n", i, array[i]);
		}
	}
	
	if (d.has_key("emdata")) {
		EMData * img = d["emdata"];
		if (img) {
			printf("image size = (%d, %d, %d)\n",
				   img->get_xsize(), img->get_ysize(), img->get_zsize());
		}
	}
	
	if (d.has_key("int")) {
		int n = d["int"];
		printf("int n = %d\n", n);
	}
	
	if (d.has_key("float")) {
		float f = d["float"];
		printf("float f = %f\n", f);
	}
	
	if (d.has_key("long")) {
		int l = (int)d["long"];
		printf("long l = %d\n", l);
	}
	
}


void TestUtil::to_IntPoint(const IntPoint & p)
{
	printf("IntPoint p = (%d, %d, %d)\n", p[0], p[1], p[2]);
}


IntPoint TestUtil::from_IntPoint()
{
	return IntPoint(10, 20, 30);
}


void TestUtil::to_FloatPoint(const FloatPoint & p)
{
	printf("FloatPoint p = (%f, %f, %f)\n", p[0], p[1], p[2]);
}

FloatPoint TestUtil::from_FloatPoint()
{
	return FloatPoint(1.1, 2.2, 3.3);
}
	
void TestUtil::to_IntSize(const IntSize & p)
{
	printf("IntSize p = (%d, %d, %d)\n", p[0], p[1], p[2]);
}

IntSize TestUtil::from_IntSize()
{
	return IntSize(10, 20, 30);
}

void TestUtil::to_FloatSize(const FloatSize & p)
{
	printf("FloatSize p = (%f, %f, %f)\n", p[0], p[1], p[2]);
}

FloatSize TestUtil::from_FloatSize()
{
	return FloatSize(1.1, 2.2, 3.3);
}

void TestUtil::to_Vec3f(const Vec3f & p)
{
	printf("Vec3f p = (%f, %f, %f)\n", p[0], p[1], p[2]);
}

Vec3f TestUtil::from_Vec3f()
{
	return Vec3f(1.1, 2.2, 3.3);
}

void TestUtil::to_Vec3i(const Vec3i & p)
{
	printf("Vec3i p = (%d, %d, %d)\n", p[0], p[1], p[2]);
}

Vec3i TestUtil::from_Vec3i()
{
	return Vec3i(10, 20, 30);
}



