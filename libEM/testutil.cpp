#include "testutil.h"
#include "emdata.h"
#include "emobject.h"

#include <vector>
#include <string>
#include <map>
#include <assert.h>

using std::vector;
using std::string;
using std::map;

using namespace EMAN;

int TestUtil::ti[] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
float TestUtil::tf[] = {1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5};
const char *TestUtil::ts[] = {"aa", "bb", "cc", "dd", "ee", "ff", "gg", "hh", "ii", "jj"};


int TestUtil::get_debug_int(int i)
{
	return ti[i];
}

float TestUtil::get_debug_float(int i)
{
	return tf[i];
}

const char* TestUtil::get_debug_string(int i)
{
	return ts[i];
}

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


IntPoint TestUtil::test_IntPoint(const IntPoint & p)
{
	assert(p[0] == ti[0]);
	assert(p[1] == ti[1]);
	assert(p[2] == ti[2]);
	printf("IntPoint p = (%d, %d, %d)\n", p[0], p[1], p[2]);
	return IntPoint(ti[0], ti[1], ti[2]);
}

FloatPoint TestUtil::test_FloatPoint(const FloatPoint & p)
{
	assert(p[0] == tf[0]);
	assert(p[1] == tf[1]);
	assert(p[2] == tf[2]);
	printf("FloatPoint p = (%f, %f, %f)\n", p[0], p[1], p[2]);
	return FloatPoint(tf[0], tf[1], tf[2]);
}

	
IntSize TestUtil::test_IntSize(const IntSize & p)
{
	assert(p[0] == ti[0]);
	assert(p[1] == ti[1]);
	assert(p[2] == ti[2]);
	printf("IntSize p = (%d, %d, %d)\n", p[0], p[1], p[2]);
	return IntSize(ti[0], ti[1], ti[2]);
}


FloatSize TestUtil::test_FloatSize(const FloatSize & p)
{
	assert(p[0] == tf[0]);
	assert(p[1] == tf[1]);
	assert(p[2] == tf[2]);
	printf("FloatSize p = (%f, %f, %f)\n", p[0], p[1], p[2]);
	return FloatSize(tf[0], tf[1], tf[2]);
}


Vec3i TestUtil::test_Vec3i(const Vec3i & p)
{
	assert(p[0] == ti[0]);
	assert(p[1] == ti[1]);
	assert(p[2] == ti[2]);
	printf("Vec3i p = (%d, %d, %d)\n", p[0], p[1], p[2]);
	return Vec3i(ti[0], ti[1], ti[2]);
}

Vec3f TestUtil::test_Vec3f(const Vec3f & p)
{
	assert(p[0] == tf[0]);
	assert(p[1] == tf[1]);
	assert(p[2] == tf[2]);
	printf("Vec3f p = (%f, %f, %f)\n", p[0], p[1], p[2]);
	return Vec3f(tf[0], tf[1], tf[2]);
}
