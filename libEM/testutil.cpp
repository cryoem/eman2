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

int TestUtil::get_debug_int(int i)
{
	return ti[i];
}

float TestUtil::get_debug_float(int i)
{
	return tf[i];
}

string TestUtil::get_debug_string(int i)
{
	char a[32];
	sprintf(a, "%d%d", i+1, i+1);
	return string(a);
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

#if 1
map<string, int> TestUtil::test_map_int(const map<string, int>& d)
{
	map<string, int> r;
	map<string, int>::const_iterator p;
	for (p = d.begin(); p != d.end(); p++) {
		printf("map[\"%s\"] = %d; ", p->first.c_str(), p->second);
		r[p->first] = p->second;
	}
	printf("\n");
	return r;
}

map<string, long> TestUtil::test_map_long(const map<string, long>& d)
{
	map<string, long> r;
	map<string, long>::const_iterator p;
	for (p = d.begin(); p != d.end(); p++) {
		printf("map[\"%s\"] = %d; ", p->first.c_str(), p->second);
		r[p->first] = p->second;
	}
	printf("\n");
	return r;
}

map<string, float> TestUtil::test_map_float(const map<string, float>& d)
{
	map<string, float> r;
	map<string, float>::const_iterator p;
	for (p = d.begin(); p != d.end(); p++) {
		printf("map[\"%s\"] = %f; ", p->first.c_str(), p->second);
		r[p->first] = p->second;
	}
	printf("\n");
	return r;
}

map<string, string> TestUtil::test_map_string(const map<string, string>& d)
{
	map<string, string> r;
	map<string, string>::const_iterator p;
	for (p = d.begin(); p != d.end(); p++) {
		printf("map[\"%s\"] = %s; ", p->first.c_str(), p->second.c_str());
		r[p->first] = p->second;
	}
	printf("\n");
	return r;
}

map<string, EMObject> TestUtil::test_map_emobject(const map<string, EMObject>& d)
{
	map<string, EMObject> r;
	map<string, EMObject>::const_iterator p;
	for (p = d.begin(); p != d.end(); p++) {
		printf("map[\"%s\"] = %f; ", p->first.c_str(), (float)(p->second));
		r[p->first] = EMObject((float)p->second);
	}
	printf("\n");
	return r;
}

map<string, vector<string> > TestUtil::test_map_vecstring(const map<string,
														  vector<string> >& d)
{
	map<string, vector<string> > r;
	return r;
}

#endif

vector<int> TestUtil::test_vector_int(const vector<int> & v)
{
	vector<int> r;
	for (size_t i = 0; i < v.size(); i++) {
		printf("v[%d]=%d; ", i, v[i]);
		assert(v[i] == ti[i]);
		r.push_back(v[i]);
	}
	printf("\n");
	return r;
}

vector<float> TestUtil::test_vector_float(const vector<float> & v)
{
	vector<float> r;
	for (size_t i = 0; i < v.size(); i++) {
		printf("v[%d]=%f; ", i, v[i]);
		assert(v[i] == tf[i]);
		r.push_back(v[i]);
	}
	printf("\n");
	return r;
}

vector<long> TestUtil::test_vector_long(const vector<long> & v)
{
	vector<long> r;
	for (size_t i = 0; i < v.size(); i++) {
		printf("v[%d]=%d; ", i, (int)v[i]);
		assert((int)v[i] == ti[i]);
		r.push_back(v[i]);
	}
	printf("\n");
	return r;
}

vector<string> TestUtil::test_vector_string(const vector<string> & v)
{
	vector<string> r;
	for (size_t i = 0; i < v.size(); i++) {
		printf("v[%d]=%s; ", i, v[i].c_str());
		r.push_back(v[i]);
	}
	printf("\n");
	return r;
}

vector<EMData*> TestUtil::test_vector_emdata(const vector<EMData*> & v)
{
	vector<EMData*> r;
	for (size_t i = 0; i < v.size(); i++) {
		EMData * e = v[i];
		printf("Image(%d,%d,%d);", e->get_xsize(), e->get_ysize(), e->get_zsize());
		r.push_back(v[i]);
	}
	printf("\n");
	return r;
}

vector<Pixel> TestUtil::test_vector_pixel(const vector<Pixel> & v)
{
	vector<Pixel> r;
	return r;
}

#if 0
void  TestUtil::test_map_int(const map<string, int>& d)
{


}

void TestUtil::test_map_long(const map<string, long>& d)
{


}

void  TestUtil::test_map_float(const map<string, float>& d)
{

}

void  TestUtil::test_map_string(const map<string, string>& d)
{
}
#endif
