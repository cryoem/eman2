/**
 * $Id$
 */

/*
 * Author: Liwei Peng, 12/16/2004 (sludtke@bcm.edu)
 * Copyright (c) 2000-2006 Baylor College of Medicine
 *
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * */

#ifdef _WIN32
#include <windows.h>
	#define MAXPATHLEN MAX_PATH
#else
#include <sys/param.h>
#endif	//_WIN32

#include "testutil.h"
#include "xydata.h"

#include <algorithm>
#include "emassert.h"
//#include <stdlib.h>

using std::vector;
using std::string;
using std::map;

using namespace EMAN;

int TestUtil::ti[] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
float TestUtil::tf[] = {1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5};

const char * TestUtil::EMDATA_HEADER_EXT = ".head";
const char * TestUtil::EMDATA_DATA_EXT = ".data";

string TestUtil::progname = "";

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

Transform TestUtil::get_debug_transform(int i)
{
	vector<Transform> v(3);
	for (int j=0; j<3; j++) {
		Transform t;
		t.set_trans(j, j+1, j+2);
		v[j] = t;
	}
	return v[i];
}

string TestUtil::get_debug_image(const string & imagename)
{
	char imgpath[MAXPATHLEN];
	char * path_env = getenv("DEBUG_IMAGE_PATH");
	if (path_env) {
		sprintf(imgpath, "%s/%s", path_env, imagename.c_str());
	}
	else {
		sprintf(imgpath, "%s/images/%s", getenv("HOME"), imagename.c_str());
	}
	return string(imgpath);
}

string TestUtil::get_golden_image(const string & imagename)
{
	char imgpath[MAXPATHLEN];
	char * path_env = getenv("DEBUG_IMAGE_PATH");
	if (path_env) {
		sprintf(imgpath, "%s/testdata/%s", path_env, imagename.c_str());
	}
	else {
		sprintf(imgpath, "%s/images/testdata/%s", getenv("HOME"), imagename.c_str());
	}
	return string(imgpath);
}

void TestUtil::to_emobject(const Dict& d)
{
	if (d.has_key("floatarray")) {
		vector<float> array = d["floatarray"];
		for (size_t i = 0; i < array.size(); i++) {
			Assert(array[i] == tf[i]);
			LOGDEBUG("floatarray[%d] = %f\n", i, array[i]);
		}
	}

	if (d.has_key("emdata")) {
		EMData * img = d["emdata"];
		if (img) {
			int nx = img->get_xsize();
			int ny = img->get_ysize();
			int nz = img->get_zsize();
			Assert(nx == ti[0]);
			Assert(ny == ti[1]);
			Assert(nz == ti[2]);
			LOGDEBUG("image size = (%d, %d, %d)\n", nx, ny, nz);
		}
	}

	if (d.has_key("int")) {
		int n = d["int"];
		Assert(n == ti[0]);
		LOGDEBUG("int n = %d\n", n);
	}

	if (d.has_key("float")) {
		float f = d["float"];
		Assert(f == tf[0]);
		LOGDEBUG("float f = %f\n", f);
	}

	if (d.has_key("long")) {
		int l = (int)d["long"];
		Assert(l == ti[0]);
		LOGDEBUG("long l = %d\n", l);
	}

    if (d.has_key("string")) {
        string s = (const char*)d["string"];
        string s2 = get_debug_string(0);
        Assert(s == s2);
    }


	if (d.has_key("xydata")) {
		XYData *xyd = d["xydata"];
		size_t nitems = xyd->get_size();
		for (size_t i = 0; i < nitems; i++) {
			float xi = xyd->get_x(i);
			float yi = xyd->get_y(i);
			LOGDEBUG("xydata[%d] = (%f,%f)\n", i, xi, yi);
			Assert(xi == tf[i]);
			Assert(yi == tf[i]);
		}
	}

	if (d.has_key("stringarray")) {
		vector<string> array = d["stringarray"];
		for (size_t i = 0; i < array.size(); i++) {
			Assert(array[i] == get_debug_string(i));
			LOGDEBUG("stringarray[%d] = %s\n", i, array[i].c_str());
		}
	}

	if (d.has_key("transformarray")) {
		vector<Transform> array = d["transformarray"];
		for (size_t i = 0; i < array.size(); i++) {
//			array[i].printme();
			Assert(array[i] == get_debug_transform(i));
//			LOGDEBUG("transformarray[%d] = %s\n", i, array[i].to_str());
		}
	}
}


IntPoint TestUtil::test_IntPoint(const IntPoint & p)
{
	Assert(p[0] == ti[0]);
	Assert(p[1] == ti[1]);
	Assert(p[2] == ti[2]);
	LOGDEBUG("IntPoint p = (%d, %d, %d)\n", p[0], p[1], p[2]);
	return IntPoint(ti[0], ti[1], ti[2]);
}

FloatPoint TestUtil::test_FloatPoint(const FloatPoint & p)
{
	Assert(p[0] == tf[0]);
	Assert(p[1] == tf[1]);
	Assert(p[2] == tf[2]);
	LOGDEBUG("FloatPoint p = (%f, %f, %f)\n", p[0], p[1], p[2]);
	return FloatPoint(tf[0], tf[1], tf[2]);
}


IntSize TestUtil::test_IntSize(const IntSize & p)
{
	Assert(p[0] == ti[0]);
	Assert(p[1] == ti[1]);
	Assert(p[2] == ti[2]);
	LOGDEBUG("IntSize p = (%d, %d, %d)\n", p[0], p[1], p[2]);
	return IntSize(ti[0], ti[1], ti[2]);
}


FloatSize TestUtil::test_FloatSize(const FloatSize & p)
{
	Assert(p[0] == tf[0]);
	Assert(p[1] == tf[1]);
	Assert(p[2] == tf[2]);
	LOGDEBUG("FloatSize p = (%f, %f, %f)\n", p[0], p[1], p[2]);
	return FloatSize(tf[0], tf[1], tf[2]);
}


Vec3i TestUtil::test_Vec3i(const Vec3i & p)
{
	Assert(p[0] == ti[0]);
	Assert(p[1] == ti[1]);
	Assert(p[2] == ti[2]);
	LOGDEBUG("Vec3i p = (%d, %d, %d)\n", p[0], p[1], p[2]);
	return Vec3i(ti[0], ti[1], ti[2]);
}

Vec3f TestUtil::test_Vec3f(const Vec3f & p)
{
	Assert(p[0] == tf[0]);
	Assert(p[1] == tf[1]);
	Assert(p[2] == tf[2]);
	LOGDEBUG("Vec3f p = (%f, %f, %f)\n", p[0], p[1], p[2]);
	return Vec3f(tf[0], tf[1], tf[2]);
}


map<string, int> TestUtil::test_map_int(const map<string, int>& d)
{
	map<string, int> r;
	map<string, int>::const_iterator p;
	for (p = d.begin(); p != d.end(); p++) {
		LOGDEBUG("map[\"%s\"] = %d; ", p->first.c_str(), p->second);
		r[p->first] = p->second;
	}
	LOGDEBUG("\n");
	return r;
}

map<string, long> TestUtil::test_map_long(const map<string, long>& d)
{
	map<string, long> r;
	map<string, long>::const_iterator p;
	for (p = d.begin(); p != d.end(); p++) {
		LOGDEBUG("map[\"%s\"] = %d; ", p->first.c_str(), p->second);
		r[p->first] = p->second;
	}
	LOGDEBUG("\n");
	return r;
}

map<string, float> TestUtil::test_map_float(const map<string, float>& d)
{
	map<string, float> r;
	map<string, float>::const_iterator p;
	for (p = d.begin(); p != d.end(); p++) {
		LOGDEBUG("map[\"%s\"] = %f; ", p->first.c_str(), p->second);
		r[p->first] = p->second;
	}
	LOGDEBUG("\n");
	return r;
}

map<string, string> TestUtil::test_map_string(const map<string, string>& d)
{
	map<string, string> r;
	map<string, string>::const_iterator p;
	for (p = d.begin(); p != d.end(); p++) {
		LOGDEBUG("map[\"%s\"] = %s; ", p->first.c_str(), p->second.c_str());
		r[p->first] = p->second;
	}
	LOGDEBUG("\n");
	return r;
}

map<string, EMObject> TestUtil::test_map_emobject(const map<string, EMObject>& d)
{
	map<string, EMObject> r;
	map<string, EMObject>::const_iterator p;
	for (p = d.begin(); p != d.end(); p++) {
		LOGDEBUG("map[\"%s\"] = %f; ", p->first.c_str(), (float)(p->second));
		r[p->first] = EMObject(p->second);
	}
	LOGDEBUG("\n");
	return r;
}

map<string, vector<string> > TestUtil::test_map_vecstring(const map<string,
														  vector<string> >&)
{
	map<string, vector<string> > r;
	return r;
}


vector<int> TestUtil::test_vector_int(const vector<int> & v)
{
	vector<int> r;
	for (size_t i = 0; i < v.size(); i++) {
		LOGDEBUG("v[%d]=%d; ", i, v[i]);
		Assert(v[i] == ti[i]);
		r.push_back(v[i]);
	}
	LOGDEBUG("\n");
	return r;
}

vector<float> TestUtil::test_vector_float(const vector<float> & v)
{
	vector<float> r;
	for (size_t i = 0; i < v.size(); i++) {
		LOGDEBUG("v[%d]=%f; ", i, v[i]);
		Assert(v[i] == tf[i]);
		r.push_back(v[i]);
	}
	LOGDEBUG("\n");
	return r;
}

vector<long> TestUtil::test_vector_long(const vector<long> & v)
{
	vector<long> r;
	for (size_t i = 0; i < v.size(); i++) {
		LOGDEBUG("v[%d]=%d; ", i, (int)v[i]);
		Assert((int)v[i] == ti[i]);
		r.push_back(v[i]);
	}
	LOGDEBUG("\n");
	return r;
}

vector<string> TestUtil::test_vector_string(const vector<string> & v)
{
	vector<string> r;
	for (size_t i = 0; i < v.size(); i++) {
		LOGDEBUG("v[%d]=%s; ", i, v[i].c_str());
		r.push_back(v[i]);
	}
	LOGDEBUG("\n");
	return r;
}

vector<EMData*> TestUtil::test_vector_emdata(const vector<EMData*> & v)
{
	vector<EMData*> r;
	for (size_t i = 0; i < v.size(); i++) {
		EMData * e = v[i];
		LOGDEBUG("Image(%d,%d,%d); ", e->get_xsize(), e->get_ysize(), e->get_zsize());
		r.push_back(v[i]);
	}
	LOGDEBUG("\n");
	return r;
}

vector<Pixel> TestUtil::test_vector_pixel(const vector<Pixel> & v)
{
	vector<Pixel> r;
	for (size_t i = 0; i < v.size(); i++) {
		Pixel p = v[i];
		LOGDEBUG("Pixel(%d,%d,%d)=%4.2f; ", p.x, p.y, p.z, p.value);
		Pixel p2(p.x, p.y, p.z, p.value);
		r.push_back(p2);
	}

	return r;
}

Dict TestUtil::test_dict(const Dict & d)
{
	Dict r;

	vector<string> keys = d.keys();
	sort(keys.begin(), keys.end());

	for (size_t i = 0; i < keys.size(); i++) {
		LOGDEBUG("keys[%s] = %f\n", keys[i].c_str(), (float)d[keys[i]]);
		Assert(keys[i] == get_debug_string(i));
		Assert(((float)d[keys[i]]) == tf[i]);
		r[keys[i]] = d[keys[i]];
	}

	return r;
}


int TestUtil::check_image(const string& imagefile, EMData * image)
{
#if DEBUG
	string headerfile1 = Util::sbasename(imagefile) + EMDATA_HEADER_EXT;
	string datafile1 = Util::sbasename(imagefile) + EMDATA_DATA_EXT;

	char imgpath[MAXPATHLEN];
	char * path_env = getenv("DEBUG_IMAGE_PATH");
	if (path_env) {
		sprintf(imgpath, "%s/testdata/%s/", path_env, progname.c_str());
	}
	else {
		sprintf(imgpath, "%s/images/testdata/%s/", getenv("HOME"), progname.c_str());
	}

	string headerfile2 = string(imgpath) + headerfile1;
	string datafile2 = string(imgpath) + datafile1;


	if (image) {
		dump_emdata(image, imagefile);
	}
	else {
		dump_image_from_file(imagefile);
	}

    if (!Util::is_file_exist(headerfile2) ||
        !Util::is_file_exist(datafile2)) {
        return 0;
    }

	string diffcmd1 = "diff " + headerfile1 + " " + headerfile2;

	int err = system(diffcmd1.c_str());
	if (!err) {
		string diffcmd2 = "diff " + datafile1 + " " + datafile2;
		err = system(diffcmd2.c_str());
	}
	if (err) {
		LOGERR("check_image on %s FAILED\n", imagefile.c_str());
	}

	return err;
#endif
    return 0;
}

void TestUtil::dump_image_from_file(const string & filename)
{
	EMData * e = new EMData();
	e->read_image(filename);
	dump_emdata(e, filename);
	if( e )
	{
		delete e;
		e = 0;
	}
}

void TestUtil::dump_emdata(EMData * image, const string & filename)
{
	string filebase = Util::sbasename(filename);
	string headerfile = filebase + EMDATA_HEADER_EXT;
	string datafile = filebase + EMDATA_DATA_EXT;

	FILE *hfile = fopen(headerfile.c_str(), "wb");
	if (!hfile) {
		throw FileAccessException(headerfile);
	}
#if 0
	vector<string> excl_keys;
	excl_keys.push_back("MRC.label");
	excl_keys.push_back("IMAGIC.minute");
	excl_keys.push_back("IMAGIC.sec");

	Dict attr_dict = image->get_attr_dict();
	vector < string > keys = attr_dict.keys();



	for (size_t i = 0; i < keys.size(); i++) {

		bool is_exclude = false;
		for (size_t j = 0; j < excl_keys.size(); j++) {
			if (Util::sstrncmp(keys[i].c_str(), excl_keys[j].c_str())) {
				is_exclude = true;
				break;
			}
		}
		if (!is_exclude) {
			fprintf(hfile, "%s = %s\n", keys[i].c_str(),
					attr_dict[keys[i]].to_str().c_str());
		}
	}
#endif

	fprintf(hfile, "nx = %d\n", image->get_xsize());
	fprintf(hfile, "ny = %d\n", image->get_ysize());
	fprintf(hfile, "nz = %d\n", image->get_zsize());

	fclose(hfile);
	hfile = 0;

	FILE *dfile = fopen(datafile.c_str(), "wb");
	if (!dfile) {
		throw FileAccessException(datafile);
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	int nz = image->get_zsize();

	size_t row_size = nx * sizeof(float);
	size_t nxy = nx * ny;
	float * rdata = image->get_data();

	for (int i = 0; i < nz; i++) {
		for (int j = 0; j < ny; j++) {
			fwrite(&rdata[i * nxy + j * nx], row_size, 1, dfile);
		}
	}
	fclose(dfile);
	dfile = 0;
}

void TestUtil::set_progname(const string & cur_progname)
{
	progname = Util::sbasename(cur_progname);
}


void TestUtil::make_image_file_by_mode(const string & filename,
									   EMUtil::ImageType image_type, int mode,
									   EMUtil::EMDataType datatype,
									   int nx, int ny, int nz)
{
    EMData * e = new EMData();
    e->set_size(nx, ny, nz);
	bool is_complex = EMUtil::is_complex_type(datatype);

	e->set_attr("is_complex", (int)is_complex);
    e->set_attr("datatype", (int)datatype);
    float * data = e->get_data();

	size_t l = 0;
    for (int i = 0; i < nz; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nx; k++) {
                if (mode == 1) {
                    data[l] = get_pixel_value_by_dist1(nx, ny, nz, k, j, i);
                }
                else if (mode == 2) {
                    data[l] = get_pixel_value_by_dist2(nx, ny, nz, k, j, i);
                }
				l++;
            }
        }
    }

    if (!is_complex) {
        e->write_image(filename, 0, image_type, false, 0, datatype, true);
    }
    else {
    	e->update();
        e->set_attr("is_complex", false);
        EMData * fft = e->do_fft();
        fft->write_image(filename, 0, image_type, false, 0, datatype, true);
        if( fft )
        {
        	delete fft;
        	fft = 0;
        }
    }

	if( e )
	{
    	delete e;
    	e = 0;
	}
}

int TestUtil::verify_image_file_by_mode(const string & filename,
										EMUtil::ImageType, int mode,
										EMUtil::EMDataType datatype,
										int nx, int ny, int nz)
{
	int err = 0;

	EMData * e = new EMData();
	e->read_image(filename);

	Dict attr_dict = e->get_attr_dict();
	bool is_complex = EMUtil::is_complex_type(datatype);

	if (is_complex) {
		nx = (nx+2);
	}

	if (nx != (int) attr_dict["nx"]) {
        LOGERR("nx: %d != %d\n", nx, (int) attr_dict["nx"]);
        return 1;
    }

    if (ny != (int) attr_dict["ny"]) {
        LOGERR("ny: %d != %d\n", ny, (int) attr_dict["ny"]);
        return 1;
    }

    if (nz != (int) attr_dict["nz"]) {
        LOGERR("nz: %d != %d\n", nz, (int) attr_dict["nz"]);
        return 1;
    }

	if (datatype != (int) attr_dict["datatype"]) {
        LOGERR("datatype: %d != %d\n", datatype, (int) attr_dict["datatype"]);
        return 1;
    }


	if ((int)is_complex != (int) attr_dict["is_complex"]) {
        LOGERR("is_complex: %d != %d\n", is_complex, (int) attr_dict["is_complex"]);
        return 1;
    }


	if (!is_complex) {
		float * data = e->get_data();
		size_t l = 0;
		for (int i = 0; i < nz; i++) {
			for (int j = 0; j < ny; j++) {
				for (int k = 0; k < nx; k++) {

                    int d2 = 0;
                    if (mode == 1) {
                        d2 = (int)get_pixel_value_by_dist1(nx,ny,nz,k,j,i);
                    }
                    else if (mode == 2) {
                        d2 = (int)get_pixel_value_by_dist2(nx,ny,nz,k,j,i);
                    }

					if ((int)data[l] != d2) {
                        LOGERR("(%d,%d,%d): %d != %d\n", i,j,k,(int)data[l], d2);
						break;
						err = 1;
					}
					l++;
				}
			}
		}
	}

	return err;
}

EMObject TestUtil::emobject_to_py(bool b)
{
	return EMObject(b);
}

EMObject TestUtil::emobject_to_py(unsigned int un)
{
	return EMObject(un);
}

EMObject TestUtil::emobject_to_py(int n)
{
	return EMObject(n);
}

EMObject TestUtil::emobject_to_py(float f)
{
	return EMObject(f);
}

EMObject TestUtil::emobject_to_py(double f)
{
	return EMObject(f);
}

EMObject TestUtil::emobject_to_py(const string& str)
{
	return EMObject(str);
}


EMObject TestUtil::emobject_to_py(EMData * emdata)
{
	return EMObject(emdata);
}


EMObject TestUtil::emobject_to_py(XYData * xydata)
{
	return EMObject(xydata);
}


EMObject TestUtil::emobject_farray_to_py()
{
	vector<float> v(3);
	for (int i = 0; i < 3; i++) {
		v[i] = tf[i];
	}
	return EMObject(v);
}


EMObject TestUtil::emobject_strarray_to_py()
{
	vector<string> v(3);
	for (int i = 0; i < 3; i++) {
		v[i] = get_debug_string(i);
	}
	return EMObject(v);
}

EMObject TestUtil::emobject_transformarray_to_py()
{
	vector<Transform> v(3);
	for (int i=0; i<3; i++) {
		Transform t;
		t.set_trans(i, i+1, i+2);
		v[i] = t;
	}
	return EMObject(v);
}

EMObject TestUtil::emobject_to_py(Transform * t)
{
	return EMObject(t);
}

EMObject TestUtil::emobject_to_py(Ctf * ctf_)
{
	return EMObject(ctf_);
}




