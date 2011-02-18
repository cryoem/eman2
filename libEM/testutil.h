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

#ifndef eman__testutil_h__
#define eman__testutil_h__

#include "emdata.h"

using std::map;
using std::vector;
using std::string;

namespace EMAN
{
	/**TestUtil defines function assisting testing of EMAN2.
	 */

	class Dict;

	class TestUtil
	{
	public:
		static const char * EMDATA_HEADER_EXT;
		static const char * EMDATA_DATA_EXT;

		static int get_debug_int(int i);
		static float get_debug_float(int i);
		static string get_debug_string(int i);
		static Transform get_debug_transform(int i);

		static string get_debug_image(const string & imagename);
		static string get_golden_image(const string & imagename);

		static void to_emobject(const Dict & d);

		static EMObject emobject_to_py(bool b);
		static EMObject emobject_to_py(unsigned int n);
		static EMObject emobject_to_py(int n);
		static EMObject emobject_to_py(float f);
		static EMObject emobject_to_py(double f);
		static EMObject emobject_to_py(const string& str);
		static EMObject emobject_to_py(EMData * emdata);
		static EMObject emobject_to_py(XYData * xydata);
		static EMObject emobject_farray_to_py();
		static EMObject emobject_strarray_to_py();
		static EMObject emobject_transformarray_to_py();
		static EMObject emobject_to_py(Transform * t);
		static EMObject emobject_to_py(Ctf * ctf_);

		static IntPoint test_IntPoint(const IntPoint & p);
		static FloatPoint test_FloatPoint(const FloatPoint & p);
		static IntSize test_IntSize(const IntSize & p);
		static FloatSize test_FloatSize(const FloatSize & p);
		static Vec3i test_Vec3i(const Vec3i & p);
		static Vec3f test_Vec3f(const Vec3f & p);

		static vector<int> test_vector_int(const vector<int> & v);
		static vector<float> test_vector_float(const vector<float> & v);
		static vector<long> test_vector_long(const vector<long> & v);
		static vector<string> test_vector_string(const vector<string> & v);
		static vector<EMData*> test_vector_emdata(const vector<EMData*> & v);
		static vector<Pixel> test_vector_pixel(const vector<Pixel> & v);

		static map<string, int> test_map_int(const map<string, int>& d);
		static map<string, long> test_map_long(const map<string, long>& d);
		static map<string, float> test_map_float(const map<string, float>& d);
		static map<string, string> test_map_string(const map<string, string>& d);
		static map<string, EMObject> test_map_emobject(const map<string, EMObject>& d);
		static map<string, vector<string> > test_map_vecstring(const map<string,
															   vector<string> >& d);

		static Dict test_dict(const Dict & d);

		static void dump_image_from_file(const string & filename);
		static void dump_emdata(EMData * image, const string & filename);
		static int check_image(const string& imagefile, EMData * image = 0);
        static void set_progname(const string & cur_progname);

        static void make_image_file(const string & filename,
									EMUtil::ImageType image_type,
									EMUtil::EMDataType datatype = EMUtil::EM_FLOAT,
									int nx = 16, int ny = 16, int nz = 1)
        {
            make_image_file_by_mode(filename, image_type, 1, datatype, nx, ny, nz);
        }

		static int verify_image_file(const string & filename,
									 EMUtil::ImageType image_type,
									 EMUtil::EMDataType datatype = EMUtil::EM_FLOAT,
									 int nx = 16, int ny = 16, int nz = 1)
        {
            return verify_image_file_by_mode(filename, image_type, 1, datatype, nx, ny, nz);
        }

        static void make_image_file2(const string & filename,
									 EMUtil::ImageType image_type,
									 EMUtil::EMDataType datatype = EMUtil::EM_FLOAT,
									 int nx = 16, int ny = 16, int nz = 1)
        {
            make_image_file_by_mode(filename, image_type, 2, datatype,nx, ny, nz);
        }

		static int verify_image_file2(const string & filename,
									  EMUtil::ImageType image_type,
									  EMUtil::EMDataType datatype = EMUtil::EM_FLOAT,
									  int nx = 16, int ny = 16, int nz = 1)
        {
            return verify_image_file_by_mode(filename, image_type, 2,
											 datatype, nx, ny, nz);
        }




	private:
		static float tf[10];
		static int ti[10];
		static string progname;

        static void make_image_file_by_mode(const string & filename,
											EMUtil::ImageType image_type, int mode,
											EMUtil::EMDataType datatype = EMUtil::EM_FLOAT,
											int nx = 16, int ny = 16, int nz = 1);

		static int verify_image_file_by_mode(const string & filename,
											 EMUtil::ImageType image_type, int mode,
											 EMUtil::EMDataType datatype = EMUtil::EM_FLOAT,
											 int nx = 16, int ny = 16, int nz = 1);


		static float get_pixel_value_by_dist1(int nx, int ny, int nz, int x, int y, int z)
		{
            int x2 = x;
            int y2 = y;
            int z2 = z;

            x2 = abs(nx/2-x);
            y2 = abs(ny/2-y);

            if (z > nz/2) {
                z2 = nz-z;
            }

            if (nz == 1) {
                return (float)(x2*x2 + y2*y2);
            }
            else {
                int areax = (int)((float)nx * z2 / nz);
                int areay = (int)((float)ny * z2 / nz);
                if ((abs(x-nx/2) <= areax) && (abs(y-ny/2) <= areay)) {
                    return (float)(x2*x2 + y2*y2);
                }
                else {
                    return 0;
                }
            }
		}

        static float get_pixel_value_by_dist2(int nx, int ny, int nz, int x, int y, int z)
		{
            int x2 = x;
            int y2 = y;
            int z2 = z;


            if (x > nx/2) {
                x2 = nx-x;
            }
            if (y > ny/2) {
                y2 = ny-y;
            }

            if (z > nz/2) {
                z2 = nz-z;
            }

            if (nz == 1) {
                return (float)(x2*x2 + y2*y2);
            }
            else {
                int areax = (int)((float)nx * z2 / nz);
                int areay = (int)((float)ny * z2 / nz);
                if ((abs(x-nx/2) <= areax) && (abs(y-ny/2) <= areay)) {
                    return (float)(x2*x2 + y2*y2);
                }
                else {
                    return 0;
                }
            }
        }
	};





}

#endif	//eman__testutil_h__


