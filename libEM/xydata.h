/**
 * $Id$
 */

/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
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

#ifndef eman__xydata_h__
#define eman__xydata_h__ 1

#include <string>
#include <vector>

using std::string;
using std::vector;

namespace EMAN
{
	/** XYData defines a 1D (x,y) data set.
     */
	class XYData
	{
	  public:
		struct Pair
		{
			Pair(float xx, float yy)
			:x(xx), y(yy)
			{
			}

			bool operator<(const Pair & p) const
			{
				return (x < p.x);
			}

			float x;
			float y;
		};

	  public:
		XYData();
		virtual ~ XYData()
		{
		}

		int read_file(const string & filename);

		int write_file(const string & filename) const;

		float calc_correlation(XYData * xy, float minx, float maxx) const;

		void update();

		float get_yatx(float x);

		float get_x(size_t i) const
		{
			return data[i].x;
		}

		void set_x(size_t i, float x)
		{
			if (i>=data.size()) data.resize(i+1,Pair(0,0));
			data[i].x = x;
		}

		float get_y(size_t i) const
		{
			return data[i].y;
		}

		void set_y(size_t i, float y)
		{
			if (i>=data.size()) data.resize(i+1,Pair(0,0));
			data[i].y = y;
		}

		vector<float> get_xlist() const;

		vector<float> get_ylist() const;

		void set_xy_list(const vector<float>& xlist, const vector<float>& ylist);

		size_t get_size() const
		{
			return data.size();
		}

		void set_size(size_t n);

		float get_miny()
		{
			update();
			return ymin;
		}

		float get_maxy()
		{
			update();
			return ymax;
		}

		bool is_validx(float x) const
		{
			if (x < data[0].x || x > data[data.size() - 1].x)
			{
				return false;
			}
			return true;
		}

	  private:
		vector < Pair > data;
		float ymin;
		float ymax;
		float mean_x_spacing;
	};
}


#endif
