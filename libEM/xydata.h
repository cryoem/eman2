/**
 * $Id$
 */
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

		float get_yatx(float x) const;

		float get_x(size_t i) const
		{
			return data[i].x;
		}

		void set_x(size_t i, float x)
		{
			data[i].x = x;
		}

		float get_y(size_t i) const
		{
			return data[i].y;
		}

		void set_y(size_t i, float y)
		{
			data[i].y = y;
		}

		size_t get_size() const
		{
			return data.size();
		}

		float get_miny() const
		{
			return ymin;
		}

		float get_maxy() const
		{
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
