#ifndef eman__xydata_h__
#define eman__xydata_h__ 1

#include <string>
#include <vector>

using std::string;
using std::vector;

namespace EMAN
{
    class XYData
    {
    public:
	struct Pair
	{
	    Pair(float xx, float yy)
		: x(xx), y(yy)
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
	virtual ~XYData()
	{
	}

	int read_file(string filename);

	int write_file(string filename) const;

	float calc_correlation(XYData * xy, float minx, float maxx) const;

	float get_yatx(float x) const;

	float get_x(size_t i) const
	{
	    if (i < 0 || i > data.size()) {
		return 0;
	    }
	    return data[i].x;
	}

	float get_y(size_t i) const
	{
	    if (i < 0 || i > data.size()) {
		return 0;
	    }
	    return data[i].y;
	}

	int get_size() const
	{
	    return (int)data.size();
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
	    if (x < data[0].x || x > data[data.size() - 1].x) {
		return false;
	    }
	    return true;
	}

    private:
	vector<Pair> data;
	float ymin;
	float ymax;
	float mean_x_spacing;

	void update();
    };
}


#endif
