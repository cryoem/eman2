#ifndef __xydata_h__
#define __xydata_h__

#include <string>
using std::string;


namespace EMAN {

    class XYData {
    public:
	XYData();
	virtual ~XYData();

	int read_file(string filename);
	int write_file(string filename) const;
	
	void set_size(int new_size);
	void update();
	float calc_correlation(XYData* xy, float minx, float maxx) const;

	int write_x_data(float xdata[]);
	int write_y_data(float ydata[]);
	
	float* get_x_data();
	float* get_y_data();
	
	void get_x_data(float xdata[]);
	void get_y_data(float ydata[]);
	
	
	float get_yatx(float x) const;
	
	float get_x(int i) const;
	float get_y(int i) const;
	
	int get_size() const;
	
	float get_miny() const;
	float get_maxy() const;

	float get_minx() const;
	float get_maxx() const;
	
	bool is_validx(float x) const;
    };
}


#endif
