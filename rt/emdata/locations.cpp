#include "emdata.h"
#include "log.h"
#include <assert.h>

using namespace EMAN;


const char* get_test_image()
{
	static char filename[256];
	static bool done = false;
	if (!done) {
		sprintf(filename, "%s/images/groel2d.mrc", getenv("HOME"));
		done = true;
	}
	return filename;
}

int main()
{
	int err = 0;
	
	EMData *image = new EMData();
	const char* imagefile = get_test_image();
	
	image->read_image(imagefile);

	float * data = image->get_data();

	int max_index = image->calc_max_index();
	Point<int> max_location = image->calc_max_location();

	int nx = image->get_xsize();
	int ny = image->get_ysize();
	//int nz = image->get_zsize();
	
	int max_index2 = max_location.z * nx * ny + max_location.y * nx + max_location.x;

	if (max_index2 != max_index) {
		LOGERR(" calc_max_index() and calc_max_location() have diff results: %d!=%d\n",
			   max_index, max_index2);
		err = 1;
	}
	
	LOGDEBUG("%s:  max index = %d\n", imagefile, max_index);
	LOGDEBUG("%s:  max location = (%d,%d,%d)\n",
			 imagefile, max_location.x,
			 max_location.y, max_location.z);

	float max = (float) image->get_attr("maximum");
	float max_loc_val = data[max_index];

	if (max != max_loc_val) {
		LOGERR("max_value: %f != max_location_value: %f\n",
			   max, max_loc_val);
		err = 1;
	}
	
	float threshold = 3.25f;
	const int num_high_pixels = 6;
	
	vector<Pixel> high_pixels = image->calc_highest_locations(threshold);
	int num_high_pixels2 = (int) high_pixels.size();

	if (num_high_pixels != num_high_pixels2) {
		LOGERR("# of highest pixels (value > %f) should be %d. You got %d\n",
			   threshold, num_high_pixels, num_high_pixels2);
		err = 1;
	}

	for (int i = 0; i < num_high_pixels2; i++) {
		LOGDEBUG("    (%d,%d,%d) = %f\n", high_pixels[i].x, high_pixels[i].y,
				 high_pixels[i].z, high_pixels[i].value);
	}

	
	return err;
}
