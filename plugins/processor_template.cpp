#include "processor_template.h"

using namespace EMAN;


/** define your Processor operation
 */
void XYZProcessor::process(EMData * image)
{
	if (!image) {
		return;
	}


	// The following are the sample code to get your parameters. Then
	// go through the image data pixel.
#if 0
	int value1 = params["value1"];
	float value2 = params["value2"];

	float *data = image->get_data();
	int size = image->get_xsize() * image->get_ysize() * image->get_zsize();
	for (int i = 0; i < size; i++) {
		if (data[i] <= value1 && data[i] >= value2) {
			data[i] = 0;
		}
	}
#endif

}
