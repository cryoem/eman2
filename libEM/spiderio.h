/**
 * $Id$
 */
#ifndef eman__spiderio_h__
#define eman__spiderio_h__ 1


#include "imageio.h"
#include <stdio.h>

namespace EMAN
{
    class SpiderIO : public ImageIO
    {
    public:
	SpiderIO(string filename, IOMode rw_mode = READ_ONLY);
	~SpiderIO();

	DEFINE_IMAGEIO_FUNC;
	static bool is_valid(const void *first_block);
	
    protected:
	struct SpiderHeader
	{
	    float nslice;	// 1 for images
	    float ny;
	    float u1;		// unused
	    float u2;
	    float type;		// type of data 1=2d image,3=volume,-1=2d fft,-3=3d fft
	    float mmvalid;	// 1 if max and min are valid
	    float max;
	    float min;
	    float mean;
	    float sigma;	// std dev, -1=unknown
	    float u3;
	    float nx;
	    float headrec;	// # records in header
	    float angvalid;	// 1 if tilt angles have been computed
	    float phi;
	    float theta;
	    float gamma;
	    float dx;
	    float dy;
	    float dz;		// translations
	    float scale;
	    float headlen;	// in bytes
	    float reclen;
	    float istack;	// file contains a stack of images (2 in overall header, -1 in images)
	    float inuse;	// indicates that the current image is present in the stack
	    float maxim;
	    float imgnum;
	    float u6;
	    float u7;
	    float ang2valid;	// additional angles
	    float phi1;
	    float theta1;
	    float psi1;
	    float phi2;
	    float theta2;
	    float psi2;
	    float xf[14];	// Jose Maria's transforms
	    float u8[161];
	    char date[12];
	    char time[8];
	    char name[160];
	};

	enum SpiderType {
	    IMAGE_2D = 1,
	    IMAGE_3D = 3,
	    IMAGE_2D_FFT = -1,
	    IMAGE_3D_FFT = -3
	};

	enum {
	    SINGLE_IMAGE_HEADER = 0,
	    STACK_IMAGE_HEADER = -1,
	    STACK_OVERALL_HEADER = 2,
	    NUM_FLOATS_IN_HEADER = 211
	};

	int write_single_header(const Dict & dict);
	int write_single_data(float *data);
	virtual bool is_valid_spider(const void *first_block);

    private:
	bool need_swap() const;
	void swap_data(float *data, int nitems);
	void swap_header(SpiderHeader * header);

    private:
	string filename;
	IOMode rw_mode;

	FILE *spider_file;
	SpiderHeader *first_h;
	SpiderHeader *cur_h;

	bool is_big_endian;
	bool initialized;
	bool is_new_file;

    };

}

#endif
