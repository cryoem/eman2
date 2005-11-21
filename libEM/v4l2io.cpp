/**
 * $Id$
 */
#ifdef ENABLE_V4L2
#include "v4l2io.h"
#include "geometry.h"
#include "util.h"
#include "portable_fileio.h"

#include <limits.h>
#include <ctype.h>
#include <sys/ioctl.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

using namespace EMAN;

static int xioctl (int v4l_file, int request, void *arg)
{
        int r;

        do r = ioctl (v4l_file, request, arg);
        while (-1 == r && EINTR == errno);

        return r;
}

void errno_exit (const char *s)
{
        fprintf (stderr, "%s error %d, %s\n",s, errno, strerror (errno));
        exit (EXIT_FAILURE);
}


V4L2IO::V4L2IO(const string & file, IOMode rw)
: initialized(false)
{
	filename=(char *)malloc(file.length()+1);
	strcpy(filename,file.c_str());
	v4l_file=0;
	nx = 0;
	ny = 0;

}

V4L2IO::~V4L2IO()
{
	if (v4l_file) {
		close(v4l_file);
		v4l_file = 0;
	}
	if (filename) free(filename);
}

void V4L2IO::init()
{
	ENTERFUNC;
	
	if (initialized) {
		return;
	}

	initialized = true;

	v4l_file = open (filename, O_RDWR | O_NONBLOCK, 0);

	unsigned int min;
	
	if (-1 == xioctl (v4l_file, VIDIOC_QUERYCAP, &cap)) {
			if (EINVAL == errno) {
					fprintf (stderr, "%s is not a V4L2 device, try /dev/vbi*\n",filename);
					exit (-1);
			} else {
					errno_exit ("VIDIOC_QUERYCAP");
			}
	}
	
	printf("driver: %s\ncard %s\n",cap.driver,cap.card);
	
	int input=1;
	if (-1 == xioctl (v4l_file, VIDIOC_S_INPUT, &input)) errno_exit ("VIDIOC_S_INPUT");
	
	int std=V4L2_STD_NTSC_M;
	if (-1 == xioctl (v4l_file, VIDIOC_S_STD, &std)) errno_exit ("VIDIOC_S_STD");
	
	
	if (!(cap.capabilities & V4L2_CAP_VIDEO_CAPTURE)) {
			fprintf (stderr, "%s is not a video capture device\n",filename);
			exit (EXIT_FAILURE);
	}
	
	if (!(cap.capabilities & V4L2_CAP_READWRITE)) {
			fprintf (stderr, "%s does not support read i/o\n",filename);
			exit (EXIT_FAILURE);
	}
	
/*	cropcap.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
	if (-1 == xioctl (v4l_file, VIDIOC_CROPCAP, &cropcap)) {
		fprintf(stderr,"VIDIOC_CROPCAP failed %d %d %d %d\n",cropcap.bounds.left,cropcap.bounds.top,cropcap.bounds.width,cropcap.bounds.height);
//		exit(EXIT_FAILURE);
	}
	
	crop.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
	crop.c = cropcap.defrect;
	
	if (-1 == xioctl (v4l_file, VIDIOC_S_CROP, &crop)) {
			switch (errno) {
			case EINVAL:
					//Cropping not supported.
					break;
			default:
					// Errors ignored.
					break;
			}
	}
*/

	printf("Supported formats:\n");
	struct v4l2_fmtdesc fdq;
	fdq.type=V4L2_BUF_TYPE_VIDEO_CAPTURE;
	for (int idx=0; idx<100; idx++) {
		fdq.index=idx;
		if (xioctl(v4l_file,VIDIOC_ENUM_FMT,&fdq)==-1) break;
		printf("%4s %s\n",&fdq.pixelformat,fdq.description);
	}

//	memset(&fmt,0,sizeof(fmt));
	fmt.type=V4L2_BUF_TYPE_VIDEO_CAPTURE;
	if (-1 == xioctl (v4l_file, VIDIOC_G_FMT, &fmt)) errno_exit ("VIDIOC_G_FMT");
	
	fmt.type                = V4L2_BUF_TYPE_VIDEO_CAPTURE;
//	fmt.fmt.pix.width       = 640; 
//	fmt.fmt.pix.height      = 480;
	fmt.fmt.pix.pixelformat = V4L2_PIX_FMT_GREY;
//	fmt.fmt.pix.pixelformat = V4L2_PIX_FMT_YUYV;
//	fmt.fmt.pix.field       = V4L2_FIELD_INTERLACED;
	
	if (-1 == xioctl (v4l_file, VIDIOC_S_FMT, &fmt)) errno_exit ("VIDIOC_S_FMT");
	
	/* Note VIDIOC_S_FMT may change width and height. */
	
	/* Buggy driver paranoia. */
	min = fmt.fmt.pix.width * 2;
	if (fmt.fmt.pix.bytesperline < min)
			fmt.fmt.pix.bytesperline = min;
	min = fmt.fmt.pix.bytesperline * fmt.fmt.pix.height;
	if (fmt.fmt.pix.sizeimage < min)
			fmt.fmt.pix.sizeimage = min;
	
	printf("fmt.fmt.pix.width = %d\n",fmt.fmt.pix.width);
	printf("fmt.fmt.pix.height = %d\n",fmt.fmt.pix.height);
	printf("fmt.fmt.pix.pixelformat = %4s\n",&fmt.fmt.pix.pixelformat);
	printf("fmt.fmt.pix.bytesperline = %d\n",fmt.fmt.pix.bytesperline);
	printf("fmt.fmt.pix.sizeimage = %d\n",fmt.fmt.pix.sizeimage);
			
	EXITFUNC;
}

bool V4L2IO::is_valid(const void *first_block)
{
	ENTERFUNC;
	return true;
}

int V4L2IO::read_header(Dict & dict, int image_index, const Region * area, bool)
{
	ENTERFUNC;
	if (!initialized) init();
	
	dict["nx"] = (int)(fmt.fmt.pix.width);
	dict["ny"] = (int)(fmt.fmt.pix.height);
	dict["nz"] = 1;
			
	dict["datatype"] = EMUtil::EM_UCHAR;
	
	EXITFUNC;
	return 0;
}

int V4L2IO::write_header(const Dict & dict, int image_index, const Region*,
						EMUtil::EMDataType, bool)
{
	ENTERFUNC;	
	EXITFUNC;
	return 1;	// No write capability
}

int V4L2IO::read_data(float *data, int image_index, const Region * area, bool)
{
	if (!initialized) init();
	
	Region tmp(0,0,fmt.fmt.pix.width,fmt.fmt.pix.height);
	if (!area) area=&tmp;
	
	int x,y;
	ENTERFUNC;
	unsigned char *dbuf = (unsigned char *)malloc(fmt.fmt.pix.sizeimage);
	read(v4l_file,dbuf,fmt.fmt.pix.sizeimage);
	
	printf("buf.length=%d  area.x=%d area.y=%d\n",fmt.fmt.pix.sizeimage,(int)area->size[0],(int)area->size[1]);
	
	for (y=0; y<area->size[1]; y++) {
		for (x=0; x<area->size[0]; x++) {
			data[(int)area->size[0]*y+x]=dbuf[fmt.fmt.pix.bytesperline*y+x+1]/256.0;
		}
	}

	int sum=0;
	for (x=0; x<fmt.fmt.pix.sizeimage; x++) sum+=dbuf[x];
	printf("%d\n",sum);
	
	free(dbuf);
	EXITFUNC;
	return 0;
}

int V4L2IO::write_data(float *data, int image_index, const Region* ,
					  EMUtil::EMDataType, bool)
{
	ENTERFUNC;

	EXITFUNC;
	return 1;
}

void V4L2IO::flush()
{

}


bool V4L2IO::is_complex_mode()
{
	return false;
}

bool V4L2IO::is_image_big_endian()
{
	return 1;
}

#endif
