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
struct v4l2_format v4l2_fmt;	// should be a static member, but doesn't seem to work right

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
{
	filename=(char *)malloc(file.length()+1);
	strcpy(filename,file.c_str());
	v4l_file=0;
	nx = 0;
	ny = 0;
	initialized=false;

}

V4L2IO::~V4L2IO()
{
	if (v4l_file) {
		close(v4l_file);
		v4l_file = 0;
	}
	if (filename) free(filename);
}

void V4L2IO::init() {
	ENTERFUNC;
	static int ginit=0;
	
	if (!ginit) globalinit(filename);
	ginit=1;
	
	if (initialized) {
		return;
	}

	initialized = true;

	v4l_file = open (filename, O_RDWR, 0);

	EXITFUNC;
}

void V4L2IO::globalinit(const char *fsp)
{
	ENTERFUNC;
	
//	v4l_file = open (filename, O_RDWR | O_NONBLOCK, 0);
	int vfile = open (fsp, O_RDWR, 0);

	unsigned int min;
	struct v4l2_capability cap;
	struct v4l2_cropcap cropcap;
	struct v4l2_crop crop;
	
	if (-1 == xioctl (vfile, VIDIOC_QUERYCAP, &cap)) {
			if (EINVAL == errno) {
					fprintf (stderr, "%s is not a V4L2 device, try /dev/vbi*\n",fsp);
					exit (-1);
			} else {
					errno_exit ("VIDIOC_QUERYCAP");
			}
	}
	
	printf("driver: %s\ncard %s\n",cap.driver,cap.card);
	
	int input=0;
	if (-1 == xioctl (vfile, VIDIOC_S_INPUT, &input)) errno_exit ("VIDIOC_S_INPUT");
	
	int std=V4L2_STD_NTSC_M;
	if (-1 == xioctl (vfile, VIDIOC_S_STD, &std)) errno_exit ("VIDIOC_S_STD");
	
	
	if (!(cap.capabilities & V4L2_CAP_VIDEO_CAPTURE)) {
			fprintf (stderr, "%s is not a video capture device\n",fsp);
			exit (EXIT_FAILURE);
	}
	
	if (!(cap.capabilities & V4L2_CAP_READWRITE)) {
			fprintf (stderr, "%s does not support read i/o\n",fsp);
			exit (EXIT_FAILURE);
	}
	
/*	cropcap.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
	if (-1 == xioctl (vfile, VIDIOC_CROPCAP, &cropcap)) {
		fprintf(stderr,"VIDIOC_CROPCAP failed %d %d %d %d\n",cropcap.bounds.left,cropcap.bounds.top,cropcap.bounds.width,cropcap.bounds.height);
//		exit(EXIT_FAILURE);
	}
	
	crop.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
	crop.c = cropcap.defrect;
	
	if (-1 == xioctl (vfile, VIDIOC_S_CROP, &crop)) {
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

	/*
	printf("Supported formats:\n");
	struct v4l2_fmtdesc fdq;
	fdq.type=V4L2_BUF_TYPE_VIDEO_CAPTURE;
	for (int idx=0; idx<100; idx++) {
		fdq.index=idx;
		if (xioctl(vfile,VIDIOC_ENUM_FMT,&fdq)==-1) break;
		printf("%4s %s\n",&fdq.pixelformat,fdq.description);
	}
	*/

//	memset(&v4l2_fmt,0,sizeof(v4l2_fmt));
	v4l2_fmt.type=V4L2_BUF_TYPE_VIDEO_CAPTURE;
	if (-1 == xioctl (vfile, VIDIOC_G_FMT, &v4l2_fmt)) errno_exit ("VIDIOC_G_FMT");
	
	v4l2_fmt.type                = V4L2_BUF_TYPE_VIDEO_CAPTURE;
//	v4l2_fmt.fmt.pix.width       = 640; 
//	v4l2_fmt.fmt.pix.height      = 480;
//	v4l2_fmt.fmt.pix.pixelformat = V4L2_PIX_FMT_GREY;
	v4l2_fmt.fmt.pix.pixelformat = V4L2_PIX_FMT_YUYV;
	v4l2_fmt.fmt.pix.field       = V4L2_FIELD_NONE;
//	v4l2_fmt.fmt.pix.field       = V4L2_FIELD_INTERLACED;
	
	if (-1 == xioctl (vfile, VIDIOC_S_FMT, &v4l2_fmt)) errno_exit ("VIDIOC_S_FMT");
	
	/* Note VIDIOC_S_FMT may change width and height. */
	
	/* Buggy driver paranoia. */
	min = v4l2_fmt.fmt.pix.width * 2;
	if (v4l2_fmt.fmt.pix.bytesperline < min)
			v4l2_fmt.fmt.pix.bytesperline = min;
	min = v4l2_fmt.fmt.pix.bytesperline * v4l2_fmt.fmt.pix.height;
	if (v4l2_fmt.fmt.pix.sizeimage < min)
			v4l2_fmt.fmt.pix.sizeimage = min;
	
/*	printf("fmt.fmt.pix.width = %d\n",v4l2_fmt.fmt.pix.width);
	printf("fmt.fmt.pix.height = %d\n",v4l2_fmt.fmt.pix.height);
	printf("fmt.fmt.pix.pixelformat = %4s\n",&v4l2_fmt.fmt.pix.pixelformat);
	printf("fmt.fmt.pix.bytesperline = %d\n",v4l2_fmt.fmt.pix.bytesperline);
	printf("fmt.fmt.pix.sizeimage = %d\n",v4l2_fmt.fmt.pix.sizeimage);
	printf("fmt.fmt.pix.field = %d\n",v4l2_fmt.fmt.pix.field); */
	
	close(vfile);
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
	
	dict["nx"] = (int)(v4l2_fmt.fmt.pix.width);
	dict["ny"] = (int)(v4l2_fmt.fmt.pix.height);
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
	
	Region tmp(0,0,v4l2_fmt.fmt.pix.width,v4l2_fmt.fmt.pix.height);
	if (!area) area=&tmp;
	
	int x,y;
	ENTERFUNC;
	unsigned char *dbuf = (unsigned char *)malloc(v4l2_fmt.fmt.pix.sizeimage);

	read(v4l_file,dbuf,v4l2_fmt.fmt.pix.sizeimage);

/*	int count=v4l2_fmt.fmt.pix.sizeimage;
	int pos=0;
	while (count) {
		int c=read(v4l_file,dbuf+pos,count);
		printf("# %d\n",c);
		if (c==-1) { printf("ern %d\n",errno); break; }
		count-=c;
		pos+=c;
	}
	
	printf("buf.length=%d  area.x=%d area.y=%d\n",v4l2_fmt.fmt.pix.sizeimage,(int)area->size[0],(int)area->size[1]);
	*/
	
	for (y=0; y<area->size[1]; y++) {
		for (x=0; x<area->size[0]; x++) {
			data[(int)(area->size[0]*(area->size[1]-y-1)+x)]=dbuf[v4l2_fmt.fmt.pix.bytesperline*y+x*2]/256.0;
		}
	}

//	int sum=0;
//	for (x=0; x<v4l2_fmt.fmt.pix.sizeimage; x++) sum+=dbuf[x];
//	printf("%d\n",sum);
	
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
