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

#ifdef ENABLE_V4L2
#include "v4l2io.h"
#include "geometry.h"
#include "util.h"
#include "portable_fileio.h"

#include <climits>
#include <cctype>
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
//	if (v4l_file) {
		//close(v4l_file);
		//v4l_file = 0;
	//}
	if (filename) free(filename);
}

void V4L2IO::init() {
	ENTERFUNC;
	static int ginit=-1;
	
	if (ginit==-1) ginit=globalinit(filename,0,90,48,8,200,210); // planets on webcam
//	if (ginit==-1) ginit=globalinit(filename,0,80,48,0,240,220); // good for deep sky on webcam
	//if (ginit==-1) ginit=globalinit(filename,0,12,-1,60);
	v4l_file=ginit;
	
	initialized=true;

	EXITFUNC;
}

int V4L2IO::globalinit(const char *fsp,int input,int brt,int cont,int gamma,int expos,int gain)
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
	
	if (-1 == xioctl (vfile, VIDIOC_S_INPUT, &input)) errno_exit ("VIDIOC_S_INPUT");
	
	int std=V4L2_STD_NTSC_M;
	if (-1 == xioctl (vfile, VIDIOC_S_STD, &std)) printf("Can't set NTSC standard\n");
	
	
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
	v4l2_fmt.fmt.pix.width       = 640; 
	v4l2_fmt.fmt.pix.height      = 480;
//	v4l2_fmt.fmt.pix.pixelformat = V4L2_PIX_FMT_GREY;
	v4l2_fmt.fmt.pix.pixelformat = V4L2_PIX_FMT_YUV420;
//	v4l2_fmt.fmt.pix.pixelformat = V4L2_PIX_FMT_YUYV;
//	v4l2_fmt.fmt.pix.field       = V4L2_FIELD_NONE;
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
	
	printf("fmt.fmt.pix.width = %d\n",v4l2_fmt.fmt.pix.width);
	printf("fmt.fmt.pix.height = %d\n",v4l2_fmt.fmt.pix.height);
	printf("fmt.fmt.pix.pixelformat = %4s\n",&v4l2_fmt.fmt.pix.pixelformat);
	printf("fmt.fmt.pix.bytesperline = %d\n",v4l2_fmt.fmt.pix.bytesperline);
	printf("fmt.fmt.pix.sizeimage = %d\n",v4l2_fmt.fmt.pix.sizeimage);
//	printf("fmt.fmt.pix.field = %d\n",v4l2_fmt.fmt.pix.field); 
	
	struct v4l2_queryctrl qc;
	struct v4l2_control con;

	qc.id=V4L2_CID_BRIGHTNESS;
	ioctl(vfile,VIDIOC_QUERYCTRL,&qc);
	printf("brightness = %d - %d by %d %d\n",qc.minimum,qc.maximum,qc.step,qc.default_value);

	qc.id=V4L2_CID_CONTRAST;
	ioctl(vfile,VIDIOC_QUERYCTRL,&qc);
	printf("contrast = %d - %d by %d %d\n",qc.minimum,qc.maximum,qc.step,qc.default_value);

	qc.id=V4L2_CID_GAMMA;
	ioctl(vfile,VIDIOC_QUERYCTRL,&qc);
	printf("gamma = %d - %d by %d %d\n",qc.minimum,qc.maximum,qc.step,qc.default_value);

	qc.id=V4L2_CID_EXPOSURE;
	ioctl(vfile,VIDIOC_QUERYCTRL,&qc);
	printf("exposure = %d - %d by %d %d\n",qc.minimum,qc.maximum,qc.step,qc.default_value);

	qc.id=V4L2_CID_GAIN;
	ioctl(vfile,VIDIOC_QUERYCTRL,&qc);
	printf("gain = %d - %d by %d %d\n",qc.minimum,qc.maximum,qc.step,qc.default_value);

	con.id=V4L2_CID_AUTOGAIN;
	con.value=0;
	ioctl(vfile,VIDIOC_S_CTRL,&con);
	
	if (brt!=-1) {
	con.id=V4L2_CID_BRIGHTNESS;
	con.value=brt;
	ioctl(vfile,VIDIOC_S_CTRL,&con);
	}

	if (cont!=-1) {
	con.id=V4L2_CID_CONTRAST;
	con.value=cont;
	ioctl(vfile,VIDIOC_S_CTRL,&con);
	}

	if (gamma!=-1) {
	con.id=V4L2_CID_GAMMA;
	con.value=gamma;
	ioctl(vfile,VIDIOC_S_CTRL,&con);
	}

	if (expos!=-1) {
	con.id=V4L2_CID_EXPOSURE;
	con.value=expos;
	ioctl(vfile,VIDIOC_S_CTRL,&con);
	}

	if (gain!=-1) {
	con.id=V4L2_CID_GAIN;
	con.value=gain;
	ioctl(vfile,VIDIOC_S_CTRL,&con);
	}

//	close(vfile);
	EXITFUNC;
	return vfile;
}

bool V4L2IO::is_valid(const void *first_block)
{
	ENTERFUNC;
	return true;
}

int V4L2IO::read_header(Dict & dict, int image_index, const Region * area, bool)
{
	ENTERFUNC;
	
	if(image_index == -1) {
		image_index = 0;
	}

	if(image_index != 0) {
		throw ImageReadException(filename, "no stack allowed for MRC image. For take 2D slice out of 3D image, read the 3D image first, then use get_clip().");
	}
	
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
	LOGWARN("V4L write is not supported.");
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
// This is for YUYV
//			data[(int)(area->size[0]*(area->size[1]-y-1)+x)]=dbuf[v4l2_fmt.fmt.pix.bytesperline*y+x*2]/256.0;
// This is for YUV420
//			data[(int)(area->size[0]*(area->size[1]-y-1)+x)]=dbuf[v4l2_fmt.fmt.pix.bytesperline*y+x]/256.0;
			data[(int)(area->size[0]*(area->size[1]-y-1)+x)]=dbuf[(int)area->size[0]*y+x]/256.0;
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
