from libpyAligner2 import *
from libpyAverager2 import *
from libpyCmp2 import *
from libpyFilter2 import *
from libpyReconstructor2 import * 
from libpyProjector2 import *
from libpyEMObject2 import * 
from libpyEMData2 import *
from libpyGeometry2 import *
from libpyTransform2 import *
from libpyUtils2 import * 
from libpyPointArray2 import *
from libpyTypeConverter2 import *
from bisect import bisect_left
from pyemtbx.imagetypes import *
from pyemtbx.box import *


EMANVERSION="EMAN2 v1.90"

Vec3f.__str__=lambda x:"Vec3f"+str(x.as_list())

Transform.__str__=lambda x:"Transform(\t%7.4g\t%7.4g\t%7.4g\n\t\t%7.4g\t%7.4g\t%7.4g\n\t\t%7.4g\t%7.4g\t%7.4g)\nPretrans:%s\nPosttrans:%s"%(x.at(0,0),x.at(0,1),x.at(0,2),x.at(1,0),x.at(1,1),x.at(1,2),x.at(2,0),x.at(2,1),x.at(2,2),str(x.get_pretrans()),str(x.get_posttrans()))

def display(img):
	"""This will use 'v2', and EMAN1 program to view an image
	or a list/tuple of images. This is basically a hack."""
	os.unlink("/tmp/img.spi")
	if isinstance(img,list) or isinstance(img,tuple) :
		for i in img: i.write_image("/tmp/img.spi",-1)
	else:
		img.write_image("/tmp/img.hdf")
	from os import system
	system("v2 /tmp/img.spi")


__doc__ = \
"EMAN classes and routines for image/volume processing in \n\
single particle reconstructions.\n\
\n\
The following classes are defined: \n\
  EMData - the primary class to process electronic microscopy images. \n\
  Transform - defines a transformation. \n\
"


