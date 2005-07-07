from EMAN2_cppwrap import *
from bisect import bisect_left
from pyemtbx.imagetypes import *
from pyemtbx.box import *
#from Sparx import *
from sys import exit
import os

EMANVERSION="EMAN2 v1.90"

Vec3f.__str__=lambda x:"Vec3f"+str(x.as_list())

Transform3D.__str__=lambda x:"Transform3D(\t%7.4g\t%7.4g\t%7.4g\n\t\t%7.4g\t%7.4g\t%7.4g\n\t\t%7.4g\t%7.4g\t%7.4g)\nPretrans:%s\nPosttrans:%s"%(x.at(0,0),x.at(0,1),x.at(0,2),x.at(1,0),x.at(1,1),x.at(1,2),x.at(2,0),x.at(2,1),x.at(2,2),str(x.get_pretrans()),str(x.get_posttrans()))

def display(img):
	"""This will use 'v2', and EMAN1 program to view an image
	or a list/tuple of images. This is basically a hack."""
	try: os.unlink("/tmp/img.hdf")
	except: pass
	if isinstance(img,list) or isinstance(img,tuple) :
		for i in img: i.write_image("/tmp/img.hdf",-1)
	else:
		img.write_image("/tmp/img.hdf")
	os.system("v2 /tmp/img.hdf")

def error_exit(s) :
	"""A quick hack until I can figure out the logging stuff. This function
	should still remain as a shortcut"""
	print s
	exit(1)
	

__doc__ = \
"EMAN classes and routines for image/volume processing in \n\
single particle reconstructions.\n\
\n\
The following classes are defined: \n\
  EMData - the primary class to process electronic microscopy images. \n\
 \n\
  Quaternion - implements mathematical quaternion. \n\
  Region - defines a rectangular 2D/3D region. \n\
  Transform3D - defines a transformation including rotation, translation, and different Euler angles. \n\
  Vec3i - a 3-element integer vector. \n\
  Vec3f - a 3-element floating number vector. \n\
\n\
  EMObject - A wrapper class for int, float, double, string, EMData, XYData, list etc. \n\
  Pixel - defines a pixel's 3D coordinates and its value. \n\
  SimpleCtf - defines EMAN CTF parameter model. \n\
  XYData - implements operations on a series of (x y) data pair. \n\
\n\
  Aligners - Aligner factory. Each Aligner alignes 2D/3D images. \n\
  Averagers - Averager factory. Each Averager averages a set of images. \n\
  Cmps  - Cmp factory. Each Cmp defines a way to do image comparison. \n\
  Filters - Filter factory. Each processor implements an image-processing algorithm. \n\
  Projectors - Projector factory. Each Projector implements an algorithm for 3D image projection. \n\
  Reconstructors - Reconstructor factory. Each Reconstructor implements an algorithm for image reconstruction. \n\
\n\
  EMUtil - Defines utility functions for EMData-related operations. \n\
  TestUtil - Unit test utility functions. \n\
  Util - Generic utility functions. \n\
\n\
  EMNumPy - Defines functions for conversions between EMData and Numeric python array. \n\
  Log - Log output at different verbose level. \n\
  PointArray - Point array. \n\
\n\
  dump_aligners() - Print out all Aligners and their parameters. \n\
  dump_averagers() - Print out all Averagers and their parameters. \n\
  dump_cmps() - Print out all Cmps and their parameters. \n\
  dump_processors() - Print out all Filters and their parameters. \n\
  dump_projectors() - Print out all Projectors and their parameters. \n\
  dump_reconstructors() - Print out all Reconstructors and their parameters. \n\
"
  
  
