#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
#
#

from EMAN2_cppwrap import *
from bisect import bisect_left
from pyemtbx.imagetypes import *
from pyemtbx.box import *
#from Sparx import *
from sys import exit
import os
import time
import shelve
import re

EMANVERSION="EMAN2 v1.90"

Vec3f.__str__=lambda x:"Vec3f"+str(x.as_list())

Transform3D.__str__=lambda x:"Transform3D(\t%7.4g\t%7.4g\t%7.4g\n\t\t%7.4g\t%7.4g\t%7.4g\n\t\t%7.4g\t%7.4g\t%7.4g)\nPretrans:%s\nPosttrans:%s"%(x.at(0,0),x.at(0,1),x.at(0,2),x.at(1,0),x.at(1,1),x.at(1,2),x.at(2,0),x.at(2,1),x.at(2,2),str(x.get_pretrans()),str(x.get_posttrans()))

GUIMode=0
GUIbeingdragged=None

def timer(fn,n=1):
	a=time.time()
	for i in range(n): fn()
	print time.time()-a

def E2init(argv) :
        from OpenGL import GLUT
	"""E2init(argv)
This function is called to log information about the current job to the local logfile"""
	try:
		db=shelve.open(".eman2log")
	except:
		return -1
		
	try:
		n=db["count"]
		db["count"]=n+1
	except:
		n=1
		db["count"]=n
	db[str(n)]={"pid":os.getpid(),"start":time.time(),"args":argv}
	db.close()
	GLUT.glutInit( argv )

	return n

def E2end(n):
	"""E2end(n)
This function is called to log the end of the current job. n is returned by E2init"""
	db=shelve.open(".eman2log")
	d=db[str(n)]
	d["end"]=time.time()
	db[str(n)]=d
	db.close()
	
	return n

parseparmobj1=re.compile("([^\(]*)\(([^\)]*)\)")	# This parses test(n=v,n2=v2) into ("test","n=v,n2=v2")
parseparmobj2=re.compile("([^=,]*)=([^,]*)")		# This parses "n=v,n2=v2" into [("n","v"),("n2","v2")]
def parsemodopt(optstr):
	"""This is used so the user can provide the name of a comparitor, processor, etc. with options
	in a convenient form. It will parse "dot(normalize=1,negative=0)" and return
	("dot",{"normalize":1,"negative":0})"""
	
	if not optstr or len(optstr)==0 : return (None,{})
	
	p1=re.findall(parseparmobj1,optstr)
	if len(p1)==0 : return (optstr,{})
	
	p2=re.findall(parseparmobj2,p1[0][1])
	r2={}
	# this will convert values to integers, floats or strings. Note that an integer '1' can be forced to be 
	# considered a string with test(name="1")
	for i,j in enumerate(p2):
		v=j[1]
		try: v=int(v)
		except:
			try: v=float(v)
			except:
				if v[0]=='"' and v[-1]=='"' : v=j[1][1:-1]
		r2[j[0]]=v
		
	return (p1[0][0],r2)

def display(img):
	
	if GUIMode:
		import emimage
		return emimage.EMImage(img)
	else:
		# In non interactive GUI mode, this will display an image or list of images with e2display
		try: os.unlink("/tmp/img.hed")
		except: pass
		try: os.unlink("/tmp/img.img")
		except: pass
		if isinstance(img,list) or isinstance(img,tuple) :
			for i in img: i.write_image("/tmp/img.hed",-1)
		else:
			img.write_image("/tmp/img.hed")
	#	os.system("v2 /tmp/img.hed")
		os.system("e2display.py /tmp/img.hed")

def plot(data,show=1,size=(800,600),path="plot.png"):
	"""plots an image or an array using the matplotlib library"""
	import matplotlib
	matplotlib.use('Agg')
	import pylab
	pylab.figure(figsize=(size[0]/72.0,size[1]/72.0),dpi=72)
	if isinstance(data,EMData) :
		a=[]
		for i in range(data.get_xsize()): 
			a.append(data.get_value_at(i,0,0))
		pylab.plot(a)
	elif isinstance(data,list) or isinstance(data,tuple):
		if isinstance(data[0],list) or isinstance(data[0],tuple) :
			pylab.plot(data[0],data[1])
		else:
			try:
				a=float(data[0])
				pylab.plot(data)
			except:
				print "List, but data isn't floats"
				return
	else :
		print "I don't know how to plot that type"
		return
	
	pylab.savefig(path)
	if show:
		try: os.system("display "+path)
		except: pass

	if show: os.system("display "+path)

def qplot(img):
	"""This will plot a 1D image using qplot"""
	out=file("/tmp/plt.txt","w")
	for i in range(img.get_xsize()):
		out.write("%d\t%f\n"%(i,img.get_value_at(i,0)))
	out.close()
	os.system("qplot /tmp/plt.txt")

def error_exit(s) :
	"""A quick hack until I can figure out the logging stuff. This function
	should still remain as a shortcut"""
	print s
	exit(1)
	
def test_image(type=0,size=(128,128)):
	"""Returns a simple standard test image
	type=0  scurve
	type=1	gaussian noise, 0 mean, sigma 1
	type=2  square
	type=3  hollow square
	type=4  circular sinewave
	size=(128,128) """
	ret=EMData()
	ret.set_size(*size)
	if type==0 :
		ret.process_inplace("testimage.scurve")
	elif type==1 :
		ret.process_inplace("testimage.noise.gauss")
	elif type==2:
		ret.process_inplace("testimage.squarecube",{"edge_length":size[0]/2})
	elif type==3:
		ret.process_inplace("testimage.squarecube",{"fill":1,"edge_length":size[0]/2})
	elif type==4:
		ret.process_inplace("testimage.sinewave.circular")
	
	return ret

def isosurface(marchingcubes, threshhold, smooth=False):
	"""Return the Isosurface points, triangles, normals(smooth=True), normalsSm(smooth=False)"""
	marchingcubes.set_surface_value(threshhold)
	d = marchingcubes.get_isosurface(smooth)
	return d['points'], d['faces'], d['normals']

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
  Processors - Processor factory. Each processor implements an image-processing algorithm. \n\
  Projectors - Projector factory. Each Projector implements an algorithm for 3D image projection. \n\
  Reconstructors - Reconstructor factory. Each Reconstructor implements an algorithm for image reconstruction. \n\
\n\
  EMUtil - Defines utility functions for EMData-related operations. \n\
  TestUtil - Unit test utility functions. \n\
  Util - Generic utility functions. \n\
\n\
  EMNumPy - Defines functions for conversions between EMData and numpy array. \n\
  Log - Log output at different verbose level. \n\
  PointArray - Point array. \n\
\n\
  dump_aligners() - Print out all Aligners and their parameters. \n\
  dump_averagers() - Print out all Averagers and their sparameters. \n\
  dump_cmps() - Print out all Cmps and their parameters. \n\
  dump_processors() - Print out all Processor`s and their parameters. \n\
  dump_projectors() - Print out all Projectors and their parameters. \n\
  dump_reconstructors() - Print out all Reconstructors and their parameters. \n\
  dump_analyzers() - Print out all Analyzers anf their parameters. \n\
"
  
  
