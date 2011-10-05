#!/usr/bin/env python
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

from EMAN2 import *
from math import *
import os
import pickle
import sys


def main():
	global tdim,pdim
	global cmp_probe,cmp_target
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <tile file>
	
	Operates on files containing sets of tiled JPEG images representing larger images. Used for 
	interactive web browsing. Not generally useful to end-users."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--build", type=str, help="Build a new tile file from the specified image")
	parser.add_argument("--buildpspec",type=float,help="Builds 1D and 2D power spectra for the images when building, Value is A/pix for image.")
	parser.add_argument("--tilesize", type=int,default=256, help="Build a new tile file from this image")
	parser.add_argument("--dump",action="store_true",default=False,help="Dump the tile dictionary from the file")
	parser.add_argument("--display",type=str,default="",help="Displays a specific tile (level,x,y))")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
	if len(args)<1 : parser.error("Tile file required")
	try: chains=options.chains
	except: chains=None
		
	if options.build:
		# read the target and probe
		orig=EMData()
		orig.read_image(options.build)
		orig.process_inplace("normalize")
		opt={}
		try: opt["pspec"]=options.buildpspec
		except: pass
	
		build_tiles(orig,args[0],256,opt)

	if options.dump:
		td=tile_list(args[0])
		print td
		k=td.keys()
		k.sort()
		for i in k: print i,td[i]

	if options.display!="" :
		l,x,y=options.display.split(',')
		try: img=get_tile(args[0],int(l),int(x),int(y))
		except:
			print "Tile not present in file"
			sys.exit(1)
		o=file("/tmp/tile.jpg","w")
		o.write(img)
		o.close()
		os.system("display /tmp/tile.jpg")
		

def tile_list(tilefile):
	"""tile_list(tilefile)
	Extract dictionary of tiles from a tilefile"""
	
	tf=file(tilefile,"r")
	
	td=pickle.load(tf)

	tf.close()
	return td

def get_tile(tilefile,level,x,y):
	"""get_tile(tilefile,level,x,y)
	retrieve a tile from the file"""
	
	tf=file(tilefile,"r")
	
	td=pickle.load(tf)
	a=td[(level,x,y)]
	
	tf.seek(a[0],1)
	ret=tf.read(a[1])
	
	tf.close()
	return ret
	

def build_tiles(img,tilefile,tilesize,options=[]):
	"""build_tiles(img,tilefile,tilesize)
	This will construct a new tile file from a source image.
	options may include : pspec """
	levels=ceil(log(max(img.get_xsize(),img.get_ysize())/tilesize)/log(2.0))
	
	tf=file(tilefile,"w")
	
	tile_dict={}
	pos=0
	img2=img.copy()
	xs,ys=img2.get_xsize(),img2.get_ysize()
	for l in range(int(levels)):
		rmin=img2.get_attr("mean")-img2.get_attr("sigma")*3.0
		rmax=img2.get_attr("mean")+img2.get_attr("sigma")*3.0
		for x in range(0,img2.get_xsize(),tilesize):
			for y in range(0,img2.get_ysize(),tilesize):
				i=img2.get_clip(Region(x,y,tilesize,tilesize))
				i.set_attr("render_min",rmin)
				i.set_attr("render_max",rmax)
				i.set_attr("jpeg_quality",70)
				fsp="tmpimg.%d.%03d.%03d.jpg"%(l,x/tilesize,y/tilesize)
				i.write_image(fsp)
				sz=os.stat(fsp).st_size
				tile_dict[(l,x/tilesize,y/tilesize)]=(pos,sz)
				pos+=sz
		img2.process_inplace("math.meanshrink",{"n":2})
	
	# This will produce 2 power spectrum images in the tile file
	# with scale factors -1 and -2
	if "pspec" in options :
		nx,ny=img.get_xsize()/512,img.get_ysize()/512
		a=EMData()
		a.set_size(512,512)
		if (ny>2 and nx>2) :
			for y in range(1,ny-1):
				for x in range(1,nx-1):
					c=img.get_clip(Region(x*512,y*512,512,512))
					c.process_inplace("normalize")
					c.process_inplace("math.realtofft")
					c.process_inplace("math.squared")
					a+=c
			a.set_value_at(256,256,0,.01)
			a-=a.get_attr("minimum")-a.get_attr("sigma")*.01
			a.process_inplace("math.log")
			a-=a.get_attr("minimum")
			a.set_attr("render_min",a.get_attr("minimum")-a.get_attr("sigma")*.1)
			a.set_attr("render_max",a.get_attr("mean")+a.get_attr("sigma")*4.0)
			a.set_attr("jepg_quality",80)
			a.write_image("/tmp/tmpimg.mrc")
			fsp="tmpimg.jpg"
			a.write_image(fsp)
			sz=os.stat(fsp).st_size
			tile_dict[(-1,0,0)]=(pos,sz)
			pos+=sz
	
#		try:
			import matplotlib
			matplotlib.use('Agg')
			import pylab
			manager = pylab.get_current_fig_manager()
			apix=options["pspec"]
			dx=1.0/(2.0*apix*256.0)
			x=pylab.arange(dx,dx*255.9,dx)
			y=a.calc_radial_dist(255,1,1,0)	# radial power spectrum (log)
			pylab.figure(figsize=(8,6),dpi=96)
			pylab.axes([.08,.08,.9,.9], axisbg='w')
			pylab.plot(x,y)
			pylab.axis([0,dx*256,min(y),max(y)])
			pylab.xlabel("Spatial Freq. (1/A)")
			pylab.ylabel("Log Intensity (10^x)")
#			print y
			
			fsp="tmpimg2.png"
			pylab.savefig(fsp,dpi=96)
			sz=os.stat(fsp).st_size
			tile_dict[(-2,0,0)]=(pos,sz)
			pos+=sz

#		except:
#			print "Unable to generate plot (need matplotlib)"
			
	
	pickle.dump(tile_dict,tf)
	
	for l in range(int(levels)):
		for x in range(0,xs,tilesize):
			for y in range(0,ys,tilesize):
				fsp="tmpimg.%d.%03d.%03d.jpg"%(l,x/tilesize,y/tilesize)
				a=file(fsp,"r")
				b=a.read()
				a.close()
				tf.write(b)
				os.remove(fsp)
		xs/=2
		ys/=2
	
	if "pspec" in options :
		for fsp in ["tmpimg.jpg","tmpimg2.png"] :
			a=file(fsp,"r")
			b=a.read()
			a.close()
			tf.write(b)
#			os.remove(fsp)
	
	tf.close()

if __name__ == "__main__":
    main()
