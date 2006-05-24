#!/usr/bin/env python

from EMAN2 import *
from optparse import OptionParser
from math import *
import os
import pickle
import sys

def main():
	global tdim,pdim
	global cmp_probe,cmp_target
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <tile file>
	
Operates on files containing sets of tiled JPEG images representing larger images. Used for 
interactive web browsing."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--build", type="string", help="Build a new tile file from this image")
	parser.add_option("--tilesize", type="int",default=256, help="Build a new tile file from this image")
	parser.add_option("--dump",action="store_true",default=False,help="Dump the tile dictionary from the file")
	parser.add_option("--display",type='string',default="",help="Displays a specific tile (level,x,y))")
	
	(options, args) = parser.parse_args()
	if len(args)<1 : parser.error("Tile file required")
	try: chains=options.chains
	except: chains=None
		
	if options.build:
		# read the target and probe
		orig=EMData()
		orig.read_image(options.build)
	
		build_tiles(orig,args[0],256)

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
	

def build_tiles(img,tilefile,tilesize):
	"""build_tiles(img,tilefile,tilesize)
	This will construct a new tile file from a source image."""
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
				i.set_attr("jpeg_quality",60)
				fsp="tmpimg.%d.%03d.%03d.jpg"%(l,x/tilesize,y/tilesize)
				i.write_image(fsp)
				sz=os.stat(fsp).st_size
				tile_dict[(l,x/tilesize,y/tilesize)]=(pos,sz)
				pos+=sz
		img2.mean_shrink(2)
	
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
	
	tf.close()

if __name__ == "__main__":
    main()
