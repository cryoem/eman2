#!/bin/env python
# e2boxer.py  07/27/2004  Steven Ludtke
# This program is used to box out particles from micrographs/CCD frames

from EMAN2 import *
from optparse import OptionParser
from math import *
import os
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage: %prog [options] <image>
	
Automatic and manual particle selection. This version is specifically aimed at square boxes
for single particle analysis."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--gui",action="store_true",help="Start the GUI for interactive boxing",default=False)
	parser.add_option("--box","-B",type="int",help="Box size in pixels",default=-1)
	parser.add_option("--ptclsize","-P",type="int",help="Approximate size (diameter) of the particle in pixels. Not required if reference particles are provided.",default=-1)
	parser.add_option("--refptcl","-R",type="string",help="A stack of reference images. Must have the same scale as the image being boxed.",default=None)
	parser.add_option("--refvol","-V",type="string",help="A 3D model to use as a reference for autoboxing",default=None)
	parser.add_option("--sym","-S",type="string",help="Symmetry of the 3D model",default=None)
	parser.add_option("--auto","-A",type="string",action="append",help="Autobox using specified method: circle, ref",default=[])
			
	(options, args) = parser.parse_args()
	if len(args)<1 : parser.error("Input image required")

	image=EMData()
	image.read_image(args[0])
	
	refptcl=None
	if options.refptcl :
		refptcl=EMData.read_images(options.refptcl)
		refbox=refptcl[0].get_xsize()
		print "%d reference particles read (%d x %d)"%(len(refptcl),refbox,refbox)
	
	if options.box<5 :
		if options.refptcl : options.box=refptcl[0].get_xsize()
		elif options.ptclsize : 
			options.box=good_boxsize(options.ptclsize*1.2)
		else : parser.error("Please specify a box size")
	else:
		if not options.box in good_box_sizes:
			print "Note: EMAN2 processing would be more efficient with a boxsize of %d"%good_boxsize(options.box)
	
	shrinkfactor=int(ceil(options.box/16))
	print "Shrink factor = ",shrinkfactor
	#shrinkfactor=int(ceil(image.get_ysize()/1024.0))
	#if options.box/shrinkfactor<12 : shrinkfactor/=2
	
	image.filter("normalize")
	shrink=image
	shrink.mean_shrink(shrinkfactor)		# shrunken original image
	
	# This confusing line insures the shrunken images has even dimensions
	# since odd FFTs haven't been fixed yet
	if shrink.get_xsize()&1 or shrink.get_ysize()&1 :
		shrink=shrink.get_clip(Region(0,0,(shrink.get_xsize()|1)^1,(shrink.get_ysize()|1)^1))

	# now we try to clean up long range density variations in the image
	flt=shrink.copy_head()
	flt.to_one()
	flt.filter("mask.sharp",{"outer_radius":options.box*2/shrinkfactor})
	flt/=(float(flt.get_attr("mean"))*flt.get_xsize()*flt.get_ysize())
	flt.filter("xform.phaseorigin")
	a=shrink.convolute(flt)
	a*=a.get_xsize()*a.get_ysize()
	shrink-=a
	a=None
	
	shrink2=shrink.copy(0)
	shrink2.filter("math.squared")
#	image=EMData()
#	image.read_image(args[0])
#	shrink.write_image("e.mrc")
#	shrink2.write_image("f.mrc")
		
	print "Autobox mode ",options.auto
	
	if "ref" in options.auto:
		if not refptcl: error_exit("Reference particles required")
		
		# refptcls will contain shrunken normalized reference particles
		refptcls=[]
		for n,i in enumerate(refptcl):
			ic=i.copy(0)
			refptcls.append(ic)
			ic.mean_shrink(shrinkfactor)
			# first a circular mask
			ic.filter("mask.sharp",{"outer_radius":ic.get_xsize()/2-1})
			
			# make the unmasked portion mean -> 0
			ic.add(-float(ic.get_attr("mean_nonzero")),1)
			ic.filter("normalize.unitlen")
#			ic.write_image("scaled_refs.hdf",-1)

		# prepare a mask to use for local sigma calculaton
		circle=shrink.copy_head()
		circle.to_one()
		circle.filter("mask.sharp",{"outer_radius":options.box/(shrinkfactor*2)-1})
		circle/=(float(circle.get_attr("mean"))*circle.get_xsize()*circle.get_ysize())
		
		ccfmean=shrink.calc_ccf(circle,True,None)
		ccfsig=shrink2.calc_ccf(circle,True,None)
		ccfmean.filter("math.squared")
		ccfsig-=ccfmean		# ccfsig is now pointwise standard deviation of local mean
		ccfsig.filter("math.sqrt")
#		shrink.write_image("z0.mrc")
#		ccfsig.write_image("z1.mrc")
		
		xs=shrink.get_xsize()
		ys=shrink.get_ysize()
		pks=[]
		for n,i in enumerate(refptcls):
			print n
			j=i.get_clip(Region(-(xs-i.get_xsize())/2,-(ys-i.get_ysize())/2,xs,ys))
#			j.write_image("0.%0d.mrc"%n)
			ccfone=shrink.calc_ccf(j,True,None)
#			ccfone.write_image("a.%0d.mrc"%n)
			ccfone/=ccfsig
#			ccfone.write_image("b.%0d.mrc"%n)
			sig=float(ccfone.get_attr("sigma"))
			ccfone.filter("mask.onlypeaks",{"npeaks":0})
#			ccfone.write_image("c.%0d.mrc"%n)
			pk=ccfone.calc_highest_locations(sig*4.0)
			for m,p in enumerate(pk):
				pk[m]=(-p.value,n,p.x,p.y)
			pks+=pk
			
		pks.sort()		# an ordered list of the best particle locations
		
		# this will produce a new list excluding any lower valued boxes within
		# 1/2 a box size of a higher one. It also rescales the boxes.
		# (ok, you could do this with clever syntax, but this is more readable)
		goodpks=[]
		bf=options.box/(shrinkfactor*2)
		for n,i in enumerate(pks):
			for nn,ii in enumerate(pks[:n]):
				if i[2]<bf or i[3]<bf or i[2]>xs-bf-1 or i[3]>ys-bf-1 : break
				if hypot(i[2]-ii[2],i[3]-ii[3])<bf : break
			else: goodpks.append((i[0],i[1],i[2]*shrinkfactor-options.box/2,i[3]*shrinkfactor-options.box/2))
		
		# This will optimize the center location of each particle and improve
		# the similarity calculation
		for n,i in enumerate(goodpks):
			b=EMData()
			b.read_image(args[0],0,0,Region(i[2],i[3],options.box,options.box))
			ba=refptcl[i[1].align("RotateTranslateFlip",{"to":b})]
			
			
		# Write EMAN1 style box database
		out=open(args[0][:-3]+"box","w")
		for i in goodpks[:500]:
			out.write("%d\t%d\t%d\t%d\t-3\n"%(i[2],i[3],options.box,options.box))
		
		out.close()
		
	if "circle" in options.auto:
		shrinksq=shrink.copy(0)
		shrinksq*=shrinksq			# shrunken original image squared
		
		# outer and inner ring mask
		outer=EMData()
		sbox=int(options.box/shrinkfactor)
		outer.set_size(shrink.get_xsize(),shrink.get_ysize(),1)
		outer.to_one()
		inner=outer.copy(0)
		
		outer.filter("mask.sharp",{"inner_radius":sbox*2/5,"outer_radius":sbox/2})
		inner.filter("mask.sharp",{"outer_radius":sbox*2/5})
		
		outer.write_image("b_outer.mrc")
		inner.write_image("b_inner.mrc")

		ccf1=shrinksq.calc_ccf(inner,True,None)
		ccf2=shrinksq.calc_ccf(outer,True,None)
		
		ccf1.write_image("b_ccf1.mrc")
		ccf2.write_image("b_ccf2.mrc")
		
		ccf1/=ccf2
		
		ccf1.write_image("b_result.mrc")
	
if __name__ == "__main__":
	main()
