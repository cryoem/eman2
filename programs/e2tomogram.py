#!/bin/env python
# e2tomogram.py  11/21/2004  Steven Ludtke
# This program is used to analyze tilt series


from EMAN2 import *
from optparse import OptionParser
from math import *
import os
import sys

def matrixalign(im1,im2,box,padbox,debug=0) :
	"""This will calculate a set of alignment vectors between two images
box is the size of the area used to find the alignment vectors, padbox
is the size of the expanded box to search for the smaller box in.
returns a list of tuples :  (I,x,y,dx,dy)  where I is the alignment peak
intensity."""
	sx=im1.get_xsize()
	sy=im1.get_ysize()
	bigpad=(padbox-box)/2
	nx=int((sx-2*bigpad)/box)*2
	ny=int((sy-2*bigpad)/box)*2
	
	ret=[]
	for y in range(ny):
		for x in range(nx):
			clip2=im2.get_clip(Region(x*box/2,y*box/2,padbox,padbox))
			# if there are zeroes in the clipped area we don't want to do correlations here
#			if float(clip2.get_attr("mean"))!=float(clip2.get_attr("mean_nonzero")) # continue
			
			clip1=im1.get_clip(Region(x*box/2,y*box/2,padbox,padbox))
			
			# mask out the center of im1 and find it within im2
			
			# note that this normalization, making the masked out region mean value
			# exactly 0, is critical to obtaining correct alignments, since it insures
			# that
			mask=clip1.copy_head()
			mask.to_one()
			mask.filter("MaskSharp",{"outer_radius":box/2})
			clip1*=mask
			clip1-=float(clip1.get_attr("mean_nonzero"))
			clip1*=mask
			
			clip2.filter("NormalizeStd")
			clip2s=clip2.copy(0)
			clip2s.filter("ValueSquared")
			
			ccf=clip1.calc_ccf(clip2,1)
			ccfs=mask.calc_ccf(clip2s,1)	# this is the sum of the masked values^2 for each pixel center
			ccfs.filter("ValueSqrt")
			ccf/=ccfs
	
			ccf.filter("NormalizeStd")		# peaks relative to 1 std-dev
			ccf.filter("MaskSharp",{"outer_radius":bigpad})		# max translation
			
			if (debug):
				clip1.write_image("dbug.hed",-1)
				clip2.write_image("dbug.hed",-1)
				ccf.write_image("dbug.hed",-1)
			maxloc=ccf.calc_max_location()
			
			ret.append((float(ccf.get_attr("maximum")),x*box+bigpad/2-sx/2,y*box+bigpad/2-sy/2,maxloc[0]-padbox/2,maxloc[1]-padbox/2))

	return ret
	
def mode(vals):
	"""calculates the most common value in a list, if no value is more
	common, returns the median"""
	d={}
	for i in vals:
		try: d[i]=d[i]+1
		except: d[i]=1

	cnt=[(i[1],i[0]) for i in d.items()]
	
	cnt.sort()
	
	try:
		if cnt[-1][0]==cnt[-2][0] :
			vals.sort()
			return vals[len(vals)/2]
	except:
		pass
		
	return cnt[-1][1]
	
def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage: %prog [options] input_stack.hed output.hed
	
Processes a tomographic tilt series"""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--tilt", "-T", type="float", help="Angular spacing between tilts (fixed)",default=0.0)
	parser.add_option("--maxshift","-M", type="int", help="Maximum translational error between images (pixels), default=64",default=64.0)
	parser.add_option("--box","-B", type="int", help="Box size for alignment probe (pixels), default=96",default=96.0)
	parser.add_option("--highpass",type="float",help="Highpass Gaussian filter radius (pixels), default none", default=-1.0)
	parser.add_option("--lowpass",type="float",help="Lowpass Gaussian filter radius (pixels), default none",default=-1.0)
#	parser.add_option("--apix", "-A", type="float", help="A/voxel", default=1.0)
#	parser.add_option("--res", "-R", type="float", help="Resolution in A, equivalent to Gaussian lowpass with 1/e width at 1/res",default=2.8)
#	parser.add_option("--box", "-B", type="string", help="Box size in pixels, <xyz> or <x>,<y>,<z>")
#	parser.add_option("--het", action="store_true", help="Include HET atoms in the map", default=False)
#	parser.add_option("--chains",type="string",help="String list of chain identifiers to include, eg 'ABEFG'")
#	parser.add_option("--quiet",action="store_true",default=False,help="Verbose is the default")
	
	(options, args) = parser.parse_args()
	if len(args)<2 : parser.error("Input and output files required")
	if options.tilt<=0 : parser.error("--tilt must be specified")
	
	nimg=EMUtil.get_image_count(args[0])
	if (nimg<3) : parser.error("Input file must contain at least 3 images")
	
	# copy the file with possible format conversion
#	Log.logger().set_level(Log.LogLevel.VARIABLE_LOG)
	for i in range(nimg):
		a=EMData()
		a.read_image(args[0],i)
		a.write_image(args[1],i)
		
	
	cmplist=[(x,x+1) for x in range(nimg/2,nimg-1)]+[(x,x-1) for x in range(nimg/2,0,-1)]
	for i in cmplist:
		im1=EMData()
		im1.read_image(args[1],i[0])
		im1.filter("NormalizeEdgeMean")
		if options.highpass>0 :im1.filter("HighpassGauss",{"highpass":options.highpass})
		if (options.lowpass>0) : im1.filter("LowpassGauss",{"lowpass":options.lowpass})
		im2=EMData()
		im2.read_image(args[1],i[1])
		im2.filter("NormalizeEdgeMean")
		if options.highpass>0 : im2.filter("HighpassGauss",{"highpass":options.highpass})
		if (options.lowpass>0) : im2.filter("LowpassGauss",{"lowpass":options.lowpass})
		
		vec=matrixalign(im1,im2,options.box,options.box+options.maxshift*2,i[0]==63)
		
		vec.sort()			# sort in order of peak height
		vec2=vec[-len(vec)/4:]		# take the 25% strongest correlation peaks
		
		vec3=[(hypot(x[1],x[2]),x[0],x[1],x[2],x[3],x[4]) for x in vec2]
		vec3.sort()					# sort in order of distance from center
		vec4=vec3[:len(vec3)/2]		# take the 1/2 closest to the center
#		vec4=vec3
#		for x in vec4: print x
		
		dxs=[int(x[4]) for x in vec4]
		dys=[int(x[5]) for x in vec4]
		
		best=(mode(dxs),mode(dys))
	
		print i,best
		im2.rotate_translate(0,0,0,best[0],best[1],0)
		im2.write_image(args[1],i[1])
			
if __name__ == "__main__":
    main()
