#!/bin/env python
# e2tomogram.py  11/21/2004  Steven Ludtke
# This program is used to analyze tilt series


from EMAN2 import *
from optparse import OptionParser
from math import *
import os
import sys

def matrixalign(im1,im2,box,padbox,maxrange=64,debug=0) :
	"""This will calculate a set of alignment vectors between two images
box is the size of the area used to find the alignment vectors, padbox
is the size of the expanded box to search for the smaller box in.
returns a dict of tuples :  (I,x,y,dx,dy)  where I is the alignment peak
intensity. The key is the box number centered on the origin, ie - (0,0)
is the alignment of the center of the image (1,0) is one box to the right
of the center. maxrange allows calculating a limited distance from the center"""
	sx=im1.get_xsize()
	sy=im1.get_ysize()
	bigpad=(padbox-box)/2
	nx=int((sx-2*bigpad)/box)*2-1	# this insures that we have a box at the origin
	ny=int((sy-2*bigpad)/box)*2-1
	dx=(sx-padbox)/float(nx-1)
	dy=(sy-padbox)/float(ny-1)
	ret={}
	for y in range(ny):
		for x in range(nx):
			if (abs(x-nx/2-1)>maxrange or abs(y-ny/2-1)>maxrange) : continue
			
			clip2=im2.get_clip(Region(int(x*dx),int(y*dy),padbox,padbox))
			# if there are zeroes in the clipped area we don't want to do correlations here
#			if float(clip2.get_attr("mean"))!=float(clip2.get_attr("mean_nonzero")) # continue
			
			clip1=im1.get_clip(Region(int(x*dx),int(y*dy),padbox,padbox))
			
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
			
			ret[(x-nx/2-1,y-ny/2-1)]=((float(ccf.get_attr("maximum")),x*box+bigpad/2-sx/2,y*box+bigpad/2-sy/2,maxloc[0]-padbox/2,maxloc[1]-padbox/2))

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
	parser.add_option("--mode",type="string",help="centering mode 'modeshift' or 'censym'",default="censym")
#	parser.add_option("--het", action="store_true", help="Include HET atoms in the map", default=False)
	
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
		
		
		if options.mode=="modeshift" :
			dct=matrixalign(im1,im2,options.box,options.box+options.maxshift*2,debug=i[0]==63)
			vec=dct.values()
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
		elif options.mode=="censym" :
			dct=matrixalign(im1,im2,options.box,options.box+options.maxshift*2,2)
			pairs=[]
			for x in range(3):
				for y in range(-2,3):
					if y<=0 and x==0 : continue
					a=dct[(x,y)]
					b=dct[(-x,-y)]
					if hypot(a[3]-b[3],a[4]-b[4])>6.0 : continue
					pairs.append((x,y,(a[3]+b[3])/2.0,(a[4]+b[4])/2.0,hypot(a[3]-b[3],a[4]-b[4])))
#					print "%d,%d\t%5.2f %5.2f\t%5.2f"%(pairs[-1][0],pairs[-1][1],pairs[-1][2],pairs[-1][3],pairs[-1][4])
			
			if len(pairs)==0 : 
				print "Alignment failed on image %d (%d)"%(i[1],i[0])
				if (i[0]>0) : cmplist.append((i[0]-1,i[1]))
				else : cmplist.append((i[0]+1,i[1]))
				best=(0,0)
			else :
				# start by finding the average pair-matched shift
				sum=[0,0]
				norm=0
				for p in pairs:
					sum[0]+=p[2]*(1.0/(1.0+p[4]))
					sum[1]+=p[3]*(1.0/(1.0+p[4]))
					norm+=(1.0/(1.0+p[4]))
				best=(sum[0]/norm,sum[1]/norm)

				# now do it again, but exclude any outliers from the average
				sum=[0,0]
				norm=0
				for p in pairs:
					if hypot(p[2]-best[0],p[3]-best[1])>5.0 :continue
					sum[0]+=p[2]*(1.0/(1.0+p[4]))
					sum[1]+=p[3]*(1.0/(1.0+p[4]))
					norm+=(1.0/(1.0+p[4]))
					best=(sum[0]/norm,sum[1]/norm)
					
		print "%d.\t%5.2f\t%5.2f"%(i[1],best[0],best[1])
		im2.rotate_translate(0,0,0,best[0],best[1],0)
		im2.filter("NormalizeStd")
		im2.write_image(args[1],i[1])
	
	return
	# now we look for the common-line in the aligned images
	sum=im1.do_fft()
	sum.to_zero()
	for i in range(nimg):
		a=EMData()
		a.read_image(args[1],i)
		b=a.do_fft()
		sum+=b
	
	sum=sum.do_ift()
	sum.write_image("fft.mrc",0)
	
			
if __name__ == "__main__":
    main()
