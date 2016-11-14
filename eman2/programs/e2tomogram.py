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

# e2tomogram.py  11/21/2004  Steven Ludtke
# This program is used to analyze tilt series


from EMAN2 import *
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
	bigpad=(padbox-box)/2			# ostensibly this is the max translation we should
									# allow, but we actually allow 2x this
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
			mask.process_inplace("mask.sharp",{"outer_radius":box/2})
			clip1*=mask
			clip1-=float(clip1.get_attr("mean_nonzero"))
			clip1*=mask
			
			clip2.process_inplace("normalize")
			clip2s=clip2.copy()
			clip2s.process_inplace("math.squared")
			
			ccf=clip1.calc_ccf(clip2,fp_flag.CIRCULANT)
			ccfs=mask.calc_ccf(clip2s,fp_flag.CIRCULANT)	# this is the sum of the masked values^2 for each pixel center
			ccfs.process_inplace("math.sqrt")
			ccf/=ccfs
	
			ccf.process_inplace("normalize")		# peaks relative to 1 std-dev
			if bigpad*2>padbox/2 : ccf.process_inplace("mask.sharp",{"outer_radius":padbox/2-1})
			else : ccf.process_inplace("mask.sharp",{"outer_radius":bigpad*2})		# max translation
			ccf.set_value_at(int(padbox/2),int(padbox/2),0,0)		# remove 0 shift artifacts
			
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
	usage = """prog [options] input_stack.hed output.hed
	
	Fiducial-less alignment of tomograms. This program has many limitations, and is still being developed.
	Not yet recommended for routine use.
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--tilt", "-T", type=float, help="Angular spacing between tilts (fixed)",default=0.0)
	parser.add_argument("--maxshift","-M", type=int, help="Maximum translational error between images (pixels), default=64",default=64.0)
	parser.add_argument("--box","-B", type=int, help="Box size for alignment probe (pixels), default=96",default=96.0)
	parser.add_argument("--highpass",type=float,help="Highpass Gaussian processor radius (pixels), default none", default=-1.0)
	parser.add_argument("--lowpass",type=float,help="Lowpass Gaussian processor radius (pixels), default none",default=-1.0)
	parser.add_argument("--mode",type=str,help="centering mode 'modeshift', 'censym' or 'region,<x>,<y>,<clipsize>,<alisize>",default="censym")
	parser.add_argument("--localavg",type=int,help="Average several images for the alignment",default=1)
	parser.add_argument("--tiltaxis",type=float,help="Skip automatic tilt axis location, use fixed angle from x",default=400.0)
	parser.add_argument("--twopass",action="store_true",default=False,help="Skip automatic tilt axis location, use fixed angle from x")
	parser.add_argument("--nozero",action="store_true",default=False,help="Do not allow 0-translations between images")
	#parser.add_argument("--het", action="store_true", help="Include HET atoms in the map", default=False)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
	if len(args)<2 : parser.error("Input and output files required")
	if options.tilt<=0 : parser.error("--tilt must be specified")
	
	nimg=EMUtil.get_image_count(args[0])
	if (nimg<3) : parser.error("Input file must contain at least 3 images")
	
#	Log.logger().set_level(Log.LogLevel.VARIABLE_LOG)
	
	if options.mode[:6]=="region":
		rgnp=[int(x) for x in options.mode[7:].split(',')]
		cen=(rgnp[0],rgnp[1])
		for i in range(nimg):
			a=EMData()
			a.read_image(args[0],i)
			b=a.get_clip(Region(rgnp[0]+a.get_xsize()/2-rgnp[2]/2,rgnp[1]+a.get_ysize()/2-rgnp[2]/2,rgnp[2],rgnp[2]))
			b.write_image(args[1],i)
	else:
		# copy the file with possible format conversion
		for i in range(nimg):
			a=EMData()
			a.read_image(args[0],i)
			a.write_image(args[1],i)
		

		
	cmplist=[(x,x+1) for x in range(nimg/2,nimg-1)]+[(x,x-1) for x in range(nimg/2,0,-1)]
	if options.twopass : cmplist+=cmplist
	ii=-1
	while ii< len(cmplist)-1:
		ii+=1
		i=cmplist[ii]
		if options.mode[:6]=="region" : inn=0
		else : inn=1
		
		# read local set of images to average for alignment
		if i[1]>i[0] :
			iml=EMData.read_images(args[inn],range(i[0]-options.localavg+1,i[0]+1))
		else :
			iml=EMData.read_images(args[inn],range(i[0],i[0]+options.localavg))
		for img in iml:
			img.process_inplace("normalize.edgemean")
		im1=iml[0].copy()
		for img in iml[1:]:
			im1+=img
		iml=None
		im1.process_inplace("normalize.edgemean")
		if options.highpass>0 :im1.process_inplace("filter.highpass.gauss",{"cutoff_abs":options.highpass})
		if (options.lowpass>0) : im1.process_inplace("filter.lowpass.gauss",{"cutoff_abs":options.lowpass})
		if options.localavg>1: im1.write_image("aliref.hed",i[0])
		
		im2=EMData()
		im2.read_image(args[inn],i[1])
		im2.process_inplace("normalize.edgemean")
		if options.highpass>0 : im2.process_inplace("filter.highpass.gauss",{"cutoff_abs":options.highpass})
		if (options.lowpass>0) : im2.process_inplace("filter.lowpass.gauss",{"cutoff_abs":options.lowpass})
		
		
		if options.mode=="modeshift" :
			dct=matrixalign(im1,im2,options.box,options.box+options.maxshift,debug=i[0]==63)
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
		elif options.mode[:6]=="region":
			if i[0]==nimg/2 :
				cen=(rgnp[0],rgnp[1])

			ref=im1.copy()
			
			mask=ref.copy_head()
			mask.to_one()
			mask.process_inplace("mask.sharp",{"outer_radius":rgnp[3]/2,"dx":cen[0],"dy":cen[1]})
			ref*=mask
			ref-=float(ref.get_attr("mean_nonzero"))
			ref*=mask
			
			im2.process_inplace("normalize")
			im2s=im2.copy()
			im2s.process_inplace("math.squared")
			
#			ref.write_image("dbug.hed",-1)
#			im2.write_image("dbug.hed",-1)
			ccf=ref.calc_ccf(im2,fp_flag.CIRCULANT)
			ccfs=mask.calc_ccf(im2s,fp_flag.CIRCULANT)	# this is the sum of the masked values^2 for each pixel center
			ccfs.process_inplace("math.sqrt")
			ccf/=ccfs
#			ccf.process_inplace("mask.sharp",{"outer_radius":(im1.get_xsize()-rgnp[2])/2})
			ccf.process_inplace("normalize")		# peaks relative to 1 std-dev
			ccf.process_inplace("mask.onlypeaks",{"npeaks":0})
			ccf.process_inplace("mask.sharp",{"outer_radius":options.maxshift})
			if options.nozero : ccf.set_value_at(ccf.get_xsize()/2,ccf.get_ysize()/2,0,0)

			if i[1] in range(72,77) : ccf.write_image("dbug.hed",-1)
			maxloc=ccf.calc_max_location()
			maxloc=(maxloc[0]-im1.get_xsize()/2,maxloc[1]-im1.get_ysize()/2)
			print maxloc
			
			out=im2.get_clip(Region(cen[0]-maxloc[0]-rgnp[2]/2+im2.get_xsize()/2,cen[1]-maxloc[1]-rgnp[2]/2+im2.get_ysize()/2,rgnp[2],rgnp[2]))
			cen=(cen[0]-maxloc[0],cen[1]-maxloc[1])
			out.write_image(args[1],i[1])
			print "%d.\t%d\t%d"%(i[1],cen[0],cen[1])
			continue			
		elif options.mode=="censym" :
			dct=matrixalign(im1,im2,options.box,options.box+options.maxshift,2)
			pairs=[]
			for x in range(3):
				for y in range(-2,3):
					if y<=0 and x==0 : continue
					a=dct[(x,y)]
					b=dct[(-x,-y)]
					if hypot(a[3]-b[3],a[4]-b[4])>7.0 : continue
					pairs.append((x,y,(a[3]+b[3])/2.0,(a[4]+b[4])/2.0,hypot(a[3]-b[3],a[4]-b[4])))
#					print "%d,%d\t%5.2f %5.2f\t%5.2f"%(pairs[-1][0],pairs[-1][1],pairs[-1][2],pairs[-1][3],pairs[-1][4])
			
			if len(pairs)==0 : 
				print "Alignment failed on image %d (%d)"%(i[1],i[0])
				if (i[1]>nimg/2) : cmplist[ii]=(i[0]-1,i[1])
				else : cmplist[ii]=(i[0]+1,i[1])
				if abs(i[0]-i[1])<5 : ii-=1
				continue
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
		im2.process_inplace("normalize")
		im2.write_image(args[1],i[1])
	
	print "Alignment Stage Complete"
	
	if options.tiltaxis!=400.0 :
		tiltaxis=(0,options.tiltaxis)
	else:
		print "Begin tilt axis search"
		
		# now we look for the common-line in the aligned images
		im1.read_image(args[1],0)
		sum=im1.do_fft()
		sum.to_zero()
		for i in range(nimg):
			a=EMData()
			a.read_image(args[1],i)
			a.process_inplace("mask.dampedzeroedgefill")
			a.process_inplace("normalize")
			a.process_inplace("mask.gaussian",{"outer_radius":a.get_xsize()/4})
			b=a.do_fft()
			b.process_inplace("complex.normpixels")
			sum+=b
		print "Phase average calculated"
			
		sum.ri2ap()
		curve=[]
		for angi in range(-90*4,90*4+1):
				ang=angi*pi/(180.0*4.0)
				v=0
				for r in range(a.get_xsize()/64,a.get_xsize()/2):
						v+=sum.get_value_at(2*int(r*cos(ang)),a.get_ysize()/2+int(r*sin(ang)))
				curve.append((v,ang*180.0/pi))
		tiltaxis=max(curve)
	
	print "Tilt axis at %1.2f degrees from x"%tiltaxis[1]

	for i in range(nimg):
		a=EMData()
		a.read_image(args[1],i)
		a.set_rotation(options.tilt*(i-nimg/2-1)*pi/180.0,0.0,-tiltaxis[1]*pi/180.0)
		a.set_attr("ptcl_repr",1)
		a.write_image(args[1],i)
		
#	sum=sum.do_ift()
#	sum.write_image("fft.mrc",0)
	
			
if __name__ == "__main__":
    main()
