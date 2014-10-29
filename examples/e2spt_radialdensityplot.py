#!/usr/bin/env python
#
# Author: Jesus Galaz, 06/05/2012 - Last change 12/17/2012
# Copyright (c) 2011 Baylor College of Medicine
#
# ******** CHANGES FOR DONGHUA INCLUDED *********
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
	
import os, sys, commands
from EMAN2 import *
from pylab import figure, show	
import math

def main():

	progname = os.path.basename(sys.argv[0])
	usage = """This program allows you to examine density variations along one or more given volume (at the same time).
				It calculates the mean intensity either for slices (planes) along any the three cartesian axes (X, Y or Z), or for radial consecutive shells of increasing radius, 
				or for cylindrical shells of varying or fixed height, starting from the center of the volume. 
				All mean density values are saved to .txt files, and plots are produced with them and saved as .png images. In fact, to compare different volumes you can plot all curves in a single plot."""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--vols", type=str, help="""Volume whose radial density plot you 
		want to compute. For multiple volumes, either provide them as an .hdf stack, or 
		separate them by commas --vols=first.hdf,second.hdf,etc...""", default=None)
	
	parser.add_argument("--mode", type=str, help="""provide --mode=x, y, or z to get the 
		average density per slice in the indicated direction. 
		--mode=cylinder for concentric cylindrical shell; default is --mode=sphere.
		For MULTIPLE modes, separate them by commas, for example --mode=x,y,z,cylinder""", default='sphere')
	
	parser.add_argument("--fixedcylinderheight", type=int, help="""Works only if --mode=cylinder, 
		and keeps the height of the cylinder at a constant value, while varying the radius.""", default=0)

	
	parser.add_argument("--mask",type=str,help="""Mask processor applied to volumes before alignment. 
		Default is None.""", default=None)
	
	parser.add_argument("--normproc",type=str,help="""Normalization processor applied to 
		volumes before computing density values. Default is None.
		If normalize.mask is used, results of the mask option will be passed in automatically.""", default=None)
	
	parser.add_argument("--preprocess",type=str,help="""Any processor (as in e2proc3d.py) 
		to be applied to each volume prior to radial density plot computation.""", default=None)
	
	parser.add_argument("--lowpass",type=str,help="""A lowpass filtering processor (as in e2proc3d.py) 
		to be applied to each volume prior to radial density plot computation.""", default=None)
	
	parser.add_argument("--highpass",type=str,help="""A highpass filtering processor (as in e2proc3d.py) 
		to be applied to each volume prior to radial density plot computation.""", default=None)	
	
	parser.add_argument("--shrink", type=int,default=1,help="""Optionally shrink the input 
		volumes by an integer amount.""")	
	
	#parser.add_argument("--apix", type=float, help="Provide --apix to overrride the value found in the volumes' header paramter.", default=0)
	
	parser.add_argument("--singleplotperfile", action="store_true",default=False,help="""
		Plot all the Radial Density Profiles of the volumes provided in each .hdf stack in 
		one single plot.""")	
	
	parser.add_argument("--singlefinalplot", action="store_true",default=False,help="""Plot 
		all the Radial Density Profiles of the volumes provided in all .hdf stacks in one 
		FINAL single 'master' plot.""")	
	
	parser.add_argument("--threshold", action="store_true",default=False,help="""If on, 
		this will turn all negative pixel values into 0.""")	
	
	parser.add_argument("--normalizeplot", action="store_true",default=False,help="""This 
		will make the maximum density in each plot or curve equal to 1.""")
	
	parser.add_argument("--ppid", type=int, help="""Set the PID of the parent process, 
		used for cross platform PPID""",default=-1)
	
	parser.add_argument("--verbose", "-v", default=0, help="""Verbose level [0-9], higner 
		number means higher level of verboseness""",dest="verbose", action="store", metavar="n",type=int)
	
	(options, args) = parser.parse_args()
	
	import matplotlib.pyplot as plt
	from matplotlib.ticker import MaxNLocator
	
	logger = E2init(sys.argv, options.ppid)

	if options.normproc: 
		options.normproc=parsemodopt(options.normproc)
	
	if options.mask: 
		options.mask=parsemodopt(options.mask)
	
	if options.preprocess: 
		options.preprocess=parsemodopt(options.preprocess)
		
	if options.lowpass: 
		options.lowpass=parsemodopt(options.lowpass)
	
	if options.highpass: 
		options.highpass=parsemodopt(options.highpass)
			
	files = options.vols
	files = files.split(',')
	
	for i in xrange(0,len(files)):
		for j in range(i+1,len(files)):
			if files[i] == files[j]:
				print "ERROR: You have supplied a file twice, see", files[i],files[j]
				sys.exit()
	modes=options.mode.split(',')
	
	for m in modes:
		options.mode=m
		modetag = '_MODE' + m
		finalvalues = {}
		imgsperstack=1
		
		names=[]
		
		finalvalues=[]
		for i in files:	
			n = EMUtil.get_image_count(i)
			
			print "The stack %s has %d images in it" % ( i, n ) 
			
			kk=0
			stack={}
			stackvalues = []
			suffix = modetag
			for j in range(n):
				ptcl = EMData(i,j)
				if n > 1:
					suffix = modetag + str(kk).zfill(len(str(n)))
				
				values = calcvalues(ptcl,options)
					
				print "For file %s img number %d the max is %f" %(i,j,max(values))
				
				if options.normalizeplot:
					
					minv = min(values)
					for v in range(len(values)):
						values[v] = values[v] - minv
					
					maxv = max(values)
					for v in range(len(values)):	
						values[v] = values[v]/maxv	
					print "Therefore, max is", max(values)
				
				id=i.replace('.',suffix + '.')
				stackvalues.append([id,values])
				kk+=1
			stack.update({i:stackvalues})
			finalvalues.append(stack)
			 		
		plotname = 'finalplot_MODE' + m + '.png'
		fileid=''
		
		cc=0
		for ele in finalvalues:
			i = ele.keys()[0]
			key = i
			
			n = EMUtil.get_image_count(i)
			
			if options.singleplotperfile:
				fileid=i.split('.')[0]	
				plotname = fileid + modetag + '.png'
				
			kk=0
			for f in range(n):
				apix = EMData(i,f,True)['apix_x']
				
				values = ele[key][f][1]
				id = ele[key][f][0]
				
				x = range(len(values))				
					
				for j in range(len(x)):
					x[j] = int(round(x[j] * apix))				
				txtname = i.split('.')[0] + modetag + str(kk).zfill(len(str(n))) + '.txt'
				txtf = open(txtname,'w')
				lines = []
					
				for v in range(len(values)):
					line = str(v) +  ' ' + str(values[v]) + '\n'
					lines.append(line)
				txtf.writelines(lines)
				txtf.close()
				
				plt.plot(x,values,linewidth=2)
				
				if not options.singleplotperfile and not options.singlefinalplot:
					#plotname=i.split('.')[0]+str(kk).zfill(len(str(n))) + '.png'
					plotname=id.split('.')[0] + '.png'
					fileid = plotname.split('.')[0]
				
				if options.mode == 'sphere':
					plt.title("Spherical radial density plot " + fileid)
					plt.xlabel("Radius (angstroms)")
					plt.ylabel("Density (arbitrary units)")
				
				if options.mode == 'x':
					plt.title("Density plot of slices along x-axis "+ fileid)
					plt.xlabel("X (angstroms)")
					plt.ylabel("Density (arbitrary units)")

				if options.mode == 'y':
					plt.title("Density plot of slices along y-axis "+ fileid)
					plt.xlabel("Y (angstroms)")
					plt.ylabel("Density (arbitrary units)")

				if options.mode == 'z':
					plt.title("Density plot of slices along z-axis "+ fileid)
					plt.xlabel("Z (angstroms)")
					plt.ylabel("Density (arbitrary units)")
				
				if options.mode == 'cylinder':
					plt.title("Density plot of concentric cylyndrical shells "+ fileid)
					plt.xlabel("Radius (angstroms)")
					plt.ylabel("Density (arbitrary units)")			
								
				if not options.singleplotperfile and not options.singlefinalplot:
					plt.savefig(plotname)
					plt.clf()
				else:
					pass
				kk+=1
				cc+=1
			if options.singleplotperfile and not options.singlefinalplot:
				plt.savefig(plotname)
				plt.clf()
			else:
				pass
	
		if options.singlefinalplot:
			plt.savefig(plotname)
			plt.clf()
	return()				
				

def calcvalues(a,options):
	# Make the mask first, use it to normalize (optionally), then apply it 
	mask=EMData(a["nx"],a["ny"],a["nz"])
	mask.to_one()

	if options.mask != None:
		mask.process_inplace(options.mask[0],options.mask[1])

	# normalize
	if options.normproc != None:
		if options.normproc[0]=="normalize.mask": 
			options.normproc[1]["mask"]=mask
		a.process_inplace(options.normproc[0],options.normproc[1])

	a.mult(mask)

	if options.normproc != None:
		if options.normproc[0]=="normalize.mask": 
			options.normproc[1]["mask"]=mask
		a.process_inplace(options.normproc[0],options.normproc[1])

	a.mult(mask)

	# preprocess
	if options.preprocess != None:
		a.process_inplace(options.preprocess[0],options.preprocess[1])

	# lowpass
	if options.lowpass != None:
		a.process_inplace(options.lowpass[0],options.lowpass[1])

	# highpass
	if options.highpass != None:
		a.process_inplace(options.highpass[0],options.highpass[1])

	# Shrink
	if options.shrink>1 :
		shrinkfactor = options.shrink
		x = a['nx']
		y = a['ny']
		z = a['nz']
		f = min(x,y,z)
		if shrinkfactor > f:
			print """You have supplied a shrinkfactor that is larger than the smallest dimensions of your image.
				Therefore, the shrinkfactor will be changed to""", f
			shrinkfactor = f
			
		a=a.process("math.meanshrink",{"n":shrinkfactor})
	
	if options.threshold:
		a=a.process("threshold.belowtozero",{"minval":0.0})

	if options.mode == 'sphere':
		print "I will calculate the radial density"
		values = a.calc_radial_dist(a['nx']/2, 0, 1, 1)
		return(values)
	
	elif options.mode == 'cylinder':
		values = cylinder(a,options)
		return(values)
		
	elif options.mode == 'x' or options.mode == 'y' or options.mode == 'z':
		values = direction(a,options)
		return(values)


def cylinder(a,options):
	values = []
	mask = EMData(a['nx'],a['ny'],a['nz'])
	mask.to_one()
	
	for i in xrange(1,a['nx']/2):
		heightout = i
		heightin = heightout-1
		radius = i
		if options.fixedcylinderheight:
			heightout = options.fixedcylinderheight
			heightin = heightout
		#print "Heighout, heightin and radius are", heightout, heightin, radius
		maskout = mask.process("testimage.cylinder",{'height':heightout,'radius':radius})
		maskin = mask.process("testimage.cylinder",{'height':heightin,'radius':radius-1})
		
		finalmask = maskout - maskin
		
		b = a.copy()
		#b.mult(maskout)
		#b.mult(maskin)
		b.mult(finalmask)
		value = b ['mean_nonzero']
		values.append(value)
		
	return(values)


def calc_dimension(a):
	dimensionality = 0
	nx=a['nx']
	ny=a['ny']
	nz=a['nz']
	if nx == 1 and ny == 1 and nz == 1:
		pass
	if (nx == 1 and ny == 1 and nz > 1) or (nx == 1 and nz == 1 and ny > 1) or (ny == 1 and nz == 1 and nx > 1): 
		dimensionality = 1
	if (nx == 1 and ny > 1 and nz > 1) or (ny == 1 and nx > 1 and ny > 1) or (nz == 1 and nx > 1 and ny > 1):
		dimensionality = 2
	if nx >1 and ny > 1 and nz > 1:
		dimensionality = 3
	
	return (dimensionality)


def direction(a,options):
	values = []
	dimensionality = calc_dimension(a)
	if dimensionality < 2:
		if dimensionality == 0:
			print "Your image is a point. There's nothing to analyze. Perhaps you over shrunk"
			sys.exit()
	x = a['nx']
	y = a['ny']
	z = a['nz']
	mask = EMData(a['nx'],a['ny'],a['nz'])
	mask.to_one()
	#print "The size of the image is", x,y,z
	rng = 0
	if options.mode == 'x':
		rng = x
	if options.mode == 'y':
		rng = y
	if options.mode == 'z':
		rng = z
	
	#print "The mode is", options.mode
	#print "and it should be equal to y see", 'y' == options.mode
	#print "And the range for values calculation is", rng
	
	for i in xrange(0,rng):
		#print "\nFor slice", i
		maskslice = mask
		#print "I will mask the image, whose dimensionality is", dimensionality
		if options.mode == 'x' and x > 1:
			if dimensionality == 3:
				maskslice = mask.process("mask.zeroedge3d",{'x0':i,'x1':a['nx'] -i -1,'y0':0,'y1':0,'z0':0,'z1':0})
			elif dimensionality == 2:
				maskslice = mask.process("mask.zeroedge2d",{'x0':i,'x1':a['nx'] -i -1,'y0':0,'y1':0})
				
		if options.mode == 'y' and y > 1:
			if dimensionality == 3:
				maskslice = mask.process("mask.zeroedge3d",{'x0':0,'x1':0,'y0':i,'y1':a['ny'] -i -1,'z0':0,'z1':0})
			elif dimensionality == 2:
				maskslice = mask.process("mask.zeroedge2d",{'x0':0,'x1':0,'y0':i,'y1':a['ny'] -i -1})
			#print "The mask in Y was computed."
		if options.mode == 'z' and z > 1:
			if dimensionality == 3:
				maskslice = mask.process("mask.zeroedge3d",{'x0':0,'x1':0,'y0':0,'y1':0,'z0':i,'z1':a['nz'] -i -1})
			else:
				print "ERROR: It makes no sense to look for density variations across z y a 2D image"
				sys.exit()	
		b = a.copy()
		b.mult(maskslice)
		#print "The mask was applied.\n"
		value = b ['mean_nonzero']
		values.append(value)
		
	return(values)

if __name__ == '__main__':
	main()
