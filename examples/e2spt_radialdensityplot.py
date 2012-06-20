#!/usr/bin/env python
#
# Author: Jesus Galaz, 06/05/2012
# Copyright (c) 2011 Baylor College of Medicine
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
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from pylab import figure, show	

def main():

	progname = os.path.basename(sys.argv[0])
	usage = """Aligns a 3d volume to another by executing e2classaverage3d.py and then calculates the FSC between them by calling e2proc3d.py . It returns both a number for the resolution based on the FSC0.5 
	criterion(on the screen) and a plot as an image in .png format."""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--vols", type=str, help="Volume whose radial density plot you want to compute. For multiple volumes, either provide them as an .hdf stack, or separate them by commas --vols=first.hdf,second.hdf,etc...", default=None)
	parser.add_argument("--output", type=str, help="Name for the output .png and .txt files that contain the plots and the numeric values for them. Must be specified if --singleplot is on.", default=None)
	parser.add_argument("--mode", type=str, help="provide --mode=x, y, or z to get the average density per slice in the indicated direction. --mode=cylinder for concentric cylindrical shell; default is --mode=sphere", default='sphere')
	parser.add_argument("--fixedcylinderheight", type=int, help="Works only if --mode=cylinder, and keeps the height of the cylinder at a constant value, while varying the radius.", default=0)

	parser.add_argument("--mask",type=str,help="Mask processor applied to volumes before alignment. Default is mask.sharp:outer_radius=-2", default="mask.sharp:outer_radius=-2")
	parser.add_argument("--normproc",type=str,help="Normalization processor applied to volumes before alignment. Default is to use normalize.mask. If normalize.mask is used, results of the mask option will be passed in automatically. If you want to turn this option off specify 'None' ", default=None)
	parser.add_argument("--preprocess",type=str,help="Any processor (as in e2proc3d.py) to be applied to each volume prior to radial density plot computation.", default=None)
	parser.add_argument("--lowpass",type=str,help="A lowpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to radial density plot computation.", default=None)
	parser.add_argument("--highpass",type=str,help="A highpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to radial density plot computation.", default=None)	
	parser.add_argument("--shrink", type=int,default=1,help="Optionally shrink the input volumes by an integer amount.")	
	#parser.add_argument("--apix", type=float, help="Provide --apix to overrride the value found in the volumes' header paramter.", default=0)
	parser.add_argument("--singleplot", action="store_true",default=False,help="Plot all the Radial Density Profiles of the volumes provided in one single plot.")	
	parser.add_argument("--threshold", action="store_true",default=False,help="If on, this will turn all negative pixel values into 0.")	

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
	
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
		
	names = options.vols
	names = names.split(',')
	
	finalvalues = {}
	for i in names:
		
		n = EMUtil.get_image_count(i)
		
		stackvalues = []
		print "The stack %s has %d images in it" % ( i, n ) 
		for j in range(n):
			ptcl = EMData(i,j)
			#radius = a.get_xsize()/2
			values = calcvalues(ptcl,options)	
			stackvalues.append(values)
		finalvalues.update({i:stackvalues})	
	
	print "\n\nfinal values are", finalvalues
	
	if options.singleplot and len(names) > 1:
		for i in names:
			if len(finalvalues[i]) > 1:
				print "ERROR: You can plot RD profiles for multiple particles IN ONE PLOT if each individual .hdf file has one particle only, or if you supply ONE stack only with all the particles in it."
				print "In this case, you've supplied %d files, and the file %s has %d particles in it" % (len(names), i, len(finalvalues[i]))
				sys.exit()

		if not options.output:
			print "ERROR: You must supply and output name if you want to plot multiple RD profiels from different .hdf files into one plot."
			sys.exit()
		else:	
			plotname = options.output
			if '.png' not in plotname:
				plotname += '.png'

			plt.title("Spherical radial density plot")
			plt.ylabel("Density (arbitrary units)")
			plt.xlabel("Radius (angstroms)")

			if options.mode == 'x':
				plt.title("Density plot of slices along x-axis")
				plt.xlabel("X (angstroms)")

			if options.mode == 'y':
				plt.title("Density plot of slices along y-axis")
				plt.xlabel("Y (angstroms)")

			if options.mode == 'z':
				plt.title("Density plot of slices along z-axis")
				plt.xlabel("Z (angstroms)")
				
			if options.mode == 'cylinder':
				plt.title("Density plot of concentric cylyndrical shells")
				plt.xlabel("Radius (angstroms)")
			
			for i in names:
				apix = EMData(i,0,True)['apix_x']
	
				x = range(len(values))
				for j in range(len(x)):
					x[j] = int(round(x[j] * apix))		
			
				values = finalvalues[i][0]
				txtname = plotname.replace('.png', '_' + str(j).zfill(len(names)) + '.txt') 
				f = open(txtname,'w')
				lines = []
				for value in values:
					line = str(value) + '\n'
					lines.append(line)
				f.writelines(lines)
				f.close()
				
				plt.plot(x,values, linewidth=2)

			p = plt.gca()
			plt.savefig(plotname)
			plt.clf()

	else:
		for i in names:
			plotname = ''
			if options.output:
				plotname = options.output
			else:
				plotname = i.replace('.hdf','.png')
			
			if plotname and '.png' not in plotname:
				plotname += '.png'
				
			for j in range(len(finalvalues[i])):
				apix = EMData(i,j,True)['apix_x']

				plotname = plotname.replace('.png','_' + str(j).zfill(len(finalvalues[i])) + '.png')
				txtname = plotname.replace('.png','.txt')
				values = finalvalues[i][j]
				
				x=range(len(values))
				for i in range(len(x)):
					x[i] = int(round(x[i] * apix))
				
				plt.title("Spherical radial density plot")
				plt.ylabel("Density (arbitrary units)")
				plt.xlabel("Radius (angstroms)")

				if options.mode == 'x':
					plt.title("Density plot of slices along x-axis")
					plt.xlabel("X (angstroms)")
				
				if options.mode == 'y':
					plt.title("Density plot of slices along y-axis")
					plt.xlabel("Y (angstroms)")
				
				if options.mode == 'z':
					plt.title("Density plot of slices along z-axis")
					plt.xlabel("Z (angstroms)")
					
				if options.mode == 'cylinder':
					plt.title("Density plot of concentric cylyndrical shells")
					plt.xlabel("Radius (angstroms)")

				plt.plot(x,values,linewidth=2)
				p = plt.gca()
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
	if options.shrink !=None and options.shrink>1 :
		a=a.process("math.meanshrink",{"n":options.shrink})
	
	if options.threshold:
		a=a.process("threshold.belowtozero",{"minval":0.0})

	if options.mode == 'sphere':
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


def direction(a,options):
	values = []
	mask = EMData(a['nx'],a['ny'],a['nz'])
	mask.to_one()
	
	rng = a['nx']
	if options.mode == 'y':
		rng == a['ny']
	
	if options.mode == 'z':
		rng == a['nz']
	
	print "The mode is", options.mode
	print "And the range for values calculation is", rng
	
	for i in xrange(0,rng):
		maskslice = mask
		if options.mode == 'x':
			maskslice = mask.process("mask.zeroedge3d",{'x0':i,'x1':a['nx'] -i -1,'y0':0,'y1':0,'z0':0,'z1':0})
		
		if options.mode == 'y':
			maskslice = mask.process("mask.zeroedge3d",{'x0':0,'x1':0,'y0':i,'y1':a['ny'] -i -1,'z0':0,'z1':0})
		
		if options.mode == 'z':
			maskslice = mask.process("mask.zeroedge3d",{'x0':0,'x1':0,'y0':0,'y1':0,'z0':i,'z1':a['nz'] -i -1})
		
		b = a.copy()
		b.mult(maskslice)
		value = b ['mean_nonzero']
		values.append(value)
		
	return(values)

if __name__ == '__main__':
	main()
