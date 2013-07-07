#!/usr/bin/env python

# Author: Jesus Galaz, 02/Feb/2013, last update 26/June/2013
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


import os
from EMAN2 import *
import sys
import time
import numpy
import pylab
#from operator import itemgetter
from matplotlib.ticker import MaxNLocator
from pylab import figure, show	
import matplotlib.pyplot as plt
import math


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """UNDER DEVELOPOMENT. Extracts particles from each image in an aligned tilt series based on their position in the reconstructed tomogram."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	'''
	Parameters for adding ctf and noise
	'''
	parser.add_argument('--tiltseries',type=str,default='',help='File in .ali, .mrc or .hdf format of the aligned tiltseries.')
	parser.add_argument('--tiltangles',type=str,default='',help='File in .tlt or .txt format containing the tilt angle of each tilt image in the tiltseries.')
	parser.add_argument('--coords',type=str,default='',help='File in .txt format containing the coordinates of particles determined from the reconstructed tomogram of the supplied tiltseries.')
	parser.add_argument('--path',type=str,default='spt_subtilt',help='Directory to save the results.')
	parser.add_argument('--boxsize',type=int,default=128,help='Size of the 2D "tiles" or images for each particle from each image in the tiltseries.')
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument('--subset', type=int, default=0, help='''Specify how many sub-tiltseries (or particles) from the coordinates file you want to extract; e.g, if you specify 10, the first 10 particles will be boxed.\n0 means "box them all" because it makes no sense to box none''')
	#parser.add_argument('--tomogramthickness',type=int,default=None,help='Z dimension of the reconstructed tomogram.')
	parser.add_argument('--tomosides',type=str,default='',help='Comma separated values for the tomogram dimensions. Alternatively, provide the path to the tomogram itself through --tomogram.')
	parser.add_argument('--tomgoram',type=str,default='',help='Path to the tomogram.')
	parser.add_argument("--shrink", type=int,default=0,help="Optionally shrink the coordinates by a factor of --shrink=N to speed up the process. Might compromise accuracy if two points in the coordinates file are veyr close to eachother.")
	
	(options, args) = parser.parse_args()
	logger = E2init(sys.argv, options.ppid)
	return
	
	#def subtiltextractor(parameters):
	
	if not options.tiltseries or not options.tiltangles or not options.coords or not options.tomogramthickness:
		print "ERROR: You must provide ALL of the following options: --tiltseries, --tiltangles, --coords and --tomogramthickness."
		sys.exit()
		
	if options.subset:
		if options.subset > nptcls:
			print "WARNING: The total amount of lines in the coordinates files is LESS than the subset of particles to box you specified; therefore, ALL particles will be extracted."
		else:
			nptcls=options.subset
	
	print "The size of the set of sub-tiltseries to extract is", nptcls
	
	#k=-1
	#name = options.output
	
	serieshdr = EMData(options.tiltseries,0,True)
	nslices = serieshdr['nz']
	nx = serieshdr['nx']
	ny = serieshdr['ny']
	
	anglesfile = open(options.tiltangles,'r')				#Open tilt angles file
	alines = anglesfile.readlines()							#Read its lines
	anglesfile.close()										#Close the file
	
	tiltangles = [ alines[i].replace('\n','') for i in range(len(alines)) ]	#Eliminate trailing return character, '\n', for each line in the tiltangles file
	
	#for i in range(len(alines)):
	#	alines[i]=alines[i].replace('\n')					#Eliminate
	
	ntiltangles = len(tiltangles)
	
	if int(nslices) != int(ntiltanges):
		print """ERROR: The tiltangles file doesn't seem to correspond to the tiltseries provided.
				The number of images in --tiltseries (z dimension of MRC stack) must be equal to the number
				of lines in --tiltangles."""
		sys.exit()
	
	
	#"""
	#You do not need to keep track of the mathematics of tilting and rotating and find the correspondence between
	#tomogram and each image in the tilt series.
	#Instead, let's make a simple 3D model, containing a bright dot at the position of each particle.
	#Then, rotate that using the known tilt angles, generate a projection for each, and find the dots (maxima) in each.
	#Use the location of these maxima, to extract each particle from each image in the tilt series.
	#"""
	
	if options.tomosides:							#Read tomogram dimensions.
		sides=options.tomosides.split(',')
		tomox = sides[0]
		tomoy = sides[1]
		tomoz = sides[2]
	
	if options.tomogram:
		tomohdr = EMData(options.tomogram,0,True)	#Read tomogram dimensions from tomogram header, if the tomogram is provided.
		tomox = tomohdr['nx']
		tomoy = tomohdr['ny']
		tomoz = tomohdr['nz']

	if options.shrink:								#The 'MODEL' to build for the coordinates need not use the full size of the tomogram.
		tomox = tomox/4.0
		tomoy = tomoy/4.0
		tomoz = tomoz/4.0
	
	tomovol = EMData(tomox,tomoy,tomoz)				#Create empty volume for the MODEL to build
	tomovol.to_zero()								#Make sure it's empty
	
	cfile = open(options.coords,'r')				#Open coordinates file
	clines = cfile.readlines()						#Read its lines
	cfile.close()									#Close the file
	
	'''
	Iterate over the lines of the coordinates file.
	Some people might manually make ABERRANT coordinates files with commas, tabs, or more than once space in between coordinates.
	Each line needs to be parsed.
	'''
	p=1
	for line in clines:
		line =line.replace(", ",' ')	
		line = line.replace(",",' ')
		line = line.replace("x",'')
		line = line.replace("y",'')
		line = line.replace("z",'')
		line = line.replace("=",'')
		line = line.replace("_",' ')
		line = line.replace("\n",' ')
		line = line.replace("\t",' ')
		line = line.replace("  ",' ')
		
		line = line.split()		
	
		xc = int(clines[i][0])				#Determine x y z coordinates for each line
		yc = int(clines[i][1])
		zc = int(clines[i][2])
		
		for k in range(len(tiltangles)):
			
			xcshift = xc - tomox/2.0
			
			#RY = math.array( [ [math.cos(tilt),0, -1 * math.sin(tilt), 0], [0,1,0,0], [math.sin(tilt), 0, math.cos(tilt), 0], [0,0,0,1] ] )
			#zt = zc* math.cos(tiltangles[k]) - xcshift * math.sin(tiltangles[k])
			
			xtshift = zc* math.sin(tiltangles[k) + xcshift * math.cos(tiltangles[k])
			yt = yc
			
			xt = xtshift + tomox/2.0
			
			r = Region( (2*xt-options.boxsize)/2, (2*yt-options.boxsize)/2, k, options.boxsize, options.boxsize, k+1)
			e = EMData()
			e.read_image(options.tiltseries,0,False,r)
			e.write_image('subtilt_' + str(p) + '.hdf',p-1)
						
		p+=1
		
		#		(cos q  0  -sin q   0)
		#Ry(q) = (0      1    0      0)
		#        (sin q  0  cos q    0)
		#        (0      0    0     1) 

		
		
		
		if options.shrink:					#Shrink them if the model will be smaller than the original tomgoram
			xc = xc/4.o
			yc = yc/4.0
			zc = zc/4.0
			
		tomovol.set_value_at(xc,yc,zc,1)	#Set the value of the center pixel where any particles were picked to 1.
		
	
	for tilt in tiltangles:
		t=Transform({'type':'eman','az':0,'phi':0,'alt':tilt})
		prj = tomovo.project('standard',t)
		
		for i in range(len(clines)):
			
		
	
	
	
	
	
	
	
	'''
	
	"""
	Old mathematical approach
	"""	
	for j in range(nslices):
		
		2Dregion = Region(0,0,j,nx,ny,j+1)
		#a = EMData(s)
		slice=EMData()
		slice.read_image(options.tiltseries,0,False,2Dregion)
		tiltaxisx = int(nx)/2
		
		alpha = tiltangles[j]
		#k = 1
		for i in range(nptcls):
			#Some people might manually make ABERRANT coordinates files with commas, tabs, or more than once space in between coordinates
			clines[i] = clines[i].replace(", ",' ')	
			clines[i] = clines[i].replace(",",' ')
			clines[i] = clines[i].replace("x",'')
			clines[i] = clines[i].replace("y",'')
			clines[i] = clines[i].replace("z",'')
			clines[i] = clines[i].replace("=",'')
			clines[i] = clines[i].replace("_",' ')
			clines[i] = clines[i].replace("\n",' ')
			clines[i] = clines[i].replace("\t",' ')
			clines[i] = clines[i].replace("  ",' ')
			clines[i] = clines[i].split()		
		
			xc = int(clines[i][0])
			yc = int(clines[i][1])
			zc = int(clines[i][2])
		
			y2 = y
			
			xp2ta = tiltaxisx - xc
			zp = zc - options.thickness/2
			
			if alpha > 0:
				zp = -1 * zp
			
			dxa = (xp2ta + zp) * cos(alpha)
			
			x2 = tiltaxis - la
			
			r = Region((x2 - boxsize)/2,(y2 - boxsize)/2, boxsize, boxsize)
        		e = EMData()
			e.read_image(s,0,False,r)
			
			name = 'particle#' + str(k).zfill(len(pcoords)) + '_slice' + str(j).zfill(len(pcoords)) + '.mrc'
	'''
		
if __name__ == '__main__':
	
	main()
