#!/usr/bin/env python

# Author: Jesus Galaz, 02/02/2013, last update 02/01/2013
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


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Extracts particles from each image in an aligned tilt series based on their position in the reconstructed tomogram."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	'''
	Parameters for adding ctf and noise
	'''
	parser.add_argument('--tiltseries',type=str,default=None,help='File in .ali, .mrc or .hdf format of the aligned tiltseries.')
	parser.add_argument('--tiltangles',type=str,default=None,help='File in .tlt or .txt format containing the tilt angle of each tilt image in the tiltseries.')
	parser.add_argument('--coords',type=str,default=None,help='File in .txt format containing the coordinates of particles determined from the reconstructed tomogram of the supplied tiltseries.')
	parser.add_argument('--path',type=str,default='spt_subtilt',help='Directory to save the results.')
	parser.add_argument('--boxsize',type=int,default=128,help='Size of the 2D "tiles" or images for each particle from each image in the tiltseries.')
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument('--subset', type=int, default=0, help='''Specify how many sub-tiltseries (or particles) from the coordinates file you want to extract; e.g, if you specify 10, the first 10 particles will be boxed.\n0 means "box them all" because it makes no sense to box none''')
	parser.add_argument('--tomogramthickness',type=int,default=None,help='Z dimension of the reconstructed tomogram.')

	
	(options, args) = parser.parse_args()
	logger = E2init(sys.argv, options.ppid)
	return
	
	#def subtiltextractor(parameters):
	
	if not options.tiltseries or not options.tiltangles or not options.coords or not options.tomogramthickness:
		print "ERROR: You must provide ALL of the following options: --tiltseries, --tiltangles, --coords and --tomogramthickness."
		sys.exit()
		
	if options.subset:
		if options.subset > nptcls:
			print "WARNING: The total amount of lines in the coordinates files is LESS that the subset of particles to box you specified; therefore, ALL particles will be extracted"
		else:
			nptcls=options.subset
	
	print "The size of the set of sub-tiltseries to extract is", nptcls
	
	#k=-1
	#name = options.output
	
	clines = open(options.coords,'r').readlines()
	nptcls = len(clines)
	
	serieshdr = EMData(options.tiltseries,0,True)
	nslices = serieshdr['nz']
	nx = serieshdr['nx']
	ny = serieshdr['ny']
	
	tiltangles = open(options.tiltangles,'r').readlines()
	ntiltangles = len(tiltanges)
	
	if int(nslices) != int(ntiltanges):
		print """ERROR: The tiltangles file doesn't seem to correspond to the tiltseries provided.
				The number of images in --tiltseries, or the its z dimension, must be equal to the number
				of lines in --tiltangles."""
		sys.exit()
		
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
		
if __name__ == '__main__':
	
	main()
