#!/usr/bin/env python

from EMAN2 import *
import os, sys
from sys import argv

stem = argv[1]
process =  argv[2]



c = os.getcwd()
findir = os.listdir(c)
imgextensions = ['mrcs','MRCS','hdf','HDF','ali','ALI','st','ST','mrc','MRC','tif','TIF','DM3','dm3']

for f in findir:
	dimension=2
	program='e2proc2d.py'

	if stem in f:
		format = f.split('.')[-1]
		if format in imgextensions:
			print "\nprocessing file", f
			hdr = EMData(f,0,True)
			nz = hdr['nz']

			if nz>1:
				dimension=3
			
			print "dimensionality of image is",dimension

			if format=='mrcs' or format=='MRCS':
				dimension=2
				print "but this is an mrcs file, so dimensionality", dimension
			
			if dimension == 3:
				program='e2proc3d.py'
			
			output = f.split('.' + format)[0] + '_proc.' + format
			cmd = program + ' ' + f + ' ' + output + ' --process ' + process
			print "\ncmd to run is",cmd
			os.popen( cmd )
		else:
			print "\nERROR: skipping invalid file", f