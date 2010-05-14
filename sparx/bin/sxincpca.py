#!/usr/bin/env python
from EMAN2 import *
from sparx import *
from string import atoi, replace, split, atof
from sys import argv, exit
from user_functions import *
import  global_def
from global_def import *
from optparse import OptionParser
import sys

if len(argv) != 8:
	print "incpca.py prefix start end average.hdf output.hdf maskfile nvec "
	exit(-1)

if global_def.CACHE_DISABLE:
	from utilities import disable_bdb_cache
	disable_bdb_cache()
prefix = argv[1]
start  = atoi(argv[2])
end    = atoi(argv[3])
avg    = get_im(argv[4])
output = argv[5]
mask   = getImage(argv[6])
nvec   = atof(argv[7])

ana = Analyzers.get("pca_large", {"mask":mask, "nvec":nvec})
totnimg = 0

for i in xrange(start, end):
	filename = prefix + ('%04d.hdf'%i)
	print 'Loading file: ', filename
	nimg = EMUtil.get_image_count(filename)
	for j in xrange(nimg):
		img = get_im(filename, j)
		img -= avg

		ana.insert_image( img )
		totnimg += 1

		print totnimg, ' inserted'

		'''
		if totnimg%100==0:
			odd_var = odd_varer.getvar()
			eve_var = eve_varer.getvar()
			all_var = all_varer.getvar()
		
			odd_var.write_image( 'odd_' + output, iwrite )
			eve_var.write_image( 'eve_' + output, iwrite )
			all_var.write_image( output, iwrite )
			iwrite += 1  
			print 'ntot, ccc: %6d %10.3f' % ( totnimg, ccc(odd_var, eve_var,m) )  
		'''    

vecs = ana.analyze()
for iout in  xrange( len(vecs) ):
	vecs[iout].write_image( output, iout )

