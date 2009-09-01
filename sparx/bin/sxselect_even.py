#!/usr/bin/env python

from EMAN2 import *
from sparx import *
from string import atoi, split, atof
from sys import argv, exit
import global_def

if len(argv) != 5:
	print "Usage: select_even.py prj_stack stat.txt each_prj prj_out"
	exit(-1)

if global_def.CACHE_DISABLE:
	from utilities import disable_bdb_cache
	disable_bdb_cache()

stack_in = argv[1]
fstat    = open( argv[2], 'r' )
each_prj = atoi( argv[3] )
stack_ot = argv[4]

def cmp_defocus(i1, i2):
	return i1[1] > i2[1]

itotal = 0

line = fstat.readline()
while len(line)>0:
	items = split(line)
	nprj  = atoi( items[5] )

	prjinfo = []

	for i in xrange(nprj):
		prj_line = fstat.readline()
		prj_item = split( prj_line )
		prj_id   = atoi( prj_item[0] )
		defocus  = atof( prj_item[1] )
		prjinfo.append( [prj_id, defocus] )

	prjinfo.sort( lambda x,y : cmp(y[1],x[1]) )

	if each_prj > nprj : curt_prj = nprj
	else:                curt_prj = each_prj    


	for i in xrange(curt_prj):
		prj_id = prjinfo[i][0]
		prj    =  get_im( stack_in, prj_id )

		defocus1 = prjinfo[i][1]
		defocus2 = prj.get_attr( 'defocus' )
		assert defocus1==defocus2

		print 'prj %6d wrote to %6d, defocus %10.3f' % (prj_id, itotal, defocus1)

		prj.write_image( stack_ot, itotal )
		itotal += 1

	line = fstat.readline()
