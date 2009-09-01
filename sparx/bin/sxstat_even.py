#!/usr/bin/env python
from EMAN2 import *
from sparx import *
import global_def
from math import pi,sin,cos

angle_to_rad = pi/180.0

def getvec( phi, tht ):

	if tht > 180.0:
		tht = tht - 180.0
		phi = phi + 180.0
	elif tht > 90.0:
		tht = 180.0 - tht
		phi = phi + 180.0

	assert tht <=90.0

	tht *= angle_to_rad
	phi *= angle_to_rad

	x = sin(tht)*cos(phi) 
	y = sin(tht)*sin(phi)
	z = cos(tht)

	return (x,y,z)

def nearest_ang( vecs, phi, tht ) :
	vec = getvec( phi, tht )

	best_s = -1.0
	best_i = -1

	for i in xrange( len(vecs) ):
		s = vecs[i][0]*vec[0] + vecs[i][1]*vec[1] + vecs[i][2]*vec[2]
		if s > best_s:
			best_s = s
			best_i = i

	return best_i

from sys import argv, exit
from string import atof

if len(argv)!=4:
	print "Usage: sxstat_even.py prj_stack delta out.txt"
	exit(-1)
if global_def.CACHE_DISABLE:
	from utilities import disable_bdb_cache
	disable_bdb_cache()

prj_stack = argv[1]
eve_angls = even_angles( atof(argv[2]), 0.0, 90.0, method='P' )
of = open( argv[3], 'w' )

eve_vecs = [None]*len(eve_angls)

for i in xrange( len(eve_angls) ):
	eve_vecs[i] = getvec( eve_angls[i][0], eve_angls[i][1] )


nimage = EMUtil.get_image_count( prj_stack )

hist = [None]*len(eve_vecs)

phis = []
thts = []
prj = EMData()
for i in xrange(nimage):

	prj.read_image( prj_stack, i, True )

	phi = prj.get_attr( 'phi' )
	tht = prj.get_attr( 'theta' )

	phis.append( phi )
	thts.append( tht )

	angid  = nearest_ang( eve_vecs, phi, tht )

	if hist[angid] is None:
		hist[angid] = []

	defocus = prj.get_attr( 'defocus' )
	hist[angid].append( [i, defocus] )

	print "prj #%6d %10.3f %10.3f close to ori #%6d %10.3f %10.3f" % (i, phi, tht, angid, eve_angls[angid][0], eve_angls[angid][1])


for i in xrange( len(hist) ):
	if hist[i] is None:
		hist[i] = []

	of.write( "ori %6d %10.3f %10.3f nprj %6d\n" % (i, eve_angls[i][0], eve_angls[i][1], len(hist[i]) ) )

	for j in xrange( len(hist[i]) ):
		imgid = hist[i][j][0]
		defocus = hist[i][j][1]
		of.write( "    %6d %10.3f %10.3f %10.3f\n" % (imgid, defocus, phis[imgid], thts[imgid]) )
