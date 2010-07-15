#!/usr/local/EMAN2/python/Python-2.6.5-ucs4/bin/python

# it is likely that the above line would have to be changed depending on python location

from EMAN2 import *
from sparx import *

from sys import  *
from string import *
from math import sqrt

doc_home = "/home/ryan/"
prj_home = "/home/ryan/Data_In/"

i_ctfs = open( doc_home + "ctfs.ext" )
ignore = i_ctfs.readline()

pixel = 2.38

Cs = 2.0
voltage = 200.0
amp_contrast = 10
bfactor = 0.0

total_high_proj = 0

proj_out = "bdb:/home/ryan/HOSPITALIS320"


for ii in xrange(1) :

	proj_in = prj_home + "ssspdbox320dcs2.dat"
	#proj_in = replace( proj_in, ' ', '0' )

	#fn_angl = doc_home + "FINAL/defgrp%3d/angles298.ext" % ii
	#fn_angl = replace( fn_angl, ' ', '0' )

	#fn_tran = doc_home + "FINAL/defgrp%3d/trans298-SPI13.ext" % ii
	#fn_tran = replace( fn_tran, ' ', '0' )

	#fn_high = doc_home + "FINAL/defgrp%3d/select_298_pop01.ext" % ii
	#fn_high = replace( fn_high, ' ', '0' )

	#selected_particles = readSpiderDoc(fn_high)

	nimage  = EMUtil.get_image_count( proj_in )
	#trans = read_txt_col(fn_tran)
	#angl = read_txt_col(fn_angl)

	print "Converting ", proj_in, " Number of my great particles is ",nimage

	prev_defocus = 0.

	for iq in xrange(nimage):
		defocus = atof( split( i_ctfs.readline() )[2] )
		if defocus != prev_defocus:
			print "Defocus is ", defocus
		prev_defocus = defocus
		i = iq  #int(selected_particles[iq][0])-1

		# Here one can have a problem because sometimes in trans doc two numbers
		# are connected together like "8.1753E-02-3.1129", readSpiderDoc can not handle this.
		#psi	= angl[i][2]
		#theta   = angl[i][3]
		#phi	= angl[i][4]

		#alpha   = trans[i][2]
		#sx	= trans[i][3]
		#sy	= trans[i][4]
		#mir	= trans[i][5]
		active = 1

		#assb, sxnb, synb, sc = compose_transform2(0.0, sx, sy, 1,  -alpha, 0.,0.,1)
		#if( mir > 0.5 ) :
		#	phi   = (540.0 + phi)%360.0
		#	theta = 180.0 - theta
		#	psi   = (540.0 - psi + assb)%360.0
		#else :
		#	psi = (psi + assb)%360.0

		data = EMData()
		data.read_image( proj_in, i )
		#  One could decimate the data here as spider images are usually oversampled.  This would require setting scale and adjusting
		#  pixel size appropriately
		#deci = Util.window(resample(data,scale),128,128,1,0,0,0)
		deci = data
		#set_params_proj(deci, [phi, theta, psi, sxnb, synb])
		set_params_proj(deci, [0.,0.,0.,0.,0.])
		deci.set_attr_dict({'active':1, 'ctf_applied':0})

		# Here, we convert the amp_contrast into the new convention
		set_ctf(deci, [defocus, Cs, voltage, pixel, bfactor, amp_contrast])

		deci.write_image( proj_out, total_high_proj )

		total_high_proj += 1



print "Total of ", total_high_proj, "  projections written"

