#!/usr/bin/env python

'''
====================
Author: Jesus Galaz - May/2017, Last update: may/2016
====================

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
'''
from __future__ import print_function
from __future__ import division

from past.utils import old_div
from builtins import range
from EMAN2 import *
from EMAN2jsondb import JSTask,jsonclasses

import sys
import numpy
import math
import random


def main():

	usage = """e2spt_transformplot.py file <options>.
	'file' can be a .json file or a stack of 3D subvolumes. In the latter case, the transformation (angles and translations) will be obtained from the header parameter 'xform.align3d'.
	Supplying the file directly is redundant with supplying the parameter --input.
	"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)	
	
	parser.add_argument("--input", type=str, default='', help="""A .json file with subtomogram alignment information, or a .hdf stack of aligned particles with correct alignment transformations in their header parameter 'xform.align3d'.""")
	parser.add_argument("--inversetoo",action="store_true",default=False,help="""Also plots the angles for the inverse of a transform.""")

	parser.add_argument("--path",type=str,default='spttransformplot',help="""Directory to store results in. The default is a numbered series of directories containing the prefix 'spttransformplot'; for example, spttransformplot_02 will be the directory by default if 'spttransformplot_01' already exists.""")
	parser.add_argument("--ppid", type=int, default=-1, help="Set the PID of the parent process, used for cross platform PPID")
	
	parser.add_argument("--subset",type=int,default=0,help="""Default=0 (not used). Plot only this substet of transforms from the hdf stack or json file provided.""")
	
	parser.add_argument("--verbose", "-v", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness", dest="verbose", action="store", metavar="n")
		
	(options, args) = parser.parse_args()	
	
	if not options.input:
		options.input = sys.argv[1]

	logger = E2init(sys.argv, options.ppid)

	orientations = {}

	if '.json' in options.input:				
		jsonfile= options.input
		jsonfileopen = js_open_dict( jsonfile )
		n=len(jsonfileopen)
		originaln=n
		print("\nthe number of transformations to plot is %d" %(n))
		if options.subset:
			n = options.subset
			print("\nplotting only a subset, n=%d" %(n))

		for j in range(n):
			xformslabel = 'subtomo_' + str( j ).zfill( len( str( originaln ) ) )
			t=jsonfileopen[xformslabel][0]
			print("\nread transform from .json file", t)
			orientations.update( {j:t} )
			#jsA.setval( xformslabel, [ t , score ] )
		jsonfileopen.close()

	elif '.hdf' in options.input:
		n=EMUtil.get_image_count(options.input)
		originaln=n
		print("\nthe number of transformations to plot is %d" %(n))
		if options.subset:
			n = options.subset
			print("\nplotting only a subset, n=%d" %(n))

		for j in range(n):
			t=EMData( options.input, j, True)['xform.align3d']
			orientations.update( {j:t} )


	azs=[]
	alts=[]
	phis=[]
	xs=[]
	ys=[]
	zs=[]

	azsi=[]
	altsi=[]
	phisi=[]
	xsi=[]
	ysi=[]
	zsi=[]

	if len( orientations ) > 2:

		from EMAN2_utils import makepath
		options = makepath(options,'spttransformplot')
		
		for i in orientations:
			t = orientations[i]
			print("\n t to get rotations and translations from transform number %d is" %(i))
			print(t)
			
			rots=t.get_rotation()
			
			az=rots['az']
			azs.append(az)

			alt=rots['alt']
			alts.append(alt)
			
			phi=rots['phi']
			phis.append(phi)
			
		
			trans=t.get_trans()
			
			x=trans[0]
			xs.append(x)

			y=trans[1]
			ys.append(y)

			z=trans[2]
			zs.append(z)

			if options.inversetoo:
				ti = t.inverse()
				print("\n t inverse to get rotations and translations from transform number %d is" %(i))
				print(ti)
				
				rotsi=ti.get_rotation()
				
				azi=rotsi['az']
				azsi.append(azi)

				alti=rotsi['alt']
				altsi.append(alti)
				
				phii=rots['phi']
				phisi.append(phii)
				
			
				transi=ti.get_trans()
				
				xi=transi[0]
				xsi.append(xi)

				yi=transi[1]
				ysi.append(yi)

				zi=transi[2]
				zsi.append(zi)


		textwriter(options, azs,'az')
		plotvals( options, azs, 'az' )

		textwriter(options, alts,'alt')
		plotvals( options, alts, 'alt' )

		textwriter(options, phis,'phi')
		plotvals( options, phis, 'phi' )

	
		textwriter(options, xs, 'x' )
		plotvals( options, xs, 'x', 1.0 )

		textwriter(options, ys, 'y' )
		plotvals( options,  ys,'y', 1.0 )

		textwriter(options, zs, 'z' )
		plotvals( options, zs, 'z', 1.0 )


		if options.inversetoo:
			textwriter(options, azsi,'az_inverse')
			plotvals( options, azsi, 'az_inverse' )

			textwriter(options, altsi,'alt_inverse')
			plotvals( options, altsi, 'alt_inverse' )

			textwriter(options, phisi,'phi_inverse')
			plotvals( options, phisi, 'phi_inverse' )

		
			textwriter(options, xsi, 'x_inverse' )
			plotvals( options, xsi, 'x_inverse', 1.0 )

			textwriter(options, ysi, 'y_inverse' )
			plotvals( options,  ysi, 'y_inverse', 1.0 )

			textwriter(options, zsi, 'z_inverse' )
			plotvals( options, zsi, 'z_inverse', 1.0 )

	else:
		print("\nthere's fewer than 2 transforms. no point in plotting a single (or null) value.")
		sys.exit(1)

	E2end(logger)

	return


def plotvals( options, vals, tag, binwidth=0 ):
	import matplotlib.pyplot as plt

	sigmavals= numpy.std(vals)
	meanvals = numpy.mean(vals)

	cuberoot = numpy.power(len(vals),old_div(1.0,3.0))
	width = old_div((3.5*sigmavals),cuberoot)

	if binwidth:
		width=binwidth
	
	#print "Therefore, according to Scott's normal reference rule, width = (3.5*std)/cuberoot(n), the width of the histogram bins will be", width
	
	minvals = min(vals)
	maxvals = max(vals)
	
	#if options.bins:
	#	calcbins = options.bins

	if 'az' in tag or 'phi' in tag:
		minvals = 0
		maxvals = 360
	elif 'alt' in tag:
		minvals = 0
		maxvals = 180
	elif 'x' in tag or 'y' in tag or 'z' in tag:
		pass

	calcbins = int(round( old_div((maxvals - minvals ), width) ))

	#count, bins, ignored = plt.hist(vals, 30, normed=True)
	ignored = plt.hist(vals, calcbins)
		
	#plt.plot(bins, 1/(sigmavals * numpy.sqrt(2 * numpy.pi)) * numpy.exp( - (calcbins - meanvals)**2 / (2 * sigmavals**2) ), linewidth=2, color='r')

	plt.title( tag + ' distribution' )
	plt.ylabel("n")
	plt.xlabel(tag)
	
	plt.savefig(options.path + '/' + tag + '.png')
	plt.clf()

	return


def textwriter(options,data,tag):
	
	#if options.path not in name:
	name = options.path + '/' + tag + '.txt'
	
	print("I am in the text writer for this file", name)
	
	f=open(name,'w')
	lines=[]
	for i in range(len(data)):
			
		line2write = str(i) + '\t' + str(data[i]) + '\n'
		#print "THe line to write is"
		lines.append(line2write)
	
	f.writelines(lines)
	f.close()

	return

if __name__ == '__main__':
	main()

