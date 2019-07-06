#!/usr/bin/env python

'''
====================
Author: Jesus Galaz - oct/2017, Last update: nov/2017
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

from builtins import range
from EMAN2 import *
from EMAN2jsondb import JSTask,jsonclasses

import sys
import numpy
import math
import random


def main():

	usage = """e2spt_transform.py file <options>.
	This program applies a transform (rotations and translations) to an image stack, or (translations only) to a coordinates file. The transformations will be read from a .json file with transforms in it for each volume in the image stack or each row in the coordinates file.
	"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--input", type=str, default=None, help="""default=None. stack of image in .hdf format, presumably raw, without any alignment.""")

	parser.add_argument("--coords", type=str, default=None, help="""default=None. coordinates file to apply translations from --transformfile.""")
	
	parser.add_argument("--transformfile", type=str, default=None, help="""default=None. .json file with alignment parameters produced by other e2spt programs""") 
	
	parser.add_argument("--translateonly", action='store_true', default=False, help="""default=False (not used). requieres --input. If supplied, this option will ensure that only translations are applied to the --input stack.""")

	#parser.add_argument("--path",type=str,default='spttransform',help="""Directory to store results in. The default is a numbered series of directories containing the prefix 'spttransform'; for example, spttransform_02 will be the directory by default if 'spttransform_01' already exists.""")
	parser.add_argument("--ppid", type=int, default=-1, help="set the PID of the parent process, used for cross platform PPID")
	
	#parser.add_argument("--subset",type=int,default=0,help="""Default=0 (not used). Plot only this substet of transforms from the hdf stack or json file provided.""")
	
	parser.add_argument("--verbose", "-v", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness", dest="verbose", action="store", metavar="n")
	
	#arser.add_argument("--inversetoo",action="store_true",default=False,help="""Also plots the angles for the inverse of a transform.""")
	
	(options, args) = parser.parse_args()	

	if options.input and options.coords:
		print("\nto minimize confusion, apply --transformfile to --input and to --coords on separate runs of the program.")
		sys.exit(1)
	#elif if not coords and not options.input:
	#	options.input = sys.argv[1]

	logger = E2init(sys.argv, options.ppid)

	orientations = {}

	#nptcls = 0
	#ncoords=0
	n = 0
	lines = []
	newlines = []
	coordlines={}
	if options.coords:
		with open( options.coords, 'r') as coordsfile:
			lines = coordsfile.readlines()
			n = int(len(lines))
			k=0
			for line in lines:
				#print "\nadding line {}".format(line)
				coordlines.update({k:line})
				k+=1

	elif options.input:
		n = int(EMUtil.get_image_count(options.input))

	#n=nptcls
	if options.transformfile:				
		jsonfile = options.transformfile
		jsonfileopen = js_open_dict( jsonfile )
		#print "\njsonfile to open = {}".format(jsonfile)
		nt = int(len(jsonfileopen))
		
		#originaln=n
		#print "\nthe number of transformations to plot is %d" %(n)
		#print "\n n={}, type={}, nt={}, type={}".format(n,type(n),nt,type(nt))
		if int(n) != int(nt):
			print("int(n)={}, int(nt)={}".format(int(n),int(nt)))
			if int(n) < int(nt):
				print("\nWARNING: fewer particles/coordinates n={} than transform parameters nt={}. Only transforming as many particles from {} as transforms in {}".format(n,nt,options.input,options.transformfile))
			elif int(n) > int(nt):
				n = int(nt)
				print("\nWARNING: more particles/coordinates n={} than transform parameters nt={}".format(n,nt))

		keys = list(jsonfileopen.keys())
		keys.sort()
		for j in range(n):
			label = keys[j]
			t = jsonfileopen[label][0]
			
			if options.input:
				ptcl = EMData(options.input,j)
				outname = options.input.replace('.hdf','_transformed.hdf')
				#translations are in the frame of the rotated particle.
				#to apply translations only, you need to convert them to the unrotated frame by compute the inverse transform Ti
				#since Ti will be 'untranslated', the negative of the translations from Ti gives the translation that 'centers' the particles translationally (assuming they were correctly aligned)
				if options.translateonly:
					ttransform = translationtransform(t)
					ptcl.transform(ttransform)
					outname = options.input.replace('.hdf','_translated.hdf')				
				else:
					ptcl.transform(t)

				ptcl.write_image(outname,j)
	
			elif options.coords:
				coord = coordlines[j].split()
				x = float(coord[0])
				y = float(coord[1])
				z = float(coord[2])
				tc = Transform({'type':'eman','tx':x,'ty':y,'tz':z})

				ttransform = translationtransform(t)
				
				tfinal = ttransform*tc
				
				newcoord = tfinal.get_trans()
				newx = int(round(float(newcoord[0])))
				newy = int(round(float(newcoord[1])))
				newz = int(round(float(newcoord[2])))
				newlines.append(str(newx)+'\t'+str(newy)+'\t'+str(newz)+'\n')
			
			if options.verbose > 9:
				print("\ntransformed particle/coordinate n={}/{}, using transform={} from .json file {}".format( j, n, t, options.transformfile ))

		if options.coords:
			name,ext = os.path.splitext(options.coords)
			if newlines:
				if int(len(newlines)) < int(len(coordlines)):
					dif = len(coordlines) - len(newlines)
					print("\nlen(newlines)={} , len(coordlines)={} , dif={}".format(len(newlines),len(coordlines),dif))
					for kk in range(dif):
						indx = len(coordlines) -1*kk -1
						print("\nWARNING: appending untranslated cooridates from --coords {}".format(coordlines[indx]))
						newlines.append(coordlines[indx])	#if the coordinates file grew with respect to the alignment parameters file because new particles were boxed (but not aligned), this adds untranslated coordinates from the end of the original --coords file 

				with open( name + '_trans' + ext, 'w') as newcoordsfile:
					newcoordsfile.writelines(newlines)
			elif not newlines:
				print("\nERROR: no transformed coordinates to write. something went terribly -terribly- wrong! (sorry).")
				sys.exit(1)


			#orientations.update( {j:t} )
			#jsA.setval( xformslabel, [ t , score ] )
		jsonfileopen.close()

	E2end(logger)

	return


def translationtransform(t):
	
	ti = t.inverse()
	solutiontrans = -1* ti.get_trans()
	solutiont = Transform({'type':'eman','tx':solutiontrans[0],'ty':solutiontrans[1],'tz':solutiontrans[2]})

	return solutiont


if __name__ == '__main__':
	main()

