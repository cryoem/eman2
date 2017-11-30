#!/usr/bin/env python

'''
====================
Author: Jesus Galaz - nov/2017, Last update: nov/2017
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

from optparse import OptionParser
from EMAN2 import *
from EMAN2jsondb import JSTask,jsonclasses

import sys
import numpy
import math
import random


def main():

	usage = """e2spt_translate.py file <options>.
	This program translationally aligns a stack of images to a reference.
	"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--input", type=str, default=None, help="""default=None. stack of image in .hdf format, presumably raw, without any alignment.""")

	parser.add_argument("--ref", type=str, default=None, help="""default=None. volume to be used as a reference.""")
	
	#parser.add_argument("--transformfile", type=str, default=None, help="""default=None. .json file with alignment parameters produced by other e2spt programs""") 
	
	#parser.add_argument("--translateonly", action='store_true', default=False, help="""default=False (not used). requieres --input. If supplied, this option will ensure that only translations are applied to the --input stack.""")

	parser.add_argument("--path",type=str,default='spttranslate',help="""Directory to store results in. The default is a numbered series of directories containing the prefix 'spttranslate'; for example, spttranslate_02 will be the directory by default if 'spttranslate_01' already exists.""")
	
	parser.add_argument("--ppid", type=int, default=-1, help="set the PID of the parent process, used for cross platform PPID.")
	
	#parser.add_argument("--subset",type=int,default=0,help="""Default=0 (not used). Plot only this substet of transforms from the hdf stack or json file provided.""")
	parser.add_argument("intonly", action='store_true', default=False, help="default=Flase. If on, this will allow integer translations only.")

	parser.add_argument("masked", action='store_true', default=False, help="default=False. treat zero pixels in --ref as a mask for normalization.")
	
	parser.add_argument("maxshift", type=int, default=0, help="default=0 (not used). Maximum allowable translation in pixels.")
	
	parser.add_argument("nozero", action='store_true', default=False, help="default=False. Zero translations not permitted (useful for CCD images)")
	
	parser.add_argument("useflcf", action='store_true', default=False, help="default=False. Use Fast Local Correlation Function rather than CCF.")

	parser.add_argument("--verbose", "-v", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness", dest="verbose", action="store", metavar="n")
	
	#arser.add_argument("--inversetoo",action="store_true",default=False,help="""Also plots the angles for the inverse of a transform.""")
	
	(options, args) = parser.parse_args()	

	logger = E2init(sys.argv, options.ppid)

	alignerparams = {}
	
	if options.intonly:
		alignerparams.update({"intonly":1})
	
	if options.masked:
		alignerparams.update({"masked":1})
	
	if options.maxshift:
		alignerparams.update({"maxshift":options.maxshift})
	
	if options.nozero:
		alignerparams.update({"nozero":1})
	
	if options.useflcf:
		alignerparams.update({"useflcf":1})
	
	from EMAN2_utils import makepath
	options = makepath(options,'spttranslate')

	n = EMUtil.get_image_count(options.input)

	ref = EMData(options.ref,0)
	
	jsAliParamsPath = options.path + '/sptali_trans.json'

	jsA = js_open_dict( jsAliParamsPath )
	
	basename = os.path.basename(options.input)
	stem,extension = os.splitext(basename)
	outfile = options.path + '/' + stem + '_translated.hdf'

	for i in xrange(0,n):

		img = EMData(options.input,i)
		
		ali = img.align('translational',ref,alignerparams)
		
		print('\nali={}'.format(ali))

		#bestcoarse = img.xform_align_nbest(options.align[0],sfixedimage,options.align[1],options.npeakstorefine,options.aligncmp[0],options.aligncmp[1])

		xformslabel = 'subtomo_' + str( 0 ).zfill( len( str( nptcl ) ) )
		t = ali['xform.align3d']
		score = ali['score']			
		jsA.setval( xformslabel, [ t , score ] )

		imgt=img.transform(t)

		imgt.write_image(outfile,i)


	jsA.close()
	
	E2end(logger)

	return


if __name__ == '__main__':
	
	main()
