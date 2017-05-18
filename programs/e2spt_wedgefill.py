#!/usr/bin/env python

'''
====================
Author: Steven Ludkte - 2016, Last update: May/2017 (Jesus)
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

from EMAN2 import *
from sys import argv,stdout,exit


def main():

	#if len(argv)<4 : 
	#	print """usage: e2spt_wedgefill <input stack> --options
	#This program will identify and fill in the missing wedge in a stack of subtomograms (with a missing wedge) which have been aligned to an average (with no missing wedge). While in a sense this creates 'Frankenstien particles' it is assumed that having values from an average of many particles is better than having zero in the same areas, at least for certain purposes.
	#
	#You must provide an average volume from which to draw the missing values and it is critical that this volume be an average of the particles being corrected. The particles should also have been properly normalized prior to averaging to achieve the desired effect.
	#
	#Specifically, this process is designed to make it possible to run PCA or other classification processes on the particles. Any missing values will all be at the same point near the center of the cloud defined by the entire population, rather than having zero values distorting the cloud. Once particles have been classified, the originals wihout filled wedge can be used for averaging. 
	#
	#	exit(1)
	
	usage = """usage: e2spt_wedgefill <input stack> --options
	This program will identify and fill in the missing wedge in a stack of subtomograms (with a missing wedge) which have been aligned to an average (with no missing wedge). While in a sense this creates 'Frankenstien particles' it is assumed that having values from an average of many particles is better than having zero in the same areas, at least for certain purposes.
	
	You must provide an average volume from which to draw the missing values and it is critical that this volume be an average of the particles being corrected. The particles should also have been properly normalized prior to averaging to achieve the desired effect.
	
	#Specifically, this process is designed to make it possible to run PCA or other classification processes on the particles. Any missing values will all be at the same point near the center of the cloud defined by the entire population, rather than having zero values distorting the cloud. Once particles have been classified, the originals wihout filled wedge can be used for averaging. 
	#.
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)	
	parser.add_argument("--input", type=str, default='', help="""this is redundant with supplying the image stack directly after the program name.""")
	parser.add_argument("--output", type=str, default='', help="""optional file name to save the stack with a filled wedge.""")

	parser.add_argument("--fillimage", type=str, default='', help="""the iamge to use to fill in the missing wedge of the images in --input. Ideally, this is the average of the aligned images in --input, or the reference the images in --input were aligned to.""")
	
	parser.add_argument("--fillwithnoise",action="store_true",default=False,help="""this will fill the missing wedge with gaussian noise. --matchto will be turned on by default if this option is supplied.""")

	parser.add_argument("--matchto",action="store_true",default=False,help="""this will match the power spectrum of each image in --input to that --fillimage so that things are properly normalized.""")

	#parser.add_argument("--path",type=str,default='spttransformplot',help="""Directory to store results in. The default is a numbered series of directories containing the prefix 'spttransformplot'; for example, spttransformplot_02 will be the directory by default if 'spttransformplot_01' already exists.""")
	parser.add_argument("--ppid", type=int, default=-1, help="Set the PID of the parent process, used for cross platform PPID")
	parser.add_argument("--subset",type=int,default=0,help="""Default=0 (not used). Plot only this substet of transforms from the hdf stack or json file provided.""")
	parser.add_argument("--verbose", "-v", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness", dest="verbose", action="store", metavar="n")
	
	(options, args) = parser.parse_args()	

	logger = E2init(sys.argv, options.ppid)

	if not options.input:
		options.input = argv[1]

	if '.hdf' not in options.input[-4:]:
		print "\nerror. --input must be in .hdf format"
		sys.exit(1)

	if not options.fillimage and not options.fillwithnoise:
		print "\nERROR: you must provide either --fillimage or --fillwithnoise."
		sys.exit(1)
	elif options.fillimage and options.fillwithnoise:
		print "\nERROR: you must provide only one of either --fillimage or --fillwithnoise."
		sys.exit(1)
	else:

		if not options.output:
			options.output = options.input.replace('.hdf','_filled.hdf')

		fill = None
		fillf = None

		if options.fillimage:
			fill=EMData(options.fillimage,0)
			
			if not options.matchto:
				fillf=fill.do_fft()

		elif options.fillwithnoise:
			inputhdr=EMData(options.input,0,True)
			fill=EMData(inputhdr['nx'],inputhdr['ny'],inputhdr['nz'])
			fill.to_one()
			radius = min(inputhdr['nx'],inputhdr['ny'],inputhdr['nz'])/2 -4
			fill.process_inplace('mask.soft',{'outer_radius':radius})

			options.matchto = True
		
		#n=EMUtil.get_image_count(argv[1])
		
		n=EMUtil.get_image_count(options.input)

		if options.subset:
			n = options.subset

		fillp = fill.copy()
		for i in xrange(n):
			im=EMData(options.input,i)
			
			if options.fillwithnoise:
				filln = fill.copy()
				filln.process_inplace('testimage.noise.gauss', {'mean':0,'sigma':1})
				fillp = filln.copy()

			if options.matchto:
				fillmatch = fillp.process('filter.matchto',{'to':im})
				fillp = fillmatch.copy()
			
			if options.fillwithnoise or options.matchto:
				fillf = fillp.do_fft()

			imf=im.do_fft()

			imf.process_inplace("mask.wedgefill",{"fillsource":fillf})
		#	imf.process_inplace("filter.lowpass.gauss",{"cutoff_freq":.05})
		#	imf.process_inplace("filter.highpass.gauss",{"cutoff_pixels":2})
			im=imf.do_ift()
			im.process_inplace("normalize")
			im.write_image( options.output, i )

			if i%10==0:
				print "  %d/%d\r"%(i+1,n),
				sys.stdout.flush()

		print "\ndone\n"


		E2end(logger)

	return


if __name__ == '__main__':
	
	main()
