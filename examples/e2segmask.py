#!/usr/bin/env python

'''
====================
Author: Jesus Galaz - Sep/2017, Last update: Sep/2017
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

from EMAN2 import *
from EMAN2jsondb import JSTask,jsonclasses

import os


def main():

	usage = """e2segmask.py <options>.
	Program to produce volumes that are masked with a segmentation file.
	"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--data", type=str, default='', help="""Tomgoram file.""")

	parser.add_argument("--invert",action="store_true",default=False,help="""Reverse the contrast of the output file (needed if --data contains a raw tomogram with protein corresponding to 'black' densities).""")

	parser.add_argument("--mask", type=str, default='', help="""Segmentation mask. If you want to supply multiple masks, separate them by commas.""")
	
	parser.add_argument("--overwrite_size",action="store_true",default=False,help="""If the size of --mask does not match that of --data, --mask will be clipped; the --overwrite_size option means that the original mask will be deleted and replaced by the correctly sized one.""")

	#parser.add_argument("--path",type=str,default='segmask',help="""Directory to store results in. The default is a numbered series of directories containing the prefix 'segmask'; for example, segmas_02 will be the directory by default if 'segmask_01' already exists.""")
	
	parser.add_argument("--ppid", type=int, default=-1, help="Set the PID of the parent process, used for cross platform PPID")

	parser.add_argument("--tag", type=str, default='', help="""String to identify the masked feature, by appending it to the output file name. For example 'mt' for microtubules, or 'ribo' for ribosomes. Works only when 1 file is provided through --mask. For multiple masks, the output files will be numbered.""")

	parser.add_argument("--threshold",type=float,default=0.0,help="""Default=0.0. All negative densities will be zeroed out from --mask.""")

	parser.add_argument("--verbose", "-v", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness", dest="verbose", action="store", metavar="n")
	

	(options, args) = parser.parse_args()	
	
	if not options.data:
		options.input = sys.argv[1]

	logger = E2init(sys.argv, options.ppid)

	indata = options.data
	validextensions = ['.hdf','.HDF','.mrc','.mrcs','.MRC','.MRCS','.st','.ST','.ali','.ALI','.rec','.REC']
	extension=''
	for v in validextensions:
		if v in indata:
			extension=v

	masks = options.mask.split(',')
	
	outdata=''
	if extension:
		if options.tag and len( masks) < 2:
			outdata = indata.replace(extension, '_' + options.tag + '.hdf')
		else:
			outdata = indata.replace(extension, '.hdf')

	else:
		print ("\nERROR: invalid extension. Format must be " + ' or '.join([v for v in validextensions]) )
		print ("the extension of your data seems to be {} or {}".format(indata[-5:],indata[-4:]) )
		sys.exit(1)


	datahdr = EMData( indata, 0, True )
	
	datax = datahdr['nx']
	datay = datahdr['ny']
	dataz = datahdr['nz']

	count = 0
	for m in masks:
		maskhdr = EMData( m, 0, True )
		mx = maskhdr['nx']
		my = maskhdr['ny']
		mz = maskhdr['nz']

		cmd=''
		mask2use = m
		if datax != mx or datay != my or dataz != mz:
			if options.overwrite_size:
				cmd = 'e2proc3d.py ' + m + ' ' + m + ' --clip ' + str(datax) + ',' + str(datay) + ',' + str(dataz) #+ ' --process threshold.binary:value' + str(options.threshold)
			else:
				newmask = m.replace( m[-4:], '.hdf')
	
				cmd = 'e2proc3d.py ' + m + ' ' + newmask + ' --clip ' + str(datax) + ',' + str(datay) + ',' + str(dataz) #+ ' --process threshold.binary:value' + str(options.threshold)

				masks[count] = newmask

				mask2use = newmask

			print ("\nWARNING: --mask has a different size than --data; running the following command to clip --mask to the size of --data {}".format(cmd) )
			runcmd( options, cmd)
		
		newmaskhdr = EMData( mask2use, 0, True )
		newmx = maskhdr['nx']
		newmy = maskhdr['ny']
		newmz = maskhdr['nz']

		if datax != newmx or datay != newmy or dataz != newmz:
			print ("\nERROR: cannot fix size discrepancy. clip --mask manually. size of --data is x={} y={} z={}; size of --mask is x={} y={} z={}".format(datax,datay,dataz,newmx,newmy,newmz))
			sys.exit(1)

		count+=1

		thresholdvalstring = str(options.threshold).replace('.','p')

		mask2useth = mask2use.replace( mask2use[ -4:], '_th' + thresholdvalstring + '.hdf' )

		cmd = 'e2proc3d.py ' + mask2use + ' ' + mask2useth +  ' --process threshold.binary:value' + str(options.threshold)

		print ("\napplying threshold to mask")
		runcmd( options, cmd )

		if len(masks) > 1:
			num = str(count).zfill( len(str(len(masks))))
			outdata = outdata.replace('.hdf','_' + num + '.hdf')

		cmd  = 'e2proc3d.py ' + indata +  ' ' + outdata + ' --multfile ' + mask2useth

		print("\nmultiplying data by mask {}".format(mask2useth))
		runcmd( options, cmd )


		if options.invert:
			cmd = 'e2proc3d.py ' + outdata + ' ' + outdata.replace('.hdf','_inv.hdf') + ' --mult -1 ' + ' --process threshold.belowtozero:minval=' + str(options.threshold)
			
			print ("\ninverting contrast")
			runcmd( options, cmd )


	E2end(logger)
	
	return


def runcmd(options,cmd):
	if options.verbose > 9:
		print ("\n(e2segmask)(runcmd) running command = {}".format(cmd))
	
	p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	text=p.communicate()	
	p.stdout.close()
	
	if options.verbose > 8:
		print ("\n(e2segmask)(runcmd) done")
	
	return 1


if __name__ == '__main__':
	
	main()
