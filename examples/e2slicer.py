#!/usr/bin/env python
from past.utils import old_div
from builtins import range
'''
====================
Author: Jesus Galaz-Montoya - 02/March/2013, Last update: 01/2023
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
from EMAN2_utils import *

import sys
from sys import argv
import operator
import random
import numpy

from sys import argv

import time
from datetime import timedelta

def main():
	start = time.perf_counter()

	usage = """
			This program produces orthogonal slices of an EM volume of --nslices thickness.
			e2slicer.py <input_file1> <input_file2> ... <input_fileN> <options> . 
			The options should be supplied in "--option=value" format (or --options=value:parameter1=value:parameter2=value... etc, 
			replacing "option" for a valid option name, and "value", "parameter" for a acceptable entries for that option. 
			This program extracts slices from a volume or multiple volumes in a stack. By default, the program will extract three orthogonal slices 
			(one from each direction X, Y and Z) that go through the middle of the volume(s). Other options also provide all the slices along any of the 3 cartesian axes. 
			"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)	
	
	parser.add_argument("--allx",action='store_true',default=False,help="Get ALL the slices in a volume along the x axis.")
	parser.add_argument("--ally",action='store_true',default=False,help="Get ALL the slices in a volume along the y axis.")
	parser.add_argument("--allz",action='store_true',default=False,help="Get ALL the slices in a volume along the z axis.")
	
	parser.add_argument("--input", type=str, default='', help="""Default=None. This is redundant with supplying input files directly. The file name containing a volume or stack of volumes from which you want to generate slices. You can supply more than one file either by providing an .hdf stack of volumes, or by providing multiple files separated by commas.""")
		
	parser.add_argument("--nslices", type=int, default=1,help="""default=1. Number of slices to average around the central sections (not compatible with --allx, --ally, or --allz)""")

	parser.add_argument("--onlymidx",action='store_true',default=False,help="Only extract the middle slice of the volume parallel to the YZ plane.")
	parser.add_argument("--onlymidy",action='store_true',default=False,help="Only extract the middle slice of the volume parallel to the XZ plane.")
	parser.add_argument("--onlymidz",action='store_true',default=False,help="Only extract the middle slice of the volume parallel to the XY plane.")
	
	parser.add_argument("--path",type=str,default="slices",help="""Defautl=slices. Directory to store results in. The default is a numbered series of directories containing the prefix 'slices'; for example, slices_02 will be the directory by default if 'slices_01' already exists.""")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	parser.add_argument("--singlestack",action='store_true',default=False,help="""This option will save slices from all particles into a single .hdf stack file if --onlymidz or --onlymidy or --onlymidx are provided, instead of one slice file per volume.""")
	parser.add_argument("--shrink",type=int,default=0,help="""Integer factor to use for shrinking volumes prior to extracting slices.""")
	
	parser.add_argument("--orthogonaloff",action='store_true',default=False,help="""By default, the program will extract three orthogonal slices through the middle of the input volume(s). If this parameter is specified, it will not, and only the other options that are supplied will be valid.""")
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="default=0. verbose level [0-9], higher number means higher level of verboseness.")


	(options, args) = parser.parse_args()	
	
	#c:this checks whether input is supplied directly via sys.argv, or --input, counts valid input files, and puts them in a list for future use.
	options,inputs = checkinput( options )

	if not inputs:
		print("\nERROR: failed to load any valid input files. Make sure you're in the correct directory and filenames are correct. Exiting.")
		sys.exit(1)

	if options.verbose: print("\ninputs are {}".format(inputs))

	logger = E2init(sys.argv, options.ppid)

	
	'''#
	#Check for sanity of some supplied parameters
	'''#
	if options.onlymidz:
		if options.onlymidx or options.onlymidy:
			print("ERROR: You can only supply one of --onlymidx, --onlymidy or --onlymidz at a time.")
			sys.exit(1)
	
	if options.onlymidx:
		if options.onlymidy or options.onlymidz:
			print("ERROR: You can only supply one of --onlymidx, --onlymidy or --onlymidz at a time.")
			sys.exit(1)
			
	if options.onlymidy:
		if options.onlymidx or options.onlymidz:
			print("ERROR: You can only supply one of --onlymidx, --onlymidy or --onlymidz at a time.")
			sys.exit(1)

	if options.allx or options.ally or options.allz:
		if options.onlymidx or options.onlymidy or options.onlymidz:
			print("ERROR: Cannot supply --allx, --ally or --allz at the same time than any of --onlymidx={}, --onlymidy={} or --onlymidz={}".format(options.onlymidx,options.onlymidy,otions.onlymidz))
			sys.exit(1)

	'''#
	#Make a directory where to store the results
	'''#
	from EMAN2_utils import makepath
	options = makepath(options,'slices')
	
	
	'''#
	#Generate orthogonal slice regions, using the transforms to get x and y axis perspectives
	'''#

	t = Transform()
	ty = Transform({'type':'eman','az':0,'alt':-90,'phi':0})
	tx = Transform({'type':'eman','az':-90,'alt':90,'phi':90})
	
	axes_dict = {}
	if options.onlymidx:
		axes_dict.update({'x':tx})
	elif options.onlymidy:
		axes_dict.update({'y':ty})
	elif options.onlymidz:
		axes_dict.update({'z':t})													

	elif options.allz or options.allx or options.ally:
		if options.allx: 
			axes_dict.update({'x':tx})
		if options.ally: 
			axes_dict.update({'y':ty})
		if options.allz: 
			axes_dict.update({'z':t})

	if not options.orthogonaloff:
		axes_dict.update({'x':tx,'y':ty,'z':t})


	for f in inputs:
		print("\nprocessing file {}".format(f))
		ext = os.path.splitext(f)[-1]
		try:
			hdr=EMData(f,0,True)
		except:
			print("\ninvalid file f={}".format(f))

		n = EMUtil.get_image_count( f )
		
		for i in range(n):
			if options.verbose: print("\nprocessing particle {}".format(i))
			a = EMData(f,i)
			
			ap = a.copy()
			if options.shrink:
				nyquist = ap['apix_x'] * 2.0
				shnyquist = nyquist * options.shrink
				ap.process_inplace('filter.lowpass.tanh',{'cutoff_freq':old_div(1.0,shnyquist)})
				ap.process_inplace('math.fft.resample',{'n':options.shrink})

			ptcltag = ''
			slicestag = ''
			
			if n > 1:
				slicestag = '_SLICESmid'
				
				if options.allz or options.ally or options.allx or not options.singlestack:
					ptcltag = '_ptcl' + str(i).zfill( len( str(n) ))
			
	
			for axis in axes_dict:
				slicef=None
				if options.nslices > 1 or options.allx or options.ally or options.allz:
					t=axes_dict[axis]

					aprot = ap.copy()
					if options.verbose>5: 
						print("\n(e2slicer)(main) --nslices={} > 1, therefore rotating ptcl with transform={} because axis={}".format(options.nslices,t,axis))
					aprot.transform(t)

					if options.nslices >1:
						slicef=get_slices(options,aprot,axis=axis)

					if options.allx and axis == 'x':
						if options.verbose > 5: print("Generating all x slices")		
						rotvolxname = options.path + '/tmp_rotx.hdf'
						aprot.write_image(rotvolxname,0)
						outname = options.path + '/' + os.path.basename( f ).replace(ext,ptcltag+'_all_x.hdf')
						runcmd(options,'e2proc2d.py ' + rotvolxname + ' ' + outname + ' --threed2twod && cd ' + options.path + ' && rm tmp*')
					if options.ally and axis == 'y':
						if options.verbose > 5: print("Generating all y slices")
						rotvolyname = options.path + '/tmp_roty.hdf'
						aprot.write_image(rotvolyname,0)
						outname = options.path + '/' + os.path.basename( f ).replace(ext,ptcltag+'_all_y.hdf')
						runcmd(options,'e2proc2d.py ' + rotvolyname + ' ' + outname + ' --threed2twod && cd ' + options.path + ' && rm tmp*')
					if options.allz and axis == 'z':
						if options.verbose > 5: print("Generating all z slices")
						rotvolzname = options.path + '/tmp_rotz.hdf'
						aprot.write_image(rotvolzname,0)
						outname = options.path + '/' + os.path.basename( f ).replace(ext,ptcltag+'_all_z.hdf')
						runcmd(options,'e2proc2d.py ' + f + ' ' + outname + ' --threed2twod && cd ' + options.path + ' && rm tmp*')

				else:
					if options.verbose>5: print("\n(e2slicer)(main) getting single slice from axis={}".format(axis))
					slicef=get_single_slice(options,ap,axis=axis)

				if slicef:
					if options.onlymidx or options.onlymidy or options.onlymidz:
						slices_mid_output = options.path + '/' + os.path.basename( f ).replace(ext,ptcltag + '_' + axis + '.hdf')
						slicef.write_image(slices_mid_output,i)

					elif not options.orthogonaloff:
						slices_ortho_output = options.path + '/' + os.path.basename( f ).replace(ext,ptcltag + '_ortho.hdf')
						#print("\noutput name for ortho stack is {}".format(slices_ortho_output))
						#print("\noptions.nslices={}".format(options.nslices))
						if options.nslices > 1:
							slices_ortho_output = slices_ortho_output.replace(".hdf","_nthick" + str(options.nslices) +".hdf" )
							#print("\nupdated output name for ortho stack is {} becasue --nslices={}".format(slices_ortho_output, options.nslices))
						slicef.write_image(slices_ortho_output,-1)
	
	E2end(logger)
	
	elapsed = time.perf_counter() - start
	print(str(timedelta(seconds=elapsed)))

	return


def get_single_slice(options,img,axis=''):

	nx=img['nx']
	ny=img['ny']
	nz=img['nz']

	rmid=None

	if axis=='x':
		rmid = Region(old_div(nx,2), 0, 0, 1, ny, nz)
		x=nz
		y=ny
	elif axis=='y': 
		rmid = Region(0, old_div(ny,2), 0, nx, 1, nz)
		x=nx
		y=nz
	elif axis=='z': 
		rmid = Region(0, 0, old_div(nz,2), nx, ny, 1)
		x=nx
		y=ny

	if options.verbose > 6: print("The region for the orthogonal z slice is", rmid)

	print("\n(e2slicer)(get_single_slice) getting slice from region rmid={} because axis={}".format(rmid,axis))
	slicemid = img.get_clip( rmid )

	#c: apparently, setting the size of the output is essential for it to show as a plane in xy as opposed to viewing the slice "from the side"
	slicemid.set_size(x,y,1)

	return slicemid


def get_slices(options,img,axis=''):

	nx=img['nx']
	ny=img['ny']
	nz=img['nz']
	
	rmid=None

	rmid = Region(0, 0, old_div(nz,2)-int(ceil(options.nslices/2.0)), nx, ny, options.nslices)
	if options.verbose > 6: print("The region for the orthogonal z nslices is", rmid)

	slicemid = img.get_clip( rmid )
	slicemid.set_size(nx,ny,options.nslices)

	prj = slicemid.project("standard",Transform())
	prj.set_attr('xform.projection',Transform())
	prj.set_size(nx,ny,1)

	return prj


if __name__ == '__main__':
	main()
