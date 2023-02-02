#!/usr/bin/env python
from past.utils import old_div
from builtins import range
'''
====================
Author: Jesus Galaz-Montoya - 02/March/2013, Last update: 07/Nov/2017
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
import sys
from sys import argv
import EMAN2
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
	
	parser.add_argument("--lowpass",type=str, default='', help="Default=None. A lowpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.")

	parser.add_argument("--mask",type=str,default='',help="Default=None. Mask processor applied to particles before alignment." )
	
	parser.add_argument("--normproc",type=str,default='',help="""Default=None (not used). Normalization processor applied to particles before computing slices.""")
	parser.add_argument("--nslices", type=int, default=1,help="""default=1. Number of slices to average around the central sections (not compatible with --allx, --ally, or --allz)""")

	parser.add_argument("--onlymidx",action='store_true',default=False,help="Only extract the middle slice of the volume parallel to the YZ plane.")
	parser.add_argument("--onlymidy",action='store_true',default=False,help="Only extract the middle slice of the volume parallel to the XZ plane.")
	parser.add_argument("--onlymidz",action='store_true',default=False,help="Only extract the middle slice of the volume parallel to the XY plane.")
	
	parser.add_argument("--path",type=str,default="slices",help="""Defautl=slices. Directory to store results in. The default is a numbered series of directories containing the prefix 'slices'; for example, slices_02 will be the directory by default if 'slices_01' already exists.""")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	parser.add_argument("--singlestack",action='store_true',default=False,help="""This option will save slices from all particles into a single .hdf stack file if --onlymidz or --onlymidy or --onlymidx are provided, instead of one slice file per volume.""")
	parser.add_argument("--shrink",type=int,default=0,help="""Integer factor to use for shrinking volumes prior to extracting slices.""")
	
	parser.add_argument("--orthogonaloff",action='store_true',default=False,help="""By default, the program will extract three orthogonal slices through the middle of the input volume(s). If this parameter is specified, it will not, and only the other options that are supplied will be valid.""")
	
	#parser.add_argument("--saverotvol",action='store_true',default=False,help="Will save the volume in each rotated position used to generate a projection.")

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness.")

	#parser.add_argument("--transformsfile",type=str,help="A text files containing lines with one triplet of az,alt,phi values each, representing the transforms to use to project a single volume supplied. ", default='')
	#parser.add_argument("--angles",type=str,help="A single comma or space separated triplet of az,alt,phi values representing the particle rotation to apply before projecting it.", default='')
	#parser.add_argument("--tag",type=str,help="When supplying --angles, tag the output projection with a string provided through --tag", default='')

	(options, args) = parser.parse_args()	
	
	#c:this checks whether input is supplied directly via sys.argv, or --input, counts valid input files, and puts them in a list for future use.
	options,inputs = checkinput( options )

	print("\ninputs are {}".format(inputs))

	logger = E2init(sys.argv, options.ppid)

	if options.lowpass: 
		options.lowpass=parsemodopt(options.lowpass)

	if options.mask: 
		options.mask=parsemodopt(options.mask)
			
	if options.normproc: 
		options.normproc=parsemodopt(options.normproc)

	
	'''#
	#Check for sanity of some supplied parameters
	'''#
	if options.onlymidz:
		if options.onlymidx or options.onlymidy:
			print("ERROR: You can only supply one of --onlymidx, --onlymidy or --onlymidz at a time.")
			sys.exit()
	
	if options.onlymidx:
		if options.onlymidy or options.onlymidz:
			print("ERROR: You can only supply one of --onlymidx, --onlymidy or --onlymidz at a time.")
			sys.exit()
			
	if options.onlymidy:
		if options.onlymidx or options.onlymidz:
			print("ERROR: You can only supply one of --onlymidx, --onlymidy or --onlymidz at a time.")
			sys.exit()


	'''#
	#Make a directory where to store the results
	'''#
	from EMAN2_utils import makepath
	options = makepath(options,'slices')
	
	
	'''#
	#Generate orthogonal slice regions
	'''#
	for f in inputs:
		ext = os.path.splitext(f)[-1]

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
				
			nx=ap['nx']
			ny=ap['ny']
			nz=ap['nz']
			
			ptcltag = ''
			slicestag = ''
			
			if n > 1:
				slicestag = '_SLICESmid'
				
				if options.allz or options.ally or options.allx or not options.singlestack:
					ptcltag = '_ptcl' + str(i).zfill( len( str(n) ))
			
			rmid = None
						
			if options.onlymidz or options.onlymidy or options.onlymidx:
				t = Transform()
				
				if options.onlymidz:	
					#rmid = Region(0, 0, old_div(nz,2), nx, ny, 1)
					#if options.nslices:
					rmid = Region(0, 0, old_div(nz,2)-int(ceil(options.nslices/2.0)), nx, ny, options.nslices)
					if options.verbose > 6: print("The region for the orthogonal z slice is", rmid)
					
					slicestag += 'z'
					if n < 2:
						slicestag = '_SLICEmidz'			
										
				elif options.onlymidx:
					t = Transform({'type':'eman','az':90,'alt':90,'phi':0})

					#rmid = Region(old_div(nx,2), 0, 0, 1, ny, nz)
					#if options.nslices:
					rmid = Region(old_div(nx,2)-int(ceil(options.nslices/2.0)), 0, 0, options.nslices, ny, nz)
						
					if options.verbose: print("The region for the orthogonal x slice is", rmid)

					#slicemidx=a.get_clip(rmidx)
					#slicemidx.write_image(options.path + '/' + options.input.replace('.',ptcltag+'_SLICEmidx.'),0)
					
					slicestag += 'x'
					if n < 2:
						slicestag = '_SLICEmidx'
					
					#Tx = Transform({'type':'eman','az':0,'alt':-90,'phi':0})
					#ap.transform( Tx )
		
				elif options.onlymidy:
					t = Transform({'type':'eman','az':0,'alt':90,'phi':0})

					#rmid = Region(0, old_div(ny,2), 0, nx, 1, nz)
					#if options.nslices:
					rmid = Region(0, old_div(ny,2)-int(ceil(options.nslices/2.0)), 0, nx, options.nslices, nz)
					if options.verbose: print("The region for the orthogonal y slice is", rmid)
		
					#slicemidy=a.get_clip(rmidy)
					#slicemidy.write_image(options.path + '/' + options.input.replace('.',ptcltag+'_SLICEmidy.'),0)
					
					slicestag += 'y'
					if n < 2:
						slicestag = '_SLICEmidy'
					
					#Ty = Transform({'type':'eman','az':0,'alt':-90,'phi':-90})
					#ap.transform( Ty )
								
				slicemid = ap.get_clip( rmid )
				prj = slicemid.project("standard",t)
				prj.set_attr('xform.projection',t)

				slicemidf=prj.copy()
				slicemidf.set_size( nx, ny, 1)
				
				slicemidf.write_image(options.path + '/' + os.path.basename( f ).replace(ext,ptcltag + slicestag + '.hdf'),i)
			
			elif not options.onlymidz and not options.onlymidy and not options.onlymidx:
				app = a.copy()
				if options.shrink:
					nyquist = app['apix_x'] * 2.0
					shnyquist = nyquist * options.shrink
					app.process_inplace('filter.lowpass.gauss',{'cutoff_freq':old_div(1.0,shnyquist)})
					app.process_inplace('math.meanshrink',{'n':options.shrink})
				
				#regions={}
				if not options.orthogonaloff:
					if options.verbose > 5: print("Generating orthogonal slices")
					rmidz = Region(0, 0, old_div(nz,2), nx, ny, 1)
					rmidx = Region(old_div(nx,2), 0, 0, 1, ny, nz)
					rmidy = Region(0, old_div(ny,2), 0, nx, 1, nz)
					if options.nslices:
						rmidz = Region(0, 0, old_div(nz,2)-int(ceil(options.nslices/2.0)), nx, ny, options.nslices)
						rmidx = Region(old_div(nx,2)-int(ceil(options.nslices/2.0)), 0, 0, options.nslices, ny, nz)
						rmidy = Region(0, old_div(ny,2)-int(ceil(options.nslices/2.0)), 0, nx, options.nslices, nz)

					#tz = Transform({'type':'eman','az':0,'alt':0,'phi':0})

					regions={0:rmidz,1:rmidx,2:rmidy}
					#k=0
					for kk in regions:
						z=1
						if kk == 0:
							x=nx
							y=ny
					
						elif kk == 1:
							x=nz
							y=ny
						
						elif kk == 2:
							x=nx
							y=nz
					
						#if options.threed2threed or options.threed2twod:
						#d = EMData()
						#d.read_image(options.input, 0, False, regions[tag])
					
						if options.verbose > 6: print("I have extracted this orthogonal region".format(regions[kk]))
						slices = app.get_clip(regions[kk])

						prj = slices.project("standard",t)
						prj.set_attr('xform.projection',t)

						slicesf=prj.copy()

						slicesf.set_size(x,y,1)
						if options.verbose > 6: print("ptcltag={}, slicesf={}, type(slicesf)={}".format(ptcltag,slicesf,type(slicesf)))
					
						outname = options.path + '/' + os.path.basename( f ).replace(ext,ptcltag+'_SLICESortho.hdf')
						
						slicesf.write_image( outname,kk)
												
						if options.verbose > 6: print("mean={}, index={}".format( slicesf['mean'],kk) )
					
				if options.allz:
					if options.verbose > 5: print("Generating all z slices")
					outname = options.path + '/' + os.path.basename( f ).replace(ext,ptcltag+'_SLICESz.hdf')
					os.system('e2proc2d.py ' + f + ' ' + outname + ' --threed2twod')
			
				if options.allx:
					if options.verbose > 5: print("Generating all x slices")
					Tx = Transform({'type':'eman','az':0,'alt':-90,'phi':0})
					volx = app.copy()
					volx.transform(Tx)
					rotvolxname = options.path + '/' + os.path.basename( ft ).replace(ext, ptcltag+'rotx.hdf')
					volx.write_image(rotvolxname,0)
				
					outname = options.path + '/' + os.path.basename( f ).replace(ext,ptcltag+'_SLICESx.hdf')
				
					os.system('e2proc2d.py ' + rotvolxname + ' ' + outname + ' --threed2twod')
			
				if options.ally:	
					if options.verbose > 5: print("Generating all y slices")
					Ty = Transform({'type':'eman','az':0,'alt':-90,'phi':-90})
					voly = app.copy()
					voly.transform(Ty)
					rotvolyname = options.path + '/' + os.path.basename( f ).replace(ext, ptcltag+'roty.hdf')
					voly.write_image(rotvolyname,0)
				
					outname = options.path + '/' + os.path.basename( f ).replace(ext,ptcltag+'_SLICESy.hdf')
				
					os.system('e2proc2d.py ' + rotvolyname + ' ' + outname + ' --threed2twod')
	
	E2end(logger)
	
	elapsed = time.perf_counter() - start
	print(str(timedelta(seconds=elapsed)))

	return()	

if __name__ == '__main__':
	main()
