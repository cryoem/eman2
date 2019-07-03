#!/usr/bin/env python
# Author: Jesus Galaz-Montoya 2018 (jgalaz@gmail.com); last update Feb/12
# Copyright (c) 2000-2011 Baylor College of Medicine
#
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
from __future__ import print_function
from __future__ import division
from EMAN2_utils import *
from EMAN2 import *
import os
import sys
from sys import argv
from EMAN2jsondb import JSTask,jsonclasses

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] 
	This program stacks DDD frames saved as independent files but collected as part of a tomographic tiltseries, provided that the tilt angle
	information is in the filename between '[' and ']'; e.g., "afsdas_img1234_[-60.0].mrc", "afsdasd_img1235_[-60.0].mrc", "asdfasd_img1236_[-60.0].mrc", 
	... "afsdasfad_imgN_[-58.0].mrc", etc.
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument('--apix', type=float, default=0.0, help="""default=0.0 (not used). Reset the apix on the header of each image to this value.""")
	
	parser.add_argument('--extension', type=str, default='.mrc', help="""default=.mrc. Valid choices are .mrcs, .hdf, .tiff""")

	parser.add_argument('--inputstring', type=str, default='.mrc', help="""default=.mrc. String common to all images to be processed. E.g., with '.mrc' as the default, all files in the directory ending in '.mrc' (or containing this as part of the filename) will be analyzed.""")

	parser.add_argument('--nframes', type=int, default=0, help="""default=0 (not used). Number of expected frames per tilt angle, used to check whether errors occur in the stacking of images.""")

	parser.add_argument('--outputstem', type=str, default='stack', help="""default=stack. Output filename stem/root to write stacked frames out to; e.g., "stack_01", "stack_02", "stack_03", etc. Explicitly choose an extension/format by supplying --extension; e.g., --extension=.hdf, and the output will be 'stack_01_angle+02.hdf','stack_02_angle+00.hdf','stack_03_angle-02.hdf', etc, etc.""")
	#parser.add_argument("--parallel","-P",type=str,help="Run in parallel, specify type:<option>=<value>:<option>:<value>",default=None, guitype='strbox', row=8, col=0, rowspan=1, colspan=2, mode="align")

	parser.add_argument("--ppid", type=int, default=-1, help="Set the PID of the parent process, used for cross platform PPID")

	parser.add_argument('--twodstack', action='store_true', default=False, help="""default=False (not used). Produces a stack of 2D images (Z number of images of X,Y size under the same "virutal stack" file name). By default, the output stack will be a 3D stack (a single image of size X,Y,Z, where Z is equal to the number of stacked frames).""")

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness.")	

	(options, args) = parser.parse_args()
	
	logid=E2init(sys.argv,options.ppid)
	
	extensions=['.hdf','.mrc','.mrcs','.st','.ali','.rec','.tiff','.dm3']

	if options.extension not in extensions:
		print("\nERROR: invalid extension={}, must match one of the following extensions={}".format(extension,extensions))
		sys.exit(1)

	if options.twodstack:
		if options.extension not in ['.mrcs','.hdf']:
			print("\nERROR: --twodstack requires .mrcs or .hdf as --extension. EXITING.")
			sys.exit(1)

	c=os.getcwd()
	findir=os.listdir(c)

	outputroot = options.outputstem + options.extension

	frames=[]
	angles=set([])
	masterdict = {}
	for f in findir:
		fextension = os.path.splitext(f)[-1]
		if options.inputstring in f and fextension in extensions:
			imgfile = f
			angle = float(f.split('[')[1].split(']')[0])
			angles.add(angle)

			try:
				frames = masterdict[ angle ]
			except:
				masterdict.update({ angle:[] })
				frames = []

			framenumber = int(f.split('-')[-1].replace(fextension,''))
			frames.append([framenumber,f])

			masterdict.update({ angle: frames })

	angles = list(angles)
	angles.sort()

	angletag=''
	ntilts = len(angles)
	if options.verbose:
		print("\nthere are these many tiltangles, n={}".format(ntilts))

	if options.verbose > 8:	
		print("\nwhich are:\n{}".format(angles))

	for angle in angles:
		frames = masterdict[angle]
		nframes = len(frames)
		if options.nframes:

			if nframes != options.nframes:
				print("\nWARNING: the number of frames={} for tiltangle={} does not match the expected number of frames={}".format(nframes,angle,options.nframes))
				if nframes > options.nframes:
					print("\nin fact, there are MORE frames than expected")
				elif nframes < options.nframes:
					print("\nin fact, there are FEWER frames than expected")			
				
				print("\nERROR: exiting")
				sys.exit(1)

	count=0
	for angle in angles:
		print("\nstacking frames for angle {}/{}, angle={}".format(count+1,ntilts,angle))
		if float(angle) < 0.0:
			angletag = 'angle_' + str(angle).zfill(5) + '_stacked'
		elif float(angle) >= 0.0:
			angletag = 'angle_+' + str(angle).zfill(4) + '_stacked'

		output = outputroot.replace(options.extension,'_' + angletag + options.extension)
		
		if options.verbose > 5:
			print("\noutput image will be {}".format(output))
		
		frames = masterdict[angle]
		frames.sort()

		nframes = len(frames)
		#nimgs = len(frames)
		if options.nframes:
			if nframes != options.nframes:
				print("\nWARNING: nframes={} != options.nframes={}".format(nframes,options.nframes))
				print("\nbefore writing out images, frames are {}",format(frames))
				print("\nERROR: exiting")
				sys.exit(1)

		#if not options.twodstack:
		vol = EMData()
		hdr = EMData(frames[0][-1],0,True)
		vol.set_size(hdr['nx'], hdr['ny'], nframes)

		#print("\noutput will be {}".format(output))

		k=0
		for frame in frames:
			framenumber = int(frame[0])
			imgfile = frame[-1]
			if options.verbose:
				print("\ninput imgfile is {}".format(imgfile))
			
			img = EMData(imgfile,0)
			if options.apix:
				img['apix_x']=options.apix
				img['apix_y']=options.apix
				img['apix_z']=options.apix

			#img.write_image(output, framenumber)
			
			if options.twodstack:
				img.write_image(output, -1)
				if options.verbose:
					print("\nwrote frame {}/{} to output 2D stack".format(k+1,nframes))
			else:
				vol.insert_clip(img, (0, 0, k))
				if options.verbose:
					print("\ninserted frame {}/{} into output 3D volume".format(k+1,nframes))
		
			#cmd = 'e2proc2d.py ' + imgfile + ' ' + output
			#runcmd(options, cmd)

			print("\nadded frame {}".format(framenumber))
			k+=1

		if not options.twodstack:
			if options.apix:
				vol['apix_x']=options.apix
				vol['apix_y']=options.apix
				vol['apix_z']=options.apix
			vol.write_image(output,0)

		count+=1

	E2end(logid)

	return


if __name__ == '__main__':
	main()

