#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

# LAST update: July/2016
# Author: Steven Ludtke  2/8/2011 (rewritten)
# Author: Jesus Galaz-Montoya, all command line functionality + updates/enhancements/fixes.
# Author: John Flanagan  9/7/2011 (helixboxer)
# Copyright (c) 2011- Baylor College of Medicine
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
#


from past.utils import old_div
from builtins import range
import sys
import os
import weakref
import threading
import math
from EMAN2 import *

'''
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import Qt
from eman2_gui.emapplication import get_application, EMApp
from eman2_gui.emimage2d import EMImage2DWidget
from eman2_gui.emimagemx import EMImageMXWidget
from eman2_gui.emimage3d import EMImage3DWidget
from eman2_gui.emscene3d import EMScene3D
from eman2_gui.emdataitem3d import EMDataItem3D, EMIsosurface
from eman2_gui.emshape import EMShape
'''
#from eman2_gui.valslider import *


from sys import argv
from EMAN2jsondb import JSTask,jsonclasses


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <Volume file>

	WARNING: This program still under development.

	Tomography 3-D particle picker and annotation tool. Still under development."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_header(name="tbheader", help='Options below this label are specific to e2spt_boxer', title="### e2spt_boxer options ###", row=1, col=0, rowspan=1, colspan=3, mode="boxing")
	parser.add_pos_argument(name="tomogram",help="The tomogram to use for boxing.", default="", guitype='filebox', browser="EMTomoDataTable(withmodal=True,multiselect=False)",  row=0, col=0, rowspan=1, colspan=3, mode="boxing")
	parser.add_argument("--boxsize","-B",type=int,help="Box size in pixels",default=0)

	parser.add_argument("--centerbox", action="store_true", default=False, help="""DEPRECATED.
		Will apply xform.centerofmass to the boxed subvolumes before the final extraction.""")


	parser.add_argument("--autocenter", type=str, default='', help="""Options are
		--autocenter=xform.centerofmass (applies center of mass after --autocentermask if specified)
		or --autocenter=xform.centeracf (applies autoconvolution after --autocentermask if sepcified).
		This is intended to yield raw subtomograms that are better centered in the box.
		Mostly aplicable for freestanding, soluble, globular single particles.""")

	parser.add_argument("--autocentermask",type=str,default='',help="""Mask used to autocenter particles.
		Default is --autocentermask=None""")


	parser.add_argument("--path",default='',type=str,help="Pathname to save data to")
	parser.add_argument("--inmemory",action="store_true",default=False,help="This will read the entire tomogram into memory. Much faster, but you must have enough ram !", guitype='boolbox', row=2, col=1, rowspan=1, colspan=1, mode="boxing")
	parser.add_argument("--yshort",action="store_true",default=False,help="This means you have a file where y is the short axis", guitype='boolbox', row=2, col=0, rowspan=1, colspan=1, mode="boxing")
	parser.add_argument("--apix",type=float,default=0.0,help="Use THIS A/pix value to display the tomogram (if filtering) and to write to the header of the extracted subvolumes, instead of using the apix value one stored in the tomogram's header.", guitype='floatbox', row=3, col=0, rowspan=1, colspan=1, mode="boxing['self.pm().getAPIX()']")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--helixboxer",action="store_true",default=False,help="Helix Boxer Mode", guitype='boolbox', row=2, col=2, rowspan=1, colspan=1, mode="boxing")
	parser.add_argument("--normproc",type=str,help="""Normalization processor applied to particles before extraction. Default is None.
													Use --normproc=normalize, --normproc=normalize.edgemean or --normalize.mask, depending on your specimen and purposes.
													If using the latter, you must provide --masknorm, otherwise, a default --masknorm=mask.sharp:outer_radius=-2 will be used.""", default='')
	parser.add_argument("--masknorm",type=str,help="Mask used to normalize particles before extraction. Default is mask.sharp:outer_radius=-2", default='')
	parser.add_argument("--thresh",type=float,help="Threshold particles before writing them out to get rid of too high and/or too low pixel values.", default=0.0)


	#parser.add_argument('--bin', type=int, default=1, help="""Specify the binning/shrinking factor you want to use (for X,Y and Z) when opening the tomogram for boxing. \nDon't worry, the sub-volumes will be extracted from the UNBINNED tomogram. \nIf binx, biny or binz are also specified, they will override the general bin value for the corresponding X, Y or Z directions""", guitype='intbox', row=3, col=1, rowspan=1, colspan=1, mode="boxing")

	parser.add_argument('--shrink', type=int, default=1, help="""Specify the binning/shrinking factor you want to use (for X,Y and Z) when opening the tomogram for boxing. \nDon't worry, the sub-volumes will be extracted from the UNBINNED tomogram.""", guitype='intbox', row=3, col=1, rowspan=1, colspan=1, mode="boxing")

	parser.add_argument("--lowpass",type=int,help="Resolution (integer, in Angstroms) at which you want to apply a gaussian lowpass filter to the tomogram prior to loading it for boxing",default=0, guitype='intbox', row=3, col=2, rowspan=1, colspan=1, mode="boxing")
	parser.add_argument("--preprocess",type=str,help="""A processor (as in e2proc3d.py) to be applied to the tomogram before opening it. \nFor example, a specific filter with specific parameters you might like. \nType 'e2proc3d.py --processors' at the commandline to see a list of the available processors and their usage""",default=None)

	#parser.add_argument('--reverse_contrast', action="store_true", default=False, help='''This means you want the contrast to me inverted while boxing, AND for the extracted sub-volumes.\nRemember that EMAN2 **MUST** work with "white" protein. You can very easily figure out what the original color\nof the protein is in your data by looking at the gold fiducials or the edge of the carbon hole in your tomogram.\nIf they look black you MUST specify this option''', guitype='boolbox', row=4, col=0, rowspan=1, colspan=1, mode="boxing")

	parser.add_argument('--invert', action="store_true", default=False, help='''This means you want the contrast to me inverted while boxing, AND for the extracted sub-volumes.\nRemember that EMAN2 **MUST** work with "white" protein. You can very easily figure out what the original color\nof the protein is in your data by looking at the gold fiducials or the edge of the carbon hole in your tomogram.\nIf they look black you MUST specify this option''', guitype='boolbox', row=4, col=0, rowspan=1, colspan=1, mode="boxing")

	#parameters for commandline boxer

	parser.add_argument('--coords', type=str, default='', help='Provide a coordinates file that contains the center coordinates of the sub-volumes you want to extract, to box from the command line.')

	parser.add_argument('--cshrink', type=int, default=1, help='''WARNING: The coordinates file gets written on the scale of the input tomogram. This means that if you supply the raw, unbinned tomogram to e2pt_boxer but also say --shrink at the command line, a temporary tomogram will be created for particle localization, but the coordinates file will contain the full-size coordinates. If, on the other hand, you pre-shrink the tomogram with another tool (e2proc3d.py or binvol in IMOD), the coordinates in the coordinates file will be shrunk. It is in such case that you would use --cshrink to specifies the factor by which to multiply the coordinates in the coordinates file, so that they can be at the same scale as the RAW tomogram (or whatever tomogram you intend for the particles to be extracted from).\nFor example, provide --cshrink=2 if the coordinates are on a 2K x 2K scale because you used a 2K x 2K tomogram to find the particles,\nbut you want to extract the subvolumes from a UN-shrunk 4K x 4K tomogram.''')

	parser.add_argument('--subset', type=int, default=0, help='''Specify how many sub-volumes 
		from the coordinates file you want to extract; e.g, if you specify 10, the first 10 
		particles will be boxed.\n0 means "box them all" because it makes no sense to box none''')
		
	parser.add_argument('--output', type=str, default='stack.hdf', help="Specify the name of the stack file where to write the extracted sub-volumes")
	parser.add_argument('--output_format', type=str, default='stack', help='''Specify 'single' if you want the sub-volumes to be written to individual files. You MUST still provide an output name in the regular way.\nFor example, if you specify --output=myparticles.hdf\nbut also specify --output_format=single\nthen the particles will be written as individual files named myparticles_000.hdf myparticles_001.hdf...etc''')

	parser.add_argument('--bruteaverage', action="store_true", default=False, help='Will generate an average of all the subvolumes (no alignment done). This is useful to see if, on average, particles are the desired specimen and reasonably centered')

	#parser.add_argument('--swapyz', action="store_true", default=False, help='''This means that the coordinates file and the actual tomogram do not agree regarding which is the "short" direction.\nFor example, the coordinates file migh thave a line like this:\n1243 3412 45\nwhere clearly the "short" direction is Z; yet, if in the actual tomogram the short direction is Y, as they come out fromIMOD by default, then the line should have been:\n1243 45 3412\n''')
	#parser.add_argument('--normalize', action="store_true", default=False, help='Will normalize each subvolume so that the mean is zero and standard deviation one').
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	parser.add_argument("--parallel",type=str, default=None, help="""Used only when extracting subtomograms from the commandline (NOT from the GUI interface) if --coords is provided. Default=None (not used). Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel""")


	global options,args
	(options, args) = parser.parse_args()


	if options.yshort:
		print("\nERROR: --yshort is not supported. Rotate your tomogram so that Z is the shortest side (e.g., using IMOD from the commandline, type 'clip rotx tomogram.rec tomogram_ZSHORT.rec'")
		sys.exit()
	else:
		tomohdr=EMData(args[0],0,True)
		nz=tomohdr['nz']
		ny=tomohdr['ny']
		if nz>ny:
			print("\nERROR: nz=%d is larger than ny=%d in the tomgram. Rotate the tomogram so that Z is the shortest side (e.g., using IMOD from the commandline, type 'clip rotx tomogram.rec tomogram_ZSHORT.rec'" %(nz,ny))



	if options.normproc:
		options.normproc=parsemodopt(options.normproc)

	if options.masknorm:
		options.masknorm=parsemodopt(options.masknorm)

	if options.thresh:
		options.thresh=parsemodopt(options.thresh)

	if options.autocentermask:
		options.autocentermask=parsemodopt(options.autocentermask)

	options.centerbox = options.autocenter

	#print "\n\n\n\n\nAUTOCENTER IS", options.autocenter
	#print "Therefore center is", options.centerbox


	'''
	Create the path where subtomograms will be saved
	'''

	from e2spt_classaverage import sptmakepath

	if options.path:
		options = sptmakepath( options, 'spt_boxer')

	if len(args) != 1:
		parser.error("You must specify a single volume data file on the command-line.")
	if not file_exists(args[0]):
		parser.error("%s does not exist" %args[0])

	if options.coords:
		#commandline_tomoboxer(args[0],options.coords,options.subset,options.boxsize,options.cshrink,options.output,options.output_format,options.swapyz,options.invert,options.centerbox)

		logger = E2init(sys.argv, options.ppid)

		commandline_tomoboxer(args[0],options)

		cleanstack(options)

		print("\ncomputing bruteaverage")
		cmdavg = 'e2proc3d.py ' + options.output + ' ' + options.output.replace('.hdf','__bruteavg.hdf') + ' --average'
		
		retavg=runcmd( options, cmdavg )
		if retavg:
			print("done")


		E2end(logger)

	else:
		sptboxergui(options,args)
		'''
		from PyQt4 import QtCore, QtGui
		from PyQt4.QtCore import Qt
		from eman2_gui.emapplication import get_application, EMApp
		from eman2_gui.emimage2d import EMImage2DWidget
		from eman2_gui.emimagemx import EMImageMXWidget
		from eman2_gui.emimage3d import EMImage3DWidget
		from eman2_gui.emscene3d import EMScene3D
		from eman2_gui.emdataitem3d import EMDataItem3D, EMIsosurface
		from eman2_gui.emshape import EMShape


		if options.path and options.verbose:
			print "OPTIONS.PATH IS!!!!\n\n\n", options.path

		img = args[0]

		app = EMApp()
		if options.inmemory:


			imghdr = EMData(img,0,True)
			thisapix = imghdr['apix_x']

			if options.shrink and int(options.shrink) > 1:

				# The new shrinking scheme
				print "Shrinking, please wait :)"
				img = EMData()
				img.read_binedimage(args[0],0,options.shrink)
			else:
				print "You better have A LOT of memory (more than 8GB) if this is an un-shrunk 4k x 4k x 0.5k tomogram, because I'm loading it un-shrunk. It is wise to wait patiently."
				img = EMData(img,0)

			if options.invert:
				img = img*(-1)
				print "The contrast of the tomogram has been reversed"

			if options.lowpass:
				filt=1.0/options.lowpass
				print "The tomogram is being low pass filtered to %d Angstroms resolution" %(options.lowpass)
				img = img.process('filter.lowpass.gauss',{'cutoff_freq':filt,"apix":img['apix_x']})

			print "Done !"

			if options.apix:
				thisapix = options.apix
			
			box = 32
			if options.boxsize:
				box = options.boxsize
			
			boxer = EMTomoBoxer(app,data=img,datafile=None,yshort=options.yshort,apix=thisapix,boxsize=box,shrink=options.shrink,contrast=options.invert,center=options.centerbox,mod=None,normalize=options.normproc)
		else :
	#		boxer=EMTomoBoxer(app,datafile=args[0],yshort=options.yshort,apix=options.apix,boxsize=options.boxsize)
			img=args[0]

			imghdr = EMData(img,0,True)
			thisapix = imghdr['apix_x']

			"""The modd variable is used as a check that determines whether a "modified"
			tomogram needs to be deleted prior to extracting boxes from disk.
			Because boxing from disk does NOT open the whole tomogram, a modified copy
			needs to be generated and written to file when you want to find
			particles in a shrunk or pre-low pass filtered tomogram, BUT still extract
			from the raw tomogram"""
			modd = False
			if options.shrink > 1:

				imgnew = img
				if '_editedtemp.' not in img:
					imgnew = img.replace('.','_editedtemp.')
				print "Shrinking, please wait :)"
				imgfile = EMData()
				imgfile.read_binedimage(img,0,options.shrink)
				imgfile.write_image(imgnew)

				img = imgnew
				modd = True

			if options.invert:
				imgnew = img
				if '_editedtemp.' not in img:
					os.system('rm *_editedtemp*')
					imgnew = img.split('/')[-1].replace('.','_editedtemp.')

				cmd = 'e2proc3d.py ' + img + ' ' + imgnew + ' --mult=-1'
				os.system(cmd)
				img = imgnew
				modd = True

			if options.lowpass:
				imgnew = img
				if '_editedtemp.' not in img:
					os.system('rm *_editedtemp*')
					imgnew = img.split('/')[-1].replace('.','_editedtemp.')
				filt=1.0/options.lowpass

				imghdr = EMData(img,0,True)
				tapix = imghdr['apix_x']
				if options.apix:
					tapix = options.apix

				cmd = 'e2proc3d.py ' + img + ' ' + imgnew + ' --process=filter.lowpass.gauss:cutoff_freq=' + str(filt) + ':apix=' + str(tapix)
				os.system(cmd)
				img = imgnew
				modd = True

			#print "The shrink factor default is", options.shrink


			print "\nDatabfile and type of datafile are", img, type(img)



			if options.apix:
				thisapix = options.apix
			box = 32
			if options.boxsize:
				box = options.boxsize
			
			boxer=EMTomoBoxer(app,data=None,datafile=img,yshort=options.yshort,apix=thisapix,boxsize=box,shrink=options.shrink,contrast=options.invert,center=options.centerbox,mod=modd,normalize=options.normproc)

			#boxer=EMTomoBoxer(app,datafile=img,yshort=options.yshort,apix=options.apix,boxsize=options.boxsize,shrink=options.shrink,contrast=options.invert,center=options.centerbox,mod=None,normalize=options.normproc)

		boxer.show()
		app.execute()
		'''
	return


"""
Function to check whether all boxes in the extracted stack seem "healthy" and if not eliminate the bad ones
"""
def cleanstack(options):
	n=EMUtil.get_image_count(options.output)
	badptcls = []
	print("\n(e2spt_boxer)(cleanstack) checking for sanity of output stack %s" %(options.output))
	for i in range(n):
		try:
			ptcl=EMData(options.output,i)
			if float(ptcl['sigma']) == 0.0:
				badptcls.append(i)
				print("WARNING: bad particle %d will be removed from output stack %s" %(i,options.output))

		except:
			print("WARNING: bad particle %d will be removed from output stack %s" %(i,options.output))
			badptcls.append(i)

	if badptcls:
		print("\n%d bad particles identified"%(len(badptcls)))

		tmpoutput = options.output.replace('.hdf','_tmp.hdf')
		badptcls.sort()
		first=0
		#kl=0
		nbad=len(badptcls)

		for kl in range(len(badptcls)):
			print("pruning bad ptcl", badptcls[kl])

			if kl > 0:
				first=badptcls[kl-1]+1
			
			last=badptcls[kl]-1
			
			cmd ='e2proc3d.py ' + options.output + ' ' + tmpoutput + ' --first ' + str(first) + ' --last ' + str(last)
			
			if kl == nbad-1:
				
				if last != n-1 and last < n-1:
					cmd+= ' && e2proc3d.py ' + options.output + ' ' + tmpoutput + ' --first ' + str(last+2) + ' --last ' + str(n-1) + ' --append'
	
			if kl>0 and 'append' not in cmd:
				cmd += ' --append'

			print("cmd is", cmd)
			runcmd(options,cmd)

			kl+=1

		os.rename(tmpoutput,options.output)

	return 


"""
This function is called to extract sub-volumes from the RAW tomogram, regardless
of where their coordinates are being found (the tomogram to find the coordinates might be
shrunk and/or lowpass filtered).
It is also called when boxing from the commandline, without GUI usage, as when you already have
a coordinates file
"""
def unbinned_extractor(options,x,y,z,tomogram,coordindx=None):
	#print "inside unbinned extractor for particle", coordindx
	boxsize = options.boxsize
	#print "boxsize",boxsize

	cshrink = options.cshrink
	#print "cshrink",cshrink

	invert = options.invert
	#print "invert",invert

	center = options.centerbox
	#print "center",center


	#if options.verbose:
	#	print "\n\nUnbinned extractor received this center", center

	print("\n(e2spt_boxer.py)(unbinned_extractor) reading tomogram header from %s" %(tomogram))
	tomo_header=EMData(tomogram,0,True)
	apix = tomo_header['apix_x']
	#print "done"
	
	#if options.verbose:
	#	print "Which has a size of", tomo_header['nx'],tomo_header['ny'],tomo_header['nz']
	
	x=round(x*cshrink)
	y=round(y*cshrink)
	z=round(z*cshrink)

	#if options.verbose:
	#	print "The actual coordinates used for extraction are", x, y, z

	#print "defining region"
	r = Region(old_div((2*x-boxsize),2),old_div((2*y-boxsize),2), old_div((2*z-boxsize),2), boxsize, boxsize, boxsize)
	#print "r",r
	e = EMData()
	e.read_image(tomogram,0,False,r)
	#print "extracted particle",coordindx

	#IF the boxed out particle is NOT empty, perform BASIC RAW-PARTICLE EDITING: contrast reversal and normalization
	#Sometimes empty boxes are picked when boxing from the commandline if yshort isn't specified but should have,
	#or if erroneous binning factors are provided

	#print "read particle from tomogram %s using coords x=%d,y=%d,z=%d" %(tomogram,x,y,z)
	#print "has mean, meannonzero, and sigma",e['mean'],e['mean_nonzero'],e['sigma']

	if float(e['sigma']) != 0.0:
		#print "mean is nonzero",e['mean']
		

		e['xform.align3d'] = Transform() #Make sure the default alignment parameters are zero

		'''
		Attempt to center particles. They must be masked to ensure other particles will not cause shifts in center of mass
		'''
		if options.verbose:
			print("\n\n\n\nCENTER is",center)

		if center:
			ec = e.copy()
			ec = ec*-1

			if options.autocentermask:
				if options.verbose : print("\nMasking for autocentering")
				ec.process_inplace(options.autocentermask[0],options.autocentermask[1])

			ec.process_inplace('normalize')

			if options.autocentermask:
				if options.verbose : print("\nMasking for autocentering")
				ec.process_inplace(options.autocentermask[0],options.autocentermask[1])

			ec.process_inplace("threshold.belowtozero",{'minval':0.0})

			if center == 'xform.centerofmass':
				ec.process_inplace('xform.centerofmass')
				if options.verbose : print("\nApplying center of mass")

			elif center == 'xform.centeracf':
				ec.process_inplace('xform.centeracf')
				if options.verbose : print("\nApplying xform.centeracf")


			#rad = boxsize/2.0+1.0
			#ec = e.process('mask.sharp',{'outer_radius':rad})
			#ec.process_inplace('xform.centerofmass')

			trans = ec['xform.align3d'].get_trans()
			tx = trans[0]
			ty = trans[1]
			tz = trans[2]

			if options.verbose : 
				print("\nAutocentering translations are", tx,ty,tz)
				print("\n\nThe old coordinates were", x, y, z)

			x = x - tx
			y = y - ty
			z = z - tz
			if options.verbose : print("Thus the new ones are", x, y, z)

			r = Region(old_div((2*x - boxsize),2),old_div((2*y - boxsize),2), old_div((2*z - boxsize),2), boxsize, boxsize, boxsize)
			e = EMData()
			e.read_image(tomogram,0,False,r)

		#It IS CONVENIENT to record any processing done on the particles as header parameters
		#print "writing particle header parameters"

		e['ptcl_source_image'] = os.path.basename(tomogram)
		e['ptcl_source_coord'] = (x,y,z)
		e['spt_tomogram'] = tomogram
		e['spt_originalstack'] = options.output


		#Make sure the transform parameter on the header is "clean", so that any later processing transformations are meaningful
		e['xform.align3d'] = Transform({"type":'eman','az':0,'alt':0,'phi':0,'tx':0,'ty':0,'tz':0})
		#print "done"
		if options.verbose : 
			print("The extracted particle has this boxsize", e['nx'],e['ny'],e['nz'])
			print("And the following mean BEFORE normalization", e['mean'])

		if options.normproc:
			#print "normalizing"
			if options.verbose: 
				print("WARNING! particle being normalized!")
			
			if options.normproc[0]=="normalize.mask":
				mask=EMData(e["nx"],e["ny"],e["nz"])
				mask.to_one()
				if options.masknorm:
					mask.process_inplace(options.masknorm[0],options.masknorm[1])
				options.normproc[1]["mask"]=mask

			e.process_inplace(options.normproc[0],options.normproc[1])
			e['spt_normalization'] = str(options.normproc[0])+' '+str(options.normproc[1])
			if options.verbose : print("This is the mean AFTER normalization", e['mean'])
			#print "done"
		if invert:
			#print "inverting contrast"
			if options.verbose : print("Particle has the following mean BEFORE contrast inversion", e['mean'])
			if options.verbose : print("Inverting contrast because --invert is", invert)
			e=e*-1
			if options.verbose : print("Particle has the following mean AFTER contrast inversion", e['mean'])
			#print "done"
		if options.thresh:
			if options.verbose : print("The thresh to apply is", options.thresh)
			e.process_inplace(options.thresh[0],options.thresh[1])

		#print "projecting"
		prjT = Transform({'type':'eman','az':0,'alt':0,'phi':0})
		
		#if options.yshort:
		#	prjT = Transform({'type':'eman','az':-90,'alt':-90,'phi':0})

		prj = e.project("standard",prjT)
		#print "done"
		apix = e['apix_x']
	
		if options.verbose : 
			print("\nThere was a particle successfully returned, with the following box size and mean value")
			print(e['nx'],e['ny'],e['nz'])
			print(e['mean'])

		
		e['apix_x'] = apix
		e['apix_y'] = apix
		e['apix_z'] = apix
		
		#The origin WILL most likely be MESSED UP if you don't explicitely set it to ZERO.
		#This can create ANNOYING visualization problems in Chimera
		#print "resetting origin"
		e['origin_x'] = 0
		e['origin_y'] = 0
		e['origin_z'] = 0
		#print "done"
		if options.coords:
			#if options.verbose:
			#	print "writing particle %d with size %d,%d,%d to temporary output file %s" % (coordindx,e['nx'],e['ny'],e['nz'],'sptboxer_dummy.hdf')
			
			if options.parallel:
				e.write_image('sptboxer_dummy.hdf',coordindx)
			else:
				e.write_image(options.output,-1)


			
			print("extracted particle %d to temporary file" %(coordindx))
		
		c=os.getcwd()
		findir=os.listdir(c)
		if prj:

			prj.set_attr('xform.projection',prjT)
			prj['apix_x']=apix
			prj['apix_y']=apix
			
			if not options.coords:
				nameprjs = options.output
				if '.mrc' in options.output:
					nameprjs = nameprjs.replace('.mrc','.hdf')

				#if not options.yshort:
				nameprjs = nameprjs.replace('.hdf','__prjsz.hdf')
				#elif options.yshort:
				#	nameprjs = nameprjs.replace('.hdf','__prjsy.hdf')

				if options.path not in nameprjs:
					nameprjs = options.path + '/' + nameprjs
				prj['spt_originalstack']=nameprjs.split('/')[-1]
				prj.process_inplace('normalize')
				prj.write_image(nameprjs,-1)
			
			elif options.coords:
				if options.parallel:
					prj.write_image('sptboxer_dummy_prjs.hdf',coordindx)
				else:
					prj.write_image(options.output.replace('.hdf','__prjsz.hdf'),-1)

		return(e,prj)

	else:
		if not coordindx:
			coordindx==-1
		print("""\nWARNING! particle %d at coordinates x=%d, y=%d, z=%d. was skipped (and not boxed) because it's SIGMA was ZERO (suggesting the box was empty).
			Your coordinates file and/or the shrinking factors specified might be MESSED UP, or you might need to swap Y and Z
			""" %(coordindx,x,y,z))
		return None


"""
This function enables extracting sub-volumes from the command line, without opening the GUI.
Usually used when "re-extracting" sub-volumes (for whatever reason) from a coordinates file previously generated.
It allows for extraction of smaller sub-sets too.
"""
def commandline_tomoboxer(tomogram,options):

	if not options.boxsize:
		print("\nERROR: --boxsize required")
		sys.exit()
		
	clines = open(options.coords,'r').readlines()
	ncoords = len(clines)

	if options.subset:
		if options.subset > ncoords:
			print("WARNING: The total amount of lines in the coordinates files is LESS than the subset of particles to box you specified; therefore, ALL particles will be extracted")
		else:
			ncoords=options.subset

	print("The size of the set of sub-volumes to extract is", ncoords)

	k=-1
	#name = options.output
	if options.path and options.path not in options.output:
		options.output = options.path + '/' + options.output

	if ".hdf" not in options.output:
		print("Format ERROR: Only .hdf fomart supported.")
		sys.exit()
	
	#avgr=None
	#if options.bruteaverage:
	#	avgr=Averagers.get('mean.tomo')

	jj=0
	
	coordsdict = {}
	xs = []
	ys = []
	zs = []
	
	apix = EMData( tomogram, 0, True )['apix_x']
	if options.apix:
		apix = options.apix
	
	c = os.getcwd()
	findir = os.listdir( c )
	
	try:
		os.remove( 'sptboxer_dummy.hdf' )
		os.remove( 'sptboxer_dummy_prjs.hdf' )
	except:
		pass


	apix = EMData(tomogram,0,True)['apix_x']
	if options.apix:
		apix=options.apix

	if options.parallel:
		dimg = EMData(8,8,8)
		dimg['apix_x']=apix
		dimg['apix_y']=apix
		dimg['apix_z']=apix
		dimg.to_one()

		dimg2D =EMData(8,8)
		dimg2D['apix_x']=apix
		dimg2D['apix_y']=apix
		dimg2D.to_one()

	for i in range(ncoords):

		#Some people might manually make ABERRANT coordinates files with commas, tabs, or more than one space in between coordinates
		clines[i] = clines[i].replace(", ",' ')
		clines[i] = clines[i].replace(",",' ')
		clines[i] = clines[i].replace("x",'')
		clines[i] = clines[i].replace("y",'')
		clines[i] = clines[i].replace("z",'')
		clines[i] = clines[i].replace("=",'')
		clines[i] = clines[i].replace("_",' ')
		clines[i] = clines[i].replace("\n",'')
		clines[i] = clines[i].replace("\t",' ')
		clines[i] = clines[i].replace("  ",' ')
		clines[i] = clines[i].split()

		x = int( float(clines[i][0]) )
		xs.append( x )
		
		y = int( float(clines[i][1]) )
		ys.append( y )
		
		z = int( float(clines[i][2]) )
		zs.append( z )
		
		coordsdict.update({i:[x,y,z]})
		
		if options.verbose: 
			print("The raw coordinates from the coordinates file provided for particle#%d are x=%d, y=%d, z=%d " % (i,x,y,z))

		#if options.swapyz:
		#	if options.verbose : print "You indicated Y and Z are flipped in the coords file, respect to the tomogram's orientation; therefore, they will be swapped"
		#	aux = y
		#	y = z
		#	z = aux
		#	if options.verbose : print "Therefore, the swapped coordinates are", x, y, z

		if options.parallel:
			dimg.write_image( 'sptboxer_dummy.hdf', i )
			dimg2D.write_image( 'sptboxer_dummy_prjs.hdf', i )

	
		if options.verbose:
			print("\n(e2spt_preproc)(main) wrote dummy ptcls to %s" %( 'sptboxer_dummy' ))
	
	if options.parallel:
		dummylarge = EMData(options.boxsize,options.boxsize,options.boxsize)
		dummylarge2D = EMData(options.boxsize,options.boxsize)

		dummylarge.to_one()
		dummylarge.write_image('sptboxer_dummy.hdf',0)
		
		dummylarge2D.to_one()
		dummylarge2D.write_image('sptboxer_dummy_prjs.hdf',0)


		print("\n(e2spt_preproc)(main) - INITIALIZING PARALLELISM!\n")
	
	
		from EMAN2PAR import EMTaskCustomer
		etc=EMTaskCustomer(options.parallel)
	#pclist=[tomogram]

	#etc.precache(tomogram)
	#print "\n(e2spt_preproc)(main) - precaching tomogram"

	tasks=[]
	results=[]
	counter=1

	results=None
	for coordindx in coordsdict:

		coordx=coordsdict[coordindx][0]
		coordy=coordsdict[coordindx][1]
		coordz=coordsdict[coordindx][2]
		
		if options.parallel:
			task = TomoBoxer3DTask( options, coordindx, coordx, coordy, coordz, tomogram)
			tasks.append(task)
			counter+=1
		else:
			print("\nthe tomogram file is", tomogram)
			results=unbinned_extractor( options, coordx, coordy, coordz, tomogram, coordindx )

		
	if options.parallel:	
		if tasks:
			tids = etc.send_tasks(tasks)
			if options.verbose: 
				print("\n(e2spt_preproc)(main) preprocessing %d tasks queued" % (len(tids))) 

	
			results = get_results( etc, tids, options )

	if results:

		if options.parallel:
			os.rename('sptboxer_dummy.hdf',options.output)
			nameprjs = options.output.replace('.hdf','__prjsz.hdf')
			os.rename('sptboxer_dummy_prjs.hdf',nameprjs)

		radius = old_div(options.boxsize,4.0)	#the particle's diameter is boxsize/2
		
		#if options.cshrink:
		#	radius /= options.cshrink
		
		print("\nextraction to output file completed")

		cmd = 'e2spt_icethicknessplot.py --plotparticleradii --fit --apix ' + str( apix ) + ' --radius ' + str( int(radius) ) + ' --files ' + options.coords
		
		if options.cshrink:
			cmd += ' --cshrink ' + str( int(options.cshrink) )
		
		print("\nplotting particle distribution")

		retice=runcmd( options, cmd )
		if retice:
			print("done")
		#if options.bruteaverage:
		
		if options.path:
			c=os.getcwd()
			findir=os.listdir(c)
			stemcoords=os.path.splitext(options.coords)[0]

			for fi in findir:
				if stemcoords in fi and '.png' in fi:
					os.rename(fi,options.path+'/'+fi)

		

	return


def runcmd(options,cmd):
	if options.verbose > 8:
		print("(e2spt_classaverage)(runcmd) running command", cmd)
	
	p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	text=p.communicate()	
	p.stdout.close()
	
	if options.verbose > 8:
		print("(e2spt_classaverage)(runcmd) done")
	
	#if options.verbose > 9:
	#	print text
	
	return 1


'''
CLASS TO PARALLELIZE PARTICLE EXTRACTION FROM COMMANDLINE
'''
class TomoBoxer3DTask(JSTask):
	'''This is a task object for the parallelism system.'''


	def __init__(self, options, coordindx, coordx, coordy, coordz, tomogramname):
	
		#data={"image":image}
		
		#JSTask.__init__(self,"TomoBoxer3d",data,{},"")
		JSTask.__init__(self,"TomoBoxer3d",{},"")
		self.classoptions={"options":options,"coordindx":coordindx, "coordx":coordx,"coordy":coordy,"coordz":coordz,"tomogramname":tomogramname}
	
	def execute(self,callback=None):
		
		options = self.classoptions['options']
		coordindx = self.classoptions['coordindx']
		coordx = self.classoptions['coordx']
		coordy = self.classoptions['coordy']
		coordz = self.classoptions['coordz']
		tomogramname = self.classoptions['tomogramname']
		#counter=self.classoptions['counter']
		#ncoords=self.classoptions['ncoords']
		#print "inside class TomoBoxer3DTask"
		#print "calling unbinned_extractor"
		unbinned_extractor(options,coordx,coordy,coordz,tomogramname,coordindx)
		#print "returned from unbinned_extractor"

		return


def get_results(etc,tids,options):
	"""This will get results for a list of submitted tasks. Won't return until it has all requested results.
	aside from the use of options["ptcl"] this is fairly generalizable code. """
	
	# wait for them to finish and get the results
	# results for each will just be a list of (qual,Transform) pairs
	results=[0]*len(tids)		# storage for results
	ncomplete=0
	tidsleft=tids[:]
	while 1:
		time.sleep(5)
		proglist=etc.check_task(tidsleft)
		nwait=0
		for i,prog in enumerate(proglist):
			if prog==-1 : nwait+=1
			if prog==100 :
				r=etc.get_results(tidsleft[i])				#Results for a completed task
				
				if r:
					#print "r is", r
					#ptcl=r[0].classoptions["ptclnum"]		#Get the particle number from the task rather than trying to work back to it
					#results[ptcl] = r[1]
					results[i] = 1					
				ncomplete+=1
		
		tidsleft=[j for i,j in enumerate(tidsleft) if proglist[i]!=100]		# remove any completed tasks from the list we ask about
		if options.verbose:
			print(("\n%d tasks, %d complete, %d waiting to start \r" % (len(tids),ncomplete,nwait)))
			sys.stdout.flush()
	
		if len(tidsleft)==0: break
		
	return results


def sptboxergui(options,args):

	from PyQt4 import QtCore, QtGui
	from PyQt4.QtCore import Qt
	from eman2_gui.emapplication import get_application, EMApp
	from eman2_gui.emimage2d import EMImage2DWidget
	from eman2_gui.emimagemx import EMImageMXWidget
	from eman2_gui.emimage3d import EMImage3DWidget
	from eman2_gui.emscene3d import EMScene3D
	from eman2_gui.emdataitem3d import EMDataItem3D, EMIsosurface
	from eman2_gui.emshape import EMShape
	from eman2_gui.valslider import ValSlider, ValBox

	class EMAverageViewer(QtGui.QWidget):
		"""This is a multi-paned view showing a single boxed out particle from a larger tomogram"""

		def __init__(self,parent):
			QtGui.QWidget.__init__(self)

			self.setWindowTitle("Particle Average")

			self.parent=weakref.ref(parent)

			self.resize(300,500)

			self.gbl = QtWidgets.QGridLayout(self)
			#self.xyview = EMImage2DWidget()
			#self.gbl.addWidget(self.xyview,0,1)

			#self.xzview = EMImage2DWidget()
			#self.gbl.addWidget(self.xzview,1,1)

			#self.zyview = EMImage2DWidget()
			#self.gbl.addWidget(self.zyview,0,0)

			# This code puts an isosurface in the lower left, but we seem to get a lot of segfaults in the isosurfaces. Replacing it with a 2-D view for now
			#self.d3view = EMScene3D()
			#self.d3viewdata = EMDataItem3D(test_image_3d(3), transform=Transform())
			#isosurface = EMIsosurface(self.d3viewdata, transform=Transform())
			#self.d3view.insertNewNode('', self.d3viewdata, parentnode=self.d3view)
			#self.d3view.insertNewNode("Iso", isosurface, parentnode=self.d3viewdata)

			self.d3view = EMImage2DWidget()
			self.gbl.addWidget(self.d3view,0,0)

			self.gbl2 = QtWidgets.QGridLayout()
			self.gbl.addLayout(self.gbl2,1,0)

			self.wfilt = ValSlider(rng=(0,50),label="Filter:",value=0.0)
			self.gbl2.addWidget(self.wfilt,2,0,1,2)

			self.wmask = ValSlider(rng=(0,100),label="Mask:",value=0.0)
			self.gbl2.addWidget(self.wmask,3,0,1,2)

			self.wsymlbl=QtWidgets.QLabel("Symmetry:")
			self.gbl2.addWidget(self.wsymlbl,4,0)

			self.wsym=QtGui.QLineEdit("C1")
			self.gbl2.addWidget(self.wsym,4,1)

			self.wprog=QtGui.QProgressBar()
			self.wprog.setRange(0,100)
			self.gbl2.addWidget(self.wprog,5,0,1,2)

			self.wrestart=QtGui.QPushButton("Restart")
			self.gbl2.addWidget(self.wrestart,6,1)

			self.needupd=0					# Set by the second thread when a display update is ready, 1 means progress update, 2 means volume update
			self.threadrestart=False		# Set by the GUI thread when the second thread needs to restart from scratch
			self.threadprog=0				# Thread progress (0-100)
			self.threadprogstr=""			# String describing thread action
			self.data=None

			# These are values from the widgets, stored so the thread can get at them without making GUI calls
			self.sym="c1"
			self.filt=0.0
			self.mask=0.0

			self.wfilt.valueChanged.connect(self.event_filter)
			self.wmask.valueChanged.connect(self.event_mask)
			self.wsym.editingFinished.connect(self.event_symchange)
			self.wrestart.clicked[bool].connect(self.event_restart)


			# The timer event handles displaying the results processed by the other thread
			self.timer=QtCore.QTimer(self)
			self.timer.timeout.connect(self.event_timer)
			self.timer.start(500)

			# The processing is all done in the background by the other thread
			self.bgthread=threading.Thread(target=self.thread_process)
			self.bgthread.daemon=True
			self.bgthread.start()

		def event_timer(self):
			if self.needupd&1 :
				self.wprog.setValue(self.threadprog)
			if self.needupd&2 : self.update()
			self.needupd=0

		def event_symchange(self):
			#print "sym"
			self.sym=self.wsym.text()
			self.wrestart.setEnabled(True)

		def event_filter(self,value):
			#print "filt"
			self.filt=value
			self.wrestart.setEnabled(True)

		def event_mask(self,value):
			#print "mask"
			self.mask=value
			self.wrestart.setEnabled(True)

		def event_restart(self):
			#print "restart"
			self.threadrestart=True
			self.wrestart.setEnabled(False)

		#def set_data(self,data):
			#"""Sets the current volume to display"""

			#self.data=data

			#self.update()
			#self.show()

		def update(self):
			#if self.wfilt.getValue()!=0.0 :
				#self.fdata=self.data.process("filter.lowpass.gauss",{"cutoff_freq":1.0/self.wfilt.getValue()})

			#xyd=self.fdata.process("misc.directional_sum",{"axis":"z"})
			#xzd=self.fdata.process("misc.directional_sum",{"axis":"y"})
			#zyd=self.fdata.process("misc.directional_sum",{"axis":"x"})

			#self.xyview.set_data(xyd)
			#self.xzview.set_data(xzd)
			#self.zyview.set_data(zyd)

			#self.d3viewdata.setData(self.data)
			#self.d3view.updateSG()
			self.d3view.set_data(self.data)

		def thread_process(self):

			while 1:
				time.sleep(5)


	class EMBoxViewer(QtGui.QWidget):
		"""This is a multi-paned view showing a single boxed out particle from a larger tomogram"""

		def __init__(self):
			QtGui.QWidget.__init__(self)
			self.setWindowTitle("Single Particle View")

			self.resize(300,300)

			self.gbl = QtWidgets.QGridLayout(self)
			self.xyview = EMImage2DWidget()
			self.gbl.addWidget(self.xyview,0,1)

			self.xzview = EMImage2DWidget()
			self.gbl.addWidget(self.xzview,1,1)

			self.zyview = EMImage2DWidget()
			self.gbl.addWidget(self.zyview,0,0)
			self.data = None


			# This puts an isosurface view in the lower left corner, but was causing a lot of segfaults, so switching to 2-D slices for now
			#self.d3view = EMScene3D()
			#self.d3viewdata = EMDataItem3D(test_image_3d(3), transform=Transform())
			#isosurface = EMIsosurface(self.d3viewdata, transform=Transform())
			#self.d3view.insertNewNode('', self.d3viewdata, parentnode=self.d3view)
			#self.d3view.insertNewNode("Iso", isosurface, parentnode=self.d3viewdata )

			self.d3view = EMImage2DWidget()
			self.gbl.addWidget(self.d3view,1,0)

			self.wfilt = ValSlider(rng=(0,50),label="Filter:",value=0.0)
			self.gbl.addWidget(self.wfilt,2,0,1,2)

			self.wfilt.valueChanged.connect(self.event_filter)

			self.gbl.setRowStretch(2,1)
			self.gbl.setRowStretch(0,5)
			self.gbl.setRowStretch(1,5)
			#QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("mousedown"),self.xy_down)
			#QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("mousedrag"),self.xy_drag)
			#QtCore.QObject.connect(self.xyview,QtCore.SIGNAL("mouseup")  ,self.xy_up  )

			#QtCore.QObject.connect(self.xzview,QtCore.SIGNAL("mousedown"),self.xz_down)
			#QtCore.QObject.connect(self.xzview,QtCore.SIGNAL("mousedrag"),self.xz_drag)
			#QtCore.QObject.connect(self.xzview,QtCore.SIGNAL("mouseup")  ,self.xz_up  )

			#QtCore.QObject.connect(self.zyview,QtCore.SIGNAL("mousedown"),self.zy_down)
			#QtCore.QObject.connect(self.zyview,QtCore.SIGNAL("mousedrag"),self.zy_drag)
			#QtCore.QObject.connect(self.zyview,QtCore.SIGNAL("mouseup")  ,self.zy_up  )

	#		self.setSizeGripEnabled(True)

	#		if get_platform() == "Darwin": # because OpenGL widgets in Qt don't leave room in the bottom right hand corner for the resize tool
	#			self.status = QtGui.QStatusBar()
	#			self.gbl.addWidget(self.status,3,0,1,2)
	#			self.margin = 0

		def set_data(self,data):
			"""Sets the current volume to display"""

			self.data=data
			self.fdata=data

			self.update()
			self.show()

		def get_data(self):
			return self.data

		def update(self):
			if self.data==None:
				self.xyview.set_data(None)
				self.xzview.set_data(None)
				self.zyview.set_data(None)

				#self.d3viewdata.setData(test_image_3d(3))
				#self.d3view.updateSG()
				self.d3view.set_data(test_image_3d(3))

				return

			if self.wfilt.getValue()>4 :
				self.fdata=self.data.process("filter.lowpass.gauss",{"cutoff_freq":old_div(1.0,self.wfilt.getValue()),"apix":self.data['apix_x']}) #JESUS

			xyd=self.fdata.process("misc.directional_sum",{"axis":"z"})
			xzd=self.fdata.process("misc.directional_sum",{"axis":"y"})
			zyd=self.fdata.process("misc.directional_sum",{"axis":"x"})

			self.xyview.set_data(xyd)
			self.xzview.set_data(xzd)
			self.zyview.set_data(zyd)

			#self.d3viewdata.setData(self.fdata)
			#self.d3view.updateSG()
			self.d3view.set_data(self.fdata)


		def event_filter(self,value):
			self.update()

		def closeEvent(self, event):
			self.d3view.close()


	class EMTomoBoxer(QtGui.QMainWindow):
		"""This class represents the EMTomoBoxer application instance.  """
		module_closed = QtCore.pyqtSignal()

		def __init__(self,application,data=None,datafile=None,yshort=False,apix=0.0,boxsize=32,shrink=1,contrast=None,center=None,mod=False,normalize=False):
			QtGui.QWidget.__init__(self)

			self.app=weakref.ref(application)
			self.yshort=yshort
			self.apix=apix

			self.shrink=shrink
			self.contrast=contrast
			self.mod=mod
			self.center=center
			self.normalize=normalize
			self.setWindowTitle("Main Window (e2spt_boxer.py)")

	#		self.setWindowTitle("e2spt_boxer.py")

			# Menu Bar
			self.mfile=self.menuBar().addMenu("File")
			self.mfile_open=self.mfile.addAction("Open")
			self.mfile_read_boxloc=self.mfile.addAction("Read Box Coord")
			self.mfile_save_boxloc=self.mfile.addAction("Save Box Coord")
			self.mfile_save_boxes=self.mfile.addAction("Save Boxed Data")
			self.mfile_save_boxes_stack=self.mfile.addAction("Save Boxes as Stack")
			self.mfile_quit=self.mfile.addAction("Quit")

			self.mwin=self.menuBar().addMenu("Window")
			self.mwin_boxes=self.mwin.addAction("Particles")
			self.mwin_single=self.mwin.addAction("Single Particle")
			self.mwin_average=self.mwin.addAction("Averaging")


			self.setCentralWidget(QtGui.QWidget())
			self.gbl = QtWidgets.QGridLayout(self.centralWidget())

			# relative stretch factors
			self.gbl.setColumnStretch(0,1)
			self.gbl.setColumnStretch(1,4)
			self.gbl.setColumnStretch(2,0)
			self.gbl.setRowStretch(1,1)
			self.gbl.setRowStretch(0,4)

			# 3 orthogonal restricted projection views
			self.xyview = EMImage2DWidget()
			self.gbl.addWidget(self.xyview,0,1)

			self.xzview = EMImage2DWidget()
			self.gbl.addWidget(self.xzview,1,1)

			self.zyview = EMImage2DWidget()
			self.gbl.addWidget(self.zyview,0,0)

			# Select Z for xy view
			self.wdepth = QtGui.QSlider()
			self.gbl.addWidget(self.wdepth,1,2)

			### Control panel area in upper left corner
			self.gbl2 = QtWidgets.QGridLayout()
			self.gbl.addLayout(self.gbl2,1,0)

			# box size
			self.wboxsize=ValBox(label="Box Size:",value=boxsize)
			self.gbl2.addWidget(self.wboxsize,1,0,1,2)
			self.oldboxsize=boxsize

			# max or mean
			self.wmaxmean=QtGui.QPushButton("MaxProj")
			self.wmaxmean.setCheckable(True)
			self.gbl2.addWidget(self.wmaxmean,2,0)

			# number slices
			self.wnlayers=QtGui.QSpinBox()
			self.wnlayers.setMinimum(1)
			self.wnlayers.setMaximum(256)
			self.wnlayers.setValue(1)
			self.gbl2.addWidget(self.wnlayers,2,1)

			# Local boxes in side view
			self.wlocalbox=QtWidgets.QCheckBox("Limit Side Boxes")
			self.gbl2.addWidget(self.wlocalbox,3,0)

			# scale factor
			self.wscale=ValSlider(rng=(.1,2),label="Sca:",value=1.0)
			self.gbl2.addWidget(self.wscale,4,0,1,2)

			# 2-D filters
			self.wfilt = ValSlider(rng=(0,150),label="Filt:",value=0.0)
			self.gbl2.addWidget(self.wfilt,5,0,1,2)

			self.curbox=-1
			self.boxes=[]						# array of box info, each is (x,y,z,...)
			self.helixboxes=[]					# array of helix box info. each is (xi, yi, zi, xf, yf, zf)
			self.boxesimgs=[]					# z projection of each box
			self.xydown=None
			self.firsthbclick = None

			# coordinate display
			self.wcoords=QtWidgets.QLabel("X: " + str(self.get_x()) + "\t\t" + "Y: " + str(self.get_y()) + "\t\t" + "Z: " + str(self.get_z()))
			self.gbl2.addWidget(self.wcoords, 0, 0, 1, 2)

			# file menu
			self.mfile_open.triggered[bool].connect(self.menu_file_open)
			self.mfile_read_boxloc.triggered[bool].connect(self.menu_file_read_boxloc)
			self.mfile_save_boxloc.triggered[bool].connect(self.menu_file_save_boxloc)
			self.mfile_save_boxes.triggered[bool].connect(self.menu_file_save_boxes)
			self.mfile_save_boxes_stack.triggered[bool].connect(self.menu_file_save_boxes_stack)
			self.mfile_quit.triggered[bool].connect(self.menu_file_quit)

			# window menu
			self.mwin_boxes.triggered[bool].connect(self.menu_win_boxes)
			self.mwin_single.triggered[bool].connect(self.menu_win_single)
	#		QtCore.QObject.connect(self.mwin_average,QtCore.SIGNAL("triggered(bool)")  ,self.menu_win_average  )

			# all other widgets
			self.wdepth.valueChanged[int].connect(self.event_depth)
			self.wnlayers.valueChanged[int].connect(self.event_nlayers)
			self.wboxsize.valueChanged.connect(self.event_boxsize)
			self.wmaxmean.clicked[bool].connect(self.event_projmode)
			self.wscale.valueChanged.connect(self.event_scale)
			self.wfilt.valueChanged.connect(self.event_filter)
			self.wlocalbox.stateChanged[int].connect(self.event_localbox)

			self.xyview.mousedown.connect(self.xy_down)
			self.xyview.mousedrag.connect(self.xy_drag)
			self.xyview.mouseup.connect(self.xy_up)
			self.xyview.mousewheel.connect(self.xy_wheel)
			self.xyview.signal_set_scale.connect(self.xy_scale)
			self.xyview.origin_update.connect(self.xy_origin)

			self.xzview.mousedown.connect(self.xz_down)
			self.xzview.mousedrag.connect(self.xz_drag)
			self.xzview.mouseup.connect(self.xz_up)
			self.xzview.signal_set_scale.connect(self.xz_scale)
			self.xzview.origin_update.connect(self.xz_origin)

			self.zyview.mousedown.connect(self.zy_down)
			self.zyview.mousedrag.connect(self.zy_drag)
			self.zyview.mouseup.connect(self.zy_up)
			self.zyview.signal_set_scale.connect(self.zy_scale)
			self.zyview.origin_update.connect(self.zy_origin)

			if datafile!=None:
				print("\nIn ETomoBoxer, datafile is", datafile)
				self.set_datafile(datafile)		# This triggers a lot of things to happen, so we do it last

			if data!=None:
				self.set_data(data)

			# Boxviewer subwidget (details of a single box)
			self.boxviewer=EMBoxViewer()
			#self.app().attach_child(self.boxviewer)

			# Boxes Viewer (z projections of all boxes)
			self.boxesviewer=EMImageMXWidget()
			#self.app().attach_child(self.boxesviewer)
			self.boxesviewer.show()
			self.boxesviewer.set_mouse_mode("App")
			self.boxesviewer.setWindowTitle("Particle List")

			# Average viewer shows results of background tomographic processing
	#		self.averageviewer=EMAverageViewer(self)
			#self.averageviewer.show()

			self.boxesviewer.mx_image_selected.connect(self.img_selected)
			self.e = None

		def menu_win_boxes(self) : self.boxesviewer.show()
		def menu_win_single(self) : self.boxviewer.show()
	#	def menu_win_average(self) : self.averageviewer.show()

		def set_datafile(self,datafile):
			print("\nIn set_datafile, received datafile", datafile)
			if datafile==None :
				self.datafile=None
				self.data=None
				self.xyview.set_data(None)
				self.xzview.set_data(None)
				self.zyview.set_data(None)
				return

			self.data=None
			self.datafile=datafile

			print("\nDatafile set, see!", self.datafile, type(self.datafile))

			imgh=EMData(datafile,0,1)

			if self.yshort:
				self.datasize=(imgh["nx"],imgh["nz"],imgh["ny"])
			else:
				self.datasize=(imgh["nx"],imgh["ny"],imgh["nz"])

			self.wdepth.setRange(0,self.datasize[2]-1)
			self.boxes=[]
			self.curbox=-1

			self.wdepth.setValue(old_div(self.datasize[2],2))
			self.update_all()

		def set_data(self,data):
			if data==None :
				self.datafile=None
				self.data=None
				self.xyview.set_data(None)
				self.xzview.set_data(None)
				self.zyview.set_data(None)
				return

			self.data=data
			self.datafile=None

			if self.yshort:
				self.datasize=(data["nx"],data["nz"],data["ny"])
			else:
				self.datasize=(data["nx"],data["ny"],data["nz"])

			self.wdepth.setRange(0,self.datasize[2]-1)
			self.boxes=[]
			self.curbox=-1

			self.wdepth.setValue(old_div(self.datasize[2],2))
			self.update_all()

		def get_cube(self,x,y,z):
			"""Returns a box-sized cube at the given center location"""
			bs=self.boxsize()

			if self.yshort:
				if self.data!=None:
					r=self.data.get_clip(Region(x-old_div(bs,2),z-old_div(bs,2),y-old_div(bs,2),bs,bs,bs))
					if options.normproc:
						r.process_inplace(options.normproc)
					r.process_inplace("xform",{"transform":Transform({"type":"eman","alt":90.0})})
					r.process_inplace("xform.mirror",{"axis":"z"})
				elif self.datafile!=None:
					r=EMData(self.datafile,0,0,Region(x-old_div(bs,2),z-old_div(bs,2),y-old_div(bs,2),bs,bs,bs))
					if options.normproc:
						r.process_inplace(options.normproc)
					r.process_inplace("xform",{"transform":Transform({"type":"eman","alt":90.0})})
					r.process_inplace("xform.mirror",{"axis":"z"})
				else: return None

			else :
				if self.data!=None:
					r=self.data.get_clip(Region(x-old_div(bs,2),y-old_div(bs,2),z-old_div(bs,2),bs,bs,bs))
				elif self.datafile!=None:
					r=EMData(self.datafile,0,0,Region(x-old_div(bs,2),y-old_div(bs,2),z-old_div(bs,2),bs,bs,bs))
				else: return None

			if self.apix!=0 :
				r["apix_x"]=self.apix
				r["apix_y"]=self.apix
				r["apix_z"]=self.apix

			if options.normproc:
				r.process_inplace(options.normproc)
			return r

		def get_slice(self,n,xyz):
			"""Reads a slice either from a file or the preloaded memory array.
			xyz is the axis along which 'n' runs, 0=x (yz), 1=y (xz), 2=z (xy)"""
			if self.yshort:
				if self.data!=None :
					if xyz==0:
						r=self.data.get_clip(Region(n,0,0,1,self.datasize[2],self.datasize[1]))
						r.set_size(self.datasize[2],self.datasize[1],1)
					elif xyz==2:
						r=self.data.get_clip(Region(0,n,0,self.datasize[0],1,self.datasize[1]))
						r.set_size(self.datasize[0],self.datasize[1],1)
					else:
						r=self.data.get_clip(Region(0,0,n,self.datasize[0],self.datasize[2],1))

				elif self.datafile!=None:
					if xyz==0:
						r=EMData()
						r.read_image(self.datafile,0,0,Region(n,0,0,1,self.datasize[2],self.datasize[1]))
						r.set_size(self.datasize[2],self.datasize[1],1)

					elif xyz==2:
						r=EMData()
						r.read_image(self.datafile,0,0,Region(0,n,0,self.datasize[0],1,self.datasize[1]))
						r.set_size(self.datasize[0],self.datasize[1],1)
					else:
						r=EMData()
						r.read_image(self.datafile,0,0,Region(0,0,n,self.datasize[0],self.datasize[2],1))
				else:
					return None

			else :
				if self.data!=None :
					if xyz==0:
						r=self.data.get_clip(Region(n,0,0,1,self.datasize[1],self.datasize[2]))
						r.set_size(self.datasize[1],self.datasize[2],1)
					elif xyz==1:
						r=self.data.get_clip(Region(0,n,0,self.datasize[0],1,self.datasize[2]))
						r.set_size(self.datasize[0],self.datasize[2],1)
					else:
						r=self.data.get_clip(Region(0,0,n,self.datasize[0],self.datasize[1],1))

				elif self.datafile!=None:
					if xyz==0:
						r=EMData()
						r.read_image(self.datafile,0,0,Region(n,0,0,1,self.datasize[1],self.datasize[2]))
						r.set_size(self.datasize[1],self.datasize[2],1)
					elif xyz==1:
						r=EMData()
						r.read_image(self.datafile,0,0,Region(0,n,0,self.datasize[0],1,self.datasize[2]))
						r.set_size(self.datasize[0],self.datasize[2],1)
					else:
						r=EMData()
						r.read_image(self.datafile,0,0,Region(0,0,n,self.datasize[0],self.datasize[1],1))

				else :
					return None

			if self.apix!=0 :
				r["apix_x"]=self.apix
				r["apix_y"]=self.apix
				r["apix_z"]=self.apix
			return r

		def event_boxsize(self):
			if self.boxsize()==self.oldboxsize:
				return
			self.oldboxsize=self.boxsize()

			cb=self.curbox
			for i in range(len(self.boxes)):
				self.update_box(i)
			self.update_box(cb)

		def event_projmode(self,state):
			"""Projection mode can be simple average (state=False) or maximum projection (state=True)"""
			self.update_all()

		def event_scale(self,newscale):
			self.xyview.set_scale(newscale)
			self.xzview.set_scale(newscale)
			self.zyview.set_scale(newscale)

		def event_depth(self):
			self.update_xy()

		def event_nlayers(self):
			self.update_all()

		def event_filter(self):
			self.update_all()

		def event_localbox(self,tog):
			self.update_sides()

		def boxsize(self):
			return int(self.wboxsize.getValue())

		def nlayers(self):
			return int(self.wnlayers.value())

		def depth(self):
			return int(self.wdepth.value())

		def scale(self):
			return self.wscale.getValue()

		def get_x(self):
			return self.get_coord(0)

		def get_y(self):
			return self.get_coord(1)

		def get_z(self):
			return self.depth()

		def get_coord(self, coord_index):
			if len(self.boxes) > 1:
				if self.curbox:
					return self.boxes[self.curbox][coord_index]
				else:
					return self.boxes[-1][coord_index]
			else:
				return 0


		def menu_file_open(self,tog):
			QtGui.QMessageBox.warning(None,"Error","Sorry, in the current version, you must provide a file to open on the command-line.")

		def load_box_yshort(self, boxcoords):
			if options.yshort:
				return [boxcoords[0], boxcoords[2], boxcoords[1]]
			else:
				return boxcoords

		def menu_file_read_boxloc(self):
			fsp=str(QtWidgets.QFileDialog.getOpenFileName(self, "Select output text file"))

			f=open(fsp,"r")
			if options.helixboxer:
				for b in f:
					b2=[old_div(int(float(i)),self.shrink) for i in b.split()[:6]]
					self.boxes.append(self.load_box_yshort(b2[3:6]))
					self.update_box(len(self.boxes)-1)
					self.helixboxes.append(b2)
					self.update_helixbox(len(self.helixboxes)-1)
					self.boxes.append(self.load_box_yshort(b2[0:3]))
					self.update_box(len(self.boxes)-1)
			else:
				for b in f:
					b2=[old_div(int(float(i)),self.shrink) for i in b.split()[:3]]
					self.boxes.append(b2)
					self.update_box(len(self.boxes)-1)
			f.close()

		def menu_file_save_boxloc(self):
			shrinkf=self.shrink 								#jesus

			fsp=str(QtWidgets.QFileDialog.getSaveFileName(self, "Select output text file"))

			out=open(fsp,"w")
			if options.helixboxer:
				for b in self.helixboxes:
					out.write("%d\t%d\t%d\t%d\t%d\t%d\n"%(b[0]*shrinkf,b[1]*shrinkf,b[2]*shrinkf,b[3]*shrinkf,b[4]*shrinkf,b[5]*shrinkf))
			else:
				for b in self.boxes:
					out.write("%d\t%d\t%d\n"%(b[0]*shrinkf,b[1]*shrinkf,b[2]*shrinkf))
			out.close()

		def menu_file_save_boxes(self):
			fsp=os.path.basename(str(QtWidgets.QFileDialog.getSaveFileName(self, "Select output file (numbers added)")))
			if ".hdf" not in fsp[-4:]:
				fsp += '.hdf'
			
			fspprjs=fsp.replace('.hdf','_prjs.hdf')
			prj=EMData() #Dummy

			progress = QtGui.QProgressDialog("Saving", "Abort", 0, len(self.boxes),None)
			if options.helixboxer:
				for i,b in enumerate(self.helixboxes):
					img = self.extract_subtomo_box(self.get_extended_a_vector(b), cshrink=self.shrink)

					#img['origin_x'] = 0
					#img['origin_y'] = 0
					#img['origin_z'] = 0

					if self.normalize:
						img.process_inplace(normalize)
					#img=img.process('normalize.edgemean')

					#if fsp[:4].lower()=="bdb:":
					#	img.write_image(os.path.join(options.path,"%s_%03d"%(fsp,i)),0)

					#if "." in fsp:
					img.write_image(os.path.join(options.path,"%s_%03d.%s"%(fsp.rsplit(".hdf",1)[0],i,'hdf')))
					
					#else:
					#	QtGui.QMessageBox.warning(None,"Error","Please provide a valid image file extension. The numerical sequence will be inserted before the extension.")
					#	return

					progress.setValue(i+1)
					if progress.wasCanceled() : break
			else:
				for i,b in enumerate(self.boxes):
					#img=self.get_cube(b[0],b[1],b[2])
					bs=self.boxsize()
					shrinkf=self.shrink
					if shrinkf >1:
						bs=bs*shrinkf

					contrast=self.contrast
					center=self.center

					#if self.yshort:
					#	ret = unbinned_extractor(options,bs,b[0],b[2],b[1],shrinkf,contrast,center,args[0])
					#	img = ret[0]
					#	prj = ret[1]
					#else:
					#ret = unbinned_extractor(options,bs,b[0],b[1],b[2],shrinkf,contrast,center,args[0])
					options.boxsize = bs
					ret = unbinned_extractor(options,b[0],b[1],b[2],args[0],i)
					img = ret[0]
					prj = ret[1]

					print("extracting particle %d from coordinates x=%d, y=%d, z=%d" % (i, b[0],b[1],b[2]))

					#if "." in fsp:
						#img.write_image(os.path.join(options.path,"%s_%03d.%s"%(fsp.rsplit(".",1)[0],i,fsp.rsplit(".",1)[1])))
					img.write_image(os.path.join(options.path,"%s_%03d.%s"%(fsp.rsplit(".",1)[0],i,'hdf')))
					prj.write_image(fspprjs,-1)

					#else:
					#	QtGui.QMessageBox.warning(None,"Error","Please provide a valid image file extension. The numerical sequence will be inserted before the extension.")
					#	return

					progress.setValue(i+1)
					if progress.wasCanceled() : break

		def menu_file_save_boxes_stack(self):

			fsp=os.path.join(options.path,os.path.basename(str(QtWidgets.QFileDialog.getSaveFileName(self, "Select output file (.hdf supported only)"))))
			#if fsp[:4].lower()!="bdb:" and fsp[-4:].lower()!=".hdf" :


			if fsp[-4:].lower()!=".hdf" :
				#QtGui.QMessageBox.warning(None,"Error","3-D stacks supported only for .hdf files")
				#return
				fsp+='.hdf'

			fspprjs=fsp.replace('.hdf','_prjs.hdf')
			prj=EMData() #Dummy

			progress = QtGui.QProgressDialog("Saving", "Abort", 0, len(self.boxes),None)
			if options.helixboxer:
				for i,b in enumerate(self.helixboxes):
					img = self.extract_subtomo_box(self.get_extended_a_vector(b), cshrink=self.shrink)

					#img['origin_x'] = 0
					#img['origin_y'] = 0
					#img['origin_z'] = 0
					if self.normalize:
						e.process_inplace(normalize)
					#img=img.process('normalize.edgemean')
					print("extracting particle %d" % (i))

					img.write_image(fsp,i)

					progress.setValue(i+1)
					if progress.wasCanceled():
						break
			else:
				for i,b in enumerate(self.boxes):
					#img=self.get_cube(b[0],b[1],b[2])
					bs=self.boxsize()
					shrinkf=self.shrink
					if shrinkf >1:
						bs=bs*shrinkf

					contrast=self.contrast
					center=self.center

					#if self.yshort:
					#	ret = unbinned_extractor(options,bs,b[0],b[2],b[1],shrinkf,contrast,center,args[0])
					#	img = ret[0]
					#	prj = ret[1]
					#else:
					#ret = unbinned_extractor(options,bs,b[0],b[1],b[2],shrinkf,contrast,center,args[0])
					#print "before extraction shrinkf is", shrinkf
					print("extracting particle %d from coordinates x=%d, y=%d, z=%d" % (i, b[0],b[1],b[2]))
					options.boxsize=bs
					ret = unbinned_extractor(options,b[0],b[1],b[2],args[0],i)
					img = ret[0]
					prj = ret[1]
					#print "returned image is of size",img['nx']

					#img['origin_x'] = 0
					#img['origin_y'] = 0
					#img['origin_z'] = 0
					if self.normalize:
						img.process_inplace(normalize)
					#img=img.process('normalize.edgemean')

					img.write_image(fsp,i)
					prj.write_image(fspprjs,-1)

					progress.setValue(i+1)
					if progress.wasCanceled():
						break


		def menu_file_quit(self):
			self.close()

		def transform_coords(self, point, xform):
			xvec = xform.get_matrix()
			return [xvec[0]*point[0] + xvec[4]*point[1] + xvec[8]*point[2] + xvec[3], xvec[1]*point[0] + xvec[5]*point[1] + xvec[9]*point[2] + xvec[7], xvec[2]*point[0] + xvec[6]*point[1] + xvec[10]*point[2] + xvec[11]]

		def extract_subtomo_box(self, helixbox, cshrink=1, tomogram=argv[1]):
			""" Retruns an extracted subtomogram box"""
			# Only scale the helix boxer values transiently
			x1 = round(helixbox[0]*cshrink)
			y1 = round(helixbox[1]*cshrink)
			z1 = round(helixbox[2]*cshrink)
			x2 = round(helixbox[3]*cshrink)
			y2 = round(helixbox[4]*cshrink)
			z2 = round(helixbox[5]*cshrink)

			bs=old_div(self.boxsize(),2)
			# Get the extended vector based on boxsize
			a = Vec3f((x2-x1), (y2-y1), (z2-z1))	# Find the a, the long vector
			tcs = self.get_box_coord_system([x1,y1,z1,x2,y2,z2])							# Get the local coord system
			# Get the new coord system
			# First extract a subtomo gram bounding region from the tomogram so we do have to read the whole bloody thing in!
			rv = [self.transform_coords([0, -bs, -bs], tcs), self.transform_coords([0, bs, bs], tcs), self.transform_coords([0, bs, -bs], tcs), self.transform_coords([0, -bs, bs], tcs), self.transform_coords([a.length(), -bs, -bs], tcs), self.transform_coords([a.length(), bs, bs], tcs), self.transform_coords([a.length(), bs, -bs], tcs), self.transform_coords([a.length(), -bs, bs], tcs)]
			rvmin = [int(min([i[0] for i in rv])), int(min([i[1] for i in rv])), int(min([i[2] for i in rv]))]	# Min bounding box extension
			rvmax = [int(max([i[0] for i in rv])), int(max([i[1] for i in rv])), int(max([i[2] for i in rv]))]	# Max bounding box extension
			r = Region(rvmin[0],rvmin[1],rvmin[2],rvmax[0]-rvmin[0],rvmax[1]-rvmin[1],rvmax[2]-rvmin[2])		# Extract the region
			e = EMData()
			e.read_image(tomogram,0,False,r)
			e.set_attr("source_path", tomogram)
			e["ptcl_source_image"]=tomogram
			e["ptcl_source_coord"]=(old_div((rvmin[0]+rvmax[0]),2),old_div((rvmin[1]+rvmax[1]),2),old_div((rvmin[2]+rvmax[2]),2))
			# Next adjust the transform matrix to move it to the origin
			origin = self.transform_coords([0,0,0], tcs)
			tcs.set_trans(origin[0] - rvmin[0], origin[1] - rvmin[1], origin[2] - rvmin[2])

			return e.extract_box(tcs, Region(0, -bs, -bs, a.length(), bs, bs))

		def get_averager(self):
			"""returns an averager of the appropriate type for generating projection views"""
			if self.wmaxmean.isChecked() : return Averagers.get("minmax",{"max":1})

			return Averagers.get("mean")

		def update_sides(self):
			"""updates xz and yz views due to a new center location"""

			#print "\n\n\n\n\nIn update sides, self.datafile is", self.datafile
			#print "\n\n\n\n"

			if self.datafile==None and self.data==None:
				return

			if self.curbox==-1 :
				x=old_div(self.datasize[0],2)
				y=old_div(self.datasize[1],2)
				z=0
			else:
				x,y,z=self.boxes[self.curbox][:3]

			self.cury=y
			self.curx=x
			bs=self.boxsize()

			# update shape display
			if self.wlocalbox.isChecked():
				xzs=self.xzview.get_shapes()
				if options.helixboxer:
					for i in range(0, len(self.boxes), 2):
						if i + 1 >= len(self.boxes):
							break
						#if abs(self.boxes[i][1] - zc) < bs/2 or abs(self.boxes[i+1][2] - zc) < bs/2:
						if (self.boxes[i][1]<self.cury+old_div(bs,2) and self.boxes[i][1]>self.cury-old_div(bs,2)) or (self.boxes[i+1][1]<self.cury+old_div(bs,2) and self.boxes[i+1][1]>self.cury-old_div(bs,2)):
							xzs[i][0]="rect"
							xzs[str(old_div(i,2))+"helix"][0]="line"
							xzs[i+1][0]="rect"
						else:
							xzs[i][0]="hidden"
							xzs[str(old_div(i,2))+"helix"][0]="hidden"
							xzs[i+1][0]="hidden"
				else:
						for i in range(len(self.boxes)):
							if self.boxes[i][1]<self.cury+old_div(bs,2) and self.boxes[i][1]>self.cury-old_div(bs,2):
								xzs[i][0]="rect"
							else:
								xzs[i][0]="hidden"

				zys=self.zyview.get_shapes()
				if options.helixboxer:
					for i in range(0, len(self.boxes), 2):
						if i + 1 >= len(self.boxes):
							break
						if (self.boxes[i][0]<self.curx+old_div(bs,2) and self.boxes[i][0]>self.curx-old_div(bs,2)) or (self.boxes[i+1][0]<self.curx+old_div(bs,2) and self.boxes[i+1][0]>self.curx-old_div(bs,2)):
							zys[i][0]="rect"
							zys[str(old_div(i,2))+"helix"][0]="line"
							zys[i+1][0]="rect"
						else:
							zys[i][0]="hidden"
							zys[str(old_div(i,2))+"helix"][0]="hidden"
							zys[i+1][0]="hidden"
				else:
					for i in range(len(self.boxes)):
						if self.boxes[i][0]<self.curx+old_div(bs,2) and self.boxes[i][0]>self.curx-old_div(bs,2):
							zys[i][0]="rect"
						else:
							zys[i][0]="hidden"
			else :
				xzs=self.xzview.get_shapes()
				zys=self.zyview.get_shapes()
				if options.helixboxer:
					for i in range(0, len(self.boxes), 2):
						if i + 1 >= len(self.boxes):
							break
						xzs[i][0]="rect"
						xzs[str(old_div(i,2))+"helix"][0]="line"
						xzs[i+1][0]="rect"
						zys[i][0]="rect"
						zys[str(old_div(i,2))+"helix"][0]="line"
						zys[i+1][0]="rect"
				else:
					for i in range(len(self.boxes)):
						xzs[i][0]="rect"
						zys[i][0]="rect"

			self.xzview.shapechange=1
			self.zyview.shapechange=1

			# yz
			avgr=self.get_averager()

			for x in range(x-old_div(self.nlayers(),2),x+old_div((self.nlayers()+1),2)):
				slc=self.get_slice(x,0)
				avgr.add_image(slc)

			av=avgr.finish()
			if not self.yshort:
				av.process_inplace("xform.transpose")

			if self.wfilt.getValue()!=0.0:
				av.process_inplace("filter.lowpass.gauss",{"cutoff_freq":old_div(1.0,self.wfilt.getValue()),"apix":self.apix})

			self.zyview.set_data(av)

			# xz
			avgr=self.get_averager()

			for y in range(y-old_div(self.nlayers(),2),y+old_div((self.nlayers()+1),2)):
				slc=self.get_slice(y,1)
				avgr.add_image(slc)

			av=avgr.finish()
			if self.wfilt.getValue()!=0.0:
				av.process_inplace("filter.lowpass.gauss",{"cutoff_freq":old_div(1.0,self.wfilt.getValue()),"apix":self.apix})

			self.xzview.set_data(av)


		def update_xy(self):
			"""updates xy view due to a new slice range"""

			#print "\n\n\n\n\nIn update_xy, self.datafile is", self.datafile
			#print "\n\n\n\n"

			if self.datafile==None and self.data==None:
				return



			# Boxes should also be limited by default in the XY view
			if len(self.boxes) > 0:
				zc=self.wdepth.value()
				#print "The current depth is", self.wdepth.value()
				bs=self.boxsize()
				xys=self.xyview.get_shapes()
				if options.helixboxer:
					for i in range(0, len(self.boxes), 2):
						if i + 1 >= len(self.boxes):
							break
						if abs(self.boxes[i][2] - zc) < old_div(bs,2) or abs(self.boxes[i+1][2] - zc) < old_div(bs,2):
							xys[i][0]="rect"
							xys[str(old_div(i,2))+"helix"][0]="line"
							xys[i+1][0]="rect"
						else:
							xys[i][0]="hidden"
							xys[str(old_div(i,2))+"helix"][0]="hidden"
							xys[i+1][0]="hidden"
				else:
					for i in range(len(self.boxes)):
						#print "the z coord of box %d is %d" %(i,self.boxes[i][2])
						#print "therefore the criteria to determine whether to display it is", abs(self.boxes[i][2] - zc)
						if abs(self.boxes[i][2] - zc) < old_div(bs,2):
							#print "Which is less than half the box thus it survives"
							xys[i][0]="rect"
						else :
							xys[i][0]="hidden"
							#print "Which is more than half the box and thus it dies"

				self.xyview.shapechange=1

			if self.wmaxmean.isChecked():
				avgr=Averagers.get("minmax",{"max":1})

			else:
				avgr=Averagers.get("mean")

			slc=EMData()
			for z in range(self.wdepth.value()-old_div(self.nlayers(),2),self.wdepth.value()+old_div((self.nlayers()+1),2)):
				slc=self.get_slice(z,2)
				avgr.add_image(slc)

			av=avgr.finish()

			#print "\n\nIn update xy, av and type are", av, type(av)

			if self.wfilt.getValue()!=0.0:

				av.process_inplace("filter.lowpass.gauss",{"cutoff_freq":old_div(1.0,self.wfilt.getValue()),"apix":self.apix})
			self.xyview.set_data(av)

		def update_all(self):
			"""redisplay of all widgets"""

			#print "\n\n\n\n\nIn update all, self.datafile is", self.datafile
			#print "\n\n\n\n"
			if self.datafile==None and self.data==None:
				return



			self.update_xy()
			self.update_sides()

			#self.xyview.update()
			#self.xzview.update()
			#self.zyview.update()

		def update_coords(self):
			self.wcoords.setText("X: " + str(self.get_x()) + "\t\t" + "Y: " + str(self.get_y()) + "\t\t" + "Z: " + str(self.get_z()))

		def inside_box(self,n,x=-1,y=-1,z=-1):
			"""Checks to see if a point in image coordinates is inside box number n. If any value is negative, it will not be checked."""
			box=self.boxes[n]
			if x>=0 and (x<box[0]-old_div(self.boxsize(),2) or x>box[0]+old_div(self.boxsize(),2)) : return False
			if y>=0 and (y<box[1]-old_div(self.boxsize(),2) or y>box[1]+old_div(self.boxsize(),2)) : return False
			if z>=0 and (z<box[2]-old_div(self.boxsize(),2) or z>box[2]+old_div(self.boxsize(),2)) : return False
			return True

		def do_deletion(self, n, delimgs=True):
			""" Helper for del_box"""
			if n==len(self.boxes)-1 :
				self.boxes.pop()
				if delimgs:
					self.boxesimgs.pop()
					self.boxesviewer.set_data(self.boxesimgs)
					self.boxesviewer.update()
				self.xyview.del_shape(n)
				self.xzview.del_shape(n)
				self.zyview.del_shape(n)
				self.curbox=-1
				self.xyview.update()
				self.xzview.update()
				self.zyview.update()
				#self.update()
			else :
				a=self.boxes.pop()
				self.boxes[n]=a
				if delimgs:
					a=self.boxesimgs.pop()
					self.boxesimgs[n]=a
					self.boxesviewer.set_data(self.boxesimgs)
					self.boxesviewer.set_selected([],True)
					self.boxesviewer.update()
				self.xyview.del_shape(len(self.boxes))
				self.xzview.del_shape(len(self.boxes))
				self.zyview.del_shape(len(self.boxes))
				self.update_box(n,True)
	#			self.update()

		def do_helix_deletion(self, n):
			if n==len(self.helixboxes)-1 :
				self.helixboxes.pop()
				self.xyview.del_shape(str(n)+"helix")
				self.xzview.del_shape(str(n)+"helix")
				self.zyview.del_shape(str(n)+"helix")
			else:
				a=self.helixboxes.pop()
				self.helixboxes[n]=a
				self.xyview.del_shape(str(len(self.helixboxes))+"helix")
				self.xzview.del_shape(str(len(self.helixboxes))+"helix")
				self.zyview.del_shape(str(len(self.helixboxes))+"helix")
				self.update_helixbox(n)

		def del_box(self,n):
			"""Delete an existing box by replacing the deleted box with the last box. A bit funny, but otherwise
			update after deletion is REALLY slow."""
	#		print "del ",n
			if n<0 or n>=len(self.boxes): return

			if self.boxviewer.get_data(): self.boxviewer.set_data(None)
			self.curbox=-1
			if options.helixboxer:
				if n + 1 == len(self.boxes) and len(self.boxes) % 2 == 1: 	# Delete unpaired box
					self.do_deletion(n, delimgs=False)
				else:								# Delete box pairs
					if n % 2:
						self.do_helix_deletion(int(old_div(n,2)))
						self.do_deletion(n, delimgs=False)
						self.do_deletion(n-1, delimgs=False)
					else:
						self.do_helix_deletion(int(old_div(n,2)))
						self.do_deletion(n+1, delimgs=False)
						self.do_deletion(n, delimgs=False)
					return "DELHELIX"	# If we have deleted a pair do not reset the pair toggle/counter
			else:
				self.do_deletion(n)

		def compute_crossAB(self, a, b):
			c1 = a[1]*b[2] - a[2]*b[1]
			c2 = a[2]*b[0] - a[0]*b[2]
			c3 = a[0]*b[1] - a[1]*b[0]
			return Vec3f(c1,c2,c3)

		def compute_perpZ(self, a):
			# Z axis
			b1 = -a[1]
			b2 = a[0]
			b3 = 0
			return Vec3f(b1,b2,b3)

		def compute_perpY(self, a):
			# Y axis
			b1 = -a[2]
			b2 = 0
			b3 = a[0]

			return Vec3f(b1,b2,b3)

		def get_box_coord_system(self, helixbox):
			"""
			Compute the coordinate system for the box
			"""
			a = Vec3f((helixbox[0]-helixbox[3]), (helixbox[1]-helixbox[4]), (helixbox[2]-helixbox[5]))

			a.normalize()
			b = self.compute_perpZ(a)
			b.normalize()
			c = self.compute_crossAB(a, b)

			return Transform([a[0],a[1],a[2],helixbox[3],b[0],b[1],b[2],helixbox[4],c[0],c[1],c[2],helixbox[5]])


		def get_extended_a_vector(self, helixbox):
			"""
			Extend the A vector to the box ends
			"""
			a = Vec3f((helixbox[3]-helixbox[0]), (helixbox[4]-helixbox[1]), (helixbox[5]-helixbox[2]))
			a.normalize()
			bs = self.boxsize()
			return [(helixbox[0] - old_div(a[0]*bs,2)),(helixbox[1] - old_div(a[1]*bs,2)),(helixbox[2] - old_div(a[2]*bs,2)),(helixbox[3] + old_div(a[0]*bs,2)),(helixbox[4] + old_div(a[1]*bs,2)),(helixbox[5] + old_div(a[2]*bs,2))]


		def update_helixbox(self, n, quiet=False):
			"""
			Update a helix box
			"""
			if n > len(self.helixboxes)-1: return	# Some boxes may not be paired
			helixbox = self.get_extended_a_vector(self.helixboxes[n])

			key = str(n)+"helix"
			if options.yshort:
				self.xyview.add_shape(key,EMShape(("line",.2,.2,.8, helixbox[0], helixbox[2], helixbox[3], helixbox[5],2)))
				self.xzview.add_shape(key,EMShape(("line",.2,.2,.8, helixbox[0], helixbox[1], helixbox[3], helixbox[4],2)))
				self.zyview.add_shape(key,EMShape(("line",.2,.2,.8, helixbox[1], helixbox[2], helixbox[4], helixbox[5],2)))
			else:
				self.xyview.add_shape(key,EMShape(("line",.2,.2,.8, helixbox[0], helixbox[1], helixbox[3], helixbox[4],2)))
				self.xzview.add_shape(key,EMShape(("line",.2,.2,.8, helixbox[0], helixbox[2], helixbox[3], helixbox[5],2)))
				self.zyview.add_shape(key,EMShape(("line",.2,.2,.8, helixbox[2], helixbox[1], helixbox[5], helixbox[4],2)))
			self.xyview.update()
			self.xzview.update()
			self.zyview.update()

			if not quiet and options.helixboxer:
				hb = self.extract_subtomo_box(helixbox, cshrink=self.shrink)
				self.boxviewer.set_data(hb)

				proj=hb.process("misc.directional_sum",{"axis":"z"})
				try: self.boxesimgs[n]=proj
				except:
					for i in range(len(self.boxesimgs),n+1): self.boxesimgs.append(None)
					self.boxesimgs[n]=proj
				self.boxesviewer.set_data(self.boxesimgs)
				self.boxesviewer.update()

			if n!=self.curbox and options.helixboxer:
				self.boxesviewer.set_selected((n,),True)

		def update_box(self,n,quiet=False):
			"""After adjusting a box, call this"""
	#		print "upd ",n,quiet

			try:
				box=self.boxes[n]
			except IndexError:
				return
			bs2=old_div(self.boxsize(),2)

			#if self.curbox!=n :
				#self.xzview.scroll_to(None,box[2])
				#self.zyview.scroll_to(box[2],None)


			# Boxes may not extend outside the tomogram
			if box[0]<bs2 : box[0]=bs2
			if box[0]>self.datasize[0]-bs2 : box[0]=self.datasize[0]-bs2
			if box[1]<bs2 : box[1]=bs2
			if box[1]>self.datasize[1]-bs2 : box[1]=self.datasize[1]-bs2
			if box[2]<bs2 : box[2]=bs2
			if box[2]>self.datasize[2]-bs2 : box[2]=self.datasize[2]-bs2
	#		print self.boxes
			self.xyview.add_shape(n,EMShape(("rect",.2,.2,.8,box[0]-bs2,box[1]-bs2,box[0]+bs2,box[1]+bs2,2)))
			self.xyview.add_shape("xl",EMShape(("line",.8,.8,.1,0,box[1],self.datasize[0],box[1],1)))
			self.xyview.add_shape("yl",EMShape(("line",.8,.8,.1,box[0],0,box[0],self.datasize[1],1)))
			self.xzview.add_shape(n,EMShape(("rect",.2,.2,.8,box[0]-bs2,box[2]-bs2,box[0]+bs2,box[2]+bs2,2)))
			self.xzview.add_shape("xl",EMShape(("line",.8,.8,.1,0,box[2],self.datasize[0],box[2],1)))
			self.xzview.add_shape("zl",EMShape(("line",.8,.8,.1,box[0],0,box[0],self.datasize[2],1)))
			self.zyview.add_shape(n,EMShape(("rect",.2,.2,.8,box[2]-bs2,box[1]-bs2,box[2]+bs2,box[1]+bs2,2)))
			self.zyview.add_shape("yl",EMShape(("line",.8,.8,.1,box[2],0,box[2],self.datasize[1],1)))
			self.zyview.add_shape("zl",EMShape(("line",.8,.8,.1,0,box[1],self.datasize[2],box[1],1)))

			if self.depth()!=box[2]:
				self.wdepth.setValue(box[2])
			else:
				self.xyview.update()
			self.update_sides()

			# For speed, we turn off updates while dragging a box around. Quiet is set until the mouse-up
			if not quiet and not options.helixboxer:
				# Get the cube from the original data (normalized)
				cube=self.get_cube(*box)
				self.boxviewer.set_data(cube)

				# Make a z projection and store it in the list of all boxes
				proj=cube.process("misc.directional_sum",{"axis":"z"})
				try: self.boxesimgs[n]=proj
				except:
					for i in range(len(self.boxesimgs),n+1): self.boxesimgs.append(None)
					self.boxesimgs[n]=proj
				self.boxesviewer.set_data(self.boxesimgs)
				self.boxesviewer.update()

			if n!=self.curbox and not options.helixboxer:
				self.boxesviewer.set_selected((n,),True)

			self.curbox=n
			self.update_coords()



		def img_selected(self,event,lc):
	#		print "sel",lc[0]
			if event.modifiers()&Qt.ShiftModifier:
				self.del_box(lc[0])
			else:
				self.update_box(lc[0])
			if self.curbox>=0 :
				box=self.boxes[self.curbox]
				self.xyview.scroll_to(box[0],box[1])
				self.xzview.scroll_to(None,box[2])
				self.zyview.scroll_to(box[2],None)

		def add_helix_box(self, xf, yf, zf, xi, yi, zi):
			print(xf, yf, zf, xi, yi, zi)
			if options.yshort:
				self.helixboxes.append([xf, zf, yf, xi, zi, yi])
			else:
				self.helixboxes.append([xf, yf, zf, xi, yi, zi])

		def xy_down(self,event):
			x,y=self.xyview.scr_to_img((event.x(),event.y()))
			x,y=int(x),int(y)
			z=int(self.get_z())
			self.xydown=None
			if x<0 or y<0 : return		# no clicking outside the image (on 2 sides)

			for i in range(len(self.boxes)):
				if self.inside_box(i,x,y,z):
					if event.modifiers()&Qt.ShiftModifier:
						if self.del_box(i) != "DELHELIX": self.firsthbclick = None
					else:
						self.xydown=(i,x,y,self.boxes[i][0],self.boxes[i][1])
						if options.helixboxer: self.update_helixbox(int(old_div(i,2)))
						self.update_box(i)
					break
			else:
	#			if x>self.boxsize()/2 and x<self.datasize[0]-self.boxsize()/2 and y>self.boxsize()/2 and y<self.datasize[1]-self.boxsize()/2 and self.depth()>self.boxsize()/2 and self.depth()<self.datasize[2]-self.boxsize()/2 :
				if not event.modifiers()&Qt.ShiftModifier:
					###########
					if options.helixboxer:	# Only create a helixbox every 2 clicks
						if self.firsthbclick:
							self.add_helix_box(x, y, self.depth(), self.firsthbclick[0], self.firsthbclick[1], self.firsthbclick[2])
							self.firsthbclick = None
							self.update_helixbox(len(self.helixboxes)-1)
						else:
							self.firsthbclick = [x, y, self.depth()]
					###########
					self.boxes.append(([x,y,self.depth()]))
					self.xydown=(len(self.boxes)-1,x,y,x,y)		# box #, x down, y down, x box at down, y box at down
					self.update_box(self.xydown[0])

			if self.curbox>=0:
				box=self.boxes[self.curbox]
				self.xzview.scroll_to(None,box[2])
				self.zyview.scroll_to(box[2],None)

		def xy_drag(self,event):
			if self.xydown==None : return

			x,y=self.xyview.scr_to_img((event.x(),event.y()))
			x,y=int(x),int(y)

			dx=x-self.xydown[1]
			dy=y-self.xydown[2]
			if options.helixboxer:
				if len(self.boxes) % 2 == 0 or (self.xydown[0] != len(self.boxes)-1):	# Only update the helix boxer if it is paired, otherwise treat it as a regular box
					hb = self.helixboxes[int(old_div(self.xydown[0],2))]
					if self.xydown[0] % 2 == 0:
						hb[3] = dx+self.xydown[3]
						hb[4] = dy+self.xydown[4]
					else:
						hb[0] = dx+self.xydown[3]
						hb[1] = dy+self.xydown[4]
					self.update_helixbox(int(old_div(self.xydown[0],2)))
				else:
					self.firsthbclick[0] = x
					self.firsthbclick[1] = y

			self.boxes[self.xydown[0]][0]=dx+self.xydown[3]
			self.boxes[self.xydown[0]][1]=dy+self.xydown[4]
			self.update_box(self.curbox,True)

		def xy_up  (self,event):
			if self.xydown!=None: self.update_box(self.curbox)
			self.xydown=None

		def xy_wheel (self,event):
			if event.delta() > 0:
				#self.wdepth.setValue(self.wdepth.value()+4)
				self.wdepth.setValue(self.wdepth.value()+1) #jesus

			elif event.delta() < 0:
				#self.wdepth.setValue(self.wdepth.value()-4)
				self.wdepth.setValue(self.wdepth.value()-1) #jesus


		def xy_scale(self,news):
			"xy image view has been rescaled"
			self.wscale.setValue(news)
			#self.xzview.set_scale(news,True)
			#self.zyview.set_scale(news,True)

		def xy_origin(self,newor):
			"xy origin change"
			xzo=self.xzview.get_origin()
			self.xzview.set_origin(newor[0],xzo[1],True)

			zyo=self.zyview.get_origin()
			self.zyview.set_origin(zyo[0],newor[1],True)

		def xz_down(self,event):
			x,z=self.xzview.scr_to_img((event.x(),event.y()))
			x,z=int(x),int(z)
			y=int(self.get_y())
			self.xzdown=None
			if x<0 or z<0 : return		# no clicking outside the image (on 2 sides)

			for i in range(len(self.boxes)):
				if (not self.wlocalbox.isChecked() and self.inside_box(i,x,y,z)) or self.inside_box(i,x,self.cury,z) :
					if event.modifiers()&Qt.ShiftModifier:
						if self.del_box(i) != "DELHELIX": self.firsthbclick = None
					else :
						self.xzdown=(i,x,z,self.boxes[i][0],self.boxes[i][2])
						if options.helixboxer: self.update_helixbox(int(old_div(i,2)))
						self.update_box(i)
					break
			else:
				if not event.modifiers()&Qt.ShiftModifier:
					###########
					if options.helixboxer:	# Only create a helixbox every 2 clicks
						if self.firsthbclick:
							self.add_helix_box(x, self.cury, z, self.firsthbclick[0], self.firsthbclick[1], self.firsthbclick[2])
							self.firsthbclick = None
							self.update_helixbox(len(self.helixboxes)-1)
						else:
							self.firsthbclick = [x, self.cury, z]
					###########
					self.boxes.append(([x,self.cury,z]))
					self.xzdown=(len(self.boxes)-1,x,z,x,z)		# box #, x down, y down, x box at down, y box at down
					self.update_box(self.xzdown[0])

			if self.curbox>=0 :
				box=self.boxes[self.curbox]
				self.xyview.scroll_to(None,box[1])
				self.zyview.scroll_to(box[2],None)

		def xz_drag(self,event):
			if self.xzdown==None : return

			x,z=self.xzview.scr_to_img((event.x(),event.y()))
			x,z=int(x),int(z)

			dx=x-self.xzdown[1]
			dz=z-self.xzdown[2]
			if options.helixboxer:
				if len(self.boxes) % 2 == 0 or (self.xzdown[0] != len(self.boxes)-1):	# Only update the helix boxer if it is paired, otherwise treat it as a regular box
					hb = self.helixboxes[int(old_div(self.xzdown[0],2))]
					if self.xzdown[0] % 2 == 0:
						hb[3] = dx+self.xzdown[3]
						hb[5] = dz+self.xzdown[4]
					else:
						hb[0] = dx+self.xzdown[3]
						hb[2] = dz+self.xzdown[4]
					self.update_helixbox(int(old_div(self.xzdown[0],2)))
				else:
					self.firsthbclick[0] = x
					self.firsthbclick[2] = z

			self.boxes[self.xzdown[0]][0]=dx+self.xzdown[3]
			self.boxes[self.xzdown[0]][2]=dz+self.xzdown[4]
			self.update_box(self.curbox,True)

		def xz_up  (self,event):
			if self.xzdown!=None: self.update_box(self.curbox)
			self.xzdown=None

		def xz_scale(self,news):
			"xy image view has been rescaled"
			self.wscale.setValue(news)
			#self.xyview.set_scale(news,True)
			#self.zyview.set_scale(news,True)

		def xz_origin(self,newor):
			"xy origin change"
			xyo=self.xyview.get_origin()
			self.xyview.set_origin(newor[0],xyo[1],True)

			#zyo=self.zyview.get_origin()
			#self.zyview.set_origin(zyo[0],newor[1],True)


		def zy_down(self,event):
			z,y=self.zyview.scr_to_img((event.x(),event.y()))
			z,y=int(z),int(y)
			x=int(self.get_x())
			self.xydown=None
			if z<0 or y<0 : return		# no clicking outside the image (on 2 sides)

			for i in range(len(self.boxes)):
				if (not self.wlocalbox.isChecked() and self.inside_box(i,x,y,z)) or  self.inside_box(i,self.curx,y,z):
					if event.modifiers()&Qt.ShiftModifier:
						if self.del_box(i) != "DELHELIX": self.firsthbclick = None
					else :
						self.zydown=(i,z,y,self.boxes[i][2],self.boxes[i][1])
						if options.helixboxer: self.update_helixbox(int(old_div(i,2)))
						self.update_box(i)
					break
			else:
				if not event.modifiers()&Qt.ShiftModifier:
					###########
					if options.helixboxer:	# Only create a helixbox every 2 clicks
						if self.firsthbclick:
							self.add_helix_box(self.curx, y, z, self.firsthbclick[0], self.firsthbclick[1], self.firsthbclick[2])
							self.firsthbclick = None
							self.update_helixbox(len(self.helixboxes)-1)
						else:
							self.firsthbclick = [self.curx, y, z]
					###########
					self.boxes.append(([self.curx,y,z]))
					self.zydown=(len(self.boxes)-1,z,y,z,y)		# box #, x down, y down, x box at down, y box at down
					self.update_box(self.zydown[0])

			if self.curbox>=0 :
				box=self.boxes[self.curbox]
				self.xyview.scroll_to(box[0],None)
				self.xzview.scroll_to(None,box[2])

		def zy_drag(self,event):
			if self.zydown==None : return

			z,y=self.zyview.scr_to_img((event.x(),event.y()))
			z,y=int(z),int(y)

			dz=z-self.zydown[1]
			dy=y-self.zydown[2]
			if options.helixboxer:
				if len(self.boxes) % 2 == 0 or (self.zydown[0] != len(self.boxes)-1):	# Only update the helix boxer if it is paired, otherwise treat it as a regular box
					hb = self.helixboxes[int(old_div(self.zydown[0],2))]
					if self.zydown[0] % 2 == 0:
						hb[5] = dz+self.zydown[3]
						hb[4] = dy+self.zydown[4]
					else:
						hb[2] =  dz+self.zydown[3]
						hb[1] = dy+self.zydown[4]
					self.update_helixbox(int(old_div(self.zydown[0],2)))
				else:
					self.firsthbclick[2] = z
					self.firsthbclick[1] = y

			self.boxes[self.zydown[0]][2]=dz+self.zydown[3]
			self.boxes[self.zydown[0]][1]=dy+self.zydown[4]
			self.update_box(self.curbox,True)

		def zy_up  (self,event):
			if self.zydown!=None:
				self.update_box(self.curbox)
			self.zydown=None

		def zy_scale(self,news):
			"xy image view has been rescaled"
			self.wscale.setValue(news)
			#self.xyview.set_scale(news,True)
			#self.xzview.set_scale(news,True)

		def zy_origin(self,newor):
			"xy origin change"
			xyo=self.xyview.get_origin()
			self.xyview.set_origin(xyo[0],newor[1],True)

			#xzo=self.xzview.get_origin()
			#self.xzview.set_origin(xzo[0],newor[1],True)

		def closeEvent(self,event):
			print("Exiting")
			self.boxviewer.close()
			self.boxesviewer.close()
	#		self.averageviewer.close()
			event.accept()
			#self.app().close_specific(self)
			self.module_closed.emit() # this signal is important when e2ctf is being used by a program running its own event loop

		#def closeEvent(self,event):
			#self.target().done()


	if options.path and options.verbose:
		print("OPTIONS.PATH IS!!!!\n\n\n", options.path)

	img = args[0]

	app = EMApp()
	if options.inmemory:


		imghdr = EMData(img,0,True)
		thisapix = imghdr['apix_x']

		if options.shrink and int(options.shrink) > 1:

			# The new shrinking scheme
			print("Shrinking, please wait :)")
			img = EMData()
			img.read_binedimage(args[0],0,options.shrink)
		else:
			print("You better have A LOT of memory (more than 8GB) if this is an un-shrunk 4k x 4k x 0.5k tomogram, because I'm loading it un-shrunk. It is wise to wait patiently.")
			img = EMData(img,0)

		if options.invert:
			img = img*(-1)
			print("The contrast of the tomogram has been reversed")

		if options.lowpass:
			filt=old_div(1.0,options.lowpass)
			print("The tomogram is being low pass filtered to %d Angstroms resolution" %(options.lowpass))
			img = img.process('filter.lowpass.gauss',{'cutoff_freq':filt,"apix":img['apix_x']})

		print("Done !")

		if options.apix:
			thisapix = options.apix
		
		box = 32
		if options.boxsize:
			box = options.boxsize
		
		boxer = EMTomoBoxer(app,data=img,datafile=None,yshort=options.yshort,apix=thisapix,boxsize=box,shrink=options.shrink,contrast=options.invert,center=options.centerbox,mod=None,normalize=options.normproc)
	else :
#		boxer=EMTomoBoxer(app,datafile=args[0],yshort=options.yshort,apix=options.apix,boxsize=options.boxsize)
		img=args[0]

		imghdr = EMData(img,0,True)
		thisapix = imghdr['apix_x']

		'''The modd variable is used as a check that determines whether a "modified"
		tomogram needs to be deleted prior to extracting boxes from disk.
		Because boxing from disk does NOT open the whole tomogram, a modified copy
		needs to be generated and written to file when you want to find
		particles in a shrunk or pre-low pass filtered tomogram, BUT still extract
		from the raw tomogram'''
		modd = False
		if options.shrink > 1:

			imgnew = img
			if '_editedtemp.' not in img:
				imgnew = img.replace('.','_editedtemp.')
			print("Shrinking, please wait :)")
			imgfile = EMData()
			imgfile.read_binedimage(img,0,options.shrink)
			imgfile.write_image(imgnew)

			img = imgnew
			modd = True

		if options.invert:
			imgnew = img
			if '_editedtemp.' not in img:
				os.system('rm *_editedtemp*')
				imgnew = img.split('/')[-1].replace('.','_editedtemp.')

			cmd = 'e2proc3d.py ' + img + ' ' + imgnew + ' --mult=-1'
			os.system(cmd)
			img = imgnew
			modd = True

		if options.lowpass:
			imgnew = img
			if '_editedtemp.' not in img:
				os.system('rm *_editedtemp*')
				imgnew = img.split('/')[-1].replace('.','_editedtemp.')
			filt=old_div(1.0,options.lowpass)

			imghdr = EMData(img,0,True)
			tapix = imghdr['apix_x']
			if options.apix:
				tapix = options.apix

			cmd = 'e2proc3d.py ' + img + ' ' + imgnew + ' --process=filter.lowpass.gauss:cutoff_freq=' + str(filt) + ':apix=' + str(tapix)
			os.system(cmd)
			img = imgnew
			modd = True

		#print "The shrink factor default is", options.shrink


		print("\nDatafile and type of datafile are", img, type(img))



		if options.apix:
			thisapix = options.apix
		box = 32
		if options.boxsize:
			box = options.boxsize
		
		boxer=EMTomoBoxer(app,data=None,datafile=img,yshort=options.yshort,apix=thisapix,boxsize=box,shrink=options.shrink,contrast=options.invert,center=options.centerbox,mod=modd,normalize=options.normproc)

		#boxer=EMTomoBoxer(app,datafile=img,yshort=options.yshort,apix=options.apix,boxsize=options.boxsize,shrink=options.shrink,contrast=options.invert,center=options.centerbox,mod=None,normalize=options.normproc)

	boxer.show()
	app.execute()

	return


if __name__ == "__main__":
	main()



