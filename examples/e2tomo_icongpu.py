#!/usr/bin/env python
'''
====================
Author: Jesus Galaz-Montoya - 2017, Last update: 12/Sep/2017
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

#import matplotlib
#matplotlib.use('Agg',warn=False)

#import matplotlib.pyplot as plt
#import pylab

#import matplotlib.colors as mcol
#import matplotlib.cm as cm

#import colorsys

import sys, os

from EMAN2 import *


def main():

	progname = os.path.basename(sys.argv[0])
	usage = """Wrapper to run ICOn reconstructions using fewer parameters. The program automatically crops the tiltseries into a square, runs tests alignments
	to make sure a good XY size is picked (some sizes make ICON crash), and picks the correct iteration numbers.

	Preprocessing and reconstruction need to happen in separate runs since preprocessing acts on the raw (or ideally, X-ray corrected) tiltseries, the .st file from IMOD after X-ray correction.
	On the other hand, the reconstruction should be performed with the preprocessed tiltseries AFTER it has been aligned with IMOD, that is, the .ali file, and AFTER gold fiducials have been deleted.
	This program automatically backs up the .st and .ali IMOD files while replacing them with the ICONPreProcessed files (gold deletion needs to proceed in IMOD's pipeline) prior to

	To preprocess the raw tiltseries, run:

	e2tomo_icongpu.py --tiltseries=<.st file> --iconpreproc --thickness=<value in pixels>

	This MUST be run from the IMOD reconstruction directory, with all alignment files in it, so that the preprocessed tiltseries can be realigned automatically without going through IMOD's ETOMO pipeline again.


	For reconstruction, run:

	e2tomo_icongpu.py --tiltseries=<.ali file> --sizez=<size of output tomogram in z>

	The .tlt file should be located automatically if it has the same name (but with .tlt extension) as the --tiltseries, as should be the case for an intact IMOD directory.
	Otherwise, to supply an alternative .tlt file, add --tltfile=<.tlt file from imod> to the command above.

	(This can be run in any given directory as long as the .ali and .tlt files are there).


	If you don't need or want to delete gold fiducials and want to run BOTH steps at once (ICONPreProcess and reconstruction with ICON-GPU), add --skipgolderasing
	"""
		
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--gpus", type=str, default="-1", help="""default=-1 (all available GPUs will be used). To select specific GPUs to use, provide them as a comma-separated list of integers.""")
	
	parser.add_argument("--iconpreproc", action='store_true', default=False, help="""default=False. If on, this will trigger a call to ICONPreProcess on the .st file supplied through --tiltseries. This""")
		
	parser.add_argument("--path", type=str, default='plotfig', help="""Default=icongpu. Name of the directory where to store the output results. Only works when reconstructing the .ali file (preprocessing of the .st file will output the preprocessed tiltseries to the current directory). A numbered series of 'icongpu' directories will be created (i.e., if the program is run more than once, results will be stored in iconpu_01, icongpu_02, etc., directories, to avoid overwriting data).""")
	parser.add_argument("--ppid", type=int, default=-1, help="Default=-1. Set the PID of the parent process, used for cross platform PPID")

	parser.add_argument("--sizez",type=int, default=0, help="""Default=0. Output size in Z for the reconstructed tomogram. This should be the same as the --thickness value provided during tiltseries preprocessing, or larger (it's good to make sure the entire reconstruction will fit in the reconstruction volume without being too tight). If running a reconstruction of the .ali file and both --thickness and --sizez are provided, the latter will supersede the former.""")

	parser.add_argument("--thickness", type=int, default=0, help="""default=0. Thickness of the specimen as seen in a preliminary weighted back projection reconstruction from IMOD (through how many slices in Z are there specimen densities?).""")
	parser.add_argument("--tiltseries", type=str, default='', help="""default=None. .st file from IMOD if --iconpreproc is turned on. Otherwise, supply the .ali file from IMOD *after* X-ray correction, iconpreprocessing, and alignment with IMOD.""")
	parser.add_argument("--tltfile", type=str, default='', help="""default=None. .tlt file from IMOD.""")

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness.")

	(options, args) = parser.parse_args()	
	
	logger = E2init(sys.argv, options.ppid)

	if not options.tiltseries:
		print "\nERROR: --tiltseries required"
		sys.exit(1)

	filename, extension = os.path.splitext(options.tiltseries)
	
	alifile = options.tiltseries.replace(extension,'.ali')
	
	c=os.getcwd()
	findir=os.listdir(c)

	if options.iconpreproc and alifile not in findir:
		print "\nERROR: the aligned tiltseries must be in the same directory, and should match the name of the raw .st tiltseries, except that the extension should be .ali instead of .st; the expected file is {}".format(alifile)
		sys.exit(1)

	if options.iconpreproc:

		if '.st' not in extension:
			print "\nERROR: the extension of the --tiltseries is {} instead of .st; make sure this is the correct tiltseries, and change the extension to .st".format(extension)
			sys.exit(1)
		
		if not options.thickness:
			print "\nERROR: --thickness required for ICONPreProcess."
			sys.exit(1)

		if options.skipgolderasing:
			if not options.tltfile:

		iconpreprocfunc(options,alifile)
	
	elif not options.iconpreproc or options.skipgolderasing:
		
		if not options.thickness and not options.sizez:
			print "\nERROR: --thickness or --sizez required"
			sys.exit(1)
		elif options.thickness and options.sizez:
			print "\nWARNING: --thickness={} and --sizez={} were both provided; only --sizez={} will be used for reconstruction".format(options.thickness,options.sizez,options.sizez)
		
		if not options.tltfile:
			print "\nWARNING: --tltfile not provided. The program will attempt to find it automatically"
			tltfile = options.tiltseries.replace(extension,'.tlt')
			
			if tltfile not in findir:
				print "\nERROR: in the abscence of --tltfile, the .tlt file with tiltseries angles must be in the running (current) directory, and should match the name of the raw .st tiltseries, except that the extension should be .tlt instead of .st; the expected file is {}".format(alifile)
				sys.exit(1)

		if not options.iconpreproc and '.ali' not in extension:
			print "\nWARNING: the extension of the --tiltseries is {} instead of .ali; make sure this is the correct tiltseries".format(extension)
		
		if not options.iconpreproc:
			alifile = options.tiltseries

		icongpufunc(options,alifile,extension)

	E2end(logger)


	return


def iconpreprocfunc(options,alifile):
	
	outfile = options.tiltseries.replace(extension,'_iconpreproc_th' + str(options.thickness) + '.st')
	backupst = 'backup.' + options.tiltseries
	backupali = 'backup.' + alifile

	cmd1 = "ICONPreProcess -input " + options.tiltseries + " -tiltfile " + options.tltfile + " -thickness " + str(options.thickness) + " -output " + outfile
	runcmd(options,cmd1)
	if options.verbose:
		print "\n(e2tomo_icongpu)(iconpreprocfunc) running ICONPreProcess"
	
	cmd2 = "cp " + options.tiltseries + " " + backupst
	runcmd(options,cmd2)
	if options.verbose:
		print "\n(e2tomo_icongpu)(iconpreprocfunc) backing up raw (x-ray corrected) tiltseries {} to file {}".format(options.tiltseries,backupst)

	cmd3 = "cp " + alifile + " " + backupali
	runcmd(options,cmd3)
	if options.verbose:
		print "\n(e2tomo_icongpu)(iconpreprocfunc) backing up aligned tiltseries {} to file {}".format(alifile,backupali)

	cmd4 = "subm newst"
	runcmd(options,cmd4)
	if options.verbose:
		print "\n(e2tomo_icongpu)(iconpreprocfunc) generating new .ali file after ICONPreProcess"

	return


def icongpufunc(options,alifile,extension):
	
	if options.verbose:
		print "\n(e2tomo_icongpu)(icongpufunc) making directory {} to store reconstruction results".format(options.path)
	
	from EMAN2_utils import makepath
	options = makepath(options) 

	hdr = EMData(alifile,0,True)
	nx = hdr['nx']
	ny = hdr['ny']

	outsize = nx

	if nx != ny:
	
		if ny < nx:
			outsize = ny

		if outsize % 2:
			outsize -= 1
	
		alifile = cropper(options,alifile,extension,nx,ny,outsize)

	
	passtest = False

	while not passtest and outsize > 63:
			
		it1,i2,i3 = calciterations(outsize)
		
		iterationsstring = str(it1)+','+str(it2)+','+str(it3)
	
		tmpdir = "tmpicontest"
		
		cmdicon = "mkdir " + tmpdir + " ; ICON-GPU -input " + alifile + " -tiltfile " + options.tltfile + " -outputPath " + tmpdir + " -slice 0,1 -ICONIteration " + iterationsstring + " -dataType 1 -threshold 0 -gpu " + options.gpus

		runcmd(cmdicon)

		findirtest = os.listdir(tmpdir+'/reconstruction/')
		
		for f in findirtest:
			if '.mrc' in f[-4:]:
				imgpath = tmpdir+'/reconstruction/'+f
				img = EMData(imgpath,0)
				sigma = img['sigma']
				sigmanonzero = hdrtest['sigma_nonzero']

				if sigma and sigmanonzero:
					passtest = True
				else:
					passtest = False
					outsize -= 2
					print "\nICON-GPU failed becase the image size was bad; cropping the images in the tiltseries to this new size, nx={}, nx={}".format(outsize,outsize)
					alifile = cropper(options,alifile,extension,nx,ny,outsize)
					os.remove(tmpdir)
					break

	if passtest:
		icondir = options.path


	return


def calciterations(outsize):
	#iteration guide from ICON creators, tweaked empirically
	#512*512, -ICONIteration 10,100,60
	#1024*1024, -ICONIteration 10,140,80
	#2048*2048 -ICONIteration 10,180,110 
	#4096*4096 -ICONIteration 10,260,140
	
	it1 = 10

	if outsize < 513:
		it1,it2 = scaleiterations(outsize,512,100,60)
	
	elif outsize > 512 and outsize < 1025:
		it1,it2 = scaleiterations(outsize,1024,140,80)

	elif outsize > 1024 and outsize < 2049:
		it1,it2 = scaleiterations(outsize,2048,190,110)

	elif outsize > 2048 and outsize < 4097:
		it1,it2 = scaleiterations(outsize,4096,260,150)
	
	elif outsize > 4096:
		it1,it2 = scaleiterations(outsize,8192,350,200)

	return it1,it2,it3


def scaleiterations(outsize,upperlimit,it2base,it3base):
	sizefactor = outsize/upperlimit
	it2 = int(ceil(it2base * sizefactor))
	it3 = int(ceil(it3base * sizefactor))

	return it2,it3


def cropper(options,alifile,extension,nx,ny,outsize):

	if options.verbose:
		print "\nWARNING: the tiltseries will be resized since its images are not squares. The original size is nx={}, ny={}, and will be cropped to nx={}, ny={}".format(nx,ny,outsize,outsize)
		
	outcrop = alifile.replace(extension,'_clip' + str(outsize) + extension)
	cmdcrop = "e2proc2d.py " + alifile + ' ' + outcrop + ' --clip ' + str(outsize) + ',' + str(outsize)
	
	runcmd(options,cmdcrop)

	return outcrop,outsize


def runcmd(options,cmd):
	if options.verbose > 9:
		print "(e2tomo_icongpu)(runcmd) running command", cmd
	
	p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	text=p.communicate()	
	p.stdout.close()
	
	if options.verbose > 8:
		print "(e2segmask)(runcmd) done"
	
	return 1


if __name__ == '__main__':
	main()