#!/usr/bin/env python
'''
====================
Author: Jesus Galaz-Montoya - 2017, Last update: Jul/2018
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
import sys, os

from EMAN2 import *
from EMAN2_utils import runcmd

import shutil


def main():

	progname = os.path.basename(sys.argv[0])
	usage = """This program requires ICON-GPU and IMOD. Wrapper to run ICON-GPU reconstructions using fewer parameters. The program automatically crops the tiltseries into a square, runs test alignments
	to make sure a good XY size is picked (some sizes make ICON crash), and picks the correct iteration numbers. 
	The reconstructed tomogram is also automatically rotated around x so that z is the shortest size using IMOD's clip rotx function, and can be optionally binned if --shrink is provided.

	Preprocessing and reconstruction normally happen in separate runs since preprocessing acts on the raw (or ideally, X-ray corrected) tiltseries, the .st file from IMOD after X-ray correction.
	On the other hand, the reconstruction should be performed with the preprocessed tiltseries AFTER it has been aligned with IMOD (the .ali file) and AFTER gold fiducials have been deleted.
	This program automatically backs up the .st and .ali IMOD files while replacing them with the ICONPreProcessed files (gold deletion needs to proceed in IMOD's pipeline) prior to reconstruction.
	An option is provided through --skipgolderasing to streamline the entire pipeline with a single run of this program if the sample has no gold fiducials or the user doesn't want to delete them.

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
	
	parser.add_argument("--highpasspixels", type=int, default=4, help="""Default=4. Number of Fourier pixels to zero out during highpass filtering AFTER reconstruction (postprocessing). Provide 'None' or '0' to disactivate.""")
	
	parser.add_argument("--iconpreproc", action='store_true', default=False, help="""default=False. If on, this will trigger a call to ICONPreProcess on the .st file supplied through --tiltseries.""")

	parser.add_argument("--lowpassresolution", type=float, default=100.0, help="""Default=100. Resolution in angstroms to lowpass filter the tomogram AFTER reconstruction (postprocessing). Type 'None' or '0' to disactivate.""")

	parser.add_argument("--path", type=str, default='icongpu', help="""Default=icongpu. Name of the directory where to store the output results. Only works when reconstructing the .ali file (preprocessing of the .st file will output the preprocessed tiltseries to the current directory). A numbered series of 'icongpu' directories will be created (i.e., if the program is run more than once, results will be stored in iconpu_01, icongpu_02, etc., directories, to avoid overwriting data).""")
	parser.add_argument("--ppid", type=int, default=-1, help="Default=-1. Set the PID of the parent process, used for cross platform PPID")

	parser.add_argument("--shrink",type=int, default=0, help="""Default=0 (not used). Shrink factor to provide IMOD's binvol program with to bin/shrink the output tomogram after rotation about the x axis.""")
	parser.add_argument("--sizez",type=int, default=0, help="""Default=0 (not used). Output size in Z for the reconstructed tomogram. This should be the same as the --thickness value provided during tiltseries preprocessing, or larger (it's good to make sure the entire reconstruction will fit in the reconstruction volume without being too tight). If running a reconstruction of the .ali file and both --thickness and --sizez are provided, the latter will supersede the former.""")
	parser.add_argument("--skipgolderasing", action='store_true', default=False, help="""default=False. If on, this will call IMOD to generate a new aligned tiltseries after ICONPreProcess, and then ICON-GPU will be automatically called to perform the reconstruction. Typically, one would NOT use this option as it is useful to delete the gold fiducials prior to reconstruction.""")

	parser.add_argument("--thickness", type=int, default=0, help="""default=0. Thickness of the specimen as seen in a preliminary weighted back projection reconstruction from IMOD (through how many slices in Z are there specimen densities?).""")
	parser.add_argument("--tiltseries", type=str, default=None, help="""default=None. .st file from IMOD if --iconpreproc is turned on. Otherwise, supply the .ali file from IMOD *after* X-ray correction, iconpreprocessing, and alignment with IMOD.""")
	parser.add_argument("--tltfile", type=str, default=None, help="""default=None. .tlt file from IMOD.""")

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness.")

	(options, args) = parser.parse_args()	
	
	logger = E2init(sys.argv, options.ppid)

	if not options.tiltseries:
		print("\nERROR: --tiltseries required")
		sys.exit(1)

	filename, extension = os.path.splitext(options.tiltseries)
	
	alifile = options.tiltseries.replace(extension,'.ali')
	
	c=os.getcwd()
	findir=os.listdir(c)

	if not options.tltfile:
		anglesfile = os.path.basename(options.tiltseries.replace(extension,'.tlt'))
		if anglesfile in findir:
			options.tltfile = anglesfile
		else:
			print("\nERROR: expected tlt file = {}, (text file with the list of tilt angles) not found. Supply --tltfile explicitly.".format(anglesfile))
			sys.exit(1)

	if options.iconpreproc and alifile not in findir:
		print("\nERROR: the aligned tiltseries must be in the same directory, and should match the name of the raw .st tiltseries, except that the extension should be .ali instead of .st; the expected file is {}".format(alifile))
		sys.exit(1)

	if options.verbose:
		print("\n(e2tomo_icongpu)(icongpufunc) making directory {} to store reconstruction results".format(options.path))
	
	from EMAN2_utils import makepath
	options = makepath(options)

	cmdsfilepath = options.path + '/cmds.txt'

	if options.iconpreproc:

		if '.st' not in extension:
			print("\nERROR: the extension of the --tiltseries is {} instead of .st; make sure this is the correct tiltseries, and change the extension to .st".format(extension))
			shutil.rmtree(options.path)
			sys.exit(1)
		
		if not options.thickness:
			print("\nERROR: --thickness required for ICONPreProcess.")
			shutil.rmtree(options.path)
			sys.exit(1)

		iconpreprocfunc(options,alifile,extension,cmdsfilepath)
	
	elif not options.iconpreproc or options.skipgolderasing:
		
		if not options.thickness and not options.sizez:
			print("\nERROR: --thickness or --sizez required")
			sys.exit(1)
		elif options.thickness and options.sizez:
			print("\nWARNING: --thickness={} and --sizez={} were both provided; only --sizez={} will be used for reconstruction".format(options.thickness,options.sizez,options.sizez))
			options.thickness = options.sizez

		if not options.tltfile:
			print("\nWARNING: --tltfile not provided. The program will attempt to find it automatically")
			tltfile = options.tiltseries.replace(extension,'.tlt')
			
			if tltfile not in findir:
				print("\nERROR: in the abscence of --tltfile, the .tlt file with tiltseries angles must be in the running (current) directory, and should match the name of the raw .st tiltseries, except that the extension should be .tlt instead of .st; the expected file is {}".format(alifile))
				sys.exit(1)

		if not options.iconpreproc and '.ali' not in extension:
			print("\nWARNING: the extension of the --tiltseries is {} instead of .ali; make sure this is the correct tiltseries".format(extension))
		
		if not options.iconpreproc:
			alifile = options.tiltseries

		icongpufunc(options,alifile,cmdsfilepath)

	E2end(logger)


	return


def iconpreprocfunc(options,alifile,extension,cmdsfilepath):

	#from shutil import copyfile
	
	outfile = options.tiltseries.replace(extension,'_iconpreproc_th' + str(options.thickness) + '.st')
	backupst = 'backup.' + options.tiltseries
	backupali = 'backup.' + alifile

	if options.verbose:
		print("\n(e2tomo_icongpu)(iconpreprocfunc) running ICONPreProcess")
	cmd1 = "ICONPreProcess -input " + options.tiltseries + " -tiltfile " + options.tltfile + " -thickness " + str(options.thickness) + " -output " + outfile
	runcmd(options,cmd1,cmdsfilepath)

	
	if options.verbose:
		print("\n(e2tomo_icongpu)(iconpreprocfunc) backing up raw (x-ray corrected) tiltseries {} to file {}".format(options.tiltseries,backupst))
	shutil.copyfile(options.tiltseries, backupst)
	#cmd2 = "cp " + options.tiltseries + " " + backupst
	#runcmd(options,cmd2,cmdsfilepath)

	if options.verbose:
		print("\n(e2tomo_icongpu)(iconpreprocfunc) copying the preprocessed tiltseries {} to the original tiltseries file {}".format(outfile,options.tiltseries))
	shutil.copyfile(outfile, options.tiltseries)
	#cmd4 = "cp " + outfile + " " + options.tiltseries
	#runcmd(options,cmd4,cmdsfilepath)

	if options.verbose:
		print("\n(e2tomo_icongpu)(iconpreprocfunc) backing up aligned tiltseries {} to file {}".format(alifile,backupali))
	shutil.copyfile(alifile, backupali)
	#cmd3 = "cp " + alifile + " " + backupali
	#runcmd(options,cmd3,cmdsfilepath)

	if options.verbose:
		print("\n(e2tomo_icongpu)(iconpreprocfunc) generating new .ali file after ICONPreProcess")
	cmd2 = "subm newst"
	runcmd(options,cmd2,cmdsfilepath)

	return


def icongpufunc(options,alifile,cmdsfilepath):
	
	alifilename, aliextension = os.path.splitext(alifile)

	hdr = EMData(alifile,0,True)
	nx = hdr['nx']
	ny = hdr['ny']

	outsize = nx

	if nx != ny:
	
		if ny < nx:
			outsize = ny

		if outsize % 2:
			outsize -= 1
	
		alifile = cropper(options,alifile,aliextension,outsize,cmdsfilepath)

	alifile,outsize,iterationsstring = icontest(options,alifile,outsize,cmdsfilepath,aliextension)

	icondir = options.path

	thickness = old_div(float(outsize),4.0)

	if options.thickness:
		thickness = float(options.thickness)
		#print "\noptions.thickness has changed thickness to {}".format(thickness)
	if options.sizez:
		thickness=float(options.sizez)
		#print "\noptions.sizez has changed thickness to {}".format(thickness)

	thickenss=int(round(thickness))
	print("\naaaafter options, thickness is {}".format(thickness))

	outtomogram = alifile.replace(aliextension,'_icongpu.mrc')
 	
	if options.verbose:
		print("\n(e2tomo_icongpu)(icongpufunc) calling ICON-GPU.")
	cmdicon1 = 'ICON-GPU -input ' + alifile + ' -tiltfile ' + options.tltfile + ' -outputPath ' + icondir + ' -slice 0,' + str(outsize-1) + ' -ICONIteration '+ iterationsstring + ' -dataType 1 -threshold 0 -gpu ' + options.gpus 
	runcmd(options,cmdicon1,cmdsfilepath)

	findir=os.listdir(icondir+'/crossValidation')

	if 'crossV.frc' in findir:
		sigma=1000000000
		sigmanonzero=100000000
		try:
			hdr = EMData( icondir + '/crossValidation/crossV.frc' )
			sigma = hdr['sigma']
			sigmanonzero = hdr['sigma_nonzero']
		except:
			with open(icondir + '/crossValidation/crossV.frc', 'r') as crossVfrcfile: 
				lines = crossVfrcfile.readlines()

				if not lines: sigma=sigmanonzero=0.0
				elif int(len(lines)) < int(old_div(nx,2.0) -1.0):
					sigma=0.0
					sigmanonzero=0.0

		if float(sigma) == 0.0 or float(sigmanonzero) == 0.0:
			print("\nWARNING: the crossV.frc file generated by ICON-GPU is empty; generating it 'manually' with e2proc3d")
			cmdcrossvfrc = 'cd ' + icondir + '/crossValidation && e2proc3d.py crossV_reprojection.mrc --calcfsc=GroundTruth.mrc crossV.frc'
			runcmd(options,cmdcrossvfrc,cmdsfilepath)

			#hdr = EMData( icondir + '/crossValidation/crossV.frc' )
			#sigma = hdr['sigma']
			#sigmanonzero = hdr['sigma_nonzero']

			#if float(sigma) == 0.0 or float(sigmanonzero) == 0.0:
			#	print("\nERROR: manual computation of crossV.frc with e2proc3d failed")
			#	sys.exit(1)


	elif 'crossV.frc' not in findir:
		print("\nERROR: no crossV.frc file found inside {}/crossValidation directory".format(icondir))
		sys.exit(1)	
		
	if options.verbose:
		print("\n(e2tomo_icongpu)(icongpufunc) calling ICONMask3.")
	cmdicon2 = 'ICONMask3 -inputPath ' + icondir + '/reconstruction -tiltfile ' + options.tltfile + ' -output ' + outtomogram + ' -slice 0,' + str(outsize-1) + ' -thickness ' + str(thickness) + ' -crossVfrc ' + icondir + '/crossValidation/crossV.frc -fullRecfrc ' + icondir + '/crossValidation/fullRec.frc' 
	runcmd(options,cmdicon2,cmdsfilepath)

	outtomogramzshort = outtomogram.replace('.mrc','_ZSHORT.mrc') 
	if options.verbose:
		print("\n(e2tomo_icongpu)(icongpufunc) calling IMOD to rotate the reconstructed volume around x and shrink it if --shrink > 1 was specified.")
	cmdimod1 = 'clip rotx ' + outtomogram + ' ' + outtomogramzshort
	runcmd(options,cmdimod1,cmdsfilepath)
	os.remove(outtomogram)
	sourcetomo = outtomogramzshort

	outtomogramzshortbinned = outtomogramzshort.replace('.mrc','_bin'+str(options.shrink)+'.mrc') 
	if options.shrink and int(options.shrink) > 1:
		cmdimod2 = 'binvol ' + outtomogramzshort + ' ' + outtomogramzshortbinned + ' --binning ' + str(options.shrink) + ' --antialias -1'
		runcmd(options,cmdimod2,cmdsfilepath)
		sourcetomo = outtomogramzshortbinned

	outtomogramzshortpostproc = sourcetomo
	if options.lowpassresolution and not options.highpasspixels:
		outtomogramzshortpostproc = sourcetomo.replace('.mrc', '_lp' + str(options.lowpassresolution) )
		cmdeman2 = 'e2proc3d.py ' + sourcetomo + ' ' + outtomogramzshortpostproc + ' --process filter.lowpass.tanh:cutoff_freq=' + str(old_div(1.0,options.lowpassresolution)) 
		runcmd(options,cmdeman2,cmdsfilepath)
	
	elif options.highpasspixels and not options.lowpassresolution:
		outtomogramzshortpostproc = sourcetomo.replace('.mrc', '_hpp' + str(options.highpasspixels) + '.mrc')
		cmdeman2 = 'e2proc3d.py ' + sourcetomo + ' ' + outtomogramzshortpostproc + ' --process filter.highpass.gauss:cutoff_pixels=' + str(options.highpasspixels)
		runcmd(options,cmdeman2,cmdsfilepath)

	elif options.lowpassresolution and options.highpasspixels:
		outtomogramzshortpostproc = sourcetomo.replace('.mrc', '_lp' + str(int(round(options.lowpassresolution))) + '_hpp' + str(int(options.highpasspixels)) + '.mrc')
		cmdeman2 = 'e2proc3d.py ' + sourcetomo + ' ' + outtomogramzshortpostproc + ' --process filter.lowpass.tanh:cutoff_freq=' + str(old_div(1.0,options.lowpassresolution)) + ' --process filter.highpass.gauss:cutoff_pixels:' + str(options.highpasspixels)
		runcmd(options,cmdeman2,cmdsfilepath)
	
	os.rename(outtomogramzshort, options.path + '/' + outtomogramzshort)
	
	if options.shrink and options.shrink > 1:
		os.rename(outtomogramzshortbinned, options.path + '/' + outtomogramzshortbinned)
	
	if options.lowpassresolution or options.highpasspixels:
		os.rename(outtomogramzshortpostproc, options.path + '/' + outtomogramzshortpostproc)

	return


def icontest(options,alifile,outsize,cmdsfilepath,aliextension):
	
	passtest = False

	while not passtest and outsize > 63:
			
		it1,it2,it3 = calciterations(outsize)
		
		iterationsstring = str(it1)+','+str(it2)+','+str(it3)
		print("iterationsstring is {} of type {} ".format(iterationsstring,type(iterationsstring)))

		tmpdir = "tmpicontest"
		
		#print "options.gpus={}, type={}".format(options.gpus,type(options.gpus))
		#print "tmpdir={}, type={}".format(tmpdir,type(tmpdir))
		#print "alifile={}, type={}".format(alifile,type(alifile))
		#print "options.tltfile={}, type={}".format(options.tltfile,type(options.tltfile))
		
		icongpucmd = "ICON-GPU -input " + alifile + " -tiltfile " + options.tltfile + " -outputPath " + tmpdir + " -slice 0,1 -ICONIteration " + iterationsstring + " -dataType 1 -threshold 0 -gpu " + options.gpus
		cmdicontest = "mkdir " + tmpdir + " ; " + icongpucmd
		
		try:
			runcmd(options,cmdicontest,cmdsfilepath)

			findirtest = os.listdir(tmpdir+'/reconstruction/')
		
			for f in findirtest:
				if '.mrc' in f[-4:]:
					imgpath = tmpdir+'/reconstruction/'+f
					img = EMData(imgpath,0)
					sigma = img['sigma']
					sigmanonzero = img['sigma_nonzero']

					shutil.rmtree(tmpdir)
					if sigma and sigmanonzero:
						passtest = True
						print("\nthe test passed; the tiltseries has a good size now nx={}, ny={}".format(img['nx'],img['ny']))
						return alifile,outsize,iterationsstring
					else:
						passtest = False
						outsize -= 2
						print("\nICON-GPU failed becase the image size was bad; cropping the images in the tiltseries to this new size, nx={}, nx={}".format(outsize,outsize))
						alifile = cropper(options,alifile,aliextension,outsize,cmdsfilepath)	
						break
		except:
			print("\nWARNING: ICON-GPU command failed; the command run was cmd={}".format(icongpucmd))
		
	return


def calciterations(outsize):
	#iteration guide from ICON creators, tweaked empirically
	#512*512, -ICONIteration 10,100,60
	#1024*1024, -ICONIteration 10,140,80
	#2048*2048 -ICONIteration 10,180,110 
	#4096*4096 -ICONIteration 10,260,140
	
	it1 = 10

	if outsize < 513:
		it2,it3 = scaleiterations(outsize,512,100,60)
	
	elif outsize > 512 and outsize < 1025:
		it2,it3 = scaleiterations(outsize,1024,140,80)

	elif outsize > 1024 and outsize < 2049:
		it2,it3 = scaleiterations(outsize,2048,190,110)

	elif outsize > 2048 and outsize < 4097:
		it2,it3 = scaleiterations(outsize,4096,260,150)
	
	elif outsize > 4096:
		it2,it3 = scaleiterations(outsize,8192,350,200)

	return it1,it2,it3


def scaleiterations(outsize,upperlimit,it2base,it3base):
	sizefactor = old_div(float(outsize),float(upperlimit))
	it2 = int(ceil(it2base * sizefactor))
	it3 = int(ceil(it3base * sizefactor))

	return it2,it3


def cropper(options,alifile,extension,outsize,cmdsfilepath):
	
	hdr = EMData(alifile,0,True)
	nx=hdr['nx']
	ny=hdr['ny']

	if options.verbose:
		print("\nWARNING: the tiltseries will be resized since its images are not squares. The original size is nx={}, ny={}, and will be cropped to nx={}, ny={}".format(nx,ny,outsize,outsize))
	
	outcrop = alifile.replace(extension,'_clip' + str(outsize) + extension)
	if '_clip' in alifile:
		outcrop = alifile.split('_clip')[0]+extension 
		outcrop = outcrop.replace(extension,'_clip' + str(outsize) + extension)

	cmdcrop = "e2proc2d.py " + alifile + ' ' + outcrop + ' --clip ' + str(outsize) + ',' + str(outsize)
	
	runcmd(options,cmdcrop,cmdsfilepath)

	return outcrop


#def runcmd(options,cmd,cmdsfilepath=''):
#	if options.verbose > 9:
#		print("(e2tomo_icongpu)(runcmd) running command", cmd)
#	
#	p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#	text=p.communicate()	
#	p.stdout.close()
#	
#	if options.verbose > 8:
#		print("(e2segmask)(runcmd) done")
#	
#	if cmdsfilepath:
#		with open(cmdsfilepath,'a') as cmdfile: cmdfile.write( cmd + '\n')
#
#	return 1


if __name__ == '__main__':
	main()
