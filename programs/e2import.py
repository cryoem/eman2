#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
#
# Author: John Flanagan (jfflanag@bcm.edu)
# Copyright (c) 2000-2011 Baylor College of Medicine


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
#

from past.utils import old_div
from builtins import range
import os, shutil, glob
from EMAN2 import *
from EMAN2star import StarFile
import numpy as np

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] files
	This program performs a variety of tasks for getting data or metadata from other programs into an EMAN2 project.

	import_movies - imports DDD movie data for a cryoEM project (--importation copy recommended)

	import_particles - will simply copy a set of per-micrograph particle files into EMAN2.1's preferred HDF format in particles/
	import_boxes - will read EMAN1 '.box' files (text files containing coordinates) into appropriate info/*json files (see --box_type)
	import_eman1 - will convert a typical EMAN1 phase-flipped start.hed/img file into an EMAN2 project (converting files, fixing CTF, splitting, ...)

	import_tomos - flag to handle imported files as subtomogams for a SPT project (see also --importation)
	import_serialem - flag to handle imported files as SerialEM mdoc files
	import_fei - flag to handle imported files as FEI tomography metadata files
	import_ucsftomo - flag to handle imported files as UCSF tomo metadata files
	import_tiltseries - imports tilt series for a tomography project (--importation copy recommended)
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_pos_argument(name="files",help="List the files to import here.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=True)",  row=0, col=0, rowspan=1, colspan=3, nosharedb=True, mode='coords,parts,tomos,eman1,movies,rawtilts,meta,tiltseries')

	parser.add_header(name="filterheader", help='Options below this label are specific to e2import', title="### e2import options ###", row=2, col=0, rowspan=1, colspan=2, mode='coords,parts')

	# Type Flags
	parser.add_argument("--import_movies",action="store_true",help="Import DDD movies",default=False, guitype='boolbox', row=3, col=0, rowspan=1, colspan=1, mode='movies[True]')

	parser.add_argument("--darkrefs",help="Specify a comma separated list of dark refereence stacks/images to import. Files will be placed in movierefs_raw. See --importation for additional options.",default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=True)", row=4, col=0, rowspan=1, colspan=2, mode='movies')
	parser.add_argument("--gainrefs",help="Specify a comma separated list of gain refereence stacks/images to import. Files will be placed in movierefs_raw. See --importation for additional options.",default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=True)", row=5, col=0, rowspan=1, colspan=2, mode='movies')

	#parser.add_argument("--import_rawtilts",action="store_true",help="Import tilt images",default=False, guitype='boolbox', row=4, col=2, rowspan=1, colspan=1, mode='rawtilts[True]')
	parser.add_argument("--apix",help="Specify the apix of the tiltseries you are importing. If -1 (default), the apix in the header will not be changed.",type=float,default=-1,guitype='floatbox', row=5, col=1, rowspan=1, colspan=1,mode='tiltseries[-1]')

	parser.add_argument("--import_tiltseries",action="store_true",help="Import tiltseries",default=False, guitype='boolbox', row=5, col=2, rowspan=1, colspan=1, mode='tiltseries[True]')
	parser.add_argument("--import_tomos",action="store_true",help="Import tomograms for segmentation and/or subtomogram averaging",default=False, guitype='boolbox', row=4, col=2, rowspan=1, colspan=1, mode='tomos[True]')

	#parser.add_pos_argument(name="tilt_angles",help="Specify a file containing tilt angles corresponding to the input tilt images.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=True)",  row=0, col=0, rowspan=1, colspan=2, nosharedb=True, mode='rawtilts')

	#parser.add_argument(name="--rawtlt",help="List the text file containing tilt angles for the tiltseries to be imported.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)",  row=3, col=0, rowspan=1, colspan=3, nosharedb=True, mode='tiltseries')
	#parser.add_argument(name="--start",help="First tilt angle. Increment determined by number of tilts. Custom tilt angles can be specified by a tilt angles text file.", default="", guitype='floatbox', row=4, col=0, rowspan=1, colspan=3, nosharedb=True, mode='tiltseries')
	#parser.add_argument(name="--stop",help="Final tilt angle. Increment determined by number of tilts. Custom tilt angles can be specified by a tilt angles text file.", default="", guitype='floatbox', row=4, col=1, rowspan=1, colspan=3, nosharedb=True, mode='tiltseries')

	# parser.add_argument("--serialem_mdoc",action="store_true",help="Import metadata from corresponding SerialEM '.mdoc' files.",default=False, guitype='boolbox', row=1, col=0, rowspan=1, colspan=3, mode='meta')
	# parser.add_argument("--fei_tomo",action="store_true",help="Import metadata from corresponding FEI tomography files.",default=False, guitype='boolbox', row=2, col=0, rowspan=1, colspan=3, mode='meta')
	# parser.add_argument("--ucsf_tomo",action="store_true",help="Import metadata from corresponding UCSF Tomo files.",default=False, guitype='boolbox', row=3, col=0, rowspan=1, colspan=3, mode='meta')

	parser.add_argument("--import_particles",action="store_true",help="Import particles",default=False, guitype='boolbox', row=3, col=0, rowspan=1, colspan=1, mode='parts[True]')
	parser.add_argument("--import_eman1",action="store_true",help="This will import a phase-flipped particle stack from EMAN1",default=False, guitype='boolbox', row=3, col=0, rowspan=1, colspan=1, mode='eman1[True]')

	parser.add_argument("--importation",help="Specify import mode: move, copy or link",default='copy',guitype='combobox',choicelist='["move","copy","link"]',row=9,col=0,rowspan=1,colspan=2, mode='tomos["copy"],rawtilts["copy"],movies["move"],tiltseries["copy"]',choices=["move","copy","link"])

	parser.add_argument("--invert",action="store_true",help="Invert the contrast before importing tomograms",default=False, guitype='boolbox', row=5, col=0, rowspan=1, colspan=1, mode='tomos,rawtilts,tiltseries')
	#parser.add_argument("--tomoseg_auto",action="store_true",help="Default process for tomogram segmentation, including lowpass, highpass, normalize, clampminmax.",default=True, guitype='boolbox', row=4, col=1, rowspan=1, colspan=1, mode='tomos,rawtilts,tiltseries')
	parser.add_argument("--shrink",type=int,help="Shrink tomograms before importing. Does not work while not copying.",default=1, guitype='intbox', row=6, col=0, rowspan=1, colspan=1, mode='tomos')
	#parser.add_argument("--preprocess",type=str,help="Other pre-processing operation before importing tomograms. Dose not work while not copying.",default="", guitype='strbox', row=6, col=0, rowspan=1, colspan=2, mode='tomos,rawtilts,tiltseries')

	parser.add_argument("--import_boxes",action="store_true",help="Import boxes",default=False, guitype='boolbox', row=3, col=0, rowspan=1, colspan=1, mode='coords[True]')
	parser.add_argument("--extension",type=str,help="Extension of the micrographs that the boxes match", default='dm3')
	parser.add_argument("--box_type",help="Type of boxes to import, normally boxes, but for tilted data use tiltedboxes, and untiltedboxes for the tilted  particle partner",default="boxes",guitype='combobox',choicelist='["boxes","coords","relion_star","tiltedboxes","untiltedboxes"]',row=3,col=1,rowspan=1,colspan=1, mode="coords['boxes']")
	parser.add_argument("--boxsize",help="Specify the boxsize for each particle.",type=int,default=256)
	parser.add_argument("--curdefocushint",action="store_true",help="Used with import_eman1, will use EMAN1 defocus as starting point",default=False, guitype='boolbox', row=5, col=0, rowspan=1, colspan=1, mode='eman1[True]')
	parser.add_argument("--curdefocusfix",action="store_true",help="Used with import_eman1, will use EMAN1 defocus unchanged (+-.001 um)",default=False, guitype='boolbox', row=5, col=1, rowspan=1, colspan=1, mode='eman1[False]')

	parser.add_argument("--threads", default=1,type=int,help="Number of threads to run in parallel on a single computer when multi-computer parallelism isn't useful",guitype='intbox', row=7, col=0, rowspan=1, colspan=1, mode='eman1[1]')
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv,options.ppid)

	# Import EMAN1
	# will read start.hed/img, split by micrograph (based on defocus), and reprocess CTF in EMAN2 style
	if options.import_eman1 :
		try:
			n=EMUtil.get_image_count(args[0])
		except:
			print("Error, couldn't read images from: ",args[0])
			sys.exit(1)

		try:
			img=EMData(args[0],0)
			ctf=img["ctf"]
		except:
			print("Error, start.hed/img must be phase-flipped to import")
			sys.exit(1)

		db=js_open_dict("info/project.json")
		db["global.apix"]=ctf.apix
		db["global.cs"]=ctf.cs
		db["global.voltage"]=ctf.voltage

		try: os.mkdir("particles")
		except: pass

		imgnum=0
		lastdf=-1.0
		for i in range(n):
			img=EMData(args[0],i)
			ctf=img["ctf"]
			img.del_attr("ctf")
			fft1=img.do_fft()
			if ctf.defocus!=lastdf :
				imgnum+=1
				if options.verbose>0: print("Defocus {:4.2f} particles{:03d}".format(ctf.defocus,imgnum))
				db=js_open_dict("info/particles{:03d}_info.json".format(imgnum))
				ctf2=EMAN2Ctf()
				ctf2.defocus=ctf.defocus
				ctf2.cs=ctf.cs
				ctf2.apix=ctf.apix
				ctf2.voltage=ctf.voltage
				ctf2.ampcont=ctf.ampcont
				ctf2.dfdiff=0
				ctf2.dfang=0
				db["ctf"]=[ctf2]
				db.close()

				flipim=fft1.copy()
				ctf2.compute_2d_complex(flipim,Ctf.CtfType.CTF_SIGN)

			lastdf=ctf.defocus

			# unflip the EMAN1 phases (hopefully accurate enough)
			fft1.mult(flipim)
			img=fft1.do_ift()
			img.write_image("particles/particles{:03d}.hdf".format(imgnum),-1)		# append particle to stack

		if options.curdefocusfix:
			flag="--curdefocusfix"
			rbysnr=" "
		elif options.curdefocushint:
			flag="--curdefocushint"
			rbysnr="--refinebysnr"
		else:
			flag=""
			rbysnr="--refinebysnr"

		# fill in the needed CTF info
		launch_childprocess("e2ctf.py --autofit {} --allparticles --threads {} --voltage {} --cs {} --ac {} --apix {} --computesf".format(flag,options.threads,ctf.voltage,ctf.cs,ctf.ampcont,ctf.apix))
		launch_childprocess("e2ctf.py --autofit {} --allparticles --threads {} --voltage {} --cs {} --ac {} --apix {}".format(flag,options.threads,ctf.voltage,ctf.cs,ctf.ampcont,ctf.apix))

		# reflip flip the phases, and make "proc" images
		launch_childprocess("e2ctf.py {} --phaseflip --allparticles --phaseflipproc filter.highpass.gauss:cutoff_freq=0.005 --phaseflipproc2 filter.lowpass.gauss:cutoff_freq=0.08 --phaseflipproc3 math.meanshrink:n=2".format(rbysnr))

		# build sets
		launch_childprocess("e2buildsets.py --allparticles --setname all")

	# Import boxes
	if options.import_boxes:
		# Check to make sure there are micrographs
		if not os.access("info", os.R_OK):
			os.mkdir("info")
		# Do imports
		# we add boxsize/2 to the coords since box files are stored with origin being the lower left side of the box, but in EMAN2 origin is in the center
		if options.box_type == 'boxes':
			micros=os.listdir("micrographs")
			for filename in args:
				boxlist = []
				fh = open(filename, 'r')
				for line in fh.readlines():
					if line[0]=="#" : continue
					fields = line.split()
					if len(fields)<4 : continue		# skip lines that don't work
					boxlist.append([float(fields[0])+old_div(float(fields[3]),2), float(fields[1])+old_div(float(fields[3]),2), 'manual'])

				js=js_open_dict(info_name(filename,nodir=True))
				js["boxes"]=boxlist
				js.close()
				if not "{}.hdf".format(base_name(filename,nodir=True)) in micros:
					print("Warning: Imported boxes for {}, but micrographs/{}.hdf does not exist".format(base_name(filename),base_name(filename,True)))

		elif options.box_type == 'coords':
			micros=os.listdir("micrographs")
			for filename in args:
				boxlist = []
				fh = open(filename, 'r')
				for line in fh.readlines():
					if line[0]=="#" : continue
					fields = line.split()
					if len(fields)<2 : continue		# skip lines that don't work
					boxlist.append([float(fields[0]), float(fields[1]), 'manual'])
				js_open_dict(info_name(filename,nodir=True))["boxes"]=boxlist
				if not "{}.hdf".format(base_name(filename,nodir=True)) in micros:
					print("Warning: Imported boxes for {}, but micrographs/{}.hdf does not exist".format(base_name(filename),base_name(filename,True)))


		elif options.box_type == 'tiltedboxes':

			for filename in args:
				boxlist = []
				fh = open(filename, 'r')
				for line in fh.readlines():
					if line[0]=="#" : continue
					fields = line.split()
					if len(fields)<4 : continue		# skip lines that don't work
					boxlist.append([float(fields[0])+old_div(float(fields[3]),2), float(fields[1])+old_div(float(fields[3]),2), 'tilted'])
				js_open_dict(info_name(filename,nodir=True))["boxes_rct"]=boxlist

		elif options.box_type == 'untiltedboxes':
			for filename in args:
				boxlist = []
				fh = open(filename, 'r')
				for line in fh.readlines():
					if line[0]=="#" : continue
					fields = line.split()
					if len(fields)<4 : continue		# skip lines that don't work
					boxlist.append([float(fields[0])+old_div(float(fields[3]),2), float(fields[1])+old_div(float(fields[3]),2), 'untilted'])
				js_open_dict(info_name(filename,nodir=True))["boxes_rct"]=boxlist

		elif options.box_type == 'relion_star':
			bs = options.boxsize
			starfs = [f for f in args if '.star' in f]
			if len(starfs) < 1:
				print("You must specify at least one .star file containing particle coordinates")
				exit(1)
			for filename in starfs:
				print(("Importing from {}.star".format(base_name(filename,nodir=True))))
				sf = StarFile(filename)
				hdr = list(sf.keys())
				if len(hdr) < 3:
					print(("Could not parse {}".format(filename)))
					continue
				mk = "rlnMicrographName"
				yk = "rlnCoordinateY"
				xk = "rlnCoordinateX"
				project_micros = os.listdir('micrographs')
				if mk not in hdr or yk not in hdr or xk not in hdr:
					possible = "{}.hdf".format(base_name(filename.replace('_autopick.star',''),nodir=True))
					if possible in project_micros:
						micros = [possible]
					else:
						print(("{} does not follow the RELION header convention for single particle data. To use this program".format(filename)))
						if mk not in hdr: print("Micrograph names should be listed under _rlnMicrographName")
						if yk not in hdr: print("Y coordinates must be listed under _rlnCoordinateY")
						if xk not in hdr: print("X coordinates must be listed under _rlnCoordinateX")
						continue
				else: micros=[i.split('/')[-1] for i in np.unique(sf[mk])]
				if len(micros) == 1:
					mg = micros[0]
					boxlist = []
					print(("Found {} boxes for {}".format(len(sf[xk]),mg)))
					for x,y in zip(sf[xk],sf[yk]):
						xc = int(x)
						yc = int(y)
						boxlist.append([xc,yc,'manual']) # should probably be 'relion' or 'from_star'
					js_open_dict(info_name(mg,nodir=True))["boxes"]=boxlist
					if not "{}.hdf".format(base_name(mg,nodir=True)) in project_micros:
						print("Warning: Imported boxes for {}.hdf, but micrographs/{}.hdf does not exist".format(base_name(filename),base_name(mg,nodir=True)))
				elif len(micros) > 1:
					for mg in project_micros:
						boxlist = []
						ptcls = []
						for i,name in enumerate(sf[mk]):
							mgname = name.split('/')[-1].split('.')[0]
							#print(mgname,hdf_name,mg)
							if mg[:-4] in mgname: ptcls.append(i)
						print(("Found {} boxes for {}".format(len(ptcls),mg)))
						for p in ptcls:
							xc = int(sf[xk][p])
							yc = int(sf[yk][p])
							boxlist.append([xc,yc,'manual'])
						js_open_dict(info_name(mg,nodir=True))["boxes"]=boxlist
						if not "{}.hdf".format(base_name(mg,nodir=True)) in project_micros:
							print("Warning: Imported boxes for {}, but micrographs/{}.hdf does not exist".format(base_name(mg),base_name(mg,nodir=True)))

		else : print("ERROR: Unknown box_type")

	# Import particles
	if options.import_particles:
		if not os.access("particles", os.R_OK):
			os.mkdir("particles")

		fset=set([base_name(i) for i in args])
		if len(fset)!=len(args):
			print("ERROR: You specified multiple files to import with the same base name, eg - a10/abc123.spi and a12/abc123.spi. If you have multiple images with the same \
name, you will need to modify your naming convention (perhaps by prefixing the date) before importing. If the input files are in IMAGIC format, so you have .hed and .img files \
with the same name, you should specify only the .hed files (no renaming is necessary).")
			sys.exit(1)

		for i,fsp in enumerate(args):
			E2progress(logid,old_div(float(i),len(args)))
			if EMData(fsp,0,True)["nz"]>1 :
				run("e2proc2d.py {} particles/{}.hdf --threed2twod --inplace".format(fsp,base_name(fsp)))
			else: run("e2proc2d.py {} particles/{}.hdf --inplace".format(fsp,base_name(fsp)))

	if options.gainrefs != "" or options.darkrefs != "":
		refsdir = os.path.join(".","movierefs")
		if not os.access(refsdir, os.R_OK):
			os.mkdir(refsdir)

	if options.gainrefs != "":
		for ref in options.gainrefs.split(","):
			refname=os.path.join(refsdir,os.path.basename(ref))
			if refname[-4:] == ".mrc": refname+="s"
			if not os.path.isfile(refname):
				if options.importation == "move":
					os.rename(ref,refname)
				elif options.importation == "link":
					print("Movie references must be moved or copied. Linking is not supported.")
					sys.exit(1)
				elif options.importation == "copy":
					run("e2proc2d.py {} {} ".format(ref, refname))

	if options.darkrefs != "":
		for ref in options.darkrefs.split(","):
			refname=os.path.join(refsdir,os.path.basename(ref))
			if refname[-4:] == ".mrc": refname+="s"
			if not os.path.isfile(refname):
				if options.importation == "move":
					os.rename(ref,refname)
				elif options.importation == "link":
					print("Movie references must be moved or copied. Linking is not supported.")
					sys.exit(1)
				elif options.importation == "copy":
					run("e2proc2d.py {} {} ".format(ref, refname))


	if options.import_movies:

		moviesdir = os.path.join(".","movies")
		if not os.access(moviesdir, os.R_OK):
			os.mkdir(moviesdir)

		for filename in args:
			newname=os.path.join(moviesdir,os.path.basename(filename))
			if not os.path.isfile(newname):
				if newname[-4:] == ".mrc": newname+="s"
				if options.importation == "move":
					os.rename(filename,newname)
				if options.importation == "copy":
					run("e2proc2d.py {} {} ".format(filename, newname))
				if options.importation == "link":
					os.symlink(filename,newname)
		print("Done.")

	# Import tilts
	# if options.import_rawtilts:

	# 	stdir = os.path.join(".","raw_tilts")
	# 	if not os.access(stdir, os.R_OK):
	# 		os.mkdir("tilts")

	# 	for filename in args:
	# 		newname=os.path.join(stdir,os.path.basename(filename))
	# 		if options.importation == "move":
	# 			os.rename(filename,newname)
	# 		if options.importation == "copy":
	# 			tpos=filename.rfind('.')
	# 			if tpos>0: newname=os.path.join(stdir,os.path.basename(filename[:tpos]+'.hdf'))
	# 			else: newname=os.path.join(stdir,os.path.basename(filename))
	# 			cmd="e2proc2d.py {} {} ".format(filename, newname)
	# 			if options.invert: cmd+=" --mult -1 --process normalize "
	# 			#if options.tomoseg_auto:
	# 			#	cmd+=" --process filter.lowpass.gauss:cutoff_abs=.25 --process filter.highpass.gauss:cutoff_pixels=5 --process threshold.clampminmax.nsigma:nsigma=3 "
	# 			cmd+=options.preprocess
	# 			run(cmd)
	# 			print("Done.")
	# 		if options.importation == "link":
	# 			os.symlink(filename,newname)

	# Import tilt series
	if options.import_tiltseries:
		# try:
		# 	db=js_open_dict("info/project.json")
		# 	if options.apix == -1: 
		# 		options.apix = db["global.apix"]
		# 		print("Using global apix: {}".format(db["global.apix"]))
		# except: pass

		stdir = os.path.join(".","tiltseries")
		if not os.access(stdir, os.R_OK) : os.mkdir("tiltseries")

		for filename in args:
			newname=os.path.join(stdir,os.path.basename(filename))
			if options.importation == "move":
				os.rename(filename,newname)
			if options.importation == "copy":
				if os.path.isfile(newname): os.remove(newname)
				tpos=filename.rfind('.')
				if tpos>0: newname=os.path.join(stdir,os.path.basename(filename[:tpos]+'.hdf'))
				else: newname=os.path.join(stdir,os.path.basename(filename))
				cmd="e2proc2d.py {} {} --inplace ".format(filename, newname)
				if options.invert: cmd+=" --mult -1 --process normalize "
				if options.apix != -1: cmd += " --apix {} ".format(options.apix)
				#if options.tomoseg_auto:
				#	cmd+=" --process filter.lowpass.gauss:cutoff_abs=.25 --process filter.highpass.gauss:cutoff_pixels=5 --process threshold.clampminmax.nsigma:nsigma=3 "
				#cmd+=options.preprocess
				run(cmd)
				print("Done.")
			if options.importation == "link": os.symlink(filename,newname)
			if (options.importation == "link" or options.importation == "move") and options.apix != -1:
				run("e2proc3d.py {} {} --apix {} --threed2twod".format(newname, newname, options.apix))
	
	# Import tomograms
	if options.import_tomos:
		tomosdir = os.path.join(".","tomograms")
		if not os.access(tomosdir, os.R_OK):
			os.mkdir("tomograms")
		for filename in args:
			if options.importation == "move":
				os.rename(filename,os.path.join(tomosdir,os.path.basename(filename)))
			if options.importation == "copy":
				### use hdf file as output
				if options.shrink>1:
					shrinkstr="_bin{:d}".format(options.shrink)
				else:
					shrinkstr=""
				tpos=filename.rfind('.')
				if tpos>0:
					newname=os.path.join(tomosdir,os.path.basename(filename[:tpos]+shrinkstr+'.hdf'))
				else:
					newname=os.path.join(tomosdir,os.path.basename(filename))
				hdr=EMData(filename,0,True)
				cmd="e2proc3d.py {} {} ".format(filename, newname)
				
				# shrink, and clip from the origin to give us a good box size. Done from origin so segmentation results will be positioned well for scaling
				if options.shrink>1:
					nx=good_size(old_div(hdr["nx"],options.shrink))*options.shrink
					ny=good_size(old_div(hdr["ny"],options.shrink))*options.shrink
					nz=good_size(old_div(hdr["nz"],options.shrink))*options.shrink
#					cmd+="--clip {},{},{},{},{},{} --meanshrink {:d} ".format(nx,ny,nz,nx/2,ny/2,nz/2,options.shrink)
					cmd+="--clip {},{},{} --meanshrink {:d} ".format(nx,ny,nz,options.shrink)		# no origin shift. If result is scaled up and clipped to original size this may work better
				else:
					nx=good_size(hdr["nx"])
					ny=good_size(hdr["ny"])
					nz=good_size(hdr["nz"])
#					cmd+="--clip {},{},{},{},{},{} ".format(nx,ny,nz,nx/2,ny/2,nz/2,options.shrink)
					cmd+="--clip {},{},{} ".format(nx,ny,nz,options.shrink)
					
				if options.invert:
					cmd+=" --mult -1 --process normalize "
					
				cmd+=" --process normalize "
				#if options.tomoseg_auto:
				#	cmd+=" --process filter.lowpass.gauss:cutoff_abs=.25 --process filter.highpass.gauss:cutoff_pixels=5 --process normalize --process threshold.clampminmax.nsigma:nsigma=3 "
				#cmd+=options.preprocess
				run(cmd)
				print("Done.")
			#shutil.copy(filename,os.path.join(tomosdir,os.path.basename(filename)))
			if options.importation == "link":
				os.symlink(filename,os.path.join(tomosdir,os.path.basename(filename)))

	# # Import serialEM metadata
	# if options.serialem_mdoc:
	# 	for fn in args:
	# 		mdoc = read_mdoc(fn)
	# 		# check and correct project parameters from MDOC file contents
	# 		d = js_open_dict("info/project.json")
	# 		try: d.setval("global.apix",mdoc["PixelSpacing"],deferupdate=True)
	# 		except: pass
	# 		try: d.setval("global.microscope_voltage",mdoc["Voltage"],deferupdate=True)
	# 		except: pass
	# 		d.close()
	# 		# for each referenced image, append pertinent keys/values to corresponding info.json
	# 		for z in range(mdoc["zval"]+1):
	# 			tlt = mdoc[z]["SubFramePath"].rsplit("\\")[-1]+"_RawImages"
	# 			d = js_open_dict(info_name(tlt))
	# 			for k in mdoc[z].keys():
	# 				d.setval(k,mdoc[z][k],deferupdate=True)
	# 			d.close()
	#
	# # Import FEI metadata
	# if options.fei_tomo:
	# 	print("FEI tomography metadata not yet handled by this program.")
	# 	sys.exit(1)
	#
	# # Import UCSF tomo metadata
	# if options.ucsf_tomo:
	# 	print("UCSF tomography metadata not yet handled by this program.")
	# 	sys.exit(1)

	E2end(logid)

# def read_mdoc(mdoc):
# 	movie = {}
# 	frames = {}
# 	zval = -1
# 	frames[zval] = {}
# 	frames["misc"] = []
# 	frames["labels"] = []
# 	with open(mdoc) as mdocf:
# 		for l in mdocf.readlines():
# 			p = l.strip()
# 			if "ZValue" in p:
# 				zval+=1
# 				frames[zval] = {}
# 			elif p != "":
# 				x,y = p.split("=")[:2]
# 				x = x.strip()
# 				if x == "TiltAngle": frames[zval]["tilt_angle"]=float(y)
# 				elif x == "Magnification": frames[zval]["magnification"] = int(y)
# 				elif x == "Intensity": frames[zval]["intensity"]=float(y)
# 				elif x == "SpotSize": frames[zval]["spot_size"]=int(y)
# 				elif x == "Defocus": frames[zval]["defocus"]=float(y)
# 				elif x == "ExposureTime": frames[zval]["exposure_time"]=float(y)
# 				elif x == "Binning": frames[zval]["binning"]=float(y)
# 				elif x == "ExposureDose": frames[zval]["exposure_dose"]=float(y)
# 				elif x == "RotationAngle": frames[zval]["rotation_angle"]=float(y)
# 				elif x == "StageZ": frames[zval]["stage_z"]=float(y)
# 				elif x == "CameraIndex": frames[zval]["camera_index"]=int(y)
# 				elif x == "DividedBy2": frames[zval]["divide_by_2"]=int(y)
# 				elif x == "MagIndex": frames[zval]["mag_index"]=int(y)
# 				elif x == "TargetDefocus": frames[zval]["target_defocus"]=float(y)
# 				elif x == "ImageShift": frames[zval]["image_shift"]=map(float,y.split())
# 				elif x == "StagePosition": frames[zval]["stage_position"]=map(float,y.split())
# 				elif x == "MinMaxMean":
# 					vals = map(float,y.split())
# 					frames[zval]["min"]=vals[0]
# 					frames[zval]["max"]=vals[1]
# 					frames[zval]["mean"]=vals[2]
# 				elif x == "SubFramePath":
# 					sfp = base_name(y).split("-")[-1]
# 					frames[zval]["SubFramePath"]=sfp
# 				elif x == "DateTime": frames[zval]["date_time"]=y
# 				elif x == "PixelSpacing": frames["global.apix"] = float(y)
# 				elif x == "Voltage": frames["global.microscope_voltage"] = float(y)
# 				elif x == "ImageFile": frames["image_file"] = str(y)
# 				elif x == "ImageSize": frames["image_size"] = y.split()
# 				elif x == "DataMode": frames["data_mode"] = y
# 				elif x == "PriorRecordDose": frames["prior_record_dose"] = y
# 				elif x == "FrameDosesAndNumber": frames["frame_doses_and_number"] = y
# 				elif x == "[T": frames["labels"].append(y.replace("]",""))
# 				elif "PreexposureTime" in x: frames[zval]["preexposure_time"] = float(y)
# 				elif "TotalNumberOfFrames" in x: frames[zval]["frame_count"] = int(y)
# 				elif "FramesPerSecond" in x: frames[zval]["frames_per_second"] = float(y)
# 				elif "ProtectionCoverMode" in x: frames[zval]["protection_cover_mode"] = y
# 				elif "ProtectionCoverOpenDelay" in x: frames[zval]["protection_cover_open_delay"] = float(y)
# 				elif "TemperatureDetector" in x: frames[zval]["detector_temperature"] = float(y)
# 				elif "FaradayPlatePeakReading" in x: frames[zval]["faraday_plate_peak_reading"] = float(y)
# 				elif "SensorModuleSerialNumber" in x: frames[zval]["sensor_module_serial_number"] = y
# 				elif "ServerSoftwareVersion" in x: frames[zval]["server_software_version"] = y
# 				elif "SensorReadoutDelay" in x: frames[zval]["sensor_readout_delay"] = y
# 				else: frames["misc"].append(y) # catches any missed parameters
#
# 	frames["zval"] = zval
# 	return frames

def run(command):
	"Mostly here for debugging, allows you to control how commands are executed (os.system is normal)"

	print("{}: {}".format(time.ctime(time.time()),command))
	ret=launch_childprocess(command)

	# We put the exit here since this is what we'd do in every case anyway. Saves replication of error detection code above.
	if ret !=0 :
		print("Error running: ",command)
		sys.exit(1)

	return

if __name__ == "__main__":
	main()