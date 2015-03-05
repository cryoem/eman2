#!/usr/bin/env python
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
#
import os, shutil, glob
from EMAN2 import *

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] files
	This program performs a variety of tasks for getting data or metadata from other programs into an EMAN2 project. 
	
	import_particles - will simply copy a set of per-micrograph particle files into EMAN2.1's preferred HDF format in particles/
	import_boxes - will read EMAN1 '.box' files (text files containing coordinates) into appropriate info/*json files (see --box_type)
	import_tomos - imports subtomogams for a SPT project (see also --importation)
	import_eman1 - will convert a typical EMAN1 phase-flipped start.hed/img file into an EMAN2 project (converting files, fixing CTF, splitting, ...)
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_pos_argument(name="import_files",help="List the files to import here.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=True)",  row=0, col=0, rowspan=1, colspan=2, nosharedb=True, mode='coords,parts,tomos,eman1')
	parser.add_header(name="filterheader", help='Options below this label are specific to e2import', title="### e2import options ###", row=1, col=0, rowspan=1, colspan=2, mode='coords,parts,tomos')
	parser.add_argument("--import_particles",action="store_true",help="Import particles",default=False, guitype='boolbox', row=2, col=0, rowspan=1, colspan=1, mode='parts[True]')
	parser.add_argument("--import_eman1",action="store_true",help="This will import a phase-flipped particle stack from EMAN1",default=False, guitype='boolbox', row=2, col=0, rowspan=1, colspan=1, mode='eman1[True]')
	parser.add_argument("--import_tomos",action="store_true",help="Import tomograms",default=False, guitype='boolbox', row=2, col=0, rowspan=1, colspan=1, mode='tomos[True]')
	parser.add_argument("--importation",help="Specify mode move, copy or link, for importing tomograms only",default='copy',guitype='combobox',choicelist='["move","copy","link"]',row=2,col=1,rowspan=1,colspan=1, mode='tomos')
	parser.add_argument("--import_boxes",action="store_true",help="Import boxes",default=False, guitype='boolbox', row=2, col=0, rowspan=1, colspan=1, mode='coords[True]')
	parser.add_argument("--extension",type=str,help="Extension of the micrographs that the boxes match", default='dm3')
	parser.add_argument("--box_type",help="Type of boxes to import, normally boxes, but for tilted data use tiltedboxes, and untiltedboxes for the tilted  particle partner",default="boxes",guitype='combobox',choicelist='["boxes","coords","tiltedboxes","untiltedboxes"]',row=2,col=1,rowspan=1,colspan=1, mode="coords['boxes']")
	parser.add_argument("--curdefocushint",action="store_true",help="Used with import_eman1, will use EMAN1 defocus as starting point",default=False, guitype='boolbox', row=4, col=0, rowspan=1, colspan=1, mode='eman1[True]')
	parser.add_argument("--curdefocusfix",action="store_true",help="Used with import_eman1, will use EMAN1 defocus unchanged (+-.001 um)",default=False, guitype='boolbox', row=4, col=1, rowspan=1, colspan=1, mode='eman1[False]')
	parser.add_argument("--threads", default=1,type=int,help="Number of threads to run in parallel on a single computer when multi-computer parallelism isn't useful",guitype='intbox', row=6, col=0, rowspan=1, colspan=1, mode='eman1[1]')
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv,options.ppid)

	# Import EMAN1
	# will read start.hed/img, split by micrograph (based on defocus), and reprocess CTF in EMAN2 style
	if options.import_eman1 :
		try:
			n=EMUtil.get_image_count(args[0])
		except:
			print "Error, couldn't read images from: ",args[0]
			sys.exit(1)
		
		try:
			img=EMData(args[0],0)
			ctf=img["ctf"]
		except:
			print "Error, start.hed/img must be phase-flipped to import"
			sys.exit(1)
		
		db=js_open_dict("info/project.json")
		db["global.apix"]=ctf.apix
		db["global.cs"]=ctf.cs
		db["global.voltage"]=ctf.voltage
		
		try: os.mkdir("particles")
		except: pass
		
		imgnum=0
		lastdf=-1.0
		for i in xrange(n):
			img=EMData(args[0],i)
			ctf=img["ctf"]
			img.del_attr("ctf")
			fft1=img.do_fft()
			if ctf.defocus!=lastdf :
				imgnum+=1
				if options.verbose>0: print "Defocus {:4.2f} particles{:03d}".format(ctf.defocus,imgnum)
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
					boxlist.append([float(fields[0])+float(fields[3])/2, float(fields[1])+float(fields[3])/2, 'manual'])

				js_open_dict(info_name(filename,nodir=True))["boxes"]=boxlist
				if not "{}.hdf".format(base_name(filename,nodir=True)) in micros:
					print "Warning: Imported boxes for {}, but micrographs/{}.hdf does not exist".format(base_name(filename),base_name(filename,True))

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
					print "Warning: Imported boxes for {}, but micrographs/{}.hdf does not exist".format(base_name(filename),base_name(filename,True))


		elif options.box_type == 'tiltedboxes':

			for filename in args:
				boxlist = []
				fh = open(filename, 'r')
				for line in fh.readlines():
					if line[0]=="#" : continue
					fields = line.split()
					if len(fields)<4 : continue		# skip lines that don't work
					boxlist.append([float(fields[0])+float(fields[3])/2, float(fields[1])+float(fields[3])/2, 'tilted'])
				js_open_dict(info_name(filename,nodir=True))["boxes_rct"]=boxlist

		elif options.box_type == 'untiltedboxes':
			for filename in args:
				boxlist = []
				fh = open(filename, 'r')
				for line in fh.readlines():
					if line[0]=="#" : continue
					fields = line.split()
					if len(fields)<4 : continue		# skip lines that don't work
					boxlist.append([float(fields[0])+float(fields[3])/2, float(fields[1])+float(fields[3])/2, 'untilted'])
				js_open_dict(info_name(filename,nodir=True))["boxes_rct"]=boxlist

		else : print "ERROR: Unknown box_type"

	# Import particles
	if options.import_particles:
		if not os.access("particles", os.R_OK):
			os.mkdir("particles")

		fset=set([base_name(i) for i in args])
		if len(fset)!=len(args):
			print "ERROR: You specified multiple files to import with the same base name, eg - a10/abc123.spi and a12/abc123.spi. If you have multiple images with the same \
name, you will need to modify your naming convention (perhaps by prefixing the date) before importing. If the input files are in IMAGIC format, so you have .hed and .img files \
with the same name, you should specify only the .hed files (no renaming is necessary)."
			sys.exit(1)

		for i,fsp in enumerate(args):
			E2progress(logid,float(i)/len(args))
			if EMData(fsp,0,True)["nz"]>1 :
				run("e2proc2d.py {} particles/{}.hdf --threed2twod --inplace".format(fsp,base_name(fsp)))
			else: run("e2proc2d.py {} particles/{}.hdf --inplace".format(fsp,base_name(fsp)))

	# Import tomograms
	if options.import_tomos:
		tomosdir = os.path.join(".","rawtomograms")
		if not os.access(tomosdir, os.R_OK):
			os.mkdir("rawtomograms")
		for filename in args:
			if options.importation == "move":
				os.rename(filename,os.path.join(tomosdir,os.path.basename(filename)))
			if options.importation == "copy":
				shutil.copy(filename,os.path.join(tomosdir,os.path.basename(filename)))
			if options.importation == "link":
				os.symlink(filename,os.path.join(tomosdir,os.path.basename(filename)))
	E2end(logid)

def run(command):
	"Mostly here for debugging, allows you to control how commands are executed (os.system is normal)"

	print "{}: {}".format(time.ctime(time.time()),command)
	ret=launch_childprocess(command)

	# We put the exit here since this is what we'd do in every case anyway. Saves replication of error detection code above.
	if ret !=0 :
		print "Error running: ",command
		sys.exit(1)

	return

if __name__ == "__main__":
	main()
