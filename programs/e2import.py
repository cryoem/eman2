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
	usage = """prog [options] <micrgrpah1, microgaph2....
	Import particles coordinats from box files. To work the box file name must be the same as the microgrpah name save the extension.
	>"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_pos_argument(name="import_files",help="List the files to import here.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=True)",  row=0, col=0, rowspan=1, colspan=2, nosharedb=True, mode='coords,parts,tomos')
	parser.add_header(name="filterheader", help='Options below this label are specific to e2import', title="### e2import options ###", row=1, col=0, rowspan=1, colspan=2, mode='coords,parts,tomos')
	parser.add_argument("--import_particles",action="store_true",help="Import particles",default=False, guitype='boolbox', row=2, col=0, rowspan=1, colspan=1, mode='parts[True]')
	parser.add_argument("--import_tomos",action="store_true",help="Import tomograms",default=False, guitype='boolbox', row=2, col=0, rowspan=1, colspan=1, mode='tomos[True]')
	parser.add_argument("--importation",help="import particles",default='copy',guitype='combobox',choicelist='["move","copy","link"]',row=2,col=1,rowspan=1,colspan=1, mode='tomos')
	parser.add_argument("--import_boxes",action="store_true",help="Import boxes",default=False, guitype='boolbox', row=2, col=0, rowspan=1, colspan=1, mode='coords[True]')
	parser.add_argument("--extension",type=str,help="Extension of the micrographs that the boxes match", default='dm3', guitype='strbox',row=3, col=0, rowspan=1, colspan=1, mode='coords')
	parser.add_argument("--box_type",help="Type of boxes to import, normally boxes, but for tilted data use tiltedboxes, and untiltedboxes for the tilted  particle partner",default=None,guitype='combobox',choicelist='["boxes","tiltedboxes","untiltedboxes"]',row=2,col=1,rowspan=1,colspan=1, mode="coords['boxes']")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv,options.ppid)
	
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
					fields = line.split()
					boxlist.append([float(fields[0])+float(fields[3])/2, float(fields[1])+float(fields[3])/2, 'manual'])
					
				js_open_dict(info_name(filename))["boxes"]=boxlist
				if not "{}.hdf".format(base_name(filename)) in micros:
					print "Warning: Imported boxes for {}, but micrographs/{}.hdf does not exist".format(base_name(filename),base_name(filename))
				
		elif options.box_type == 'tiltedboxes':
			
			for filename in args:
				boxlist = []
				fh = open(filename, 'r')
				for line in fh.readlines():
					fields = line.split()
					boxlist.append([float(fields[0])+float(fields[3])/2, float(fields[1])+float(fields[3])/2, 'tilted'])
				js_open_dict(info_name(filename))["boxes_tilted"]=boxlist
				
		elif options.box_type == 'untiltedboxes':		
			for filename in args:
				boxlist = []
				fh = open(filename, 'r')
				for line in fh.readlines():
					fields = line.split()
					boxlist.append([float(fields[0])+float(fields[3])/2, float(fields[1])+float(fields[3])/2, 'untilted'])
				js_open_dict(info_name(filename))["boxes_untilted"]=boxlist

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
			
		for fsp in args:
			if EMData(fsp,0,True)["nz"]>1 :
				run("e2proc2d.py {} particles/{}.hdf --threed2twod".format(fsp,base_name(fsp)))
			else: run("e2proc2d.py {} particles/{}.hdf".format(fsp,base_name(fsp)))
	   
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
