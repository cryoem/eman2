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
import re, os
from EMAN2 import *

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <micrgrpah1, microgaph2....>
	Use this program to import and filter raw micrographs. If you choose to filter and/or convert format, this program will process each micrograph
	and dump them into the directory './micrographs.', otherwise the micrographs will simply be moved into './micrographs'. If you select the option
	--moverawdata AND you filter or change format, your original micrographs will be moved into the directory './raw_micrographs' and your 
	filtered micrographs will be in './micrographs as usual. BDB files are not moved, but they can be processed."""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_pos_argument(name="micrographs",help="List the micrographs to filter here.", default="", guitype='filebox', browser="EMRawDataTable(withmodal=True,multiselect=True)", positional=True, row=0, col=0,rowspan=1, colspan=2, mode='filter')
	parser.add_pos_argument(name="import_files",help="List the files to import/filter here.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=True)", positional=True, row=0, col=0,rowspan=1, colspan=2, mode='import')
	parser.add_header(name="filterheader", help='Options below this label are specific to filtering', title="### filtering options ###", row=1, col=0, rowspan=1, colspan=2, mode='import,filter')
	parser.add_argument("--invert",action="store_true",help="Invert contrast",default=False, guitype='boolbox', row=2, col=0, rowspan=1, colspan=1, mode='filter[True]')
	parser.add_argument("--edgenorm",action="store_true",help="Edge normalize",default=False, guitype='boolbox', row=2, col=1, rowspan=1, colspan=1, mode='filter[True]')
	parser.add_argument("--xraypixel",action="store_true",help="Filter X-ray pixels",default=False, guitype='boolbox', row=3, col=0, rowspan=1, colspan=1, mode='filter[True]')
	parser.add_argument("--inplace",action="store_true",help="Do processing inplace",default=False, guitype='boolbox', row=4, col=0, rowspan=1, colspan=1)
	parser.add_argument("--moverawdata",action="store_true",help="Move raw data to directory ./raw_micrographs after filtration",default=False, guitype='boolbox', row=3, col=1, rowspan=1, colspan=1, mode='filter[True]')
	parser.add_argument("--suffix",type=str,help="Filename suffix",default="_filt",guitype='strbox', row=5, col=0, rowspan=1, colspan=2, mode='filter')
	parser.add_argument("--format", help="Format of the output particles images, should be bdb,img,spi or hdf", default="hdf", guitype='combobox', choicelist="['bdb','hdf','img','spi']", row=6, col=0, rowspan=1, colspan=2, mode="filter")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	(options, args) = parser.parse_args()
	
	microdir = os.path.join(".","micrographs")
	if not os.access(microdir, os.R_OK):
		os.mkdir("micrographs")
	
	# If we do not do any filtering and keep the micorgpahs in the same format, only need to move micros.
	if not options.invert and not options.edgenorm and not options.xraypixel and (os.path.splitext(arg[0])[1] == ("."+options.format)):
		for arg in args:
			os.rename(arg,os.path.join(microdir,os.path.basename(arg)))
		print "exiting...."
		exit(0)
	
	
	logid=E2init(sys.argv,options.ppid)
	
	# After filtration we move micrographs to a directory 'raw_micrographs', if desired
	if options.moverawdata:
		originalsdir = os.path.join(".","raw_micrographs")
		if not os.access(originalsdir, os.R_OK):
			os.mkdir("raw_micrographs")
	for i,arg in enumerate(args):
		launch_childprocess(get_proc2dcmd(options, arg))
		if options.moverawdata: 
			if arg[:4].upper() != "BDB:": os.rename(arg,os.path.join(originalsdir,os.path.basename(arg)))
		E2progress(logid,(float(i)/float(len(args))))
		
	E2end(logid)
		
def get_proc2dcmd(options, filename):
	""" Return a proc2d command """
	if not options.inplace:
		base = os.path.splitext(os.path.basename(filename))[0]
		if options.format == 'bdb': 
			output = "bdb:micrographs#"+base
		else:
			output = os.path.join(os.path.join(".","micrographs"),base+options.suffix+"."+options.format)
	
	if options.inplace:
		cmd = "e2proc2d.py %s"%filename
	else:
		cmd = "e2proc2d.py %s %s"%(filename,output)
		
	if options.invert: cmd += " --mult=-1"
	if options.edgenorm: cmd += " --process=normalize.edgemean"
	if options.xraypixel: cmd += " --process=threshold.clampminmax.nsigma:nsigma=4:tomean=1"
	if options.inplace: cmd += " --inplace"
	
	return cmd
	
if __name__ == "__main__":
	main()