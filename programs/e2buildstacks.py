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
import os, re
from EMAN2 import *

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <micrgrpah1, microgaph2....>
	Use e2bdb.py or e2proc2d.py to build particle sets
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_pos_argument(name="stack_files",help="List building material (sets) here.", default="", guitype='filebox', browser="EMParticlesEditTable(withmodal=True,multiselect=True)",  row=0, col=0,rowspan=1, colspan=2, nosharedb=True)
	parser.add_header(name="buildheader", help='Options below this label are specific to e2buildstacks', title="### e2buildstacks options ###", row=1, col=0, rowspan=1, colspan=2)
	parser.add_argument("--stackname",type=str,help="Name of the stack to build", default='my_stack', guitype='strbox',row=2, col=0, rowspan=1, colspan=1)
	parser.add_argument("--filetype",help="Type of file",default='bdb',guitype='combobox',choicelist='["bdb","hdf","spi"]',row=3,col=0,rowspan=1,colspan=1)
	parser.add_argument("--excludebad",action="store_true",help="Exculde bad particles. Only works for BDB stacks",default=True, guitype='boolbox',row=4,col=0,rowspan=1,colspan=1)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	(options, args) = parser.parse_args()
	
	# Now do the work
	path = "sets"

	boxesdir = os.path.join(".",path)
	if not os.access(boxesdir, os.R_OK):
		os.mkdir(path)
	
	#If we use bdb, then make a virtual stack (but only if all input files are also BDB)
	out_name = ""
	if options.filetype == "bdb":
		out_name = "bdb:%s#%s"%(path,os.path.splitext(options.stackname)[0])
	
		# Check to see if any of the input files are not BDB format if this is so then we Cannot make a virtual stack and must make a 'regular' stack
		bdbfiles = True
		for arg in args:
			if arg[:4] != "bdb:": 
				bdbfiles = False
				break
			
		# Make a virtual stack. This is only possible if all input files and the output format are BDB. Otherwise pass over this code and make a 'real' stack
		if bdbfiles:		
			# Check exclude DB
			if db_check_dict("bdb:select"):
				select_db = db_open_dict("bdb:select")
			else:
				select_db = None
			cmd = ""
	
			for i,name in enumerate(args):
				cmd = "e2bdb.py"
				no_exc = True
				if options.excludebad and select_db != None:
					exc_db = re.sub(r'_ctf_flip$|_ctf_wiener$','',get_file_tag(name))
					if select_db.has_key(exc_db):
						cmd += " "+name+"?exclude."+exc_db
						no_exc=False
				
				if no_exc:
					cmd += " "+name
				
				cmd += " --appendvstack="+out_name
		
				success = (os.system(cmd) in (0,11,12))
			
			exit(0)
	
	#Make a 'real' stack
	if options.filetype == "hdf":
		out_name = os.path.join(path,"%s.hdf"%os.path.splitext(options.stackname)[0])
	if options.filetype == "spi":
		out_name = os.path.join(path,"%s.spi"%os.path.splitext(options.stackname)[0])
		
	for arg in args:
		launch_childprocess("e2proc2d.py %s %s"%(arg,out_name))
	
if __name__ == "__main__":
	main()