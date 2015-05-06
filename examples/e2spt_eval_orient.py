#!/usr/bin/env python
#
# Author: Steven Ludtke, March 2015
# Copyright (c) 2015- Baylor College of Medicine
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

from EMAN2 import *
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <json file>
This program reads a JSON file produced by e2spt_classaverage and analyzes the orientation paramters
too look for issues with preferred orientation, etc.
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
#	parser.add_argument("--input", default='',type=str, help="The name of the hdf stack of volumes to process.")
	parser.add_argument("--output", default="orient.txt",type=str, help="The name of the output average volume.")
	
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness.")

	(options, args) = parser.parse_args()

	try: db=js_open_dict(args[0])
	except:
		print "ERROR: could not open ",args[0]
		sys.exit(1)

	logid=E2init(sys.argv,options.ppid)

	out=file(options.output,"w")

	alts=[0]*9
	for k in db.keys():
		xf=db[k][0]
		xfd=xf.get_params("eman")
		out.write("%1.3f,\t%1.3f,\t%1.3f,\t%1.3g\n"%(xfd["az"],xfd["alt"],xfd["phi"],float(db[k][1])))
		
		try:
			an=int(floor(xfd["alt"]/20.0001))
			alts[an]+=1
		except: pass
	

	print "Altitude distribution:"

	for i in xrange(9) :
		print "%d - %d: %d"%(i*20.0,(i+1)*20.0,alts[i])


	print "See also: ",options.output
	E2end(logid)
	

if __name__ == '__main__':
	main()
