#!/usr/bin/env python

#
# Author: Steve Ludtke 03/04/2017 (sludtke@bcm.edu)
# Copyright (c) 2014- Baylor College of Medicine
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
#
#


from EMAN2 import *
from optparse import OptionParser
from math import *
import os
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <SPARX hdf stack>

	This program will read a SPARX particle HDF stack with embedded CTF parameters and produce an EMAN2 project. 

"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	#options associated with e2refine.py
	parser.add_header(name="text1", help='Important instructions', title="Use this to create an EMAN2.2 project from a SPARX stack ", row=0, col=0, rowspan=1, colspan=3)
	parser.add_pos_argument(name="sparx_file",help="Select SPARX file", default="", guitype='filebox', browser="EMParticlesEditTable(withmodal=True,multiselect=False)",  row=6, col=0,rowspan=1, colspan=3)
#	parser.add_argument("--refinedefocus",  action="store_true", help="Will refit defocus to +-0.1 micron then further optimize using SSNR",default=False,guitype='boolbox', row=10, col=2, rowspan=1, colspan=1,mode="[True]")
#	parser.add_argument("--refitdefocus",  action="store_true", help="Will use EMAN2 CTF fitting to refit the defocus values within +-0.1 micron, astigmatism unchanged",default=False,guitype='boolbox', row=10, col=1, rowspan=1, colspan=1)
#	parser.add_argument("--fixeddefocus",  action="store_true", help="Will use defocus values from original HDF file.",default=False,guitype='boolbox', row=10, col=1, rowspan=1, colspan=1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()


	logid=E2init(sys.argv,options.ppid)

	try: os.mkdir("eman2")
	except: pass
	try: os.mkdir("eman2/particles")
	except: pass
	os.chdir("eman2")	# many things need to happen with the project directory as a base

	args[0]="../"+args[0]
	img=EMData(args[0],0)
	ctf=img["ctf"]

	prj=js_open_dict("info/project.json")
	try:
		prj["global.apix"]=ctf.apix
		prj["global.microscope_cs"]=ctf.cs
		prj["global.microscope_voltage"]=ctf.voltage
	except:
		print "No ctf info found. This shouldn't happen with a SPARX input file."
		sys.exit(1)
		
	print "Processing ",args[0]

	N=EMUtil.get_image_count(args[0])
	olddf=0
	micronum=0		# number of micrograph
	for i in xrange(N):
		img=EMData(args[0],i)
		ctf=img["ctf"]
		if ctf.defocus!=olddf :
			micronum+=1
			microname="particles/img{:05d}.hdf".format(micronum)
			jdb=js_open_dict(info_name(microname))
			jdb["ctf_frame"]=[512,ctf,(256,256),tuple(),5,1]
			olddf=ctf.defocus
			if options.verbose : print "{}) ({}/{}) defocus={}".format(micronum,i,N,ctf.defocus)
		img.del_attr("ctf")
		img.write_image(microname,-1)

	print micronum," micrographs found"


#	print "Defocus processing"

	#if options.refinedefocus : 
		#dfopt="--curdefocushint --refinebysnr"
		#if options.verbose>0 : print "Defocus Refit and SSNR Refinement to +-0.1 micron"
	#elif options.refitdefocus : 
		#dfopt="--curdefocushint"
		#if options.verbose>0 : print "Defocus Refit to +-0.1 micron from Relion values"
	#elif options.fixeddefocus: 
		#dfopt="--curdefocusfix"
		#if options.verbose>0 : print "Computing particle SNRs"
	#else:
		#dfopt=" "
		#if options.verbose>0 : print "Full refit of CTF parameters"

#	launch_childprocess("e2ctf_auto.py --voltage {} --cs {} --ac {} --apix {} --allparticles --autofit {} -v {}".format(ctf.voltage,ctf.cs,ctf.ampcont,ctf.apix,dfopt,options.verbose-1))
	
	print "Done. Project in eman2/"
	print "Please run e2ctf_auto.py to complete CTF parameter determination"
	E2end(logid)

if __name__ == "__main__":
    main()
