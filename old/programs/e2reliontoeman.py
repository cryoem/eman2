#!/usr/bin/env python

#
# Author: Steve Ludtke 04/16/14 (sludtke@bcm.edu)
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
import time
import traceback
from EMAN2star import StarFile


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <relion STAR file>

This program will take data from a Relion project and convert it into a basic EMAN2 project. Provide 
the name of the Relion STAR file associated with the raw particle data. An eman2 subdirectory will be
created, and the images, and available metadata will be copied into the new project. CTF parameters
will be extracted from the STAR file and will be automatically processed through EMAN2's CTF procedure
(this takes most of the time the script will require).

"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	#options associated with e2refine.py
	parser.add_header(name="text1", help='Important instructions', title="Use this to create an EMAN2.1 project from a Relion project:", row=0, col=0, rowspan=1, colspan=3)
	parser.add_header(name="text2", help='Important instructions', title="* cd <folder with Relion STAR file>", row=1, col=0, rowspan=1, colspan=3)
	parser.add_header(name="text3", help='Important instructions', title="* run e2projectmanager, and use this tool", row=2, col=0, rowspan=1, colspan=3)
	parser.add_header(name="text4", help='Important instructions', title="* exit PM, cd eman2, run PM from new eman2 folder", row=3, col=0, rowspan=1, colspan=3)
	parser.add_pos_argument(name="star_file",help="Select STAR file", default="", guitype='filebox', browser="EMParticlesEditTable(withmodal=True,multiselect=False)",  row=6, col=0,rowspan=1, colspan=3)
	parser.add_argument("--apix", default=0, type=float,help="The angstrom per pixel of the input particles.", guitype='floatbox', row=8, col=0, rowspan=1, colspan=1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--fixeddefocus",  action="store_true", help="Defocus and astigmatism are used unchanged from STAR file",default=False,guitype='boolbox', row=10, col=0, rowspan=1, colspan=1)
	parser.add_argument("--refinedefocus",  action="store_true", help="Will refit defocus to +-0.1 micron then further optimize using SSNR",default=False,guitype='boolbox', row=10, col=2, rowspan=1, colspan=1,mode="[True]")
	parser.add_argument("--refitdefocus",  action="store_true", help="Will use EMAN2 CTF fitting to refit the defocus values within +-0.1 micron, astigmatism unchanged",default=False,guitype='boolbox', row=10, col=1, rowspan=1, colspan=1)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()


	logid=E2init(sys.argv,options.ppid)

	try: os.mkdir("eman2")
	except: pass
	try: os.mkdir("eman2/particles")
	except: pass
	os.chdir("eman2")	# many things need to happen with the project directory as a base

	if options.verbose>0 : print "Parsing STAR file"
	star=StarFile("../"+args[0])

	if options.apix<=0 :
		try:
			options.apix=star["rlnDetectorPixelSize"][0]/star["rlnMagnification"][0]*10000.0
			print "Using {} A/pix from Relion file".format(options.apix)
		except:
			print "A/pix not specified and not found in STAR file"
			sys.exit(1)

	prj=js_open_dict("info/project.json")
	try:
		prj["global.apix"]=options.apix
		prj["global.microscope_cs"]=star["rlnSphericalAberration"][0]
		if prj["global.microscope_cs"]<=0.0 : prj["global.microscope_cs"]=0.001
		prj["global.microscope_voltage"]=star["rlnVoltage"][0]
		print "V={} Cs={}".format(prj["global.microscope_voltage"],prj["global.microscope_cs"])
	except:
		print "Did not find Voltage and Cs in Relion file"
		
	oldname=""
	olddf=-1.0
	micronum=0		# number of micrograph
	fnum=0			# image number in new file
	for i in xrange(len(star["rlnImageName"])):
		name=star["rlnImageName"][i].split("@")[1]
		imgnum=int(star["rlnImageName"][i].split("@")[0])-1
	
		if name!=oldname:
			hdr=EMData("../"+name,0,True)
			nx=hdr["nx"]
			ny=hdr["ny"]
			oldname=name
			if options.verbose>0 : print "Particle dimensions: {}x{}".format(nx,ny)
		
		if i==0 or star["rlnDefocusU"][i-1]!=star["rlnDefocusU"][i]:
			if micronum>0 and options.verbose>0 : print "Image {}: {} particles processed, df={}".format(micronum,fnum,ctf.defocus)
			micronum+=1
			fnum=0
			microname="particles/{}_{:04d}.hdf".format(base_name(name),micronum)
			jdb=js_open_dict(info_name(microname))
			
			# Make a "micrograph" CTF entry for each set of different defocuses to use when fitting
			dfu=star["rlnDefocusU"][i]
			dfv=star["rlnDefocusV"][i]
			dfang=star["rlnDefocusAngle"][i]
			ctf=EMAN2Ctf()
			ctf.from_dict({"defocus":(dfu+dfv)/20000.0,"dfang":dfang,"dfdiff":(dfu-dfv)/10000.0,"voltage":star["rlnVoltage"][i],"cs":max(star["rlnSphericalAberration"][i],0.0001),"ampcont":star["rlnAmplitudeContrast"][i]*100.0,"apix":options.apix})
			jdb["ctf_frame"]=[512,ctf,(256,256),tuple(),5,1]
		
		# copy the image
		if name[-5:]==".mrcs" : img=EMData("../"+name,imgnum)		# read one slice from the MRC stack
		else: img=EMData("../"+name,0,False,Region(0,0,imgnum,nx,ny,1))		# read one slice from the MRC stack
		img.write_image(microname,fnum)
		fnum+=1

	if options.refinedefocus : 
		dfopt="--curdefocushint --refinebysnr"
		if options.verbose>0 : print "Defocus Refit and SSNR Refinement to +-0.1 micron"
	elif options.refitdefocus : 
		dfopt="--curdefocushint"
		if options.verbose>0 : print "Defocus Refit to +-0.1 micron from Relion values"
	elif options.fixeddefocus: 
		dfopt="--curdefocusfix"
		if options.verbose>0 : print "Computing particle SNRs"
	else:
		dfopt=" "
		if options.verbose>0 : print "Full refit of CTF parameters"

	launch_childprocess("e2ctf.py --voltage {} --cs {} --ac {} --apix {} --allparticles --autofit {} -v {}".format(ctf.voltage,ctf.cs,ctf.ampcont,ctf.apix,dfopt,options.verbose-1))
	
	E2end(logid)

if __name__ == "__main__":
    main()
