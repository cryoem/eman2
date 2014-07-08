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
	parser.add_argument("--apix", default=0, type=float,help="The angstrom per pixel of the input particles. This argument is required if you specify the --mass argument. If unspecified (set to 0), the convergence plot is generated using either the project apix, or if not an apix of 1.", guitype='floatbox', row=16, col=0, rowspan=1, colspan=1, mode="refinement['self.pm().getAPIX()']")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--refinedefocus",  action="store_true", help="Will use EMAN2 CTF fitting to refine the defocus by SNR optimization (+-0.1 micron from the current values, no astigmatism adjustment)")
	parser.add_argument("--refitdefocus",  action="store_true", help="Will use EMAN2 CTF fitting to refit the defocus values (+-0.1 micron, astigmatism unchanged)")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	if options.apix<=0 :
		print "A/pix must be specified"
		sys.exit(1)

	logid=E2init(sys.argv,options.ppid)

	try: os.mkdir("eman2")
	except: pass
	try: os.mkdir("eman2/particles")
	except: pass
	os.chdir("eman2")	# many things need to happen with the project directory as a base

	if options.verbose>0 : print "Parsing STAR file"
	star=StarFile("../"+args[0])
		
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
			ctf.from_dict({"defocus":(dfu+dfv)/20000.0,"dfang":dfang,"dfdiff":(dfu-dfv)/10000.0,"voltage":star["rlnVoltage"][i],"cs":star["rlnSphericalAberration"][i],"ampcont":star["rlnAmplitudeContrast"][i]*100.0,"apix":options.apix})
			jdb["ctf_frame"]=[512,ctf,(256,256),tuple(),5,1]
		
		# copy the image
		if name[-5:]==".mrcs" : img=EMData("../"+name,imgnum)		# read one slice from the MRC stack
		else: img=EMData("../"+name,0,False,Region(0,0,imgnum,nx,ny,1))		# read one slice from the MRC stack
		img.write_image(microname,fnum)
		fnum+=1

	if options.refinedefocus : 
		dfopt="--curdefocushint --refinebysnr"
		if options.verbose>0 : print "CTF Refinement"
	elif options.refitdefocus : 
		dfopt="--curdefocushint"
		if options.verbose>0 : print "CTF Refit"
	else: 
		dfopt="--curdefocusfix"
		if options.verbose>0 : print "Computing particle SNRs"

	launch_childprocess("e2ctf.py --voltage {} --cs {} --ac {} --apix {} --allparticles --autofit {} -v {}".format(ctf.voltage,ctf.cs,ctf.ampcont,ctf.apix,dfopt,options.verbose-1))
	
	E2end(logid)

if __name__ == "__main__":
    main()
