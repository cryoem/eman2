#!/usr/bin/env python

#
# Author: Steven Ludtke, 10/18/2014 (sludtke@bcm.edu)
# Copyright (c) 2000- Baylor College of Medicine
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
import socket


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """e2framestats.py [options]

Extracts some information about each frame present in particles and generates a multicolumn text file for plotting. Current output is:

n SNR(20-200) SNR(4-10) quality defocus SNR1*SNR2

"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

#	parser.add_argument("--refine",type=str,default=None,help="Automatically get parameters for a refine directory")
	parser.add_argument("--output",type=str,help="Output text file file",default="frames.txt")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	out=file(options.output,"w")

	# Some info about the project for the header... may be useful when comparing multiple projects
	db=js_open_dict("info/project.json")
	try: out.write("# {}	# {} # {}\n".format(os.getcwd(),socket.gethostname(),os.getlogin()))
	except: out.write("# # # \n")
	try: out.write("# mass={} # apix={} # voltage={} # Cs={}\n".format(db["global.particle_mass"],db["global.apix"],db["global.microscope_voltage"],db["global.microscope_cs"]))
	except: out.write("# mass=-1 # apix=-1 # voltage=-1 # Cs=-1\n")
	out.write("#%N %SNR (20-200) %SNR(4-10) %Quality %Defocus %lohiprod %Defocus Diff %Defocus Ang %SNR(10-20)\n")

	files=["particles/"+i for i in os.listdir("particles") if "__" not in i and i[0]!="." and ".hed" not in i ]
	args.sort()

	for i,f in enumerate(files):
		inf=info_name(f)
		db=js_open_dict(inf)
		try: ctf=db["ctf"][0]
		except: continue
		r1=int(floor(1.0/(200.0*ctf.dsbg)))
		r2=int(ceil(1.0/(20.0*ctf.dsbg)))
		r3=int(floor(1.0/(10.0*ctf.dsbg)))
		r4=int(ceil(1.0/(4.0*ctf.dsbg)))
		losnr=sum(ctf.snr[r1:r2])/(r2-r1)
		medsnr=sum(ctf.snr[r2:r3])/(r3-r2)
		hisnr=sum(ctf.snr[r3:r4])/(r4-r3)
		qual=db["quality"]

		out.write("{:5d}\t{:6.4f}\t{:6.4f}\t{:1d}\t{:6.4f}\t{:7.5f}\t{:6.4f}\t{:6.1f}\t{:6.4f}\n".format(i,losnr,hisnr,qual,ctf.defocus,sqrt(losnr*hisnr),ctf.dfdiff,ctf.dfang,medsnr))
		db.close()

if __name__ == "__main__":  main()
