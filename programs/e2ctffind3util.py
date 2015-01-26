#!/usr/bin/env python

#
# Author: Stephen Murray, 11/5/2014 (scmurray@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
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
# e2ctffind3.py  11/5/2014 Stephen Murray
# This is a program for interacting with ctffind3 and ctffind3 results from within EMAN2

from EMAN2 import *
from EMAN2db import db_open_dict, db_close_dict, db_check_dict, db_list_dicts
from optparse import OptionParser
from math import *
import os
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] 
For more information on ctffind3 please see: Mindell, JA, Grigorieff N.  2003.  Accurate determination of local defocus and specimen tilt in electron microscopy. J Struct Biol. 142:334-47.

"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_header(name="ctffind3header", help='Options below this label are specific to e2ctffind3util.py', title="### e2ctffind3util options (requires that ctffind3 be installed)###", row=0, col=0, rowspan=1, colspan=2, mode="import,run")

	#options associated with e2ctffind3.py
	parser.add_argument("--apix", default=0.0, type=float,help="The angstrom per pixel of the micrographs", guitype='floatbox', row=3, col=0, rowspan=1, colspan=1, mode="import,run")
	parser.add_argument("--cs", default=0.0, type=float,help="The spherical aberration of the microscope", guitype='floatbox', row=3, col=1, rowspan=1, colspan=1, mode="import,run")
	parser.add_argument("--voltage", default=0.0, type=float,help="The voltage (in kV) of the microscope", guitype='floatbox', row=4, col=0, rowspan=1, colspan=1, mode="import,run")
	parser.add_argument("--ac", default=0.0, type=float,help="The amplitude contrast of the micrographs", guitype='floatbox', row=4, col=1, rowspan=1, colspan=1, mode="import,run")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness", guitype='intbox', row=5, col=0, rowspan=1, colspan=1, mode="import,run")
	parser.add_argument("--importctffind3", default=False, action="store_true",help="Import ctffind3 data?", guitype='boolbox', row=6, col=0, rowspan=1, colspan=1, mode='import[True]')
	parser.add_argument("--importctffind4", default=False, action="store_true",help="Import ctffind4 data?", guitype='boolbox', row=6, col=1, rowspan=1, colspan=1, mode='import[False]')
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	parser.add_pos_argument(name="micrographs",help="List the micrographs to run ctffind3 on here.", default="", guitype='filebox', browser="EMRawDataTable(withmodal=True,multiselect=True)", filecheck=False, row=1, col=0,rowspan=1, colspan=2, mode='run')
	parser.add_argument("--allmicrographs", default=False, action="store_true",help="Run Ctffind3 on all micrographs in the micrographs directory?", guitype='boolbox', row=2, col=0, rowspan=1, colspan=1, mode="run")
	parser.add_argument("--maxres", default=0.0, type=float,help="The highest resolution to be fitted (Angstroms)", guitype='floatbox', row=6, col=0, rowspan=1, colspan=1, mode="run")
	parser.add_argument("--minres", default=0.0, type=float,help="The lowest resolution to be fitted (Angstroms)", guitype='floatbox', row=6, col=1, rowspan=1, colspan=1, mode="run")
	parser.add_argument("--defocusmin", default=0.0, type=float,help="The starting defocus value for grid search (microns)", guitype='floatbox', row=7, col=0, rowspan=1, colspan=1, mode="run")
	parser.add_argument("--defocusmax", default=0.0, type=float,help="The end defocus value for grid search (microns)", guitype='floatbox', row=7, col=1, rowspan=1, colspan=1, mode="run")
	parser.add_argument("--defocusstep", default=0.0, type=float,help="The step width for grid search (microns)", guitype='floatbox', row=8, col=0, rowspan=1, colspan=1, mode="run")
	parser.add_argument("--runctffind3", default=False, action="store_true",help="Run ctffind3 on the selected micrographs?", guitype='boolbox', row=9, col=0, rowspan=1, colspan=1, mode='run[True]')
	parser.add_argument("--runctffind4", default=False, action="store_true",help="Run ctffind4 on the selected micrographs?", guitype='boolbox', row=9, col=1, rowspan=1, colspan=1, mode='run[False]')
	parser.add_argument("--windowsize", default=0, type=int,help="The amplitude contrast of the micrographs", guitype='intbox', row=8, col=1, rowspan=1, colspan=1, mode="run")
	
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv,options.ppid)
	if options.apix<= 0 :
		print "Angstrom per pixel (apix) must be specified!"
		exit(-1)
	if options.cs<= 0 :
		print "Spherical Aberration (cs) must be specified!"
		exit(-2)
	if options.voltage<= 0 :
		print "Voltage must be specified!"
		exit(-3)
	if options.ac<= 0 :
		print "Amplitude Contrast must be specified!"
		exit(-4)
		
	if options.allmicrographs:
		args = []
		for item in os.listdir("micrographs"):
			try: 
				if item.split(".")[1] == "mrc" or item.split(".")[1] == "hdf":
					args.append("micrographs/" + item)
			except:
				pass
	
	if options.runctffind4:
		version = "ctffind4"
		run_ctffind(options.apix, args, options.cs, options.voltage, options.ac, options.windowsize, options.minres, options.maxres, options.defocusmin*10000, options.defocusmax*10000, options.defocusstep*10000, options.verbose,version)
	elif options.runctffind3:
		version = "ctffind3"
		run_ctffind(options.apix, args, options.cs, options.voltage, options.ac, options.windowsize, options.minres, options.maxres, options.defocusmin*10000, options.defocusmax*10000, options.defocusstep*10000, options.verbose, version)

	if options.importctffind4 or options.runctffind4:
		version = "ctffind4"
	else:
		version = "ctffind3"
	import_ctf(options.voltage, options.cs, options.ac, options.apix, options.verbose, version)

	print "e2ctffind3util.py complete!"
	E2end(logid)

def import_ctf(voltage, cs, ac, apix, verbose, version):
	if not os.path.exists(version):
		print "no " + version + " directory found. Please see usage instructions!"
		exit(-5)
		
	for filename in os.listdir("micrographs"):
		if os.path.exists(version + "/" + base_name(filename)+"_" + version + ".log"):
			f_log = open(version + "/" + base_name(filename) + "_" + version + ".log")
			for line in f_log.readlines():
				if len(line) > 1:
					if len(line.split()) == 6:
						if line.split()[5] == "Values":
							defocusu = float(line.split()[0])
							defocusv = float(line.split()[1])
							dfang =  float(line.split()[2])
							cc = float(line.split()[3])
							e2defocus = (defocusu + defocusv) / 20000.0
							e2dfdiff = abs(defocusu - defocusv) / 10000.0
							e2dfang = dfang
							if not os.path.exists(os.getcwd() + "/info"):
								os.mkdir(os.getcwd() + "/info")
							jdb = js_open_dict(os.getcwd() + "/info/" + base_name(filename) + "_info.json")
							if "ctf" in jdb.keys():
								jdb.delete('ctf')
							ctf = EMAN2Ctf()
							ctf.from_dict({"defocus":e2defocus,"dfang":e2dfang,"dfdiff":e2dfdiff,"voltage":voltage,"cs":cs,"ampcont":ac,"apix":apix})
							jdb['ctf_frame'] = [512,ctf,(256,256),tuple(),5,1]
							#launch_childprocess("e2ctf.py --voltage {} --cs {} --ac {} --apix {} --autofit --curdefocusfix --verbose {} {}".format(voltage,cs,ac,apix,verbose-1,))
	launch_childprocess("e2ctf.py --voltage {} --cs {} --ac {} --apix {} --allparticles --autofit --curdefocusfix --astigmatism --verbose {}".format(voltage,cs,ac,apix,verbose-1))

def run_ctffind(apix, args, cs, voltage, ac, windowsize, minres, maxres, defocusmin, defocusmax, defocusstep,verbose, version):
	print "Running " + version
	dstep = 10.0
	mag = dstep / apix * 10000
	
	try: os.mkdir(version)
	except: pass
	created = False
	for image in args:
		if not os.path.exists(image):
			print "Image Does not exist: " + image
			exit(-6)
		card = open("card.txt",'w')
		if image.split(".")[1] != "mrc":
			launch_childprocess("e2proc2d.py {} {}.mrc --verbose={}".format(image,image.split(".")[0],verbose-1))
			created = True
		card.write(image.split(".")[0] + ".mrc\n" + version + "/" + base_name(image) + "_" + version + ".ctf\n" + str(cs) + "," + str(voltage) + "," + str(ac) + "," + str(mag) + "," + str(dstep) + "\n" + str(windowsize) + "," + str(minres) + "," + str(maxres) + "," + str(defocusmin) + "," + str(defocusmax) + "," + str(defocusstep))
		card.close()
		print "running " + version + " on: " + image
		if version == "ctffind3":
			s = "`which ctffind3.exe` < card.txt >ctffind3/" + base_name(image) + "_ctffind3.log"
		else:
			s = "`which ctffind` --old-school-input < card.txt >ctffind4/" + base_name(image) + "_ctffind4.log"
		call(s,shell=True)
		if created:
			os.remove(image.split(".")[0] + ".mrc")
			created = False
	#os.remove("card.txt")
	
if __name__=="__main__":
	main()
