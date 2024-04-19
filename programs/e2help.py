#!/usr/bin/env python
#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
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
#

# e2help.py  07/23/3004  Steven Ludtke
# This program will provide a variety of EMAN2 help

from builtins import range
from EMAN2 import *
from EMAN2_meta import EMANVERSION, FULLVERSIONSTRING
from math import *
import os
import sys
from sys import exit


def main():
	progname = os.path.basename(sys.argv[0])
	helpstring =  """Help is available on the following topics:
boxsizes, processors, cmps, aligners, averagers, projectors, reconstructors, analyzers, symmetries, orientgens, rotationtypes, imagetypes"""
	usage = """prog <topic> [contains]
	
Interactive help on a variety of the eman2 library's modular functions. The optional 'contains' argument will
act as a filter on the names of the algorithms."""
	usage += " "+helpstring

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	#parser.add_argument("--res", "-R", type=float, help="Resolution in A, equivalent to Gaussian lowpass with 1/e width at 1/res",default=2.8)
	#parser.add_argument("--box", "-B", type=str, help="Box size in pixels, <xyz> or <x>,<y>,<z>")
	parser.add_argument("--gui", action="store_true", help="Use the GUI for display help", default=False)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
		
	if options.gui:
		from e2projectmanager import TheHelp
		from eman2_gui.emapplication import EMApp
		app = EMApp()
		thehelp = TheHelp()
		thehelp.show()
		if args:
			print(args[0])
			if args[0] in ("aligner","aligners"):
				thehelp._helpchange(0)
			elif args[0] in ("analyzer","analyzers"):
				thehelp._helpchange(1)
			elif args[0] in ("averager","averagers"):
				thehelp._helpchange(2)
			elif args[0] in ("cmp","cmps"):
				thehelp._helpchange(3)
			elif args[0] in ("orientgen","orientationgen","orientgens","orientationgens","orientationgenerators"):
				thehelp._helpchange(4)
			elif args[0] in ("processor","processors"):
				thehelp._helpchange(5)
			elif args[0] in ("projector","projectors"):
				thehelp._helpchange(6)
			elif args[0] in ("reconstructor","reconstructors"):
				thehelp._helpchange(7)
			elif args[0] in ("sym","symmetry","symmetries"):
				thehelp._helpchange(8)
		app.exec_()
		exit(0)

	if len(args)<1 : 
		print(helpstring)
		exit(0)
		
	l=None
	if args[0] in ("cmp","cmps") :
		print("Available comparators:")
		l=dump_cmps_list()
	elif args[0] in ("analyzer","analyzers") :
		print("Available analysers:")
		l=dump_analyzers_list()
	elif args[0] in ("averager","averagers") :
		print("Available averagers:")
		l=dump_averagers_list()
	elif args[0] in ("processor","processors") :
		print("Available processors:")
		l=dump_processors_list()
	elif args[0] in ("projector","projectors") :
		print("Available projectors:")
		l=dump_projectors_list()
	elif args[0] in ("reconstructor","reconstructors") :
		print("Available reconstructors:")
		l=dump_reconstructors_list()
	elif args[0] in ("aligner","aligners") :
		print("Available aligners:")
		l=dump_aligners_list()
	elif args[0] in ("sym","symmetry","symmetries") :
		print("Available symmetries:")
		l=dump_symmetries_list()
	elif args[0] in ("orientgen","orientationgen","orientgens","orientationgens","orientationgenerators") :
		print("Available orientation generators:")
		l=dump_orientgens_list()
	elif args[0][:3]=="box":
		if len(args)>1 :
			print(good_size(int(args[1])))
		else:
			print(good_box_sizes)
		exit(0)
	elif args[0][:8]=="rotation" :
		print("Available rotation conventions:")
		l={"eman":["EMAN convention, az(Z),alt(X),phi(Z') Eulers","alt","FLOAT","Altitude, X-axis","az","FLOAT","Azimuth, Z-axis","phi","FLOAT","Z' Axis. in-plane rotation in 2-D"],
		"imagic":["IMAGIC convention","alpha","FLOAT","alpha","beta","FLOAT","beta","gamma","FLOAT","gamma"],
		"spider":["SPIDER convention","phi","FLOAT","phi","theta","FLOAT","theta","psi","FLOAT","psi"],
		"mrc":["MRC/CCP4 convention","omega","FLOAT","omega","theta","FLOAT","theta","psi","FLOAT","psi"],
		"xyz":["XYZ convention (Chimera)","x","FLOAT","X-axis","y","FLOAT","Y-axis","z","FLOAT","Z-axis"],
		"spin":["Spin-Axis (n1,n2,n3) vector with angle omega","n1","FLOAT","X vector component","n2","FLOAT","Y vector component","n3","FLOAT","Z vector component","omega","FLOAT","Angle of rotation in degrees"],
		"spinvec":["Spin vector, length determines rotation amount, 0.5 -> 180 deg","v1","FLOAT","X vector component","v2","FLOAT","Y vector component","v3","FLOAT","Z vector component"],
		"sgirot":["SGI Spin-Axis (n1,n2,n3) vector with angle q","n1","FLOAT","X vector component","n2","FLOAT","Y vector component","n3","FLOAT","Z vector component","q","FLOAT","Angle of rotation in degrees"],
		"quaternion":["Standard 4 component quaternion (e0,e1,e2,e3)","e0","FLOAT","e0","e1","FLOAT","e1","e2","FLOAT","e2","e3","FLOAT","e3"]}
	
	elif args[0] in ("imagetypes", "imagetype", "image", "images", "imageformats", "imageformat"):
		print("Available image types:")
		
		header_line = ["Type", "Extension", "Read", "Write", "3D", "ImgStack", "VolStack", "RgnI/O"]
		img_types = [
			 ["HDF5",          "hdf",      "Y", "Y", "Y", "Y", "Y",     "Y"],
			 ["DM2 (Gatan)",   "dm2",      "Y", "N", "N", "N", "N",     "N"],
			 ["DM3 (Gatan)",   "dm3",      "Y", "N", "N", "N", "",      "N"],
			 ["DM4 (Gatan)",   "dm4",      "Y", "N", "Y", "Y", "",      "N"],
			 ["SER (FEI)",     "ser",      "Y", "N", "N", "Y", "",      "N"],
			 ["EER (TF)",      "eer",      "Y", "N", "N", "Y", "N",     "N "],
			 ["EM",            "em",       "Y", "Y", "Y", "N", "",      "Y"],
			 ["ICOS",          "icos",     "Y", "Y", "Y", "N", "",      "Y"],
			 ["Imagic",        "img/hed",  "Y", "Y", "Y", "Y", "Y",     "Y"],
			 ["MRC",           "mrc",      "Y", "Y", "N", "Y", "",      "Y"],
			 ["MRCS",          "mrcs",     "Y", "Y", "N", "Y", "soon?", "Y"],
			 ["Spider Stack",  "spi",      "Y", "Y", "Y", "Y", "",      "Y"],
			 ["Spider Single", "spi",      "Y", "Y", "Y", "N", "",      "Y"],
			 ["SER",           "ser",      "Y", "N", "N", "Y", "",      "N"],
			 ["BDB",           "N/A",      "Y", "Y", "Y", "Y", "",      "Y"],
			 ["Amira",         "am",       "Y", "Y", "Y", "N", "",      "N"],
			 ["DF3",           "df3",      "Y", "Y", "Y", "N", "",      "N"],
			 ["FITS",          "fts",      "Y", "N", "Y", "N", "",      "N"],
			 ["JPEG",          "jpg/jpeg", "N", "Y", "N", "N", "",      "N"],
			 ["LST",           "lst",      "Y", "Y", "Y", "Y", "",      "N"],
			 ["LSTFAST",       "lsx/lst",  "Y", "Y", "Y", "Y", "",      "N"],
			 ["OMAP",          "omap",     "Y", "N", "Y", "N", "",      "N"],
			 ["PGM",           "pgm",      "Y", "Y", "N", "N", "",      "N"],
			 ["PIF",           "pif",      "Y", "Y", "Y", "Y", "",      "N"],
			 ["PNG",           "png",      "Y", "Y", "N", "N", "",      "N"],
			 ["SAL",           "hdr/img",  "Y", "N", "N", "N", "",      "N"],
			 ["SITUS",         "situs",    "Y", "Y", "Y", "N", "",      "N"],
			 ["TIFF",          "tiff/tif", "Y", "Y", "N", "N", "",      "N"],
			 ["V4L",           "v4l",      "Y", "N", "N", "N", "",      "N"],
			 ["VTK",           "vtk",      "Y", "Y", "Y", "N", "",      "N"],
			 ["XPLOR",         "xplor",    "Y", "Y", "Y", "N", "",      "N"]]
		
		type_len_max = max([len(img[0]) for img in img_types])
		
		print(f"{header_line[0]:>{type_len_max}}: " + " ".join(header_line[1:]))
		
		for img in img_types:
			print(  f"{img[0]:<{type_len_max}}: " 
				  + f"{img[1]:{len(header_line[1])}}"
				  + f" {img[2]:^{len(header_line[2])}}"
				  + f" {img[3]:^{len(header_line[3])}}"
				  + f" {img[4]:^{len(header_line[4])}}"
				  + f" {img[5]:^{len(header_line[5])}}"
				  + f" {img[6]:^{len(header_line[6])}}"
				  + f" {img[7]:^{len(header_line[7])}}"
				  )


	elif args[0] in ("version"):
		print(FULLVERSIONSTRING) 

	elif args[0] in ("tophat"):
		print("There are multiple ways to filter the 3D maps in SPA or SPT refinements:")
		print("    wiener: wiener filter based on FSC curve. default mode in most programs.")
		print("    global: tophat filter across the map at the resolution cutoff 0.143 from fsc_masked_xx.txt.") 
		print("    localwiener: wiener filter based on the fsc curve of local regions from the even/odd maps. see e2fsc_real_local.py")
		print("    local: tophat filter based on local resolution calculated from the even/odd maps at 0.143 cutoff. see e2fsc_real_local.py") 
		print("    localold: an old version of local resolution based filter. see e2fsc_local.py") 
	else:
		print(helpstring)
		print("unknown option:",args[0])
		
	if l:
		if options.verbose>0:
			if len(args)>1 : k=[i for i in list(l.keys()) if args[1] in i]
			else: k=list(l.keys())
			k.sort()
			for i in k:
				print("%s : %s"%(i, l[i][0]))
				for j in range(1,len(l[i]),3): 
					print("\t%s(%s) - %s"%(l[i][j],l[i][j+1],l[i][j+2]))
		else :
			if len(args)>1 : k=[i for i in list(l.keys()) if args[1] in i]
			else: k=list(l.keys())
			if len(k)==0 :
				print("Empty list - no items met search criteria")
				sys.exit(0)
			maxk=max([len(ii) for ii in k])
			fmt="%%-%0ds : "%maxk
			k.sort()
			for i in k:
				print(fmt%i, end=' ')
				for j in range(1,len(l[i]),3): 
					print("%s(%s)  "%(l[i][j],l[i][j+1]), end=' ')
				if len(k)>1: print("")

if __name__ == "__main__":
	main()
