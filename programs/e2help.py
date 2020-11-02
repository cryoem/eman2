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
		"sgirot":["SGI Spin-Axis (n1,n2,n3) vector with angle q","n1","FLOAT","X vector component","n2","FLOAT","Y vector component","n3","FLOAT","Z vector component","q","FLOAT","Angle of rotation in degrees"],
		"quaternion":["Standard 4 component quaternion (e0,e1,e2,e3)","e0","FLOAT","e0","e1","FLOAT","e1","e2","FLOAT","e2","e3","FLOAT","e3"]}
	
	elif args[0] in ("imagetypes", "imagetype", "image", "images", "imageformats", "imageformat"):
		print("Available image types:")
		
		header_line = ["Type", "Extension", "Read", "Write", "3D", "ImgStack", "VolStack", "RgnI/O", "Comments"]
		img_types = [
			 ["HDF5", "hdf", "Y", "Y", "Y", "Y", "Y", "Y", "HDF5 is an international standard for scientific data (http://www.hdfgroup.org/HDF5/). It supports arbitrary metadata (header info) and is very portable. This is the standard interchange format for EMAN2. Chimera can read EMAN2 style HDF files."],
			 ["DM2 (Gatan)", "dm2", "Y", "N", "N", "N", "N", "N", "Proprietary Gatan format (older version)"],
			 ["DM3 (Gatan)", "dm3", "Y", "N", "N", "N", "", "N", "Proprietary Gatan format from Digital Micrograph"],
			 ["DM4 (Gatan)", "dm4", "Y", "N", "Y", "Y", "", "N", "Proprietary Gatan format from Digital Micrograph, used with K2 cameras"],
			 ["SER (FEI)", "ser", "Y", "N", "N", "Y", "", "N", "Proprietary FEI format (Falcon camera ?)"],
			 ["EER (TF)", "eer", "Y", "N", "N", "Y", "N", "N ", "Falcon 4 camera counting mode format. Extremely large frame count with RLE compression to make frames very small. Supports up to 4x oversampling of counting data. Default reader is without oversampling. See below for details."],
			 ["EM", "em", "Y", "Y", "Y", "N", "", "Y", "As produced by the EM software package"],
			 ["ICOS", "icos", "Y", "Y", "Y", "N", "", "Y", "Old icosahedral format"],
			 ["Imagic", "img/hed", "Y", "Y", "Y", "Y", "Y", "Y", "This format stores header and image data in 2 separate files. Region I/O is only available for 2D. The Imagic format in EMAN2 is fully compatible with Imagic4D standard since the 2.0 release."],
			 ["MRC", "mrc", "Y", "Y", "N", "Y", "", "Y", "Largely compatible with CCP4. Note that some programs will treat 3D MRC files as stacks of 2D imagess (like IMOD). This behavior is partially supported in EMAN, but be aware that it is impossible to store metadata about each image in the stack when doing this, so it is not suitable as an export format for single particle work. EMAN2 supports reading of FEI MRC, which is an extended MRC format for tomography. The extra header information will be read into the header. All FEI MRC images will be 2-byte integer."],
			 ["MRCS", "mrcs", "Y", "Y", "N", "Y", "soon?", "Y", "Identical to MRC format above. If the filename is .mrcs, then a 3-D volume file will automatically be treated as a stack of 2-D images. If any other extension is used, it will appear to be a single 3-D volume."],
			 ["Spider Stack", "spi", "Y", "Y", "Y", "Y", "", "Y", "To read the overall image header in a stacked spider file, use image_index = -1."],
			 ["Spider Single", "spi", "Y", "Y", "Y", "N", "", "Y", "Specify \"--outtype=spidersingle\" to use with e2proc2d/3d"],
			 ["SER", "ser", "Y", "N", "N", "Y", "", "N", "Also known as TIA (Emospec) file format, used by FEI Tecnai and Titan microscope for acquiring and displaying scanned images and spectra"],
			 ["BDB", "N/A", "Y", "Y", "Y", "Y", "", "Y", "This entry is for EMAN2's (retired) embedded database system. While it is still possible to read/write BDB's for backwards compatibility, we do not suggest any new use of this format in EMAN2 (SPARX still uses it for many operations)"],
			 ["Amira", "am", "Y", "Y", "Y", "N", "", "N", "A native format for the Amira visualization package"],
			 ["DF3", "df3", "Y", "Y", "Y", "N", "", "N", "File format for POV-Ray, support 8,16,32 bit integer per pixel"],
			 ["FITS", "fts", "Y", "N", "Y", "N", "", "N", "Widely used file format in astronomy"],
			 ["JPEG", "jpg/jpeg", "N", "Y", "N", "N", "", "N", "Note that JPEG images use lossy compression and are NOT suitable for quantitative analysis. PNG (lossless compression) is a better alternative unless file size is of critical importance."],
			 ["LST", "lst", "Y", "Y", "Y", "Y", "", "N", "ASCII file contains a list of image file names and numbers. Two variants, LST and LSX. LSX is normally used in EMAN2 and has the additional restraint that all lines have the same length."],
			 ["LSTFAST", "lsx/lst", "Y", "Y", "Y", "Y", "", "N", "Optomized version of LST"],
			 ["OMAP", "omap", "Y", "N", "Y", "N", "", "N", "Also called DSN6 map, 1 byte integer per pixel"],
			 ["PGM", "pgm", "Y", "Y", "N", "N", "", "N", "Standard graphics format with 8 bit greyscale images. No compression."],
			 ["PIF", "pif", "Y", "Y", "Y", "Y", "", "N", "Purdue Image Format. This will read most, but not all PIF images. Recent support added for mode 40 and 46 (boxed particles). Some of the FFT formats cannot be read by EMAN2. PIF writing is normally done in FLOAT mode, which is not used very often in PIF. PIF technically permits only images with odd dimensions, EMAN does not enforce this."],
			 ["PNG", "png", "Y", "Y", "N", "N", "", "N", "Excellent format for presentations. Lossless data compression, 8 bit or 16 bit per pixel"],
			 ["SAL", "hdr/img", "Y", "N", "N", "N", "", "N", "Scans-A-Lot. Old proprietary scanner format. Separate header and data file"],
			 ["SITUS", "situs", "Y", "Y", "Y", "N", "", "N", "Situs-specific ASCII format on a cubic lattice. Used by Situs programs"],
			 ["TIFF", "tiff/tif", "Y", "Y", "N", "N", "", "N", "Good format for use with programs like photoshop. Some variants are good for quantitative analysis, but JPEG compression should be avoided."],
			 ["V4L", "v4l", "Y", "N", "N", "N", "", "N", "Used by some video-capture boards in Linux. Acquires images from the V4L2 interface in real-time(video4linux)."],
			 ["VTK", "vtk", "Y", "Y", "Y", "N", "", "N", "Native format from Visualization Toolkit"],
			 ["XPLOR", "xplor", "Y", "Y", "Y", "N", "", "N", "8 bytes integer, 12.5E float ASCII format"]]
		
		type_len_max = max([len(img[0]) for img in img_types])
		
		print(f"{header_line[0]:>{type_len_max}}: " + " ".join(header_line[1:-1]))
		
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
