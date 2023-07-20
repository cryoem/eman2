#!/usr/bin/env python
#
# Author: Steven Ludtke  07/19/2023
# Copyright (c) 2023- Baylor College of Medicine
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
#

from EMAN3 import *
import os
import xml.etree.ElementTree as ET
import sys

def main():

	usage="""e3iflm.py <iflm_xml_file>

Turn an iflm focal series into a 3-D volume
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
#	parser.add_argument("--est_gain", type=str,help="specify output file for gain image. Estimates a gain image when given a set of many movies via hierarchical median estimation", default=None)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")

	(options, args) = parser.parse_args()

	tree=ET.parse(args[0])
	el_root=tree.getroot()

	el_images=el_root.find("ImageMatrix").find("Images").findall("Image")

	umx=float(el_root.find("ImageMatrix").find("TileWidth").text)		# X image size in microns
	umy=float(el_root.find("ImageMatrix").find("TileHeight").text) 	# Y image size in microns
	nx=int(el_root.find("ImageMatrix").find("TilePixelWidth").text)
	ny=int(el_root.find("ImageMatrix").find("TilePixelHeight").text)
	nz=len(el_images)

	focus1=float(el_images[1].find("Position").find("{http://www.thermofisher.com/schemas/sem}Position").find("Focus").text)
	focus0=float(el_images[0].find("Position").find("{http://www.thermofisher.com/schemas/sem}Position").find("Focus").text)

	map=EMData(nx,ny,nz)
	map["apix_x"]=umx*10000.0/nx
	map["apix_y"]=umy*10000.0/ny
	map["apix_z"]=(focus1-focus0)*10000

	for i,el_image in enumerate(el_images):
		tiff_path=el_image.find("RelativePath").text.replace("\\","/")
		focus=float(el_image.find("Position").find("{http://www.thermofisher.com/schemas/sem}Position").find("Focus").text)
		plane=int(el_image.find("Index").find("Plane").text)

		img=EMData(tiff_path)
	#	img.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.4})
		map.insert_clip(img,(0,0,i))
		
		print(".",end="")
		sys.stdout.flush()
	
	print("")

#	map.process_inplace("mask.onlypeaks",{"xsize":0,"ysize":0,"zsize":4})
	map.write_image(args[0].replace(".xml",".hdf:12"),0)

if __name__ == '__main__':
	main()
