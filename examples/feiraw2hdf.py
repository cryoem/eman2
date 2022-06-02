#!/usr/bin/env python
# Read's FEI's 'RAW' format and writes to hdf
# feiraw2hdf.py <input> <input> ... <output> 
# also writes the sum of all inputs/1024 to output position 0

# This version modified so it only keeps the last 2 input files in the output (as well as the overall average)

from struct import unpack
from sys import argv
import sys,os
from EMAN2 import *

if argv[-1][-4:]!=".hdf" : 
	print("output should be an HDF file")
	sys.exit(1)

# iterate over input files
for i,fi in enumerate(argv[1:-1]):
	fin=open(fi,"rb")

	# check magic number and skip
	a=fin.read(13)
	if a[:7]!=b"FEI Raw" :
		print("Image magic number '%s', not 'FEI RawImage'"%a)
		sys.exit(1)

	# read the simple header
	(one,nx,ny,chan,bits,encoding,offset,stridex,stridey)=unpack("9I",fin.read(36))

	img=EMData(nx,ny)
	print("%s: %d x %d  %d "%(fi,nx,ny,encoding),end="")

	fin=None			# close the file for the process_region_io below

	if encoding==0 :		# signed int
		EMUtil.read_raw_emdata(img,fi,49,1,0,2,None)
		print("Signed Int")
	elif encoding==1 :		# unsigned int
		EMUtil.read_raw_emdata(img,fi,49,1,0,1,None)
		print("Unsigned Int")
	elif encoding==2 :		# float
		EMUtil.read_raw_emdata(img,fi,49,1,0,0,None)
		print("Float")
	else :
		print(f"Unknown encoding {encoding}")


	if i==0: imga=img.copy()
	else: imga.add(img)

	if i>4 : 
#		img.process_inplace("normalize.edgemean")	# not sure why someone thought this was a good idea
		img.write_image(argv[-1],i-4)	# write to output file

#imga.mult(1/1024.0)
#imga.process_inplace("normalize.edgemean")	# destroys the value of the image when used on things like gain references
imga.write_image(argv[-1],0)

#process_region_io(EMAN::EMData *ths, const char* path, size_t offset, int rw_mode, int image_index, size_t mode_size, const EMAN::Region * area=0
