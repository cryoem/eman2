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

from past.utils import old_div
from builtins import range
from EMAN3 import *
from EMAN3jax import *
import sys
import os.path
import math
import random
import os
import datetime
import time
import traceback

# constants

HEADER_ONLY=True
HEADER_AND_DATA=False

xyplanes = ['xy', 'yx']
xzplanes = ['xz', 'zx']
yzplanes = ['yz', 'yz']
threedplanes = xyplanes + xzplanes + yzplanes


def changed_file_name(input_name, output_pattern, input_number, multiple_inputs):
	# convert an input file name to an output file name
	# by replacing every @ or * in output_pattern with
	# the input file base name (no extension)

	outname = output_pattern

	if multiple_inputs or 'FILL_ONE_FILE' in os.environ :
		base_name = os.path.basename(os.path.splitext(input_name)[0])

		if '@' in outname :
			outname = outname.replace('@', base_name)

		if '*' in outname :
			outname = outname.replace('*', base_name)

		if '%i' in outname :
			outname = outname.replace('%i', str(input_number))

		if '%j' in outname :
			outname = outname.replace('%j', str(input_number-1))

		if '%d' in outname :
			dt = datetime.datetime.now()
			date_str = str(dt.year) + '_' + str(dt.month) + '_' + str(dt.day)
			outname = outname.replace('%d', date_str)

		if '%t' in outname :
			dt = datetime.datetime.now()
			time_str = str(dt.hour) + '_' + str(dt.minute) + '_' + str(dt.second)
			outname = outname.replace('%t', time_str)

		if '%%' in outname :
			outname = outname.replace('%%', '%')

	return outname

def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + """ [options] <inputfile> ... <inputfile>

	NOTE to EMAN2 users: Output must be specified explicitly with --output option. It is no longer the last argument.

	Generic 2-D image processing and file format conversion program. Acts on stacks of 2-D images
	(multiple images in one file). All EMAN2 recognized file formats accepted (see Wiki for list).

	MRC stack files MUST use the .mrcs extension to be treated as a set of 2-D images

	--output can use {} values to make individual output paths for each input:
	n - file number in command-line order starting at 0
	base - the 'base name' of the input, excludes path, extension and anything following "__"


	"""

	parser = EMArgumentParser(usage=usage,allow_abbrev=False,version=EMANVERSION)

	parser.add_argument("--output", type=str, help="Output filename or pattern. Appends output images to existing image files. ",default=None)
	parser.add_argument("--apix", type=float, help="Override A/pixel from header on all inputs",default=0)
	parser.add_argument("--fouriershrink", metavar="n", type=float, action="append", help="Reduce an image size by an arbitrary scaling factor by clipping in Fourier space. eg - 2 will reduce image size to 1/2.")
	parser.add_argument("--outmode", type=str, default="float", help="All EMAN3 programs write images with 4-byte floating point values by default. This allows specifying an alternate format when supported (float, int8, int16, int32, uint8, uint16, uint32). Values are rescaled to fill MIN-MAX range. The ':bits' filename specification is the preferred approach for this.")
	parser.add_argument("--createnew", type=str, default=None, help="<nx>,<ny>[,<nimg=1>[,<value=0>]]  If specified with no input filenames, will create nimg nx x ny images with constant value to use as input.")
	parser.add_argument("--normalize", type=str, default=None, help="<mode>  Normalize images: standard, edgemean")

	# parser.add_argument("--average", action="store_true", help="Averages all input images (without alignment) and writes a single output image")
	# parser.add_argument("--avgseq", type=int,default=0, help="Averages sets of N sequential frames. eg - if N=4 and the input contains 100 images, the output would be 25 images")
	# parser.add_argument("--averager",type=str,help="If --average is specified, this is the averager to use (e2help.py averager). Default=mean",default="mean")
	# parser.add_argument("--calcsf", metavar="outputfile", type=str, help="calculate a radial structure factor for the image and write it to the output file, must specify apix. divide into <n> angular bins")
	# parser.add_argument("--calccont", action="store_true", help="Compute the low resolution azimuthal contrast of each image and put it in the header as eval_contrast_lowres. Larger values imply more 'interesting' images.")
	# parser.add_argument("--clip", metavar="xsize,ysize", type=str, action="append", help="Specify the output size in pixels xsize,ysize[,xcenter,ycenter], images can be made larger or smaller.")
	# parser.add_argument("--exclude", metavar="exclude-list-file", type=str, help="Excludes image numbers, either a list of comma separated values, or a filename with one number per line, first image == 0")
	# parser.add_argument("--fftavg", metavar="filename", type=str, help="Incoherent Fourier average of all images and write a single power spectrum image")
	# parser.add_argument("--process", metavar="processor_name:param1=value1:param2=value2", type=str, action="append", help="apply a processor named 'processorname' with all its parameters/values.")
	# parser.add_argument("--mult", metavar="k", type=float, help="Multiply image by a constant. mult=-1 to invert contrast.")
	# parser.add_argument("--add", metavar="f", type=float,action="append",help="Adds a constant 'f' to the densities")
	# parser.add_argument("--addfile", type=str, action="append",help="Adds the volume to another volume of identical size")
	# parser.add_argument("--first", metavar="n", type=int, default=0, help="the first image in the input to process [0 - n-1])")
	# parser.add_argument("--last", metavar="n", type=int, default=-1, help="the last image in the input to process")
	# parser.add_argument("--list", metavar="listfile", type=str, help="Works only on the image numbers in LIST file")
	# parser.add_argument("--randomn", metavar="n", type=int, default=0, help="Selects a random subset of N particles from the file to operate on.")
	# parser.add_argument("--inplace", action="store_true", help="Output overwrites input, USE SAME FILENAME, DO NOT 'clip' images.")
	# parser.add_argument("--interlv", metavar="interleave-file", type=str, help="Specifies a 2nd input file. Output will be 2 files interleaved.")
	# parser.add_argument("--extractboxes",default=False, action="store_true",help="Extracts box locations from the image header to produce a set of .box files for only the particles in the .lst files")
	# parser.add_argument("--meanshrink", metavar="n", type=float, action="append", help="Reduce an image size by an integral (1.5 also allowed) scaling factor using average. eg - 2 will reduce image size to 1/2. Clip is not required.")
	# parser.add_argument("--medianshrink", metavar="n", type=int, action="append", help="Reduce an image size by an integral scaling factor, uses median filter. eg - 2 will reduce image size to 1/2. Clip is not required.")
	# parser.add_argument("--fouriershrink", metavar="n", type=float, action="append", help="Reduce an image size by an arbitrary scaling factor by clipping in Fourier space. eg - 2 will reduce image size to 1/2.")
	# parser.add_argument("--mraprep",  action="store_true", help="this is an experimental option")
	# parser.add_argument("--compressbits", type=int,help="HDF only. Bits to keep for compression. -1 for no compression",default=-1)
	# parser.add_argument("--outnorescale", action="store_true", default=False, help="If specified, floating point values will not be rescaled when writing data as integers. Values outside of range are truncated.")
	# parser.add_argument("--mrc16bit",  action="store_true", help="(deprecated, use --outmode instead) output as 16 bit MRC file")
	# parser.add_argument("--mrc8bit",  action="store_true", help="(deprecated, use --outmode instead) output as 8 bit MRC file")
	# parser.add_argument("--fixintscaling", type=str, default=None, help="When writing to an 8 or 16 bit integer format the data must be scaled. 'noscale' will assume the pixel values are already correct, 'full' will insure the full range of values are included in the output, 'sane' will pick a good range, a number will set the range to mean+=sigma*number")
	# parser.add_argument("--multfile", type=str, action="append", help="Multiplies the volume by another volume of identical size. This can be used to apply masks, etc.")
 #
	# parser.add_argument("--norefs", action="store_true", help="Skip any input images which are marked as references (usually used with classes.*)")
	# parser.add_argument("--outtype", metavar="image-type", type=str, help="output image format, 'mrc', 'imagic', 'hdf', etc. if specify spidersingle will output single 2D image rather than 2D stack.")
	# parser.add_argument("--radon",  action="store_true", help="Do Radon transform")
	# parser.add_argument("--randomize", type=str, action="append",help="Randomly rotate/translate the image. Specify: da,dxy,flip  da is a uniform distribution over +-da degrees, dxy is a uniform distribution on x/y, if flip is 1, random handedness changes will occur")
	# parser.add_argument("--rotavg", action="store_true", help="Compute the 1-D rotational average of each image as a final step before writing the output", default=False)
	# parser.add_argument("--rotate", type=float, action="append", help="Rotate clockwise (in degrees)")
	# parser.add_argument("--rfp",  action="store_true", help="this is an experimental option")
	# parser.add_argument("--fp",  type=int, help="This generates rotational/translational 'footprints' for each input particle, the number indicates which algorithm to use (0-6)")
	# parser.add_argument("--scale", metavar="f", type=float, action="append", help="Scale by specified scaling factor. Clip must also be specified to change the dimensions of the output map.")
	# parser.add_argument("--anisotropic", type=str,action="append", help="Anisotropic scaling, stretches on one axis and compresses the orthogonal axis. Specify amount,angle. See e2evalrefine")
	# parser.add_argument("--selfcl", metavar="steps mode", type=int, nargs=2, help="Output file will be a 180x180 self-common lines map for each image.")
	# parser.add_argument("--setsfpairs",  action="store_true", help="Applies the radial structure factor of the 1st image to the 2nd, the 3rd to the 4th, etc")
	# parser.add_argument("--split", metavar="n", type=int, help="Splits the input file into a set of n output files")
	# parser.add_argument("--translate", type=str, action="append", help="Translate by x,y pixels")
	# parser.add_argument("--headertransform", type=int, action="append", help="This will take the xform.align2d header value from each particle, and apply it. Pass 0 to perform the transform or 1 to perform the inverse.")
	# parser.add_argument("--plane", metavar=threedplanes, type=str, default='xy', help="Change the plane of image processing, useful for processing 3D mrcs as 2D images.")
	# parser.add_argument("--writejunk", action="store_true", help="Writes the image even if its sigma is 0.", default=False)
	# parser.add_argument("--swap", action="store_true", help="Swap the byte order", default=False)
	# parser.add_argument("--threed2threed", action="store_true", help="Process 3D image as a stack of 2D slices, then output as a 3D image", default=False)
	# parser.add_argument("--threed2twod", action="store_true", help="Process 3D image as a stack of 2D slices, then output as a 2D stack", default=False)
	# parser.add_argument("--twod2threed", action="store_true", help="Process a stack of 2D images, then output as a 3D image.", default=False)
	# parser.add_argument("--unstacking", action="store_true", help="Process a stack of 2D images, then output as a series of numbered single image files", default=False)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, help="verbose level [0-9], higher number means higher level of verboseness",default=1)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)

	# eer_input_group = parser.add_mutually_exclusive_group()
	# eer_input_group.add_argument("--eer2x", action="store_true", help="Render EER file on 8k grid.")
	# eer_input_group.add_argument("--eer4x", action="store_true", help="Render EER file on 16k grid.")


	append_options = ["anisotropic","clip", "process", "meanshrink", "medianshrink", "fouriershrink", "scale", "randomize", "rotate", "translate", "multfile","addfile","add", "headertransform"]

#	optionlist = get_optionlist(sys.argv[1:])

	(options, args) = parser.parse_args()

	if options.outmode not in file_mode_map:
		error_exit("Invalid output mode, please specify one of :\n",str(list(file_mode_map.keys())).translate(None,'"[]'))

	# user requested creation of new empty image stack to use instead of input files
	if options.createnew is not None:
		if len(args)>0: error_exit("--createnew cannot be used with input filenames")
		if os.path.exists("tmp_input_p2d.hdf"): error_exit("--createnew makes a temporary file 'tmp_input_p2d.hdf', but this file already exists. Please remove and try again.")
		parm=options.createnew.split(",")
		if len(parm)<2: error_exit("--createnew <nx>,<ny>,<nimg>,<value>")
		img=EMData(int(parm[0]),int(parm[1]),1)
		img.to_zero()
		try: img.add(float(parm[3]))
		except: pass
		try: n=int(parm[2])
		except: n=1
		if options.apix>0 : img["apix_x"]=img["apix_y"]=img["apix_z"]=options.apix
		for i in range(n): img.write_image("tmp_input_p2d.hdf",-1)
		args=["tmp_input_p2d.hdf"]

	if len(args) < 1: error_exit("e3proc2d.py --help  for usage information")

	### Validation of input files, and collecting statistics
	if options.output is None: error_exit("--output requred")

	nperfile=[]
	dimperfile=[]
	outperfile=[]
	for i,fsp in enumerate(args):
		if options.verbose>0 and len(args)>100:
			tlast=print_progress(tlast,"checking files",i+1,len(args))

		try: nperfile.append(file_image_count(fsp))
		except:
			error_exit(f"{fsp} does not appear to be a valid/accessible image file")

		hdr=get_header(fsp,0)
		if hdr["nz"]!=1 and base_name(fsp)[1]=="mrc": error_exit(f"ERROR: {fsp} is an MRC with 3-D data. If this is actually a stack file, please rename the file with a .mrcs extension")
		if hdr["nz"]!=1 : error_exit(f"ERROR: {fsp} contains 3-D data. Please use e3proc3d.py to work with 3-D files.")

		dimperfile.append((hdr["nx"],hdr["ny"],max(options.apix,hdr["apix_x"])))

		outperfile.append(options.output.format(n=i,base=base_name(fsp)[0]))

	if options.verbose>0: print(f"{sum(nperfile)} particles found in {len(args)} files")



	logid = E3init(sys.argv,options.ppid)

	tlast=0
	I=0
	# The actual data processing loop
	for FSPIN,N,DIM,FSPOUT in zip(args,nperfile,dimperfile,outperfile):
		I+=1
		if options.verbose>0 and len(args)>100:
			tlast=print_progress(tlast,"processing files",I-1,len(args))

		# if ":" in outfile:
		# 	outfile, bits, mode, rendermin_abs, rendermax_abs, rendermin_s, rendermax_s = parse_outfile_arg(outfile)

		chunksize=max(1000000000//(DIM[0]*DIM[1]*4),1)		# target no more than ~1 GB of input at a time
		for chunk in range(0,N,chunksize):
			if options.verbose>0 and len(args)<=100 :
				cp=sum(nperfile[:I-1])+chunk
				tlast=print_progress(tlast,f"{FSPIN} -> {FSPOUT} : {chunk}",cp,sum(nperfile))

			### read image chunk
			imgs=EMStack2D(EMData.read_images(FSPIN,range(chunk,min(chunk+chunksize,N))))
			hdrs=[im.get_attr_dict_cache() for im in imgs.emdata]		# do this fast before the data gets converted to JAX
			if options.apix>0:
				for hdr in hdrs: hdr["apix_x"]=hdr["apix_y"]=hdr["apix_z"]=options.apix

			if options.fouriershrink :
				fac=imgs.shape[1]/good_size(imgs.shape[1]/options.fouriershrink)
				imgs=imgs.downsample(good_size(imgs.shape[1]/options.fouriershrink)).do_ift()
				for hdr in hdrs:
					hdr["apix_x"]*=fac
					hdr["apix_y"]*=fac
					hdr["apix_z"]*=fac

			if options.normalize is not None:
				if options.normalize=="standard": imgs=imgs.normalize_standard()
				elif options.normalize=="edgemean": imgs=imgs.normalize_edgemean()

			for im,hd in zip(imgs.emdata,hdrs): im.set_attr_dict(hd)
			EMData.write_images(FSPOUT,imgs.emdata,chunk)		# there may be a bit specification in OUTFSP, so we use emdata.write_images directly

	if options.verbose>0 : tlast=print_progress(tlast,f"Processing complete",len(args),len(args))


	try: os.unlink("tmp_input_p2d.hdf")		# remove temporary file generated by --create
	except: pass
	E3end(logid)

def doparallel(argv,parallel,args):
	print("Parallelism not supported. Please use e2proc2dpar.py")

if __name__ == "__main__":
	main()
