#!/usr/bin/env python
#
# Author: T. Durmaz 08/29/2014 (tunay.durmaz@uth.tmc.edu)
# Copyright (c) 2014 The University of Texas - Houston Medical School
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

import os, sys
import json

from optparse import *
# from argparse import *
from EMAN2 import *
from EMAN2db import *
from EMAN2jsondb import *
# from emboxerbase import *

from sparx import *

"""
This program is used to window particles from a micrograph. The coordinates of the particles are given as input.
"""

def baseroot(path):
	return os.path.splitext(os.path.basename(path))[0]

def build_micnames(options, args):
	# 	Build micrograph basename list
	extension = options.micsuffix
	if len(args) > 0:
		micnames = args
	else:
		import glob
		micnames = glob.glob(os.path.join(options.indir, options.nameroot + "*" + extension))

	micnames = [ baseroot(m) for m in micnames]

	return micnames

def check_options(options, progname):
	if options.outdir == None:
		print "\nOutput directory must be specified with option --outdir. Type %s -h for help.\n" % progname
		sys.exit()
		
	if  options.coords_dir == None:
		print "\nCoordinates directory must be specified with option --coords_dir. Type %s -h for help.\n" % progname
		sys.exit()

	if  options.coords_extension == None:
		print "\nExtension of coordinates file must be specified with option --coords_extension. Type %s -h for help.\n" % progname
		sys.exit()

	if options.coords_format == None:
		print "\nCoordinate file format must be specified with option --coords_format. Type %s -h for help.\n" % progname
		sys.exit()
		
	if not(options.coords_format.lower() == 'sparx' or options.coords_format.lower() == 'eman1' or options.coords_format.lower() == 'eman2' or options.coords_format.lower() == 'spider'):
		print "\nInvalid option value: --coords_format=%s. Type %s -h for help.\n" % (options.coords_format.lower(), progname)
		sys.exit()

	if options.limitctf == True:
		if options.importctf == None:
			print "\nCTF parameters (--importctf) must be specified with option --limitctf. Type %s -h for help.\n" % progname
			sys.exit()
		
	if (options.resample_ratio <= 0.0 or options.resample_ratio > 1.0):
		print "\nInvalid option value: --resample_ratio=%s. Type %s -h for help.\n" % (options.resample_ratio, progname)
		sys.exit()
	
def main(args):
# 	parser1 = argparse.ArgumentParser(description='This program is used to window particles from a micrograph. The coordinates of the particles are given as input.')
# 	parser1.add_argument()

	progname = os.path.basename(sys.argv[0])
	usage = progname + " micrographs_list  --coords_dir=coords_dir  --coords_suffix=coords_suffix" + \
	                                          "  --coords_extension=coords_extension  --coords_format=coords_format" + \
	                                          "  --indir=input_dir  --importctf=ctf_file  --limitctf" + \
	                                          "  --resample_ratio=resample_ratio  --box_size=box_size" + \
	                                          "  --outdir=outdir  --outsuffix=outsuffix  --micsuffix=micsuffix" + \
	                                          "  --nameroot=nameroot  --invert" + \
	                                          "  --defocuserror=defocuserror  --astigmatismerror=astigmatismerror"
	parser = OptionParser(usage, version=SPARXVERSION)

	parser.add_option("--coords_dir",       type="string",        default=".",       help="<Coordinates Directory> Directory containing files with particle coordinates. (Default: current directory)")
	parser.add_option("--coords_suffix",    type="string",        default="",        help="<Coordinates File Suffix> Suffix of coordinate files. For example '_ptcls'. ")
	parser.add_option("--coords_extension", type="string",        default="box",     help="<Coordinates File Extension> File extension of coordinate files. e.g 'box' for eman1, 'json' for eman2, ...")
	parser.add_option("--coords_format",    type="string",        default="eman1",   help="<Coordinates File Format> Format of coordinates file: 'sparx', 'eman1', 'eman2', or 'spider'. The coordinates of sparx, eman2, and spider format is particle center. The coordinates of eman1 format is particle box conner associated with the original box size.")
	parser.add_option("--indir",            type="string",        default=".",       help="<Micrograph Directory> Directory containing micrographs to be processed. (Default: current directory)")
	parser.add_option("--nameroot",         type="string",        default="",        help="<Micrograph Root Name> Root name (Prefix) of micrographs to be processed.")
	parser.add_option("--micsuffix",        type="string",        default="hdf",     help="<Micrograph Extension > A string denoting micrograph type. (Default 'hdf')")
	parser.add_option("--outdir",           type="string",        default=".",       help="<Output Directory> Output directory (Default: current directory)")
	parser.add_option("--outsuffix",        type="string",        default="_ptcls",  help="<Output File Suffix> Suffix for output stack. (Default '_ptcls')")
	parser.add_option("--importctf",        type="string",        default="",        help="<CTER CTF File> File name with CTF parameters produced by sxcter.") 
	parser.add_option("--box_size",         type="int",           default=256,       help="<Box Size> x and y dimension in pixels of square area to be windowed. Pixel size after resampling is assumed when resample_ratio < 1.0 (Default 256)")
	parser.add_option("--invert",           action="store_true",  default=False,     help="<Invert Contrast> Invert image contrast (recommended for cryo data) (Default, no contrast inversion)")
	parser.add_option("--resample_ratio",   type="float",         default=1.0,       help="<Resample Ratio> Ratio of new to old image size (or old to new pixel size) for resampling. Valid range is 0.0 < resample_ratio <= 1.0. (Default: 1.0)  (advanced)")
	parser.add_option("--limitctf",         action="store_true",  default=False,     help="<Apply CTF-Limit Filter> Filter micrographs based on the CTF limit. It requires --importctf. (Default: no filter) (advanced)")	
	parser.add_option("--defocuserror",     type="float",         default=1000000.0, help="<Defocus Error Limit> Exclude micrographs whose relative defocus error as estimated by sxcter is larger than defocuserror percent.  The error is computed as (std dev defocus)/defocus*100%. (Default: include all irrespective of error values.) (advanced)" )
	parser.add_option("--astigmatismerror", type="float",         default=360.0,     help="<Astigmatism Error Limit> Set to zero astigmatism for micrographs whose astigmatism angular error as estimated by sxcter is larger than astigmatismerror degrees. (Default: include all irrespective of error values.)  (advanced)")
	
	# must be switched off in production
	# parser.add_option("--use_latest_master_directory", action="store_true", dest="use_latest_master_directory", default=False)
	# 
	# parser.add_option("--restart_section", type="string", default="", help="<restart section name> (no spaces) followed immediately by comma, followed immediately by generation to restart, example: \n--restart_section=candidate_class_averages,1         (Sections: restart, candidate_class_averages, reproducible_class_averages)")
	# parser.add_option("--stop_after_candidates",          action="store_true", default=False,   help="<stop_after_candidates> stops after the 'candidate_class_averages' section")
	# 
	parser.add_option("--return_options", action="store_true", dest="return_options", default=False, help = SUPPRESS_HELP)

	(options, args) = parser.parse_args(args)

	if options.return_options:
		return parser
	
# 	Set local constants
	box_size = options.box_size
	box_half = box_size // 2
	options.micsuffix = "." + options.micsuffix
	cterr = [options.defocuserror/100.0, options.astigmatismerror]
		
	check_options(options, progname)
	
	extension_coord = options.coords_suffix + "." + options.coords_extension
	
# 	Build micrograph basename list
	micnames = build_micnames(options, args)
	print_msg('Detected micrographs : %6d ...\n' % (len(micnames)))
	
# 	If there is no micrographs, exit
	if len(micnames) == 0:
		print usage
		sys.exit()
	
# 	Load CTFs
	n_reject_defocus_error = 0
	if options.importctf != None:
		ctfs0 = read_text_row(options.importctf)
		print_msg('Detected CTF entries : %6d ...\n' % (len(ctfs0)))

		ctfs={}
		for i in xrange(len(ctfs0)):
			ctf=ctfs0[i]
			basemic = baseroot(ctf[-1])
			if(ctf[8]/ctf[0] > cterr[0]):
				print_msg('Defocus error %f exceeds the threshold. Micrograph %s rejected.\n' % (ctf[8]/ctf[0], basemic))
				n_reject_defocus_error += 1
			else:
				if(ctf[10] > cterr[1] ):
					ctf[6] = 0.0
					ctf[7] = 0.0
				ctfs[basemic] = ctf
		print_msg('Rejected micrographs by defocus error  : %6d ...\n' % (n_reject_defocus_error))
	
# 	Create circular 2D mask for ...
	mask = model_circle(box_size//2, box_size, box_size)

# 	Prepare loop variables
	n_micrographs_process = 0
	n_micrographs_reject_no_micrograph = 0
	n_micrographs_reject_no_coordinates = 0
	n_micrographs_reject_no_cter_entry = 0
	n_total_coordinates_detect = 0
	n_total_coordinates_process = 0
	n_total_coordinates_reject_out_of_boundary = 0
	cutoffhistogram = []		#@ming compute the histogram for micrographs cut of by ctf limit.
# 	Loop over micrographs
	for k in range(len(micnames)):
		# basename is name of micrograph minus the path and extension
		# Here, assuming micrograph and coordinates have the same file basename
		basename = micnames[k]
		f_mic    = os.path.join(os.path.abspath(options.indir), basename + options.micsuffix)
		f_info   = os.path.join(options.coords_dir, basename + extension_coord)

# 		CHECKS: BEGIN
# 		IF micrograph exists
		if not os.path.exists(f_mic):
			print_msg('    Cannot read %s. Skipping %s ...\n' % (f_mic, basename))
			n_micrographs_reject_no_micrograph += 1
			continue
		
# 		IF coordinates file exists
		if not os.path.exists(f_info):
			print_msg('    Cannot read %s. Skipping %s ...\n' % (f_info, basename))
			n_micrographs_reject_no_coordinates += 1
			continue
		
# 		IF micrograph is in CTER results
		if options.importctf != None:
			if basename not in ctfs:
				print_msg('    Is not listed in CTER results, skipping %s...\n' % (basename))
				n_micrographs_reject_no_cter_entry += 1
				continue
			else:
				ctf = ctfs[basename]		
# 		CHECKS: END

		n_micrographs_process += 1
		
		print_msg('\n')
		print_msg('Processing micrograph %s... Path: %s... Coordinates file %s\n' % (basename, f_mic, f_info))
	
# 		Read coordinates according to the specified format and 
# 		make the coordinates the center of particle image 
		if options.coords_format.lower() == 'sparx' :
			coords = read_text_row(f_info)
		elif options.coords_format.lower() == 'eman1':
			coords = read_text_row(f_info)
			for i in range(len(coords)):
				coords[i] = [coords[i][0] + coords[i][2]//2  ,coords[i][1] + coords[i][3]//2]
		elif options.coords_format.lower() == 'eman2':
			coords = js_open_dict(f_info)["boxes"]
			for i in range(len(coords)):
				coords[i] = [coords[i][0],coords[i][1]]
		elif options.coords_format.lower() == 'spider':
			coords = read_text_row(f_info)
			for i in range(len(coords)):
				coords[i] = [coords[i][2] ,coords[i][3]]
		else:
			assert(False) # Unreachable code
		
# 		Load micrograph from the file
		immic = get_im(f_mic)
		
# 		Calculate the new pixel size
		resample_ratio = options.resample_ratio
		if options.importctf != None:		
			pixel_size_orig = ctf[3]
			
			if resample_ratio < 1.0:
				assert(resample_ratio > 0.0)
				new_pixel_size = pixel_size_orig / resample_ratio
				print_msg('Resample micrograph to pixel size %6.4f and window segments from resampled micrograph\n' % new_pixel_size)
			else:
				# assert(resample_ratio == 1.0)
				new_pixel_size = pixel_size_orig
		
# 			Set ctf along with new pixel size in resampled micrograph
			ctf[3] = new_pixel_size	
		else:
			assert(options.importctf == None)
			if resample_ratio < 1.0:
				assert(resample_ratio > 0.0)
				print_msg('Resample micrograph with ratio %6.4f and window segments from resampled micrograph\n' % resample_ratio)
			# else:
			#	assert(resample_ratio == 1.0)
			
# 		Apply filters to micrograph
		fftip(immic)
		if options.limitctf:
			assert(options.importctf != None)
# 			Cut off frequency components higher than CTF limit 
			q1, q2 = ctflimit(box_size,ctf[0],ctf[1],ctf[2],new_pixel_size)
			
# 			This is absolute frequency of the CTF limit in the scale of original micrograph
			if resample_ratio < 1.0:
				assert(resample_ratio > 0.0)
				q1 = resample_ratio * q1 / float(box_size) # q1 = (pixel_size_orig / new_pixel_size) * q1/float(box_size)
			else:
				# assert(resample_ratio == 1.0) -> pixel_size_orig == new_pixel_size -> pixel_size_orig / new_pixel_size == 1.0
				q1 = q1 / float(box_size)
			
			if q1 < 0.5:          #@ming
				immic = filt_tanl(immic, q1, 0.01)
				cutoffhistogram.append(q1)
		
# 		Cut off frequency components lower than the box size can express 
		immic = fft(filt_gaussh( immic, resample_ratio/box_size ))
		
# 		Resample micrograph, map coordinates, and window segments from resampled micrograph using new coordinates
# 		after resampling by resample_ratio, new pixel size will be pixel_size/resample_ratio = new_pixel_size
#		NOTE: 2015/04/13 Toshio Moriya
#		resample() efficiently takes care of the case resample_ratio = 1.0 but
#		it does not set apix_*. Even though it sets apix_* when resample_ratio < 1.0 ...
		immic = resample(immic, resample_ratio)
				
		if options.invert:
			stt = Util.infomask(immic, None, True)
			Util.mul_scalar(immic, -1.0)
			immic += 2*stt[0]
		
		if options.importctf != None:
			from utilities import generate_ctf
			ctf = generate_ctf(ctf)

# 		Prepare loop variables
		nx = immic.get_xsize() 
		ny = immic.get_ysize()
		x0 = nx//2
		y0 = ny//2
		print_msg('\n')
		print_msg('Micrograph size := (%6d, %6d)\n' % (nx, ny))

		otcl_images  = "bdb:%s/" % options.outdir + basename + options.outsuffix
		ind = 0
		
		n_coordinates_reject_out_of_boundary = 0
		
# 		Loop over coordinates
		for i in range(len(coords)):
		
			source_x = int(coords[i][0])
			source_y = int(coords[i][1])
					
			x = source_x		
			y = source_y
			
			if resample_ratio < 1.0:
				assert(resample_ratio > 0.0)
				x = int(x * resample_ratio)	
				y = int(y * resample_ratio)
			# else:
			# 	assert(resample_ratio == 1.0)
				
			if( (0 <= x - box_half) and ( x + box_half <= nx ) and (0 <= y - box_half) and ( y + box_half <= ny ) ):
				imw = Util.window(immic, box_size, box_size, 1, x-x0, y-y0)
			else:
				print_msg('Coordinates ID = %04d (x = %4d, y = %4d, box_size = %4d) is out of micrograph bound, skipping ....\n' % (i, x, y, box_size))
				n_coordinates_reject_out_of_boundary += 1
				continue
			
			imw = ramp(imw)
			stat = Util.infomask( imw, mask, False )
			imw -= stat[0]
			imw /= stat[1]

#			NOTE: 2015/04/09 Toshio Moriya
#		    ptcl_source_image might be redundant information ...
#		    Consider re-organizing header entries...
			imw.set_attr("ptcl_source_image", f_mic)
			imw.set_attr("ptcl_source_coord_id", i)
			imw.set_attr("ptcl_source_coord", [source_x, source_y])
			imw.set_attr("resample_ratio", resample_ratio)
						
#			NOTE: 2015/04/13 Toshio Moriya
#			apix_* attributes are updated by resample() only when resample_ratio != 1.0
# 			Let's make sure header info is consistent by setting apix_* = 1.0 
#			regardless of options, so it is not passed down the processing line
			imw.set_attr("apix_x", 1.0)
			imw.set_attr("apix_y", 1.0)
			imw.set_attr("apix_z", 1.0)			
			if options.importctf != None:
				imw.set_attr("ctf",ctf)
				imw.set_attr("ctf_applied", 0)
				imw.set_attr("pixel_size_orig", pixel_size_orig)
				# imw.set_attr("apix_x", new_pixel_size)
				# imw.set_attr("apix_y", new_pixel_size)
				# imw.set_attr("apix_z", new_pixel_size)
#			NOTE: 2015/04/13 Toshio Moriya 
#			Pawel Comment: Micrograph is not supposed to have CTF header info.
#			So, let's assume it does not exist & ignore its presence.
#           Note that resample() "correctly" updates pixel size of CTF header info if it exists
			# elif (imw.has_ctff()):
			# 	assert(options.importctf == None)
			# 	ctf_origin = imw.get_attr("ctf")
			# 	pixel_size_origin = round(ctf_origin.apix, 5) # Because SXCTER ouputs up to 5 digits 
			# 	imw.set_attr("apix_x",pixel_size_origin)
			# 	imw.set_attr("apix_y",pixel_size_origin)
			# 	imw.set_attr("apix_z",pixel_size_origin)	
			
			imw.write_image(otcl_images, ind)
			ind += 1
		
		n_total_coordinates_detect += len(coords)
		n_total_coordinates_process += ind
		n_total_coordinates_reject_out_of_boundary += n_coordinates_reject_out_of_boundary
		
#		Print out the summary of this micrograph
		print_msg('\n')
		print_msg('Micrograph summary of coordinates...\n')
		print_msg('Detected                        : %4d\n' % (len(coords)))
		print_msg('Processed                       : %4d\n' % (ind))
		print_msg('Rejected by out of boundary     : %4d\n' % (n_coordinates_reject_out_of_boundary))
	
	if options.limitctf:
#		Print out the summary of CTF-limit filtering
		print_msg('\n')
		print_msg('Global summary of CTF-limit filtering (--limitctf) ...\n')
		print_msg('Percentage of filtered micrographs: %8.2f\n' % (len(cutoffhistogram) * 100.0 / len(micnames)))
		
		lhist = 10
		if len(cutoffhistogram) >= lhist:
			from statistics import hist_list
			region,hist = hist_list(cutoffhistogram, lhist)	
			print_msg("      Histogram of cut off frequency\n")
			print_msg("      ERROR       number of frequencies\n")
			for lhx in xrange(lhist):
				print_msg(" %14.7f     %7d\n" % (region[lhx], hist[lhx]))  # print_msg(" %10.3f     %7d\n" % (region[lhx], hist[lhx]))  
		else:
			print_msg("The number of filtered micrographs (%d) is less than the number of bins (%d). No histogram is produced.\n" % (len(cutoffhistogram), lhist))
		
#	Print out the summary of all micrographs
	print_msg('\n')
	print_msg('Global summary of micrographs ...\n')
	print_msg('Detected                        : %6d\n' % (len(micnames)))
	print_msg('Processed                       : %6d\n' % (n_micrographs_process))
	print_msg('Rejected by no micrograph file  : %6d\n' % (n_micrographs_reject_no_micrograph))
	print_msg('Rejected by no coordinates file : %6d\n' % (n_micrographs_reject_no_coordinates))
	print_msg('Rejected by no CTER entry       : %6d\n' % (n_micrographs_reject_no_cter_entry))
	print_msg('\n')
	print_msg('Global summary of coordinates ...\n')
	print_msg('Detected                        : %6d\n' % (n_total_coordinates_detect))
	print_msg('Processed                       : %6d\n' % (n_total_coordinates_process))
	print_msg('Rejected by out of boundary     : %6d\n' % (n_total_coordinates_reject_out_of_boundary))
	print_msg('\n')

						
if __name__=='__main__':
	main(sys.argv[1:])
