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
from EMAN2 import *
from EMAN2db import *
from EMAN2jsondb import *
from emboxerbase import *

from utilities import *
from fundamentals import *
from filter import *
from global_def import *
from distutils.extension import Extension

def get_suffix_and_extension(db_dir, db_file):
	db = js_open_dict(os.path.join(db_dir, db_file))
	suffix    = str(db['suffix'])    # db['suffix'] is unicode, so need to call str()
	extension = str(db['extension']) # db['extension'] is unicode, so need to call str()
	
	return suffix, extension

def get_mic_base_names(options):
	micnames = []
	import glob
	
	for f in glob.glob(os.path.join(options.topdir, '*.hdf')):   # currently handles only hdf formatted micrographs
			micnames.append(base_name(f))
	return micnames

def ctf(options):
	"""
	Read ctf information.
	"""

	ctfs = read_text_row(options.importctf)
	cterr = [options.defocuserror/100.0, options.astigmatismerror]

	for i in xrange(len(ctfs)):
		smic = ctfs[i][-1].split('/')
		ctfilename = (smic[-1].split('.'))[0]
# 		if(ctfs[i][8]/ctfs[i][0] > cterr[0]):
# 			print_msg('Defocus error %f exceeds the threshold. Micrograph %s rejected.\n'%(ctfs[i][8]/ctfs[i][0], ctfilename))
# 			ctfs[i][0]=10.0
# 			continue
		if(ctfs[i][10] > cterr[1] ):
			ctfs[i][6] = 0.0
			ctfs[i][7] = 0.0
		
	return ctfs

def window_micrograph(options, basename, ctf, suffix, extension):
	f_mic = os.path.join(options.topdir, basename + extension)
	f_info = info_name(f_mic)
			
	otcl_images  = "bdb:%s/"%options.outdir + basename + suffix
# 	print basename, f_info, f_mic

	box_size = options.box_size
	mask = pad(model_circle(box_size//2, box_size, box_size), box_size, box_size, 1, 0.0)

	im = get_im(f_mic)
	x0 = im.get_xsize()//2  #  Floor division or integer division
	y0 = im.get_ysize()//2
		
	coords = js_open_dict(f_info)["boxes"]
	for j in range(len(coords)):

		x = int(coords[j][0])
		y = int(coords[j][1])

		imn=Util.window(im, box_size, box_size, 1, x-x0, y-y0)
		imn.set_attr('ptcl_source_image',f_mic)
		imn.set_attr('ptcl_source_coord',[x,y])
		stat = Util.infomask(imn, mask, False)   

		imn = ramp(imn)
		imn -= stat[0]
		Util.mul_scalar(imn, 1.0/stat[1])
			
		ctff=generate_ctf(ctf)
# 		print ctff
		imn.set_attr("ctf",ctff)
		imn.set_attr("ctf_applied", 0)
		imn.write_image(otcl_images, j)
# 		imn.write_image(otcl_images, iImg)
# 		iImg = iImg + 1

def window(data):
	"""
	Using coordinates window particles, and add ctf information to it.
	"""
	for k, info in data.items():
		print 'Processing {0}'.format(k)
		box_size = data[k]['box_size']
		pixel_ratio = float(data[k]['input_pixel'])/float(data[k]['output_pixel'])
		img = EMData()
		img.read_image(k)
		img_filt = filt_gaussh(img, pixel_ratio/box_size)


		if pixel_ratio != 1.0:
			print "Generating downsampled image\n"
			sb = Util.sincBlackman(15, .5 * pixel_ratio,1999) # 1999 taken directly from util_sparx.h
			img_filt = img_filt.downsample(sb, pixel_ratio)
			box_size = box_size / pixel_ratio

		output_file_name = 'out_' + os.path.basename(k)
		clip = EMData()
		for i, (x, y) in enumerate(data[k]['coordinates']):
			reg = Region((x * pixel_ratio)-box_size//2, (y * pixel_ratio)-box_size//2, box_size, box_size)
			clip = img_filt.get_clip(reg)
			clip.write_image(output_file_name, i)
		# Set ctf
		set_ctf(clip, data[k]['ctf'])
		print 'Windowed prticles for {0} -> {1}'.format(k, output_file_name)


def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " --coords_dir=coords_dir  --importctf=ctf_file  --topdir=topdir  --input_pixel=input_pixel  --output_pixel=output_pixel"
	
	parser = OptionParser(usage, version=SPARXVERSION)

	parser.add_option('--coords_dir',   dest='coordsdir',                help='Directory containing particle coordinates')
# 	parser.add_option('--importctf',    dest='ctffile',                  help='File name with CTF parameters produced by sxcter.')
	parser.add_option('--topdir',       dest='topdir',       default='./', help='Path name of directory containing relevant micrograph directories')
# 	parser.add_option('--input_pixel',  dest='input_pixel',  default=1,  help='input pixel size')
# 	parser.add_option('--output_pixel', dest='output_pixel', default=1,  help='output pixel size')
	parser.add_option('--box_size',     dest='box_size',     type=int,   help='box size')
	parser.add_option('--outdir',     dest='outdir',      help='Output directory')
	parser.add_option('--outstack',     dest='outstack',      help='Output stack name')

	# import ctf estimates done using cter
	parser.add_option("--input",              type="string",	default= None,     		  help="Input particles.")
	parser.add_option("--importctf",          type="string",	default= None,     		  help="Name of the file containing CTF parameters produced by sxcter.")
	parser.add_option("--defocuserror",       type="float",  	default=1000000.0,        help="Exclude micrographs whose relative defocus error as estimated by sxcter is larger than defocuserror percent.  The error is computed as (std dev defocus)/defocus*100%")
	parser.add_option("--astigmatismerror",   type="float",  	default=360.0,            help="Set to zero astigmatism for micrographs whose astigmatism angular error as estimated by sxcter is larger than astigmatismerror degrees.")

	(options, args) = parser.parse_args()
	
	if len(args) > 0:
		print "\nusage: " + usage
		print "Please run '" + progname + " -h' for detailed options\n"
	else:
		ctfs = ctf(options)
		
# 		otcl_images  = "bdb:%s/"%options.outdir + options.outstack + suffix
# 		iImg=0
		suffix, extension = get_suffix_and_extension("e2boxercache","quality.json")
		micnames = get_mic_base_names(options)
# 		print micnames		
		for i in range(len(micnames)):
			window_micrograph(options, micnames[i], ctfs[i], suffix, extension)
# 			basename = micnames[i]
# 			f_mic = os.path.join(options.topdir, basename + extension)
# 			f_info = info_name(f_mic)
# 			
# 			otcl_images  = "bdb:%s/"%options.outdir + basename + suffix
# # 			print basename, f_info, f_mic
# 				
# 			im = get_im(f_mic)
# 			x0 = im.get_xsize()//2  #  Floor division or integer division
# 			y0 = im.get_ysize()//2
# 				
# 			coords = js_open_dict(f_info)["boxes"]
# 			for j in range(len(coords)):
# 
# 				x = int(coords[j][0])
# 				y = int(coords[j][1])
# 
# 				imn=Util.window(im, box_size, box_size, 1, x-x0, y-y0)
# 				imn.set_attr('ptcl_source_image',f_mic)
# 				imn.set_attr('ptcl_source_coord',[x,y])
# 				stat = Util.infomask(imn, mask, False)   
# 
# 				imn = ramp(imn)
# 				imn -= stat[0]
# 				Util.mul_scalar(imn, 1.0/stat[1])
# 					
# 				ctff=generate_ctf(ctfs[i])
# # 				print ctff
# 				imn.set_attr("ctf",ctff)
# 				imn.set_attr("ctf_applied", 0)
# 				imn.write_image(otcl_images, j)
# # 				imn.write_image(otcl_images, iImg)
# # 				iImg = iImg + 1

if __name__=='__main__':
	main()
