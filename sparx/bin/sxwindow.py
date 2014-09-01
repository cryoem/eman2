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
from filter import *
from global_def import *


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

# 	parser.add_option('--coords_dir',   dest='coordsdir',                help='Directory containing particle coordinates')
# 	parser.add_option('--importctf',    dest='ctffile',                  help='File name with CTF parameters produced by sxcter.')
# 	parser.add_option('--topdir',       dest='topdir',       default='', help='Path name of directory containing relevant micrograph directories')
# 	parser.add_option('--input_pixel',  dest='input_pixel',  default=1,  help='input pixel size')
# 	parser.add_option('--output_pixel', dest='output_pixel', default=1,  help='output pixel size')

	(options, args) = parser.parse_args()
	
	if len(args) < 1:
		print "\nusage: " + usage
		print "Please run '" + progname + " -h' for detailed options\n"
	else:
		params = {}  # In this case you do not need dictionary.  It just confuses the code.
		# the next line seems to be incorrect.  args is a list of arguments, the first one is the name of the program itself
		#  In addition, the way you are trying to do it the user would have to give the list of actual name
		#  Please check how this is done in sxhelixoboxer!!
		params["filenames"] = args
		params["format"] = 'json'
# 		params["coordsdir"] = options.coordsdir
# 		params["ctffile"]   = options.ctffile
		database = "e2boxercache"
		db = js_open_dict(os.path.join(database,"quality.json"))
		box_size = # box size is a user-defined parameter

		for f in params["filenames"]:
			im = EMData(f)
			x0=im.get_xsize()/2
			y0=im.get_ysize()/2

			coords = js_open_dict(info_name(f))["boxes"]
			for i in range(len(coords)):
				x = coords[i][0]
				y = coords[i][1]
				imn=Util.window(im, box_size, box_size, 1, int(x-x0),int(y-y0))
				# next line - why str?  the variables you use are laready strings.
				#  where suffix and extension are defined?
				imn.write_image(str(base_name(f) + db['suffix'] + db['extension']), -1) # -1: appending to the image stack


if __name__=='__main__':
	main()
