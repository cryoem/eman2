#!/usr/bin/env python

#
# Author: Ross Coleman (racolema@gmail.com)
# Copyright (c) 2009- Baylor College of Medicine
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
#

from EMAN2 import get_image_directory, Transform, Region, EMANVERSION, EMData, E2init, E2end, EMArgumentParser
from EMAN2db import db_open_dict, db_check_dict, db_close_dict
from math import *
import sys
import os

try:
	from PyQt4 import QtGui, QtCore
	from emapplication import EMApp, get_application
	from emimage2d import EMImage2DWidget
	from emselector import EMSelectorDialog
	from emshape import EMShape, EMShapeDict
	
	ENABLE_GUI = True
	
except ImportError, e:
	print "Importing GUI libraries failed!"
	print e
	print "GUI features are disabled."
	ENABLE_GUI = False

"""
This program is used to box helices from a micrograph and extract segments from the helices.
A helical region of a macromolecule can appear as a rectangle on a 2D micrograph; thus a "helix" in this program is a rectangle. 
Different helices from the same micrograph will generally have different lengths but the same width.
A "particle" in this program is a square or rectangular region taken from a "helix." Particles are chosen so that they overlap each other
within a helix. Usually, all segments from a micrograph will have the same dimensions.
"""

SXHELIXBOXER_DB = "bdb:" # used to use "bdb:e2helixboxercache#"


def main():
	usage = """sxhelixboxer.py --gui <micrograph1> <<micrograph2> <micrograph3> ...
	sxhelixboxer.py --gui --helix-width=<width> <micrograph1> <<micrograph2> <micrograph3> ...
	sxhelixboxer.py <options (not --gui)> <micrograph>    
	sxhelixboxer.py <outdir> --window --dirid=mic --micid=mymic --micsuffix=hdf --invert_contrast --apix=1.84 --boxsize='200,45' --hcoords_suffix=_boxes.txt --ptcl-dst=30 --rmax=92.0 --importctf=ctfestimates.txt

	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--gui", action="store_true", help="Start the graphic user interface for boxing helices.")

	parser.add_argument("--helix-coords", "-X", type=str, help="Save coordinates for helices to the file specified, which will have the EMAN1 *.box format:\t\t\tx1-w/2        y1-w/2        w        w        -1                    x2-w/2        y2-w/2        w        w        -2")
	parser.add_argument("--helix-images", "-x", type=str, help="Save images of the helices. The file name specified will have helix numbers added to it.")
	parser.add_argument("--ptcl-coords",  "-P", type=str, help="Save coordinates of the centers of segments to the specified formatted text file")
	parser.add_argument("--ptcl-images",  "-p", type=str, help="Save images of the segments. The file name specified will have helix numbers (and segment numbers if the file type does not support image stacks) added to it.")
	parser.add_argument("--ptcl-images-stack-mode", type=str, default="multiple", help="Options for saving segment images to stack files. 'single' uses one stack file, 'multiple' (default) uses one stack file per helix, 'none' uses a file for each segment and is always used when the output file format does not support image stacks.")
	
	parser.add_argument("--db-add-hcoords",     type=str, help="Append any unique helix coordinates to the database from the specified file (in EMAN1 *.box format). Use --helix-width to specify a width for all boxes.")
	parser.add_argument("--db-set-hcoords",     type=str, help="Replaces the helix coordinates in the database with the coordinates from the specified file (in EMAN1 *.box format). Use --helix-width to specify a width for all boxes.")
	parser.add_argument("--ptcl-dst", 	        type=int, 	         dest="ptcl_dst", 			  help="Distance between windowed squares in pixels", default=-1)
	parser.add_argument("--helix-width", "-w",  type=int,            dest="helix_width", help="Helix width in pixels. Overrides widths saved in the database or in an input file.", default=-1)
	parser.add_argument("--ptcl-length",        type=int,            dest="ptcl_length", help="Segment length in pixels", default=-1)
	parser.add_argument("--ptcl-width",         type=int,            dest="ptcl_width", help="Segment width in pixels", default=-1)
	parser.add_argument("--ptcl-not-rotated",   action="store_true", dest="ptcl_not_rotated", help="Segments are oriented as on the micrograph. They are square with length max(ptcl_length, ptcl_width).")
	parser.add_argument("--ptcl-norm-edge-mean",action="store_true", help="Apply the normalize.edgemean processor to each segment.")
	parser.add_argument("--gridding",           action="store_true", default=False, help="Use a gridding method for rotation operations on segments. Requires segments to be square. (Used by default. If flag present, gridding will be turned off)")
	parser.add_argument("--save-ext",           type=str,default="hdf",dest="save_ext",help="The default file extension to use when saving 'segment' images. This is simply a convenience for improved workflow. If a format other than HDF is used, metadata will be lost when saving.")
	parser.add_argument("--ppid",               type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	parser.add_argument("--invert_contrast",     action="store_true", default=False, help="Invert contrast of micrograph of contrast of boxed filaments and segments are inverted")
	
	# window segments from boxed filaments 
	parser.add_argument("--window",             action="store_true",	 default=False,               help="window segments from boxed filaments (or helices in sxhelixboxer lingo)")
	parser.add_argument("--dirid",              type=str,				 default="",                  help="A string for identifying directories containing relevant micrographs. Any directory containing dirid as a contiguous string will be searched for micrographs. These micrographs are assumed to be those which were used to box the helices which are to be windowed.")
	parser.add_argument("--micid",              type=str,				 default="",                  help="A string for identifying the name (minus extension) of relevant micrographs.")
	parser.add_argument("--micsuffix",          type=str,				 default="",                  help="A string denoting micrograph type. Currently only handles suffix types, e.g. 'hdf', 'ser' etc.")
	parser.add_argument("--boxsize",            type=int,				 default=256,            	  help="x and y dimension in pixels of square area to be windowed out from filament. Pixel size is assumed to be new_pixel_size.")
	#parser.add_argument("--outstacknameall",    type=str,				 default="bdb:data",          help="File name plus path and type (only handles bdb and hdf right now) under which ALL windowed segments from ALL micrograph directories will be saved, e.g. 'bdb:adata' or 'adata.hdf'")
	parser.add_argument("--hcoords_dir",        type=str,				 default="",        		  help="Directory containing helix box coordinates")
	parser.add_argument("--hcoords_suffix",     type=str,				 default="_boxes.txt",        help="String identifier which when concatenated with a micrograph name (minus extension) gives the name of the text file containing coordinates of ALL helices boxed from the micrograph. If there is no such file, helices boxed from the micrograph will not be windowed. Default is '_boxes.txt', so if mic0.hdf is a micrograph, then the text file containing coordinates of the helices boxed in it is mic0_boxes.txt. The coordinate file is assumed to be in the format used by sxhelixboxer.")
	parser.add_argument("--new_apix",           type=float,  			 default=-1.0,                help="New target pixel size to which the micrograph should be resampled. Default is -1, in which case there is no resampling.")
	parser.add_argument("--freq",               type=float, 			 default=-1.0,                help="Cut-off frequency at which to high-pass filter micrographs before windowing. Default is -1, in which case, the micrographs will be high-pass filtered with cut-off frequency 1.0/segnx, where segnx is the target x dimension of the segments.") 
	parser.add_argument("--apix",               type=float,				 default= -1.0,               help="pixel size in Angstroms")   
	parser.add_argument("--rmax",               type=float, 			 default= 80.0,               help="maximal radius for hsearch (Angstroms)")
	parser.add_argument("--minseg",             type=int,	 			 default=6,                   help="Skip filaments that yield fewer then minseg segments. Default is 6.")
	parser.add_argument("--dbg",                type=int,	 			 default=1,                   help="If 1, then intermediate files in output directory will not be deleted; if 0, then all output directories where intermediate files were stored will be deleted. Default is 1.")
	parser.add_argument("--topdir",             type=str,				 default="",                  help="Path name of directory containing relevant micrograph directories")
	
	# import ctf estimates done using cter
	parser.add_argument("--importctf",          type=str,				 default= None,     		  help="File name with CTF parameters produced by sxcter.")
	parser.add_argument("--defocuserror",       type=float,  			 default=1000000.0,           help="Exclude micrographs whose relative defocus error as estimated by sxcter is larger than defocuserror percent.  The error is computed as (std dev defocus)/defocus")
	parser.add_argument("--astigmatismerror",   type=float,  			 default=360.0,               help="Set to zero astigmatism for micrographs whose astigmatism angular error as estimated by sxcter is larger than astigmatismerror degrees.")
	parser.add_argument("--limitctf",             action="store_true",     default=False,               help="Filter micrographs based on the CTF limit. (Default: no filter)")

	(options, args) = parser.parse_args()
	
	# invert meaning of gridding
	options.gridding = not options.gridding
	# Window filaments 
	if options.window:
		if len(args) < 1:
			print "Must specify name of output directory where intermediate files are to be deposited."
			return
		outdir = args[0]
		tdir = options.topdir
		if len(tdir) < 1:  tdir = None
		if options.importctf:  cterr = [options.defocuserror/100.0, options.astigmatismerror]
		else:                  cterr = None
		options.outstacknameall = "bdb:data"  # this was disabled, to be removed later PAP
 		windowallmic(options.dirid, options.micid, options.micsuffix, outdir, pixel_size=options.apix, boxsize=options.boxsize, minseg=options.minseg,\
				outstacknameall=options.outstacknameall, hcoords_dir = options.hcoords_dir, hcoords_suffix = options.hcoords_suffix, ptcl_dst=options.ptcl_dst, \
				inv_contrast=options.invert_contrast, new_pixel_size=options.new_apix, rmax = options.rmax, freq=options.freq, \
				debug = options.dbg, do_rotation = True, do_gridding=options.gridding, topdir=tdir, importctf=options.importctf, limitctf=options.limitctf, cterr = cterr)
		return

	if options.helix_width < 1:
		helix_width = None
	else:
		helix_width = options.helix_width

	if options.ptcl_width < 1:
		px_width = None
	else:
		px_width = options.ptcl_width

	if options.ptcl_length < 1:
		px_length = None
	else:
		px_length = options.ptcl_length
	
	if options.ptcl_dst < 0:
		px_dst = None
	else:
		px_dst = options.ptcl_dst	
	
	logid=E2init(sys.argv)

	if options.gui:
		if ENABLE_GUI:
			logid=E2init(sys.argv,options.ppid)
			app = EMApp()
			helixboxer = EMHelixBoxerWidget(args, app, options.helix_width,options.save_ext, options.invert_contrast)
			helixboxer.show()
			app.execute()
			E2end(logid)
		else:
			return
	else:
		if len(args) == 1:
			logid=E2init(sys.argv)
			micrograph_filepath = args[0]
			if options.db_add_hcoords:
				db_load_helix_coords(micrograph_filepath, options.db_add_hcoords, True, helix_width)
			if options.db_set_hcoords:
				db_load_helix_coords(micrograph_filepath, options.db_set_hcoords, False, helix_width)
			if options.helix_coords:
				db_save_helix_coords(micrograph_filepath, options.helix_coords, helix_width)
			if options.helix_images:
				db_save_helices(micrograph_filepath, options.helix_images, helix_width)
			if options.ptcl_coords:
				db_save_particle_coords(micrograph_filepath, options.ptcl_coords, px_dst, px_length, px_width)
			if options.ptcl_images:
				db_save_particles(micrograph_filepath, options.ptcl_images, px_dst, px_length, px_width, not options.ptcl_not_rotated, options.ptcl_norm_edge_mean, options.gridding, options.ptcl_images_stack_mode)
				
			E2end(logid)
		elif len(args) == 0:
			print 'You must specify a micrograph file or use the "--gui" option.'
			return
		elif len(args) > 1:
			print 'Multiple micrographs can only be specified with the "--gui" option'
			return

def counterGen():
	"""
	Calling this function will create a counter object.
	Ex: "counter = counterGen()"
	Then calling "counter.next()" will return the next number in {1, 2, 3, ...}
	"""
	i = 0
	while True:
		i += 1
		yield i

def get_helix_from_coords(micrograph, x1, y1, x2, y2, width):
	"""
	Gets the rectangular helix of the image specified by the coordinates.
	@param micrograph: the EMData object which holds the micrograph from which helices and particles are chosen
	@param x1: x-coordinate in pixels of the first end-point along the long axis of symmetry of the rectangle
	@param y1: y-coordinate in pixels of the first end-point along the long axis of symmetry of the rectangle
	@param x2: x-coordinate in pixels of the second end-point along the long axis of symmetry of the rectangle
	@param y2: y-coordinate in pixels of the second end-point along the long axis of symmetry of the rectangle
	@param width: the width in pixels of the rectangle
	@return: the rectangular EMData helix specified by the coordinates and width 
	"""
	rot_angle = get_helix_rotation_angle(x1, y1, x2, y2, width)     ## rot_angle   @ming
	centroid = ( (x1+x2)/2.0,(y1+y2)/2.0 )
	l_vect = (x2-x1, y2-y1)
	length = sqrt(l_vect[0]**2+l_vect[1]**2)
	tr = Transform()
	tr.set_trans(centroid)
	tr.set_rotation({"type":"2d", "alpha":rot_angle})
	helix_dimensions = ( int(round(width)), int(round(length)), 1 )
	helix = micrograph.get_rotated_clip( tr, helix_dimensions )
	helix["ptcl_helix_coords"] = (x1, y1, x2, y2, width)
	helix["xform.align2d"]     = Transform()
	helix["xform.projection"]  = Transform()   
	helix["astig_jiao"]        = rot_angle                ##  add the rotation angle to astigmatism angle @ming
	return helix

def get_helix_rotation_angle(x1, y1, x2, y2, width):
	l_vect = (x2-x1, y2-y1)
	length = sqrt(l_vect[0]**2+l_vect[1]**2)
	assert length != 0
	l_uvect = (l_vect[0]/length, l_vect[1]/length)
	
	#Rotate so that the length is parallel to the y-axis
	#Angle between l_uvect and y-axis: l_uvect (dot) j_hat = cos (rot_angle)
	rot_angle = degrees(acos( l_uvect[1]))       #180/pi*acos( l_uvect[1] )                  ## positive angle for counterclockwise rotation @ming
	#Whether we rotate clockwise or counterclockwise depends on the sign of l_uvect[0] (the x-component)
	if l_uvect[0] < 0:
		rot_angle *= -1 #rotate the box clockwise
	return rot_angle
	
def get_particle_centroids(helix_coords, px_overlap, px_length, px_width, is_rotated = False):
	"""
	Finds the centroids for each particle in a helix.
	@param helix_coords: (x1, y1, x2, y2, width) -- the helix endpoints on the micrograph and width of the helix in pixels
	@param px_overlap: the overlap of consecutive particles in pixels, defaults to 90% of particle length
	@param px_length: the length in pixels of a particle along the axis of the helix
	@param px_width: the width in pixels of the particle
	@return: a list of coordinates on the micrograph for each particle's centroid: [(x0,y0),(x1,y1),(x2,y2),...]
	"""
	(x1,y1,x2,y2,w) = helix_coords
	l_vect = (x2-x1,y2-y1)
	helix_length = sqrt(l_vect[0]**2+l_vect[1]**2)
	l_uvect = (l_vect[0]/helix_length, l_vect[1]/helix_length)
	w_uvect = (-l_uvect[1],l_uvect[0])
	
	assert px_length > px_overlap, "The overlap must be smaller than the segment length"
	px_step = px_length - px_overlap
	l = px_length/2.0 #distance from (x1, y1) on the helix axis
	ptcl_coords = []
	while l < helix_length - px_length/2.0 + px_step:
		(x,y) = (x1 + l*l_uvect[0], y1 + l*l_uvect[1])
		if is_rotated:
			dst1 = ((x - x1)**2 + (y - y1)**2)**0.5
			dst2 = ((x - x2)**2 + (y - y2)**2)**0.5
			if 2*min(dst1, dst2) >= (px_length - 0.01):
				ptcl_coords.append((x,y))	
		else:
			ptcl_coords.append((x,y))	
		l += px_step
	
	return ptcl_coords

def get_rotated_particles( micrograph, helix_coords, px_dst = None, px_length = None, px_width = None, gridding = False , mic_name=""):
	"""
	Gets the overlapping square/rectangular "particles" with "lengths" (could be less than "widths") 
	parallel to the helical axis. They are then rotated so the "lengths" are vertical.
	@param micrograph: EMData object for the micrograph
	@param helix_coords: (x1,y1,x2,y2,width) tuple for a helix
	@param px_dst:    length of overlap in pixels of the rectangular particles, measured along the line connecting particle centroids, defaults to 0.1*px_length
	@param px_length: distance in pixels between the centroids of two adjacent particles, defaults to the width of the helix
	@param px_width:  width of the particles in pixels, defaults to px_length
	@param gridding:  use gridding method for rotation operation, requires square particles (px_length == px_width)
	@return: a list of EMData particles
	"""
	(x1,y1,x2,y2,w) = helix_coords
	l_vect = (x2-x1, y2-y1)
	length = sqrt(l_vect[0]**2+l_vect[1]**2)
	assert length != 0
	l_uvect = (l_vect[0]/length, l_vect[1]/length)
		
	#Rotate so that the length is parallel to the y-axis
	#Angle between l_uvect and y-axis: l_uvect (dot) j_hat = cos (rot_angle)
	rot_angle = 180/pi*acos( l_uvect[1] )
	#Whether we rotate clockwise or counterclockwise depends on the sign of l_uvect[0] (the x-component)
	if l_uvect[0] < 0: rot_angle *= -1 #rotate the box clockwise
	
	if not px_width:  px_width = w
	if not px_length: px_length = px_width
	if not px_dst:    px_overlap = 0.9*px_length
	else:             px_overlap = px_length - px_dst
	assert px_length > px_overlap, "The overlap must be smaller than the segment length"
	
	centroids = get_particle_centroids(helix_coords, px_overlap, px_length, px_width, is_rotated = True)
	particles = []
	if gridding == True:
		assert int(round(px_width)) == int(round(px_length)), "The segment must be square for gridding method"
		side1 = int(px_width + 0.5)
		sidepadded = int(side1*1.42 + 0.5)
		enth  = (sidepadded - side1)//2
		nxc   = micrograph.get_xsize()//2
		nyc   = micrograph.get_ysize()//2
		#print  px_width,px_length,px_dst,side1,sidepadded,enth,nxc,nyc
	else:
		ptcl_dimensions = ( int(round(px_width)), int(round(px_length)), 1 )


	for centroid in centroids:

		# liner interpolation for rotation
		if gridding == False:
			tr = Transform()
			tr.set_trans(centroid)
			tr.set_rotation({"type":"2d", "alpha":rot_angle})
			ptcl = micrograph.get_rotated_clip( tr, ptcl_dimensions )
			ptcl["ptcl_helix_coords"] = tuple(helix_coords)
			ptcl["ptcl_source_coord"] = tuple(centroid)
			ptcl["xform.align2d"]     = Transform()
			ptcl["xform.projection"]  = Transform()
			ptcl["ptcl_source_image"] = mic_name
			particles.append(ptcl)

		else:
			try:
				from EMAN2 import Util
				from fundamentals import rot_shift2D
				ptcl = Util.window( micrograph, sidepadded, sidepadded, 1, int(round(centroid[0] - nxc)), int(round(centroid[1] - nyc)), 0) 
				ptcl = Util.window( rot_shift2D(ptcl, -rot_angle, interpolation_method="gridding"), side1, side1, 1, 0, 0, 0)
				ptcl["ptcl_helix_coords"] = tuple(helix_coords)
				ptcl["ptcl_source_coord"] = tuple(centroid)
				ptcl["xform.align2d"]     = Transform()
				ptcl["xform.projection"]  = Transform()
				ptcl["ptcl_source_image"] = mic_name
				particles.append(ptcl)
			except:
				continue
	
	return particles

def get_unrotated_particles(micrograph, helix_coords, px_dst = None, px_length = None, px_width = None, mic_name=""):
	"""
	Gets the image data for each particle, without first rotating the helix and its corresponding particles to be vertical.
	@param micrograph: EMData object that holds the image data for the helix
	@param helix_coords: the (x1, y1, x2, y2, width) tuple (in pixels) that specifies the helix on the micrograph
	@param px_dst: the number of pixels between consecutive particles,measured along the long-axis of the helix, defaults to 10% of max(px_length, px_width)
	@param px_length: the distance between consecutive particle midpoints in pixels, chosen to correspond with rotated case, defaults to helix width
	@param px_width: corresponds to particle width in the rotated case, in the unrotated case only used to set length of the square to max(px_width, px_length), defaults to px_length
	@return: a list of EMData objects for each unrotated particle in the helix 
	"""
	#Will be square
	if px_length and px_width:
		side = max(px_length, px_width)
	elif not px_length and not px_width:
		side = helix_coords[4]
	elif px_width:
		side = px_width
	elif px_length:
		side = px_length
	
	if not px_dst:
		px_overlap = 0.9*side
	else:
		px_overlap = side - px_dst	
	centroids = get_particle_centroids(helix_coords, px_overlap, side, side)
	particles = []
	rot_angle = get_helix_rotation_angle(*helix_coords)
	#tr = Transform({"type":"eman","alt":90,"phi":rot_angle}) #How to rotate a cylinder that is oriented along the z-axis to align it along the boxed helix
	for centroid in centroids:
		ptcl= micrograph.get_clip( Region(centroid[0]-side/2.0, centroid[1]-side/2.0, side, side) )
		ptcl["ptcl_helix_coords"] = tuple(helix_coords)
		ptcl["ptcl_source_coord"] = tuple(centroid)
		ptcl["xform.projection"]  = Transform()
		ptcl["xform.align2d"]     = Transform()
		ptcl["ptcl_source_image"] = mic_name
		particles.append(ptcl)
	return particles

def load_helix_coords(coords_filepath, specified_width=None):
	"""
	load coordinates from a tab-delimited text file specifying helix coordinates as in EMAN1 *.box format    
	Uses the EMAN1 *.box file format (r is half the width (w) of the boxes)
		x1-r    y1-r    w    w    -1
		x2-r    y2-r    w    w    -2
	@param coords_filepath: file path to a tab-delimited text file specifying helix coordinates as in the EMAN1 *.box format
	@param specified_width: force all helix coordinates to have the specified width
	@return a list of tuples [(x0, x1, y1, y2, width), ...] 
	"""
	data = []
	datum = [None]*5
	for line in open(coords_filepath):
		line = line.split()
		for i in range(len(line)):
			line[i] = int(line[i])
		if line[4] == -1:
			w = line[2]
			r = w / 2.0
			datum[0] = line[0] + r
			datum[1] = line[1] + r
			if specified_width:
				datum[4] = specified_width
			else:
				datum[4] = w
		elif line[4] == -2:
			assert line[2] == w
			r = w / 2.0
			datum[2] = line[0] + r
			datum[3] = line[1] + r
			w = None
			r = None
			data.append(tuple(datum))
			datum = [None]*5
	
	return data
def save_helix_coords(coords_list, output_filepath, helix_width = None):
	"""
	Saves coordinates and widths of the boxed helices to a file.
	Uses the EMAN1 *.box file format (r is half the width (w) of the boxes)
		x1-r    y1-r    w    w    -1
		x2-r    y2-r    w    w    -2
	@param coords_list: a list of tuples (x1, y1, x2, y2, width), with each tuple corresponding to a helix
	@param output_filepath: the directory and file name in which to save the coordinates
	@param helix_width: if specified, it replaces the widths in coords_list as the width for each helix
	"""
	out_file = open(output_filepath, "w")
	
	for coords in coords_list:
		(x1, y1) = (coords[0], coords[1])
		(x2, y2) = (coords[2], coords[3])
		if helix_width:
			width = helix_width
		else:
			width = coords[4]
		r = width / 2.0
					
		#For some reason, EMAN1 subtracts half the box width from each coordinate
		#EMAN1 uses <cstdio> fprintf() and "%1.0f", which rounds half-integers away from zero
		#the string format operator works the same in Python as it does in C for decimal floats
		out_file.write( "%1.0f\t%1.0f\t%1.0f\t%1.0f\t-1\n" % (x1 - r, y1 - r, width, width) )
		out_file.write( "%1.0f\t%1.0f\t%1.0f\t%1.0f\t-2\n" % (x2 - r, y2 - r, width, width) )
	out_file.close()

def save_helix(helix_emdata, helix_filepath, helix_num):
	"""
	Saves a boxed helix to an image file.
	@param helix_emdata: the EMData object that holds the image data for the helix
	@param helix_filepath: a template for the output file path -- the value of helix_num will be added before the file extension
	@param helix_num: the number that identifies this helix among those boxed from this micrograph
	"""
	(path, ext) = os.path.splitext(helix_filepath)
	helix_filepath = "%s_%i%s" % (path, helix_num, ext)
	if os.access(helix_filepath, os.F_OK):
		os.remove(helix_filepath) #in case it's a stack file, otherwise helix_emdata.write_image(helix_filepath) only overwrites the first in the stack
	helix_emdata.write_image( helix_filepath )

def save_particle_coords(helix_particle_coords_dict, output_filepath, micrograph_filepath, ptcl_length, ptcl_width):
	"""
	Saves the coordinates on the micrograph for the center of each particle, with comment lines identifying the micrograph and helices for the particles
	@param helix_particle_coords_dict: {(h0_x1,h0_y1,h0_x2,h0_y2,h0_w):[(x0,y0),(x1,y1),(x2,y2),...],(h1_x1,h1_y1,h1_x2,h1_y2,h1_w):[(x0,y0),(x1,y1),(x2,y2),...],...}
	@param output_filepath: the directory and file name in which to save the coordinates
	@param micrograph_filepath: the directory and file name of the micrograph
	@param ptcl_length: the length of each particle
	@param ptcl_width: the width of each particle    
	"""
	out_file = open(output_filepath, "w")
	out_file.write("#micrograph: " + micrograph_filepath + "\n")
	out_file.write("#segment length: " + str(ptcl_length) + "\n")
	out_file.write("#segment width: " + str(ptcl_width) + "\n")
	for helix_coords in helix_particle_coords_dict.keys():
		out_file.write("#helix: " + str(tuple(helix_coords[0:2])) + "," + str(tuple(helix_coords[2:4])) + "," + str(helix_coords[4]) + "\n")
		particle_list = helix_particle_coords_dict[helix_coords]
		for ptcl_center in particle_list:
			out_file.write(str(ptcl_center[0]) + "\t" + str(ptcl_center[1]) + "\n")
	out_file.close()

def save_particles(particles, ptcl_filepath, do_edge_norm=False, stack_file_mode = "multiple"):
	"""
	saves the particles in a helix to a stack file
	@param particles: [ [helix0_particle0, helix0_particle1, ...], 
						[helix1_particle0, helix1_particle1, ...],
						...
					]
	@param ptcl_filepath: a template for the output file path, will be modified for modes that save multiple files
	@param do_edge_norm: Apply the processor "normalize.edgemean" to each particle before saving
	@param stack_file_mode: "single" -- saves all particles to a single stack file -- path/filename.extension
							"multiple" -- saves a stack file for each helix -- path/filename_0.extension for helix 0
							"none" -- saves a file for each particle -- path/filename_0_1.extension for helix 0 particle 1
							NOTE: if the file format does not support stack files, the "none" mode is used
	"""
	(path, ext) = os.path.splitext(ptcl_filepath)
	
	
	#Testing file writing support
	#If file type doesn't support writing, use HDF instead
	#If file type doesn't support stack files, store each particle individually regardless of stack_file_mode
	testdata = EMData(100,100)
	testdata.to_value(1)
	testfilename = ".HelixBoxerTestFile%s" % ext    
	try:
		testdata.write_image(testfilename, 0) #Test for write support
	except RuntimeError, e:
		ext = ".hdf"
	try:
		testdata.write_image(testfilename, 1) #Test for stack file support
	except RuntimeError, e:
		stack_file_mode = "none"
	finally:
		if os.access(testfilename, os.F_OK):
			os.remove(testfilename)


	
	if stack_file_mode == "single":
		ptcl_filepath = "%s%s" % (path, ext) #Note: ext may have been changed to a default type, so this is needed
		if os.access(ptcl_filepath, os.F_OK):
			os.remove(ptcl_filepath)
		for ptcl_lst in particles:
			for ptcl in ptcl_lst:
				if do_edge_norm:
					ptcl = ptcl.process("normalize.edgemean")
				ptcl.write_image(ptcl_filepath, -1) #appending to the image stack

	elif stack_file_mode == "multiple":
		for helix_num in range(len(particles)):
			ptcl_filepath = "%s_%i%s" % (path, helix_num, ext)
			if os.access(ptcl_filepath, os.F_OK):
				os.remove(ptcl_filepath)
			for ptcl_num in range(len(particles[helix_num])):
				ptcl = particles[helix_num][ptcl_num]
				if do_edge_norm:
					ptcl = ptcl.process("normalize.edgemean")
				ptcl.write_image(ptcl_filepath, ptcl_num) #appending to the image stack

	elif stack_file_mode == "none":
		for helix_num in range(len(particles)):
			ptcl_list = particles[helix_num]
			for ptcl_num in range(len(ptcl_list)):
				ptcl_filepath = "%s_%i_%i%s" % (path, helix_num, ptcl_num, ext)
				if os.access(ptcl_filepath, os.F_OK):
					os.remove(ptcl_filepath)
				ptcl = ptcl_list[ptcl_num]
				if do_edge_norm:
					ptcl = ptcl.process("normalize.edgemean")
				ptcl.write_image(ptcl_filepath)

def db_get_item(micrograph_filepath, key):
	"""
	gets the value stored in the sxhelixboxer database for the specified micrograph and key 
	"""
	db_name = SXHELIXBOXER_DB + key
	db = db_open_dict(db_name)
	val = db[micrograph_filepath]
	db_close_dict(db_name)
	return val

def db_set_item(micrograph_filepath, key, value):
	"""
	sets the value stored in the sxhelixboxer database for the specified micrograph and key 
	"""
	db_name = SXHELIXBOXER_DB + key
	db = db_open_dict(db_name)
	db[micrograph_filepath] = value
	db_close_dict(db_name)

def db_get_helices_dict(micrograph_filepath, helix_width = None):
	"""
	gets a dictionary of helices
	@param micrograph_filepath: the path to the image file for the micrograph
	@param helix_width: if specified, it replaces the widths in the database as the width for each helix
	@return: a dictionary formed like {(x1, y1, x2, y2, width): particle_EMData_object, ...}
	"""
	micrograph = EMData(micrograph_filepath)
	db = db_open_dict(SXHELIXBOXER_DB + "helixboxes")
	box_coords_list = db[micrograph_filepath]
	if not box_coords_list:
		return {}
	helices_dict = {}
	for coords in box_coords_list:
		if helix_width:
			coords[4] = helix_width
		helix = get_helix_from_coords(micrograph, *coords)
		helix["ptcl_source_image"] = micrograph_filepath
		helices_dict[tuple(coords)] = helix
	return helices_dict

def win_get_helices_dict(micrograph_filepath, helix_width = None):     #new function ##@ming
	"""
	gets a dictionary of helices
	@param micrograph_filepath: the path to the image file for the micrograph
	@param helix_width: if specified, it replaces the widths in the database as the width for each helix
	@return: a dictionary formed like {(x1, y1, x2, y2, width): particle_EMData_object, ...}
	"""
	micrograph = EMData(micrograph_filepath)
	db = db_open_dict(SXHELIXBOXER_DB + "helixboxes")
	box_coords_list = db[micrograph_filepath]
	if not box_coords_list:
		return {}
	helices_dict = {}
	#if helix_width:
	#	coords = [tuple( list(coords[:4]) + [helix_width] ) for coords in box_coords_list]
	for coords in box_coords_list:
		if helix_width:
			#coords = [tuple( list(coords[:4]) + [helix_width] )]
			coords=list(coords)
			coords[4]=helix_width
			coords=tuple(coords)
		#	coords[4] = helix_width
		helix = get_helix_from_coords(micrograph, *coords)
		helix["ptcl_source_image"] = micrograph_filepath
		helices_dict[tuple(coords)] = helix
	return helices_dict
	
	
def db_load_helix_coords(micrograph_filepath, coords_filepath, keep_current_boxes = True, specified_width=None):
	"""
	@param micrograph_filepath: the path to the image file for the micrograph
	@param coords_filepath: file path to a tab-delimited text file specifying helix coordinates as in the EMAN1 *.box format
	@param keep_current_boxes: whether to add to or replace the helix coordinates currently in the database
	@param specified_width: force all helix coordinates to have the specified width
	"""
	coords_list = load_helix_coords(coords_filepath, specified_width)
	if keep_current_boxes:
		db_coords = db_get_item(micrograph_filepath, "helixboxes")
		if specified_width:
			db_coords = [tuple( list(coords[:4]) + [specified_width] ) for coords in db_coords]
		for coords in coords_list:
			if not coords in db_coords:
				db_coords.append(coords)
		db_set_item(micrograph_filepath, "helixboxes", db_coords)
	else:
		db_set_item(micrograph_filepath, "helixboxes", coords_list)

def db_save_helix_coords(micrograph_filepath, output_filepath=None, helix_width = None):
	"""
	@param helix_width: if specified, it replaces the widths in coords_list as the width for each helix
	"""
	micrograph_filename = os.path.basename(micrograph_filepath)
	micrograph_name = os.path.splitext( micrograph_filename )[0]
	if not output_filepath:
		output_filepath = os.getcwd() + micrograph_name + "_boxes.txt"
	db = db_open_dict(SXHELIXBOXER_DB + "helixboxes")
	box_coords_list = db[micrograph_filepath]
	save_helix_coords(box_coords_list, output_filepath, helix_width)

def db_save_helices(micrograph_filepath, helix_filepath=None, helix_width = None):
	"""
	@param helix_width: if specified, it replaces the widths in coords_list as the width for each helix
	"""
	micrograph_filename = os.path.basename(micrograph_filepath)
	micrograph_name = os.path.splitext( micrograph_filename )[0]
	if not helix_filepath:
		helix_filepath = "%s_helix.hdf" % ( os.path.join(os.getcwd(), micrograph_name) )
	helices_dict = db_get_helices_dict(micrograph_filepath, helix_width)
	i = 0
	for coords in helices_dict:
		helix = helices_dict[coords]
		save_helix(helix, helix_filepath, i)
		i+=1

def db_save_particle_coords(micrograph_filepath, output_filepath = None, px_dst = None, px_length = None, px_width = None, do_rotation = False):
	base = os.path.splitext( micrograph_filepath )[0]
	if not output_filepath:
		output_filepath = "%s_helix_ptcl_coords.txt" % (base)
	helix_coords_list = db_get_item(micrograph_filepath, "helixboxes")
	helix_width = helix_coords_list[0][4]
	if not px_width:
		px_width = helix_width
	if not px_length:
		px_length = px_width
	if px_dst:
		px_overlap = px_length - px_dst
	else:
		px_overlap = 0.9*px_length
	assert px_overlap < px_length
	helix_particle_coords_dict = {}
	for helix_coords in helix_coords_list:
		particle_centroids = get_particle_centroids(helix_coords, px_overlap, px_length, px_width, is_rotated = do_rotation)
		helix_particle_coords_dict[tuple(helix_coords)] = particle_centroids
	
	save_particle_coords(helix_particle_coords_dict, output_filepath, micrograph_filepath, px_length, px_width)
	
def db_save_particles(micrograph_filepath, ptcl_filepath = None, px_dst = None, px_length = None, px_width = None, rotated = True, do_edge_norm = False, gridding = False, stack_file_mode = "multiple", do_filt = True, filt_freq = -1):
	from filter 		import filt_gaussh
	
	micrograph_filename = os.path.basename(micrograph_filepath)
	micrograph_name = os.path.splitext( micrograph_filename )[0]
	if not ptcl_filepath:
		ptcl_filepath = "%s_helix_ptcl.hdf" % ( os.path.join(os.getcwd(), micrograph_name) )
	(base, ext) = os.path.splitext(ptcl_filepath)
	helices_dict = db_get_helices_dict(micrograph_filepath)
	
	all_particles = []
	micrograph = EMData(micrograph_filepath)
	if do_filt:
		if filt_freq <= 0:
			filt_freq = 1.0/px_length
		micrograph = filt_gaussh(micrograph, filt_freq)		
	nhelix = 0
	for coords in helices_dict:
		helix = helices_dict[coords]
		if rotated:
			helix_particles = get_rotated_particles(micrograph, coords, px_dst, px_length, px_width, gridding, mic_name = micrograph_filename)
		else:
			helix_particles = get_unrotated_particles(micrograph, coords, px_dst, px_length, px_width,mic_name = micrograph_filename)
		for ii in xrange(len(helix_particles)):
			(helix_particles[ii]).set_attr("filament", micrograph_filename+"%04d"%nhelix)
		nhelix = nhelix + 1
		all_particles.append(helix_particles)

	save_particles(all_particles, ptcl_filepath, do_edge_norm, stack_file_mode)

if ENABLE_GUI:
	class EMWriteHelixFilesDialog(QtGui.QDialog):
		"""
		options for writing helices and particles to files
		"""
		def __init__(self, qparent=None,saveext="hdf"):
			QtGui.QWidget.__init__(self, qparent)
			self.setWindowTitle(self.tr("Write Helix and Particle Files"))
			self.__create_ui()
	#        self.helices_file_extension_dict = {"MRC":"mrc", "Spider":"spi", "Imagic": "img", "HDF5": "hdf"}
	#        self.ptcls_file_extension_dict = {"Spider":"spi", "Imagic": "img", "HDF5": "hdf"}

			width = self.parentWidget().get_width()
			self.ptcls_width_spinbox.setValue( width )
			self.ptcls_length_spinbox.setValue( width )
			self.ptcls_distance_spinbox.setValue( int(0.1*width) )


			micrograph_filepath = self.parentWidget().micrograph_filepath
			(micrograph_dir, micrograph_filename) = os.path.split(micrograph_filepath)
			self.micrograph_filename = micrograph_filename
			self.micrograph_name = os.path.splitext(micrograph_filename)[0]
			self.default_dir = os.getcwd()

			self.helices_coords_line_edit.setText( os.path.join(self.default_dir, self.micrograph_name + "_boxes.txt") )
			self.helices_images_line_edit.setText( os.path.join(self.default_dir, self.micrograph_name + "_helix."+saveext) )
			self.ptcls_coords_line_edit.setText( os.path.join(self.default_dir, self.micrograph_name + "_helix_ptcl_coords.txt") )
			self.ptcls_images_line_edit.setText( os.path.join(self.default_dir, self.micrograph_name + "_helix_ptcl."+saveext) )

			self.connect(self.helices_coords_browse_button, QtCore.SIGNAL("clicked()"), self.browse_helix_coords)
			self.connect(self.helices_images_browse_button, QtCore.SIGNAL("clicked()"), self.browse_helix_images)
			self.connect(self.ptcls_coords_browse_button, QtCore.SIGNAL("clicked()"), self.browse_ptcl_coords)
			self.connect(self.ptcls_images_browse_button, QtCore.SIGNAL("clicked()"), self.browse_ptcl_images)
			self.connect(self.button_box, QtCore.SIGNAL("accepted()"), self.save)
			self.connect(self.button_box, QtCore.SIGNAL("rejected()"), self.cancel)

		def __create_ui(self):
			self.helices_groupbox = QtGui.QGroupBox(self.tr("Write &Helices:"))
			self.helices_groupbox.setCheckable(True)
			
			self.helices_coords_groupbox = QtGui.QGroupBox(self.tr("Helix Coordinates (EMAN1 format)"))
			self.helices_coords_groupbox.setCheckable(True)
			helices_coords_label = QtGui.QLabel(self.tr("Path:"))
			self.helices_coords_line_edit = QtGui.QLineEdit()
			self.helices_coords_line_edit.setMinimumWidth(300)
			self.helices_coords_browse_button = QtGui.QPushButton(self.tr("Browse"))
			
			self.helices_images_groupbox = QtGui.QGroupBox(self.tr("Helix Images"))
			self.helices_images_groupbox.setCheckable(True)
			helices_images_label = QtGui.QLabel(self.tr("Path:"))
			self.helices_images_line_edit = QtGui.QLineEdit()
			self.helices_images_browse_button = QtGui.QPushButton(self.tr("Browse"))

			self.ptcls_groupbox = QtGui.QGroupBox(self.tr("Write &Segments:"))
			self.ptcls_groupbox.setCheckable(True)

			ptcls_distance_label = QtGui.QLabel(self.tr("&Distance:"))
			self.ptcls_distance_spinbox = QtGui.QSpinBox()
			self.ptcls_distance_spinbox.setMaximum(10000)
			ptcls_distance_label.setBuddy(self.ptcls_distance_spinbox)
			ptcls_length_label = QtGui.QLabel(self.tr("&Length:"))
			self.ptcls_length_spinbox = QtGui.QSpinBox()
			self.ptcls_length_spinbox.setMaximum(10000)
			ptcls_length_label.setBuddy(self.ptcls_length_spinbox)
			ptcls_width_label = QtGui.QLabel(self.tr("W&idth:"))
			self.ptcls_width_spinbox = QtGui.QSpinBox()
			self.ptcls_width_spinbox.setMaximum(10000)
			ptcls_width_label.setBuddy(self.ptcls_width_spinbox)

			self.ptcls_coords_groupbox = QtGui.QGroupBox(self.tr("Segment Coordinates"))
			self.ptcls_coords_groupbox.setCheckable(True)
			ptcls_coords_label = QtGui.QLabel(self.tr("Path:"))
			self.ptcls_coords_line_edit = QtGui.QLineEdit()
			self.ptcls_coords_browse_button = QtGui.QPushButton(self.tr("Browse"))
				
			self.ptcls_images_groupbox = QtGui.QGroupBox(self.tr("Segment Images"))
			self.ptcls_images_groupbox.setCheckable(True)
			self.ptcls_edgenorm_checkbox = QtGui.QCheckBox(self.tr("&Normalize Edge-Mean"))
			self.ptcls_edgenorm_checkbox.setChecked(False)
			self.ptcls_edgenorm_checkbox.setToolTip("Uses normalize.edgemean processor on each segment: pixel-value -> (pixel-value - edge-mean) / standard deviation")
			
			self.ptcls_rotation_groupbox = QtGui.QGroupBox(self.tr("Rotation"))
			self.ptcls_bilinear_rotation_radiobutton = QtGui.QRadioButton(self.tr("Bilinear Rotation"))
			self.ptcls_bilinear_rotation_radiobutton.setToolTip("Rectangular segments. Rotation angle is the one that makes associated helix vertical. Bilinear rotation algorithm.")
			self.ptcls_bilinear_rotation_radiobutton.setChecked(True)
			self.ptcls_gridding_rotation_radiobutton = QtGui.QRadioButton(self.tr("Gridding Rotation"))
			self.ptcls_gridding_rotation_radiobutton.setToolTip("Square segments with sides = max(Length, Width). Rotation angle is the one that makes associated helix vertical. Gridding rotation algorithm.")
			self.ptcls_no_rotation_radiobutton = QtGui.QRadioButton(self.tr("No Rotation"))
			self.ptcls_no_rotation_radiobutton.setToolTip("Segments are not rotated from the micrograph. Square segments with sides = max(length, width)")
			
			self.ptcls_stack_groupbox = QtGui.QGroupBox(self.tr("Image Stacks"))
			self.ptcls_single_stack_radiobutton = QtGui.QRadioButton(self.tr("Single image stack"))
			self.ptcls_single_stack_radiobutton.setChecked(True)
			self.ptcls_single_stack_radiobutton.setToolTip("Saves a single image stack file for all the helices. Fails for incompatible file formats.")
			self.ptcls_multiple_stack_radiobutton = QtGui.QRadioButton(self.tr("Image stack per helix"))
			self.ptcls_multiple_stack_radiobutton.setToolTip("Saves an image stack file for each helix. Fails for incompatible file formats.")
			self.ptcls_no_stack_radiobutton = QtGui.QRadioButton(self.tr("File for each segment"))

			ptcls_images_label = QtGui.QLabel(self.tr("Path:"))
			self.ptcls_images_line_edit = QtGui.QLineEdit()
			self.ptcls_images_browse_button = QtGui.QPushButton(self.tr("Browse"))

			self.button_box = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Save | QtGui.QDialogButtonBox.Cancel)


			
			helices_coords_layout = QtGui.QHBoxLayout()
			helices_coords_layout.addWidget(helices_coords_label)
			helices_coords_layout.addWidget(self.helices_coords_line_edit)
			helices_coords_layout.addWidget(self.helices_coords_browse_button)
			self.helices_coords_groupbox.setLayout(helices_coords_layout)
			
			helices_images_layout = QtGui.QHBoxLayout()
			helices_images_layout.addWidget(helices_images_label)
			helices_images_layout.addWidget(self.helices_images_line_edit)
			helices_images_layout.addWidget(self.helices_images_browse_button)
			self.helices_images_groupbox.setLayout(helices_images_layout)
			
			helices_layout = QtGui.QVBoxLayout()
			helices_layout.addWidget(self.helices_coords_groupbox)
			helices_layout.addWidget(self.helices_images_groupbox)
			
			self.helices_groupbox.setLayout(helices_layout)
			
			ptcls_distance_layout = QtGui.QHBoxLayout()
			ptcls_distance_layout.addWidget(ptcls_distance_label)
			ptcls_distance_layout.addWidget(self.ptcls_distance_spinbox)
			
			ptcls_length_layout = QtGui.QHBoxLayout()
			ptcls_length_layout.addWidget(ptcls_length_label)
			ptcls_length_layout.addWidget(self.ptcls_length_spinbox)
			
			ptcls_width_layout = QtGui.QHBoxLayout()
			ptcls_width_layout.addWidget(ptcls_width_label)
			ptcls_width_layout.addWidget(self.ptcls_width_spinbox)
			
			ptcls_coords_layout = QtGui.QHBoxLayout()
			ptcls_coords_layout.addWidget(ptcls_coords_label)
			ptcls_coords_layout.addWidget(self.ptcls_coords_line_edit)
			ptcls_coords_layout.addWidget(self.ptcls_coords_browse_button)
			self.ptcls_coords_groupbox.setLayout(ptcls_coords_layout)
			
			ptcls_images_path_layout = QtGui.QHBoxLayout()
			ptcls_images_path_layout.addWidget(ptcls_images_label)
			ptcls_images_path_layout.addWidget(self.ptcls_images_line_edit)
			ptcls_images_path_layout.addWidget(self.ptcls_images_browse_button)
			
			ptcls_images_rotation_layout = QtGui.QVBoxLayout()
			ptcls_images_rotation_layout.addWidget(self.ptcls_bilinear_rotation_radiobutton)
			ptcls_images_rotation_layout.addWidget(self.ptcls_gridding_rotation_radiobutton)
			ptcls_images_rotation_layout.addWidget(self.ptcls_no_rotation_radiobutton)
			self.ptcls_rotation_groupbox.setLayout(ptcls_images_rotation_layout)
			
			ptcls_imagestack_layout = QtGui.QVBoxLayout()
			ptcls_imagestack_layout.addWidget(self.ptcls_single_stack_radiobutton)
			ptcls_imagestack_layout.addWidget(self.ptcls_multiple_stack_radiobutton)
			ptcls_imagestack_layout.addWidget(self.ptcls_no_stack_radiobutton)
			self.ptcls_stack_groupbox.setLayout(ptcls_imagestack_layout)
			
			ptcls_rotation_stack_layout = QtGui.QHBoxLayout()
			ptcls_rotation_stack_layout.addWidget(self.ptcls_rotation_groupbox)
			ptcls_rotation_stack_layout.addWidget(self.ptcls_stack_groupbox)
			
			ptcls_images_layout = QtGui.QVBoxLayout()
			ptcls_images_layout.addLayout(ptcls_rotation_stack_layout)
			ptcls_images_layout.addWidget(self.ptcls_edgenorm_checkbox)
			ptcls_images_layout.addLayout(ptcls_images_path_layout)
			self.ptcls_images_groupbox.setLayout(ptcls_images_layout)
			
			ptcls_opts_layout = QtGui.QVBoxLayout()
			ptcls_opts_layout.addLayout(ptcls_distance_layout)
			ptcls_opts_layout.addLayout(ptcls_length_layout)
			ptcls_opts_layout.addLayout(ptcls_width_layout)
			ptcls_opts_layout.addWidget(self.ptcls_coords_groupbox)
			ptcls_opts_layout.addWidget(self.ptcls_images_groupbox)
			self.ptcls_groupbox.setLayout(ptcls_opts_layout)
	
			self.vbl = QtGui.QVBoxLayout(self)
			self.vbl.setMargin(0)
			self.vbl.setSpacing(6)
			self.vbl.setObjectName("vbl")
			self.vbl.addWidget(self.helices_groupbox)
			self.vbl.addWidget(self.ptcls_groupbox)
			self.vbl.addWidget(self.button_box)
		def browse_helix_coords(self):
			file_dlg = QtGui.QFileDialog(self,self.tr("Save Helix Coordinates"))
			file_dlg.setAcceptMode(QtGui.QFileDialog.AcceptSave)
			file_dlg.selectFile( os.path.join(self.default_dir, self.micrograph_name + "_boxes.txt") )
			if file_dlg.exec_():
				file_path = file_dlg.selectedFiles()[0]
				file_path = str(file_path)
				self.helices_coords_line_edit.setText(file_path)
		def browse_helix_images(self):
			file_dlg = QtGui.QFileDialog(self,self.tr("Save Helix Images"))
			file_dlg.setAcceptMode(QtGui.QFileDialog.AcceptSave)
			file_dlg.selectFile(self.helices_images_line_edit.text())
			if file_dlg.exec_():
				file_path = file_dlg.selectedFiles()[0]
				file_path = str(file_path)
				self.helices_images_line_edit.setText(file_path)
	
	# EMSelector isn't working well: it sets the filename to the innermost directory and the file filter doesn't work        
	#        selector = EMSelectorDialog(single_selection=True,save_as_mode=False)
	#        path = selector.exec_()
	#        if path:
	#            path = os.path.dirname(path)
	#            if os.path.isdir(path):
	#                path = os.path.join( path, os.path.basename(self.helices_images_line_edit.text()) )
	#            self.helices_images_line_edit.setText(path)
		def browse_ptcl_coords(self):
			file_dlg = QtGui.QFileDialog(self,self.tr("Save Helix Coordinates"))
			file_dlg.setAcceptMode(QtGui.QFileDialog.AcceptSave)
			file_dlg.selectFile(self.ptcls_coords_line_edit.text())
			if file_dlg.exec_():
				file_path = file_dlg.selectedFiles()[0]
				file_path = str(file_path)
				self.ptcls_coords_line_edit.setText(file_path)
		def browse_ptcl_images(self):
			file_dlg = QtGui.QFileDialog(self,self.tr("Save Helix Images"))
			file_dlg.setAcceptMode(QtGui.QFileDialog.AcceptSave)
			file_dlg.selectFile(self.ptcls_images_line_edit.text())
			if file_dlg.exec_():
				file_path = file_dlg.selectedFiles()[0]
				file_path = str(file_path)
				self.ptcls_images_line_edit.setText(file_path)
	# EMSelector isn't working well: it sets the filename to the innermost directory and the file filter doesn't work
	#        selector = EMSelectorDialog(single_selection=True,save_as_mode=False)
	#        path = selector.exec_()
	#        if path:
	#            if os.path.isdir(path):
	#                path = os.path.join( path, os.path.basename(self.ptcls_images_line_edit.text()) )
	#            self.ptcls_images_line_edit.setText(path)
		def cancel(self):
			self.hide()
		def save(self):
			"""
			writes the image data for the helices and segments to files if each of those options are checked
			"""
			helices_dict = self.parentWidget().helices_dict
			micrograph = self.parentWidget().main_image.get_data()
			
			if self.helices_groupbox.isChecked():
				if self.helices_coords_groupbox.isChecked():
					path = str( self.helices_coords_line_edit.text() )
					save_helix_coords(helices_dict.keys(), path)
				if self.helices_images_groupbox.isChecked():
					helix_filepath = str(self.helices_images_line_edit.text())
					i = 0
					for coords_key in helices_dict:
						helix = helices_dict[coords_key]
						save_helix(helix, helix_filepath, i)
						i += 1
			if self.ptcls_groupbox.isChecked():
				px_dst = self.ptcls_distance_spinbox.value()
				px_length = self.ptcls_length_spinbox.value()
				px_width = self.ptcls_width_spinbox.value()
				px_overlap = None
				if self.ptcls_images_groupbox.isChecked():
					do_edge_norm = self.ptcls_edgenorm_checkbox.isChecked()
					ptcl_filepath = str(self.ptcls_images_line_edit.text())

					if self.ptcls_single_stack_radiobutton.isChecked():
						stack_mode = "single"
					elif self.ptcls_multiple_stack_radiobutton.isChecked():
						stack_mode = "multiple"
					elif self.ptcls_no_stack_radiobutton.isChecked():
						stack_mode = "none"

					all_particles = []
					nhelix = 0
					for coords_key in helices_dict:
						if self.ptcls_bilinear_rotation_radiobutton.isChecked():
							helix_particles = get_rotated_particles(micrograph, coords_key, px_dst, px_length, px_width, gridding=False, mic_name=self.micrograph_filename)
							px_overlap = px_length - px_dst
						elif self.ptcls_gridding_rotation_radiobutton.isChecked():
							side = max(px_length, px_width)
							#need to prepare the fft volume for gridding at here
							helix_particles = get_rotated_particles(micrograph, coords_key, px_dst, side, side, gridding=True , mic_name=self.micrograph_filename)
							px_overlap = side - px_dst
						elif self.ptcls_no_rotation_radiobutton.isChecked():
							side = max(px_length, px_width)
							helix_particles = get_unrotated_particles(micrograph, coords_key, px_dst, side, side, mic_name=self.micrograph_filename)
							px_overlap = side - px_dst
						for ii in xrange(len(helix_particles)):
							(helix_particles[ii]).set_attr("filament", self.micrograph_filename+"%04d"%nhelix)
						nhelix = nhelix + 1
						all_particles.append(helix_particles)
					save_particles(all_particles, ptcl_filepath, do_edge_norm, stack_mode)
				if self.ptcls_coords_groupbox.isChecked():
					if not px_overlap:
						px_overlap = px_length - px_dst
					ptcl_coords_filepath = str(self.ptcls_coords_line_edit.text())
					micrograph_filepath = self.parentWidget().micrograph_filepath
					helix_particle_coords_dict = {}
					for coords_key in helices_dict:
						ptcl_coords_list = get_particle_centroids(coords_key, px_overlap, px_length, px_width, is_rotated=not(self.ptcls_no_rotation_radiobutton.isChecked()))
						helix_particle_coords_dict[coords_key] = ptcl_coords_list
					save_particle_coords(helix_particle_coords_dict, ptcl_coords_filepath, micrograph_filepath, px_length, px_width)
			
			self.hide()
			
if ENABLE_GUI:
	class EMHelixBoxerWidget(QtGui.QWidget):
		"""
		the GUI widget which contains the settings for boxing helices and writing results to files
		"""
		def __init__(self, micrograph_filepaths, app, box_width=100, saveext="hdf", invert_contrast=False):
			"""
			@param micrograph_filepath: the path to the image file for the micrograph
			@param app: the application to which this widget belongs
			"""
			QtGui.QWidget.__init__(self)

			self.doctf = False
			self.winsize = 512
			self.cs = 2.0
			self.edge = 0.0
			self.volt = 300.0
			self.kboot = 16
			self.ov = 0
			self.ac = 10.0
			self.pixelsize = 1.0

			self.invert_contrast = invert_contrast

			if box_width<1 : box_width=100
			self.box_width=box_width

			self.saveext=saveext
			self.app = app
			self.setWindowIcon(QtGui.QIcon(get_image_directory() +"green_boxes.png"))
			self.setWindowTitle("sxhelixboxer")

			self.main_image = None #Will be an EMImage2DWidget instance
			self.helix_viewer = None #Will be an EMImage2DWidget instance
			self.write_helix_files_dlg = None
			
			self.color = (0, 0, 1)
			self.selected_color = (0, 1, 0)
			self.counter = counterGen()
			self.micrograph_filepath = None
	
			self.__create_ui()

			if sys.version_info >= (2, 6):
				self.micrograph_filepath_set = set([ os.path.relpath(path) for path in micrograph_filepaths ]) #os.path.relpath is new in Python 2.6
			else:
				self.micrograph_filepath_set = set(micrograph_filepaths) # [micrograph1_filepath, micrograph2_filepath, ...]
			self.update_micrograph_table()

			self.connect(self.box_width_spinbox, QtCore.SIGNAL("valueChanged(int)"), self.width_changed)
			self.connect( self.img_quality_combobox, QtCore.SIGNAL("currentIndexChanged(int)"), self.set_image_quality )
			self.connect(self.load_boxes_action, QtCore.SIGNAL("triggered()"), self.load_boxes)
			self.connect(self.load_micrograph_action, QtCore.SIGNAL("triggered()"), self.open_micrograph)
	#        self.connect(self.write_coords_action, QtCore.SIGNAL("triggered()"), self.write_coords)
			self.connect(self.write_images_action, QtCore.SIGNAL("triggered()"), self.write_images)
			self.connect(self.quit_action, QtCore.SIGNAL("triggered()"), self.close)
			self.connect( self.micrograph_table, QtCore.SIGNAL("currentCellChanged (int,int,int,int)"), self.micrograph_table_selection)
			
			self.micrograph_table.setCurrentCell(0,0) #self.micrograph_table_selection() will display this micrograph
			
			if box_width > 0:
				self.box_width_spinbox.setValue(box_width)
				self.width_changed(box_width)
			
		def __create_ui(self):
			
			self.menu_bar = QtGui.QMenuBar(self)
			self.file_menu = QtGui.QMenu(self.tr("&File"))
			self.load_micrograph_action = QtGui.QAction(self.tr("&Open Micrographs"), self)
	#        self.write_coords_action = QtGui.QAction(self.tr("Save &Coordinates"), self)
			self.write_images_action = QtGui.QAction(self.tr("&Save"), self)
			self.load_boxes_action = QtGui.QAction(self.tr("&Load Coordinates"), self)
			self.quit_action = QtGui.QAction(self.tr("&Quit"), self)
			self.file_menu.addAction(self.load_micrograph_action)
			self.file_menu.addAction(self.load_boxes_action)
	#        self.file_menu.addAction(self.write_coords_action)
			self.file_menu.addAction(self.write_images_action)
			self.file_menu.addSeparator()
			self.file_menu.addAction(self.quit_action)
			self.menu_bar.addMenu(self.file_menu)
			
			self.box_width_label = QtGui.QLabel(self.tr("Box &Width:"))
			self.box_width_spinbox = QtGui.QSpinBox()
			self.box_width_spinbox.setMaximum(10000)
			self.box_width_label.setBuddy(self.box_width_spinbox)
			
			self.img_quality_label = QtGui.QLabel(self.tr("Image &Quality:"))
			self.img_quality_combobox = QtGui.QComboBox()
			qualities = [str(i) for i in range(5)]
			self.img_quality_combobox.addItems(qualities)
			self.img_quality_combobox.setCurrentIndex(2)
			self.img_quality_label.setBuddy(self.img_quality_combobox)
			
			self.micrograph_table = QtGui.QTableWidget(1,2)
			self.micrograph_table.setHorizontalHeaderLabels(["Micrograph", "Boxed Helices"])
			
			self.status_bar = QtGui.QStatusBar()
			#self.status_bar.showMessage("Ready",10000)
			
			widthLayout = QtGui.QHBoxLayout()
			widthLayout.addWidget(self.box_width_label)
			widthLayout.addWidget(self.box_width_spinbox)
			
			qualityLayout = QtGui.QHBoxLayout()
			qualityLayout.addWidget(self.img_quality_label)
			qualityLayout.addWidget(self.img_quality_combobox)
				
			self.vbl = QtGui.QVBoxLayout(self)
			self.vbl.setMargin(30)
			self.vbl.setSpacing(6)
			self.vbl.setObjectName("vbl")
			self.vbl.addWidget(self.menu_bar)
			self.vbl.addLayout(widthLayout)
			self.vbl.addLayout(qualityLayout)
			self.vbl.addWidget(self.micrograph_table)
			self.vbl.addWidget(self.status_bar)
			
			# add input fields for CTF estimation
			hbl_doctf = QtGui.QHBoxLayout()
			self.doctf_chk = QtGui.QCheckBox("To estimate CTF use sxcter.py instead!!!                               ")
			self.doctf_chk.setToolTip("CTF Estimation using CTER")
			self.doctf_chk.setChecked(self.doctf)
			hbl_doctf.addWidget(self.doctf_chk)
			self.vbl.addLayout(hbl_doctf)
			
			QtCore.QObject.connect(self.doctf_chk,QtCore.SIGNAL("clicked(bool)"),self.doctf_checked)

		def doctf_checked(self,val):
			if not(self.doctf):
				self.doctf = val
			
				if val:
					hgctf = QtGui.QHBoxLayout()
					ctftitle = QtGui.QLabel("<b>Parameters of CTF estimation</b>")
					hgctf.addWidget(ctftitle)
					self.vbl.addLayout(hgctf)
					
					hbl_wscs = QtGui.QHBoxLayout()
					window_size_label = QtGui.QLabel("Window size:")
					hbl_wscs.addWidget(window_size_label)
					self.ctf_window_size = QtGui.QLineEdit('512')
					hbl_wscs.addWidget(self.ctf_window_size)
					
					cs_label = QtGui.QLabel("Cs:")
					hbl_wscs.addWidget(cs_label)
					self.ctf_cs = QtGui.QLineEdit('2.0')
					hbl_wscs.addWidget(self.ctf_cs)
					self.vbl.addLayout(hbl_wscs)
					
					
					hbl_esv = QtGui.QHBoxLayout()
					edge_size_label = QtGui.QLabel("Edge size:")
					hbl_esv.addWidget(edge_size_label)
					self.ctf_edge_size = QtGui.QLineEdit('0')
					hbl_esv.addWidget(self.ctf_edge_size)
					
					voltage_label = QtGui.QLabel("Voltage:")
					hbl_esv.addWidget(voltage_label)
					self.ctf_volt = QtGui.QLineEdit('200.0')
					hbl_esv.addWidget(self.ctf_volt)
					self.vbl.addLayout(hbl_esv)
					
					hbl_oac = QtGui.QHBoxLayout()
					overlap_label = QtGui.QLabel("Overlap:")
					hbl_oac.addWidget(overlap_label)
					self.ctf_overlap_size = QtGui.QLineEdit('50')
					hbl_oac.addWidget(self.ctf_overlap_size)
					
					amplitude_contrast_label = QtGui.QLabel("Amplitude Contrast:")
					hbl_oac.addWidget(amplitude_contrast_label)
					self.ctf_ampcont = QtGui.QLineEdit('10.0')
					hbl_oac.addWidget(self.ctf_ampcont)
					self.vbl.addLayout(hbl_oac)
					
					hbl_kboot = QtGui.QHBoxLayout()
					kboot_label = QtGui.QLabel("kboot (only for CTER):")
					hbl_kboot.addWidget(kboot_label)
					self.ctf_kboot = QtGui.QLineEdit('16')
					hbl_kboot.addWidget(self.ctf_kboot)
					
					pixel_label = QtGui.QLabel("Pixel size:")
					hbl_kboot.addWidget(pixel_label)
					self.ctf_pixel = QtGui.QLineEdit('1.0')
					hbl_kboot.addWidget(self.ctf_pixel)
					
					self.vbl.addLayout(hbl_kboot)
					
					hbl_estdef = QtGui.QHBoxLayout()
					estimated_defocus_label = QtGui.QLabel("Estimated defocus:")
					hbl_estdef.addWidget(estimated_defocus_label)
					self.estdef = QtGui.QLineEdit('')
					hbl_estdef.addWidget(self.estdef)
					self.vbl.addLayout(hbl_estdef)
					
					hbl_astamp = QtGui.QHBoxLayout()
					astig_amp_label = QtGui.QLabel("Estimated astigmatism amplitude\n (only for CTER):")
					hbl_astamp.addWidget(astig_amp_label)
					self.astamp = QtGui.QLineEdit('')
					hbl_astamp.addWidget(self.astamp)
					self.vbl.addLayout(hbl_astamp)
					
					hbl_astagl = QtGui.QHBoxLayout()
					astig_angle_label = QtGui.QLabel("Estimated astigmatism angle \n(only for CTER)")
					hbl_astagl.addWidget(astig_angle_label)
					self.astagl = QtGui.QLineEdit('')
					hbl_astagl.addWidget(self.astagl)
					self.vbl.addLayout(hbl_astagl)
					
					hbl_deferr = QtGui.QHBoxLayout()
					deferr_label = QtGui.QLabel("Estimated defocus error \n(only for CTER):")
					hbl_deferr.addWidget(deferr_label)
					self.deferr = QtGui.QLineEdit('')
					hbl_deferr.addWidget(self.deferr)
					self.vbl.addLayout(hbl_deferr)
					
					hbl_astaglerr = QtGui.QHBoxLayout()
					astaglerr_label = QtGui.QLabel("Estimated astigmatism angle error \n(only for CTER):")
					hbl_astaglerr.addWidget(astaglerr_label)
					self.astaglerr = QtGui.QLineEdit('')
					hbl_astaglerr.addWidget(self.astaglerr)
					self.vbl.addLayout(hbl_astaglerr)
					
					hbl_astamperr = QtGui.QHBoxLayout()
					astamperr_label = QtGui.QLabel("Estimated astigmatism amplitude error \n (only for CTER):")
					hbl_astamperr.addWidget(astamperr_label)
					self.astamperr = QtGui.QLineEdit('')
					hbl_astamperr.addWidget(self.astamperr)
					self.vbl.addLayout(hbl_astamperr)
					
					hbl_ctf_cter = QtGui.QHBoxLayout()
					self.estimate_ctf_cter =QtGui.QPushButton("Estimate CTF using CTER")
					hbl_ctf_cter.addWidget(self.estimate_ctf_cter)
					self.vbl.addLayout(hbl_ctf_cter)
					
					QtCore.QObject.connect(self.ctf_window_size,QtCore.SIGNAL("editingFinished()"),self.new_ctf_window)
					QtCore.QObject.connect(self.ctf_cs,QtCore.SIGNAL("editingFinished()"),self.new_ctf_cs)
					QtCore.QObject.connect(self.ctf_edge_size,QtCore.SIGNAL("editingFinished()"),self.new_ctf_edge)
					QtCore.QObject.connect(self.ctf_volt,QtCore.SIGNAL("editingFinished()"),self.new_ctf_volt)
					QtCore.QObject.connect(self.ctf_overlap_size,QtCore.SIGNAL("editingFinished()"),self.new_ctf_overlap_size)
					QtCore.QObject.connect(self.ctf_ampcont,QtCore.SIGNAL("editingFinished()"),self.new_ctf_ampcont)
					QtCore.QObject.connect(self.ctf_kboot,QtCore.SIGNAL("editingFinished()"),self.new_ctf_kboot)
					QtCore.QObject.connect(self.ctf_pixel,QtCore.SIGNAL("editingFinished()"),self.new_ctf_pixel)
					QtCore.QObject.connect(self.estimate_ctf_cter,QtCore.SIGNAL("clicked(bool)"), self.calc_ctf_cter)
	
		def new_ctf_pixel(self):
			self.pixelsize=self.ctf_pixel.text()
			
		def new_ctf_window(self):
			self.winsize=self.ctf_window_size.text()

		def new_ctf_cs(self):
			self.cs=self.ctf_cs.text()
			
		def new_ctf_edge(self):
			self.edge=self.ctf_edge_size.text()
		
		def new_ctf_volt(self):
			self.volt=self.ctf_volt.text()
			
		def new_ctf_kboot(self):
			self.kboot=self.ctf_kboot.text()
			
		def new_ctf_overlap_size(self):
			self.ov=self.ctf_overlap_size.text()
		
		def new_ctf_ampcont(self):
			self.ac=self.ctf_ampcont.text()
		
		def calc_ctf_cter(self):
			# calculate ctf of ORIGINAL micrograph using cter in gui mode
			# this must mean cter is being calculated on a single micrograph!
			from utilities import get_im
			
			print "starting cter"
			
			# get the current micrograph
			image_name = self.micrograph_filepath
			img = get_im(image_name)
			
			try:
				ctf_window_size  = int(self.ctf_window_size.text())
				input_pixel_size = float(self.ctf_pixel.text())
				ctf_edge_size    = int(self.ctf_edge_size.text())
				ctf_overlap_size = int(self.ctf_overlap_size.text())
				ctf_volt         = float(self.ctf_volt.text())
				ctf_cs           = float(self.ctf_cs.text())
				ctf_ampcont      = float(self.ctf_ampcont.text())
				ctf_kboot        = int(self.ctf_kboot.text())
				
			except ValueError,extras:
				# conversion of a value failed!
				print "integer conversion failed."
				if not(extras.args is None):
					print extras.args[0]
				return
			except:
				print "error"
				return
			
			fname, fext = os.path.splitext(image_name)
			outpwrot = 'pwrot_%s'%fname
			outpartres = 'partres_%s'%fname
			
			if os.path.exists(outpwrot) or os.path.exists(outpartres):
				print "Please remove or rename %s and or %s"%(outpwrot,outpartres)
				return
			
			from morphology import cter
			defocus, ast_amp, ast_agl, error_defocus, error_astamp, error_astagl = cter(None, outpwrot, outpartres, None, None, ctf_window_size, voltage=ctf_volt, Pixel_size=input_pixel_size, Cs = ctf_cs, wgh=ctf_ampcont, kboot=ctf_kboot, MPI=False, DEBug= False, overlap_x = ctf_overlap_size, overlap_y = ctf_overlap_size, edge_x = ctf_edge_size, edge_y = ctf_edge_size, guimic=image_name)
		
			self.estdef.setText(str(defocus))
			self.estdef.setEnabled(False)
			
			self.astamp.setText(str(ast_amp))
			self.astamp.setEnabled(False)
			
			self.astagl.setText(str(ast_agl))
			self.astagl.setEnabled(False)
			
			self.deferr.setText(str(error_defocus))
			self.deferr.setEnabled(False)
			
			self.astamperr.setText(str(error_astamp))
			self.astamperr.setEnabled(False)
			
			self.astaglerr.setText(str(error_astagl))
			self.astaglerr.setEnabled(False)
			
			# XXX: wgh?? amp_cont static to 0?
			# set image properties, in order to save ctf values
			from utilities import set_ctf
			set_ctf(img, [defocus, ctf_cs, ctf_volt, input_pixel_size, 0, ctf_ampcont, ast_amp, ast_agl])
			# and rewrite image 
			img.write_image(image_name)
			print [defocus, ctf_cs, ctf_volt, input_pixel_size, 0, ctf_ampcont, ast_amp, ast_agl]
			
			print [defocus, ctf_cs, ctf_volt, input_pixel_size, 0, ctf_ampcont, ast_amp, ast_agl]
			print "CTF estimation using CTER done."
						
		def color_boxes(self):
			"""
			Sets the colors of the boxes, with the current box being colored differently from the other boxes.
			"""
			emshapes_dict = self.main_image.get_shapes()
			for key in emshapes_dict:
				shape = emshapes_dict.get(key).shape
				for i in range(3):
					shape[i+1] = self.color[i]
			current_shape = emshapes_dict.get(self.current_boxkey)
			if current_shape:
				for i in range(3):
					current_shape.shape[i+1] = self.selected_color[i]
			self.main_image.shapechange=1
			self.main_image.updateGL()
		def display_helix(self, helix_emdata):
			"""
			launches or updates an EMImage2DWidget to display helix_emdata
			@param helix_emdata: an EMData object that stores the image data for a helix
			"""
			self.color_boxes()
			if not self.helix_viewer:
				self.helix_viewer = EMImage2DWidget(application=get_application())
				#self.helix_viewer.setWindowTitle("Current Helix")
				self.helix_viewer.resize(300,800)
				self.helix_viewer.set_scale(1)
			QtCore.QObject.connect(self.helix_viewer, QtCore.SIGNAL("module_closed"), self.helix_viewer_closed)
			self.helix_viewer.set_data(helix_emdata)
			self.helix_viewer.setWindowTitle("Current Helix: %d x %d pixels" % (helix_emdata["nx"], helix_emdata["ny"]) )
			get_application().show_specific(self.helix_viewer)
			self.helix_viewer.updateGL()
		def closeEvent(self, event):
			"""
			overwriting the default so this will close all windows
			"""
			if self.main_image:
				self.main_image.closeEvent(event)
			if self.helix_viewer:
				self.helix_viewer.closeEvent(event)
			event.accept()
	
		def generate_emshape_key(self):
			"""
			creates a unique key for a new "rectline" EMShape, which is used for boxing a helix
			@return: a string that is the key for the new "rectline" EMShape
			"""
			i = self.counter.next()
			return "rectline%i" % i
		def get_width(self):
			"""
			returns the current width for the helices
			"""
			return self.box_width_spinbox.value()
		def helix_viewer_closed(self):
			"""
			This should execute when self.helix_viewer is closed.
			"""
			self.helix_viewer = None
		def load_boxes(self):
			"""
			load boxes from a file selected in a file browser dialog
			"""
			path = QtGui.QFileDialog.getOpenFileName(self, self.tr("Open Box Coordinates File"), "", self.tr("Boxes (*.txt *.box)"))
			path = str(path)
			coords_list = load_helix_coords(path)
			
			if self.main_image.shapes!=None and len(self.main_image.shapes)>0 :
				keep_boxes_msgbox = QtGui.QMessageBox()
				keep_boxes_msgbox.setText(self.tr("Keep current boxes?"))
				keep_boxes_msgbox.setInformativeText(self.tr("Do you want to keep your current boxes?"))
				keep_boxes_msgbox.setStandardButtons(QtGui.QMessageBox.No | QtGui.QMessageBox.Yes)
				keep_boxes_msgbox.setDefaultButton(QtGui.QMessageBox.Yes)
				keep_current_boxes = keep_boxes_msgbox.exec_()
		
				if keep_current_boxes == QtGui.QMessageBox.No:
					self.main_image.shapes = EMShapeDict()
					self.set_db_item("helixboxes", [])
					self.helices_dict = {}
					if self.helix_viewer:
						self.display_helix(EMData(10,10))
			
			for coords in coords_list:
				emshape = EMShape(["rectline", self.color[0], self.color[1], self.color[2], coords[0], coords[1], coords[2], coords[3], coords[4], 2])
				key = self.generate_emshape_key()
				self.main_image.add_shape(key, emshape)
				helix = get_helix_from_coords(self.main_image.get_data(), *coords)
				helix["ptcl_source_image"] = self.micrograph_filepath
				self.helices_dict[coords] = helix
				self.add_box_to_db(coords)
	
			self.main_image.updateGL()
			self.update_micrograph_table()
		
		def load_micrograph(self, micrograph_emdata):
			"""
			This displays a micrograph in self.main_image, resetting member variables that refer to other micrographs.
			If self.main_image == None, a new EMImage2DWidget is instatiated.
			@param micrograph_emdata: the EMData object that holds the micrograph to display
			"""
			self.edit_mode = None #Values are in {None, "new", "move", "2nd_point", "1st_point", "delete"}
			self.current_boxkey = None
			self.initial_helix_box_data_tuple = None
			self.click_loc = None #Will be (x,y) tuple
			
			if not self.main_image:
				self.main_image = EMImage2DWidget(application=self.app)
				QtCore.QObject.connect(self.main_image,QtCore.SIGNAL("module_closed"), self.main_image_closed)
				QtCore.QObject.connect( self.main_image, QtCore.SIGNAL("mousedown"), self.mouse_down)
				QtCore.QObject.connect( self.main_image, QtCore.SIGNAL("mousedrag"), self.mouse_drag)
				QtCore.QObject.connect( self.main_image, QtCore.SIGNAL("mouseup"), self.mouse_up)
			self.main_image.set_data( micrograph_emdata, self.micrograph_filepath )
			self.main_image.shapes = EMShapeDict()
			self.main_image.shapechange=1
			get_application().show_specific(self.main_image)
			self.helices_dict = db_get_helices_dict(self.micrograph_filepath) #Will be like {(x1,y1,x2,y2,width): emdata}
			
			if self.get_db_item("helixboxes") == None:
				self.set_db_item("helixboxes", [])
			else:
				boxList = self.get_db_item("helixboxes")
				for box_coords in boxList:
					key = self.generate_emshape_key()
					emshape_list = ["rectline"]
					emshape_list.extend(list(self.color))
					emshape_list.extend(list(box_coords))
					emshape_list.append(2)
					emshape = EMShape( emshape_list )
					self.main_image.add_shape(key, emshape)
				self.main_image.updateGL()
			
			qual = self.get_image_quality()
			if qual:
				self.img_quality_combobox.setCurrentIndex( qual )
			else:
				self.img_quality_combobox.setCurrentIndex( 2 )
			
			self.main_image.optimally_resize()
			
			width = self.box_width
			if self.helices_dict:
				first_coords = self.helices_dict.keys()[0]
				width = first_coords[4]
			self.box_width_spinbox.setValue(width)
		def main_image_closed(self):
			"""
			This should execute when self.main_image is closed.
			"""
			if self.helix_viewer:
				self.helix_viewer.close()
			self.main_image = None
		def micrograph_table_selection(self, row, column):
			"""
			When a new cell in the micrograph table is selected that is for a different 
			micrograph than the current one, the new micrograph is loaded.
			"""
			if row >= 0: #will be -1 when all rows are removed as first step of self.update_micrograph_table()
				new_filepath = str( self.micrograph_table.item(row,0).toolTip() )
				assert new_filepath in self.micrograph_filepath_set
				if new_filepath != self.micrograph_filepath:
					self.micrograph_filepath = new_filepath
					micrograph = EMData(self.micrograph_filepath)
					
					if self.invert_contrast:
						from utilities import info, model_blank
						from EMAN2 	   import Util
						# invert contrast of micrograph so average remains unchanged
						print "Inverting contrast of micrograph"
						mnx = micrograph.get_xsize()
						mny = micrograph.get_ysize()
						sttt = info(micrograph)
						avgimg = model_blank(mnx, ny=mny, bckg=sttt[0])
	
						Util.sub_img(micrograph, avgimg) # subtract average
						Util.mul_scalar(micrograph, -1.0) # multiply by -1
						Util.add_img(micrograph, avgimg) # add back average
						sttt2 = info(micrograph)
						assert(abs(sttt[0] - sttt2[0])<0.0001), "Assert failed: average of micrograph should remain same after contrast inversion!"
						
					self.load_micrograph(micrograph)
		def open_micrograph(self):
			"""
			loads a file browser to select a micrograph or multiple microgrpahs to add to the micrograph table
			"""
			selector = EMSelectorDialog(single_selection=False,save_as_mode=False)
			new_micrographs = selector.exec_()
			if isinstance(new_micrographs, str): #Just one file was selected
				if sys.version_info >= (2, 6):
					new_micrographs = os.path.relpath(new_micrographs) #os.path.relpath is new in Python 2.6
				self.micrograph_filepath_set.add(new_micrographs)
			else: #Multiple files were selected
				if sys.version_info >= (2, 6):
					new_micrographs = [os.path.relpath(path) for path in new_micrographs] #os.path.relpath is new in Python 2.6
				self.micrograph_filepath_set.update(set(new_micrographs))
			self.update_micrograph_table()
		def update_micrograph_table(self):
			"""
			sets the micrograph table cells to the data from self.micrograph_filepath_set
			"""
			self.micrograph_table.setRowCount( 0 )
			i = 0
			for micrograph_filepath in self.micrograph_filepath_set:
				file = os.path.basename(micrograph_filepath)
				micrograph = os.path.splitext(file)[0]
				boxes = db_get_item(micrograph_filepath, "helixboxes")
				if boxes:
					num_boxes = len(boxes)
				else:
					num_boxes = 0
				micrograph_item = QtGui.QTableWidgetItem(micrograph)
				micrograph_item.setToolTip(micrograph_filepath)
				num_boxes_item = QtGui.QTableWidgetItem(str(num_boxes))
				self.micrograph_table.insertRow(i)
				self.micrograph_table.setItem(i,0, micrograph_item)
				self.micrograph_table.setItem(i,1, num_boxes_item)
				if micrograph_filepath == self.micrograph_filepath:
					self.micrograph_table.setCurrentCell(i,0)
				i+=1
			self.micrograph_table.sortItems(0)
				
		def width_changed(self, width):
			"""
			updates the widths of the boxed helices when the user changes the width to use for helices
			"""
			if width < 1:
				return
			
			#resize current boxes
			#TODO: this is similar to part of self.mouse_up ==> make both methods call a function with common code
			shapes = self.main_image.get_shapes() #an EMShapeDict of EMShapes
			for box_key in shapes.keys():
				old_emshape = shapes.get(box_key)
				old_coords = old_emshape.getShape()[4:9]
				new_coords = (old_coords[0], old_coords[1], old_coords[2], old_coords[3], width)
				helix = get_helix_from_coords( self.main_image.get_data(), *new_coords )
				helix["ptcl_source_image"] = self.micrograph_filepath
							
				self.remove_box_from_db(old_coords)
				self.add_box_to_db(new_coords)
				self.helices_dict.pop(tuple(old_coords))
				self.helices_dict[new_coords] = helix
							
				new_emshape = EMShape( ["rectline", self.color[0], self.color[1], self.color[2], new_coords[0], new_coords[1], new_coords[2], new_coords[3], new_coords[4], 2] )
				shapes[box_key] = new_emshape
				
			self.main_image.shapechange=1
			self.main_image.updateGL()
			
			if self.helix_viewer:
				self.display_helix(EMData(10,10))
	#    def write_coords(self):
	#        """
	#        Save boxed helix coordinates to tab separated text file.
	#        """
	#        (micrograph_dir, micrograph_filename) = os.path.split(self.micrograph_filepath)
	#        default_filename = os.path.splitext(micrograph_filename)[0] + "_boxes.txt"
	#        file_dlg = QtGui.QFileDialog(self,self.tr("Save Helix Coordinates"), micrograph_dir)
	#        file_dlg.setAcceptMode(QtGui.QFileDialog.AcceptSave)
	#        file_dlg.selectFile(default_filename)
	#        if file_dlg.exec_():
	#            file_path = file_dlg.selectedFiles()[0]
	#            file_path = str(file_path)
	#            print file_path
	#            save_helix_coords(self.helices_dict.keys(), file_path)
		def write_images(self):
			"""
			Load EMWriteHelixFilesDialog to save helices, and particles to image files. 
			"""
			self.write_helix_files_dlg = EMWriteHelixFilesDialog(self,self.saveext)
			self.write_helix_files_dlg.setModal(True)
			self.write_helix_files_dlg.show()
			
		def get_db_item(self, key):
			"""
			gets the value stored in the sxhelixboxer database for the specified key and the current micrograph 
			"""
			db_name = SXHELIXBOXER_DB + key
			db = db_open_dict(db_name)
			val = db[self.micrograph_filepath]
			db_close_dict(db_name)
			return val
		def remove_db_item(self, key):
			"""
			removes the key and its value from the sxhelixboxer database for the current micrograph
			"""
			db_name = SXHELIXBOXER_DB + key
			db = db_open_dict(db_name)
			db.pop(key)
		def set_db_item(self, key, value):
			"""
			sets the value stored in the sxhelixboxer database for the specified key and the current micrograph 
			"""
			db_name = SXHELIXBOXER_DB + key
			db = db_open_dict(db_name)
			db[self.micrograph_filepath] = value
			db_close_dict(db_name)
		def get_image_quality(self):
			"""
			gets the value stored in the sxhelixboxer database for image quality, which is the user's subjective
			evaluation of how good the micrograph is
			"""
			return self.get_db_item("quality")
		def set_image_quality(self, quality):
			"""
			sets the value stored in the sxhelixboxer database for image quality, which is the user's subjective
			evaluation of how good the micrograph is
			"""
			self.set_db_item("quality", quality)
		def add_box_to_db(self, box_coords):
			"""
			adds the coordinates for a helix to the sxhelixboxer database for the current micrograph
			"""
			assert len(box_coords) == 5, "box_coords must have 5 items"
			db = db_open_dict(SXHELIXBOXER_DB + "helixboxes")
			boxList = db[self.micrograph_filepath] #Get a copy of the db in memory
			boxList.append(tuple(box_coords))
			db[self.micrograph_filepath] = boxList #Needed to save changes to disk
		def remove_box_from_db(self, box_coords):
			"""
			removes the coordinates for a helix in the sxhelixboxer database for the current micrograph
			"""
			assert len(box_coords) == 5, "box_coords must have 5 items"
			db = db_open_dict(SXHELIXBOXER_DB + "helixboxes")
			boxList = db[self.micrograph_filepath] #Get a copy of the db in memory
			boxList.remove(tuple(box_coords))
			db[self.micrograph_filepath] = boxList #Needed to save changes to disk
	
		def mouse_down(self, event, click_loc):
			"""
			If the shift key is pressed and the click is inside a box, delete it.
			Otherwise, either create a new box or edit an existing one depending on click location.
			Imagine drawing two (infinite) lines through the long sides of each box.
			If the click is not between two of the lines for a box, we will create a new box.
			Then the behavior depends on distance from the shorter axis of symmetry--in other
			words, how far up or down the length of the box. Clicking in the middle 3/4 of the box
			(3/8 L from the shorter axis of symmetry) will result in moving the entire box.
			Clicking on a point betwen 3/8 L and 5/8 L from the shorter axis of symmetry
			results in moving that end of the box while keeping the midpoint of the other end fixed.
			
			@param event: the mouse click event that causes a box to be added, removed, or modified
			@param click_loc: the coordinates in image (not screen) pixels of the mouse click on the image
			"""
	
			self.click_loc = click_loc
			box_key = None
			
			if self.main_image.get_shapes(): #helix boxes already exist
				box_key = self.main_image.get_shapes().closest_collision(click_loc[0], click_loc[1], fuzzy=True)
				if event.modifiers()&QtCore.Qt.ShiftModifier:
					if not box_key:
						self.edit_mode = None #Nothing to delete
					else:
						box_key = self.main_image.get_shapes().closest_collision(click_loc[0], click_loc[1], fuzzy=False)
						self.edit_mode = "delete"
				else:
					if not box_key:
						self.edit_mode = "new"
					else:
						control_points = self.main_image.get_shapes().get(box_key).control_pts()
						closest_pt_ix = 0
						point = control_points[0]
						min_squared_dist = (click_loc[0] - point[0])**2 + (click_loc[1] - point[1])**2
						for i in (1,2):
							point = control_points[i]
							dist_squared = (click_loc[0] - point[0])**2 + (click_loc[1] - point[1])**2
							if dist_squared < min_squared_dist:
								min_squared_dist = dist_squared
								closest_pt_ix = i
						if closest_pt_ix == 0: #first endpoint
							self.edit_mode = "1st_point"
						elif closest_pt_ix == 1: #second endpoint
							self.edit_mode = "2nd_point"
						elif closest_pt_ix == 2: #midpoint
							self.edit_mode = "move"
						else:
							self.edit_mode = "error"
				
			else: #no boxes exist
				if event.modifiers()&QtCore.Qt.ShiftModifier: #nothing to delete
					self.edit_mode = None
				else:
					self.edit_mode = "new" #create new box
		
					
			
			if self.edit_mode == "new" or not self.edit_mode:
				self.current_boxkey = None
				self.initial_helix_box_data_tuple = None
			elif self.edit_mode == "delete":
				box_coords = self.main_image.get_shapes().get(box_key).getShape()[4:9]
				self.remove_box_from_db(box_coords)
				self.helices_dict.pop(tuple(box_coords))
				self.main_image.del_shape(box_key)
				self.main_image.updateGL()
				self.current_boxkey = None
				(row, col) = (self.micrograph_table.currentRow(), 1)
				num_boxes = int(str( self.micrograph_table.item(row,col).text() ))
				self.micrograph_table.item(row,col).setText(str(num_boxes-1))
			else:
				self.current_boxkey = box_key
				self.initial_helix_box_data_tuple = tuple( self.main_image.get_shapes().get(box_key).getShape()[4:9] )
	
		def mouse_drag(self, event, cursor_loc):
			"""
			Boxes are deleted in mouse_down, and the decision of how to edit is made there.
			However, new boxes are made and existing boxes are edited here.
			@param event: the mouse click event that causes a box to be added, removed, or modified
			@param click_loc: the coordinates in image (not screen) pixels of the mouse click on the image
			"""
			
			if self.click_loc and self.edit_mode: #self.click_loc and self.edit_mode are set in mouse_down
				if self.edit_mode == "new":
					if self.click_loc[0] != cursor_loc[0] or self.click_loc[1] != cursor_loc[1]: #Don't make a zero-sized box
						self.current_boxkey = self.generate_emshape_key()                
						emshape_tuple = ( "rectline",self.color[0], self.color[1], self.color[2], 
											self.click_loc[0], self.click_loc[1], cursor_loc[0], cursor_loc[1], self.get_width(), 2 )
						
						emshape_box = EMShape(emshape_tuple)
						self.main_image.add_shape(self.current_boxkey, emshape_box)
						self.main_image.updateGL()
						self.initial_helix_box_data_tuple = emshape_tuple[4:9]
						self.edit_mode = "2nd_point"
						
						helix = get_helix_from_coords( self.main_image.get_data(), *self.initial_helix_box_data_tuple )
						helix["ptcl_source_image"] = self.micrograph_filepath
						self.display_helix(helix)
						(row, col) = (self.micrograph_table.currentRow(), 1)
						num_boxes = int(str( self.micrograph_table.item(row,col).text() ))
						self.micrograph_table.item(row,col).setText(str(num_boxes+1))
					
				elif self.edit_mode == "delete":
					pass
				else:
					first = self.initial_helix_box_data_tuple[:2]
					second = self.initial_helix_box_data_tuple[2:4]
					width = self.initial_helix_box_data_tuple[4]
					move = (cursor_loc[0] - self.click_loc[0], cursor_loc[1]-self.click_loc[1])
	
					if self.edit_mode == "move":
						first = (move[0]+first[0], move[1]+first[1])
						second = (move[0]+second[0], move[1]+second[1])
					elif self.edit_mode == '1st_point': #move first point
						first = (move[0]+first[0], move[1]+first[1])
					elif self.edit_mode == "2nd_point":
						second = (move[0]+second[0], move[1]+second[1])
					
					box = self.main_image.get_shapes().get(self.current_boxkey)
					box.getShape()[4] = first[0]
					box.getShape()[5] = first[1]
					box.getShape()[6] = second[0]
					box.getShape()[7] = second[1]
					self.main_image.shapechange=1
					self.main_image.updateGL()
					
					box_coords = tuple( box.getShape()[4:9] )
					helix = get_helix_from_coords( self.main_image.get_data(), *box_coords )
					helix["ptcl_source_image"] = self.micrograph_filepath
					self.display_helix(helix)
	
		def mouse_up(self, event, cursor_loc):
			"""
			Once the mouse button comes back up, creating a new box, or editing
			an existing box is complete, so we need only clear variables relevant
			to creating or editing boxes, and get the image data from the boxed area.
			@param event: the mouse click event that causes a box to be added, removed, or modified
			@param click_loc: the coordinates in image (not screen) pixels of the mouse click on the image
			"""
	
			if self.current_boxkey and self.edit_mode != "delete":
				if self.helices_dict.has_key(self.initial_helix_box_data_tuple):
					self.helices_dict.pop(self.initial_helix_box_data_tuple)
				if self.initial_helix_box_data_tuple in self.get_db_item("helixboxes"):
					self.remove_box_from_db(self.initial_helix_box_data_tuple)
				box = self.main_image.get_shapes().get(self.current_boxkey)
				box_coords = tuple( box.getShape()[4:9] )
				helix = get_helix_from_coords( self.main_image.get_data(), *box_coords )
				helix["ptcl_source_image"] = self.micrograph_filepath
				self.helices_dict[box_coords] = helix
				
				self.add_box_to_db(box_coords)
				self.display_helix(helix)
			
			self.click_loc = None
			self.edit_mode = None
			self.current_boxkey = None #We are done editing the box
			self.initial_helix_box_data_tuple = None

def windowallmic(dirid, micid, micsuffix, outdir, pixel_size, boxsize=256, minseg = 6, outstacknameall='bdb:data', \
				hcoords_dir = "", hcoords_suffix = "_boxes.txt", ptcl_dst=-1, inv_contrast=False, new_pixel_size=-1, rmax = -1.0, freq = -1, debug = 1, \
				do_rotation = True, do_gridding=True, topdir = None, importctf=None, limitctf=None, cterr = None):
	'''
	
	Windows segments from helices boxed from micrographs. 
	
	Input
		
		dirid: A string for identifying directories containing relevant micrographs.
		 
			   Any directory containing dirid as a contiguous string will be searched
		       for micrographs. 
		       
		       These micrographs are assumed to be those which were used to box the helices
		       which are to be windowed.
		       
		       The pixel size of the micrographs should be pixel_size.
		       
		micid: A string for identifying the name (minus extension) of relevant micrographs.
		       
		micsuffix: A string denoting micrograph type. Currently only handles suffix types, i.e. 'hdf', 'ser' etc.
		
		outdir: Output directory to be created in EACH micrograph directory.
		        
		        The segments windowed from the helices boxed from the micrographs in the directory,
		        and possibly resampled micrographs (if pixel size changed), will be put here.
		
		pixel_size: The pixel size of the micrographs which were used to box the helices.
		
		boxsize: Dimension of square window size. 
				 
		outstacknameall: File name with full path and type (only handles bdb and hdf right now) 
						 under which ALL windowed segments from ALL micrograph directories 
						 will be saved, e.g. 'bdb:/home/project/adata' or '/home/project/adata.hdf'
		
		hcoords_suffix: String identifier which when concatenated with a micrograph name (minus extension) gives the name of the text file 
						containing coordinates of ALL helices boxed from the micrograph.
						
					    If there is no such file, helices boxed from the micrograph will not be windowed.
					    
						Default is '_boxes.txt', so if mic0.hdf is a micrograph, then the text file containing 
						coordinates of the helices boxed in it is mic0_boxes.txt.
						
						The coordinate file is assumed to be in the format used by sxhelixboxer, e.g.:
						
								x1-w/2           y1-w/2           w           w           -1
							    x2-w/2           y2-w/2           w           w           -2
								...
								
						where (x1, y1) and (x2, y2) are the coordinates on the micrograph for the helical axis endpoints, and w is the width of the helix boxes
		
		hcoords_dir: Full path name of directory containing coordinates of ALL helices boxed from the micrograph.
						
		ptcl_dst: Integer. Distance in pixels between adjacent squares windowed from a single boxed helix.
			      If ptcl_dst < 0, then the program will set it so ptcl_dst is ~ one rise in pixels: int( (dp/new_pixel_size) + 0.5)
		
		inv_contrast: True/False, default is False. If cryo, then set to true to invert contrast so particles show up bright against dark background. 
		
		new_pixel_size: Float. New target pixel size to which the micrograph should be resampled. 
		
		 			    Default is -1, in which case new_pixel_size is assumed to be same as pixel_size.
	
		rmax: Float. Radius of filament in Angstroms. 
		
		freq: Cut-off frequency at which to high-pass filter micrographs before windowing. 
		
		      Default is -1, in which case, the micrographs will be high-pass filtered with cut-off frequency 1.0/segnx, where segnx is the target x dimension of the segments.
		      
		debug: If 1, then do NOT delete output directories where intermediate files are stored. If 0, then delete the output directories. Default is 1.
		
		topdir: Directory containing all the relevant micrograph directories. If not specified, it will be taken as the directory from which windowing program is invoked.
		
	Output
	
		outdir: In each micrograph directory, the program will write the stack of segments windowed from all micrographs in 
		the directory to suddirectory outdir under the file name 'bdb:data'.
		
		outstacknameall: The program will concatenate segment stacks from ALL micrograph directories and write them 
						 to the file name outstacknameall, e.g. outstacknameall='bdb:adata' or outstacknameall='adata.hdf'
						 in the directory where windowallmic is invoked.
						 
						 Only hdf and bdb file types are currently handled.
		
	Example of use: 
		
		In directory mic, there is a micrograph mic0.hdf.
		
		The full path for the directory is: /Users/project/mic
		
		If topdir is not specified, then windowallmic must be invoked in the directory that contains the relevant
		micrograph directories, i.e., /Users/project
		
		Otherwise, if topdir is specified, e.g. as /Users/project, then the windowing program can be invoked from anywhere.
		
		The original pixel size is 1.2, and the radius of the helical filaments 
		is 50*1.2 Angstroms.
		
		The micrograph has already been boxed and mic0_boxes.txt contains the box coordinates
		as output by sxhelixboxer. 
		
		
		The following call to windowallmic will
		write bdb:adata to the directory where windowallmic is invoked, where bdb:adata is a stack
		of segments windowed from the boxed helix in mic0.hdf.
		
		The pixel size of the segments is 1.84 and the box size is 200 by 200.
		
		windowallmic(dirid='mic', micid='mic', micsuffix='hdf', outdir='out',  pixel_size=1.2, boxsize=200, minseg = 6, outstacknameall='bdb:adata', hcoords_suffix = "_boxes.txt", ptcl_dst=15, inv_contrast=False, new_pixel_size=1.84, rmax = 60.0, topdir='/Users/project')
	'''
	import os
	from utilities      import print_begin_msg, print_end_msg, print_msg
	from sxhelixboxer	import windowmic
	from EMAN2 	        import EMUtil, Util
	
	print_begin_msg("windowallmic\n")
	
	if not(do_rotation):
		do_gridding = False
		
	if not(outstacknameall[0:4] == 'bdb:' or outstacknameall[-3:] == 'hdf'):
		print "%s must be in bdb or hdf format"%outstacknameall
		return

	if freq < 0: freq = 1.0/boxsize
		
	if micsuffix[0] == '.': micsuffix = micsuffix[1:]
	
	if new_pixel_size < 0: new_pixel_size = pixel_size
	
	if rmax < 0:  rmaxp = int(boxsize/2 - 2)
	else:         rmaxp = int( (rmax/new_pixel_size)  + 0.5)
	
	# set rmax in pixels to default if user input rmax is greater than half the box size
	if 2*rmaxp > boxsize:
		print "ERROR...The segment size should be no less than twice rmax in pixels. Current segment size is %d and twice rmax in pixels is %d."%(boxsize, 2*rmaxp)
		return
			
	# Calculate distance between adjacent squares as ~1 rise in pixels if not set by user
	if ptcl_dst < 0:
		print "ERROR...Distance between segments (ptcl_dst) not specified."
		return

	if topdir == None:
		topdir = os.getcwd()
	flist = os.listdir(topdir)
	outdirlist = [] # List of output directories to create (after checking they do not already exist)
	micdirlist = [] # List of all directories with dirid
	# Create output directory in each micrograph directory, exit with error if one already exists
	for i1, v1 in enumerate(flist):
		topv1 = os.path.join(topdir,v1)
		if not(os.path.isdir(topv1)):
			continue
		if v1.find(dirid) < 0:
			continue
		micdirlist.append(topv1)
		# v1 is a micrograph directory, create directory named outdir in v1
		coutdir = os.path.join(topv1, outdir)
		if os.path.exists(coutdir):
			print 'Output directory %s  exists, please change the name and restart the program'%coutdir
			return
		outdirlist.append(coutdir)
	for coutdir in outdirlist:
		print_msg("Creating output directory %s\n"%coutdir)
		os.mkdir(coutdir)
	cutoffhistogram = []       #@ming
	lenmicnames = 0	
	for v1 in micdirlist:
		# window all micrographs in directory v1 with micid
		flist2 = os.listdir(v1)
		coutdir = os.path.join(v1, outdir)
		# sort flist2 using case insensitive string comparison
		flist2.sort(key=str.lower)
		nfiles = len(flist2)
		print_msg('Sorted file list in %s:\n'%v1)
		for iii in xrange(nfiles):
			print_msg('%s,'%flist2[iii])
		print_msg('\n')
		for i2, v2 in enumerate(flist2):
			filename, fext = os.path.splitext(v2)
			if fext[1:] != micsuffix:
				continue
			if filename.find(micid)>-1:
				# v2 is a micrograph to window IF text file containing box coordinates exists
				hcoordsname = filename + hcoords_suffix
				if len(hcoords_dir) > 0: hcoordsname = os.path.join(hcoords_dir, hcoordsname)
				else:                    hcoordsname = os.path.join(v1, hcoordsname)
				# If any helices were boxed from this micrograph, say mic0, then ALL the helix coordinates should be saved under mic0 + hcoords_suffix
				# For example, if using default sxhelixboxer naming convention, then coordinates of all helices boxed in mic0 would be in mic0_boxes.txt
				if( os.path.exists(hcoordsname) ):
					micname = os.path.join(v1, v2)
					lenmicnames += 1
					print_msg("\n\nPreparing to window helices from micrograph %s with box coordinate file %s\n\n"%(micname, hcoordsname))
					#windowmic(outstacknameall, coutdir, micname, hcoordsname, pixel_size, boxsize, ptcl_dst, minseg, inv_contrast, new_pixel_size, rmaxp, freq, do_rotation, do_gridding, importctf, cterr)
					windowmic(outstacknameall, v1, coutdir, micname, hcoordsname, pixel_size, boxsize, ptcl_dst, minseg, inv_contrast, new_pixel_size, rmaxp, freq, do_rotation, do_gridding, importctf, limitctf, cterr, cutoffhistogram)  ##changed by @ming 
	
	if len(cutoffhistogram) > 0:		#@ming
		lhist = 3
		if len(cutoffhistogram) >= lhist:
			from statistics import hist_list
			region,hist = hist_list(cutoffhistogram,lhist)	
			msg = "      Histogram of cut off frequencies\n      ERROR       number of frequencies\n"
			print_msg(msg)
			for lhx in xrange(len(lhist)):
				msg = " %10.3f     %7d\n"%(region[lhx], hist[lhx])
				print_msg(msg)
		print_msg('The percentage of micrographs filtered by the cutoff frequency: %6f\n' % (len(cutoffhistogram)*1.0/lenmicnames))		
	# If not debug mode, then remove all output directories 				
	if debug == 0:
		from subprocess import call
		for coutdir in outdirlist:
			cmd = "rm -ir %s"%coutdir
			print_msg("cmd: %s"%cmd)
			call(cmd, shell=True)

def windowmic(outstacknameall, micpath, outdir, micname, hcoordsname, pixel_size, boxsize, ptcl_dst, minseg, inv_contrast, new_pixel_size, rmaxp, freq, do_rotation, do_gridding, importctf, limitctf, cterr, cutoffhistogram):
	'''
	
	INPUT
			outstacknameall: File name with full path and type (only handles bdb and hdf right now) 
						 under which ALL windowed segments from ALL micrograph directories 
						 will be saved, e.g. 'bdb:adata' or 'adata.hdf'
						 
			outdir: Full path name of output directory in which to put segment stack. 
			
			micname: String. Full path name of micrograph 
			
			hcoordsname: String. Full path name of file contaning coordinates of boxed helices (file assumed to be in format used by sxhelixboxer).
			
			pixel_size: The pixel size of the micrographs in which the helices were boxed.
			
			boxsize: Integer. x and y-dimension of segments to be windowed.
			
			ptcl_dst: Integer. Distance in pixels between adjacent squares windowed from a single boxed helix.
				          
			inv_contrast: True/False, default is False. If cryo, then set to true to invert contrast so particles show up bright against dark background. 
			
			new_pixel_size: Float. New target pixel size to which the micrograph should be resampled. 
			
			rmaxp: Integer. Radius of filament in pixels. 
			
			freq: Cut-off frequency at which to high-pass filter micrographs before windowing. 
	
	OUTPUT
	
	       A stack of windowed segments (in file format bdb) corresponding to EACH boxed helix in the input micrograph.
	       
	       Example: If the name of the micrograph is mic0.hdf, and there are two filaments boxed in mic0.hdf, then the program will write
	       		    two stacks of segments to outdir: mic0_abox_0.hdf	and mic0_abox_1.hdf	 
	       	    
	'''
	from utilities    import pad, model_blank, read_text_row, get_im, print_msg
	from fundamentals import ramp, resample
	from filter	  	  import filt_gaussh,filt_tanl 
	from pixel_error  import getnewhelixcoords
	from EMAN2 	      import EMUtil, Util
	from subprocess   import call
	
	# micname is full path name
	# smic[-1] is micrograph name minus path
	smic = micname.split('/')
	# filename is name of micrograph minus the path and extension
	filename = (smic[-1].split('.'))[0]

	if importctf:	
		ctfs = read_text_row(importctf)
		nx = True
		for i in xrange(len(ctfs)):
			smic = ctfs[i][-1].split('/')
			ctfilename = (smic[-1].split('.'))[0]
			if(ctfilename == filename):
				ctfs = ctfs[i]
				nx = False
				break
		if nx:
			print "Micrograph %s"%filename,"  not listed in CTER results, skipping ...."
			return
		if(ctfs[8]/ctfs[0] > cterr[0]):
			print_msg('Defocus error %f exceeds the threshold. Micrograph %s rejected.\n'%(ctfs[8]/ctfs[0], filename))
			return
		if(ctfs[10] > cterr[1] ):
			ctfs[6] = 0.0      ##astigmatism amplitude     comment by@ming
			ctfs[7] = 0.0      ##astigmatism angle         comment by@ming

	dummy = EMData()
	dummy.read_image(micname, 0, True)
	l = dummy.get_attr_dict()

	## Cut off frequency components higher than CTF limit   @ming
	img = get_im(micname)
	from morphology import ctflimit
	if limitctf:
# 			Cut off frequency components higher than CTF limit 
		q1, q2 = ctflimit(boxsize,ctfs[0],ctfs[1],ctfs[2],new_pixel_size)
		# This is absolute frequency of the CTF limit in the scale of original micrograph
		q1 = (ctfs[3] / new_pixel_size) * q1/float(boxsize)
		if q1 < 0.5:          #@ming
			img = filt_tanl(img, q1, 0.01)
			cutoffhistogram.append(q1)
	
	if new_pixel_size != pixel_size:
		# Resample micrograph, map coordinates, and window segments from resampled micrograph using new coordinates
		# Set ctf along with new pixel size in resampled micrograph
		# set hcoordsname to name of new coordinates file, and set micname to resampled micrograph name
		print_msg('Resample micrograph to pixel size %f and window segments from resampled micrograph\n'%new_pixel_size)
		resample_ratio = pixel_size/new_pixel_size
		# after resampling by resample_ratio, new pixel size will be pixel_size/resample_ratio = new_pixel_size
		nx = img.get_xsize()
		ny = img.get_ysize()
		img = resample(img, resample_ratio)
		micname = os.path.join(outdir,'resampled_%s'%smic[-1])
		img.write_image(micname)
		if importctf: ctfs[3] = new_pixel_size
		smic = micname.split('/')
		# filename is name of micrograph minus the path and extension
		filename = (smic[-1].split('.'))[0]
		# now need to get new coordinates file and set hcoordsname to that
		hcoordsname = getnewhelixcoords(hcoordsname, outdir, resample_ratio,nx,ny, newpref="resampled_", boxsize=boxsize)

	imgs_0 = filename+"_abox" # Base name for the segment stack corresponding to each boxed helix. If imgs_0='mic0_abox', then windowed segments from first windowed helix would be 'mic0_abox_0.hdf', the second 'mic0_abox_1.hdf' etc
	fimgs_0 = os.path.join(outdir, imgs_0 + ".hdf") # This has to be hdf, sxhelixboxer cannot write windowed out segments to bdb

	ptcl_images  = "   --ptcl-images="+ fimgs_0

	# Name of file under which coordinates of segments windowed from ALL helices boxed from the micrograph will be saved
	fptcl_coords = os.path.join(outdir, filename + "_helix_ptcl_coords.txt") 
	ptcl_coords  = "   --ptcl-coords="+ fptcl_coords
	
	#smic[-1] is micrograph name minus path
	
	tmpfile = os.path.join(outdir,'filt_%s.hdf'%filename)
	imgmic  = get_im(micname)
	
	if inv_contrast:
		stt = Util.infomask(imgmic, None, True)
		Util.mul_scalar(imgmic, -1.0) # multiply by -1
		imgmic += 2*stt[0]

	imgmic.write_image(tmpfile)

	# Set box coordinates in sxhelixboxer database
	db_load_helix_coords(tmpfile, hcoordsname, False, boxsize)
	
	micrograph_filename = os.path.join(micpath,'%s.hdf'%filename)                     ##@ming
	print "mic file name=%s"%micrograph_filename
	print "mic file path=%s"%micpath
	micrograph_name = os.path.splitext( micrograph_filename )[0]        ##@ming 
	helix_filepath = "%s_helix.hdf" % ( os.path.join(micpath, filename) )   ##@ming
	helices_dict = win_get_helices_dict(tmpfile, boxsize)                ##@ming

    ## set a tag to save or unsave filaments  @ming 
	#i = 0
	#for coords in helices_dict:
	#	helix = helices_dict[coords]
	#	save_helix(helix, helix_filepath, i)
	#	i+=1
		
	db_save_particle_coords(tmpfile, fptcl_coords, ptcl_dst, boxsize, boxsize, do_rotation)
	db_save_particles(tmpfile, fimgs_0, ptcl_dst, boxsize, boxsize, do_rotation, True, do_gridding, "multiple", do_filt = True, filt_freq = freq)
	os.remove(tmpfile)

	mask = pad(model_blank(rmaxp*2, boxsize, 1, 1.0), boxsize, boxsize, 1, 0.0)

	a = read_text_row(hcoordsname)
	if len(a)%2 != 0:
		print "Number of rows in helix coordinates file %s should be even!"%hcoordsname
		return
	nhelices = len(a)/2
	#try:      iseg = EMUtil.get_image_count(outstacknameall)
	#except:   iseg = 0
	if importctf:
		from utilities import generate_ctf
		ctfs = generate_ctf(ctfs)
	#for h in xrange(nhelices):
	h=0                                                           ## added by@ming
	for coords in helices_dict:                                   ## added by@ming
		helix = helices_dict[coords]
		if importctf:			                         ## added by@ming
			ctfs.dfang -= helix["astig_jiao"]          ## added by@ming

		ptcl_images  = imgs_0+"_%i.hdf"%h                         # This is what sxhelixboxer outputs, only 'hdf' format is handled.
		otcl_images  = "bdb:%s/QT"%outdir+ ptcl_images[:-4]
		ptcl_images  = os.path.join(outdir,ptcl_images)
		if( os.path.exists(ptcl_images) ):
			n1 = EMUtil.get_image_count(ptcl_images)
			if( n1 < minseg ):
				print_msg( "Filament %s has too few segments and was skipped\n"%otcl_images)
			else:
				print_msg( "otcl_images: %s\n"%otcl_images)
				print_msg( "ptcl_images: %s\n"%ptcl_images)
				for j in xrange(n1):
					prj = get_im(ptcl_images, j)
					prj = ramp(prj)
					stat = Util.infomask( prj, mask, False )
					prj -= stat[0]
					if importctf:
						prj.set_attr("ctf",ctfs)						
						prj.set_attr("ctf_applied", 0)
					prj.write_image(otcl_images, j)
					#prj.write_image(outstacknameall, iseg)
					#iseg += 1
		h+=1			
if __name__ == '__main__':
	main()
