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

# $Id$

# todo: verify the processors who have the same names in proc3d
#	   and proc2d have the same implementation
#
# todo: lp, hp, tlp vs apix

from EMAN2 import *
from optparse import OptionParser
import sys
from math import *
import os.path
import pyemtbx.options
from pyemtbx.options import intvararg_callback
from pyemtbx.options import floatvararg_callback

def print_iminfo(data, label):
	print "%s image : %dx%dx%d Mean=%1.3g Sigma=%1.3g Min=%1.3g Max=%1.3g" % \
	(label, data.get_xsize(), data.get_ysize(), data.get_zsize(),
	data.get_attr("mean"), data.get_attr("sigma"),
	data.get_attr("minimum"), data.get_attr("maximum"))

def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + """ [options] <inputfile> <outputfile>
	Generic 3-D image processing and file format conversion program. 
        All EMAN2 recognized file formats accepted (see Wiki for list).

        Examples:

        convert MRC format to HDF format:
        e2proc3d.py test.mrc test.hdf                   

        apply a 10 A low-pass filter to a volume and write output to a new file. 
        e2proc3d.py threed.hdf threed.filt.hdf --process=filter.lowpass.gauss:cutoff_freq=0.1

	extract a reconstruction from a refinement directory as an HDF file usable with Chimera
	e2proc3d.py bdb:refine_02#threed_filt_04 map_02_04.hdf

        'e2help.py processors -v 2' for a detailed list of available procesors

"""
	parser = OptionParser(usage)
	
	parser.add_option("--medianshrink", metavar="n", type="int", action="append", 
								help="Shrinks the image by integer n using median filter")

	parser.add_option("--meanshrink", metavar="n", type="int", action="append", 
								help="Shrinks the image by integer n using mean filter")
	
#    parser.add_option("--tomoshrink", metavar="n", type="int", action="append", 
#                                help="Mean shrinks the image but is careful of memory - reads small pixel blocks from disk and slowly builds up the result")

	parser.add_option("--scale", metavar="n", type="float", action="append",
								help="Rescales the image by 'n', generally used with clip option.")

	parser.add_option("--sym", dest = "sym", action="append", help = "Symmetry to impose - choices are: c<n>, d<n>, h<n>, tet, oct, icos")

	parser.add_option("--clip", metavar="x,y,z[,xc,yc,zc]", type='string', action="callback", callback=intvararg_callback, 
							help="Make the output have this size, no scaling. ")
	
	parser.add_option("--fftclip", metavar="x, y, z", type="string", action="callback", callback=floatvararg_callback,
								help="Make the output have this size, rescaling by padding FFT.")

	parser.add_option("--process", metavar="processor_name:param1=value1:param2=value2", type="string",
								action="append", help="apply a processor named 'processorname' with all its parameters/values.")

	parser.add_option("--apix", type="float", help="A/pixel for S scaling")
	
	parser.add_option("--origin", metavar="x, y, z", type="string", action="callback", callback=floatvararg_callback,
								help="Set the coordinates for the pixel (0,0,0).")

	parser.add_option("--mult", metavar="f", type="float", 
								help="Scales the densities by a fixed number in the output")
	
	parser.add_option("--multfile", type="string", action="append",
								help="Multiplies the volume by another volume of identical size. This can be used to apply masks, etc.")
	
	parser.add_option("--mrc16bit",  action="store_true", help="output as 16 bit MRC file")
	parser.add_option("--mrc8bit",  action="store_true", help="output as 8 bit MRC file")
	
	parser.add_option("--add", metavar="f", type="float", 
								help="Adds a constant 'f' to the densities")

	parser.add_option("--calcfsc", type="string", metavar="with input",
								help="Calculate a FSC curve between two models. Output is a txt file. This option is the name of the second volume.")


	parser.add_option("--calcsf", type="string", metavar="outputfile",
								help="Calculate a radial structure factor. Must specify apix.")

	parser.add_option("--setsf", type="string", metavar="inputfile",
								help="Set the radial structure factor. Must specify apix.")
	
	parser.add_option("--tophalf", action="store_true",
								help="The output only keeps the top half map")
	
	parser.add_option("--icos5fhalfmap", action="store_true",
								help="The input is the icos 5f top half map generated by the 'tophalf' option")
	
	parser.add_option("--outtype", metavar="image-type", type="string", 
								help="Set output image format, mrc, imagic, hdf, etc")

	parser.add_option("--first", metavar="n", type="int", default=0, 
								help="the first image in the input to process [0 - n-1])")

	parser.add_option("--trans", metavar="dx,dy,dz", type="string", default=0, help="Translate map by dx,dy,dz ")

	parser.add_option("--rot", metavar="az,alt,phi or convention,a1,a2,a3,a4", type="string", default=0, help="Rotate map using EMAN Euler angles z,x,z' or an arbitrary convention. NOTE, at the moment users may only specify az,alt,phi - this is a bug that will be resolved")

	parser.add_option("--last", metavar="n", type="int", default=-1, 
								help="the last image in the input to process")
	
	parser.add_option("--swap", action="store_true", help="Swap the byte order", default=False)	
	
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	append_options = ["clip", "fftclip", "process", "filter", "meanshrink", "medianshrink", "scale", "sym", "multfile"]
	
	optionlist = pyemtbx.options.get_optionlist(sys.argv[1:])
	
	(options, args) = parser.parse_args()
	
	if len(args) != 2:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
		sys.exit(1)
	
	infile = args[0]
	outfile = args[1]
	
	n0 = options.first
	n1 = options.last
	nimg = EMUtil.get_image_count(infile)
	
	index_d = {}
	for append_option in append_options:
		index_d[append_option] = 0
	
	if(n0 < 0 or n0 > nimg):
		print "Your first index is out of range, changed to zero"
		n0 = 0
	
	if(n1 == -1): 
		n1 = nimg-1
	elif(n1 > nimg-1):
		print "Your last index is out of range, changed to %d" % (nimg-1)
		n1 = nimg-1
	
	datlst = parse_infile(infile, n0, n1)
	
	logid=E2init(sys.argv)

	x = datlst[0].get_xsize()
	y = datlst[0].get_ysize()
	z = datlst[0].get_zsize()
	
	xc = x/2
	yc = y/2
	zc = z/2
	
	nx = x
	ny = y
	nz = z

	#print_iminfo(datlst[0], "Original")
	
	apix = datlst[0]["apix_x"]	# default to apix_x from file
	if options.apix:
		apix = options.apix
		for data in datlst:
			data.set_attr('apix_x', apix)
			data.set_attr('apix_y', apix)
			data.set_attr('apix_z', apix)

	if not "outtype" in optionlist:
		optionlist.append("outtype")

	print "%d images, processing %d-%d......"%(nimg, n0, n1)
	#print 'start.....'
	for data in datlst:
		for option1 in optionlist:				
			if option1 == "origin":
				if(len(options.origin)==3):
					(originx, originy, originz) = options.origin
				else:
					print ''
					return
					
				data.set_xyz_origin(originx, originy, originz)

			elif option1 == "calcfsc" :
				datafsc=EMData(options.calcfsc)
				fsc = data.calc_fourier_shell_correlation(datafsc)
				third = len(fsc)/3
				xaxis = fsc[0:third]
				fsc = fsc[third:2*third]
				saxis = [x/apix for x in xaxis]
				Util.save_data(saxis[0],saxis[1]-saxis[0],fsc,args[1])
				print "Exiting after FSC calculation"
				sys.exit(0)

			elif option1 == "calcsf":
				dataf = data.do_fft()
				curve = dataf.calc_radial_dist(ny, 0, 0.5,True)
				curve=[i/(dataf["nx"]*dataf["ny"]*dataf["nz"]) for i in curve]
				Util.save_data(0, 1.0/(apix*2.0*ny), curve, options.calcsf);

			elif option1 == "setsf":
				sf=XYData()
				sf.read_file(options.setsf)
				
				data.process_inplace("filter.setstrucfac",{"apix":apix,"strucfac":sf})
				
				#dataf = data.do_fft()
				#curve = dataf.calc_radial_dist(ny, 0, 0.5,True)
				#filt=[]
				#norm=data["nx"]*data["ny"]*data["nz"]
				#for i,j in enumerate(curve):
					#if j==0 : filt.append(0.0)
					#else : filt.append(sqrt(norm*sf.get_yatx(i/(apix*2.0*dataf["ny"]))/j))
				#dataf.apply_radial_func(0,0.5/len(filt),filt)
				#data=dataf.do_ift()
				
				#curve2 = dataf.calc_radial_dist(ny, 0, 0.5,True)
				#plot((curve,curve2))

			elif option1 == "process":
				fi = index_d[option1]
				(filtername, param_dict) = parsemodopt(options.process[fi])
				data.process_inplace(filtername, param_dict)
				index_d[option1] += 1

			elif option1 == "mult":
				data.mult(options.mult)
				
			elif option1 == "multfile":
				mf=EMData(options.multfile[index_d[option1]],0)
				data.mult(mf)
				mf=None
				index_d[option1] += 1
				
#            elif option1 == "tomoshrink":
#                from e2tomoboxer import ShrunkenTomogram
#                st = ShrunkenTomogram(args[0])
#                st.set_cache_to_db(False)
#                tmp = st.get_image()

			elif option1 == "add":
				data.add(options.add)
			
			elif option1 == "trans":
				dx,dy,dz=options.trans.split(",")
				data.translate(float(dx),float(dy),float(dz))
			
			elif option1 == "rot":
				daz,dalt,dphi=options.rot.split(",")
				data.rotate(float(daz),float(dalt),float(dphi))
			elif option1 == "clip":
				if(len(options.clip) == 6):
					(nx, ny, nz, xc, yc, zc) = options.clip
				elif(len(options.clip) == 3):
					(nx, ny, nz) = options.clip
					xc = x/2
					yc = y/2
					zc = z/2
				else:
					print 'clip option takes either 3 or 6 arguments. --clip=x,y,z[,xc,yc,zc]'
					return

				if not (xc>=0 and yc>=0 and zc>=0 and xc<x and yc<y and zc<z):
					xc = x/2
					yc = y/2
					zc = z/2
					
				if x != nx or y != ny or z != nz:
					data.clip_inplace(Region(xc-nx/2, yc-ny/2, zc-nz/2, nx, ny, nz))
					index_d[option1] += 1

			elif option1 == "sym":
				sym = options.sym[index_d[option1]]
				xf = Transform()
				xf.to_identity()
				nsym=xf.get_nsym(sym)
				ref=data.copy()
				for i in range(1,nsym):
					dc=ref.copy()
					dc.transform(xf.get_sym(sym,i))
					data.add(dc)
				data.mult(1.0/nsym)	

			elif option1 == "scale":
				scale_f = options.scale[index_d[option1]]
				if scale_f != 1.0:
					data.scale(scale_f)
				index_d[option1] += 1

			elif option1 == "medianshrink":
				shrink_f = options.medianshrink[index_d[option1]]
				if shrink_f > 1:
					data.process_inplace("math.medianshrink",{"n":shrink_f})
					nx = data.get_xsize()
					ny = data.get_ysize()
					nz = data.get_zsize()
				index_d[option1] += 1

			elif option1 == "meanshrink":
				shrink_f = options.meanshrink[index_d[option1]]
				if shrink_f > 1:
					data.process_inplace("math.meanshrink",{"n":shrink_f})
					nx = data.get_xsize()
					ny = data.get_ysize()
					nz = data.get_zsize()
				index_d[option1] += 1

			elif option1 == "fftclip":
				if(len(options.fftclip)==3):
					(fnx, fny, fnz) = options.fftclip
				else:
					print 'fftclip option takes either 3 arguments. --fftclip=x,y,z'
					return
				
				fft = data.do_fft()
				padfft = fft.get_clip(Region(0, 0, 0, fnx+2, fny, fnz))
				data = padfft.do_ift()
				index_d[option1] += 1

			elif option1 == "icos5fhalfmap":
				print "not implemented yet"

			elif option1 == "tophalf":
				half = data.get_top_half()
				data = half
				
			elif option1 == "outtype":
				if not options.outtype:
					options.outtype = "unknown"
		
		#print_iminfo(data, "Final")
		
		if 'mrc8bit' in optionlist:
			data.write_image(outfile.split('.')[0]+'.mrc', -1, EMUtil.ImageType.IMAGE_MRC, False, None, EMUtil.EMDataType.EM_UCHAR, not(options.swap))
		elif 'mrc16bit' in optionlist:
			data.write_image(outfile.split('.')[0]+'.mrc', -1, EMUtil.ImageType.IMAGE_MRC, False, None, EMUtil.EMDataType.EM_SHORT, not(options.swap))
		else:
			data.write_image(outfile, -1, EMUtil.get_image_ext_type(options.outtype), False, None, EMUtil.EMDataType.EM_FLOAT, not(options.swap))
		
		for append_option in append_options:	#clean up the multi-option counter for next image
			index_d[append_option] = 0
		
	E2end(logid)


#parse_file() wil read the input image file and return a list of EMData() object
def parse_infile(infile, first, last):
	nimg = EMUtil.get_image_count(infile)
	
	if (nimg > 1):
		#print "it appears %s contains %d image" % (infile, nimg)
		d = EMData()
		d.read_image(infile, 0)
	
		x = d.get_xsize()
		y = d.get_ysize()
		z = d.get_zsize()
		if (z == 1):
			print "the images are 2D - I will now make a 3D image out of the 2D images"
			data = []
			return_data = EMData()
			return_data.set_size(x, y, nimg)
			for i in range(0, nimg):
				d.read_image(infile, i)
				return_data.insert_clip(d, (0, 0, i))
			data.append(return_data)	
			return data
		else:
			print "the image is a 3D stack - I will process images from %d to %d" % (first, last)
			data = []
			for i in xrange(first, last+1):
				d = EMData()
				d.read_image(infile, i)
				data.append(d)
			return data
	else:
		data = []
		d = EMData()
		d.read_image(infile,0)
		data.append(d)
		return data
	

if __name__ == "__main__":
	main()
