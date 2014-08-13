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

# $Id$

from EMAN2 import *
import sys
import os.path
import math
import random
import pyemtbx.options
import time

#constants

HEADER_ONLY=True
HEADER_AND_DATA=False

xyplanes = ['xy', 'yx']
xzplanes = ['xz', 'zx']
yzplanes = ['yz', 'yz']
threedplanes = xyplanes + xzplanes + yzplanes

# usage: e2proc2d.py [options] input output

def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + """ [options] <inputfile> <outputfile>

	Generic 2-D image processing and file format conversion program. Acts on stacks of 2-D images
	(multiple images in one file). All EMAN2 recognized file formats accepted (see Wiki for list).

	Examples:

	convert IMAGIC format test.hed to HDF format:
	e2proc2d.py test.hed test.hdf		

	apply a 10 A low-pass filter to a stack of particles and write output to a new file. 
	e2proc2d.py ptcl.hdf ptcl.filt.hdf --process=filter.lowpass.gauss:cutoff_freq=0.1

	invert the contrast in a BDB database. Overwrite the original images
	e2proc2d.py bdb:particles#set1 bdb:particles#set1 --inplace --mult=-1

	'e2help.py processors -v 2' for a detailed list of available procesors
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--apix", type=float, help="A/pixel for S scaling")
	parser.add_argument("--average", action="store_true", help="Averages all input images (without alignment) and writes a single output image")
	parser.add_argument("--averager",type=str,help="If --average is specified, this is the averager to use (e2help.py averager). Default=mean",default="mean")
	parser.add_argument("--calcsf", metavar="n outputfile", type=str, nargs=2, help="calculate a radial structure factor for the image and write it to the output file, must specify apix. divide into <n> angular bins")    
	parser.add_argument("--calccont", action="store_true", help="Compute the low resolution azimuthal contrast of each image and put it in the header as eval_contrast_lowres. Larger values imply more 'interesting' images.")
	parser.add_argument("--clip", metavar="xsize,ysize", type=str, action="append", help="Specify the output size in pixels xsize,ysize[,xcenter,ycenter], images can be made larger or smaller.")
	parser.add_argument("--exclude", metavar="exclude-list-file", type=str, help="Excludes image numbers in EXCLUDE file")
	parser.add_argument("--fftavg", metavar="filename", type=str, help="Incoherent Fourier average of all images and write a single power spectrum image")
	parser.add_argument("--process", metavar="processor_name:param1=value1:param2=value2", type=str, action="append", help="apply a processor named 'processorname' with all its parameters/values.")
	parser.add_argument("--mult", metavar="k", type=float, help="Multiply image by a constant. mult=-1 to invert contrast.")
	parser.add_argument("--first", metavar="n", type=int, default=0, help="the first image in the input to process [0 - n-1])")
	parser.add_argument("--last", metavar="n", type=int, default=-1, help="the last image in the input to process")
	parser.add_argument("--list", metavar="listfile", type=str, help="Works only on the image numbers in LIST file")
	parser.add_argument("--select", metavar="selectname", type=str, help="Works only on the images in named selection set from bdb:select")
	parser.add_argument("--inplace", action="store_true", help="Output overwrites input, USE SAME FILENAME, DO NOT 'clip' images.")
	parser.add_argument("--interlv", metavar="interleave-file", type=str, help="Specifies a 2nd input file. Output will be 2 files interleaved.")
	parser.add_argument("--meanshrink", metavar="n", type=int, action="append", help="Reduce an image size by an integral scaling factor using average. Clip is not required.")
	parser.add_argument("--medianshrink", metavar="n", type=int, action="append", help="Reduce an image size by an integral scaling factor, uses median filter. Clip is not required.")
	parser.add_argument("--mraprep",  action="store_true", help="this is an experimental option")
	parser.add_argument("--outmode", type=str, default="float", help="All EMAN2 programs write images with 4-byte floating point values when possible by default. This allows specifying an alternate format when supported (float, int8, int16, int32, uint8, uint16, uint32). Values are rescaled to fill MIN-MAX range.")
	parser.add_argument("--outnorescale", action="store_true", default=False, help="If specified, floating point values will not be rescaled when writing data as integers. Values outside of range are truncated.")
	parser.add_argument("--mrc16bit",  action="store_true", help="(deprecated, use --outmode instead) output as 16 bit MRC file")
	parser.add_argument("--mrc8bit",  action="store_true", help="(deprecated, use --outmode instead) output as 8 bit MRC file")
	parser.add_argument("--fixintscaling", type=str, default=None, help="When writing to an 8 or 16 bit integer format the data must be scaled. 'noscale' will assume the pixel values are already correct, 'sane' will pick a good range, a number will set the range to mean+=sigma*number")
	parser.add_argument("--multfile", type=str, action="append", help="Multiplies the volume by another volume of identical size. This can be used to apply masks, etc.")
	
	parser.add_argument("--norefs", action="store_true", help="Skip any input images which are marked as references (usually used with classes.*)")
	parser.add_argument("--outtype", metavar="image-type", type=str, help="output image format, 'mrc', 'imagic', 'hdf', etc. if specify spidersingle will output single 2D image rather than 2D stack.")
	parser.add_argument("--radon",  action="store_true", help="Do Radon transform")
	parser.add_argument("--randomize", type=str, action="append",help="Randomly rotate/translate the image. Specify: da,dxy,flip  da is a uniform distribution over +-da degrees, dxy is a uniform distribution on x/y, if flip is 1, random handedness changes will occur")
	parser.add_argument("--rotavg", action="store_true", help="Compute the 1-D rotational average of each image as a final step before writing the output", default=False)
	parser.add_argument("--rotate", type=float, action="append", help="Rotate clockwise (in degrees)")
	parser.add_argument("--rfp",  action="store_true", help="this is an experimental option")
	parser.add_argument("--fp",  type=int, help="This generates rotational/translational 'footprints' for each input particle, the number indicates which algorithm to use (0-6)")
	parser.add_argument("--scale", metavar="f", type=float, action="append", help="Scale by specified scaling factor. Clip must also be specified to change the dimensions of the output map.")
	parser.add_argument("--selfcl", metavar="steps mode", type=int, nargs=2, help="Output file will be a 180x180 self-common lines map for each image.")
	parser.add_argument("--setsfpairs",  action="store_true", help="Applies the radial structure factor of the 1st image to the 2nd, the 3rd to the 4th, etc") 
	parser.add_argument("--split", metavar="n", type=int, help="Splits the input file into a set of n output files")
	parser.add_argument("--translate", type=str, action="append", help="Translate by x,y pixels")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, help="verbose level [0-9], higner number means higher level of verboseness",default=1)
	parser.add_argument("--plane", metavar=threedplanes, type=str, default='xy', help="Change the plane of image processing, useful for processing 3D mrcs as 2D images.")
	parser.add_argument("--writejunk", action="store_true", help="Writes the image even if its sigma is 0.", default=False)
	parser.add_argument("--swap", action="store_true", help="Swap the byte order", default=False)
	parser.add_argument("--threed2threed", action="store_true", help="Process 3D image as a stack of 2D slices, then output as a 3D image", default=False)	
	parser.add_argument("--threed2twod", action="store_true", help="Process 3D image as a stack of 2D slices, then output as a 2D stack", default=False)
	parser.add_argument("--twod2threed", action="store_true", help="Process a stack of 2D images, then output as a 3D image.", default=False)
	parser.add_argument("--unstacking", action="store_true", help="Process a stack of 2D images, then output a a series of numbered single image files", default=False)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	parser.add_argument("--step",type=str,default="0,1",help="Specify <init>,<step>. Processes only a subset of the input data. For example, 0,2 would process only the even numbered particles")
	
	# Parallelism

	parser.add_argument("--parallel","-P",type=str,help="Run in parallel, specify type:n=<proc>:option:option",default=None)
	
	append_options = ["clip", "process", "meanshrink", "medianshrink", "scale", "randomize", "rotate","translate","multfile"]

	optionlist = pyemtbx.options.get_optionlist(sys.argv[1:])
	
	(options, args) = parser.parse_args()
	
	if len(args) != 2:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
		sys.exit(1)
	
	logid=E2init(sys.argv,options.ppid)
	
	try : options.step=int(options.step.split(",")[0]),int(options.step.split(",")[1])		# convert strings to tuple
	except:
		print "Invalid --step specification"
		sys.exit(1)
		
	if options.step!=(0,1) and options.first:
		print 'Invalid options. You used --first and --step. The --step option contains both a step size and the first image to step from. Please use only the --step option rather than --step and --first'
		sys.exit(1)
		
	if options.mrc16bit:
		print "Deprecated option mrc16bit, please use outmode=int16"
		options.outmode="int16"

	if options.mrc8bit:
		print "Deprecated option mrc8bit, please use outmode=int8|uint8"
		options.outmode="int8"

	if not file_mode_map.has_key(options.outmode) :
		print "Invalid output mode, please specify one of :\n",str(file_mode_map.keys()).translate(None,'"[]')
		sys.exit(1)

	infile = args[0]
	outfile = args[1]

	if options.parallel : 
		r=doparallel(sys.argv,options.parallel,args)
		E2end(logid)
		sys.exit(r)
		
	if options.average:
		averager=parsemodopt(options.averager)
		average=Averagers.get(averager[0], averager[1])
	else : average = None

	fftavg = None
	
	n0 = options.first
	n1 = options.last
		
	sfout_n = 0
	sfout = None
	sf_amwid = 0
		
	MAXMICROCTF = 1000
	defocus_val = [0] * MAXMICROCTF
	bfactor_val = [0] * MAXMICROCTF
	
	if options.verbose>2:
		Log.logger().set_level(options.verbose-2)
	
	d = EMData()
	threed_xsize = 0
	threed_ysize = 0
	nimg = 1

	if options.threed2threed or options.threed2twod:
		d.read_image(infile, 0, True)
		if(d.get_zsize() == 1):
			print 'Error: need 3D image to use this option'
			return
		else:
			if options.verbose>0:
				print "Process 3D as a stack of %d 2D images" % d.get_zsize()
			nimg = d.get_zsize()
			if n1 > nimg:
				print 'The value for --last is greater than the number of images in the input stack. Exiting'
				n1=options.last
			if options.step[0] > n0:
				n0 = options.step[0]
			threed_xsize = d.get_xsize()
			threed_ysize = d.get_ysize()
			isthreed = False
	else:
		nimg = EMUtil.get_image_count(infile)
		# reads header only
		isthreed = False
		plane = options.plane
		[tomo_nx, tomo_ny, tomo_nz] = gimme_image_dimensions3D(infile)
		if tomo_nz != 1:
			isthreed = True
			if not plane in threedplanes:
				parser.error("the plane (%s) you specified is invalid" %plane)
		
	if not isthreed:
		if nimg <= n1 or n1 < 0:
			n1 = nimg - 1
	else:
		if plane in xyplanes:
			n1 = tomo_nz-1
		elif plane in xzplanes:
			n1 = tomo_ny-1
		elif plane in yzplanes:
			n1 = tomo_nx-1
		if options.last >= 0 and options.last < n1:
			n1 = options.last
		elif options.last > n1:
			print 'The value for --last is greater than the number of images in the input stack. It is being set to the maximum length of the images'
			n1 = tomo_nz-1
		threed = EMData()
		threed.read_image(infile)
		
	ld = EMData()
	if options.step[0] > options.first:
		n0=options.step[0]

	if options.verbose>0:
		print "%d images, processing %d-%d stepping by %d"%(nimg,n0,n1,options.step[1])

	# Now we deal with inclusion/exclusion lists

	if options.list or options.select :
		imagelist=[0]*nimg
		
		if options.list:
			for i in read_number_file(options.list) : imagelist[i]=1
			
		if options.select:
			db=db_open_dict("bdb:.#select",ro=True)
			for i in db[options.select]: imagelist[i]=1
			
	else:
		imagelist=[1]*nimg
		
	if options.exclude :
		for i in read_number_file(options.exclude) : imagelist[i]=0
	
	sfcurve1 = None

	lasttime=time.time()
	outfilename_no_ext = outfile[:-4]
	outfilename_ext = outfile[-3:]
	if outfilename_ext == "rcs" and outfile[-4:] == "mrcs":
		outfilename_ext = outfile[-4:]
	dummy=0										#JESUS

	for i in range(n0, n1+1, options.step[1]):
		if options.verbose >= 1:
			
			if time.time()-lasttime>3 or options.verbose>2 :
				sys.stdout.write(" %7d\r" %i)
				sys.stdout.flush()
				lasttime=time.time()
		if imagelist and (not imagelist[i]):
			continue

		if options.split and options.split > 1:
			outfile = outfilename_no_ext + ".%02d." % (i % options.split) + outfilename_ext

		if not isthreed:
			if options.threed2threed or options.threed2twod:
				d = EMData()
				d.read_image(infile, 0, False, Region(0,0,i,threed_xsize,threed_ysize,1))
			else:
				d = EMData()
				d.read_image(infile, i)
		else:
			if plane in xyplanes:
				roi=Region(0,0,i,tomo_nx,tomo_ny,1)
				d = threed.get_clip(roi)
				#d.read_image(infile,0, HEADER_AND_DATA, roi)
				d.set_size(tomo_nx,tomo_ny,1)
			elif plane in xzplanes:
				roi=Region(0,i,0,tomo_nx,1,tomo_nz)
				d = threed.get_clip(roi)
				#d.read_image(infile,0, HEADER_AND_DATA, roi)
				d.set_size(tomo_nx,tomo_nz,1)
			elif plane in yzplanes:
				roi=Region(i,0,0,1,tomo_ny,tomo_nz)
				d = threed.get_clip(roi)
				#d.read_image(infile,0, HEADER_AND_DATA, roi)
				d.set_size(tomo_ny,tomo_nz,1)

		sigma = d.get_attr("sigma").__float__()
		if sigma == 0:
			if options.threed2threed or options.threed2twod:
				pass
			else:
				if options.verbose>0:
					print "Warning: sigma = 0 for image ",i
				if options.writejunk == False:
					if options.verbose>0:
						print "Use the writejunk option to force writing this image to disk"
					continue

		if not "outtype" in optionlist:
			optionlist.append("outtype")
		
		index_d = {}
		for append_option in append_options:
			index_d[append_option] = 0
		for option1 in optionlist:
			
			nx = d.get_xsize()
			ny = d.get_ysize()
		
			if option1 == "apix":
				apix = options.apix
				d.set_attr('apix_x', apix)
				d.set_attr('apix_y', apix)
				d.set_attr('apix_z', apix)
				try:
					if i==n0 and d["ctf"].apix!=apix :
						if options.verbose>0:
							print "Warning: A/pix value in CTF was %1.2f, changing to %1.2f. May impact CTF parameters."%(d["ctf"].apix,apix)
					d["ctf"].apix=apix
				except: pass
			
			if option1 == "process":
				fi = index_d[option1]
				(processorname, param_dict) = parsemodopt(options.process[fi])
				if not param_dict : param_dict={}
								
				#Parse the options to convert the image file name to EMData object(for both plain image file and bdb file)
				for key in param_dict.keys():
					#print str(param_dict[key])
					if str(param_dict[key]).find('bdb:')!=-1 or not str(param_dict[key]).isdigit():
						try:
							param_dict[key] = EMData(param_dict[key])			
						except:
							pass

				d.process_inplace(processorname, param_dict)
				index_d[option1] += 1

			elif option1 == "mult" :
				d.mult(options.mult)

			elif option1 == "multfile":
				mf=EMData(options.multfile[index_d[option1]],0)
				d.mult(mf)
				mf=None
				index_d[option1] += 1

			elif option1 == "calccont":
				dd=d.process("math.rotationalsubtract")
				f=dd.do_fft()
				#f=d.do_fft()
				if d["apix_x"]<=0 : raise Exception,"Error: 'calccont' requires an A/pix value, which is missing in the input images"
				lopix=int(d["nx"]*d["apix_x"]/200.0)
				hipix=int(d["nx"]*d["apix_x"]/25.0)
				if lopix==hipix : lopix,hipix=3,d["nx"]/5		# in case the A/pix value is drastically out of range
				r=f.calc_radial_dist(d["ny"]/2,0,1.0,1)
				lo=sum(r[lopix:hipix])/(hipix-lopix)
				hi=sum(r[hipix+1:-1])/(len(r)-hipix-2)
				print lopix, hipix, lo, hi
				d["eval_contrast_lowres"]=lo/hi
#				print lopix,hipix,lo,hi,lo/hi
				
			elif option1 == "norefs" and d["ptcl_repr"] <= 0:
				continue
				
			elif option1 == "setsfpairs":
				dataf = d.do_fft()
#				d.gimme_fft()
				x0 = 0
				step = 0.5
				
				if i%2 == 0:
					sfcurve1 = dataf.calc_radial_dist(nx, x0, step)
				else:
					sfcurve2 = dataf.calc_radial_dist(nx, x0, step)
					for j in range(nx):
						if sfcurve1[j] > 0 and sfcurve2[j] > 0:
							sfcurve2[j] = sqrt(sfcurve1[j] / sfcurve2[j])
						else:
							sfcurve2[j] = 0;

						dataf.apply_radial_func(x0, step, sfcurve2);
						d = dataf.do_ift();
#						dataf.gimme_fft();
					
			elif option1 == "rfp":
				d = d.make_rotational_footprint()

			elif option1 == "fp":
				d = d.make_footprint(options.fp)

			elif option1 == "scale":
				scale_f = options.scale[index_d[option1]]
				if scale_f != 1.0:
					d.scale(scale_f)
				index_d[option1] += 1

			elif option1 == "rotate":
				rotatef = options.rotate[index_d[option1]]
				if rotatef!=0.0 : d.rotate(rotatef,0,0)
				index_d[option1] += 1
				
			elif option1 == "translate":
				tdx,tdy=options.translate[index_d[option1]].split(",")
				tdx,tdy=float(tdx),float(tdy)
				if tdx !=0.0 or tdy != 0.0 :
					d.translate(tdx,tdy,0.0)
				index_d[option1] += 1

			elif option1 == "clip":
				ci = index_d[option1]
				clipcx=nx/2
				clipcy=ny/2
				try: clipx,clipy,clipcx,clipcy = options.clip[ci].split(",")
				except: clipx, clipy = options.clip[ci].split(",")
				clipx, clipy = int(clipx),int(clipy)
				clipcx, clipcy = int(clipcx),int(clipcy)
				
				e = d.get_clip(Region(clipcx-clipx/2, clipcy-clipy/2, clipx, clipy))
				try: e.set_attr("avgnimg", d.get_attr("avgnimg"))
				except: pass
				d = e
				index_d[option1] += 1
			
			elif option1 == "randomize" :
				ci = index_d[option1]
				rnd = options.randomize[ci].split(",")
				rnd[0]=float(rnd[0])
				rnd[1]=float(rnd[1])
				rnd[2]=int(rnd[2])
				t=Transform()
				t.set_params({"type":"2d","alpha":random.uniform(-rnd[0],rnd[0]),"mirror":random.randint(0,rnd[2]),"tx":random.uniform(-rnd[1],rnd[1]),"ty":random.uniform(-rnd[1],rnd[1])})
				d.transform(t)

			elif option1 == "medianshrink":
				shrink_f = options.medianshrink[index_d[option1]]
				if shrink_f > 1:
					d.process_inplace("math.medianshrink",{"n":shrink_f})
				index_d[option1] += 1

			elif option1 == "meanshrink":
				mshrink = options.meanshrink[index_d[option1]]
				if mshrink > 1:
					d.process_inplace("math.meanshrink",{"n":mshrink})
				index_d[option1] += 1
		
			elif option1 == "selfcl":
				scl = options.selfcl[0] / 2
				sclmd = options.selfcl[1]
				sc = EMData()
			
				if sclmd == 0:
					sc.common_lines_real(d, d, scl, true)
				else:
					e = d.copy()
					e.process_inplace("xform.phaseorigin")
				
					if sclmd == 1:
						sc.common_lines(e, e, sclmd, scl, true)
						sc.process_inplace("math.linear", Dict("shift", EMObject(-90.0), "scale", EMObject(-1.0)))		
					elif sclmd == 2:
						sc.common_lines(e, e, sclmd, scl, true)
					else:
						if options.verbose>0:
							print "Error: invalid common-line mode '" + sclmd + "'"
						sys.exit(1)
				
			elif option1 == "radon":
				r = d.do_radon()
				d = r
			
			elif option1 == "average":
				average.add_image(d)
				continue

			elif option1 == "fftavg":
				if not fftavg:
					fftavg = EMData()
					fftavg.set_size(nx+2, ny)
					fftavg.set_complex(1)
					fftavg.to_zero()
				d.process_inplace("mask.ringmean")
				d.process_inplace("normalize")
				df = d.do_fft()
				df.mult(df.get_ysize())
				fftavg.add_incoherent(df)
#				d.gimme_fft
				continue

			elif option1 == "calcsf":
				sfout_n = int(options.calcsf[0])
				sfout = options.calcsf[1]
				sf_amwid = 2 * math.pi / sfout_n
					
				dataf = d.do_fft()
#				d.gimme_fft()
				curve = dataf.calc_radial_dist(ny, 0, 0.5)
				outfile2 = sfout
				
				if n1 != 0:
					outfile2 = sfout + ".%03d" % (i+100)
					
				sf_dx = 1.0 / (apix * 2.0 * ny)
				Util.save_data(0, sf_dx, curve, outfile2)
		
				if sfout_n > 0:
					for j in range(0, sfout_n):
						curve = dataf.calc_radial_dist(ny, 0, 0.5, j * sf_amwid, sf_amwid)                    
						outfile2 = os.path.basename(sfout) + "-" + str(j) + "-" + str(sfout_n) + ".pwr"
						if n1 != 0:
							outfile2 = outfile2 + ".%03d" % (i+100)
						Util.save_data(0, sf_dx, curve, outfile2)
				
			elif option1 == "interlv":
				if not options.outtype:
					options.outtype = "unknown"
				d.append_image(outfile, EMUtil.get_image_ext_type(options.outtype))
				d.read_image(options.interlv, i)
				
			elif option1 == "outtype":
				if not options.outtype:
					options.outtype = "unknown"
				if i==0:
					original_outfile = outfile	
				if options.outtype in ["mrc", "pif", "png", "pgm"]:
					if n1 != 0:
						outfile = "%03d." % (i + 100) + original_outfile
				elif options.outtype == "spidersingle":
					if n1 != 0:
						if i==0:
							if outfile.find('.spi')>0:
								nameprefix = outfile[:-4]
							else:
								nameprefix = outfile
						spiderformat = "%s%%0%dd.spi" % (nameprefix, int(math.log10(n1+1-n0))+1)						
						outfile = spiderformat % i

				#if options.inplace:
						#d.write_image(outfile, i)
				#elif options.mraprep:
						#outfile = outfile + "%04d" % i + ".lst"
						#options.outtype = "lst"
				
				if options.fixintscaling!=None :
					if options.fixintscaling=="sane" : 
						d["render_min"]=d["mean"]-d["sigma"]*2.5
						d["render_max"]=d["mean"]+d["sigma"]*2.5
					elif options.fixintscaling=="noscale":
						d["render_min"]=0.0
						if "mrc16bit" in optionlist : d["render_max"]=65535.0
						else : d["render_max"]=255.0
					else:
						try: 
							sca=int(options.fixintscaling)
							d["render_min"]=d["mean"]-d["sigma"]*sca
							d["render_max"]=d["mean"]+d["sigma"]*sca
						except:
							print "Warning: bad fixintscaling option"

				#print_iminfo(data, "Final")
				if options.outmode!="float":
#					if outfile[-4:]!=".hdf" :
#						print "WARNING: outmode is not working correctly for non HDF images in 2.1beta3. We expect to have this fixed in the next few days."
					if options.outnorescale :
						# This sets the minimum and maximum values to the range for the specified type, which should result in no rescaling
						outmode=file_mode_map[options.outmode]
						d["render_min"]=file_mode_range[outmode][0]
						d["render_max"]=file_mode_range[outmode][1]
					else:
						d["render_min"]=d["minimum"]
						d["render_max"]=d["maximum"]

				if not options.average:	#skip writing the input image to output file 
					#write processed image to file
					if options.threed2threed or options.twod2threed:    #output a single 3D image
						#shift = 0
						
						if dummy==0:										#The "dummy" volume, as termed by Grant, only needs to be created once
																			#This dummy variable changes to dummy=1 after that happens.
							#print "\n\n\n\nTOMO OPTION IS ON!"
							#print "And this is i", i
							#print "\n\n\n\nI will create DUMMY!!!"
							z = nimg
							#print "Z should be equal to nimg, lets see: nimg,z",nimg,z
							if options.list:
								f=open(options.list,'r')
								lines = f.read().split(',')
								#print "Lines are", lines
								f.close()
								z = len( lines )
								#shift = nimg - z
								#print "If --list, z should be the length of lines in list; lines,z", len(lines),z
							elif options.exclude:
								f=open(options.exclude,'r')
								lines=f.read().split(',')
								#print "lines are", lines
								f.close()
								z = nimg - len(lines)
								#shift = len(lines)
								#print "If --exclude, z should be z size of input minus lines in exclude; nimg, len lines,zout", nimg, len(lines), z 
								
							out3d_img = EMData(d.get_xsize(), d.get_ysize(), z)
							if 'mrc8bit' in optionlist:
								#print "Writing dummy mrc8bit"
								out3d_img.write_image(outfile, 0, EMUtil.get_image_ext_type(options.outtype), False, None, EMUtil.EMDataType.EM_UCHAR, not(options.swap))
								dummy=1		
								#print "Wrote dummy mrc8bit"			
							elif 'mrc16bit' in optionlist:
								#print "Writting dummy mrc16bit"
								out3d_img.write_image(outfile, 0, EMUtil.get_image_ext_type(options.outtype), False, None, EMUtil.EMDataType.EM_SHORT, not(options.swap))
								dummy=1
								#print "Wrote dummy mrc16bit"
							else:
								#print "Writting dummy float"
#								out3d_img.write_image(outfile, 0, EMUtil.get_image_ext_type(options.outtype), False, None, EMUtil.EMDataType.EM_FLOAT, not(options.swap))
								out3d_img.write_image(outfile, 0, EMUtil.get_image_ext_type(options.outtype), False, None, file_mode_map[options.outmode], not(options.swap))
								dummy=1
								#print "Wrote dummy float"
							
						#print "imagelist is", imagelist
						if options.list or options.exclude:
							if imagelist[i] != 0:
								#if options.twod2threed:
								#	region = Region(0, 0, imagelist[0:i].count(1), d.get_xsize(), d.get_ysize(), 1)
								#elif options.threed2threed:
								region = Region(0, 0, imagelist[0:i].count(1), d.get_xsize(), d.get_ysize(), 1)
						else:
							#if options.twod2threed:
							#	region = Region(0, 0, i, d.get_xsize(), d.get_ysize(), 1)
							#elif options.threed2threed:
							region = Region(0, 0, i, d.get_xsize(), d.get_ysize(), 1)
						#print "outtype is", options.outtype
						if 'mrc8bit' in optionlist:
								d.write_image(outfile, 0, EMUtil.get_image_ext_type(options.outtype), False, region, EMUtil.EMDataType.EM_UCHAR, not(options.swap))
						elif 'mrc16bit' in optionlist:
								d.write_image(outfile, 0, EMUtil.get_image_ext_type(options.outtype), False, region, EMUtil.EMDataType.EM_SHORT, not(options.swap))
						else:
#								d.write_image(outfile, 0, EMUtil.get_image_ext_type(options.outtype), False, region, EMUtil.EMDataType.EM_FLOAT, not(options.swap))
								d.write_image(outfile, 0, EMUtil.get_image_ext_type(options.outtype), False, region, file_mode_map[options.outmode], not(options.swap))
					
					elif options.unstacking:	#output a series numbered single image files
						if 'mrc8bit' in optionlist:
							d.write_image(outfile.split('.')[0]+'-'+str(i+1).zfill(len(str(nimg)))+'.mrc', -1, EMUtil.ImageType.IMAGE_MRC, False, None, EMUtil.EMDataType.EM_UCHAR, not(options.swap))
						elif 'mrc16bit' in optionlist:
							d.write_image(outfile.split('.')[0]+'-'+str(i+1).zfill(len(str(nimg)))+'.mrc', -1, EMUtil.ImageType.IMAGE_MRC, False, None, EMUtil.EMDataType.EM_SHORT, not(options.swap))
						else:
							#d.write_image(outfile.split('.')[0]+'-'+str(i+1).zfill(len(str(nimg)))+'.'+outfile.split('.')[-1])
#							d.write_image(outfile.split('.')[0]+'-'+str(i+1).zfill(len(str(nimg)))+'.'+outfile.split('.')[-1], 0, EMUtil.get_image_ext_type(options.outtype), False, None, EMUtil.EMDataType.EM_FLOAT, not(options.swap))
							d.write_image(outfile.split('.')[0]+'-'+str(i+1).zfill(len(str(nimg)))+'.'+outfile.split('.')[-1], 0, EMUtil.get_image_ext_type(options.outtype), False, None, file_mode_map[options.outmode], not(options.swap))
							#print "I will unstack to HDF" #JESUS
							
					else:   #output a single 2D image or a 2D stack			
						# optionally replace the output image with its rotational average
						if options.rotavg:
							rd=d.calc_radial_dist(d["nx"],0,0.5,0)
							d=EMData(len(rd),1,1)
							for x in xrange(len(rd)): d[x]=rd[x]
						
						if 'mrc8bit' in optionlist:
							d.write_image(outfile.split('.')[0]+'.mrc', -1, EMUtil.ImageType.IMAGE_MRC, False, None, EMUtil.EMDataType.EM_UCHAR, not(options.swap))
						elif 'mrc16bit' in optionlist:
							d.write_image(outfile.split('.')[0]+'.mrc', -1, EMUtil.ImageType.IMAGE_MRC, False, None, EMUtil.EMDataType.EM_SHORT, not(options.swap))
						else:
							if options.inplace:
#								d.write_image(outfile, i, EMUtil.get_image_ext_type(options.outtype), False, None, EMUtil.EMDataType.EM_FLOAT, not(options.swap))
								d.write_image(outfile, i, EMUtil.get_image_ext_type(options.outtype), False, None, file_mode_map[options.outmode], not(options.swap))
							else: 
#								d.write_image(outfile, -1, EMUtil.get_image_ext_type(options.outtype), False, None, EMUtil.EMDataType.EM_FLOAT, not(options.swap))
								d.write_image(outfile, -1, EMUtil.get_image_ext_type(options.outtype), False, None, file_mode_map[options.outmode], not(options.swap))
				
	#end of image loop
		
	if average:
		avg=average.finish()
		avg["ptcl_repr"]=(n1-n0+1)
#		avg.mult(1.0/(n1-n0+1.0))
#		average.process_inplace("normalize");
		if options.inplace: avg.write_image(outfile,0)
		else : avg.write_image(outfile,-1)

	if options.fftavg:
		fftavg.mult(1.0 / sqrt(n1 - n0 + 1))
		fftavg.write_image(options.fftavg, 0)

		curve = fftavg.calc_radial_dist(ny, 0, 0.5,1)
		outfile2 = options.fftavg+".txt"
		
		sf_dx = 1.0 / (apix * 2.0 * ny)
		Util.save_data(0, sf_dx, curve, outfile2)
	
	try:
		n_outimg = EMUtil.get_image_count(outfile)
		if options.verbose>0:
			print str(n_outimg) + " images"
	except:
		pass	
	
	E2end(logid)

def doparallel(argv,parallel,args):
	pass

if __name__ == "__main__":
	main()
