#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/05/2011 (sludtke@bcm.edu)
# Copyright (c) 2011 Baylor College of Medicine
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

# $Id: e2proc2dmulti.py

from EMAN2 import *
from optparse import OptionParser
import sys
import os.path
import math
import random
import pyemtbx.options
import time

def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + """ [options] <inputfiles>

	A generic 2D image processing program which acts on multiple input stacks. Default operation is similar to e2proc2d.py --inplace, but operating
	on many files.
"""

	parser = OptionParser(usage,version=EMANVERSION)

	
	parser.add_option("--postfix", type="string", help="Adds this string to each input filename to produce output filename (avoid overwriting)",Default=None)
	parser.add_option("--allparticles",action="store_true",help="Will process all particle sets stored in BDB in the particles subdirectory",default=False)
	parser.add_option("--apix", type="float", help="A/pixel for S scaling")
	parser.add_option("--clip", metavar="xsize,ysize", type="string", action="append",
					help="Specify the output size in pixels xsize,ysize[,xcenter,ycenter], images can be made larger or smaller.")
	parser.add_option("--fftavg", metavar="filename", type="string",
					help="Incoherent Fourier average of all images and write a single power spectrum image")
	parser.add_option("--process", metavar="processor_name:param1=value1:param2=value2", type="string",
					action="append", help="apply a processor named 'processorname' with all its parameters/values.")
	parser.add_option("--autoinvert", action="store_true",help="Automatically decides whether to invert each stack of images to make particles white (EMAN2 convention). Decision is made for an entire stack.",default=False)
	parser.add_option("--mult", metavar="k", type="float", help="Multiply image by a constant. mult=-1 to invert contrast.")
	parser.add_option("--meanshrink", metavar="n", type="int", action="append",
					help="Reduce an image size by an integral scaling factor using average. Clip is not required.")
	parser.add_option("--medianshrink", metavar="n", type="int", action="append",
					help="Reduce an image size by an integral scaling factor, uses median filter. Clip is not required.")
	parser.add_option("--multfile", type="string", action="append",
								help="Multiplies the volume by another volume of identical size. This can be used to apply masks, etc.")
	
	parser.add_option("--norefs", action="store_true", help="Skip any input images which are marked as references (usually used with classes.*)")
	parser.add_option("--radon",  action="store_true", help="Do Radon transform")
	parser.add_option("--randomize", type="string", action="append",help="Randomly rotate/translate the image. Specify: da,dxy,flip  da is a uniform distribution over +-da degrees, dxy is a uniform distribution on x/y, if flip is 1, random handedness changes will occur")
	parser.add_option("--rotate", type="float", action="append", help="Rotate clockwise (in degrees)")
	parser.add_option("--fp",  type="int", help="This generates rotational/translational 'footprints' for each input particle, the number indicates which algorithm to use (0-6)")
	parser.add_option("--scale", metavar="f", type="float", action="append",
					help="Scale by specified scaling factor. Clip must also be specified to change the dimensions of the output map.")
	parser.add_option("--selfcl", metavar="steps mode", type="int", nargs=2,
					help="Output file will be a 180x180 self-common lines map for each image.")
	parser.add_option("--translate", type="string", action="append", help="Translate by x,y pixels")
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n", type="int", help="verbose level [0-9], higner number means higher level of verboseness",default=1)
	parser.add_option("--writejunk", action="store_true", help="Writes the image even if its sigma is 0.", default=False)
	
	append_options = ["clip", "process", "meanshrink", "medianshrink", "scale", "randomize", "rotate", "translate", "multfile"]

	optionlist = pyemtbx.options.get_optionlist(sys.argv[1:])
	
	(options, args) = parser.parse_args()
	
	if options.allparticles:
		args=["bdb:particles#"+i for i in db_list_dicts("bdb:particles")]
		args.sort()
		if options.verbose : print "%d particle files identified"%len(args)

	if len(args)<1 : 
		print "Please specify at least one input file to operate on"
		sys.exit(1)
	
	logid=E2init(sys.argv)

	for infile in args:
		if options.postfix!=None :
			if infile[:4].lower()=="bdb:" : outfile=infile+postfix
			elif "." in infile : outfile=infile.rsplit(".",1)[0]+postfix+"."+infile.rsplit(".",1)[1]
			else : outfile=infile+postfix
		else : outfile=infile
		
		if options.verbose>2:
			Log.logger().set_level(options.verbose-2)
		
		d = EMData()
		nimg = EMUtil.get_image_count(infile)

		ld = EMData()
		print "%s : processing %d images"%(infile,nimg)

		lasttime=time.time()
		for i in xrange(nimg):
			if options.verbose >= 1:
				
				if time.time()-lasttime>3 or options.verbose>2 :
					sys.stdout.write(" %7d\r" %i)
					sys.stdout.flush()
					lasttime=time.time()

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
					print "Warning: sigma = 0 for image ",i
					if options.writejunk == False:
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
							print "Warning: A/pix value in CTF was %1.2f, changing to %1.2f. May impact CTF parameters."%(d["ctf"].apix,apix)
						d["ctf"].apix=apix
					except: pass
				
				if option1 == "process":
					fi = index_d[option1]
					(processorname, param_dict) = parsemodopt(options.process[fi])
					if not param_dict : param_dict={}
					d.process_inplace(processorname, param_dict)
					index_d[option1] += 1

				elif option1 == "mult" :
					d.mult(options.mult)

				elif option1 == "multfile":
					mf=EMData(options.multfile[index_d[option1]],0)
					d.mult(mf)
					mf=None
					index_d[option1] += 1
					
				elif option1 == "norefs" and d.get_average_nimg() <= 0:
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
							print "Error: invalid common-line mode '" + sclmd + "'"
							sys.exit(1)
					
				elif option1 == "radon":
					r = d.do_radon()
					d = r
				
				elif option1 == "average":
					if not average:
						average = d.copy()
					else:
						average.add(d)
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
					d.read_image(options.interlv, i)
					d.append_image(outfile, IMAGIC)
				
				elif option1 == "outtype":
					if not options.outtype:
						options.outtype = "unknown"
						
					if options.outtype in ["mrc", "pif", "png", "pgm"]:
						if n1 != 0:
							outfile = "%03d." % (i + 100) + outfile
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
					
					#write processed image to file
					if options.threed2threed or options.twod2threed:    #output a single 3D image
						if i==0:
							out3d_img = EMData(d.get_xsize(), d.get_ysize(), nimg)		
							
						out3d_img.insert_clip(d, (0,0,i))
						
						if i==n1:
							if 'mrc8bit' in optionlist:
								out3d_img.write_image(outfile.split('.')[0]+'.mrc', 0, EMUtil.ImageType.IMAGE_MRC, False, None, EMUtil.EMDataType.EM_UCHAR, not(options.swap))
							elif 'mrc16bit' in optionlist:
								out3d_img.write_image(outfile.split('.')[0]+'.mrc', 0, EMUtil.ImageType.IMAGE_MRC, False, None, EMUtil.EMDataType.EM_SHORT, not(options.swap))
							else:
								out3d_img.write_image(outfile)
					else:   #output a single 2D image or a 2D stack			
						if 'mrc8bit' in optionlist:
							d.write_image(outfile.split('.')[0]+'.mrc', -1, EMUtil.ImageType.IMAGE_MRC, False, None, EMUtil.EMDataType.EM_UCHAR, not(options.swap))
						elif 'mrc16bit' in optionlist:
							d.write_image(outfile.split('.')[0]+'.mrc', -1, EMUtil.ImageType.IMAGE_MRC, False, None, EMUtil.EMDataType.EM_SHORT, not(options.swap))
						else:
							if options.inplace:
								d.write_image(outfile, i, EMUtil.get_image_ext_type(options.outtype), False, None, EMUtil.EMDataType.EM_FLOAT, not(options.swap))
							else: 
								d.write_image(outfile, -1, EMUtil.get_image_ext_type(options.outtype), False, None, EMUtil.EMDataType.EM_FLOAT, not(options.swap))
					
			
	E2end(logid)

def doparallel(argv,parallel,args):
	pass

if __name__ == "__main__":
	main()
