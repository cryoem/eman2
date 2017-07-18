#!/usr/bin/env python

#
# Author: Steven Ludtke, 07/18/2017 (sludtke@bcm.edu)
# Copyright (c) 2017- Baylor College of Medicine
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
import os
import datetime
import time
import traceback
from collections import Counter

# usage: e2proc2d.py [options] input ... input output

def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + """ [options] <inputfile> <outputfile>

	Limited version of e2proc2d.py which can operate in parallel using threads. Works on a single input/output file, lacks wildcards and many other options.
	
	Note that --inplace is implied. Output image number is always the same as input image number!
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--outtype", metavar="image-type", type=str, help="output image format, 'mrc', 'imagic', 'hdf', etc. if specify spidersingle will output single 2D image rather than 2D stack.")
	parser.add_argument("--outmode", type=str, default="float", help="All EMAN2 programs write images with 4-byte floating point values when possible by default. This allows specifying an alternate format when supported (float, int8, int16, int32, uint8, uint16, uint32). Values are rescaled to fill MIN-MAX range.")
	parser.add_argument("--fixintscaling", type=str, default=None, help="When writing to an 8 or 16 bit integer format the data must be scaled. 'noscale' will assume the pixel values are already correct, 'sane' will pick a good range, a number will set the range to mean+=sigma*number")

	parser.add_argument("--apix", type=float, help="A/pixel for S scaling")
	parser.add_argument("--clip", metavar="xsize,ysize", type=str, action="append", help="Specify the output size in pixels xsize,ysize[,xcenter,ycenter], images can be made larger or smaller.")
	parser.add_argument("--process", metavar="processor_name:param1=value1:param2=value2", type=str, action="append", help="apply a processor named 'processorname' with all its parameters/values.")
	parser.add_argument("--mult", metavar="k", type=float, help="Multiply image by a constant. mult=-1 to invert contrast.")
	parser.add_argument("--add", metavar="f", type=float,action="append",help="Adds a constant 'f' to the densities")
	parser.add_argument("--meanshrink", metavar="n", type=float, action="append", help="Reduce an image size by an integral (1.5 also allowed) scaling factor using average. eg - 2 will reduce image size to 1/2. Clip is not required.")
	parser.add_argument("--medianshrink", metavar="n", type=int, action="append", help="Reduce an image size by an integral scaling factor, uses median filter. eg - 2 will reduce image size to 1/2. Clip is not required.")
	parser.add_argument("--fouriershrink", metavar="n", type=float, action="append", help="Reduce an image size by an arbitrary scaling factor by clipping in Fourier space. eg - 2 will reduce image size to 1/2.")
	parser.add_argument("--multfile", type=str, action="append", help="Multiplies the image by another image of identical size. This can be used to apply masks, etc.")
	parser.add_argument("--randomize", type=str, action="append",help="Randomly rotate/translate the image. Specify: da,dxy,flip  da is a uniform distribution over +-da degrees, dxy is a uniform distribution on x/y, if flip is 1, random handedness changes will occur")
	parser.add_argument("--rotate", type=float, action="append", help="Rotate clockwise (in degrees)")
	parser.add_argument("--rfp",  action="store_true", help="this is an experimental option")
	parser.add_argument("--fp",  type=int, help="This generates rotational/translational 'footprints' for each input particle, the number indicates which algorithm to use (0-6)")
	parser.add_argument("--scale", metavar="f", type=float, action="append", help="Scale by specified scaling factor. Clip must also be specified to change the dimensions of the output map.")
	parser.add_argument("--anisotropic", type=str,action="append", help="Anisotropic scaling, stretches on one axis and compresses the orthogonal axis. Specify amount,angle. See e2evalrefine")
	parser.add_argument("--selfcl", metavar="steps mode", type=int, nargs=2, help="Output file will be a 180x180 self-common lines map for each image.")
	parser.add_argument("--translate", type=str, action="append", help="Translate by x,y pixels")

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, help="verbose level [0-9], higner number means higher level of verboseness",default=1)
	parser.add_argument("--parallel","-P",type=str,help="Run in parallel, only thread:n supported",default=None)
	parser.add_argument("--threads", default=4,type=int,help="Number of threads to run in parallel on a single computer when multi-computer parallelism isn't useful", guitype='intbox', row=30, col=2, rowspan=1, colspan=1, mode="refinement[4]")


	optionlist = pyemtbx.options.get_optionlist(sys.argv[1:])

	(options, args) = parser.parse_args()

	if len(args) != 2:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
		sys.exit(1)

	if options.parallel!=None:
		if options.parallel[:7]=="thread:" : options.threads=int(options.parallel[7:])

	if options.threads<1 : options.threads=1
	options.threads+=1		# one thread does the I/O, so we add an extra

	logid = E2init(sys.argv,options.ppid)

	if not file_mode_map.has_key(options.outmode) :
		print "Invalid output mode, please specify one of :\n",str(file_mode_map.keys()).translate(None,'"[]')
		sys.exit(1)
	
	N=EMUtil.get_image_count(args[0])
	npt=max(min(100,N/options.threads),1)
	
	jsd=Queue.Queue(0)
	# these start as arguments, but get replaced with actual threads
	thrds=[(jsd,args,options,i,i*npt,min(i+npt,N)) for i in xrange(N/npt+1)]
	
	thrtolaunch=0
	while thrtolaunch<len(thrds) or threading.active.count>1:
		if thrtolaunch<len(thrds):
			while (threading.active_count()==options.threads) : time.sleep(0.1)
			if options.verbose>1 : 
				print "\r Starting thread {}/{}      ".format(thrtolaunch,len(thrds)),
				sys.stdout.flush()
			thrds[thrtolaunch]=threading.Thread(target=procfn,args=thrds[thrtolaunch])		# replace args
			thrds[thrtolaunch].start()
			thrtolaunch+=1
		else: time.sleep(0.1)
		
		# return is [N,dict] a dict of image# keyed processed images
		while not jsd.empty():
			rd=jsd.get()
			for k in rd[1].keys():
				writeimage(rd[1][k],args[1],k,options)
			
			thrds[rd[0]].join()
			thrds[rd[0]]=None

	logid = E2end(sys.argv,options.ppid)
	
	
def procfn(jsd,args,options,thrn,n0,n1):
	optionlist = pyemtbx.options.get_optionlist(sys.argv[1:])

	ret=[thrn,{}]
	for n in xrange(n0, n1):
		d=EMData(args[0],n)

		index_d = Counter()

		for option1 in optionlist:
			if options.verbose > 1 :
				print "option in option list =", option1
			nx = d.get_xsize()
			ny = d.get_ysize()

			if option1 == "apix":
				apix = options.apix
				d.set_attr('apix_x', apix)
				d.set_attr('apix_y', apix)
				d.set_attr('apix_z', apix)

				try:
					if i == n0 and d["ctf"].apix != apix :
						if options.verbose > 0:
							print "Warning: A/pix value in CTF was %1.2f, changing to %1.2f. May impact CTF parameters."%(d["ctf"].apix,apix)

					d["ctf"].apix = apix
				except: pass

			if option1 == "process":
				fi = index_d[option1]
				(processorname, param_dict) = parsemodopt(options.process[fi])

				if not param_dict : param_dict = {}

				# Parse the options to convert the image file name to EMData object
				# (for both plain image file and bdb file)

				for key in param_dict.keys():
					#print str(param_dict[key])

					if str(param_dict[key]).find('bdb:') != -1 or not str(param_dict[key]).isdigit():
						try:
							param_dict[key] = EMData(param_dict[key])			
						except:
							pass

				# these processors don't work _inplace
				if processorname in ["math.bispectrum.slice"]:
					d=d.process(processorname, param_dict)
				else: d.process_inplace(processorname, param_dict)
				index_d[option1] += 1

			elif option1 == "add":
				d.add(options.add[index_d[option1]])
				af=None
				index_d[option1] += 1

			elif option1 == "mult" :
				d.mult(options.mult)

			elif option1 == "multfile":
				mf = EMData(options.multfile[index_d[option1]],0)
				d.mult(mf)
				mf = None
				index_d[option1] += 1

			elif option1 == "rfp":
				d = d.make_rotational_footprint()

			elif option1 == "fp":
				d = d.make_footprint(options.fp)

			elif option1 == "anisotropic":
				try: 
					amount,angle = (options.anisotropic[index_d[option1]]).split(",")
					amount=float(amount)
					angle=float(angle)
				except:
					traceback.print_exc()
					print options.anisotropic[index_d[option1]]
					print "Error: --anisotropic specify amount,angle"
					sys.exit(1)
					
				rt=Transform({"type":"2d","alpha":angle})
				xf=rt*Transform([amount,0,0,0,0,1/amount,0,0,0,0,1,0])*rt.inverse()
				d.transform(xf)

				index_d[option1] += 1


			elif option1 == "scale":
				scale_f = options.scale[index_d[option1]]

				if scale_f != 1.0:
					d.scale(scale_f)

				index_d[option1] += 1

			elif option1 == "rotate":
				rotatef = options.rotate[index_d[option1]]

				if rotatef != 0.0 : d.rotate(rotatef,0,0)

				index_d[option1] += 1

			elif option1 == "translate":
				tdx,tdy = options.translate[index_d[option1]].split(",")
				tdx,tdy = float(tdx),float(tdy)

				if tdx != 0.0 or tdy != 0.0 :
					d.translate(tdx,tdy,0.0)

				index_d[option1] += 1

			elif option1 == "clip":
				ci = index_d[option1]
				clipcx = nx/2
				clipcy = ny/2

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
				rnd[0] = float(rnd[0])
				rnd[1] = float(rnd[1])
				rnd[2] = int(rnd[2])

				t = Transform()
				t.set_params({"type":"2d", "alpha":random.uniform(-rnd[0],rnd[0]), \
								"mirror":random.randint(0,rnd[2]), "tx":random.uniform(-rnd[1],rnd[1]), \
								"ty":random.uniform(-rnd[1],rnd[1])})
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

			elif option1 == "fouriershrink":
				fshrink = options.fouriershrink[index_d[option1]]

				if fshrink > 1:
					d.process_inplace("math.fft.resample",{"n":fshrink})

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
						if options.verbose > 0:
							print "Error: invalid common-line mode '" + sclmd + "'"

						sys.exit(1)
		ret[1][n]=d

	jsd.put(ret)

def writeimage(d,fsp,n,options):

	if not options.outtype:
		options.outtype = "unknown"

	if options.outtype in ["mrc", "pif", "png", "pgm"]:
		if n != 0:
			fsp = "{:03d}.{}".format(n,fsp)
			n=0

	if options.outmode!="float":
		if options.fixintscaling!=None:
			if options.fixintscaling=="noscale":
				if   options.outmode == "int8" :
					u =   -128.0
					v =    127.0
				elif options.outmode == "uint8" :
					u =      0.0
					v =    255.0
				elif options.outmode == "int16" :
					u = -32768.0
					v =  32767.0
				elif options.outmode == "uint16" :
					u =      0.0
					v =  65535.0
					
			elif options.fixintscaling=="sane":
				u=d["mean"]-d["sigma"]*2.5
				v=d["mean"]-d["sigma"]*2.5
				
			else:
				u=d["mean"]-d["sigma"]*float(options.fixintscaling)
				v=d["mean"]-d["sigma"]*float(options.fixintscaling)
		else:
			u=d["minimum"]
			v=d["maximum"]
				

		if u>=v : v=u+.00001			# shouldn't happen
		d["render_min"] = u
		d["render_max"] = v

	out_type = EMUtil.get_image_ext_type(options.outtype)
	out_mode = file_mode_map[options.outmode]

	d.write_image(fsp, n, out_type, False, None, out_mode)


if __name__ == "__main__":
	main()
