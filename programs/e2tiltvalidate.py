#!/usr/bin/env python
#
# Author: John Flanagan (jfflanag@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine


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

from EMAN2 import *
from optparse import OptionParser
import os, math

def main():
	"""Program to validate a reconstruction by the Richard Henderson tilt validation method. A volume to validate, a small stack (~100 imgs) of untilted and ~10-15 degree
	tilted particles must be presented. The untilted and tilted particle stack must have a one-to-one relationship. In the contour plot, the Tiltaxis is along positive 'Y'
	The tiltaxis angle can be determined from e2RCTboxer.py uisng PairPicker mode. For example, if the tiltaxis is 45 degrees and the tilt angle is -15 degrees, there should
	be a peak in the -X, -Y quadrant at 225 degrees at a magnitude of 15.
	For more details see:
	Optiomal Determination of Particle Orientation, Absolute Hand, and COntrast Loss in Single-particle Electron Cryomicroscopy. Rosenthal, P.B., and Henderson, R. JMB, 333 (2003) pg 721-745
	"""
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options]"""
	parser = OptionParser(usage=usage,version=EMANVERSION)
	
	# options associated with e2tiltvalidate.py
	parser.add_option("--path", type="string",help="The folder the results are placed", default="TiltValidate")
	parser.add_option("--volume", type="string",help="3D volume to validate",default=None)
	parser.add_option("--untiltdata", type="string",help="Stack of untilted images",default=None)
	parser.add_option("--tiltdata", type="string",help="Stack of tilted images",default=None)
	parser.add_option("--align", type="string",help="The name of a aligner to be used in comparing the aligned images",default="translational")
	parser.add_option("--cmp", type="string",help="The name of a 'cmp' to be used in comparing the aligned images",default="ccc")
	parser.add_option("--tiltrange", type="int",help="The anglular tiltranger to search",default=15)
	parser.add_option("--sym",  type="string",help="The recon symmetry", default="c1")
	parser.add_option("--tiltangle", type="float", help="The stage tiltused during data collection", default=None)
	# options associated with e2projector3d.py
	parser.add_option("--delta", type="float",help="The angular step size for alingment", default=20.0)
	# options associated with e2simmx.py
	parser.add_option("--simalign",type="string",help="The name of an 'aligner' to use prior to comparing the images (default=rotate_translate)", default="rotate_translate")
	parser.add_option("--simaligncmp",type="string",help="Name of the aligner along with its construction arguments (default=ccc)",default="ccc")
	parser.add_option("--simcmp",type="string",help="The name of a 'cmp' to be used in comparing the aligned images (default=ccc)", default="ccc")
	parser.add_option("--simralign",type="string",help="The name and parameters of the second stage aligner which refines the results of the first alignment", default=None)
	parser.add_option("--simraligncmp",type="string",help="The name and parameters of the comparitor used by the second stage aligner. (default=dot).",default="dot")
	parser.add_option("--shrink", dest="shrink", type = "int", default=None, help="Optionally shrink the input particles by an integer amount prior to computing similarity scores. For speed purposes.")
	
	parser.add_option("--parallel",type="string",help="Parallelism string",default=None)
	parser.add_option("--verbose", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	global options
	(options, args) = parser.parse_args()

	if not options.volume:
		print "Error a volume to validate must be presented"
		exit(1)
	if not options.tiltdata:
		print "Error a stack of tilted images must be presented"
		#exit(1)
	if not options.untiltdata:
		print "Error a stack of untiled images must be presented"
		exit(1)
	if not options.tiltangle:
		print "Error a tiltangle must be entered"
		exit(1)
		
	options.cmp=parsemodopt(options.cmp)
	options.align=parsemodopt(options.align)
	
	#Read in the images
	tiltimgs = EMData.read_images(options.tiltdata)
	untiltimgs = EMData.read_images(options.untiltdata)
	if len(tiltimgs) != len(untiltimgs):
		print "The untilted image stack is not the same lenght as the tilted stack!!!"
		exit(1)

	# Make a new dir for each run
	dirindex = 1
	while os.path.exists("./%s_%d"%(options.path,dirindex)):
		dirindex += 1
	global workingdir
	workingdir = "%s_%d"%(options.path,dirindex)
	os.mkdir(workingdir)
	
	# Do projections
	e2projectcmd = "e2project3d.py %s --orientgen=eman:delta=%f:inc_mirror=1:perturb=0 --outfile=bdb:%s#projections --projector=standard --sym=%s" % (options.volume,options.delta,workingdir, options.sym)
	if options.parallel: e2projectcmd += " --parallel=%s" %options.parallel
	run(e2projectcmd)
	
	# Make simmx
	e2simmxcmd = "e2simmx.py bdb:%s#projections %s bdb:%s#simmx -f --saveali --cmp=%s --align=%s --aligncmp=%s --verbose=%d" % (workingdir,options.untiltdata,workingdir,options.simcmp,options.simalign,options.simaligncmp,options.verbose)
	if options.simralign: e2simmxcmd += " --ralign=%s --raligncmp=%s" %(options.simralign,options.simraligncmp)
	if options.parallel: e2simmxcmd += " --parallel=%s" %options.parallel
	if options.shrink: e2simmxcmd += " --shrink=%d" %options.shrink
	run(e2simmxcmd)
	
	e2simmxcmd = "e2simmx.py bdb:%s#projections %s bdb:%s#simmx_tilt -f --saveali --cmp=%s --align=%s --aligncmp=%s --verbose=%d" % (workingdir,options.tiltdata,workingdir,options.simcmp,options.simalign,options.simaligncmp,options.verbose)
	if options.simralign: e2simmxcmd += " --ralign=%s --raligncmp=%s" %(options.simralign,options.simraligncmp)
	if options.parallel: e2simmxcmd += " --parallel=%s" %options.parallel
	if options.shrink: e2simmxcmd += " --shrink=%d" %options.shrink
	run(e2simmxcmd)

	# Read in the data
	simmx= EMData.read_images("bdb:%s#simmx"%workingdir)
	simmx_tilt= EMData.read_images("bdb:%s#simmx_tilt"%workingdir)
	projections = EMData.read_images("bdb:%s#projections"%workingdir)
	volume = EMData() 
	volume.read_image(options.volume) # I don't knwo why I cant EMData.read_image.......
	symmeties = Symmetries.get(options.sym)
	
	# Find the differnces in alignment pars, this is an attempt to do per image validation
	for imgnum in xrange(simmx[0].get_ysize()):
		untiltbestscore = float('inf')
		tiltbestscore = float('inf')
		untiltbestrefnum = 0
		tiltbestrefnum = 0
		tiltimgnum = imgnum
		for refnum in xrange(simmx[0].get_xsize()):
			if simmx[0].get_value_at(refnum, imgnum) < untiltbestscore:
				untiltbestscore = simmx[0].get_value_at(refnum, imgnum)
				untiltbestrefnum = refnum		
			if simmx_tilt[0].get_value_at(refnum, tiltimgnum) < tiltbestscore:
				tiltbestscore = simmx_tilt[0].get_value_at(refnum, tiltimgnum)
				tiltbestrefnum = refnum
		# Untilt
		untilt_euler_xform = projections[untiltbestrefnum].get_attr('xform.projection')
		untiltrot = untilt_euler_xform.get_rotation("eman")
		untilt_euler_xform.set_rotation({"type":"eman","az":untiltrot["az"],"alt":untiltrot["alt"],"phi":-simmx[3].get_value_at(untiltbestrefnum, imgnum)})
		# Tilt
		tilt_euler_xform = projections[tiltbestrefnum].get_attr('xform.projection')
		tiltrot = tilt_euler_xform.get_rotation("eman")
		tilt_euler_xform.set_rotation({"type":"eman","az":tiltrot["az"],"alt":tiltrot["alt"],"phi":-simmx_tilt[3].get_value_at(tiltbestrefnum, tiltimgnum)})
		# Write out test results
		volume.project("standard", {"transform":untilt_euler_xform}).write_image('untilted_test.hdf', 2*imgnum)
		untiltimgs[imgnum].write_image('untilted_test.hdf', 2*imgnum+1)
		volume.project("standard", {"transform":tilt_euler_xform}).write_image('tilted_test.hdf', 2*imgnum)
		tiltimgs[imgnum].write_image('tilted_test.hdf', 2*imgnum+1)
		
		# Compute tilt transform, and account for symmetry if needed (the tilt angle must be < symmetryangle/2 for this to work
		# This is a bit of a cheat b/c it looks for the symmetric position thaty gives the best result. I wish I coould figure out how to do this w/o.....
		invuntiltxform = untilt_euler_xform.inverse()
		bestdistance = 129600	# 360**2
		besttiltangle = 0.0
		for sym in symmeties.get_syms():
			symtilt = tilt_euler_xform*sym
			tiltxform = symtilt*invuntiltxform
			currenttiltangle = tiltxform.get_rotation("eman")["alt"]
			distance = (currenttiltangle - options.tiltangle)**2
			if distance < bestdistance:
				bestdistance = distance
				besttiltangle = currenttiltangle
		#print untilt_euler_xform.get_rotation("eman"), tilt_euler_xform.get_rotation("eman"), untiltbestrefnum, tiltbestrefnum
		print "The angluar distance between tilted and untiled is: %3.2f"%besttiltangle
	exit(1)
	
	# Make contour plot to validate each particle
	ac = 0
	distplot = EMData(options.tiltrange*2+1,options.tiltrange*2+1)
	distplot.to_zero()
	for imgnum in xrange(simmx[0].get_ysize()):
		bestscore = float('inf')
		bestrefnum = 0
		for refnum in xrange(simmx[0].get_xsize()):
			if simmx[0].get_value_at(refnum, imgnum) < bestscore:
				bestscore = simmx[0].get_value_at(refnum, imgnum)
				bestrefnum = refnum
		# Get the euler angle for this particle and call compare to tilt"bdb:%s#
		euler_xform = projections[bestrefnum].get_attr('xform.projection')
		compare_to_tilt(volume, tiltimgs[imgnum], imgnum, euler_xform, simmx[3].get_value_at(bestrefnum, imgnum), distplot, options.tiltrange, 1) # For now only ints
		#Get 2D xfrom and transform, then add the image to its class avg"bdb:%s#
		xform = Transform({"type":"2d","alpha":simmx[3].get_value_at(bestrefnum, imgnum),"tx":simmx[1].get_value_at(bestrefnum, imgnum),"ty":simmx[2].get_value_at(bestrefnum, imgnum)})
		imgprocess = untiltimgs[imgnum].process("xform", {"transform":xform})
		# add an image to its CA, otherwise make a avger
		imgprocess.write_image("aligneddata.hdf",ac)
		projections[bestrefnum].write_image("aligneddata.hdf",ac+1)
		ac+=2
		
	# Make scoremx avg
	scoremxs = EMData.read_images("bdb:%s#scorematrix"%workingdir)
	avgmxavger = Averagers.get('mean')
	for mx in scoremxs:
		avgmxavger.add_image(mx)
	avgmx = avgmxavger.finish()
	avgmx.write_image("%s/contour.hdf"%workingdir)
	distplot.write_image("%s/distplot.hdf"%workingdir)
		
# Function to compute a similarity map for a given image		
def compare_to_tilt(volume, tilted, imgnum, eulerxform, zrot, distplot, tiltrange, tiltstep):
	scoremx = EMData(2*tiltrange+1,2*tiltrange+1)
	bestscore = float('inf')
	for rotx in xrange(-tiltrange, tiltrange+1,tiltstep):
		for roty in xrange(-tiltrange, tiltrange+1,tiltstep):
			# First make the projection
			tiltangle = math.sqrt(rotx*rotx + roty*roty)
			tiltaxis = math.degrees(math.atan2(roty, rotx))
			tiltxform = Transform({"type":"eman","az":tiltaxis,"alt":tiltangle,"phi":-tiltaxis})
			inplane = Transform({"type":"eman", "phi":-zrot})
			totalxform = tiltxform*inplane*eulerxform
			testprojection = volume.project("standard",totalxform)
			tiltalign = tilted.align(options.align[0],testprojection,options.align[1],options.cmp[0],options.cmp[1])
			score = tiltalign.cmp(options.cmp[0], testprojection, options.cmp[1])
			scoremx.set_value_at(rotx+tiltrange, roty+tiltrange, score)
	scoremx.write_image("bdb:%s#scorematrix"%workingdir, imgnum)
	# Denoise the contiur plot, I need to experiment around with this
	radius = 4
	scoremx_blur = scoremx.process('eman1.filter.median',{'radius':radius})
	scoremx_blur = scoremx_blur.get_clip(Region(radius, radius, scoremx_blur.get_xsize() - radius*2, scoremx_blur.get_ysize() - radius*2))
	# Find the peak
	maxpeak = scoremx_blur.calc_min_location()
	distplot.set_value_at(maxpeak[0], maxpeak[1], distplot.get_value_at(maxpeak[0], maxpeak[1]) + 1.0)

def run(command):
	"Execute a command with optional verbose output"		    
	print command
	#exit(1)
	if options.verbose>0 : print "***************",command
	error = os.system(command)
	if error==11 :
		pass		    
#		print "Segfault running %s\nNormal on some platforms, ignoring"%command
	elif error : 
		print "Error running:\n%s"%command		    
		exit(1)
		
if __name__ == "__main__":
	main()