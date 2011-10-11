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
	usage = """prog [options]
	Tiltvalidation using Richard Hendersons technique. To use a stack of untilted and tiltimages whose set relationship is one-to-one is required along with a
	volume to validate. A image whose x and y axes show the tiltvalidation result. A valid reconstruction will have a peak at the magnitude of the stagetilt and
	along the tiltaxis. The tiltaxis can be found using e2RCTboxer. The tiltaxis is computed by it and displayed in the e2RCTboxer GUI.
	
	Output is scorematrix, and perparticletilts, a list of angluar distances between tilted and untilted paricles.
	
	For more information see:
	
	Optimal determination of particle orientation, absolute hand, and contrast loss in 
	single-particle electron cryomicroscopy.
	Rosenthal PB, Henderson R.
	J Mol Biol. 2003 Oct 31;333(4):721-45 
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	# options associated with e2tiltvalidate.py
	parser.add_header(name="tvheader", help='Options below this label are specific to e2tiltvalidate', title="### e2tiltvalidate options ###", row=3, col=0, rowspan=1, colspan=2)
	parser.add_argument("--path", type=str,help="The folder the results are placed", default="TiltValidate")
	parser.add_argument("--volume", type=str,help="3D volume to validate",default=None, guitype='filebox', row=2, col=0, rowspan=1, colspan=2)
	parser.add_argument("--untiltdata", type=str,help="Stack of untilted images",default=None, guitype='filebox', row=0, col=0, rowspan=1, colspan=2)
	parser.add_argument("--tiltdata", type=str,help="Stack of tilted images",default=None, guitype='filebox', row=1, col=0, rowspan=1, colspan=2)
	parser.add_argument("--align", type=str,help="The name of a aligner to be used in comparing the aligned images",default="translational", guitype='comboparambox', choicelist='re_filter_list(dump_aligners_list(),\'refine|3d\', 1)', row=6, col=0, rowspan=1, colspan=2)
	parser.add_argument("--cmp", type=str,help="The name of a 'cmp' to be used in comparing the aligned images",default="ccc", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', True)', row=7, col=0, rowspan=1, colspan=2)
	parser.add_argument("--tiltrange", type=int,help="The angular tiltranger to search",default=15, guitype='intbox', row=4, col=0, rowspan=1, colspan=1)
	parser.add_argument("--sym",  type=str,help="The recon symmetry", default="c1", guitype='symbox', row=5, col=0, rowspan=1, colspan=2)
	parser.add_argument("--planethres", type=float, help="Maximum out of plane threshold for the tiltaxis. 0 = perfectly in plane, 1 = normal to plane", default=0.1, guitype='floatbox', row=4, col=1, rowspan=1)
	# options associated with e2projector3d.py
	parser.add_header(name="projheader", help='Options below this label are specific to e2project', title="### e2project options ###", row=9, col=0, rowspan=1, colspan=2)
	parser.add_argument("--delta", type=float,help="The angular step size for alingment", default=20.0, guitype='floatbox', row=10, col=0, rowspan=1, colspan=1)
	# options associated with e2simmx.py
	parser.add_header(name="simmxheader", help='Options below this label are specific to e2simmx', title="### e2simmx options ###", row=11, col=0, rowspan=1, colspan=2)
	parser.add_argument("--simalign",type=str,help="The name of an 'aligner' to use prior to comparing the images (default=rotate_translate)", default="rotate_translate", guitype='comboparambox', choicelist='re_filter_list(dump_aligners_list(),\'refine|3d\', 1)', row=14, col=0, rowspan=1, colspan=2)
	parser.add_argument("--simaligncmp",type=str,help="Name of the aligner along with its construction arguments (default=ccc)",default="ccc", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', True)', row=15, col=0, rowspan=1, colspan=2)
	parser.add_argument("--simcmp",type=str,help="The name of a 'cmp' to be used in comparing the aligned images (default=ccc)", default="ccc", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', True)', row=13, col=0, rowspan=1, colspan=2 )
	parser.add_argument("--simralign",type=str,help="The name and parameters of the second stage aligner which refines the results of the first alignment", default=None, guitype='comboparambox', choicelist='re_filter_list(dump_aligners_list(),\'refine|3d\', 1)', row=16, col=0, rowspan=1, colspan=2)
	parser.add_argument("--simraligncmp",type=str,help="The name and parameters of the comparitor used by the second stage aligner. (default=dot).",default="dot", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', True)', row=17, col=0, rowspan=1, colspan=2)
	parser.add_argument("--shrink", dest="shrink", type = int, default=0, help="Optionally shrink the input particles by an integer amount prior to computing similarity scores. For speed purposes. Defulat = 0, no shrinking", guitype='intbox', row=12, col=0, rowspan=1, colspan=1)
	
	parser.add_argument("--parallel",type=str,help="Parallelism string",default=None, guitype='strbox', row=8, col=0, rowspan=1, colspan=2)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
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
	
	logid=E2init(sys.argv,options.ppid)
	
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
	tdb = db_open_dict("bdb:%s#perparticletilts"%workingdir)
	particletilt_list = []
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
		
		# Compute tiltxis and tiltangle for each particel pair. For symmetric objects the tilttransform is selected as the one that has a tiltaxis
		# closest to perpidicular to the imaging plane (this is how Richard handles it). 
		bestinplane = 1.0
		for sym in symmeties.get_syms():
			tiltxform = tilt_euler_xform*sym.inverse()*untilt_euler_xform.inverse()
			# Select the symmetry solution whose tiltaxis is in plane
			if math.fabs(tiltxform.get_rotation("spin")["n3"]) < bestinplane:
				bestinplane = math.fabs(tiltxform.get_rotation("spin")["n3"])
				besttiltangle = tiltxform.get_rotation("spin")["Omega"]
				besttiltaxis = math.degrees(math.atan2(tiltxform.get_rotation("spin")["n2"],tiltxform.get_rotation("spin")["n1"]))
			#print "\t",tiltxform.get_rotation("spin")["Omega"],tiltxform.get_rotation("spin")["n1"],tiltxform.get_rotation("spin")["n2"],tiltxform.get_rotation("spin")["n3"]
		if bestinplane > options.planethres:
			#print "Rejecting solution"
			continue
		print "The best angle is %f with a tiltaxis of %f"%(besttiltangle,besttiltaxis)
		particletilt_list.append([imgnum, besttiltangle,besttiltaxis])

	tdb["particletilt_list"] = particletilt_list
	tdb.close()
	
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
	
	E2end(logid)
		
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
	error = launch_childprocess(command)
	if error==11 :
		pass		    
#		print "Segfault running %s\nNormal on some platforms, ignoring"%command
	elif error : 
		print "Error running:\n%s"%command		    
		exit(1)
		
if __name__ == "__main__":
	main()