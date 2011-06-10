#!/usr/bin/env python

#
# Author: John Flanagan, 24/09/2010 (jfflanag@bcm.edu)
# Copyright (c) 2010 Baylor College of Medicine
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
# GNU General Public License for more details../e2rct.py --tiltdata=yRiboRed_tilted.spi --simmx="bdb"simmx_06" --stagetilt=60
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#


from EMAN2 import *
from EMAN2db import db_close_dict
from optparse import OptionParser
from os import system

def main():
  
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] 
	
	This program is designed to generate a reconstruction via the random conical tilt technique.
	Starting from a tilted and untilted dataset. Thce program interfaces with e2refine2d.py to 
	find the azimuthal angles from the untilited dataset and e2make3d.py to make the 3D model
	from the untilted dataset combined with the azimuthal angles and stage tilt. A model is made
	foreach e2refine2d.py class(used for alignment)"""
	parser = OptionParser(usage=usage,version=EMANVERSION)
	
	#options associated with e2rctV2.py
	parser.add_option("--path",type="string",default=None,help="Path for the rct reconstruction, default=auto")
	parser.add_option("--tiltdata",type="string",default=None,help="Name of the tilted dataset, default=auto")
	parser.add_option("--untiltdata",type="string",default=None,help="Name of the tilted dataset, default=auto")
	parser.add_option("--classavg",type="string",default=None,help="Name of classavg file created by e2refine2d.py, default=auto")
	parser.add_option("--stagetilt",type="float",default=None,help="Amount of tiliting of the cryo stage, default=auto")
	parser.add_option("--careject",type="string",default=None, help="class averages to reject, default=None")
	parser.add_option("--minproj",type="int",default=1,help="Minimum number of projections/images in a class average, for a class average to be used for a reconstruction, default=auto")
	parser.add_option("--align", action="store_true", help="Switch on image alignment.",default=False)
	parser.add_option("--tiltaxis", action="store_true", help="Do a tiltaxis correction(Takes into account variations in tilt axis from micrograph to micrograph.",default=False)
	parser.add_option("--maxshift",type="int",help="Maximun amount to shift the images during alignment", default=2)
	# Options for averaging RCTs
	parser.add_option("--avgrcts",action="store_true", help="If set recons from each CA will be alinged and averaged.",default=False)
	parser.add_option("--reference", type="string",default=None,help="Reference used to align RCT recons to, needs to be aligned to symetry axis is --sym is specified")
	parser.add_option("--sym", dest="sym", default="c1", help="Set the symmetry; if no value is given then the model is assumed to have no symmetry.\nChoices are: i, c, d, tet, icos, or oct.")
	parser.add_option("--cuda",action="store_true", help="Use CUDA for the alignment step.",default=False)
	parser.add_option("--aligngran",type="float",default=10.0,help="Fineness of global search in e2align3d.py, default=10.0")
	parser.add_option("--weightrecons",action="store_true", help="Weight the reconstruction by particle numbers.",default=False)
	parser.add_option("--preprocess",metavar="processor_name(param1=value1:param2=value2)",type="string",default=None,action="append",help="preprocess recons before alignment")
	parser.add_option("--verbose", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	global options
	(options, args) = parser.parse_args()
	if options.careject: options.careject = options.careject.split(',')
	
	if not options.path:
		options.path = "bdb:rct/"
	else:
		options.path = "bdb:"+options.path+"/"
	    
	if not options.tiltdata:
		print "Error, tiltdata needed! Crashing!"
		exit(1)
	
	if not options.classavg:
		print "Error, classavg needed! Crashing!"
		exit(1)	    
		
	# Now get azimuthal data
	cadb = EMData.read_images(options.classavg)
	tiltimgs = EMData.read_images(options.tiltdata)

	# Check to see if we need a tilt angle
	if not options.stagetilt:
		if not tiltimgs[0].get_attr_dict().has_key("tiltangle"):	
			print "Error, stagetilt parameter not supplied and could find tiltangle as an attribute! Crashing!"
			exit(1)
	
	#initialize some stuff
	search = tiltimgs[0].get_attr('nx')/4 # the search range for the 3D aligner
	arlist = []
	totalptcls = 0
	tiltcorrection = 0.0
	reference = None
	# Read in reference if nedded
	if options.reference: 
		reference = EMData(options.reference)
	# commence with the RCT recons
	for avnum, a in enumerate(cadb):
		if a == None:
			break # stupidly, enumerate keeps on going, even after it reachs the end of the 'list'
		if options.careject and str(avnum) in options.careject:
			if options.verbose>0: print "Rejecting class average %d" % avnum
			continue # do not use suboptimal classaverages
		if a.get_attr('ptcl_repr') > options.minproj:
			totalptcls = a.get_attr('ptcl_repr') + totalptcls
			for inum, i in enumerate(a.get_attr('class_ptcl_idxs')):
				img = EMData()  
				img.read_image(options.untiltdata, i)
				imga = img.align('rotate_translate', a)
				t = imga.get_attr("xform.align2d")
				t.invert()	# Is this really necessary?
				p = t.get_params("2d")
				x = Transform()
				
				# Use the supplied tilt angle and tilt axis if availabe
				if not options.stagetilt:
					options.stagetilt = tiltimgs[i].get_attr("tiltangle")
				if options.tiltaxis:
					tiltcorrection = tiltimgs[i].get_attr_default("tiltaxis",0.0)
				x.set_rotation({"type":"eman", "az":(p["alpha"]-tiltcorrection), "alt":options.stagetilt}) # Use supplied if available
				tiltimgs[i].set_attr("xform.projection", x)
				
			# Now do the reconstruction
			if options.align:
				centered_particles = center_particles(tiltimgs, avnum, 1)
				for r, idx in enumerate(a.get_attr('class_ptcl_idxs')):
					centered_particles[idx].write_image("%srctclasses_%02d" % (options.path,avnum), r)
			else:
				for r, idx in enumerate(a.get_attr('class_ptcl_idxs')):
					tiltimgs[idx].write_image("%srctclasses_%02d" % (options.path,avnum), r)
				
			if options.verbose>0: print "Reconstructing: %srctrecon_%02d" % (options.path,avnum)
			run("e2make3d.py --input=%srctclasses_%02d --output=%srctrecon_%02d --iter=2" % (options.path,avnum,options.path,avnum))
			#now make an averaged image
			if options.avgrcts:
				currentrct = EMData("%srctrecon_%02d" % (options.path,avnum))
				processimage(currentrct, options.preprocess)
				#if options.shrink: currentrct.process_inplace('xform.scale', {'scale':options.shrink})
				if reference:
					# Do the alignment
					if options.cuda: EMData.switchoncuda()
					nbest = currentrct.xform_align_nbest('rotate_translate_3d', reference, {'delta':options.aligngran, 'dphi':options.aligngran,'sym':'c1', 'verbose':1}, 1, 'ccc.tomo',{})
					raligned = currentrct.align('refine_3d_grid', reference, {"xform.align3d":nbest[0]["xform.align3d"],"delta":1,"range":6,"verbose":1}, "ccc.tomo")
					if options.cuda: EMData.switchoffcuda()
					# write output
					arlist.append(raligned)
					raligned.write_image(("%srctrecon_align%02d" % (options.path,avnum)),0)
				else:
					reference=currentrct
					arlist.append(currentrct)
	#now make the average
	if options.avgrcts:
		if options.verbose>0: print "Making final recon using %d class average recons" % len(arlist)
		avgr = Averagers.get('mean')		    
		for recon in arlist:
			if options.weightrecons:
				weight = len(arlist)*recon.get_attr('ptcl_repr')/totalptcls
				if options.verbose>0: print "Weighting recon using %f" % weight
				recon.mult(weight)
			avgr.add_image(recon)
		avged = avgr.finish()
		if options.sym != "c1":
			symavg = avged.process('xform.applysym',{"sym":options.sym})
			symavg.write_image("%srctrecon_symavg" % (options.path), 0)
		avged.write_image("%srctrecon_avg" % (options.path), 0)
  
	db_close_dict(options.classavg)

# I am not worried about interpoloation errors b/c the RCT resolution is very low anyways (80/20 rule....)
def center_particles(particles, avnum, iterations):
	if options.verbose>0: print "Centering tilted particles"
	centeredimgs = []
	radius = particles[0].get_attr("nx")/2 # nx = ny, always.....
	for it in xrange(iterations):
		ptclavgr = Averagers.get('mean')
		# Make average
		for img in particles:
			ptclavgr.add_image(img)
		ptclavg = ptclavgr.finish()
		ptclavg.process_inplace("math.rotationalaverage")
		ptclavg.write_image("%scentref_%02d" % (options.path,avnum), it)
		#align translationally align particles
		for img in particles:
			#img.process_inplace("mask.noise", {"dx":0,"dy":0,"outer_radius":radius})
			centeredimgs.append(img.align("translational", ptclavg, {"maxshift":options.maxshift}))
		particles = centeredimgs # Change the reference (the source data will not be affected)
		centeredimgs = []
		
	return particles # Not the same as was passed in (rereferenced)
	
def run(command):
	"Execute a command with optional verbose output"
	global options		    
	print command
	if options.verbose>0 : print "***************",command
	error = system(command)
	if error==11 :
		pass		    
#		print "Segfault running %s\nNormal on some platforms, ignoring"%command
	elif error : 
		print "Error running:\n%s"%command		    
		exit(1)
		
def processimage(image, options_to_preprocess):
	if options_to_preprocess != None:
	    for p in options_to_preprocess:
	        try:
		    (processorname, param_dict) = parsemodopt(p)
		    if not param_dict : param_dict={}
		    image.process_inplace(str(processorname), param_dict)
		except:
		    print "warning - application of the pre processor",p," failed. Continuing anyway"
	else:
		return
	
if __name__ == "__main__":
	main()