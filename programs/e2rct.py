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
from EMAN2jsondb import *
from os import system
import sys

def main():
  
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] 
	
	This program is designed to generate a reconstruction via the random conical tilt technique.
	Starting from a tilted and untilted dataset. Thce program interfaces with e2refine2d.py to 
	find the azimuthal angles from the untilited dataset and e2make3d.py to make the 3D model
	from the untilted dataset combined with the azimuthal angles and stage tilt. A model is made
	foreach e2refine2d.py class(used for alignment)"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	#options associated with e2rctV2.py
	parser.add_header(name="rctheader", help='Options below this label are specific to e2rct', title="### e2rct options ###", row=3, col=0, rowspan=1, colspan=2)
	parser.add_argument("--path",type=str,default=None,help="Path for the rct reconstruction, default=auto")
	parser.add_argument("--untiltdata",type=str,default=None,help="Name of the tilted dataset", guitype='filebox', browser="EMSetsTable(withmodal=True,multiselect=False)", row=0, col=0, rowspan=1, colspan=2)
	parser.add_argument("--tiltdata",type=str,default=None,help="Name of the tilted dataset", guitype='filebox', browser="EMSetsTable(withmodal=True,multiselect=False)", row=1, col=0, rowspan=1, colspan=2)
	parser.add_argument("--classavg",type=str,default=None,help="Name of classavg file created by e2refine2d.py",browser="EMBrowserWidget(withmodal=True,multiselect=False)", guitype='filebox', row=2, col=0, rowspan=1, colspan=2)
	parser.add_argument("--stagetilt",type=float,default=0,help="Amount of tiliting of the cryo stage, default=0, get the stage tilt from particle attributes. Only possible if e2RCTboxer was used for particle picking", guitype='intbox', row=4, col=0, rowspan=1, colspan=1)
	parser.add_argument("--careject",type=str,default=None, help="class averages to reject, default=None", guitype='strbox', row=5, col=0, rowspan=1, colspan=1)
	parser.add_argument("--minproj",type=int,default=1,help="Minimum number of projections/images in a class average, for a class average to be used for a reconstruction, default=auto", guitype='intbox', row=5, col=1, rowspan=1, colspan=1)
	parser.add_argument("--align", action="store_true", help="Switch on image alignment. This is only translational alignment for the tilted images, and uses a iterative method similar to EMAN1 centalignint",default=False, guitype='boolbox', row=6, col=0, rowspan=1, colspan=1)
	parser.add_argument("--tiltaxis", action="store_true", help="Do a tiltaxis correction(Takes into account variations in tilt axis from micrograph to micrograph. Only possible if e2RCTboxer was used for particle picking",default=False, guitype='boolbox', row=4, col=1, rowspan=1, colspan=1)
	parser.add_argument("--maxshift",type=int,help="Maximun amount to shift the images during alignment", default=2, guitype='intbox', row=6, col=1, rowspan=1, colspan=1)
	parser.add_argument("--process",metavar="processor_name(param1=value1:param2=value2)",type=str,default=None,action="append",help="process RCT recons. Usually used to filter RCTS")
	# Options for averaging RCTs
	parser.add_header(name="avgheader", help='Options below this label are specific to e2align3d', title="### e2align3d options ###", row=7, col=0, rowspan=1, colspan=2)
	parser.add_argument("--avgrcts",action="store_true", help="If set recons from each CA will be alinged and averaged.",default=False, guitype='boolbox', row=8, col=0, rowspan=1, colspan=2)
	parser.add_argument("--reference", type=str,default=None,help="Reference used to align RCT recons to, needs to be aligned to symetry axis is --sym is specified", guitype='filebox', filecheck=False, row=9, col=0, rowspan=1, colspan=2)
	parser.add_argument("--sym", dest="sym", default="c1", help="Set the symmetry; if no value is given then the model is assumed to have no symmetry.\nChoices are: i, c, d, tet, icos, or oct.", guitype='symbox', row=10, col=0, rowspan=1, colspan=2)
	parser.add_argument("--cuda",action="store_true", help="Use CUDA for the alignment step.",default=False, guitype='boolbox', expert=True, row=12, col=0, rowspan=1, colspan=1)
	parser.add_argument("--aligngran",type=float,default=10.0,help="Fineness of global search in e2align3d.py, default=10.0", guitype='floatbox', row=11, col=0, rowspan=1, colspan=1)
	parser.add_argument("--weightrecons",action="store_true", help="Weight the reconstruction by particle numbers.",default=False, guitype='boolbox', row=11, col=1, rowspan=1, colspan=1)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	global options
	(options, args) = parser.parse_args()
	
	#print help 
	if options.untiltdata == None or options.tiltdata == None or options.classavg == None:
		parser.print_help()
		exit(0)
		
	logid=E2init(sys.argv,options.ppid)
	
	if options.careject: options.careject = options.careject.split(',')
	# initialize CUDA if desired
	if options.cuda: initializeCUDAdevice()
	
	if not options.path:
#		options.path = "rct"
		i=1
		found = 1
		while found == 1:
			if i < 10:
				run_dir = '0' + str(i)
			else:
				run_dir = str(i)
			found = os.path.exists("rct" + run_dir)
			i = i+1
		os.mkdir("rct_" + run_dir)
		options.path="rct_"+run_dir
	
	else:
		options.path = options.path+"/"

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
					print p["alpha"], tiltcorrection
				x.set_rotation({"type":"eman", "az":(p["alpha"]-tiltcorrection), "alt":options.stagetilt}) # Use supplied if available
				tiltimgs[i].set_attr("xform.projection", x)
				
			# Now do the reconstruction
			if options.align:
				centered_particles = center_particles(tiltimgs, avnum, 1)
				for r, idx in enumerate(a.get_attr('class_ptcl_idxs')):
					centered_particles[idx].write_image("%s/rctclasses_%02d.hdf" % (options.path,avnum), r)
			else:
				for r, idx in enumerate(a.get_attr('class_ptcl_idxs')):
					tiltimgs[idx].write_image("%s/rctclasses_%02d.hdf" % (options.path,avnum), r)
				
			if options.verbose>0: print "Reconstructing: %s/rctrecon_%02d.hdf" % (options.path,avnum)
			run("e2make3d.py --input=%s/rctclasses_%02d.hdf --output=%s/rctrecon_%02d.hdf --iter=2" % (options.path,avnum,options.path,avnum))
			
			# Process the RCTs usually for filtering
			currentrct = EMData("%s/rctrecon_%02d.hdf" % (options.path,avnum))
			processimage(currentrct, options.process)
			currentrct.write_image(("%s/rctrecon_%02d.hdf" % (options.path,avnum)),0)
			
			#now make an averaged 
			if options.avgrcts:
				#if options.shrink: currentrct.process_inplace('xform.scale', {'scale':options.shrink})
				if reference:
					# Do the alignment, perhaps I should use e2align3d.py?
					if options.cuda: EMData.switchoncuda()
					nbest = currentrct.xform_align_nbest('rotate_translate_3d', reference, {'delta':options.aligngran, 'dphi':options.aligngran,'sym':'c1', 'verbose':1}, 1, 'ccc.tomo',{})
					raligned = currentrct.align('refine_3d_grid', reference, {"xform.align3d":nbest[0]["xform.align3d"],"delta":1,"range":(options.aligngran/2 + 1),"verbose":1}, "ccc.tomo")
					if options.cuda: EMData.switchoffcuda()
					# write output
					arlist.append(raligned)
					raligned.write_image(("%s/rctrecon_align%02d.hdf" % (options.path,avnum)),0)
					if not options.reference:
						reference = average_rcts(arlist, totalptcls)	# Make a running average
				else:
					reference=currentrct
					arlist.append(currentrct)
	#now make the average
	if options.avgrcts:
		avged = average_rcts(arlist, totalptcls)
		if options.sym != "c1":
			symavg = avged.process('xform.applysym',{"sym":options.sym})
			symavg.write_image("%s/rctrecon_symavg.hdf" % (options.path), 0)
		avged.write_image("%s/rctrecon_avg.hdf" % (options.path), 0)
  
#	js_close_dict(options.classavg)
	E2end(logid)
	
def average_rcts(arlist, totalptcls):
	if options.verbose>0: print "Making final recon using %d class average recons" % len(arlist)
	avgr = Averagers.get('mean')		    
	for recon in arlist:
		if options.weightrecons:
			weight = len(arlist)*recon.get_attr('ptcl_repr')/totalptcls
			if options.verbose>0: print "Weighting recon using %f" % weight
			recon.mult(weight)
		avgr.add_image(recon)
	return avgr.finish()
		
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
		ptclavg.write_image("%s/centref_%02d.hdf" % (options.path,avnum), it)
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