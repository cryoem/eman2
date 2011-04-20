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
#	parser.add_option("--untiltdata",type="string",default=None,help="Name of the untilted dataset, default=auto")
	parser.add_option("--classavg",type="string",default=None,help="Name of classavg file created by e2refine2d.py, default=auto")
	parser.add_option("--stagetilt",type="float",default=None,help="Amount of tiliting of the cryo stage, default=auto")
	parser.add_option("--careject",type="int",default=None,action="append",help="class averages to reject, default=None")
	parser.add_option("--align",type="string",help="Switch on image alignment (set to 1 for True)", default=None) 
#	parser.add_option("--cmp",type="string",help="The name of a 'cmp' to be used in comparing the aligned images", default="frc")
	parser.add_option("--aligngran",type="float",default=30.0,help="Fineness of global search in e2align3d.py, default=30.0")
	parser.add_option("--weightrecons",type="int",default=1,help="weight reconstructions by ptcl_repr (before averaging), default=1")
	parser.add_option("--maxshift",type="int",help="Maximun amout to shift the images during alignment", default=2)
	parser.add_option("--minproj",type="int",default=1,help="Minimum number of projections/images in a class average, for a class average to be used for a reconstruction, default=auto")
	parser.add_option("--sym", dest="sym", default="c1", help="Set the symmetry; if no value is given then the model is assumed to have no symmetry.\nChoices are: i, c, d, tet, icos, or oct.")
	parser.add_option("--preprocess",metavar="processor_name(param1=value1:param2=value2)",type="string",default=None,action="append",help="preprocess recons before alignment")
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	global options
	(options, args) = parser.parse_args()
	
	if not options.path:
		options.path = "bdb:rct/"
	else:
		options.path = "bdb:"+options.path+"/"
	
	if options.stagetilt:
		tiltangle = options.stagetilt
	else:
		print "Error, stagetilt parameter not supplied! Crashing!"
		exit(1)
	    
	if not options.tiltdata:
		print "Error, tiltdata needed! Crashing!"
		exit(1)
	
	if not options.classavg:
		print "Error, classavg needed! Crashing!"
		exit(1)	    

	#set up the preprocesses (otherwise bugs occur)
	preprocess = ""
	if options.preprocess != None:
		for p in reversed(options.preprocess):
			preprocess = ("--preprocess=%s " % str(p)) + preprocess 
		
	# Now get azimuthal data
	cadb = EMData.read_images(options.classavg)
	tiltimgs = EMData.read_images(options.tiltdata)
	
	#initialize some stuff
	search = tiltimgs[0].get_attr('nx')/4 # the search range for the 3D aligner
	arlist = []
	totalptcls = 0
	firstrecon = None

	for avnum, a in enumerate(cadb):
		if a == None:
			break # stupidly, enumerate keeps on going, even after it reachs the end of the 'list'
		if options.careject and avnum in options.careject:
			if options.verbose>0: print "Rejecting class average %d" % avnum
			continue # do not use suboptimal classaverages
		if a.get_attr('ptcl_repr') > options.minproj:
			totalptcls = a.get_attr('ptcl_repr') + totalptcls
			for inum, i in enumerate(a.get_attr('class_ptcl_idxs')):
				img = EMData()  
				img.read_image(options.untiltdata, i)
				imga = img.align('rotate_translate', a)
				t = imga.get_attr("xform.align2d")
				t.invert()
				p = t.get_params("2d")
				x = Transform()                
				x.set_rotation({"type":"eman", "az":p["alpha"], "alt":tiltangle})
				tiltimgs[i].set_attr("xform.projection", x)
				tiltimgs[i].write_image("%srctclasses_%02d" % (options.path,avnum), inum)
				# Now do the reconstruction
			if options.align:
				imglist = EMData.read_images("%srctclasses_%02d" % (options.path,avnum))
				imglist.sort(lambda a, b: cmp(a.get_attr('xform.projection').get_rotation().get('az'), b.get_attr('xform.projection').get_rotation().get('az')))
				for r in xrange(len(imglist)):
					aligned = imglist[r].align("translational", imglist[r-1], {"maxshift":options.maxshift}) 
					aligned.write_image("%srctccarejectlassesali_%02d" % (options.path,avnum), r)	            
				if options.verbose>0: print "Reconstructing: %srctreconali_%02d" % (options.path,avnum)
				run("e2make3d.py --input=%srctclassesali_%02d --output=%srctrecon_%02d --sym=%s --iter=2" % (options.path,avnum,options.path,avnum,options.sym))
			else:
				if options.verbose>0: print "Reconstructing: %srctrecon_%02d" % (options.path,avnum)
				run("e2make3d.py --input=%srctclasses_%02d --output=%srctrecon_%02d --sym=%s --iter=2" % (options.path,avnum,options.path,avnum,options.sym))
	#now make an averaged image    
			if firstrecon:
				run("e2align3d.py %s %srctrecon_%02d %srctreconali_%02d --delta=%f --dphi=%f --search=%d --sym=%s %s --cmp='frc'" % (firstrecon,options.path,avnum,options.path,avnum,options.aligngran, options.aligngran, search,options.sym,preprocess))
				arlist.append("%srctreconali_%02d" % (options.path,avnum))
			else:
				tmp = EMData()
				tmp.read_image("%srctrecon_%02d" % (options.path,avnum))
				processimage(tmp, options.preprocess)
				tmp.write_image("%srctreconali_%02d" % (options.path,avnum),0)
				firstrecon = ("%srctreconali_%02d" % (options.path,avnum))
				arlist.append(firstrecon)
		    
	#now make the average
	if options.verbose>0: print "Making final recon using %d class average recons" % len(arlist)
	avgr = Averagers.get('mean')		    
	for r in arlist:
		recon = EMData()
		recon.read_image(r)
		if options.weightrecons:
			weight = len(arlist)*recon.get_attr('ptcl_repr')/totalptcls
			if options.verbose>0: print "Weighting recon using %f" % weight
			recon.mult(weight)
		avgr.add_image(recon)
	avged = avgr.finish()
	avged.write_image("%srctrecon_avg" % (options.path), 0)
  
	db_close_dict(options.classavg)    
	        
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
		print "Error running:\n%s"%comman		    
		exit(1)
		
def processimage(image, optionspreprocess):
	if optionspreprocess != None:
	    for p in options.preprocess:
	        try:
		    (processorname, param_dict) = parsemodopt(p)
		    if not param_dict : param_dict={}
		    image.process_inplace(str(processorname), param_dict)
		except:
		    print "warning - application of the pre processor",p," failed. Continuing anyway"  
	
if __name__ == "__main__":
	main()