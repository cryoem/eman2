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
from optparse import OptionParser
from os import system, unlink, getcwd
from glob import glob

def main():
  
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] 
	
	This program is designed to generate a reconstruction via the random conical tilt technique.
	Starting from a tilted and untilted dataset. The program interfaces with e2refine2d.py to 
	find the azimuthal angles from the untilited dataset and e2make3d.py to make the 3D model
	from the untilted dataset combined with the azimuthal angles and stage tilt. A model is made
	foreach e2refine2d.py class(used for alignment)"""
        parser = OptionParser(usage=usage,version=EMANVERSION)
	
	#options associated with e2rct.py
	parser.add_option("--path",type="string",default=None,help="Path for the rct reconstruction, default=auto")
	parser.add_option("--tiltdata",type="string",default=None,help="Name of the tilted dataset, default=auto")
#	parser.add_option("--untiltdata",type="string",default=None,help="Name of the untilted dataset, default=auto")
	parser.add_option("--simmx",type="string",default=None,help="Name of simmex file created by e2refine2d.py, default=auto")
	parser.add_option("--stagetilt",type="float",default=None,help="Amount of tiliting of the cryo stage, default=auto")
	parser.add_option("--align",type="string",help="Switch on image alignment (set to 1 for True)", default=None) 
	parser.add_option("--cmp",type="string",help="The name of a 'cmp' to be used in comparing the aligned images", default="frc")
	parser.add_option("--maxshift",type="int",help="Maximun amout to shift the images during alignment", default=2)
        parser.add_option("--minproj",type="int",default=1,help="Minimum number of projections in each reconstruction, default=auto")
	parser.add_option("--sym", dest="sym", default="c1", help="Set the symmetry; if no value is given then the model is assumed to have no symmetry.\nChoices are: i, c, d, tet, icos, or oct.")
	parser.add_option("--postprocess", metavar="processor_name(param1=value1:param2=value2)", type="string", action="append", help="postprocessor to be applied to the 3D volume once the reconstruction is completed. There can be more than one postprocessor, and they are applied in the order in which they are specified. See e2help.py processors for a complete list of available processors.")
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
	
	if not options.simmx:
	    print "Error, simmx needed! Crashing!"
	    exit(1)
	
	
	# Now get azimuthal data
	smxdb = db_open_dict(options.simmx)
	data = []
	
	for r in xrange(smxdb[0].get_attr('ny')):
	    bestscore = 1
	    bestclass = 0
	    for c in xrange(smxdb[0].get_attr('nx')):
	        score = smxdb[0].get_value_at(c, r)
	        if score < bestscore:
	            bestscore = score
	            bestclass = c
	    if options.verbose >1: print "The bestscore is: %f for class %d" % (bestscore, bestclass)        
	    data.append((bestclass, smxdb[3].get_value_at(bestclass, r)))

	tiltimgs = EMData.read_images(options.tiltdata)
	
	# First delete any previous dicts and images from previous runs (prevents any data 'merging')
	for c in xrange(smxdb[0].get_attr('nx')):
	    try:
	        for f in glob("%s/%sEMAN2DB/rctclasses_%02d*" % (getcwd(), options.path[4:len(options.path)],c)):
	            os.unlink(f)
	    except:
	        pass
	   
	# Now make the stacks and write the transforms
	classpop = [0] * smxdb[0].get_attr('nx')
	for r in xrange(smxdb[0].get_attr('ny')):
	    t = Transform()
	    t.set_rotation({"type":"eman", "az":data[r][1], "alt":tiltangle})
	    tiltimgs[r].set_attr("xform.projection", t)
	    #print tiltimgs[i], classpop[data[i][0]], data[i][0]
	    tiltimgs[r].write_image("%srctclasses_%02d" % (options.path,data[r][0]), classpop[data[r][0]])
	    classpop[data[r][0]] += 1
	    
	db_close_dict(options.simmx)
	    	    	
	# For each class average, sort the images by azimuthal angle, align to nearest neighbor, and reconstruct    	    	
	for c in xrange(smxdb[0].get_attr('nx')):
	    if db_check_dict("%srctclasses_%02d"  % (options.path,c)):
	        if  classpop[c] >= options.minproj:
		    # Do an alignment of nieghboring images
		    if options.align != None:
		        if options.verbose>0: print "Starting alignment for %srctrecon_%02d" % (options.path,c)
		        # Now translatrionaly align the tilted images to each other
		        imglist = EMData.read_images("%srctclasses_%02d" % (options.path,c))
		        imglist.sort(lambda a, b: cmp(a.get_attr('xform.projection').get_rotation().get('az'), b.get_attr('xform.projection').get_rotation().get('az')))
		        
		        # Align each image to its nearest neighbor (in Fourier Space) and compute comparison, also make class averages(to check sanity)
		        avgr = Averagers.get('mean')
		        iniangle = imglist[0].get_attr('xform.projection').get_rotation().get('az')
		        ai = 0
		        for r in xrange(len(imglist)):
		            aligned = imglist[r].align("translational", imglist[r-1], {"maxshift":options.maxshift}) 
		            aligned.write_image("%srctclassesali_%02d" % (options.path,c), r)
		            sim = aligned.cmp(options.cmp, imglist[r-1])
		            
		            if imglist[r].get_attr('xform.projection').get_rotation().get('az') > iniangle:
			        ref = avgr.finish()
			        if ref != None:
			            ref.set_attr('iniangle',iniangle)
			            ref.write_image("%srctclassavg_%02d" % (options.path,c), ai)
			            ai += 1
			        iniangle = imglist[r].get_attr('xform.projection').get_rotation().get('az')    
			        avgr = Averagers.get('mean')
			    else:
			        avgr.add_image(aligned)
			        
		            if options.verbose>1: 
		                print "Image %d has an %s of %f with image %d" % (r,options.cmp,sim,(r-1)) 
		                #print aligned.get_attr("xform.align2d").get_params("2d")
		        
	                # For now we will just use default parms
	                if options.verbose>0: print "Reconstructing: %srctrecon_%02d" % (options.path,c)
	                run("e2make3d.py --input=%srctclassesali_%02d --output=%srctrecon_%02d --sym=%s --iter=2" % (options.path,c,options.path,c,options.sym))
	            else:
		      	# For now we will just use default parms  
		      	if options.verbose>0: print "Reconstructing: %srctrecon_%02d" % (options.path,r)
	                run("e2make3d.py --input=%srctclasses_%02d --output=%srctrecon_%02d --sym=%s --iter=2" % (options.path,c,options.path,c,options.sym))
	
def run(command):
	"Execute a command with optional verbose output"
	global options
	if options.verbose>0 : print "***************",command
	error = system(command)
	if error==11 :
		pass
#		print "Segfault running %s\nNormal on some platforms, ignoring"%command
	elif error : 
		print "Error running:\n%s"%command
		exit(1)
	
if __name__ == "__main__":
    main()