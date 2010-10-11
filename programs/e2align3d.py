#!/usr/bin/env python

#
# Author: John Flanagan, 10/08/2010 (jfflanag@bcm.edu)
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
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

from EMAN2 import *
from optparse import OptionParser
import sys, os
	
def main():
  
	progname = os.path.basename(sys.argv[0])
	usage = """%prog fixed_model moving_model output_model [options] 
	
	This program is designed to rotationally and translationally align two 3D models 
	Usually the two models are shrunk down to speed things up, then a global exhaustive 
	search is down by the reine.3d.sphere and then this rough alignment is refined using
	refine.3d using the full size maps. The refiner is much quicker than the global aligner
	as it uses a simplex algoritm to bring the alignment downhill, but if your inital global
	alignmnet is too rough then the refiner might get stuck in a local minima."""
        parser = OptionParser(usage=usage,version=EMANVERSION)
        
        #options associated with e2refine3d.py
        parser.add_option("--shrink",type="float",default=None,help="Fractional amount to shrink the maps by, default=-1.0")
        parser.add_option("--preprocess",metavar="processor_name(param1=value1:param2=value2)",type="string",default=None,action="append",help="preprocess maps before alignment")
        #options form the sphere alinger
        parser.add_option("--delta",type="float",default=30.0,help="step size for the orrientation generator, default=10.0")
        parser.add_option("--dphi",type="float",default=30.0,help="step size for the inplane angle phi, default=10.0")
        parser.add_option("--rphi",type="float",default=180.0,help="search range for the inplane angle phi, default=180.0")
        parser.add_option("--search",type="int",default=10,help="step size for translational search, default=5")
        parser.add_option("--sym",type="string",default='c1',help="model symmetry (using sym, if present, speeds thing up a lot), default='c1'")
        parser.add_option("--dotrans",type="int",default=1,help="Do translational search, default=1")
        parser.add_option("--verbose",type="int",default=0,help="Be Verbose?, default=0")
        #options associated with  the simplex 3D refiner
        parser.add_option("--stepalt",type="float",default=5.0,help="step size for alt angle, default=5.0")
        parser.add_option("--stepaz",type="float",default=5.0,help="step size for az angle, default=5.0")
        parser.add_option("--stepphi",type="float",default=5.0,help="step size for inplane phi angle, default=5.0")
        parser.add_option("--stepx",type="float",default=1.0,help="step size for th x direction, default=1.0")
        parser.add_option("--stepy",type="float",default=1.0,help="step size for th y direction, default=1.0")
        parser.add_option("--stepz",type="float",default=1.0,help="step size for th z direction, default=1.0")
        parser.add_option("--maxshift",type="int",default=-1.0,help="maximum shift, (-1 means dim/4), default=-1.0")
        parser.add_option("--cmp",type="string",default='ccc',help="comparitor to use for the 3D refiner, default='ccc'")

        global options
	(options, args) = parser.parse_args()
        
        #read in the input maps and check for sanity
        fixed = EMData()
        moving = EMData()
        
        try:
	    fixed.read_image(sys.argv[1])
	except:
	    print "Not able to read file %s" % sys.argv[1]
	    exit(1)
	    
	try:
	    moving.read_image(sys.argv[2])
	except:
	    print "Not able to read file %s" % sys.argv[2]
	    exit(1)
	    
	if (fixed.get_attr('nx') != fixed.get_attr('ny') != fixed.get_attr('nz')):
	    print "Fixed map must have cubic dimensions!"
	    exit(1)
	    
	if (moving.get_attr('nx') != moving.get_attr('ny') != moving.get_attr('nz')):
	    print "Fixed map must have cubic dimensions!"
	    exit(1)
	    
	if (moving.get_attr('nx') != fixed.get_attr('nx')):
	    print "Fixed and model maps must have the same dimensions!"
	    exit(1)
	    
	#now shrink the maps by a sane value
	if options.shrink > 1:
	    print "The idea is to shrink, not blow up the models! Continuing anyways though...."
	if options.shrink <= 0:  
	    if fixed.get_attr('nx') > 250:
	        options.shrink = 0.1
	    if fixed.get_attr('nx') <= 250:
	        options.shrink = 0.25
	    if fixed.get_attr('nx') <= 100:
	        options.shrink = 0.5
	    if fixed.get_attr('nx') <= 50:
	        options.shrink = 1.0
	
	#preprocess maps
	if options.preprocess != None:
	    for p in options.preprocess:
	        try:
		    (processorname, param_dict) = parsemodopt(p)
		    if not param_dict : param_dict={}
		    fixed.process_inplace(str(processorname), param_dict)
		    moving.process_inplace(str(processorname), param_dict)
		except:
		    print "warning - application of the pre processor",p," failed. Continuing anyway"
	
	sfixed = fixed.process('xform.scale', {'scale':options.shrink})
	smoving = moving.process('xform.scale', {'scale':options.shrink})
	
	#globally align the two maps
	gsaligned = smoving.align('rt.3d.sphere', sfixed, {'delta':options.delta,'dotrans':options.dotrans,'dphi':options.dphi,'search':options.search, 'rphi':options.rphi, 'sym':options.sym, 'verbose':options.verbose})
	
	#now refine the alignment
	if options.cmp == "dot":
	    cmpparms = {'normalize':1}
	else:
	    cmpparms = {}
	    
	t = gsaligned.get_attr('xform.align3d')
	galigned = moving.process('xform',{'transform':t})
	galignedref = galigned.align('refine.3d', fixed, {'maxshift':options.maxshift, 'stepalt':options.stepalt,'stepaz':options.stepaz,'stepphi':options.stepphi,'stepx':options.stepx,'stepy':options.stepy,'stepz':options.stepz}, options.cmp, cmpparms) 
	
	#now write out the aligned model
        if sys.argv[3]:
	    outfile = sys.argv[3]
	else:
	    outfile = 'alignedmodel.hdf'
	    
	#if output is mrc format 
	if outfile[-4:].lower() == ".mrc":
	    galignedref.set_attr('UCSF.chimera',1)
  
	galignedref.write_image(outfile, 0)
	
if __name__ == "__main__":
    main()
