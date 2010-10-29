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
        parser.add_option("--shrink",type="float",default=-1,help="Fractional amount to shrink the maps by (-1 = auto), default=-1.0")
        parser.add_option("--preprocess",metavar="processor_name(param1=value1:param2=value2)",type="string",default=None,action="append",help="preprocess maps before alignment")
        parser.add_option("--maskrad",type="int",default=-1,help="Mask the recon using a spherical Gaussian mask (-1 = None), default=-1.0")
        parser.add_option("--maskfoff",type="float",default=0.1,help="Fall offf of the Gaussian mask, default=0.1")
        parser.add_option("--nsolns",type="int",default=1,help="number of peaks in the global search to refine, default=1.0")
        parser.add_option("--famps",type="float",default=1,help="fraction of Fourier amps to exclude from recons, default=0.0")
        parser.add_option("--prec",type="float",default=0.01,help="Precison to determine what solutions are the 'same', default=0.01")
        #options form the sphere alinger
        parser.add_option("--delta",type="float",default=30.0,help="step size for the orrientation generator, default=10.0")
        parser.add_option("--dphi",type="float",default=30.0,help="step size for the inplane angle phi, default=10.0")
        parser.add_option("--lphi",type="float",default=0.0,help="lower bound for the inplane angle phi, default=0.0")
        parser.add_option("--uphi",type="float",default=359.0,help="Upper bound for the inplane angle phi, default=359.0")
        parser.add_option("--search",type="int",default=10,help="maximum extent of the translational search, default=5")
        parser.add_option("--sym",type="string",default='c1',help="model symmetry (using sym, if present, speeds thing up a lot), default='c1'")
        parser.add_option("--cmp",type="string",default='ccc',help="comparitor to use for the 3D refiner, default='ccc'")
        parser.add_option("--cmpparms",type="string",action="append",default=None,help="comparitor paramters")
        parser.add_option("--dotrans",type="int",default=1,help="Do translational search, default=1")
        #options associated with  the simplex 3D refiner
        parser.add_option("--stepalt",type="float",default=5.0,help="step size for alt angle, default=5.0")
        parser.add_option("--stepaz",type="float",default=5.0,help="step size for az angle, default=5.0")
        parser.add_option("--stepphi",type="float",default=5.0,help="step size for inplane phi angle, default=5.0")
        parser.add_option("--stepx",type="float",default=1.0,help="step size for th x direction, default=1.0")
        parser.add_option("--stepy",type="float",default=1.0,help="step size for th y direction, default=1.0")
        parser.add_option("--stepz",type="float",default=1.0,help="step size for th z direction, default=1.0")
        parser.add_option("--maxshift",type="int",default=-1.0,help="maximum shift, (-1 means dim/4), default=-1.0")
        parser.add_option("--maxiter",type="int",default=100,help="maximum number of iterations(you'll need more for courser global searches), default=100")
        parser.add_option("--rcmp",type="string",default='ccc',help="comparitor to use for the 3D refiner, default='ccc'")
<<<<<<< e2align3d.py
        parser.add_option("--rcmpparms",type="string",action="append",default=None,help="refine comparitor paramters")
=======
        parser.add_option("--verbose","-v",type="int",default=0,help="Level of verboseness, default=0")
>>>>>>> 1.25

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
	    
	if (fixed.get_attr('nx') != fixed.get_attr('ny') != fixed.get_a+ options.precttr('nz')):
	    print "Fixed map must have cubic dimensions!"
	    exit(1)
	    
	if (moving.get_attr('nx') != moving.get_attr('ny') != moving.get_attr('nz')):
	    print "Fixed map must have cubic dimensions!"
	    exit(1)	#(zzz, cmppars) = parsemodopt(options.cmpparms)	 
	    
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
	
	sfixed = fixed.process('xform.scale', {'scale':options.shrink})
	smoving = moving.process('xform.scale', {'scale':options.shrink})
        options.maskrad = options.maskrad*options.shrink
        
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
	
	#scope is ok here
        if options.cmpparms:
            cmpparms = parsedict(options.cmpparms)
        else:
	    cmpparms = {}
	if options.rcmpparms:
            rcmpparms = parsedict(options.rcmpparms)
        else:
	    rcmpparms = {}

        #denoise recons
	if options.famps > 0:
	    tmp = sfixed.do_fft()
	    fth = tmp.get_amplitude_thres(options.famps)
	    tmp.process_inplace('threshold.binary.fourier',{'value':fth})
	    sfixed = tmp.do_ift()
	    tmp = fixed.do_fft()
	    tmp.process_inplace('threshold.binary.fourier',{'value':fth})
	    fixed = tmp.do_ift()
	    
	    tmp = smoving.do_fft()
	    mth = tmp.get_amplitude_thres(options.famps)
	    tmp.process_inplace('threshold.binary.fourier',{'value':fth})
	    smoving = tmp.do_ift()
	    tmp = moving.do_fft()
	    tmp.process_inplace('threshold.binary.fourier',{'value':fth})
	    moving = tmp.do_ift()

        #mask out all the junk
        if options.maskrad > 0:
            sfixed.process_inplace('mask.gaussian.nonuniform', {'radius_x':options.maskrad,'radius_y':options.maskrad,'rae2align3d.pydius_z':options.maskrad,'gauss_width':options.maskfoff})
            smoving.process_inplace('mask.gaussian.nonuniform', {'radius_x':options.maskrad,'radius_y':options.maskrad,'radius_z':options.maskrad,'gauss_width':options.maskfoff})

        #sort of a debugging step
        fixed.set_attr('UCSF.chimera',1)
        fixed.write_image('filtered_fixed.mrc')

	if options.rcmp == "dot":
	    rcmpparms['normalize']=1
	
        #do the global search
        bestscore = 0
        bestmodel = 0
        galignedref = []
<<<<<<< e2align3d.py
        nbest = smoving.xform_align_nbest('rt.3d.sphere', sfixed, {'delta':options.delta,'dotrans':options.dotrans,'dphi':options.dphi,'search':options.search, 'lphi':options.lphi, 'uphi':options.uphi, 'sym':options.sym, 'verbose':options.verbose}, options.nsolns, options.cmp,cmpparms)
=======
        nbest = smoving.xform_align_nbest('rt.3d.sphere', sfixed, {'delta':options.delta,'dotrans':options.dotrans,'dphi':options.dphi,'search':options.search, 'rphi':options.rphi, 'sym':options.sym, 'verbose':options.verbose}, options.nsolns, options.cmp)
	if verbose: print len(nbest), " best orientations"
>>>>>>> 1.25

        #refine each solution found are write out the best one
        for i, n in enumerate(nbest):
            galigned = moving.process('xform',{'transform':n["xform.align3d"]})
            galignedref.append(galigned.align('refine.3d', fixed, {'maxshift':options.maxshift, 'stepalt':options.stepalt,'stepaz':options.stepaz,'stepphi':options.stepphi,'stepx':options.stepx,'stepy':options.stepy,'stepz':options.stepz,'maxiter':options.maxiter}, options.rcmp, rcmpparms))
            score = galignedref[i].get_attr('score')
            if score < bestscore:
                bestscore = score
                bestmodel = i
            if options.verbose > 0: print "Peak Num: ", i, " Transform: ", n["xform.align3d"], " Ini Score: ", n["score"], " Final Score: ", score
	
        # Find out how many peaks are 'the same'
        if options.nsolns > 1:
            thesame = 0
            for m in galignedref:
                if (m.get_attr('score') < bestscore + options.prec and m.get_attr('score') > bestscore - options.prec):
                    thesame += 1
            print str(thesame)+" solns refined to the 'same' point within a precision of "+str(options.prec)
        

        #apply the transform to the original model
        moving.read_image(sys.argv[2])
        ft = galignedref[bestmodel].get_attr("xform.align3d")*nbest[bestmodel]["xform.align3d"] #composition transform
        moving.process_inplace("xform",{"transform":ft})
        
	#now write out the aligned model
        if sys.argv[3]:
	    outfile = sys.argv[3]
	else:
	    outfile = 'alignedmodel.hdf'
	    
	#if output is mrc format 
	if outfile[-4:].lower() == ".mrc":
	    moving.set_attr('UCSF.chimera',1)
	    galignedref[bestmodel].set_attr('UCSF.chimera',1)
  
	galignedref[bestmodel].write_image(("filtered_"+outfile), 0)
	moving.write_image(outfile, 0)
	
if __name__ == "__main__":
    main()
