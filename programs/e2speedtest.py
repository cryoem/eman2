#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
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

import os
import sys
import math
import time
from optparse import OptionParser
from pprint import pprint

from EMAN2 import *

def main():
    progname = os.path.basename(sys.argv[0])
    usage = """%prog [options]
    
This program runs a set of speed tests on the current machine

This program runs a set of speed tests in the current computer. It should
have at least 1 unloaded processor, and minimal i/o occuring when this
is run. It will give a single value which should be generally proportional
to how fast refine will run on a given computer. It is a 'real world' test,
in that it simulates a variety of actual operations performed by a
refinement procedure. Don't compare values given by speed test in
different versions of EMAN, since the underlying routines may be
improved with time."""
    
    parser = OptionParser(usage=usage,version=EMANVERSION)
    
    parser.add_option("--slow",action="store_true",help="rtf_slow alignment",default=False)
    parser.add_option("--best",action="store_true",help="rtf_best alignment",default=False)
    parser.add_option("--low",action="store_true",help="low level test",default=False)
    parser.add_option("--refine",action="store_true",help="refine alignment",default=False)
    parser.add_option("--big",action="store_true",help="big size test",default=False)
    
    (options, args) = parser.parse_args()
    
    SIZE = 96
    NTT = 500
    
    print 'This could take a few minutes. Please be patient.'
    if options.big:
        SIZE = 100
        NIT = 320
        
    print 'Initializing'
    
    pat = EMData();
    pat.set_size(SIZE, SIZE, 1)
    for i in xrange(-SIZE/2, SIZE/2):
        for j in xrange(-SIZE/2, SIZE/2):
            value = -3.0 * math.exp(-Util.square(math.fabs(i)+math.fabs(j)) / 10.0) + math.exp(-Util.square( (math.fabs(i)+math.fabs(j))/2.0 )/100.0) * 2.0 if abs(i)<2 else 1.0  
            pat.set_value_at(i+SIZE/2, j+SIZE/2, value)
    pat.update()
    pat.process_inplace('normalize.circlemean')
    pat.process_inplace("mask.sharp", {"outer_radius":pat.get_xsize()/2.0});
    
    data = [EMData() for i in xrange(5000)]
    
    for i in xrange(NTT):
        data[i] = pat.copy()
        
        for j in xrange(SIZE):
            for k in xrange(SIZE):
                data[i].set_value_at_fast(j, k, Util.get_gauss_rand(0, 1.0))
        data[i].update()
        data[i].process_inplace('normalize.circlemean')
        data[i].process_inplace('mask.sharp', {'outer_radius':data[i].get_xsize()/2});
        
        if i < 5 :
            data[i].write_image('speed.hed', i, EMUtil.ImageType.IMAGE_IMAGIC)
            
    if options.low:
        print 'low'
        
    tmp = EMData()
    t1 = time.clock()
    for i in xrange(3):
        for j in xrange(5, NTT):
            if options.best:
                tmp = data[i].align('rtf_best', data[j], {"flip":None, "maxshift":SIZE/8})
            elif options.slow:
                tmp = data[i].align('rtf_slow', data[j], {"flip":None, "maxshift":SIZE/8})
            elif options.refine:
                tmp = data[i].align('rotate_translate_flip', data[j], {})
            else:
                tmp = data[i].align('rotate_translate_flip', data[j], {});
            
            if j%10 == 0:
                sys.stdout.write('.')
                sys.stdout.flush()
        print
    t2 = time.clock()
    
    ti = t2 - t1;
    
    if not options.big and not options.slow and not options.refine:
        print '\nYour machines speed factor = %1.1f\n' % (25000.0 / ti)
        print '\nThis repesents %1.2f RTFAligns/sec\n' % (3.0 * (NTT - 5.0) / ti)

if __name__ == "__main__":
    main()