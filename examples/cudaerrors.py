#!/usr/bin/env python

#
# Author: David Woolford March 2nd 2009
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


# this is my list of know cuda problems




import EMAN2
from EMAN2 import *
import sys
import math
import os
from time import time

def test_main():
	
	# Problem 1 - FFT reacts in a chaotic manner to images with strane dimensions
	# something weird going on with CUDA: but this is fine:
	a = EMData(599,600)
	b = a.do_fft_cuda()
	c = b.do_ift_cuda()
	# something weird going on with CUDA: but this fails:
	a = EMData(699,700)
	b = a.do_fft_cuda()
	c = b.do_ift_cuda()
	
	# Problem 2 - something goes wrong with a malloc if you create things in this order
	a = EMData(64,64)
	a.mult_cuda(2.0) # this seams to be the problem line
	a = EMData(128,128)
	a.mult_cuda(2.0) #  - alternatively do a._copy_cpu_to_gpu_rw()
	a = EMData(256,256)
	a.mult_cuda(2.0) # ERROR - alternatively do a._copy_cpu_to_gpu_rw()
	
	
	
	
	
if __name__ == '__main__':
	if EMUtil.cuda_available(): test_main()