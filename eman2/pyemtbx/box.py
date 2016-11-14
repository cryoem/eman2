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
from bisect import bisect_left,bisect_right
from math import *

# These are particle box sizes which have only 2,3,5 and 7 as prime factors, 
# and are all divisible by 4
good_box_sizes=[4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 48, 56, 60, 64, 72, 80, 84, 96, 100, 108, 112, 120, 128, 140, 
				144, 160, 168, 180, 192, 196, 200, 216, 224, 240, 252, 256, 280, 288, 300, 320, 324, 336, 360, 384, 
				392, 400, 420, 432, 448, 480, 500, 504, 512, 540, 560, 576, 588, 600, 640, 648, 672, 700, 720, 756, 
				768, 784, 800, 840, 864, 896, 900, 960, 972, 980, 1000, 1008, 1024, 1080, 1120, 1152, 1176, 1200, 
				1260, 1280, 1296, 1344, 1372, 1400, 1440, 1500, 1512, 1536, 1568, 1600, 1620, 1680, 1728, 1764, 
				1792, 1800, 1920, 1944, 1960, 2000, 2016, 2048, 2100, 2160, 2240, 2268, 2304, 2352, 2400, 2500, 
				2520, 2560, 2592, 2688, 2700, 2744, 2800, 2880, 2916, 2940, 3000, 3024, 3072, 3136, 3200, 3240, 
				3360, 3456, 3500, 3528, 3584, 3600, 3780, 3840, 3888, 3920, 4000, 4032, 4096, 4116, 4200, 4320, 
				4480, 4500, 4536, 4608, 4704, 4800, 4860, 4900, 5000, 5040, 5120, 5184, 5292, 5376, 5400, 5488, 
				5600, 5760, 5832, 5880, 6000, 6048, 6144, 6272, 6300, 6400, 6480, 6720, 6804, 6860, 6912, 7000, 
				7056, 7168, 7200, 7500, 7560, 7680, 7776, 7840, 8000, 8064, 8100, 8192, 8232, 8400, 8640, 8748, 
				8820, 8960, 9000, 9072, 9216, 9408, 9600, 9604, 9720, 9800, 10000, 10080, 10240, 10368, 10500, 
				10584, 10752, 10800, 10976, 11200, 11340, 11520, 11664, 11760, 12000, 12096, 12288, 12348, 
				12500, 12544, 12600, 12800, 12960, 13440, 13500, 13608, 13720, 13824, 14000, 14112, 14336, 
				14400, 14580, 14700, 15000, 15120, 15360, 15552, 15680, 15876, 16000, 16128, 16200,16384]

def good_boxsize(val,larger=True):
	"This will find the next largest 'good' boxsize"
	if larger:
		if val > 16384: return int(2**ceil(log(val)/log(2.0)))
		return good_box_sizes[bisect_left(good_box_sizes,val)]
	
	if val>16384: return int(2**floor(log(val)/log(2.0)))
	return good_box_sizes[bisect_right(good_box_sizes,val)-1]

    
