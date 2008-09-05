#!/bin/env python

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

from EMAN2  import *
from sparx  import *

    
N = 2000
K = 50
d = model_blank(N*(N-1)/2)
x = []
from random import random, gauss
from math import sqrt
Kt = 3 # actual number of groups
for i in xrange(N):
	if(i%Kt == 0):  x.append(gauss(0.0,1.0))
	elif(i%Kt == 1):  x.append(gauss(1.0e-1,1.0))
	elif(i%Kt == 2):  x.append(gauss(2.0e-1,1.0))


print  ttime()
for j in xrange(1,N):
	for i in xrange(j):
		d.set_value_at(mono(i,j),sqrt((x[i]-x[j])**2))
dmin = 1.0e23
print  ttime()
for i in xrange(100):
	o = Util.cluster_pairwise(d,K)
	if(dmin > o[N+K]):
		print  i,ttime(),o[N+K:N+K+2]
		dmin = o[N+K]
		best = o
print ttime()
#print  best[0:N]  # assignments
print   best[N:N+K]  # objects that are centers
print   best[N+K:N+K+2]  # criterion (should be minimized) and number of tierations)
