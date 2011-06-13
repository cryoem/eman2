#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/22/2010 (sludtke@bcm.edu)
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

# This recomputes the convergence plot FSC curves, which prior to 4/2010 were not being
# computed properly

import sys
from math import *
from os import *

from EMAN2db import db_open_dict
from EMAN2 import *

if len(sys.argv)>1 : dirs=sys.argv[1:]
else : dirs=[i for i in os.listdir(".") if "refine_"==i[:7] and len(i)==9]

# iterates over all refine directories
for d in dirs:
	print "Updating ",d
	
	db2=db_open_dict("bdb:%s#register"%d,True)
	try :
		initm=db2["cmd_dict"]["model"]
		print " Initial model ",initm
	except:
		print "No initial model recorded"
		initm=None

	db=db_open_dict("bdb:%s#convergence.results"%d)

	tdf=[i for i in db_list_dicts("bdb:%s"%d) if "threed_filt" in i or "threed_masked" in i]
	#k=db.keys()

	# instead of using the existing list, we construct our own list from files
	k=[]
	for i in tdf:
		n=int(i.split("_")[2])
		if i=="threed_filt_00" : k.append("init_00_fsc")
		elif "masked" in i : k.append("even_odd_%s_fsc"%i.split("_")[2])
		else : k.append("%02d_%02d_fsc"%(n-1,n))
	print " %d comparisons"%len(k)
	# iterates over all iterations in a directory
	for n in k:
		ns=n.split("_")
		# init vs 01
		if ns[0]=="init" :
			if initm==None : continue
			a=EMData(initm,0)
			b=EMData("bdb:%s#threed_filt_%s"%(d,ns[1]),0)
		# e2refine, ignore
		elif ns[0]=="threed" : continue
		# even/odd tests
		elif ns[0]=="even" :
			a=EMData("bdb:%s#threed_masked_%s_even"%(d,ns[2]),0)
			b=EMData("bdb:%s#threed_masked_%s_odd"%(d,ns[2]),0)
			#try : 
				#m=EMData("bdb:%s#threed_mask_%s"%(d,ns[2]),0)	# We apply the same (soft) mask to both images. Contrary to EMAN1
				#a.mult(m)
				#b.mult(m)
			#except:
				#print "not masking even/odd"
		# general case, convergence curves
		else : 
			try :
				a=EMData("bdb:%s#threed_filt_%s"%(d,ns[0]),0)
				b=EMData("bdb:%s#threed_filt_%s"%(d,ns[1]),0)
			except:
				print "Error with ",n
				sys.exit(1)

		# compute FSC and overwrite original results
		apix=a["apix_x"]
		fsc = a.calc_fourier_shell_correlation(b)
		third = len(fsc)/3
		xaxis = fsc[0:third]
		fsc = fsc[third:2*third]
#		error = fsc[2*third:]
		saxis = [x/apix for x in xaxis]

		db[n]=[saxis,fsc]
		print "  ",n



