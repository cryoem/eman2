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

# makebigfromseq 12/25/2005  Steven Ludtke
# This will read a series of images, translationally align them, average them
# together, and optionally iterate. Translational alignment only.
# makebigfromseq.py <infile> <sizexsize> <dot threshold> <darkref> <lightref> <1st image shift>

from EMAN2 import *
import sys
from math import *
try:
	import wx
	from emimage import *
except:
	print "no wx"

def fixup(i,d,l):
	"""Takes an image, a dark reference and a light/flat reference and corrects the image in-place"""
#	i2=i.copy()
	if d:
		d.process_inplace("normalize.toimage",{"to":i})
		i-=d
	i/=l
	i.process_inplace("normalize.edgemean",{})
#	display([i2,i,d,l])

def transalign(i1,i2,m1,m2):
	c1=i1.calc_ccf(i2)
	c2=m1.calc_ccf(m2)
	c2.process_inplace("threshold.belowtozero",{"value":c2.get_attr("sigma")/5.0})
	c2+=c2.get_attr("sigma")/5.0
#	display(c2)
	c1/=c2
	l=list(c1.calc_max_location())
	v=c1.get_value_at(l[0],l[1],0)
	l[0]=i1.get_xsize()/2-l[0]
	l[1]=i1.get_ysize()/2-l[1]
#	print l,v
	
	i=i1.copy()
	i.translate(l[0],l[1],0)
	i.set_attr("translational.dx",l[0])
	i.set_attr("translational.dy",l[1])
	i.set_attr("align_score",v)
	return i
	
def main(argv,app=None) :
	if len(argv)==1 : 
		print "makebigfromseq.py <infile> <sizexsize> <dot threshold> <darkref> <lightref> <im1 shift x,y>"
		sys.exit(1)
		
	n=EMUtil.get_image_count(argv[1])
	try: sz=[int(i) for i in argv[2].split("x")]
	except: sz=[int(argv[2])]*2
	
	try:
		thr=float(argv[3])
		thr2=thr
	except:
		thr=argv[3].split(',')
		thr2=float(thr[1])
		thr=float(thr[0])
	
	if (len(argv)>4) : 
		dark=EMData()
		try: dark.read_image(argv[4])
		except: dark=None
	else: dark=None
	
	if (len(argv)>5) : 
		light=EMData()
		light.read_image(argv[5])
	
	try:
		dx0,dy0=argv[6].split(',')
		dx0=int(dx0)
		dy0=int(dy0)
	except:
		dx0,dy0=0,0
		
	avg=EMData()
	avg.read_image(argv[1],0)
	mask=avg.copy()
	fixup(avg,dark,light)
	avg=avg.get_clip(Region(-(sz[0]-avg.get_xsize())/2-dx0,-(sz[1]-avg.get_ysize())/2-dy0,sz[0],sz[1]))
	avg2=avg.copy()
	
	mask.to_one()
	mask=mask.get_clip(Region(-(sz[0]-mask.get_xsize())/2,-(sz[1]-mask.get_ysize())/2,sz[0],sz[1]))
	mask+=.0000000001
	avgn=mask.copy()			# normalization image, start with 1 for the first image, always included
	avgn.translate(dx0,dy0,0)

	if (app) :
		ai=EMImage(avg)
	
	for i in range(1,n):
		a=EMData()
		a.read_image(argv[1],i)
		fixup(a,dark,light)
		a=a.get_clip(Region(-(sz[0]-a.get_xsize())/2,-(sz[1]-a.get_ysize())/2,sz[0],sz[1]))
		if i%25==0: avg.write_image("avg.mrc")
	#	b=a.align("translational",a2,{"maxshift":sz[0]/2})
		b=transalign(a,avg,mask,avgn)
		m2=mask.copy()			# copy and translate the map for normalization map
		m2.translate(b.get_attr("translational.dx"),b.get_attr("translational.dy"),0)
	#	olap=an.cmp("dot",m2,{"negative":0,"normalize":1})		# overlap dot product for normalization of quality factors
	#	if olap==0 : continue
	#	olap=sqrt(olap)		# we don't want to upweight small overlaps too much, this is empirical
		print i,b.get_attr("align_score"),b.get_attr("translational.dx"),b.get_attr("translational.dy"),
		if b.get_attr("align_score")<thr : 
			print " *"
			continue
		else : print ""
	#	b.process_inplace("normalize.toimage",{"noisy":a2,"keepzero":1})
		avg+=b
		avgn+=m2
		if (app) : app.Yield()
		
	avg/=avgn
	avg.write_image("avg.mrc")
	avgn.write_image("avgn.mrc")
	
	avgn2=mask.copy()			# normalization image, start with 1 for the first image, always included
	avgn2.translate(dx0,dy0,0)
	avgn2+=.0000000001
	
	avg*=avgn
	for i in range(1,n):
		a=EMData()
		a.read_image(argv[1],i)
		fixup(a,dark,light)
		a=a.get_clip(Region(-(sz[0]-a.get_xsize())/2,-(sz[1]-a.get_ysize())/2,sz[0],sz[1]))
		if i%25==0: avg2.write_image("avg2i.mrc")
	#	b=a.align("translational",avg,{"maxshift":sz[0]/2})
		b=transalign(a,avg,mask,avgn)
		m2=mask.copy()			# copy and translate the map for normalization map
		m2.translate(b.get_attr("translational.dx"),b.get_attr("translational.dy"),0)
		print i,b.get_attr("align_score"),b.get_attr("translational.dx"),b.get_attr("translational.dy"),
		if b.get_attr("align_score")<thr2: 
			print " *"
			continue
		else : print ""
	#	b.process_inplace("normalize.toimage",{"noisy":avg,"keepzero":1})
		avg2+=b
		avgn2+=m2
	avg2/=avgn2
	avg2.write_image("avg2.mrc")
	avgn.write_image("avg2n.mrc")

if __name__ == "__main__":
   	main(sys.argv)
