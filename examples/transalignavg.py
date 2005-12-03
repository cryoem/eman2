#!/bin/env python
# transalignavg 12/02/2005  Steven Ludtke
# This will read a series of images, translationally align them, average them
# together, and optionally iterate. Translational alignment only.
# transalignavg.py <infile> <dot threshold>

from EMAN2 import *
import sys
from math import *


n=EMUtil.get_image_count(sys.argv[1])
thr=float(sys.argv[2])

# Pass 1, sequential alignment to average
avg=EMData()
avg.read_image(sys.argv[1],0)
avg-=avg.get_edge_mean()
sum=1
for i in range(1,n):
	a=EMData()
	a.read_image(sys.argv[1],i)
	a-=a.get_edge_mean()
	b=a.align("translational",avg,{})
	dot=b.cmp("dot",avg,{"negative":0,"normalize":1})
	print "%4d. %3d\t%3d\t%1.4f"%(i,b.get_attr("translational.dx"),b.get_attr("translational.dy"),dot)
	if dot>thr : 
		avg+=b
		sum+=1
	
print "%d/%d used"%(sum,n)
display(avg)

