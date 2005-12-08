#!/bin/env python
# transalignavg 12/02/2005  Steven Ludtke
# This will read a series of images, translationally align them, average them
# together, and optionally iterate. Translational alignment only.
# transalignavg.py <infile> <dot threshold> <iterations> <background> <gamma>

from EMAN2 import *
import sys
from math import *


n=EMUtil.get_image_count(sys.argv[1])
thr=float(sys.argv[2])

if len(sys.argv)>3: iter=int(sys.argv[3])
else: iter=1

if len(sys.argv)>4 :
	darkref=EMData()
	darkref.read_image(sys.argv[4],0)
else: darkref=None

if len(sys.argv)>5 : gamma=float(sys.argv[5])
else : gamma=0

# Pass 1, sequential alignment to average
avg=EMData()
avg.read_image(sys.argv[1],0)
darkref.process("eman1.normalize.toimage",{"noisy":avg})
avg-=darkref
avg-=avg.get_edge_mean()
#if gamma : avg.process("math.pow",{"pow":gamma})
sum=1
for i in range(1,n):
	a=EMData()
	a.read_image(sys.argv[1],i)
	darkref.process("eman1.normalize.toimage",{"noisy":a})
	a-=darkref
	a-=a.get_edge_mean()
#	if gamma : a.process("math.pow",{"pow":gamma})
	b=a.align("translational",avg,{})
	dot=b.cmp("dot",avg,{"negative":0,"normalize":1})
	print "%4d. %3d\t%3d\t%1.4f"%(i,b.get_attr("translational.dx"),b.get_attr("translational.dy"),dot)
	if dot>thr : 
		avg+=b
		sum+=1
	
print "%d/%d used"%(sum,n)
avg-=avg.get_attr("minimum")
avg/=avg.get_attr("maximum")
avg.process("math.pow",{"pow":gamma})
avg.write_image("avg.mrc")
display(avg)

