#!/usr/bin/env python

# 11/20/19  Steven Ludtke
# Script designed to test data compression code timing and amount of compression
from __future__ import print_function
from __future__ import division
from past.utils import old_div
from builtins import range
from EMAN2 import *
import os
import sys
import time

tests=[
 [EM_FLOAT,0,"Float"],
# [EM_USHORT,0,"UShort"],
 [EM_COMPRESSED,0,"Float Cmp"],
 [EM_COMPRESSED,12,"12 bit"],
 [EM_COMPRESSED,9,"9 bit"],
 [EM_COMPRESSED,8,"8 bit"],
 [EM_COMPRESSED,5,"5 bit"],
 [EM_COMPRESSED,4,"4 bit"],
 [EM_COMPRESSED,2,"2 bit"]
]

fls=[(i,i.rsplit(".",1)[0]+".hdf") for i in os.listdir(".") if i.rsplit(".",1)[-1] in ("mrc","mrcs","tif","tiff","hdf")]


print(" "*35,end="")
for i in tests: print("{:16s}".format(i[2]),end="")
print("")

for fi,fo in fls:
	n=EMUtil.get_image_count(fi)
	print("%35s"%fi,end="")
	for ti,t in enumerate(tests):
		fot="test_{}/".format(ti)+fo
		try: os.unlink(fot)
		except: pass
		tm0=time.time()
		for i in range(n):
			#print "  ",i,"  \r",
			sys.stdout.flush()
			im=EMData(fi,i)
			#compressed float
			im["render_bits"]=t[1]
			im.write_image(fot,i,IMAGE_UNKNOWN,0,None,t[0])
		open(fot,"a")	# this should effectively do a 'sync'?
		tm1=time.time()
		if ti==0 :
			base=os.stat(fot).st_size
			print("%7d(%5.2f)  "%(int(base/1024),tm1-tm0),end="")
		else:
			print("%7.1f(%5.2f)  "%(float(base)/os.stat(fot).st_size,tm1-tm0),end="")
	print("")

