#!/usr/bin/env python

# 11/20/19  Steven Ludtke
# Script designed to test HDF5 data compression code timing, amount of compression and information preservation
from builtins import range
from EMAN2 import *
import os
import sys
import time
import numpy as np

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

fls=[(i,i.rsplit(".",1)[0]+".hdf") for i in sorted(os.listdir(".")) if i.rsplit(".",1)[-1] in ("mrc","mrcs","tif","tiff","hdf")]

for i in range(len(tests)):
	try: os.mkdir("test_{}".format(i))
	except: pass
	try: os.mkdir("testfsc_{}".format(i))
	except: pass

print(" "*35,end="")
for i in tests: print("{:^16s}".format(i[2]),end="")
print("")
print(" "*35,"{:^16s}{:^16s}".format("size kb(time)","cmpr ratio(time)"))

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
			
	print("%35s"%"mean,sigma",end="")
	for ti,t in enumerate(tests):
		fot="test_{}/".format(ti)+fo
		fsca=None
		for i in range(n):
			imo=EMData(fi,i)
			imc=EMData(fot,i)
			fsc=imo.calc_fourier_shell_correlation(imc)
			try: fsca+=np.array(fsc)
			except: fsca=np.array(fsc)
			if i==0: print("%7.2f,%5.2f"%(imc["mean"],imc["sigma"]),end="")
		fsca/=n
		write_FSC_file(fsca,"testfsc_{}/".format(ti)+fo.replace(".hdf",".txt"))
			
	print("")

