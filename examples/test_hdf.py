#!/usr/bin/env python

# 11/20/19  Steven Ludtke
# Script designed to test HDF5 data compression code timing, amount of compression and information preservation
from builtins import range
from EMAN2 import *
import os
import sys
import time
import numpy as np
from traceback import print_exc

tests=[
	[EM_FLOAT,0],
	[EM_CHAR,7],
	[EM_UCHAR,8],
	[EM_SHORT,15],
	[EM_USHORT,16],
	[EM_COMPRESSED,0],
	[EM_COMPRESSED,4],
	[EM_COMPRESSED,12]
	]

for a,ab in tests:
	for b,bb in tests:
		print(f"{file_mode_imap[a]:11s} {ab}\t{file_mode_imap[b]:11s} {bb}",end="\t")
		im=test_image()
		# these will be overwritten by write_compressed
		im["render_min"]=im["minimum"]
		im["render_max"]=im["maximum"]
		im["render_bits"]=ab
		try: os.unlink("/dev/shm/a.hdf")
		except: pass
		if a==EM_COMPRESSED: 
			im.write_compressed("/dev/shm/a.hdf",0,ab)
		else: im.write_image("/dev/shm/a.hdf",0,IMAGE_HDF,False,None,a)
		
		im2=EMData("/dev/shm/a.hdf",0)
		print(f"{-im.cmp('ccc',im2):1.5f}",end="\t")
		
		# now rewrite to the same file (different image)
		im=test_image(1)
		# these will be overwritten by write_compressed
		im["render_min"]=im["minimum"]
		im["render_max"]=im["maximum"]
		im["render_bits"]=bb
		try:
			if b==EM_COMPRESSED: 
				im.write_compressed("/dev/shm/a.hdf",0,bb)
			else: im.write_image("/dev/shm/a.hdf",0,IMAGE_HDF,False,None,b)
			print(" ",end="\t")
		except:
			print_exc()
			print("ERR",end="\t")
			os.unlink("/dev/shm/a.hdf")
			if b==EM_COMPRESSED: im.write_compressed("/dev/shm/a.hdf",0,bb)
			else: im.write_image("/dev/shm/a.hdf",0,IMAGE_HDF,False,None,b)
		
		im3=EMData("/dev/shm/a.hdf",0)
		print(f"{-im.cmp('ccc',im3):1.5f}",end="\t")
		print(f"{-im2.cmp('ccc',im3):1.5f}")
		
		
im4=test_image_3d()
im4.write_image("/dev/shm/a.hdf",0)
im4.write_compressed("/dev/shm/a.hdf",0)
					 
