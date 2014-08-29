#!/usr/bin/env python
# This program will find the user-assigned quality for all of the images in "micrographs"
# and delete any with a quality less than the specified (command-line) number

from EMAN2 import *
from sys import argv
import os

minqual=int(argv[1])

imgs=["micrographs/"+i for i in os.listdir("micrographs") if i[0]!="." and ".hed" not in i ]
imgs.sort()

badimgs=[]
for im in imgs:
	try:
		j=js_open_dict(info_name(im))
		if j["quality"]<minqual : badimgs.append(im)
	except: pass

for im in badimgs:
	print im

print len(badimgs)," identified. Are you sure you want to delete (y/n)? ",
ans=raw_input()

if len(ans)==0 or ans[0].lower()=="y" :
	for im in badimgs:
		os.unlink(im)
	print len(badimgs)," deleted"

else: print "aborted"
