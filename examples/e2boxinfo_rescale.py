#!/usr/bin/env python
from EMAN2 import *
from sys import argv,exit

if len(argv)<2 : 
	print """Usage: e2boxinfo_rescale.py <scale factor>

WARNING: this will rescale the box locations in ALL info/*info.json files."""
	sys.exit(1)

sfac=float(argv[1])
print "Rescaling by %1.2fx"%sfac

il=js_list_dicts("info")

tbox=0
tdb=0
for i in il:
	db=js_open_dict("info/{}".format(i))
	if db.has_key("boxes") :
		tdb+=1
		boxes=db["boxes"]
		for j in xrange(len(boxes)):
			boxes[j][0]*=sfac
			boxes[j][1]*=sfac
		db["boxes"]=boxes
		
		print "{} with {} boxes".format(i,len(boxes))
		tbox+=len(boxes)

print "Done. {} micrographs with {} total boxes scaled by {}".format(tdb,tbox,sfac)
