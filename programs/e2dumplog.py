#!/usr/bin/env python
# e2dumplog.py   08/08/2005  Steven Ludtke
# This program will dump the local logfile of all EMAN2 programs run in the current
# directory. The file is stored as a python 'shelve' file

import shelve
import sys

db=shelve.open(".eman2log")
try:
	n=int(db["count"])
except:
	print "no logfile"
	sys.exit(0)

for i in range(n-1):
	print " ".join(db[str(i+1)]["args"])
