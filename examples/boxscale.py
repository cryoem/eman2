#!/usr/bin/env python
# This simple script will upscale .box files
# assumes names contain "shrink" which is stripped on output


import sys

scale=4.0

for f in sys.argv[1:]:
	if "shrink" not in f :
		print "filenames must contain 'shrink'"
		sys.exit(1)
	outf=file(f.replace("shrink",""),"w")
	for l in file(f,"r"):
		x,y,s,s2=[int(i) for i in l.split()]
		if x<0 or y<0 : continue		# remove any out of bound boxes on the side where we can tell
		outf.write("%d\t%d\t%d\t%d\n"%(x*scale,y*scale,s*scale,s2*scale))
	outf.close()
