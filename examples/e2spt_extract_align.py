from EMAN2 import *
import sys

if len(sys.argv)<2 : 
	print "Usage: e2spt_extract_align.py <bdb:spt_xx#class_ptcl> <type>\ntype=eman,imagic,mrc,spider,quaternion,sgirot,spi,xyz"
	sys.exit(1)

a=Transform({"type":"eman","alt":1.0})
k=list(a.get_rotation(sys.argv[2]).keys())
k.remove("type")
if len(k)==3 : print "#{},{},{}".format(*k)
else: print "#{},{},{},{}".format(k)

n=EMUtil.get_image_count(sys.argv[1])
for i in xrange(n):
	im=EMData(sys.argv[1],i)
	xf=im["spt_ali_param"]
	r=xf.get_rotation(sys.argv[2])
	print "{}".format(i),
	for j in k: 
		print ", {}".format(r[j]),
	print ""
