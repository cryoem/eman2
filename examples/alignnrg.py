#!/usr/bin/env python

# This computes the local energy surface for 2-D alignment of a particle vs a class-average


from EMAN2 import *
from sys import argv,exit,stdout
from numpy import *

if len(argv)<2 : 
	print "alignnrg <refine#> <it#> <ptcl#> <cmp>"
	exit(1)
argv[1]=int(argv[1])
argv[2]=int(argv[2])
argv[3]=int(argv[3])
comp=parsemodopt(argv[4])

# get the name of the input particles
db=db_open_dict("bdb:refine_%02d#register"%argv[1],True)
ptcl=db["cmd_dict"]["input"]

# class for the particle & orientation parameters
#clmx=EMData.read_images("bdb:refine_%02d#cls_result_%02d"%(argv[1],argv[2]))
clmx=EMData.read_images("bdb:refine_%02d#classify_%02d"%(argv[1],argv[2]))
projn=int(clmx[0][0,argv[3]])	# Projection/class number for this particle
dx=clmx[2][0,argv[3]]
dy=clmx[3][0,argv[3]]
da=clmx[4][0,argv[3]]
df=int(clmx[5][0,argv[3]])

# Read the particle data
im1=EMData(ptcl,argv[3])						# The raw paricle
im1.process_inplace("normalize.edgemean")
im2=EMData("bdb:refine_%02d#projections_%02d"%(argv[1],argv[2]),projn)		# best classified projection
im3=EMData("bdb:refine_%02d#classes_%02d"%(argv[1],argv[2]),projn)			# best classified class-average
im2.process_inplace("normalize.edgemean")
im3.process_inplace("normalize.edgemean")

############
# This section is verifying we have self-consistency among the 2-D alignments from different sources
# This is really just for debugging. They should all be virtually identical !

# get the same alignment info from the similarity matrix for verification
simmx=[]
ncls=EMData("bdb:refine_%02d#simmx_%02d"%(argv[1],argv[2]),0,1)["nx"]
simmx.append(EMData("bdb:refine_%02d#simmx_%02d"%(argv[1],argv[2]),0,0,Region(0,argv[3],ncls,1)))	# sim value
simmx.append(EMData("bdb:refine_%02d#simmx_%02d"%(argv[1],argv[2]),1,0,Region(0,argv[3],ncls,1)))	# tx
simmx.append(EMData("bdb:refine_%02d#simmx_%02d"%(argv[1],argv[2]),2,0,Region(0,argv[3],ncls,1)))	# ty
simmx.append(EMData("bdb:refine_%02d#simmx_%02d"%(argv[1],argv[2]),3,0,Region(0,argv[3],ncls,1)))	# ta
simmx.append(EMData("bdb:refine_%02d#simmx_%02d"%(argv[1],argv[2]),4,0,Region(0,argv[3],ncls,1)))	# flip
dxc=simmx[1][projn,0]
dyc=simmx[2][projn,0]
dac=simmx[3][projn,0]
daf=simmx[4][projn,0]

projnsim=simmx[0].calc_min_index()
print projn,projnsim

# Do a rough alignment ourselves as a 3rd test
im1b=im1.align("rotate_translate_flip",im2)
ali=im1b["xform.align2d"].get_params("2d")

print "These 3 alignments should be roughly the same:"
print "cls_result: %4.1f %4.1f %6.1f %d"%(dx,dy,da,df) 
print "sim_result: %4.1f %4.1f %6.1f %d  %f"%(dxc,dyc,dac,daf,simmx[0][projn,0]) 
print "ali_result: %4.1f %4.1f %6.1f %d"%(ali["tx"],ali["ty"],ali["alpha"],ali["mirror"]) 

###########
# Display the various images as a quick visual check
xfm=Transform({"type":"2d","tx":dx,"ty":dy,"alpha":da,"mirror":df})
im1a=im1.copy()
im1a.transform(xfm)
display((im1,im1a,im1b,im2,im3))

# tst.hdf contains the reference, aligned image, and unaligned image for external testing of aligners
im2["class_id"]=projn
im2.write_image("tst.hdf",0)
im1a.write_image("tst.hdf",1)
im1.write_image("tst.hdf",2)


###########
# Now we compute the quality metric in an exhaustive search around the best alignment
out=EMData(50,50,40)
x=0
for tx in arange(dx-5,dx+5,0.2):
	y=0
	for ty in arange(dy-5,dy+5,0.2):
		a=0
		for ta in arange(da-10,da+10,0.5):
			xfm=Transform({"type":"2d","tx":tx,"ty":ty,"alpha":ta,"mirror":df})
			im1a=im1.copy()
			im1a.transform(xfm)
			sim=im2.cmp(comp[0],im1a,comp[1])
			out[x,y,a]=sim
			a+=1
		y+=1
	x+=1
	print " %d/50\r"%x,
	sys.stdout.flush()

#display(out)
out.write_image("%s_%d.hdf"%(argv[4],argv[3]),0)

# Best alignment from original file
xfm=Transform({"type":"2d","tx":dx,"ty":dy,"alpha":da,"mirror":df})
im1b=im1.copy()
im1b.transform(xfm)

# Best alignment from exhaustive search
x,y,a=out.calc_min_location()
print "\n",x,y,a
x=dx-5+x*0.2
y=dy-5+y*0.2
a=da-10+a*0.5
im1a=im1.copy()
im1a.transform(Transform({"type":"2d","tx":x,"ty":y,"alpha":a,"mirror":df}))

# Redo refine aligner
im1c=im1.align("refine",im2,{"xform.align2d":Transform({"type":"2d","tx":dx+2,"ty":dy-1.5,"alpha":da+4.0,"mirror":df}),"verbose":3},comp[0],comp[1])
newref=im1c["xform.align2d"].get_params("2d")

print "\n"
print "Original: %1.2f\t%1.2f\t%1.1f"%(dx,dy,da)
print "    Best: %1.2f\t%1.2f\t%1.1f"%(x,y,a)
print "Rerefine: %1.2f\t%1.2f\t%1.1f"%(newref["tx"],newref["ty"],newref["alpha"])

display((im2,im1b,im1a,im1c,im1))