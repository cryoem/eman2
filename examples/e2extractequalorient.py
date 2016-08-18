#!/usr/bin/env python

# This program will extract an ~equal number of particles over all orientations. 

from EMAN2 import *
from sys import argv
import sys

if len(argv)<4 :
	   print """Usage: e2extractequalorient.py <refine_xx> <iter #> <output> [micrograph restrict]

Extracts particles trying to achieve a reasonably uniform orientation distribution. After launching, 
provides some interactive statistics and asks for a limiting number. It will include at most this 
number of particles from each orientation in the refinement. <ouptut> should be an existing empty 
folder in which to put particles/. [micrograph restrict] is a text file containing the names of
particle image files to include in the output. This is an additional restriction on top of the 
count restriction.
	   """
	   sys.exit(1)

even=EMData.read_images("{}/classes_{:02d}_even.hdf".format(argv[1],int(argv[2])),None,True)
odd =EMData.read_images("{}/classes_{:02d}_odd.hdf".format(argv[1],int(argv[2])),None,True)

ncls=min(len(even),len(odd))
# Number of particles in each orientation
counts=[even[i]["ptcl_repr"]+odd[i]["ptcl_repr"] for i in xrange(ncls)]

# This allows us to restrict the ouput to specific micrographs in addition to the orientation leveling
if len(argv)==5 :
	files=set([base_name(i.strip()) for i in file(argv[4],"r")])
else:
	files=None

outpath=argv[3]+"/"
try:
	os.makedirs(outpath+"particles")
except:
	pass

# Exctract the aggregate list of good particles for each class
ptcls=[[] for i in xrange(ncls)]
for i in xrange(ncls):
	if not even[i].has_attr("class_ptcl_idxs") : lst=[]
	elif isinstance(even[i]["class_ptcl_idxs"],int) : lst=[even[i]["class_ptcl_idxs"]]
	else: lst=even[i]["class_ptcl_idxs"]
	ptcls[i].append(lst)
	
	if not even[i].has_attr("class_ptcl_idxs") : lst=[]
	elif isinstance(odd[i]["class_ptcl_idxs"],int) : lst=[odd[i]["class_ptcl_idxs"]]
	else: lst=odd[i]["class_ptcl_idxs"]
	ptcls[i].append(lst)
	

print "class particle counts range from {} - {} in {} classes".format(min(counts),max(counts),len(counts))
ntk=int(raw_input("How many particles to keep per orientation (at most): "))

# even and odd .lst files referencing original particles
lsx=[LSXFile(even[0]["class_ptcl_src"]),LSXFile( odd[0]["class_ptcl_src"])]

# This is where we actually identify which particles to keep
outfiles={}
for cls in xrange(ncls):
	nk=0
	for eo in (0,1):
		for p in ptcls[cls][eo]:
			if nk>=ntk : break
			orig=lsx[eo][p]			# this becomes number,particle file[,comment]
			if files!=None and base_name(orig[1]) not in files: 
#				print "exclude ",orig
				continue
			try: outfiles[orig[1]].append(orig[0])
			except: outfiles[orig[1]]=[orig[0]]
			nk+=1
			
# copy original particles
for k in outfiles: 
	outfiles[k].sort()
	for i in outfiles[k]:
		origk="particles/{}_ptcls.hdf".format(base_name(k))
		img=EMData(origk,i)
		img.write_image(outpath+origk,-1)
		
		
