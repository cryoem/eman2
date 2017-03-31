#!/usr/bin/env python

# This program is designed to take an aligned movie stack and extract the best quality frames and average them together. This is
# similar to capabilities in e2ddd_movie, but without the alignment step.  This program was designed in the context of a specific project
# and filenaming,etc would have to be tweaked for use in other projects

from EMAN2 import *
import sys
import os
import numpy as np

# read MRC stack even when named ".mrc"
print "read"
try: os.unlink("tmp.mrcs")
except: pass 
os.symlink(sys.argv[1],"tmp.mrcs")
stk=EMData.read_images("tmp.mrcs")
info=js_open_dict(info_name(sys.argv[1].split("-")[0]+".mrc"))
ctf=info["ctf_frame"][1]

# CTF amplitude filter for (hopefully better) quality assessment
# phase flipping shouldn't actually matter, but downweighting values near zero could be important
filter=EMData(1024*3+2,1024*3,1)
filter.set_complex(True)
filter["apix_x"]=1.26
filter["apix_y"]=1.26
filter["apix_z"]=1.26
ctf.compute_2d_complex(filter,Ctf.CtfType.CTF_AMP)
#ctf.compute_2d_complex(flipim,Ctf.CtfType.CTF_SIGN)

# we work only with a 3k x 3k region from the middle of the image
print "preprocess"
stkf=[]
for img in stk:
	imgf=img.get_clip(Region(512,512,1024*3,1024*3)).process("normalize.edgemean").do_fft()
	imgf.mult(filter)
	stkf.append(imgf.do_ift())

avg=sum(stkf)
#display(filter)
#display(avg)
#display(stkf,True)

#rng=[.01,.06,.1,.25]
rng=[80,18,10,4]

# write the quality plot for different resolutions (80-18, 18-10 and 10-4 A)
out=open(sys.argv[1].split("-")[0]+"_qual.txt","w")
qlist=[]
for i in xrange(len(stkf)):
	im=stkf[i]
	out.write("{}".format(i))
	qlist.append([i])
	for f in (0,1,2):
#		im2=im.process("filter.lowpass.tophat",{"cutoff_freq":rng[f+1]}).process("filter.highpass.tophat",{"cutoff_freq":rng[f]})
		c=-im.cmp("frc",avg,{"minres":rng[f],"maxres":rng[f+1]})
		qlist[-1].append(c)
		out.write("\t{}".format(c))

	out.write("\n")

# extract the best frames
stat=np.array(qlist)
stathires=stat[:,3]		# 3rd column
mean=stat.transpose()[3].mean()
sigma=stat.transpose()[3].std()

statgood=stat[stathires>mean-sigma*.5]		# extract the indices of the particles where the high resolution quality is > the average - sigma

print statgood[:,0], len(statgood),len(stat)
print "Write Output"

avr=Averagers.get("mean")
avr.add_image_list(stk)
av=avr.finish()
av["class_ptcl_idxs"]=range(len(stk))
av.write_image("micrographs/"+base_name(sys.argv[1].split("-")[0]+".mrc")+"__ali.hdf",0)

avr=Averagers.get("mean")
avr.add_image_list(stk[2:19])
av=avr.finish()
av["class_ptcl_idxs"]=range(len(2,19))
av.write_image("micrographs/"+base_name(sys.argv[1].split("-")[0]+".mrc")+"__218.hdf",0)


avr=Averagers.get("mean")
lst=[]
for i in statgood[:,0]:
	avr.add_image(stk[int(i)])
	lst.append(int(i))
av=avr.finish()
av["class_ptcl_idxs"]=lst
av.write_image("micrographs/"+base_name(sys.argv[1].split("-")[0]+".mrc")+"__goodali.hdf",0)

