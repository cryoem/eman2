#!/usr/bin/env python

import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
pltcolors=["k","b","g","r","m","c","darkblue","darkgreen","darkred","darkmagenta","darkcyan","0.5"]

from EMAN2 import *
from math import *
import os
import sys
import numpy as np

d=np.loadtxt(sys.argv[1]).transpose()
maxs=max(d[0])

xticklocs=[i for i in (.01,.05,.0833,.125,.1667,.2,.25,.3333,.4,.5) if i<maxs]
xticklbl=["1/100","1/20","1/12","1/8","1/6","1/5","1/4","1/3","1/2.5","1/2"][:len(xticklocs)]
yticklocs=(0.0,.125,.143,.25,.375,.5,.625,.75,.875,1.0)
yticklbl=("0"," ","0.143","0.25"," ","0.5"," ","0.75"," ","1.0")
yticklocs2=(0.0,.125,.25,.375,.5,.625,.75,.875,1.0)
yticklbl2=("0"," ","0.25"," ","0.5"," ","0.75"," ","1.0")


try:
#	plt.title("Gold Standard Resolution (tight mask)")
	plt.xlabel(r"Spatial Frequency (1/$\AA$)")
	plt.ylabel("FSC")
except:
	pass

#fscs=[i for i in os.listdir(options.path) if "fsc_maskedtight" in i and i[-4:]==".txt"]
#fscs.sort(reverse=True)
fscs=sys.argv[1:]
#nummx=int(fscs[0].split("_")[2][:2])
maxx=0.01

# iterate over fsc curves
for i,f in enumerate(fscs):
#	num=int(f.split("_")[2][:2])

	# read the fsc curve
	d=np.loadtxt(f).transpose()

	# plot the curve
	try: plt.plot(d[0],d[1],label=f.split("/")[-1],color=pltcolors[(i)%12])
#	try: plt.plot(d[0],d[1],color=pltcolors[(i)%12])
	except: pass
	maxx=max(maxx,max(d[0]))

	# find the resolution from the first curve (the highest numbered one)
	if f==fscs[0]:
#		lastnum=num
		# find the 0.143 crossing
		for si in xrange(2,len(d[0])-2):
			if d[1][si-1]>0.143 and d[1][si]<=0.143 :
				frac=(0.143-d[1][si])/(d[1][si-1]-d[1][si])		# 1.0 if 0.143 at si-1, 0.0 if .143 at si
				lastres=d[0][si]*(1.0-frac)+d[0][si-1]*frac
				try:
					plt.annotate(r"{:1.1f} $\AA$".format(1.0/lastres),xy=(lastres,0.143),
						xytext=((lastres*4+d[0][-1])/5.0,0.2),arrowprops={"width":1,"frac":.1,"headwidth":7,"shrink":.05})
				except: pass
				break
		else : lastres=0


plt.axhline(0.0,color="k")
plt.axhline(0.143,color="#306030",linestyle=":")
plt.axis((0,maxx,-.02,1.02))
#plt.legend(loc="lower left")
plt.legend(loc="upper right",fontsize="x-small")
#plt.minorticks_on()
plt.xticks(xticklocs,xticklbl)
plt.yticks(yticklocs,yticklbl)
plt.savefig("plot.png")
try: plt.savefig("plot.pdf")
except: pass
