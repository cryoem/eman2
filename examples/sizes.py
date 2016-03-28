#!/usr/bin/env python

from sys import argv
from EMAN2 import *

try:
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
except: pass

n=EMUtil.get_image_count(argv[1])

peaks=[]
peakvals=[]
out=file("%s.sizes.txt"%argv[1][:-4],"w")
for i in range(n):
	# Read the image
	im=EMData(argv[1],i)
	im.process_inplace("normalize.edgemean")

	# Lowpass filter
	im.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1})

	# Compute radial distribution
	d=im.calc_radial_dist(im.get_xsize()/2,0,1,0)
	d=[ta*tb for ta,tb in enumerate(d)]		# additionalradial weight


	# search for the peak
	mx=(2,-1,-1,-1)
	for x,y in enumerate(d[1:-1]):
		if y>mx[1] : mx=(x+1,y,d[x-1],d[x+1])
		
	
	# Weighted average of the 3 points around the peak. Also 'undoes' the radial weight
	try: peak=(mx[1]+mx[2]+mx[3])/(mx[1]/mx[0]+mx[2]/(mx[0]-1)+mx[3]/(mx[0]+1));
	except: 
		print "error on image %d"%i
		peaks.append(-1.0)
		peakvals.append(-1.0)
		continue

	out.write("%1.2f\n"%peak)
	peaks.append(peak)
	peakvals.append(mx[1]/mx[0])

out.close

# Sort by size
srt=[(j,i) for i,j in enumerate(peaks)]
srt.sort()

# write sorted images
n=0
for i,j in srt:
	im=EMData(argv[1],j)
	im.process_inplace("normalize.edgemean")
	im.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1})
	im["radius"]=i
	im.write_image("%s.srt.hdf"%argv[1][:-4],n)
	n+=1

# This is a plot of peak values vs peak location
plt.cla()
plt.scatter(peaks,peakvals)
plt.xlabel("Radius in pixels")
plt.ylabel("Peak height")
plt.savefig("%s.peakval.png"%argv[1][:-4])

# This plots a histogram of the sizes
plt.cla()
h=plt.hist(peaks,bins=im.get_xsize()/2,range=(0,im.get_xsize()/2))
avg=sum(peaks)/len(peaks)
plt.text(2,max(h[0])-5,"Mean size = %1.1f pixels"%avg)
plt.xlabel("Radius in pixels")
plt.savefig("%s.hist.png"%argv[1][:-4])

