#!/usr/bin/env python

from EMAN2 import *
from sys import argv
import sys
import threading
import Queue
from Simplex import Simplex
from numpy import array
from time import sleep,time

BOX=256
STEP=256
#NTHREADS=28
VERBOSE=1
DOFSC=0

if len(argv)<3 :
	print """Usage:
	movie_ccf <movie stack> <num threads> [gain norm img]

Will align a movie stack using all-vs-all CCFs with a global optimization strategy. Several outputs
including different frame subsets are produced, as well as a text file with the translation vector map.
"""

if len(argv)>3 :
	normimg=EMData(argv[3])
	print "Normalizing with ",argv[3]
else: normimg=None

NTHREADS=int(argv[2])

data=EMData.read_images(argv[1])
if normimg!=None:
	for i in data: i.mult(normimg)
n=len(data)
nx=data[0]["nx"]
ny=data[0]["ny"]
print "{} frames read {} x {}".format(n,nx,ny)


ccfs=Queue.Queue(0)

# CCF calculation
def calc_ccf(N,box,step,dataa,datab,out):
	for i in range(len(dataa)):
		c=dataa[i].calc_ccf(datab[i],fp_flag.CIRCULANT,True)
		try: csum.add(c)
		except: csum=c

#	csum.process_inplace("normalize.edgemean")
#	csum.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.15})
	out.put((N,csum))

# preprocess regions by normalizing and doing FFT
def split_fft(img,i,box,step,out):
	lst=[]
	for dx in range(box/2,nx-box,step):
		for dy in range(box/2,ny-box,step):
			lst.append(img.get_clip(Region(dx,dy,box,box)).process("normalize.edgemean").do_fft())
	out.put((i,lst))

def calcfsc(map1,map2):
	fsc=map1.calc_fourier_shell_correlation(map2)
	third=len(fsc)/3
	xaxis=fsc[0:third]
	fsc=fsc[third:2*third]

	return(xaxis,fsc)

def qsum(imlist):
	avg=Averagers.get("mean")
	avg.add_image_list(imlist)
	return avg.finish()

# prepare image data by clipping and FFT'ing all tiles
# this is threaded as well
immx=[0]*n
thds=[threading.Thread(target=split_fft,args=(data[i],i,BOX,STEP,ccfs)) for i in range(n)]
print "Precompute FFTs: {} threads".format(len(thds))
t0=time()

thrtolaunch=0
while thrtolaunch<len(thds) or threading.active_count()>1:
	if thrtolaunch<len(thds) :
		while (threading.active_count()==NTHREADS ) : sleep(.1)
#		if VERBOSE : print "Starting thread {}/{}".format(thrtolaunch,len(thds))
		thds[thrtolaunch].start()
		thrtolaunch+=1
	else: sleep(1)

	while not ccfs.empty():
		i,d=ccfs.get()
		immx[i]=d
		
for th in thds: th.join()


# create threads
thds=[]
i=0
for ima in range(n-1):
	for imb in range(ima+1,n):
		thds.append(threading.Thread(target=calc_ccf,args=((ima,imb),BOX,STEP,immx[ima],immx[imb],ccfs)))
		i+=1

print "{:1.1f} s\nCompute ccfs: {} threads".format(time()-t0,len(thds))
t0=time()

# here we run the threads and save the results, no actual alignment done here
csum2={}
thrtolaunch=0
while thrtolaunch<len(thds) or threading.active_count()>1:
	# If we haven't launched all threads yet, then we wait for an empty slot, and launch another
	# note that it's ok that we wait here forever, since there can't be new results if an existing
	# thread hasn't finished.
	if thrtolaunch<len(thds) :
		while (threading.active_count()==NTHREADS ) : sleep(.1)
		#if VERBOSE : print "Starting thread {}/{}".format(thrtolaunch,len(thds))
		thds[thrtolaunch].start()
		thrtolaunch+=1
	else:
		sleep(1)

	while not ccfs.empty():
		i,d=ccfs.get()
		csum2[i]=d

	if VERBOSE: 
		print "  {}/{} {}\r".format(thrtolaunch,len(thds),threading.active_count()),
		sys.stdout.flush()
	
		
for th in thds: th.join()

avgr=Averagers.get("minmax",{"max":0})
avgr.add_image_list(csum2.values())
csum=avgr.finish()
#csum=sum(csum2.values())
#csum.mult(1.0/len(csum2))
#csum.process_inplace("normalize.edgemean")
#display(csum)
#csum.write_image("a.hdf",0)
for i,k in enumerate(sorted(csum2.keys())): 
	im=csum2[k]
#	norm=im[BOX/2,BOX/2]/csum[BOX/2,BOX/2]
#	norm=im.get_clip(Region(BOX/2-5,BOX/2-5,11,11))["mean"]/csum.get_clip(Region(BOX/2-5,BOX/2-5,11,11))["mean"]
#	im.write_image("aa.hdf",i)

# This has been disabled since it eliminates the peak for zero shift. Instead we try the zero/zero elimination hack
	norm=1.0
	im.sub(csum*norm)

	# set the 0,0 peak to the average of neighboring pixels to reduce fixed pattern noise issues (this worked poorly)
	# im[BOX/2,BOX/2]=(im[BOX/2-1,BOX/2]+im[BOX/2+1,BOX/2]+im[BOX/2,BOX/2+1]+im[BOX/2,BOX/2-1])

#	s=im.process("math.sub.optimal",{"ref":csum,"ctfweight":0})

#	im.write_image("a.hdf",i+1)
	# This is critical. Without this, after filtering we get too many false peaks
	thr=im["mean"]+im["sigma"]*1.5
	im.process_inplace("threshold.belowtozero",{"minval":thr})


#####
# Alignment code
#####

# array of x,y locations of each frame, all relative to the last frame in the series, which will always have 0,0 shift
# we store the value for the last frame as well as a conveience
locs=[0]*(n*2)

def qual(locs,ccfs):
	"""computes the quality of the current alignment. Passed a dictionary of CCF images keyed by (i,j) tuple and
	an (x0,y0,x1,y1,...)  shift array. Smaller numbers are better since that's what the simplex does"""

	nrg=0.0
	cen=ccfs[(0,1)]["nx"]/2
	n=len(locs)/2
	for i in xrange(n-1):
		for j in xrange(i+1,n):
#			nrg-=ccfs[(i,j)].sget_value_at_interp(int(cen+locs[j*2]-locs[i*2]),int(cen+locs[j*2+1]-locs[i*2+1]))
			# This is a recognition that we will tend to get better correlation with near neighbors in the sequence
			nrg-=ccfs[(i,j)].sget_value_at_interp(int(cen+locs[j*2]-locs[i*2]),int(cen+locs[j*2+1]-locs[i*2+1]))*sqrt(float(n-fabs(i-j))/n)

#	print nrg
	return nrg
			
#print csum2.keys()

print "{:1.1f} s\nAlignment optimization".format(time()-t0)
t0=time()

# we start with a heavy filter, optimize, then repeat for successively less filtration
for scale in [0.02,0.04,0.07,0.1,0.5]:
	csum3={k:csum2[k].process("filter.lowpass.gauss",{"cutoff_abs":scale}) for k in csum2.keys()}

	incr=[16]*len(locs)
	incr[-1]=incr[-2]=4	# if step is zero for last 2, it gets stuck as an outlier, so we just make the starting step smaller
	simp=Simplex(qual,locs,incr,data=csum3)
	locs=simp.minimize(maxiters=int(100/scale),epsilon=.01)[0]
	locs=[int(floor(i*10+.5))/10.0 for i in locs]
	print locs
	if VERBOSE:
		out=file("path_{:02d}.txt".format(int(1.0/scale)),"w")
		for i in xrange(0,len(locs),2): out.write("%f\t%f\n"%(locs[i],locs[i+1]))
	


# compute the quality of each frame
quals=[0]*n			# quality of each frame based on its correlation peak summed over all images
cen=csum2[(0,1)]["nx"]/2
for i in xrange(n-1):
	for j in xrange(i+1,n):
		val=csum2[(i,j)].sget_value_at_interp(int(cen+locs[j*2]-locs[i*2]),int(cen+locs[j*2+1]-locs[i*2+1]))*sqrt(float(n-fabs(i-j))/n)
		quals[i]+=val
		quals[j]+=val

	
# round for integer only shifting
#locs=[int(floor(i+.5)) for i in locs]

print "{:1.1f}Write unaligned".format(time()-t0)
t0=time()

#write out the unaligned average movie
out=qsum(data)
out.write_image(argv[1].rsplit(".",1)[0]+"_noali.hdf",0)

print "Shift images ({})".format(time()-t0)
t0=time()
#write individual aligned frames
for i,im in enumerate(data):
	im.translate(int(floor(locs[i*2]+.5)),int(floor(locs[i*2+1]+.5)),0)
#	im.write_image("a_all_ali.hdf",i)
out=qsum(data)
out.write_image(argv[1].rsplit(".",1)[0]+"_allali.hdf",0)

#out=sum(data[5:15])	# FSC with the earlier frames instead of whole average
# compute fsc between each aligned frame and the average
# we tile this for better curves, since we don't need the detail
fscq=[0]*n
if DOFSC:
	for i in range(n):
		rgnc=0
		for x in range(64,out["nx"]-192,64):
			for y in range(64,out["ny"]-192,64):
				rgnc+=1.0
				cmpto=out.get_clip(Region(x,y,64,64))
				cscen=data[i].get_clip(Region(x,y,64,64))
				s,f=calcfsc(cmpto,cscen)
				f=array(f)
				try: fs+=f
				except: fs=f
		fs/=rgnc
		fs=list(fs)
		fscq.append(qsum(fs[2:24]))
	
		Util.save_data(s[1],s[1]-s[0],fs[1:-1],argv[1].rsplit(".",1)[0]+"_fsc_{:02d}.txt".format(i))

print "{:1.1f}\nSubsets".format(time()-t0)
t0=time()
# write translations and qualities 
out=open(argv[1].rsplit(".",1)[0]+"_info.txt","w")
out.write("#i,dx,dy,dr,rel dr,qual,(opt)fscqual\n")
for i in range(n):
	out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(i,locs[i*2],locs[i*2+1],hypot(locs[i*2],locs[i*2+1]),hypot(locs[i*2]-locs[i*2-2],locs[i*2+1]-locs[i*2-1]),quals[i],fscq[i]))

thr=max(quals)*0.6	# max correlation cutoff for inclusion
best=[im for i,im in enumerate(data) if quals[i]>thr]
out=qsum(best)
print "Keeping {}/{} frames".format(len(best),len(data))
out.write_image(argv[1].rsplit(".",1)[0]+"_goodali.hdf",0)

thr=max(quals)*0.75	# max correlation cutoff for inclusion
best=[im for i,im in enumerate(data) if quals[i]>thr]
out=qsum(best)
print "Keeping {}/{} frames".format(len(best),len(data))
out.write_image(argv[1].rsplit(".",1)[0]+"_bestali.hdf",0)

# skip the first 4 frames then keep 10
out=qsum(data[4:14])
out.write_image(argv[1].rsplit(".",1)[0]+"_4-14.hdf",0)

# Write out the translated correlation maps for debugging
#cen=csum2[(0,1)]["nx"]/2
#n=len(locs)/2
#for ii,k in enumerate(sorted(csum2.keys())): 
#	i,j=k
#	csum2[k].translate(-int(locs[j*2]-locs[i*2]),-int(locs[j*2+1]-locs[i*2+1]),0)
#	csum2[k].write_image("aa.hdf",ii)

print "{:1.1f}\nDone".format(time()-t0)
