#!/usr/bin/env python
# e2boxer.py  07/27/2004  Steven Ludtke
# This program is used to box out particles from micrographs/CCD frames

from EMAN2 import *
from optparse import OptionParser
from math import *
import time
import os
import sys

pl=()

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <image>
	
Automatic and manual particle selection. This version is specifically aimed at square boxes
for single particle analysis."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--gui",action="store_true",help="Start the GUI for interactive boxing",default=False)
	parser.add_option("--box","-B",type="int",help="Box size in pixels",default=-1)
	parser.add_option("--ptclsize","-P",type="int",help="Approximate size (diameter) of the particle in pixels. Not required if reference particles are provided.",default=-1)
	parser.add_option("--refptcl","-R",type="string",help="A stack of reference images. Must have the same scale as the image being boxed.",default=None)
	parser.add_option("--refvol","-V",type="string",help="A 3D model to use as a reference for autoboxing",default=None)
	parser.add_option("--sym","-S",type="string",help="Symmetry of the 3D model",default=None)
	parser.add_option("--auto","-A",type="string",action="append",help="Autobox using specified method: circle, ref",default=[])
	parser.add_option("--grid","-G",type="int",help="Box the entire image in a grid pattern with the specified overlap in pixels (can be negative)",default=None)
	parser.add_option("--threshold","-T",type="float",help="Threshold for keeping particles. 0-4, 0 excludes all, 4 keeps all.",default=2.0)
	parser.add_option("--nretest",type="int",help="Number of reference images (starting with the first) to use in the final test for particle quality.",default=-1)
	parser.add_option("--retestlist",type="string",help="Comma separated list of image numbers for retest cycle",default="")
	parser.add_option("--norm",action="store_true",help="Edgenormalize boxed particles",default=False)
	parser.add_option("--savealiref",action="store_true",help="Stores intermediate aligned particle images in boxali.hdf. Mainly for debugging.",default=False)
	parser.add_option("--farfocus",type="string",help="filename or 'next', name of an aligned far from focus image for preliminary boxing",default=None)
	parser.add_option("--dbout",type="string",help="filename to write EMAN1 style box database file to",default=None)
	parser.add_option("--ptclout",type="string",help="filename to write boxed out particles to",default=None)
	
	(options, args) = parser.parse_args()
	if len(args)<1 : parser.error("Input image required")

	logid=E2init(sys.argv)
	
	if (options.farfocus==None): initial=args[0]
	elif (options.farfocus=="next") : 
		initial=str(int(args[0].split('.')[0])+1)+'.'+args[0].split('.')[1]
		if not os.access(initial,os.R_OK) :
			print "No such file: ",initial
			sys.exit(1)
		print "Using initial file: ",initial
	else: initial=options.farfocus
	
	image=EMData()
	image.read_image(initial)
	image.set_attr("datatype",7)
	
	if options.grid!=None :
		step=int(options.box)+int(options.grid)
		print "Grid boxing with box size %d and step %d"%(options.box,step)
		x0=image.get_xsize()%step
		y0=image.get_ysize()%step
		outn=args[0][:-3]+"grid.hed"
		for y in range(y0,image.get_ysize()-step,step):
			for x in range(x0,image.get_xsize()-step,step):
				nim=image.get_clip(Region(x,y,options.box,options.box))
				nim.write_image(outn,-1)
		sys.exit(0)
				
	refptcl=None
	if options.refptcl :
		print options.refptcl
		refptcl=EMData.read_images(options.refptcl)
		refbox=refptcl[0].get_xsize()
		print "%d reference particles read (%d x %d)"%(len(refptcl),refbox,refbox)
	
	if options.box<5 :
		if options.refptcl : options.box=refptcl[0].get_xsize()
		elif options.ptclsize : 
			options.box=good_boxsize(options.ptclsize*1.2)
		else : parser.error("Please specify a box size")
	else:
		if not options.box in good_box_sizes:
			print "Note: EMAN2 processing would be more efficient with a boxsize of %d"%good_boxsize(options.box)
	
	try: options.retestlist=[int(i) for i in options.retestlist.split(',')]
	except: pass
			
	shrinkfactor=int(ceil(options.box/16))
	print "Shrink factor = ",shrinkfactor
	#shrinkfactor=int(ceil(image.get_ysize()/1024.0))
	#if options.box/shrinkfactor<12 : shrinkfactor/=2
	
	image.process_inplace("eman1.normalize")
	shrink=image
	shrink.mean_shrink(shrinkfactor)		# shrunken original image
	
	# This confusing line insures the shrunken images have even dimensions
	# since odd FFTs haven't been fixed yet
	if shrink.get_xsize()&1 or shrink.get_ysize()&1 :
		shrink=shrink.get_clip(Region(0,0,(shrink.get_xsize()|1)^1,(shrink.get_ysize()|1)^1))

	# now we try to clean up long range density variations in the image
	filtrad=options.box*2/shrinkfactor
	flt=EMData()
	flt.set_size(shrink.get_xsize()+filtrad*4,shrink.get_ysize()+filtrad*4,1)
	flt.to_one()
	flt.process_inplace("eman1.mask.sharp",{"outer_radius":filtrad})
	flt/=(float(flt.get_attr("mean"))*flt.get_xsize()*flt.get_ysize())
	flt.process_inplace("eman1.xform.phaseorigin")
	a=shrink.get_clip(Region(-filtrad*2,-filtrad*2,shrink.get_xsize()+filtrad*4,shrink.get_ysize()+filtrad*4))
	a.process_inplace("eman1.mask.zeroedgefill")
#	a.write_image("q0.hdf",0)
#	flt.write_image("q0.hdf",1)
	a=a.convolute(flt)
#	a.write_image("q0.hdf",2)
	a=a.get_clip(Region(filtrad*2,filtrad*2,shrink.get_xsize(),shrink.get_ysize()))
#	shrink.write_image("q1.hdf",0)
#	a.write_image("q1.hdf",1)
	shrink-=a
#	shrink.write_image("q1.hdf",2)
	a=None
	
	shrink2=shrink.copy()
	shrink2.process_inplace("eman1.math.squared")
#	image=EMData()
#	image.read_image(args[0])
#	shrink.write_image("e.mrc")
#	shrink2.write_image("f.mrc")
		
	print "Autobox mode ",options.auto[0]
	
	if "ref" in options.auto:
		if not refptcl: error_exit("Reference particles required")
		
		print "Prepare references"
		# refptcls will contain shrunken normalized reference particles
		refptcls=[]
		for n,i in enumerate(refptcl):
			# first a circular mask
			i.process_inplace("eman1.normalize.circlemean")
			i.process_inplace("eman1.mask.sharp",{"outer_radius":i.get_xsize()/2-1})
						
			ic=i.copy()
			refptcls.append(ic)
			ic.mean_shrink(shrinkfactor)
			# make the unmasked portion mean -> 0
			ic-=float(ic.get_attr("mean_nonzero"))
			ic.process_inplace("eman1.normalize.unitlen")
#			ic.write_image("scaled_refs.hdf",-1)

		# prepare a mask to use for local sigma calculaton
		circle=shrink.copy_head()
		circle.to_one()
		circle.process_inplace("eman1.mask.sharp",{"outer_radius":options.box/(shrinkfactor*2)-1})
		circle/=(float(circle.get_attr("mean"))*circle.get_xsize()*circle.get_ysize())
		
		ccfmean=shrink.calc_ccf(circle,fp_flag.CIRCULANT)
#		circle.write_image("z0a.hdf")
#		shrink2.write_image("z0b.hdf")
		ccfsig=shrink2.calc_ccf(circle,fp_flag.CIRCULANT)
		ccfmean.process_inplace("eman1.math.squared")
		ccfsig-=ccfmean		# ccfsig is now pointwise standard deviation of local mean
		ccfsig.process_inplace("eman1.math.sqrt")
#		shrink.write_image("z0.hdf")
#		ccfsig.write_image("z1.hdf")
		
		print "Locating possible particles"
		xs=shrink.get_xsize()
		ys=shrink.get_ysize()
		pks=[]
		for n,i in enumerate(refptcls):
			j=i.get_clip(Region(-(xs-i.get_xsize())/2,-(ys-i.get_ysize())/2,xs,ys))
#			j.write_image("0.%0d.hdf"%n)
			ccfone=shrink.calc_ccf(j,fp_flag.CIRCULANT)
#			ccfone.write_image("a.%0d.hdf"%n)
			ccfone/=ccfsig
#			ccfone.write_image("b.%0d.hdf"%n)
			sig=float(ccfone.get_attr("sigma"))
			ccfone.process_inplace("eman1.mask.onlypeaks",{"npeaks":0})
#			ccfone.write_image("c.%0d.hdf"%n)
			pk=ccfone.calc_highest_locations(sig*3.5)
			for m,p in enumerate(pk):
				pk[m]=(-p.value,n,p.x,p.y)
			pks+=pk
			
		pks.sort()		# an ordered list of the best particle locations
		print "%d putative particles located"%len(pks)
				
		# this will produce a new list excluding any lower valued boxes within
		# 1/2 a box size of a higher one. It also rescales the boxes.
		# (ok, you could do this with clever syntax, but this is more readable)

		# first we prepare a grid table of putative box locations for speed
		grid={}
		for n,i in enumerate(pks):
			x=int(floor(i[2]*shrinkfactor/options.box))
			y=int(floor(i[3]*shrinkfactor/options.box))
			try: grid[(x,y)].append(n)
			except: grid[(x,y)]=[n]

		
		goodpks=[]
		bf=options.box/(shrinkfactor*2)
		for n,i in enumerate(pks):
			if i[2]<bf or i[3]<bf or i[2]>xs-bf-1 or i[3]>ys-bf-1 : continue

			# local is a list of putative peaks near the current peak
                        x=int(floor(i[2]*shrinkfactor/options.box))
                        y=int(floor(i[3]*shrinkfactor/options.box))
			local=[]
			for xx in range(x-1,x+2):
				for yy in range(y-1,y+2):
					try: local+=grid[(xx,yy)]
					except: pass
			local=filter(lambda x: x>n,local)

			for nn in local:
				ii=pks[nn]
				if hypot(i[2]-ii[2],i[3]-ii[3])<bf*3/2 : break
			else: goodpks.append([i[0],i[1],i[2]*shrinkfactor-options.box/2,i[3]*shrinkfactor-options.box/2])
		
		print "%d putative particles after local exclusion"%len(goodpks)
		
		print "refine particle locations"
		# This will optimize the center location of each particle and improve
		# the similarity calculation
		goodpks2=[]
		n=0
		for i in goodpks:
			b=EMData()
			
			# on the first pass, rather than just using the best reference determined using the shrunken images
			# we also try a subset of the references to see which gives the best alignment/match
			# refns is the list of ref image #'s to test against
			if options.nretest>0 :
				refns=[i[1]]+range(options.nretest)	# list of the ref numbers of the particles to use in the realignment
			else :
				refns=[i[1]]+options.retestlist
			
			# read in the area where we think a particle exists
			try: b.read_image(initial,0,0,Region(i[2],i[3],options.box,options.box))
			except: continue
			b.process_inplace("eman1.normalize.edgemean")
			b.process_inplace("eman1.filter.lowpass.gaussian",{"lowpass":.1})
#			ba=refptcl[i[1]].align("rotate_translate",b,{},"SqEuclidean")
			
			# we iterate over each reference then find the best, eventually recentering ( i[2-3] ) and
			# updating the reference number for the second pass ( i[1] )
			tsts=[]
			for j in refns: 
				ba=b.align("rotate_translate",refptcl[j],{},"optvariance",{"matchfilt":1})
				tsts.append([ba.get_attr("align_score"),j,ba.get_attr("translational.dx"),ba.get_attr("translational.dy"),ba.get_attr("rotational")])
			tsts.sort()
#			if tsts[0][1]!=i[1] : print i[1]," -> ",tsts[0][1],"    %f,%f  %f"%(tsts[0][2],tsts[0][3],tsts[0][4])
			i[1]=tsts[0][1]
			i[2]-= cos(tsts[0][4])*tsts[0][2]+sin(tsts[0][4])*tsts[0][3]
			i[3]-=-sin(tsts[0][4])*tsts[0][2]+cos(tsts[0][4])*tsts[0][3]

# this code can be used to test alignment
#			b.write_image("cmp.hdf",-1)
#			refptcl[i[1]].write_image("cmp.hdf",-1)
#			ba.write_image("cmp.hdf",-1)
#			try: 
#				b.read_image(args[0],0,0,Region(i[2],i[3],options.box,options.box))
#				b.write_image("cmp.hdf",-1)
#			except: pass
			
			# now we refine this by doing a second pass with the best reference
			try: b.read_image(args[0],0,0,Region(i[2],i[3],options.box,options.box))
			except: continue
			b.process_inplace("eman1.normalize.edgemean")
			b.process_inplace("eman1.filter.lowpass.gaussian",{"lowpass":.1})
#			ba=refptcl[i[1]].align("rotate_translate",b,{},"SqEuclidean")
			ba=b.align("rotate_translate",refptcl[i[1]],{},"optvariance",{"matchfilt":1})
			dx=ba.get_attr("translational.dx")
			dy=ba.get_attr("translational.dy")
			da=ba.get_attr("rotational")
			i[2]-= cos(da)*dx+sin(da)*dy
			i[3]-=-sin(da)*dx+cos(da)*dy
			if hypot(dx,dy)>12.0 or ba.get_attr("ovcmp_m")<=0: 
				print '****'
				continue
			
#			refptcl[i[1]].write_image("at.hdf",-1)
#			ba.write_image("at.hdf",-1)

			# for testing, write out intermediate images 
#			b.write_image("t2.hdf",-1)
#			cc=refptcl[i[1]].rot_scale_trans2D(-da,1.0,-(cos(da)*dx+sin(da)*dy),-(-sin(da)*dx+cos(da)*dy))
#			cc-=ba.get_attr("ovcmp_b")
#			cc/=ba.get_attr("ovcmp_m")
#			cc.write_image("t2.hdf",-1)
#			b-=cc
#			b.write_image("t2.hdf",-1)
			
			# now we calculate the particle quality using the final
			# alignment parameters
			try: b.read_image(args[0],0,0,Region(i[2],i[3],options.box,options.box))
			except: continue
			b.process_inplace("eman1.normalize.edgemean")
#			b.process_inplace("eman1.filter.lowpass.gaussian",{"lowpass":.05})
#			print "%d ROT %f"%(n*2,da)
			rr=refptcl[i[1]].rot_scale_trans2D(da,1.0,0,0)
			rr.process_inplace("eman1.normalize")
			
#			b.cmp("optvariance",rr,{"keepzero":1})
#			b*=b.get_attr("ovcmp_m")
#			b+=b.get_attr("ovcmp_b")
#			rr.write_image("a.hdf",-1)
#			b.write_image("a.hdf",-1)

#			score=rr.cmp("quadmindot",b,{"normalize":1})+1.0			# This is 1.0-normalized dot product, ie 0 is best 2 is worst
#			score=rr.cmp("phase",b,{})+rr.cmp("optvariance",b,{"radweight":1,"matchamp":1})/rr.get_xsize()
			score=sqrt(rr.cmp("optvariance",b,{"matchfilt":1}))
#			score=sqrt(b.cmp("optvariance",rr,{"matchfilt":1}))
#			score=b.get_attr("ovcmp_m")*b.get_attr("sigma")
#			if (score<=0) : continue

			# now record the fixed up location
#			goodpks2.append((ba.get_attr("align_score")*ba.get_attr("ovcmp_m"),i[2],i[3],i[1],ba.get_attr("ovcmp_m"),n))
			goodpks2.append((score,i[2],i[3],i[1],ba.get_attr("ovcmp_m"),ba.get_attr("ovcmp_b"),n,score))
			print "%d\t%1.2f\t%1.2f\t%1.1f\t%1.4f\t%1.6f"%(n,ba.get_attr("translational.dx"),ba.get_attr("translational.dy"),ba.get_attr("rotational")*180.0/pi,ba.get_attr("ovcmp_m"),goodpks2[-1][0])
#			ba.write_image("ttt.hdf",-1)
			n+=1
#			display([b,ba,refptcl[i[1]]])
						
		goodpks2.sort()
		
		# now we do 1-D k-means to split the data into 3 groups
		pl=[(goodpks2[0][0],0),((goodpks2[0][0]+goodpks2[-1][0])/2.0,0),(goodpks2[-1][0],0)]
		
		for i in range(20):
			pl2=[(0,0),(0,0),(0,0)]
			for j in goodpks2:
				if j[0]<(pl[0][0]+pl[1][0])/2.0 : pl2[0]=(pl2[0][0]+j[0],pl2[0][1]+1.0)
				elif j[0]>(pl[1][0]+pl[2][0])/2.0 : pl2[2]=(pl2[2][0]+j[0],pl2[2][1]+1.0)
				else : pl2[1]=(pl2[1][0]+j[0],pl2[1][1]+1.0)
#			pl=[(pl2[0][0]/pl2[0][1],pl2[0][1]),(pl2[1][0]/pl2[1][1],pl2[1][1]),(pl2[2][0]/pl2[2][1],pl2[2][1])]
			# complexity here is necessary to deal with empty sets in the k-means
			if pl2[0][1]==0 : pl[0]=(goodpks2[0][0],0)
			else : pl[0]=(pl2[0][0]/pl2[0][1],pl2[0][1])
			if pl2[1][1]==0 : pl[1]=((goodpks2[0][0]+goodpks2[-1][0])/2.0,0)
			else : pl[1]=(pl2[1][0]/pl2[1][1],pl2[1][1])
			if pl2[2][1]==0 : pl[2]=(goodpks2[-1][0],0)
			else : (pl2[2][0]/pl2[2][1],pl2[2][1])
			print pl
		
		
		
		if   (options.threshold<1.0) : thr=pl[0][0]*options.threshold+goodpks2[0][0]*(1.0-options.threshold)
		elif (options.threshold<2.0) : thr=pl[1][0]*(options.threshold-1.0)+pl[0][0]*(2.0-options.threshold)
		elif (options.threshold<3.0) : thr=pl[2][0]*(options.threshold-2.0)+pl[1][0]*(3.0-options.threshold)
		elif (options.threshold<4.0) : thr=goodpks2[-1][0]*(options.threshold-3.0)+pl[2][0]*(4.0-options.threshold)
		else : thr=goodpks2[-1][0]
		
		print "Threshold : ",thr
		
		# put results in standard 'box' array so GUI can modify if desired
		boxes=[]
		boxthr=thr
		for i in goodpks2:
			boxes.append((i[1],1[2],options.box,options.box,i[0])		# x,y,xsize,ysize,quality
			
# 		out=open("box.stats","w")
# 		for i in goodpks2: out.write("%f\n"%i[0])
# 		out.close()
		
		
	if "circle" in options.auto:
		shrinksq=shrink.copy()
		shrinksq*=shrinksq			# shrunken original image squared
		
		# outer and inner ring mask
		outer=EMData()
		sbox=int(options.box/shrinkfactor)
		outer.set_size(shrink.get_xsize(),shrink.get_ysize(),1)
		outer.to_one()
		inner=outer.copy()
		
		outer.process_inplace("eman1.mask.sharp",{"inner_radius":sbox*2/5,"outer_radius":sbox/2})
		inner.process_inplace("eman1.mask.sharp",{"outer_radius":sbox*2/5})
		
#		outer.write_image("b_outer.hdf")
#		inner.write_image("b_inner.hdf")

		ccf1=shrinksq.calc_ccf(inner,fp_flag.CIRCULANT)
		ccf2=shrinksq.calc_ccf(outer,fp_flag.CIRCULANT)
		
#		ccf1.write_image("b_ccf1.hdf")
#		ccf2.write_image("b_ccf2.hdf")
		
		ccf1/=ccf2
		
#		ccf1.write_image("b_result.hdf")
	
	E2end(logid)

	# invoke the GUI if requested
	if options.gui:
		(boxes,boxthr)=dogui(args[0],boxes,boxthr)

	if options.dbout:
		# Write EMAN1 style box database
		out=open(options.dbout,"w")
		n=0
		for i in boxes:
			if i[4]>boxthr : continue
			out.write("%d\t%d\t%d\t%d\t-3\n"%(i[0],i[1],i[2],i[3]))		
		out.close()

	if options.ptclout
		# write boxed particles
		n=0
		for i in boxes:
			if i[4]>boxthr : continue
			try: b.read_image(args[0],0,0,Region(i[0],i[1],i[2],i[3]))
			except: continue
			if options.norm: b.process_inplace("eman1.normalize.edgemean")
#			print n,i
#			print i[4]
			b.write_image(options.ptclout,n)
			n+=1
			if options.savealiref:
				refptcl[i[3]].write_image("boxali.hdf",-1)
				b.write_image("boxali.hdf",-1)
				

def dogui(imagefsp,boxes,thr):
	"""Implements the 'boxer' GUI. image is the entire image, and boxes and thr specify current boxes
	to begin with. Modified boxes/thr are returned."""
	try:
		from PyQt4 import QtCore, QtGui, QtOpenGL
		from PyQt4.QtCore import Qt
	except:
		print "PyQt4 must be installed to use the --gui option"
		sys.exit(1)
	try:
		from emimagemx import *
		from emimage2d import *
	except:
		print "Cannot import EMAN image GUI objects (emimage,etc.)"
		sys.exit(1)

	guiapp = QtGui.QApplication([])
	
	image=EMData()
	image.read_image(imagefsp,0)
	
	boxupdate()
	
	ptcl=[]
	for i in boxes:
		if i[4]>boxthr: continue
		im=image.get_clip(i[0],i[1],i[2],i[3])
		ptcl.append(im)
		
	guiim=EMImage2D(image)
	guiim.show()
	
	guimx=EMImageMX(ptcl)
	guimx.show()
	
	app.exec_()
	
if __name__ == "__main__":
	main()
