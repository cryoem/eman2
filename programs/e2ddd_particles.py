#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/11/14 (sludtke@bcm.edu)
# Copyright (c) 2000- Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#

import pprint
from EMAN2 import *
import sys
from numpy import *
import numpy.linalg as LA

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <refine_xx folder>

	Based on a completed refinement, this will perform per-particle alignment using reference projections from the reconstruction.
	It will only perform alignments on particles used in the specified refinement run. Warning, this will replace the non ctf-corrected
	particles in the particles directory in-place, overwriting the originals. You may wish to consider making a backup of the project
	before running this program.
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	#parser.add_argument("--align_frames", action="store_true",help="Perform whole-frame alignment of the stack",default=False)
	#parser.add_argument("--save_aligned", action="store_true",help="Save aligned stack",default=False)
	#parser.add_argument("--dark",type=str,default=None,help="Perform dark image correction using the specified image file")
	#parser.add_argument("--gain",type=str,default=None,help="Perform gain image correction using the specified image file")
	#parser.add_argument("--gaink2",type=str,default=None,help="Perform gain image correction. Gatan K2 gain images are the reciprocal of DDD gain images.")
	#parser.add_argument("--step",type=str,default="1,1",help="Specify <first>,<step>,[last]. Processes only a subset of the input data. ie- 0,2 would process all even particles. Same step used for all input files. [last] is exclusive. Default= 1,1 (first image skipped)")
	#parser.add_argument("--frames",action="store_true",default=False,help="Save the dark/gain corrected frames")
	#parser.add_argument("--normalize",action="store_true",default=False,help="Apply edgenormalization to input images after dark/gain")
	#parser.add_argument("--movie", type=int,help="Display an n-frame averaged 'movie' of the stack, specify number of frames to average",default=0)
	#parser.add_argument("--simpleavg", action="store_true",help="Will save a simple average of the dark/gain corrected frames (no alignment or weighting)",default=False)
	#parser.add_argument("--avgs", action="store_true",help="Testing",default=False)
	parser.add_argument("--parallel", default=None, help="parallelism argument. This program supports only thread:<n>")
	parser.add_argument("--threads", default=1,type=int,help="Number of threads to run in parallel on a single computer when multi-computer parallelism isn't useful", guitype='intbox', row=24, col=2, rowspan=1, colspan=1, mode="refinement[4]")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
	
	if len(args)!=1:
		print usage
		parser.error("Specify refine_xx folder")

	if options.parallel!=None :
		if options.parallel[:7]!="thread:":
			print "ERROR: only thread:<n> parallelism supported by this program. It is i/o limited."
			sys.exit(1)
		threads=int(options.parallel[7:])
	else: threads=1

	pid=E2init(sys.argv)

	refineparms=js_open_dict(args[0]+"/0_refine_parms.json")
	inlst=refineparms["input"]
	if inlst[0][-4:]!=".lst" :
		print "Error: refine_xx must be run with a 'set' as --input following canonical EMAN2.1 guidelines"
		sys.exit(1)
		
	clsout=sorted([args[0]+"/"+i for i in os.listdir(args[0]) if "cls_result" in i])[-2:]

	if options.verbose: print "running on:",clsout

	newproj="../"+os.getcwd().split("/")[-1]+"_m"
	try: os.makedirs(newproj+"/particles")
	except: pass

	lastloc=None
	# eo is 0/1 for even/odd files
	for eo in xrange(2):
		if options.verbose: print "EO: ",eo
		cls=EMData.read_images(clsout[eo])		# not normally all that big, so we just read the whole darn thing
		projfsp=clsout[eo].replace("cls_result","projections")
		
		lst=LSXFile(inlst[eo])				# open a convenient object for accessing images from the input stack
		
		# i is particle number in the cls file
		for i in xrange(lst.n):
			ptloc=lst.read(i)		# ptloc is n,filename for the source image
			if options.verbose: print "{}/{}   ({})".format(i,lst.n,ptloc)
			
			# if the input particle file changed, we need some new info
			if lastloc!=ptloc[1]:
				movie="movieparticles/{}_ptcls.hdf".format(base_name(ptloc[1]))		# movie particle stack
				movieim=EMData(movie,0)		# number of frames in each movie for this stack
				movien=movieim["movie_frames"]
				
				# we construct a phase flipping image from the first particle CTF info
				ptcl=EMData(ptloc[1],ptloc[0])		# the original particle image, should have CTF info too
				ctf=ptcl["ctf"]
				flipim=movieim.do_fft()		# we're just getting a complex image of the right size, a bit stupid way to handle it 
				ctf.compute_2d_complex(flipim,Ctf.CtfType.CTF_SIGN)
				lastloc=ptloc[1]
			
			proj=EMData(projfsp,int(cls[0][0,i]))	# projection image for this particle
#			orient=Transform({"type":"2d","tx":cls[2][0,i],"ty":cls[3][0,i],"alpha":cls[4][0,i],"mirror":int(cls[5][0,i])})
			orient=Transform({"type":"2d","tx":0,"ty":0,"alpha":cls[4][0,i],"mirror":int(cls[5][0,i])})		# we want the alignment reference in the middle of the box
			proj.transform(orient.inverse())
			
			stack=EMData.read_images(movie,xrange(movien*ptloc[0],movien*(ptloc[0]+1)))
			for im in stack: im.process_inplace("normalize.edgemean")
			avg=sum(stack)
			avg.mult(1.0/len(stack))

			ptcl=EMData(ptloc[1],ptloc[0])		# the original particle image, should have CTF info too
			ptcl.process_inplace("normalize.edgemean")
			
			oldbox=proj["nx"]
			moviebox=stack[0]["nx"]
			
			# CTF phase flipping of the projection so the alignment works
			proj.process_inplace("normalize.edgemean")
			if oldbox!=moviebox : 
				proj=proj.get_clip(Region(-(moviebox-oldbox)/2,-(moviebox-oldbox)/2,moviebox,moviebox))
				ptcl=ptcl.get_clip(Region(-(moviebox-oldbox)/2,-(moviebox-oldbox)/2,moviebox,moviebox))
			pfft=proj.do_fft()
			pfft.mult(flipim)
			proj=pfft.do_ift()
			
#			if proj["nx"]!=ptcl["nx"] : proj=proj.get_clip
			if options.verbose>3 : proj.write_image("tmp.hdf",-1)
			#ptcl.process_inplace("normalize.toimage",{"to":proj})
			if options.verbose>3 :ptcl.write_image("tmp.hdf",-1)
			avg.process_inplace("normalize.edgemean")
			if options.verbose>3 :avg.write_image("tmp.hdf",-1)
			
			# This function is the actual stack alignment to the reference, producing the aligned average
			newpt=alignstack(proj,stack,options.verbose-1)
			if options.verbose>3 :
				newpt.process_inplace("normalize.edgemean")
				newpt.write_image("tmp.hdf",-1)
			if newpt["nx"]!=oldbox : newpt=newpt.window_center(oldbox)

			# write the aligned, clipped average to the output file
			newpt.process_inplace("normalize.edgemean")
			newpt.write_image("{}/particles/{}_ptcls.hdf".format(newproj,base_name(ptloc[1])),ptloc[0])
			
			if options.verbose>1 : print i,movie,ptloc[0],int(cls[0][0,i])
			
	E2end(pid)

def alignstack(refo,stack,verbose=0):
	"""aligns a movie stack to a reference image, with some constraints forcing the relative alignments to follow a pattern
	Returns the aligned average."""
	
	global cnt
	nx=refo["nx"]
	ny=refo["ny"]
	
	# a little preprocessing on the stack
	outim=[i.get_clip(Region(-nx/2,-ny/2,nx*2,ny*2)) for i in stack]
	
	for i in outim:
		i.scale(2.0)
		i.process_inplace("filter.highpass.gauss",{"cutoff_abs":.002})
		i.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.05})
	
	ref=refo.get_clip(Region(-nx/2,-ny/2,nx*2,ny*2))
	ref.scale(2.0)

	nx*=2
	ny*=2
	
	xali=XYData()			# this will contain the alignments which are hierarchically estimated and improved
	yali=XYData()			# x is time in both cases, y is x or y
	for it in xrange(1):
		step=len(outim)
		Zs=[]
		Zthr=0
		
		while step>1:
			step/=2
			i0=0
			while i0<len(outim):
				i1=min(i0+step,len(outim))
				if i1-i0<step/2 : 
					i0+=step
					continue
#				av=sum([outim[i].process("xform.translate.int",{"trans":(int(xali.get_yatx_smooth(i,1)),int(yali.get_yatx_smooth(i,1)))}) for i in xrange(i0,i1)])
				av=sum(outim[i0:i1])

				tloc=(i0+i1-1)/2.0		# the "time" of the current average
				if step>=len(outim)/2 : lrange=nx/4
				else: lrange=hypot(xali.get_yatx_smooth(i1,1)-xali.get_yatx_smooth(i0,1),yali.get_yatx_smooth(i1,1)-yali.get_yatx_smooth(i0,1))
				if lrange<6 : lrange=6	
				
				guess=(xali.get_yatx_smooth(tloc,1),yali.get_yatx_smooth(tloc,1))
				
				if verbose :print step,i0,guess,lrange,

				ccf=av.calc_ccf(ref,fp_flag.CIRCULANT,1)		# centered CCF
				ccf.process_inplace("normalize.edgemean")
				csig=ccf["sigma"]
				cmean=ccf["mean"]
#				ccf.write_image("xyz.hdf",-1)
				ccf.clip_inplace(Region(nx/2-lrange-guess[0],ny/2-lrange-guess[1],lrange*2,lrange*2))
				dx,dy,dz=ccf.calc_max_location()
				try: Z=(ccf[dx,dy]-cmean)/csig
				except: 
					print "\nBad alignment {},{} on {} to {}".format(dx,dy,i0,i1)
					i0+=step
					continue
				
				if verbose: print "###",dx-lrange,dy-lrange,
				
				if dx>0 and dy>0 and dx<lrange*2-1 and dy<lrange*2-1 :
					wx=(ccf[dx-1,dy-1]+ccf[dx-1,dy]+ccf[dx-1,dy+1],ccf[dx,dy-1]+ccf[dx,dy]+ccf[dx,dy+1],ccf[dx+1,dy-1]+ccf[dx+1,dy]+ccf[dx+1,dy+1]) # peak weights in x around peak
					wy=(ccf[dx-1,dy-1]+ccf[dx,dy-1]+ccf[dx+1,dy-1],ccf[dx-1,dy-1]+ccf[dx,dy-1]+ccf[dx+1,dy-1],ccf[dx-1,dy-1]+ccf[dx,dy-1]+ccf[dx+1,dy-1]) # peak weights in y around peak
					dx=((dx-1)*wx[0]+dx*wx[1]+(dx+1)*wx[2])/sum(wx)
					dy=((dy-1)*wy[0]+dy*wy[1]+(dy+1)*wy[2])/sum(wy)
				
				
				
				dx-=lrange
				dy-=lrange
				if verbose: print dx,dy,Z
				
				Zs.append(Z)
				
				# we only include points where the alignment seems reliable
				if step>=len(outim)/2 or Z>Zthr :
					xali.insort(tloc,guess[0]-dx)
					yali.insort(tloc,guess[1]-dy)
					xali.dedupx()
					yali.dedupx()
									
				i0+=step
				
			if Zthr==0:
				Zthr=max((Zs[0]+Zs[1])/3.0,3.0)
				if verbose : print "Z threshold: ",Zthr
			# Smoothing
			# we just do a very simplistic smoothing at each step

			for i in xrange(xali.get_size()-2):
				xali.set_y(i+1,(xali.get_y(i)+xali.get_y(i+1)*2.0+xali.get_y(i+2))/4.0)
				yali.set_y(i+1,(yali.get_y(i)+yali.get_y(i+1)*2.0+yali.get_y(i+2))/4.0)
				
			if xali.get_size()>=4:
				# smooth first point
				lp=xali.get_size()-1
				sx0=xali.get_x(1)-xali.get_x(0)
				sx1=xali.get_x(lp-1)-xali.get_x(1)
				sl=(xali.get_y(lp-1)-xali.get_y(1))/sx1
				xali.set_y(0,(xali.get_y(1)-sl*sx0+xali.get_y(0))/2.0)
				sl=(yali.get_y(lp-1)-yali.get_y(1))/sx1
				yali.set_y(0,(yali.get_y(1)-sl*sx0+yali.get_y(0))/2.0)

				# smooth last point
#				print xali.get_y(lp),yali.get_y(lp),"->",
				sx1=xali.get_x(lp-1)-xali.get_x(1)
				sx0=xali.get_x(lp)-xali.get_x(lp-1)
				sl=(xali.get_y(lp-1)-xali.get_y(1))/sx1
				xali.set_y(lp,(xali.get_y(lp-1)+sl*sx0+xali.get_y(lp))/2.0)
				sl=(yali.get_y(lp-1)-yali.get_y(1))/sx1
				yali.set_y(lp,(yali.get_y(lp-1)+sl*sx0+yali.get_y(lp))/2.0)
#				print xali.get_y(lp),yali.get_y(lp)
			   

			if verbose>1 :
				print ["%6.1f"%i for i in xali.get_xlist()]
				print ["%6.2f"%i for i in xali.get_ylist()]
				print ["%6.2f"%i for i in yali.get_ylist()]

			#out=file("path_{}_{}.txt".format(cnt,step),"w")
			#for i in xrange(len(outim)): out.write("{}\t{}\t{}\n".format(i,xali.get_yatx_smooth(i,1),yali.get_yatx_smooth(i,1)))

	#cnt+=1
			
		
	av=sum([stack[i].process("xform.translate.int",{"trans":(int(xali.get_yatx_smooth(i,1)/2.0),int(yali.get_yatx_smooth(i,1)/2.0))}) for i in xrange(len(outim))])
	av.mult(1.0/len(outim))
	
	
	if verbose :
		print ["%6.1f"%i for i in xali.get_xlist()]
		print ["%6.2f"%i for i in xali.get_ylist()]
		print ["%6.2f"%i for i in yali.get_ylist()]
	

	return av
	
cnt=1

if __name__ == "__main__":
	main()
