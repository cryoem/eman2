"""some utility functions for playing with movie-mode sequences"""
import EMAN2
from math import sqrt

def readstack(n):
	"read the data for the nth particle. Return ref,ptclframes,psref"
	ref=EMAN2.EMData("tmp2.hdf",n*33)
	ptcls=EMAN2.EMData.read_images("tmp2.hdf",range(n*33+1,(n+1)*33))
	ctfim=ref.do_fft()
	ctf=ref["ctf"]
	ctf.compute_2d_complex(ctfim,EMAN2.Ctf.CtfType.CTF_POWEVAL)
	ctfim.ri2ap()
	ctfim=ctfim.amplitude()
#	ctfim.process_inplace("normalize")

	return ref,ptcls,ctfim

def seqavg(stack):
	"return progressive averages of frames in stack"
	return [sum(stack[:i+1])*(1.0/(sqrt(i+1.0))) for i in xrange(len(stack))]

def runavg(stack,n):
	"return running averages by n "
	return [sum(stack[i:i+n])*(1.0/sqrt(n)) for i in xrange(len(stack)-n+1)]

def ccfs(ref,stack):
	"compute and center CCFs between ref and each member of stack"
	ret=[i.calc_ccf(ref) for i in stack]
	f=ref["nx"]/4
	for i in ret: i.process_inplace("xform.phaseorigin.tocenter")
	ret=[i.get_clip(EMAN2.Region(f,f,f*2,f*2)) for i in ret]

	return ret

def powspec(img):
	"compute 1D power spectrum"
	if not img.is_complex() :
		img=img.do_fft()

	ps=img.calc_radial_dist(img["ny"]/2,0,1,True)
	ps=[i/img["nx"]**2 for i in ps]

	return ps

def peaks(stack):
	"find peak locations from ccfs"

	nx=stack[0]["nx"]
	pk=[i.calc_max_location() for i in stack]
	pk=[(i[0]-nx/2,i[1]-nx/2) for i in pk]
	return pk

def localali(im1,im2,psref,maxdx):
	"""Takes an fft.amplitude() psref to use as a basis for optimal alignment. Exhaustive search +-maxdx pixels.
	returns (dx,dy) required to bring im2 into alignment with im1"""

	nx=im1["nx"]
	mask=EMAN2.EMData(nx,nx,1)
	mask.to_one()
	mask.process_inplace("mask.gaussian",{"inner_radius":nx/4,"outer_radius":nx/12})
#	EMAN2.display(mask)

	maxdx=int(maxdx)
	out=EMAN2.EMData(2*maxdx+1,2*maxdx+1,1)
	for dx in xrange(-maxdx,maxdx+1):
		for dy in xrange(-maxdx,maxdx+1):
			av=im1+im2.process("xform.translate.int",{"trans":(dx,dy,0)})
			av.mult(mask)	# to prevent edge effects and slightly smooth the powspec
			avf=av.do_fft()
			avf.ri2ap()
			avf=avf.amplitude()
			out[dx+maxdx,dy+maxdx]=avf.cmp("ccc",psref)
#			print dx,dy,out[dx+20,dy+20]

	ret=out.calc_max_location()
	return (int(ret[0]-maxdx),int(ret[1]-maxdx)),out

def stackaliloccor(ref,stack):
	"""aligns images in stack to ref based on the assumption that movements follow a path and don't jump around much"""
	nx=ref["nx"]
	ny=ref["ny"]

	# Should get better results if we filter by the CTF or SNR
	ctfim=ref.do_fft()
	ctf=ref["ctf"]
	ctf.bfactor=50
	ctf.compute_2d_complex(ctfim,EMAN2.Ctf.CtfType.CTF_INTEN)

	# prefilter ref so we don't have to do it over and over
	reff=ref.do_fft()
	reff.mult(ctfim)
	reff=reff.do_ift()

	# We compute the CCF for the unaligned average, then filter out values close to the max
	# to define a region of permissible translation for individual frames
	unaliavg=sum(stack)
	ccfmask=unaliavg.calc_ccf(reff,EMAN2.fp_flag.CIRCULANT,True)
	ccfmask.process_inplace("normalize.edgemean")
	ccfmask.process_inplace("threshold.binary",{"value":ccfmask["maximum"]*.8})
	ccfmask.process_inplace("mask.addshells",{"nshells":2})

#	EMAN2.display(ccfmask)

	for i,im in enumerate(stack):
		ccf=im.calc_ccf(reff,EMAN2.fp_flag.CIRCULANT,True)
		ccf.mult(ccfmask)
		pk=ccf.calc_max_location()
		dx=-(pk[0]-nx/2)
		dy=-(pk[1]-ny/2)
		print i,dx,dy

		try: avg.add(im.process("xform.translate.int",{"trans":(dx,dy)}))
		except: avg=im.process("xform.translate.int",{"trans":(dx,dy)})
		
	return avg

def stackali(ref,stack,psref,maxdx):
	"""aligns images in stack to ref, under constraint of maximizing psref"""

	unaliavg=sum(stack)
	ali=unaliavg.align("translational",ref)
	
	# This is the overall alignment of the average without considering psref
	# just to help us search in the correct neighborhood
	dx0,dy0=[int(i) for i in ali["xform.align2d"].get_trans_2d()]
	print "overall: ",dx0,dy0

	# to avoid a lot of gymnastics, and because the edge of the reference is "flat"
	# we shift the ref off-center for the alignments, then compensate at the end
	xfref=ref.process("xform.translate.int",{"trans":(-dx0,-dy0,0)})
	#EMAN2.display((ref,unaliavg,xfref),True)

	# align to the reference under psref
	for i,im in enumerate(stack):
		dx,dy=localali(xfref,im,psref,maxdx)
		try: avg1.add(im.process("xform.translate.int",{"trans":(dx,dy)}))
		except: avg1=im.process("xform.translate.int",{"trans":(dx,dy)})
		print i,dx,dy

	print "------"

	# now align to the initial average under psref (in case the reference is bad?)
	for i,im in enumerate(stack):
		dx,dy=localali(avg1,im,psref,maxdx)
		try: avg2.add(im.process("xform.translate.int",{"trans":(dx+dx0,dy+dy0)}))
		except: avg2=im.process("xform.translate.int",{"trans":(dx+dx0,dy+dy0)})
		print i,dx,dy
	
	return avg1,avg2
