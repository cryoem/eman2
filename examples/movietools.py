"""some utility functions for playing with movie-mode sequences"""
import EMAN2
from math import sqrt

def readstack(n):
	"read the data for the nth particle. Return ref,ptclframes"
	ref=EMAN2.EMData("tmp2.hdf",n*33)
	ptcls=EMAN2.EMData.read_images("tmp2.hdf",range(n*33+1,(n+1)*33))
	return ref,ptcls

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


