#!/usr/bin/env python

import numpy
from numpy import linalg as LA
from numpy import array, dot

from EMAN2 import *
from sparx import *
from mpi import *
from math import sqrt
from random import random
import sys

from os import system
from sys import argv
from sys import exit
from string import atof,atoi

from optparse import OptionParser

def mean_var_skew(data):
	nsize=len(data)

	if nsize<=0:
		return 0.0,0.0,0.0

	sum1=0.0; sum2=0.0; sum3=0.0;
	for i in range(nsize):
		sum1=sum1+data[i]
		sum2=sum2+data[i]**2
		sum3=sum3+data[i]**3
		
	ave1=sum1/nsize
	ave2=sum2/nsize
	ave3=sum3/nsize
	var = ave2-ave1**2
	skew = (ave3-3*ave1*var-ave1**3)/var**1.5

	return ave1, var, skew


def shift_gamma(F, gr):
	L = len(F)
	ds = 0.5/(L-1)
	fg = [0.0]*L
	for i in xrange(L):
		sp = 0.5*(2*i*ds)**gr
		ip = sp/ds; i1 = int(floor(ip))
		if i1<L-1:
			fg[i] = F[i1]*(i1+1.0-ip)+F[i1+1]*(ip-i1)
		else:
			fg[i] = F[i1]
	return fg
	

def create_CCFR_filter(pwu, pwp, fcrf, sg):
	# tpu := Re(T)/|P||U| (with T = |P|^2 or |P| or Gaussian):
	tpu = []; m=0.0
	csg = 0.5/sg**2
	lp = min(len(fcrf),len(pwp))
	for i in range(lp):
		#z = sqrt(pwp[i]/pwu[i])                 # T = |P|^2
		#z = 1.0/sqrt(pwu[i])                    # T = |P|
		z = exp(-i**2*csg)/sqrt(pwp[i]*pwu[i])   # T = Gaussian
		if z>m: m=z
		tpu.append(z)
	for i in range(lp):
		tpu[i] /= m
	
	# filter:
	wfilt = []; m=0.0
	for i in range(lp):
		z = fcrf[i]*tpu[i]
		if z>m: m=z
		wfilt.append(z)
	for i in range(lp):
		wfilt[i] /= m

	return wfilt


def create_noise_filter(pwG, FSC, q2):
	L = len(pwG)
	filt = [0.0]*L
	for i in xrange(min(L,len(FSC))):
		filt[i] = sqrt(pwG[i]*(1.0-FSC[i])/q2)
	return filt


def map_adjustment(volg, FSC):
	L = len(FSC)
	filt_adj = [0.0]*L
	for i in xrange(len(FSC)):
		filt_adj[i] = sqrt(FSC[i])
	volg_adj = filt_table(volg, filt_adj)
	return volg_adj


def norm_factor_noise(nx, ny, nz):
	pw = rops_table(model_gauss_noise(1.0, nx, ny, nz))
	start = len(pw)/10
	# old way:
	#oofactor = sqrt(sum(pw[start:])/len(pw[start:]))
	# new way, more accurate:
	s1=0.0; s2=0.0
	for i in range(start,len(pw)):
		s1 += pw[i]*i**2
		s2 += i**2
	oofactor = sqrt(s1/s2)
	return 1.0/oofactor


def create_noisy_map(vol, sigma_noise, filt, factor_noise):
	nx = vol.get_xsize()
	ny = vol.get_ysize()
	nz = vol.get_zsize()
	if sigma_noise==0:
		imgNoisy = vol.copy()
	else:
		imgNoisy = vol + filt_table(model_gauss_noise(sigma_noise*factor_noise, nx, ny, nz), filt)
	return imgNoisy


def fsc_mask(vol, sigma, level_ratio):
	vf = filt_gaussl(vol, sigma) 
	minval = level_ratio*Util.infomask(vf, None, True)[3]
	#return gauss_edge(binarize(vf, minval))
	return binarize(vf, minval)


def q2norm(fmask):
	# normalized squared L^2-norm of fmask:
	return Util.infomask(fmask*fmask, None, True)[0]

	
def compdisc(x, y, mask, discrep):
	return -x.cmp(discrep, y, {"mask":mask})


def ess_p(data,p=0.05):
	# find point in data corresponding to probability p
	# data is assumed already sorted in increasing order
	lend = len(data)

	if (p==0.0):
		ip = 0
	else:
		i0 = int(floor(p*lend))
		if(i0==p*lend):
			ip=i0-1
		else:
			ip=i0

	return data[ip]


def fitgamma(R,B,t,nu):
	# finds the gamma such that sfg(...,gamma) = 0,
	# by the bisection method.
	# t = target value (e.g., integral of fcr+)
	# nu = tilt parameter

	jmax = 40    # max number of iterations allowed
	acc = 0.001  # desired accuracy

	# intial bracket for gamma:
	g1 = 0.0
	g2 = 5.0

	fmid = sfg(R,B,t,nu,g2)
	f = sfg(R,B,t,nu,g1)

	if f*fmid>0.0:
		print "Error: inappropriate gamma bracket."
		mpi_finalize()
		exit(-1)

	if f<0.0:
		rtbis = g1
		dx = g2-g1
	else:
		rtbis = g2
		dx = g1-g2

	for j in range(jmax):
		dx *= 0.5
		xmid = rtbis+dx
		fmid = sfg(R,B,t,nu,xmid)
		if fmid<=0.0: rtbis = xmid
		if (abs(dx)<acc or fmid==0.0):
			#print "number of iterations for gamma:",j+1
			return rtbis

	print "Error: max number of iterations reached."
	mpi_finalize()
	exit(-2)


def sfg(R,B,t,nu,gamma):
	# integral(sqrt(FSC_{nu,gamma}))-t
	np = 100   # intervals for integration
	sf = 0.0
	ds = 0.5/float(np)
	for i in range(np+1):
		sf += fff(R,B,nu,gamma,i*ds)
	sf *= ds
	return sf-t


def fff(R,B,nu,gamma,s):
	# compute sqrt(FSC_{nu,gamma}(s))
	# (uses option 4 -- see notes)
	v = (2*s)**gamma
	if v<=0.5:
		u = 0.5*(2*v)**nu
	else:
		u = 1.0-0.5*(2*(1-v))**nu
	y = 0.5*u
	fs = sqrt((1.0+R)/(1.0+R*exp(B*y**2/4)))
	if fs<=0.5:
		u = 0.5*(2*fs)**nu
	else:
		u = 1.0-0.5*(2*(1-fs))**nu
	# a more exact way would be to cap u with sqrt(FSC(s))
	return u


def minfind(R,B,gamma,fcr):
	# finds the minimum (over nu) of int[(sqrt(FSC_{nu,gamma}-fcr)^2]
	# initial bracket for nu:
	nu1 = 0.5
	nu2 = 5.0
	np = 10        # number of subdivisions
	jmax = 40      # max number of iterations
	acc = 0.001    # desired accuracy
	for j in range(jmax):
		dx = (nu2-nu1)/float(np)
		mf = integ2(R,B,nu1,gamma,fcr)
		mi = 0
		for i in range(1,np+1):
			cf = integ2(R,B,nu1+dx*i,gamma,fcr)
			if cf<mf:
				mf = cf
				mi = i
		if mi==0:
			nu2 = nu1+dx
		if mi==np:
			nu1 = nu2-dx
		if mi>0 and mi<np:
			nu1 = nu1+dx*(mi-1)
			nu2 = nu1+dx*2
		if nu2-nu1<acc:
			return (nu1+nu2)/2.0

	print "Error: max number of iterations reached."
	mpi_finalize()
	exit(-3)


def integ2(R,B,nu,gamma,fcr):
	# computes int[(sqrt(FSC_{nu,gamma}-fcr)^2]
	np = len(fcr)
	ds = 0.5/float(np-1)
	sf = 0.0
	for i in range(np):
		sf += (fff(R,B,nu,gamma,i*ds)-fcr[i])**2
	sf *= ds
	return sf


def fillholes(f,t):
	# t=0: fill all holes
	# t=1: fill first hole only
	# t=2: polygonal line from max to max
	# t=3: like t=2 but omitting maxima that are < later ones
	lf = len(f)
	ff = [0.0]*lf
	if t==0 or t==1:
		tol = 0.03
		up = 0; dn = 0
		mic=0.0; mac=0.0
		ff[0] = f[0]
		for i in range(1,lf):
			ff[i] = f[i]
			if f[i]<=f[i-1]:
				dn = 1
				mic = f[i]
			if f[i]>f[i-1]:
				up = 1
				mac = f[i]
			if dn==1 and f[i]>mic+tol:
				for j in range(i):
					if ff[j]<f[i]: ff[j]=f[i]
			if t==1 and up==1 and f[i]<mac-tol:
				for j in range(i,lf):
					ff[j] = f[j]
				break
	elif t==2 or t==3:
		fp = [0]*lf    # indicator of local max of f
		fp[0] = 1; fp[lf-1] = 1
		for i in range(1,lf-1):
			if f[i]>=f[i-1] and f[i]>=f[i+1]:
				fp[i] = 1

		if t==3:     # eliminate maxima that are smaller than later ones:
			for i in range(1,lf-1):
				if fp[i]==1:
					for j in range(1,i):
						if fp[j]==1 and f[j]<f[i]:
							fp[j] = 0

		fm = [0.0]*lf  # polygonal
		fm[0] = f[0]
		i = 0                          # i=left endpoint of interval
		while i<lf-1:
			for j in range(i+1,lf):    # j=right endpoint
				if fp[j]==1:
					slope = (f[j]-f[i])/float(j-i)
					for k in range(i+1,j+1):
						fm[k] = slope*(k-i) + f[i]
					j1 = j
					break
			i = j1

		# return max between polygonal and given f:
		for i in range(lf):
			if f[i]>fm[i]:
				ff[i] = f[i]
			else:
				ff[i] = fm[i]
	return ff


def ellipsoid(piece):
	# computes the inertia ellipsoid of piece
	nx = img1_orig_0.get_xsize()
	ny = img1_orig_0.get_ysize()
	nz = img1_orig_0.get_zsize()

	mass = 0.0
	for x in range(nx):
		for y in range(ny):
			for z in range(nz):
				mass += piece[x,y,z]

	xb = 0.0
	for x in range(nx):
		ps = 0.0
		for y in range(ny):
			for z in range(nz):
				ps += piece[x,y,z]
		xb += ps*x
	xb /= mass
	
	yb = 0.0
	for y in range(ny):
		ps = 0.0
		for x in range(nx):
			for z in range(nz):
				ps += piece[x,y,z]
		yb += ps*y
	yb /= mass
	
	zb = 0.0
	for z in range(nz):
		ps = 0.0
		for x in range(nx):
			for y in range(ny):
				ps += piece[x,y,z]
		zb += ps*z
	zb /= mass
	
	c11 = 0.0
	for x in range(nx):
		ps = 0.0
		for y in range(ny):
			for z in range(nz):
				ps += piece[x,y,z]
		c11 += ps*(x-xb)**2
	c11 /= mass

	c22 = 0.0
	for y in range(ny):
		ps = 0.0
		for x in range(nx):
			for z in range(nz):
				ps += piece[x,y,z]
		c22 += ps*(y-yb)**2
	c22 /= mass

	c33 = 0.0
	for z in range(nz):
		ps = 0.0
		for x in range(nx):
			for y in range(ny):
				ps += piece[x,y,z]
		c33 += ps*(z-zb)**2
	c33 /= mass

	c12 = 0.0
	for x in range(nx):
		for y in range(ny):
			ps = 0.0
			for z in range(nz):
				ps += piece[x,y,z]
			c12 += ps*(x-xb)*(y-yb)
	c12 /= mass

	c13 = 0.0
	for x in range(nx):
		for z in range(nz):
			ps = 0.0
			for y in range(ny):
				ps += piece[x,y,z]
			c13 += ps*(x-xb)*(z-zb)
	c13 /= mass

	c23 = 0.0
	for y in range(ny):
		for z in range(nz):
			ps = 0.0
			for x in range(nx):
				ps += piece[x,y,z]
			c23 += ps*(y-yb)*(z-zb)
	c23 /= mass

	inmat = array([[c11,c12,c13],[c12,c22,c23],[c13,c23,c33]])

	comat = LA.inv(inmat)

	# inertia-ellipsoid map:
	mask2 = model_blank(nx,nx,nx)
	for x in range(nx):
		for y in range(ny):
			for z in range(nz):
				v = array([x-xb,y-yb,z-zb])
				mask2[x,y,z] = exp(-0.5*dot(v,dot(comat,v)))

	# test:
	#if myid == main_node:
	#	drop_image(mask2,"ellipsoid.spi","s")

	return mask2


def fitfcr(fcr1, R, B):
	# computes fitted FCR curve
	lf = len(fcr1)
	
	# fill up holes:
	fcr2 =  fillholes(fcr1,2)
	
	# set to 0 after first drop below 0.1:
	for i in range(lf):
		if fcr2[i]<0.1:
			for j in range(i,lf):
				fcr2[j] = 0.0
			break

	# find parameters nu,gamma to fit sqrt(FSC_{nu,gamma}) to the fcr curve:
	sfp = 0.0     # integral of max(fcr2,0) -- could also be some variation
	ds = 0.5/float(lf-1)
	for i in range(lf):
		if fcr2[i]>0.0:
			sfp += fcr2[i]
	sfp *= ds
	nu = 1.0
	for k in range(40):
		gamma = fitgamma(R,B,sfp,nu)
		nup = minfind(R,B,gamma,fcr2)
		if abs(nup-nu)<0.001:
			nu = nup
			break
		nu = nup

	# compute fitted fcr curve with these parameters:
	fcrf = []
	for i in range(lf):
		fcrf.append(fff(R,B,nu,gamma,i*ds))

	# extrapolate fcrf to 0 after first drop below 0.1:
	for i in range(lf):
		if fcrf[i]<0.1:
			k0 = int(fcrf[i]/(fcrf[i]-fcrf[i+1])+i)
			if k0>lf-1:
				k0 = lf-1
				slope = fcrf[i+1]/float(i+1-k0)
			else:
				slope = fcrf[i+1]-fcrf[i]
			for k in range(i+2,k0+1):
				fcrf[k] = slope*(k-i)+fcrf[i]
			for k in range(k0+1,lf):
				fcrf[k] = 0.0
			break

	return fcrf


def red(angle):
	# reduces angle to the interval [-180,180]
	from math import floor
	return angle-360*floor((angle+180)/360)



# read args from command line ---------------------------------
progname = os.path.basename(sys.argv[0])
usage = "mpirun [-np num_proc] [--host node[,node,...]] " + progname + "  map_file  segment_file [FSC_mask] [option_list]"

parser = OptionParser(usage,version=SPARXVERSION)
parser.add_option("--res",      type="float",   default= -1.0, metavar="resolution",   help="Resolution of map. If given, the FSC will be adjusted. If not, the resolution will be deduced from the FSC.")
parser.add_option("--R",        type="float",   default= 0.01, metavar="coeff_R",      help="R coefficient in the fitted FSC curve.")
parser.add_option("--B",        type="float",   default= 500.0,metavar="coeff_B",      help="B coefficient in the fitted FSC curve.")
parser.add_option("--u0",       type="float",   default= 1.0,  metavar="FSC_at_0",     help="Value of the experimental FSC curve at frequency 0.")
parser.add_option("--pix",      type="float",   default= 3.0,  metavar="pixel_size",   help="Pixel size of map, in A.")
parser.add_option("--sigN",   type="float",   default= 3.0,  metavar="sigma_noise",  help="Standard deviation of the Gaussian noise.")
parser.add_option("--alpha",    type="float",   default= 0.0,  metavar="cloud_factor", help="If not 0, an additional cloud will be created, corresponding to sigma_noise*cloud_factor, and the fitting parameters will be scaled by cloud_factor.")
parser.add_option("--fsc_rapid_fall", action="store_true", default= False, help="Do a rapid tapering of the FSC curve.")
parser.add_option("--j1",       type="int",     default= 0,    metavar="rapid_start",  help="Start index of the tapering factor.")
parser.add_option("--j2",       type="int",     default= 0,    metavar="rapid_end",    help="End index of the tapering factor.")
parser.add_option("--synth",          action="store_true", default= False, help="Create and save a simulated map, and exit.")
parser.add_option("--no_sim_noise",   action="store_true", default= False, help="Do not add noise when creating a simulated map.")
parser.add_option("--filtmap",        action="store_true", default= False, help="Filter the map instead of the segment.")
parser.add_option("--imask",          action="store_true", default= False, help="Use a mask to restrain the position of the segment.")
parser.add_option("--init_max",       action="store_true", default= False, help="Choose the initial fit using the max (instead of the median) of the correlation values.")

#parser.add_option("--nover",    type="int",     default= 1,    metavar="oversampling_factor", help="Number of times to oversample input volumes.")

(options, args) = parser.parse_args()    # by default it returns args = sys.argv[1:]

sys.argv = mpi_init(len(sys.argv),sys.argv)
myid = mpi_comm_rank(MPI_COMM_WORLD)
ncpu = mpi_comm_size(MPI_COMM_WORLD)
main_node = 0

if len(args) < 2 or len(args)>3:
	print "usage:  " + usage
	print "Please run '" + progname + " -h' for detailed options"
	mpi_finalize()
	exit(101)

largeFile = args[0]
segFile = args[1]
if len(args)==3:
	maskFile = args[2]
	user_mask = True
else:
	user_mask = False
segName = os.path.splitext(os.path.basename(segFile))[0]

R = options.R
B = options.B
u0 = options.u0
pix = options.pix

if options.res>0.0:
	resol = options.res
	aau = pix/resol
	gamres = 0.5*log(16/B*log((2*u0*(1+R)-1)/R))/log(2*aau)
else:
	aau = sqrt(4/B*log((2*u0*(1+R))/R))
	resol = round(pix/aau,1)
	gamres = 1.0

#nover = options.nover
fsc_rapid_fall = options.fsc_rapid_fall
j1 = options.j1
j2 = options.j2
if fsc_rapid_fall:
	if j1<=0 or j2<=j1:
		print "Warning: input j1,j2 not valid. Will use default values."
		epsf = 0.1   # level under which the FSC is essentially 0
		b1 = sqrt(4/B*log((u0*(1+R)-epsf)/R/epsf))
		b2 = b1+0.25*(0.5-b1)

sigNoise = options.sigN
alpha = options.alpha
synthetic = options.synth
sim_noise = not options.no_sim_noise
filtseg = not options.filtmap
imask = options.imask
init_median = not options.init_max


if myid == main_node:
	print "Options in effect:"
	print "  map_file =", largeFile
	print "  segment_file =", segFile
	if user_mask == True:
		print "  FSC_mask =",maskFile
	else:
		print "  FSC_mask = None"
	print "  res =",resol
	print "  R =",R
	print "  B =",B
	print "  u0 =",u0
	print "  pix =",pix
	print "  sigN =",sigNoise
	print "  alpha =",alpha
	print "  fsc_rapid =",fsc_rapid_fall
	print "  j1 =",j1
	print "  j2 =",j2
	print "  synth =",synthetic
	print "  no_sim_noise =", not sim_noise
	print "  filtmap =", not filtseg
	print "  imask =",imask
	print "  init_max =", not init_median
#-------------------------------------


### parameters not specified by the user: -------------------

# number of times to oversample input volumes:
nover = 1

# discrepancy measure:
discrep = "dot"  #"ccc"  "dot"  "lod"

# parameters for ali_vol:
ang_bracket = 5.0
shift_bracket = 5.0

# number of positions in the random initial shake:
# for the noise-free map:
nt1 = 20/ncpu
if nt1<1: nt1=1
# for the noise-corrupted maps:
nt2 = 10

# amplitudes of the random initial shakes:
# for the noise-free fitting:
shake_shift_1 = 2.0    # pixels
shake_ang_1 = 2.0      # degrees
# for the noise-corrupted fittings:
shake_shift_2 = 2.0    # pixels
shake_ang_2 = 2.0      # degrees

# max number of cycles for monitoring convergence:
max_cycles = 10

# number of terms, minus 1, to average successive sigmas
# during the convergence:
rave = 3

# relative max change in std for ending Monte Carlo:
sigma_max_change = 0.03

# multiplier for the restraint mask (imask):
lambdam = 1.0

# maximum desired number of ouput matrices:
nmatsout = 1000

#----------------------------------------------------------

# begin computations---------------------------------------


# segment map:
img1_orig_0 = get_image(segFile)

nx = img1_orig_0.get_xsize()
ny = img1_orig_0.get_ysize()
nz = img1_orig_0.get_zsize()

s = [0.0]*6
s = center_of_gravity_phase(img1_orig_0)

for i in xrange(1,nx//2):
	mask = cyclic_shift(model_circle(i, nx, ny, nz),int(s[3]),int(s[4]),int(s[5]))
	if Util.infomask(img1_orig_0,mask,False)[3] == 0.0:
		break
radius = i   #smallest sphere that contains the segment

# scale radius according to oversampling:
radius *= nover

# radius of mask:
mask_rad = 1.3*radius

# box size to window images (so that it's a multiple of nover):
nbox = 2*int(mask_rad/nover+1)*nover

# to prevent the window from going out of the boundary of the whole map:
if (nx-nbox)/2 <= abs(s[3]):
	nbox = int(nx-2*abs(s[3])-1)
if (ny-nbox)/2 <= abs(s[4]):
	nbox = int(ny-2*abs(s[4])-1)
if (nz-nbox)/2 <= abs(s[5]):
	nbox = int(nz-2*abs(s[5])-1)
if nbox%2 ==0: nbox=nbox-1


# target map:
img2_orig = get_image(largeFile)

nx0 = img2_orig.get_xsize()
ny0 = img2_orig.get_ysize()
nz0 = img2_orig.get_zsize()

if nover>1:
	img1_orig_0 = fpol(img1_orig_0, nx0*nover, ny0*nover, nz0*nover, True)

if myid == main_node:
	print "Radius of segment =",radius,"pix"

if nover>1:
	img2_orig = fpol(img2_orig, nx0*nover, ny0*nover, nz0*nover, True)

# digitize FSC curve:
sf0 = nx0//2
sf1 = sf0*nover+1
FSC = [0.0]*sf1
dn = 0.5/float(sf0)
for j in xrange(sf1):
	if j<=sf0:
		n = j*dn
		FSC[j] = u0*(1.0+R)/(1.0+R*exp(B*n**2/4))
	else:
		FSC[j] = 0.0

# multiply FSC by rapidly decreasing factor:
if fsc_rapid_fall:
	if j1<=0 or j2<=j1:    # if user-specified are invalid
		j1 = b1*nx0
		j2 = b2*nx0
	for j in xrange(len(FSC)):
		if j<=j1: ff = 1.0
		elif j>=j2: ff = 0.0
		else:
			x = float(j2-j1)*(1.0/(j2-j)-1.0/(j-j1))
			ff = 1.0/(1.0+exp(x))
		FSC[j] *= ff

# shift FSC to simulate desired resolution:
FSC = shift_gamma(FSC, gamres)

# get sqrt(FSC) (for capping):
sqrfsc = []
for j in xrange(len(FSC)):
	sqrfsc.append(sqrt(FSC[j]))

# if user didn't provide an FSC mask, create one:
if(not user_mask):
	sg = 0.02    # sigma of the Gaussian filter in abs. freq. units
	fthr = 0.01  # threshold for binarizing the filtered map (ratio with max)
	fmask = fsc_mask(img2_orig, sg, fthr)
else:
	fmask = get_image(maskFile)
q2 = q2norm(fmask)
if myid == main_node:
	print "squared L^2 norm of FSC mask =", q2


# normalization factor so that FT(Gaussian noise) has sigma=1:
factor_noise = norm_factor_noise(nx0, ny0, nz0)


########################################################################
# IMPORTANT: this block generates a realization of noise that will be
# different across the nodes, if run on several nodes. This will screw up
# the calculation. Hence, the synthetic case should be run on a single node;
# OR the map saved and then read back in as a real map (so it doesn't go
# through this block) -- this is way it's done here;
# OR this block should be changed so that the noise map
# is computed on the main node and then sent to the various nodes.
########################################################################
# if using a synthetic map, apply envelope function and
# add noise according to the given FSC:
if synthetic and sim_noise:
	img2_orig_e = filt_gaussl(img2_orig, aau)
	pwEF = rops_table(img2_orig_e)
	noise_filt_synth = create_noise_filter(pwEF, FSC, q2)
	img2_orig_e = map_adjustment(img2_orig_e, FSC)
	img2_orig = create_noisy_map(img2_orig_e, 1.0, noise_filt_synth, factor_noise)

	if myid == main_node:
		# save this map, to be read back in as a real EM map
		# (this avoids the problem mentioned above):
		drop_image(img2_orig, "map_simulated.spi","s")

		# fcr of the whole structure, for testing:
		fcrwhole = fsc(img2_orig, img2_orig_e)[1]
		for i in range(len(fcrwhole)):
			if FSC[i] == 0.0:
				fcrwhole[i] = 0.0
		filename = segName+"_"+str(resol)+"_fcrwhole.txt"
		write_text_file(fcrwhole, filename)

	mpi_finalize()
	exit(111)


# compute inertia-ellipsoid mask:
mask2 = ellipsoid(img1_orig_0)

# masked map and segment:
um = mask2*img2_orig
pm = mask2*img1_orig_0

# initial FCR curve:
fcr = fsc(um, pm)[1]

# cap fcr to sqrt(FSC):
lf = len(fcr); fcr1 = [0.0]*lf
for i in range(lf):
	if fcr[i]>sqrfsc[i]:
		fcr1[i] = sqrfsc[i]
	else:
		fcr1[i] = fcr[i]

# compute fitted FCR curve:
fcrf = fitfcr(fcr1, R,B)

# cap fcrf to sqrt(FSC):
for i in range(lf):
	if fcrf[i]>sqrfsc[i]:
		fcrf[i] = sqrfsc[i]


# get power spectra for filters:
pwu = rops_table(um)
pwp = rops_table(pm)

# create filter for initial fitting:
sg = min(nx0,ny0,nz0)*aau   # sigma of Gaussian constraint on ccf
Wfilt = create_CCFR_filter(pwu, pwp, fcrf, sg)


# filter target map or segment (to do the initial fitting):
if filtseg:
	img1_filt = filt_table(img1_orig_0, Wfilt)
	img2_filt = img2_orig.copy()
else:
	img1_filt = img1_orig_0.copy()
	img2_filt = filt_table(img2_orig, Wfilt)


# apply restraint mask if requested:
if imask:
	img1_orig_1 = get_image(segFile)
	if nover>1:
		img1_orig_1 = fpol(img1_orig_1, nx0*nover, ny0*nover, nz0*nover, True)
	blob = filt_gaussl(img1_orig_1, 0.1)
	ithres = 0.03 * Util.infomask(blob, None, True)[3]
	imaskvol = binarize(blob, ithres)
	img2_filt = img2_filt * imaskvol


# shift maps by segment's COG, and window:
img1_filt = Util.window(img1_filt, nbox, nbox, nbox, int(s[3]), int(s[4]), int(s[5]))
img2_filt = Util.window(img2_filt, nbox, nbox, nbox, int(s[3]), int(s[4]), int(s[5]))


names_params = ["phi", "theta", "psi", "s3x", "s3y", "s3z"]
lpar = 6

mask = model_circle(mask_rad, nbox, nbox, nbox)


# prepare segment for gridding
#(after this point, img1_filt MUST be used ONLY in rot_shift3D_grid):
img1_filt,kb = prepi3D(img1_filt)


# this is done with 0 to be able to compute the initial bccc:
x = rot_shift3D_grid(img1_filt, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, kb, "background", False)

bccc = compdisc(x, img2_filt, mask, discrep)

if myid == main_node:
	print  "initial",discrep,"=",bccc
bccc = -1.e10

# parameters are in the "XYZ" convention

params = [0.0]*lpar
bparams = [0.0]*lpar

# shift to the center of mass of segment:
tshift = Transform({"type":"spider","phi":0.0,"theta":0.0,"psi":0.0,"tx":s[3],"ty":s[4],"tz":s[5]})

if myid == main_node:
	ccparlist = [[0.0 for j in xrange(lpar+1)] for i in xrange(ncpu)]

ccparlist0 = [[0.0 for j in xrange(lpar+1)] for i in xrange(nt1)]


if myid == main_node:
	print "                  a_x      a_y      a_z       x        y        z         dot"

for i in xrange(nt1):
	lt = i*ncpu+myid

	for k in xrange(lpar):
		if k<lpar-3: params[k] = shake_ang_1*(2.0*random()-1.0)
		else: params[k] = shake_shift_1*(2.0*random()-1.0)

	params = ali_vol_grid(img1_filt, params, img2_filt, ang_bracket, shift_bracket, mask_rad, discrep, kb, False)

	# params are in the "xyz" convention, so get "spider" ones to do the rot:
	tr = Transform({"type":"xyz","xtilt":params[0],"ytilt":params[1],"ztilt":params[2], "tx":params[3], "ty":params[4], "tz":params[5]})
	qt = tr.get_params("spider")

	x = rot_shift3D_grid(img1_filt,qt['phi'],qt['theta'],qt['psi'],qt['tx'],qt['ty'],qt['tz'],1.0,kb,"background", False)

	tc = compdisc(x, img2_filt, mask, discrep)

	print "     fit #%3d  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %10.6f" % \
	(lt,params[0],params[1],params[2],params[3],params[4],params[5], tc)

	if(init_median):
		ccparlist0[i][0] = tc
		for j in xrange(lpar):
			ccparlist0[i][j+1] = params[j]
	else:
		if(tc>bccc):
			bccc = tc
			for k in xrange(lpar): bparams[k] = params[k]

if(init_median):
	ccparlist0.sort()
	bccc = ccparlist0[nt1/2][0]    # or (nt1-1)/2
	for j in xrange(lpar):
		bparams[j] = ccparlist0[nt1/2][j+1]    # or (nt1-1)/2

mpi_barrier(MPI_COMM_WORLD)

if myid == main_node:
	for j in xrange(ncpu):
		if(j != main_node):
			cc_params = mpi_recv(lpar+1, MPI_FLOAT, j, j*100, MPI_COMM_WORLD)
			if(init_median):
				for k in xrange(lpar+1):
					ccparlist[j][k] = float(cc_params[k])
			else:
				if (float(cc_params[0])>bccc):
					bccc = float(cc_params[0])
					for k in xrange(lpar): bparams[k] = float(cc_params[k+1])
		else:
			if(init_median):
				ccparlist[j][0] = bccc
				for k in xrange(lpar):
					ccparlist[j][k+1] = bparams[k]
else:
	mpi_send([bccc, bparams[0], bparams[1], bparams[2], bparams[3], bparams[4], bparams[5]], lpar+1, MPI_FLOAT, main_node, myid*100, MPI_COMM_WORLD)

mpi_barrier(MPI_COMM_WORLD)

if myid == main_node:
	if(init_median):
		ccparlist.sort()
		bccc = ccparlist[ncpu/2][0]    # or (ncpu-1)/2
		for j in xrange(lpar):
			bparams[j] = ccparlist[ncpu/2][j+1]    # or (ncpu-1)/2
		
	print
	if(init_median):
		print "BEST (median): %7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %10.6f" % \
		(bparams[0],bparams[1],bparams[2],bparams[3],bparams[4],bparams[5], bccc)
	else:
		print "BEST (max):    %7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %10.6f" % \
		(bparams[0],bparams[1],bparams[2],bparams[3],bparams[4],bparams[5], bccc)
	print

mpi_barrier(MPI_COMM_WORLD)

bparams = mpi_bcast(bparams, lpar, MPI_FLOAT, main_node, MPI_COMM_WORLD)


# transform with the resulting parameters:
tsx_init = Transform({"type":"xyz","xtilt":float(bparams[0]),"ytilt":float(bparams[1]),"ztilt":float(bparams[2]),"tx":float(bparams[3]),"ty":float(bparams[4]),"tz":float(bparams[5])})
# transform in the original box:
tsxp = tshift * tsx_init * tshift.inverse()
qt = tsxp.get_params("spider")
# best initial fit in original box:
img1_orig_g,kb0 = prepi3D(img1_orig_0)
x = rot_shift3D_grid(img1_orig_g, qt['phi'],qt['theta'],qt['psi'],qt['tx'],qt['ty'],qt['tz'], 1.0, kb0, "background", False)

# test - write out best initial fit:
if (myid == main_node):
	if nover>1:
		x = Util.decimate(x, nover, nover, nover)
	drop_image(x, segName+"_"+str(resol)+"_init_fit.spi","s")
#mpi_finalize(); exit(777)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# using new position of segment, recompute filter:

mask2 = ellipsoid(x)
um = mask2*img2_orig
pm = mask2*x
fcr = fsc(um, pm)[1]

lf = len(fcr); fcr1 = [0.0]*lf
for i in range(lf):
	if fcr[i]>sqrfsc[i]:
		fcr1[i] = sqrfsc[i]
	else:
		fcr1[i] = fcr[i]

fcrf = fitfcr(fcr1, R,B)

for i in range(lf):
	if fcrf[i]>sqrfsc[i]:
		fcrf[i] = sqrfsc[i]

pwu = rops_table(um)
pwp = rops_table(pm)
Wfilt = create_CCFR_filter(pwu, pwp, fcrf, sg)

if filtseg:
	img1_filt = filt_table(x, Wfilt)
else:
	img1_filt = x.copy()

# test:
#if myid == main_node:
#	filename = segName+"_"+str(resol)+"_fcr.txt"
#	if synthetic and sim_noise:
#		write_text_file([fcr, fcrf, Wfilt, sqrfsc], filename)
#	else:
#		write_text_file([fcr, fcrf, Wfilt, sqrfsc], filename)
#	drop_image(filt_table(img2_orig, Wfilt), segName+"_"+str(resol)+"_map_filt.spi", "s")
#	drop_image(filt_table(x, Wfilt), segName+"_"+str(resol)+"_init_fit_filt.spi", "s")
#mpi_finalize(); exit(999)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# window moved, filtered segment:
img1_filt = Util.window(img1_filt, nbox, nbox, nbox, int(s[3]), int(s[4]), int(s[5]))

# prepare moved segment for gridding:
#(after this point, img1_init MUST be used ONLY in rot_shift3D_grid):
img1_init,kbi = prepi3D(img1_filt)


# create filter to model noise:
pwG = rops_table(img2_orig)
noise_filt = create_noise_filter(pwG, FSC, q2)

# adjust target map so that the given FSC (and SSNR) are in
# exact accordance to the noise filter just computed (which is bounded):
img2_orig = map_adjustment(img2_orig, FSC)


sumbp  = [0.0]*(lpar+1)
sumbp2 = [0.0]*(lpar+1)
bp_mean  = [0.0]*(lpar+1)
bp_sigma = [0.0]*(lpar+1)
bp_skew  = [0.0]*(lpar+1)
bp_delta = [0.0]*(lpar+1)

for j in range(lpar+1):
	sumbp[j]  = 0.0
	sumbp2[j] = 0.0

# create the "cloud" volume(s) (average of moved segments):
nxc = img1_orig_0.get_xsize()
nyc = img1_orig_0.get_ysize()
nzc = img1_orig_0.get_zsize()
cloud = model_blank(nxc, nyc, nzc)
if alpha != 0.0:
	cloud2 = model_blank(nxc, nyc, nzc)


# array to store sigmas at the end of each cycle:
sigmas = [[0.0 for j in xrange(lpar+1)] for i in xrange(max_cycles)]

for cycle in xrange(max_cycles):
	lt = cycle*ncpu+myid
	img2_noisy = create_noisy_map(img2_orig, sigNoise, noise_filt, factor_noise)
	if (not filtseg):
		img2_noisy = filt_table(img2_noisy, Wfilt)
	if imask:
		img2_noisy = img2_noisy * imaskvol
	img2_noisy = Util.window(img2_noisy, nbox, nbox, nbox, int(s[3]), int(s[4]), int(s[5]))

	bccc = -1.e10

	bparams = [0.0]*lpar
	for i in xrange(nt2):
		for k in xrange(lpar):
			if k<lpar-3: params[k] = shake_ang_2*(2.0*random()-1.0)
			else: params[k] = shake_shift_2*(2.0*random()-1.0)

		params = ali_vol_grid(img1_init, params, img2_noisy, ang_bracket, shift_bracket, mask_rad, discrep, kbi, False)

		# params are in the "xyz" convention, so get "spider" ones to do the rot:
		tr = Transform({"type":"xyz","xtilt":params[0],"ytilt":params[1],"ztilt":params[2], "tx":params[3], "ty":params[4], "tz":params[5]})
		qt = tr.get_params("spider")

		x = rot_shift3D_grid(img1_init,qt['phi'],qt['theta'],qt['psi'],qt['tx'],qt['ty'],qt['tz'],1.0,kbi,"background",False)

		tc = compdisc(x, img2_noisy, mask, discrep)

		if(tc>bccc):
			bccc = tc
			for k in xrange(lpar): bparams[k] = params[k]

	print "BEST fit #%3d   %6.2f   %6.2f   %6.2f   %6.2f   %6.2f   %6.2f   %10.6f" % \
	(lt,bparams[0],bparams[1],bparams[2],bparams[3],bparams[4],bparams[5], bccc)

	#if lt<nmatsout:
	parfile = open(segName+str(resol)+"_params_"+str(lt)+".txt","w")
	row = "%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n" % (bparams[0],bparams[1],bparams[2],bparams[3],bparams[4],bparams[5])
	parfile.write(row)
	parfile.close()

	# pixel motion from rotation
	tsx = Transform({"type":"xyz","xtilt":bparams[0],"ytilt":bparams[1],"ztilt":bparams[2], "tx":bparams[3], "ty":bparams[4], "tz":bparams[5]})
	qts = tsx.get_params("spin")
	rotp = red(qts['Omega'])*radius*3.14159/180.

	for j in range(lpar):
		sumbp[j]  += bparams[j]
		sumbp2[j] += bparams[j]**2
	sumbp[lpar]  += rotp
	sumbp2[lpar] += rotp**2

	N = (cycle+1)*ncpu
	bp_mean  = mpi_reduce(sumbp,  lpar+1, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
	bp_sigma = mpi_reduce(sumbp2, lpar+1, MPI_FLOAT, MPI_SUM, main_node, MPI_COMM_WORLD)
	mpi_barrier(MPI_COMM_WORLD)

	if myid == main_node:
		for j in range(lpar+1):
			bp_mean[j] = bp_mean[j]/N
			bp_sigma[j] = sqrt((bp_sigma[j]- N*bp_mean[j]**2)/(N-1))
			sigmas[cycle][j] = bp_sigma[j]

		if cycle>rave:
			maxdev = 0.0
			for j in range(lpar-3,lpar+1):
				ss = 0.0
				for k in xrange(cycle-rave,cycle+1):
					ss += sigmas[k][j]
				devj = abs(sigmas[cycle][j]-sigmas[cycle-rave-1][j])/ss
				if devj > maxdev: maxdev = devj
				bp_sigma[j] = ss/float(rave+1)    # replacing by running average
			if maxdev<sigma_max_change:
				to_break = 1
			else:
				to_break = 0
			#test:
			print "cycle = %2d, sigma_rot = %6.2f, sigma_x = %6.2f, sigma_y = %6.2f, sigma_z = %6.2f, change = %6.2f%%" % \
				(cycle, bp_sigma[lpar], bp_sigma[lpar-3], bp_sigma[lpar-2], bp_sigma[lpar-1], maxdev*100.0)
		else:
			to_break = 0
			#test:
			print "cycle = %2d, sigma_rot = %6.2f, sigma_x = %6.2f, sigma_y = %6.2f, sigma_z = %6.2f" % \
				(cycle, bp_sigma[lpar], bp_sigma[lpar-3], bp_sigma[lpar-2], bp_sigma[lpar-1])
	else:
		to_break = 0
	
	to_break = mpi_bcast(to_break, 1, MPI_INT, main_node, MPI_COMM_WORLD)
	to_break = int(to_break[0])
	if to_break == 1: break


if myid == main_node:
	
	M = min(nmatsout,N)
	
	# collect params for all fits in the array parmat:
	parmat = [[0.0 for j in xrange(lpar)] for i in xrange(N)]

	for i in range(N):
		paramsrows = read_text_row(segName+str(resol)+"_params_"+str(i)+".txt")
		for j in range(lpar):
			parmat[i][j] = paramsrows[0][j]

	# remove all individual parameter files:
	for i in range(N):
		system("rm "+segName+str(resol)+"_params_"+str(i)+".txt")

	# compute some more stats of the fitting parameters:
	median = [0.0]*lpar; parlist = [0.0]*N; cp=[0.0]*lpar; wp=[0.0]*lpar
	for j in range(lpar):
		for i in range(N):
			parlist[i] = parmat[i][j]
		parlist.sort()
		median[j] = parlist[N/2]
		parlist_filt = []
		for i in range(N):
			if(parlist[i]>median[j]-3*bp_sigma[j] and parlist[i]<median[j]+3*bp_sigma[j]):
				parlist_filt.append(parlist[i])

		# get "min/max" of distribution:
		ap = ess_p(parlist_filt,0.05)
		bp = ess_p(parlist_filt,0.95)

		# center and width of distrib:
		cp[j] = (ap+bp)/2.0
		wp[j] = (bp-ap)/2.0
	
		# compute stats without outliers:
		ave,var,skew = mean_var_skew(parlist_filt)
		bp_mean[j]  = ave
		bp_sigma[j] = sqrt(var)
		bp_skew[j]  = skew

	# put together all parameters in one file (after filtering outliers
	# and composing with initial motion):
	parfile = open(segName+"_"+str(resol)+"_params_movie.txt","w")
	parmat_filt = []     # parmat with outliers excluded
	count = 0
	for i in range(N):
		keep = 1
		for j in range(lpar):
			if(parmat[i][j]<median[j]-3*bp_sigma[j] or parmat[i][j]>median[j]+3*bp_sigma[j]):
				keep = 0
				break
		if keep == 0:
			continue
		parmat_filt.append(parmat[i])
		# scale parmat[i] by alpha (if non-zero) wrt to the mean:
		if alpha != 0.0:
			for j in range(lpar):
				parmat[i][j] = (parmat[i][j]-bp_mean[j])*alpha+bp_mean[j]
		# compose with init motion and write out line to parameter file:
		tsx = Transform({"type":"xyz","xtilt":parmat[i][0],"ytilt":parmat[i][1],"ztilt":parmat[i][2],"tx":parmat[i][3],"ty":parmat[i][4],"tz":parmat[i][5]})
		tsx_comp = tsx * tsx_init
		qtc = tsx_comp.get_params("xyz")
		row = "%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n" % (red(qtc['xtilt']), red(qtc['ytilt']), red(qtc['ztilt']), qtc['tx'],qtc['ty'],qtc['tz'])
		parfile.write(row)
		count += 1
		if count >= nmatsout:
			break
	parfile.close()

	M = count   # final number of output fits

	print "          aves: %7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f" % \
		(bp_mean[0], bp_mean[1], bp_mean[2], bp_mean[3], bp_mean[4], bp_mean[5])
	print "        sigmas: %7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f" % \
		(bp_sigma[0], bp_sigma[1], bp_sigma[2], bp_sigma[3], bp_sigma[4], bp_sigma[5])
	print "         skews: %7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f" % \
		(bp_skew[0], bp_skew[1], bp_skew[2], bp_skew[3], bp_skew[4], bp_skew[5])

	sigma_tr  = sqrt(bp_sigma[lpar-3]**2+bp_sigma[lpar-2]**2+bp_sigma[lpar-1]**2)
	sigma_mot = sqrt(sigma_tr**2+bp_sigma[lpar]**2)
	skew_mot  = sqrt(bp_skew[3]**2+bp_skew[4]**2+bp_skew[5]**2)

	print "  sigma_noise =",sigNoise
	print "  sigma_rot = %7.3f pix, sigma_shift = %7.3f pix, sigma_mot = %7.3f pix" % \
		(bp_sigma[lpar], sigma_tr, sigma_mot)
	print "   skew_mot = %7.3f" % (skew_mot)

	# write out mean position:
	tsx_mean = Transform({"type":"xyz","xtilt":float(bp_mean[0]),"ytilt":float(bp_mean[1]),"ztilt":float(bp_mean[2]),"tx":float(bp_mean[3]),"ty":float(bp_mean[4]),"tz":float(bp_mean[5])})
	tsxp = tshift * tsx_mean * tsx_init * tshift.inverse()
	qt = tsxp.get_params("spider")
	x = rot_shift3D(img1_orig_0,qt['phi'],qt['theta'],qt['psi'],qt['tx'],qt['ty'],qt['tz'])
	if nover>1:
		x = Util.decimate(x, nover, nover, nover)
	drop_image(x, segName+"_"+str(resol)+"_meanx.spi","s")

	# ... I think it's probably better to use sigma instead of width to define the ellipsoid:
	for j in range(lpar):
		wp[j] = bp_sigma[j]

	# generate cloud(s):
	for i in range(M):
		for j in range(lpar):
			params[j] = parmat_filt[i][j]
		# project motion onto surface of ellipsoid:
		krot =   sqrt(((params[0]-cp[0])/wp[0])**2 + ((params[1]-cp[1])/wp[1])**2 + ((params[2]-cp[2])/wp[2])**2)
		kshift = sqrt(((params[3]-cp[3])/wp[3])**2 + ((params[4]-cp[4])/wp[4])**2 + ((params[5]-cp[5])/wp[5])**2)
		bp_delta[0] = (params[0]-cp[0])/krot;   bp_delta[1] = (params[1]-cp[1])/krot;
		bp_delta[2] = (params[2]-cp[2])/krot;   bp_delta[3] = (params[3]-cp[3])/kshift;
		bp_delta[4] = (params[4]-cp[4])/kshift; bp_delta[5] = (params[5]-cp[5])/kshift;

		# loop over both signs of each delta to generate 64 points on the ellipsoid:
		for s0 in [-1,1]:
			xte = s0*bp_delta[0]+cp[0]
			for s1 in [-1,1]:
				yte = s1*bp_delta[1]+cp[1]
				for s2 in [-1,1]:
					zte = s2*bp_delta[2]+cp[2]
					for s3 in [-1,1]:
						txe = s3*bp_delta[3]+cp[3]
						for s4 in [-1,1]:
							tye = s4*bp_delta[4]+cp[4]
							for s5 in [-1,1]:
								tze = s5*bp_delta[5]+cp[5]
								tsxe = Transform({"type":"xyz","xtilt":xte,"ytilt":yte,"ztilt":zte,"tx":txe,"ty":tye,"tz":tze})
								# total transformation in the original (non-windowed) volume:
								tsxp = tshift * tsxe * tsx_init * tshift.inverse()
								qt = tsxp.get_params("spider")
								img4 = rot_shift3D(img1_orig_0,qt['phi'],qt['theta'],qt['psi'],qt['tx'],qt['ty'],qt['tz'])
								cloud = 0.5*(square_root(square(cloud-img4))+cloud+img4)   # max(cloud,img4)
								#do scaled cloud if requested:
								if alpha != 0.0:
									xte2 = s0*bp_delta[0]*alpha+cp[0]
									yte2 = s1*bp_delta[1]*alpha+cp[1]
									zte2 = s2*bp_delta[2]*alpha+cp[2]
									txe2 = s3*bp_delta[3]*alpha+cp[3]
									tye2 = s4*bp_delta[4]*alpha+cp[4]
									tze2 = s5*bp_delta[5]*alpha+cp[5]
									tsxe2 = Transform({"type":"xyz","xtilt":xte2,"ytilt":yte2,"ztilt":zte2,"tx":txe2,"ty":tye2,"tz":tze2})
									# total transformation in the original (non-windowed) volume:
									tsxp2 = tshift * tsxe2 * tsx_init * tshift.inverse()
									qt2 = tsxp2.get_params("spider")
									img5 = rot_shift3D(img1_orig_0,qt2['phi'],qt2['theta'],qt2['psi'],qt2['tx'],qt2['ty'],qt2['tz'])
									cloud2 = 0.5*(square_root(square(cloud2-img5))+cloud2+img5)   # max(cloud2,img5)

	if nover>1:
		cloud = Util.decimate(cloud, nover, nover, nover)
		if alpha != 0.0:
			cloud2 = Util.decimate(cloud2, nover, nover, nover)

	drop_image(cloud, segName+"_"+str(resol)+"_cloudx.spi", "s")
	if alpha != 0.0:
		drop_image(cloud2, segName+"_"+str(resol)+"_cloudx_scaled.spi", "s")

mpi_finalize()
exit()
