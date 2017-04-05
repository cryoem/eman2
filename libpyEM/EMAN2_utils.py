#!/usr/bin/env python

#### python utilities. 
#### 2017-03

import numpy as np

amino_dict= {0: 'ALA', 1: 'ARG', 2: 'ASN', 3: 'ASP', 4: 'CYS', 5: 'GLU', 6: 'GLN', 7: 'GLY', 8: 'HIS', 9: 'ILE', 10: 'LEU', 11: 'LYS', 12: 'MET', 13: 'PHE', 14: 'PRO', 15: 'SER', 16: 'THR', 17: 'TRP', 18: 'TYR', 19: 'VAL', 20: 'ASX', 21:'GLX'}
amino_dict.update(dict((v, k) for k, v in amino_dict.iteritems()))
amino_dict.update({'A': 0, 'C': 4, 'E': 5, 'D': 3, 'G': 7, 'F': 13, 'I': 9, 'H': 8, 'K': 11, 'M': 12, 'L': 10, 'N': 2, 'Q': 6, 'P': 14, 'S': 15, 'R': 1, 'T': 16, 'W': 17, 'V': 19, 'Y': 18, 'X':20})

def pdb2numpy(fname, readres=False, readocc=False, readbfac=False):
	f=open(fname,'r')
	lines=f.readlines()
	f.close()
	data=[]
	for l in lines:
		if l.startswith("ATOM") or l.startswith("HETATM"):
			if l[13:15]!="CA": continue
			atom=[l[30:38],l[38:46],l[46:54]]
			a=[float(a) for a in atom]
			if readres:
				#print l[17:20].strip()
				a.append(amino_dict[l[17:20].strip()])
			if readocc:
				a.append(float(l[54:60].strip()))
			if readbfac:
				a.append(float(l[60:66].strip()))
			data.append(a)
	
	pts=np.array(data)
	return pts

def numpy2pdb(data,fname,occ=[],bfac=[],chainid=[], model=0, residue=[]):
	if model>0:
		ww='a'
		#print "Appending.."
	else:
		ww='w'

	f=open(fname,ww)
	f.write("MODEL     %4d\n"%model)
	if len(occ)!=len(data):
		if len(occ)>0: print ("warning: occ and data not same size!")
		occ=np.zeros(len(data))
	if len(bfac)!=len(data):
		if len(bfac)>0: print ("warning: bfac and data not same size!")
		bfac=np.zeros(len(data)) 
	if len(chainid)!=len(data):
		if len(chainid)>0: print ("warning: chainid and data not same size!")
		chainid=np.zeros(len(data)) 
	if len(residue)!=len(data):
		if len(residue)>0: print ("warning: residue and data not same size!")
		residue=np.zeros(len(data)) 
	atomid=1
	curchain=chainid[0]
	for i,d in enumerate(data):
		if chainid[i]!=curchain: atomid=0
		f.write("ATOM {atomid:6d}  CA  {res} {chainid}{atomid:4d}    {px:8.3f}{py:8.3f}{pz:8.3f}{occ:6.2f}{bfac:6.2f}     S_00  0\n".format(atomid=atomid, chainid=chr(int(chainid[i])+65), px=d[0], py=d[1], pz=d[2], occ=occ[i], bfac=bfac[i], res=amino_dict[residue[i]]))
		atomid+=1
		curchain=chainid[i]

	f.write("TER  {:6d}      ALA {}{:4d}\n""".format(i+1, 'A', i))
	f.write("ENDMDL\n")
	f.close()

def norm_vec(vec):
	if len(vec.shape)==1:
		return vec/np.sqrt(np.sum(vec**2))
	else:
		return (vec.T/np.sqrt(np.sum(vec**2,axis=1))).T
	
	
def get_fft(img):
	return np.fft.fftshift(np.fft.fftn(np.fft.fftshift(img)))

def get_img(fft):
	return np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(fft)).real)

### interpolate points. Same as np.interp, but the output can have >1 dimension
def interp_points(pts, npt=50, pmin=0., pmax=1.):
    pos=np.append(0,np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1)))
    fun_ax=interp1d(pos, pts.T, fill_value='extrapolate')
    mx=np.max(pos)
    rg=np.arange(npt,dtype=float)/(npt-1)*(pmax-pmin)*mx + pmin*mx
    ax=fun_ax(rg).T
    return ax

#### Distance from a point to a line segment
#### copied from stackoverflow..
def dist_line_point(A, B, P):
	""" segment line AB, point P, where each one is an array([x, y]) """
	if all(A == P) or all(B == P):
		return 0
	if np.arccos(np.dot((P - A) / numpy.linalg.norm(P - A), (B - A) / numpy.linalg.norm(B - A))) > np.pi / 2:
		return numpy.linalg.norm(P - A)
	if np.arccos(np.dot((P - B) / numpy.linalg.norm(P - B), (A - B) / numpy.linalg.norm(A - B))) > np.pi / 2:
		return numpy.linalg.norm(P - B)
	return numpy.linalg.norm(np.cross((P-A), (P-B))) / numpy.linalg.norm(A - B)

#### Distance from a set of points to a set of lines
#### Return the distance to the nearest line for each point
def dist_pts_lines(pts, lines):
	dsts=np.zeros((len(pts), len(lines)-1))
	for i,p in enumerate(pts):
		for j in range(len(lines)-1):
			dsts[i,j]=dist_line_point(lines[j], lines[j+1], p)
	return np.min(dsts, axis=1)

#### Moving average
def moving_average(a, n=3) :
	ret = np.cumsum(a, axis=0)
	ret[n:] = ret[n:] - ret[:-n]
	return ret[n - 1:] / n

#### Line to line distance and angle
def line2line_angle(a0, a1, b0, b1):
	a=a1-a0
	b=b1-b0
	c=b0-a0
	lang=np.dot(a,b)/(norm(a)*norm(b))
	return lang

	
#### Calculate the rotation matrix that rotate a given vector to [0,0,1]
def calc_rot_mat(v):

	tt=np.arctan2(v[2],v[1])
	rota=np.array([[1,0,0], [0, np.cos(tt), -np.sin(tt)], [0,np.sin(tt), np.cos(tt)]])
	vr=np.dot(v, rota)
	aa=np.arctan2(vr[1], vr[0])
	rotb=np.array([[np.cos(aa), -np.sin(aa), 0], [np.sin(aa), np.cos(aa),0],[0, 0, 1]])
	rot=np.dot(rota, rotb)
	m0=np.array([[0,0,1],[0,1,0],[1,0,0]])
	rot=np.dot(rot,m0)
	return rot

#### numpy version of EMAN2Ctf.compute_1d(). Takes vector of defocus input and output a matrix of CTF curves
def calc_ctf(defocus, bxsz=256, voltage=300, cs=4.7, apix=1. ,ampcnt=0.):
    
    
	b2=bxsz/2
	ds=1.0/(apix*bxsz)
	ns=min(int(floor(.25/ds)),bxsz/2)

	ctfout=np.zeros(b2)
	lbda = 12.2639 / np.sqrt(voltage * 1000.0 + 0.97845 * voltage * voltage)

	g1=np.pi/2.0*cs*1.0e7*pow(lbda,3.0);  
	g2=np.pi*lbda*defocus*10000.0;         
	acac=np.arccos(ampcnt/100.0);                 

	s=np.arange(b2, dtype=float)*ds
	gam=-g1*(s**4)+np.asarray(np.dot(np.asmatrix(g2).T, np.matrix(s**2)))
	ctfout = (np.cos(gam-acac))**2

	return ctfout


def make_missing_wedge(img, wedge=60):

	#img=img.transpose(0,1,2)
	ft=get_fft(img)
	ind=np.indices(ft.shape)-len(ft)/2
	tanx=np.arctan2(ind[2], ind[0])
	tanx=abs(abs(tanx)-np.pi/2)< (wedge/2)/180.*np.pi
	img2=get_img(ft*tanx)
	
	return img2


def idfft2(v,u,amp,phase,nx=256,ny=256,dtype=np.float32,usedegrees=False):
	"""
	Perform a vectorized, 2D discrete fourier transform. 
	Note that this approach scales poorly with box size.

	Author: Michael Bell (jmbell@bcm.edu)
	"""
	u = np.asarray(u).astype(dtype)
	v = np.asarray(v).astype(dtype)
	amp = np.asarray(amp).astype(dtype)
	phase = np.asarray(phase).astype(dtype)
	if usedegrees: phase *= np.pi/180.
	uu = nx*(u-u.min())/(u.max()-u.min())-nx/2.
	vv = ny*(v-v.min())/(v.max()-v.min())-ny/2.
	x,y=np.indices((nx,ny))
	xx = x-nx/2.
	yy = y-ny/2.
	o = np.ones((nx*ny))
	AA = np.multiply(amp.ravel()[:,np.newaxis],o[np.newaxis,:])
	pp = np.multiply(phase.ravel()[:,np.newaxis],o[np.newaxis,:])
	uuxx = np.multiply(uu.ravel()[:,np.newaxis],xx.ravel()[np.newaxis,:])
	vvyy = np.multiply(vv.ravel()[:,np.newaxis],yy.ravel()[np.newaxis,:])
	return np.sum(np.real(AA*np.exp(2*np.pi*1j*(uuxx+vvyy)+pp)).reshape(len(u),nx,ny),axis=0)
