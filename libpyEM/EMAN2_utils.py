#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

#### python utilities. 
#### 2017-03

from past.utils import old_div
from builtins import range
import numpy as np
import os
from EMAN2 import *
import importlib
import sys
from scipy.interpolate import interp1d

amino_dict= {0: 'ALA', 1: 'ARG', 2: 'ASN', 3: 'ASP', 4: 'CYS', 5: 'GLU', 6: 'GLN', 7: 'GLY', 8: 'HIS', 9: 'ILE', 10: 'LEU', 11: 'LYS', 12: 'MET', 13: 'PHE', 14: 'PRO', 15: 'SER', 16: 'THR', 17: 'TRP', 18: 'TYR', 19: 'VAL', 20: 'ASX', 21:'GLX', 22: 'SEC'}
amino_dict.update(dict((v, k) for k, v in list(amino_dict.items())))
amino_dict.update({'A': 0, 'C': 4, 'E': 5, 'D': 3, 'G': 7, 'F': 13, 'I': 9, 'H': 8, 'K': 11, 'M': 12, 'L': 10, 'N': 2, 'Q': 6, 'P': 14, 'S': 15, 'R': 1, 'T': 16, 'W': 17, 'V': 19, 'Y': 18, 'X':20, 'U': 22})

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
		return old_div(vec,np.sqrt(np.sum(vec**2)))
	else:
		return (old_div(vec.T,np.sqrt(np.sum(vec**2,axis=1)))).T
	
	
def get_fft(img):
	return np.fft.fftshift(np.fft.fftn(np.fft.fftshift(img)))

def get_img(fft):
	return np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(fft)).real)

### interpolate points. Same as np.interp, but the output can have >1 dimension
def interp_points(pts, npt=50, pmin=0., pmax=1.):
    pos=np.append(0,np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1)))
    fun_ax=interp1d(pos, pts.T, fill_value='extrapolate')
    mx=np.max(pos)
    rg=old_div(np.arange(npt,dtype=float),(npt-1))*(pmax-pmin)*mx + pmin*mx
    ax=fun_ax(rg).T
    return ax

#### Distance from a point to a line segment
#### copied from stackoverflow..
def dist_line_point(A, B, P):
	""" segment line AB, point P, where each one is an array([x, y]) """
	if all(A == P) or all(B == P):
		return 0
	if np.arccos(np.dot(old_div((P - A), numpy.linalg.norm(P - A)), old_div((B - A), numpy.linalg.norm(B - A)))) > old_div(np.pi, 2):
		return numpy.linalg.norm(P - A)
	if np.arccos(np.dot(old_div((P - B), numpy.linalg.norm(P - B)), old_div((A - B), numpy.linalg.norm(A - B)))) > old_div(np.pi, 2):
		return numpy.linalg.norm(P - B)
	return old_div(numpy.linalg.norm(np.cross((P-A), (P-B))), numpy.linalg.norm(A - B))

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
	return old_div(ret[n - 1:], n)

#### Line to line distance and angle
def line2line_angle(a0, a1, b0, b1):
	a=a1-a0
	b=b1-b0
	c=b0-a0
	lang=old_div(np.dot(a,b),(norm(a)*norm(b)))
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
    
    
	b2=old_div(bxsz,2)
	ds=old_div(1.0,(apix*bxsz))
	ns=min(int(np.floor(old_div(.25,ds))),old_div(bxsz,2))

	ctfout=np.zeros(b2)
	lbda = old_div(12.2639, np.sqrt(voltage * 1000.0 + 0.97845 * voltage * voltage))

	g1=np.pi/2.0*cs*1.0e7*pow(lbda,3.0);  
	g2=np.pi*lbda*defocus*10000.0;         
	acac=np.arccos(old_div(ampcnt,100.0));                 

	s=np.arange(b2, dtype=float)*ds
	gam=-g1*(s**4)+np.asarray(np.dot(np.asmatrix(g2).T, np.matrix(s**2)))
	ctfout = (np.cos(gam-acac))**2

	return ctfout


def make_missing_wedge(img, wedge=60):

	#img=img.transpose(0,1,2)
	ft=get_fft(img)
	ind=np.indices(ft.shape)-old_div(len(ft),2)
	tanx=np.arctan2(ind[2], ind[0])
	tanx=abs(abs(tanx)-old_div(np.pi,2))< (old_div(wedge,2))/180.*np.pi
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
	if usedegrees: phase *= old_div(np.pi,180.)
	uu = old_div(nx*(u-u.min()),(u.max()-u.min()))-old_div(nx,2.)
	vv = old_div(ny*(v-v.min()),(v.max()-v.min()))-old_div(ny,2.)
	x,y=np.indices((nx,ny))
	xx = x-old_div(nx,2.)
	yy = y-old_div(ny,2.)
	o = np.ones((nx*ny))
	AA = np.multiply(amp.ravel()[:,np.newaxis],o[np.newaxis,:])
	pp = np.multiply(phase.ravel()[:,np.newaxis],o[np.newaxis,:])
	uuxx = np.multiply(uu.ravel()[:,np.newaxis],xx.ravel()[np.newaxis,:])
	vvyy = np.multiply(vv.ravel()[:,np.newaxis],yy.ravel()[np.newaxis,:])
	return np.sum(np.real(AA*np.exp(2*np.pi*1j*(uuxx+vvyy)+pp)).reshape(len(u),nx,ny),axis=0)


def make_path(suffix):
	### make a suffix_xx folder and return the folder name
	for i in range(100):
		path="{}_{:02d}/".format(suffix, i)
		if os.path.exists(path):
			continue
		else:
			os.mkdir(path)
			break
	else:
		print("Too many {} folders in the project, or something odd happened....Exit.".format(suffix))
		exit()
		
	return path


def findfs(stem=''):
	"""
	"Returns a sorted list with the files in the current directory that contain the string(s) indicated by 'id'
	To find files with more than one string, use *. 
	For example stem=id1*id2 will find files with "id1" and "id2" in them, regardless of where these strings occur in the filename
	Author: Jesus Montoya, jgalaz@gmail.com, April 2019
	"""
	findir=set(fsincurrentdir())
	
	if stem: #only do this if a string has been passed in
		ids=[stem]
		if '*' in stem: ids=stem.split('*')

		for i in ids:
			if i: #ignore empty strings
				selectfiles=set([f for f in findir if i in f])
				findir=findir.intersection(selectfiles)

	findir=list(findir)
	findir.sort()
	return findir


def fsincurrentdir():
	"""Returns a sorted list with the files in the current directory
	Author: Jesus Montoya, jgalaz@gmail.com, April 2019
	"""
	import os
	c = os.getcwd()
	findir = os.listdir(c)
	findir.sort()
	return findir


def fisindir(f):
	"""Checks whether a file is in the current directory
	Author: Jesus Montoya, jgalaz@gmail.com, April 2019
	"""
	fs=fsincurrentdir()
	result=False
	if f in fs:
		result=True
	return result


def fileisimage(f):
	"""Checks whether a file 'f' is an image readable my EMAN2
	Author: Jesus Montoya, jgalaz@gmail.com, April 2019
	"""
	isimage=False
	try:
		a=EMData(f,0,True)
		isimage=True
	except:
		print("\n(EMAN2_utils)(fileisimage) file={} is NOT a valid image file".format(f))
	return isimage


def cleanfilenames():
	"""Renames all files in a directory and its subdirectories to NOT contain parentheses, brackets, commas, spaces
	Author: Jesus Montoya, jgalaz@gmail.com, April 2019
	"""
	fs=fsincurrentdir()
	badcharacters=['[',']','{','}','(',')',',','  ',' ']
	cmds = ["""find . -depth -name '*"""+b+"""*' -execdir bash -c 'for f; do mv -i "$f" "${f//"""+b+"""/_}"; done' bash {} +""" for b in badcharacters]
	for cmd in cmds:
		runcmdbasic(cmd)
	return	


def makepath(options, stem='e2dir'):
	"""
	Makes a numbered series of subdirectories to compartmentalize results 
	when the same programs are run in the same parent directory
	Author: Jesus Montoya, jgalaz@gmail.com
	"""
	if not options.path:
		if options.verbose:
			print("\n(EMAN2_utils)(makepath), stem={}".format(stem))
	
	elif options.path:
		stem=options.path
	
	i=1
	while os.path.exists("{}_{:02d}".format(stem,i)): i+=1
	
	options.path="{}_{:02d}".format(stem,i)
	try: 
		os.mkdir(options.path)
	except: 
		pass
	
	return options

#### this triggers a segmentation fault. in case we need it somewhere...
def do_segfault():
	exec'()'*7**6

def detectThreads(options):
	"""
	c:If parallelism isn't set, parallelize automatically unless disabled
	Author: Jesus Montoya, jgalaz@gmail.com
	"""
	import multiprocessing
	nparallel = multiprocessing.cpu_count()

	try:	
		if options.parallel and options.parallel != 'None' and options.parallel != 'none' and options.parallel != 'NONE':
			if 'mpi' not in options.parallel:

				options.parallel = 'thread:' + str(nparallel)
				print("\nfound cores n={}".format(nparallel))
				print("setting --parallel={}".format(options.parallel))
			else:
				print("\n(EMAN2_utils)(detectThreads) nothing to do; mpi paralellism specified; options.parallel={}".format(otpions.parallel)) 
	
		elif not options.parallel or options.parallel == 'None' or options.parallel == 'none':
			options.parallel = None
			print("\n(EMAN2_utils)(detectThreads) WARNING: parallelism disabled, options.parallel={}".format(options.parallel) )
	except:
		try:
			if options.threads and options.threads != 'None' and options.threads != 'none' and options.threads != 'NONE':
				
				options.threads = 'thread:' + str(nparallel)
				print("\nfound cores n={}".format(nparallel))
				print("setting --parallel={}".format(options.parallel))
			elif not options.threads or options.threads == 'None' or options.threads == 'none':
				options.parallel = None
				print("\n(EMAN2_utils)(detectThreads) WARNING: parallelism disabled, options.parallel={}".format(options.parallel) ) 
		except:
			print("\n(EMAN2_utils)(detectThreads) WARNING: No parallelism option detected, neither --parallel nor --threads.")

	return options


def checkinput(options):
	"""
	Checks for sanity of input whether directly as an argument or through --input. Both should be functional.
	Author: Jesus Montoya, jgalaz@gmail.com
	"""
	if not options.input:
		try:
			options.input = sys.argv[1]
			print("\ntrying to read input from sys.argv[1]={}".format(options.input))
			EMData(options.input,0,True)
		except:
			print("\n(EMAN2_utils)(checkinput) ERROR: input file {} seems to have an invalid format or doesn't exist; verify that the filename is correct.".format( options.input ))
			#parser.print_help()
			sys.exit(1)
	else:
		try:
			print("\ntrying to read input from --input={}".format(options.input))
			EMData(options.input,0,True)
		except:
			print("\n(EMAN2_utils)(checkinput) ERROR: --input file {} seems to have an invalid format or doesn't exist; verify that the filename is correct.".format( options.input ))
			sys.exit(1)
	return options
	

def runcmd(options,cmd,cmdsfilepath=''):
	"""
	Version of runcmd (below) with verbose feedback and option to save the executed command to a commands file for record keeping
	Author: Jesus Montoya, jgalaz@gmail.com
	"""
	if options.verbose > 9:
		print("\n(EMAN2_utils)(runcmd) running command {}".format(cmd))
	
	runcmdbasic(cmd)
	
	if cmdsfilepath:
		with open(cmdsfilepath,'a') as cmdfile: cmdfile.write( cmd + '\n')

	if options.verbose > 8:
		print("\n(EMAN2_utils)(runcmd) done")

	return 1


def runcmdbasic(cmd):
	"""
	Runs commands "properly" at the command line, April 2019
	Author: Jesus Montoya, jgalaz@gmail.com
	"""
	p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	text=p.communicate()	
	p.stdout.close()

	return


def origin2zero(img):
	"""
	sets the 'origin_' parameters in image headers to zero
	Author: Jesus Montoya, jgalaz@gmail.com
	"""
	#print("\n(EMAN2_utils)(origin2zero)")

	img['origin_x'] = 0
	img['origin_y'] = 0
	try:
		img['origin_z'] = 0
	except:
		pass
	return img


def clip3d( vol, size, center=None ):
	"""
	3d clipping function to avoid having to define a region in your code when you know the center and size of the box you want
	Author: Jesus Montoya, jgalaz@gmail.com
	"""

	#if options.verbose:
	#	print("\n(EMAN2_utils)(clip3d) starting")
	
	volxc = old_div(vol['nx'],2)
	volyc = old_div(vol['ny'],2)
	volzc = old_div(vol['nz'],2)
	
	if center:
		volxc = center[0]
		volyc = center[1]
		volzc = center[2]
	
	Rvol =  Region( old_div((2*volxc - size),2), old_div((2*volyc - size),2), old_div((2*volzc - size),2), size , size , size)
	volclip = vol.get_clip( Rvol )
	
	#if options.verbose:
	#	print("\n(EMAN2_utils)(clip3d) done")
	
	return volclip


def clip2d( img, size, center=None ):
	"""
	2d clipping function to avoid having to define a region in your code when you know the center and size of the box you want
	Author: Jesus Montoya, jgalaz@gmail.com
	"""

	#if options.verbose:
	#	print("\n(EMAN2_utils)(clip2d)")

	imgxc = old_div(img['nx'],2)
	imgyc = old_div(img['ny'],2)
	
	if center:
		imgxc = center[0]
		imgyc = center[1]
	
	Rimg = Region( old_div((2*imgxc - size),2), old_div((2*imgyc - size),2), size , size )
	imgclip = img.get_clip( Rimg )
	
	return imgclip


def textwriter(data,options,name,invert=0,xvals=None,onlydata=False):
	"""
	writes a list of values into a double column text file with rows of the form: "index value"
	E.g., 0 10.0 (first row), 1 21.2 (second row), 1 -18.2,..., N valueN
	This aims to make a file from which 'data' can be easily plotted with any other program
	Author: Jesus Montoya, jgalaz@gmail.com
	"""

	try:
		if options.path:
			if options.path not in name:
				name = options.path + '/' + name
	except:
		pass

	if options.verbose:
		print("(EMAN2_utils)(textwriter) writing the following file {}".format(name))
	
	with open(name,'w') as f:
		lines=[]
		for i in range(len(data)):
			val=data[i]
			if invert:
				val*=-1
				
			line2write = str(i) + '\t' + str(val) + '\n'
			if xvals:
				line2write = str(xvals[i]) + '\t' + str(val) + '\n'
			elif onlydata:
				line2write = str(val) + '\n'

			#print "THe line to write is"
			if line2write.replace('\n','').replace(' ','') != '':
				lines.append(line2write)
		

		#remove empty lines that somehow creep into files
		print('\n\n\n(EMAN2_utils)(textwriter) REMOVING EMPTY LINES')
		newlines=[]
		for line in lines:
			newline = line.replace('\n','').replace(' ','')
			if newline:
				newlines.append(newline+'\n')

		f.writelines(newlines)
	#f.close()

	return


def cmponetomany(reflist,target,align=None,alicmp=("dot",{}),cmp=("dot",{}), ralign=None, alircmp=("dot",{}),shrink=None,mask=None,subset=None,prefilt=False,verbose=0):
	"""Compares one image (target) to a list of many images (reflist). Returns """

	ret=[None for i in reflist]
#	target.write_image("dbug.hdf",-1)
	for i,r in enumerate(reflist):
		#print i,r
		if r[0]["sigma"]==0 : continue				# bad reference
		if subset!=None and i not in subset :
			ret[i]=None
			continue
		if prefilt :
			msk=r[0].process("threshold.notzero")					# mask from the projection
			r[0].process_inplace("filter.matchto",{"to":target})
			r[0].mult(msk)											# remask after filtering

#		print "Final: ",target["source_n"],",",r[0]["source_n"]

		if align[0] :
			r[0].del_attr("xform.align2d")
			ta=r[0].align(align[0],target,align[1],alicmp[0],alicmp[1])
			if verbose>3: print(ta.get_attr("xform.align2d"))
			#ta.debug_print_params()

			if ralign and ralign[0]:
				if r[1]!=None :
					#print "(single) using mask, and ",mask
					ralign[1]["xform.align2d"] = ta.get_attr("xform.align2d").inverse()
					r[0].del_attr("xform.align2d")
					ralign[1]["mask"]=r[1]
					alip = target.align(ralign[0],r[0],ralign[1],alircmp[0],alircmp[1])
					ta=r[0].copy()
					ta.transform(alip["xform.align2d"].inverse())
					ta["xform.align2d"]=alip["xform.align2d"].inverse()
				else:
					ralign[1]["xform.align2d"] = ta.get_attr("xform.align2d")
					r[0].del_attr("xform.align2d")
					ta = r[0].align(ralign[0],target,ralign[1],alircmp[0],alircmp[1])

				if verbose>3: print(ta.get_attr("xform.align2d"))


			t =  ta.get_attr("xform.align2d")
			t.invert()
			p = t.get_params("2d")

			scale_correction = 1.0
			if shrink != None: scale_correction = float(shrink)

			if mask!=None :
				ta.mult(mask)
				ptcl2=target.copy()
				ptcl2.mult(mask)
				ret[i]=(ptcl2.cmp(cmp[0],ta,cmp[1]),scale_correction*p["tx"],scale_correction*p["ty"],p["alpha"],p["mirror"],p["scale"])
			else:
				try:
					ret[i]=(target.cmp(cmp[0],ta,cmp[1]),scale_correction*p["tx"],scale_correction*p["ty"],p["alpha"],p["mirror"],p["scale"])
				except:
					print("ERROR: CMP FAILURE. See err.hdf")
					print(cmp)
					target.write_image("err.hdf",0)
					ta.write_image("err.hdf",1)
					sys.exit(1)
					
#			ta.write_image("dbug.hdf",-1)

#				print ta["source_n"],target["source_n"]
				#psub=target.process("math.sub.optimal",{"ref":ta})
				#nout=ta["source_n"]*3
				#ta.write_image("dbug_%d.hdf"%target["source_n"],nout)
				#target.write_image("dbug_%d.hdf"%target["source_n"],nout+1)
				#psub.write_image("dbug_%d.hdf"%target["source_n"],nout+2)


		else :
			ret[i]=(target.cmp(cmp[0],r[0],cmp[1]),0,0,0,1.0,False)

		if verbose==3 : print(ret[i][0], end=' ')

	if verbose==3 : print("")
	if verbose==2 :
		print("Best: ",sorted([(ret[i][0],i) for i in range(len(ret))])[0])
	return ret


def sptOptionsParser( options, program='' ):
	"""
	Used by some SPT programs (nov/2017); function might be deprecated or refactored in the near future	
	Author: Jesus Montoya, jgalaz@gmail.com
	"""
	
	print("\n(EMAN2_utils)(sptOptionsParser) parsing options")
	if program:
		print("from program {}".format(program))
	
	try:
		if options.align:
			#print "(e2spt_classaverage) --align to parse is", options.align
			options.align=parsemodopt(options.align)
		elif options.align == 'None' or  options.align == 'none':
			options.align=None
	except:
		if options.verbose > 9:
			print("\nWARNING (might not be relevant): --align not provided")
		
	try:
		if options.falign and options.falign != None and options.falign != 'None' and options.falign != 'none': 
			options.falign=parsemodopt(options.falign)
		elif options.falign == 'None' or  options.falign == 'none':
			options.falign=None
	except:
		if options.verbose > 9:
			print("\nWARNING (might not be relevant): --falign not provided")
	
	try:
		if options.aligncmp: 
			options.aligncmp=parsemodopt(options.aligncmp)
		elif options.aligncmp == 'None' or  options.aligncmp == 'none':
			options.aligncmp=None
	except:
		if options.verbose > 9:
			print("\nWARNING (might not be relevant): --aligncmp not provided")
	
	try:	
		if options.faligncmp: 
			options.faligncmp=parsemodopt(options.faligncmp)
		elif options.faligncmp == 'None' or  options.faligncmp == 'none':
			options.faligncmp=None
	except:
		if options.verbose > 9:
			print("\nWARNING (might not be relevant): --faligncmp not provided")
		
	try:
		if options.averager: 
			options.averager=parsemodopt(options.averager)
		elif options.averager == 'None' or  options.averager == 'none':
			options.averager=None
	except:
		if options.verbose > 9:
			print("\nWARNING (might not be relevant): --averager not provided")
		
	try:
		if options.autocenter:
			options.autocenter=parsemodopt(options.autocenter)
		elif options.autocenter == 'None' or  options.autocenter == 'none':
			options.autocenter=None
	except:
		if options.verbose > 9:
			print("\nWARNING (might not be relevant): --autocenter not provided")
		
	try:
		if options.autocentermask:
			options.autocentermask=parsemodopt(options.autocentermask)
		elif options.autocentermask == 'None' or  options.autocentermask == 'none':
			options.autocentermask=None
	except:
		if options.verbose > 9:
			print("\nWARNING (might not be relevant): --autocentermask not provided")
	
	try:
		if options.normproc and options.normproc != 'None' and options.normproc != 'none':
			options.normproc=parsemodopt(options.normproc)
		elif options.normproc == 'None' or  options.normproc == 'none':
			options.normproc=None
	except:
		if options.verbose > 9:
			print("\nWARNING (might not be relevant): --normproc not provided")
	
	try:
		if options.mask and options.mask != 'None' and options.mask != 'none':
			#print "\nparsing mask"
			#print "before = ".format(options.mask)
			options.mask = parsemodopt(options.mask)
			#print "after = ".format(options.mask)
		elif options.mask == 'None' or  options.mask == 'none':
			options.mask=None
	except:
		if options.verbose > 9:
			print("\nWARNING (might not be relevant): --mask not provided")
	
	try:	
		if options.preprocess and options.preprocess != 'None' and options.preprocess != 'none': 
			options.preprocess=parsemodopt(options.preprocess)
		elif options.preprocess == 'None' or  options.preprocess == 'none':
			options.preprocess=None
	except:
		if options.verbose > 9:
			print("\nWARNING (might not be relevant): --preprocess not provided")
	
	try:	
		if options.threshold and options.threshold != 'None' and options.threshold != 'none': 
			options.threshold=parsemodopt(options.threshold)
		elif options.threshold == 'None' or  options.threshold == 'none':
			options.threshold=None
	except:
		if options.verbose > 9:
			print("\nWARNING (might not be relevant): --threshold not provided")
	
	try:
		if options.preprocessfine and options.preprocessfine != 'None' and options.preprocessfine != 'none': 
			options.preprocessfine=parsemodopt(options.preprocessfine)
		elif options.preprocessfine == 'None' or  options.preprocessfine == 'none':
			options.preprocessfine=None
	except:
		if options.verbose > 9:
			print("\nWARNING (might not be relevant): --preprocessfine not provided")
	
	try:	
		if options.lowpass and options.lowpass != 'None' and options.lowpass != 'none': 
			options.lowpass=parsemodopt(options.lowpass)
		elif options.lowpass == 'None' or  options.lowpass == 'none':
			options.lowpass=None
	except:
		if options.verbose > 9:
			print("\nWARNING (might not be relevant): --lowpass not provided")
	
	try:
		if options.lowpassfine and options.lowpassfine != 'None' and options.lowpassfine != 'none': 
			options.lowpassfine=parsemodopt(options.lowpassfine)
		elif options.lowpassfine == 'None' or  options.lowpassfine == 'none':
			options.lowpassfine=None
	except:
		if options.verbose > 9:
			print("\nWARNING (might not be relevant): --lowpassfine not provided")
	
	try:
		if options.highpass and options.highpass != 'None' and options.highpass != 'none': 
			options.highpass=parsemodopt(options.highpass)
		elif options.highpass == 'None' or  options.highpass == 'none':
			options.highpass=None
	except:
		if options.verbose > 9:
			print("\nWARNING (might not be relevant): --highpass not provided")
	
	try:
		if options.highpassfine and options.highpassfine != 'None' and options.highpassfine != 'none': 
			options.highpassfine=parsemodopt(options.highpassfine)
		elif options.highpassfine == 'None' or  options.highpassfine == 'none':
			options.highpassfine=None
	except:
		if options.verbose > 9:
			print("\nWARNING (might not be relevant): --highpassfine not provided")
	try:
		if options.postprocess and options.postprocess != 'None' and options.postprocess != 'none': 
			options.postprocess=parsemodopt(options.postprocess)
		elif options.postprocess == 'None' or  options.postprocess == 'none':
			options.postprocess=None
	except:
		if options.verbose > 9:
			print("\nWARNING (might not be relevant): --postprocess not provided")
	
	try:
		if options.reconstructor and options.reconstructor != 'None' and options.reconstructor != 'none': 
			options.reconstructor=parsemodopt(options.reconstructor)
		elif options.reconstructor == 'None' or  options.reconstructor == 'none':
			options.reconstructor=None
	except:
		if options.verbose > 9:
			print("\nWARNING (might not be relevant): --reconstructor not provided")
	
	try:
		if options.preavgproc1 and options.preavgproc1 != 'None' and options.preavgproc1 != 'none': 
			options.preavgproc1=parsemodopt(options.preavgproc1)
		elif options.preavgproc1 == 'None' or  options.preavgproc1 == 'none':
			options.preavgproc1=None
	except:
		if options.verbose > 9:
			print("\nWARNING (might not be relevant): --reconstructor not provided")
		
	try:
		if options.preavgproc2 and options.preavgproc2 != 'None' and options.preavgproc2 != 'none': 
			options.preavgproc2=parsemodopt(options.preavgproc2)
		elif options.preavgproc2 == 'None' or  options.preavgproc2 == 'none':
			options.preavgproc2=None
	except:
		if options.verbose > 9:
			print("\nWARNING (might not be relevant): --reconstructor not provided")
	
	return options



def writeParameters( options, program, tag ):
	'''
	Used by many SPT programs. Function to write the parameters used for every run of the program to parameters.txt inside the path specified by --path.
	Unfortunately, the usability of the .eman2log.txt file is limited when it is overcrowded with commands; e.g., a program that iteratively runs other EMAN2 programs at the command line
	will SWARM the log file with commands that will obscure the command you wanted to log. Having a parameters file explicitly record what was the state of every parameter used by the program
	is useful, as it also explicitly records values for parameters that were used by DEFAULT and not set by the user at the commandline.

	Author: Jesus Montoya, jgalaz@gmail.com
	'''
	import datetime

	print("Tag received in writeParameters is {}".format(tag))

	names = dir(options)
	
	cmd = program
	lines = []
	now = datetime.datetime.now()
	lines.append(str(now)+'\n')
	
	#print "\nnames are", names
	optionscopy = options
	
	try:
		if options.search == 0 or options.search == 0.0:
			options.search = '0'
	except:
		pass
	try:
		if options.searchfine == 0 or options.searchfine == 0.0:
			options.searchfine = '0'
	except:
		pass
		
	#print "mask in write parameters is", optionscopy.mask, type(optionscopy.mask)
	for name in names:
				
		if getattr(options,name) and "__" not in name and "_" not in name:
		#if "__" not in name and "_" not in name:	
	
			#if "__" not in name and "_" not in name and str(getattr(options,name)) and 'path' not in name and str(getattr(options,name)) != 'False' and str(getattr(options,name)) != 'True' and str(getattr(options,name)) != 'None':			
			line = name + '=' + str(getattr(optionscopy,name))
					
			lines.append(line+'\n')
			
			if str(getattr(optionscopy,name)) != 'True' and str(getattr(optionscopy,name)) != 'False' and str(getattr(optionscopy,name)) != '':
			
				if name != 'parallel':
					if "{" in str( getattr(optionscopy,name) ) or "}" in  str(getattr(optionscopy,name)) or ")" in  str(getattr(optionscopy,name)) or ")"  in str(getattr(optionscopy,name)): 
						
						tail = str( getattr(optionscopy,name) ).replace(':','=').replace('(','').replace(')','').replace('{','').replace('}','').replace(',',':').replace(' ','').replace("'",'')
						if tail[-1] == ':':
							tail = tail[:-1] 
						cmd += ' --' + name + '=' + tail
					else:
						
						tail = str( getattr(optionscopy,name) )
						if tail[-1] == ':':
							tail = tail[:-1]
						cmd += ' --' + name + '=' + tail
						
				else:
					cmd += ' --' + name + '=' + str(getattr(optionscopy,name))
			
			elif str(getattr(optionscopy,name)) == 'True' or str(getattr(optionscopy,name)) == 'False':
				cmd += ' --' + name
	
	parmFile = 'parameters_' + tag + '.txt'
	lines.append('\n'+cmd+'\n')
	#f=open( optionscopy.path + '/' + parmFile,'w')
	pfile = optionscopy.path + '/' + parmFile
	f = open( pfile, 'w')
	f.writelines(lines)
	f.close()
	
	return cmd
