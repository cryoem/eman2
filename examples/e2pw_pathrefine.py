#!/usr/bin/env python
# Muyuan Chen 2015-03
from EMAN2 import *
import numpy as np
from cmy_utils import *
import scipy.spatial.distance as scidist
from time import time
from scipy.interpolate import interp1d
from Bio import pairwise2

mwid="GLY,ALA,SER,PRO,VAL,THR,CYS,ILE,LEU,ASN,ASP,GLU,LYS,GLN,MET,HIS,PHE,ARG,TYR,TRP,ASX".split(',')
mw=[0,14.027,30.026,40.065,42.081,44.052,46.087,56.108,56.108,57.052,58.037,71.079,71.122,72.064,74.141,80.089,90.125,99.136,106.124,129.161,132.5]
molwt=np.zeros(len(mwid))
for i in range(len(mwid)):
	molwt[amino_dict[mwid[i]]]=mw[i]

def group_pts(pts, dst0, d0):
	group=np.zeros(len(pts))-1
	dst=dst0.copy()
	rid=[[i] for i in np.argmin(d0, axis=0)]
	for i,r in enumerate(rid):
		group[r]=i
		dst[:,r]=np.inf

	for kk in range(100):
		nochange=True
		for i,r in enumerate(rid):
			mm=np.min(dst[r])
			if mm>1: continue
			ii=np.where(dst[r]==mm)[1]
			nochange=False
			group[ii]=i
			dst[:,ii]=np.inf
			rid[i].extend(ii)
		ss= np.sum(group<0)
		#print ss
		if ss==0: break
		if nochange: break
	return group



def read_fasta(fname):
	f=open(fname, 'r')
	lines=f.readlines()
	f.close()
	res=[]
	for l in lines:
		if l.startswith(">"):
			continue
		res.extend(l.strip())
	print "{:d} residues read from sequence".format(len(res))
	return [amino_dict[r] for r in res]
	


def numpy2pdb_sc(data,fname,occ=[],bfac=[],chainid=[], model=0, residue=[],sidechain=[]):
	if model>0:
		ww='a'
		#print "Appending.."
	else:
		ww='w'

	f=open(fname,ww)
	f.write("MODEL	 %4d\n"%model)
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
	resid=1
	if len(sidechain)>0:
		sidechain=np.array(sidechain)
	
	curchain=chainid[0]
	for i,d in enumerate(data):
		if chainid[i]!=curchain: atomid=0
		f.write("ATOM  {atomid:5d}  CA  {res} {chainid}{resid:4d}    {px:8.3f}{py:8.3f}{pz:8.3f}{occ:6.2f}{bfac:6.2f}     S_00  0\n".format(atomid=atomid, chainid=chr(int(chainid[i])+65),resid=resid, px=d[0], py=d[1], pz=d[2], occ=occ[i], bfac=bfac[i], res=amino_dict[residue[i]]))
		atomid+=1
		resid+=1

	if len(sidechain)>0:
		for s in sidechain:
			f.write("HETATM{atomid:5d}   S  {res} {chainid}{resid:4d}    {px:8.3f}{py:8.3f}{pz:8.3f}{occ:6.2f}{bfac:6.2f}     S_00  0\n".format(atomid=atomid, chainid=chr(int(chainid[i])+65),resid=int(s[0]),   px=s[-3], py=s[-2], pz=s[-1], occ=s[2], bfac=0, res=amino_dict[residue[i]]))
			atomid+=1

	curchain=chainid[i]
	

	f.write("TER  {:6d}      ALA {}{:4d}\n""".format(atomid+1, 'A', resid))
	f.write("ENDMDL\n")
	f.close()

def drawmap(pts, val, shape, apix=1.,fname="mapcolor_tmp0.hdf"):
	mapclr=np.zeros(shape)-1

	for i,p in enumerate(pts):
	#	 if dd0[i]<0: continue
		mapclr[p[2],p[1],p[0]]=val[i]
	e=from_numpy(mapclr)
	e["apix_x"]=e["apix_y"]=e["apix_z"]=apix
	e.write_image(fname)


### interpolate points
def interp_points(pts, npt=50, pmin=0., pmax=1.):
	pos=np.append(0,np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1)))
	fun_ax=interp1d(pos, pts.T, fill_value='extrapolate')
	mx=np.max(pos)
	rg=np.arange(npt,dtype=float)/(npt-1)*(pmax-pmin)*mx + pmin*mx
	ax=fun_ax(rg).T
	return ax



def calc_shortest_path(path, pts, dst, d0, pval, gap=5, gdsz=1.5, ngray=5., grp=[]):
# if 1:
	shrtpath=np.zeros(len(pts))
	order=np.zeros(len(pts))-1
	ordercount=np.zeros_like(order)
	for i0 in range(len(path)-1):
		realgap=min(gap,len(path)-i0-1)
		i1=i0+realgap
		p0=path[i0]
		p1=path[i1]

		j0=np.argmin(d0[:,i0])
		j1=np.argmin(d0[:,i1])


	#	 print thr
		dd0=np.zeros(len(pts))-1
		dd1=np.zeros(len(pts))-1

		dd0[j0]=0
		dd1[j1]=0
		
		if len(grp)>0:
			ingroup=(grp>=(i0-2)) * (grp<=(i1+2))
		connect=0
		for l in range(len(pts)):
				nbs0=np.any(dst[dd0==l]<=gdsz, axis=0) * (dd0<0)
				nbs1=np.any(dst[dd1==l]<=gdsz, axis=0) * (dd1<0)
				
				if len(grp)>0:
					nbs0*=ingroup
					nbs1*=ingroup
				
				if np.sum(nbs0)<1 and np.sum(nbs1)<1: break

	#			 print pval[nbs0]
				dd0[nbs0]=l+1
				dd1[nbs1]=l+1


				if dd0[j1]>0 or dd1[j0]>0:
					connect=1
					break
	#	 print connect
		mxl=l
#		 print l
		did=np.logical_and(dd0>=0, dd1>=0)
		if np.sum(did)<1:
			print "Density break between {} and {}".format(i0, i1)
#			 drawmap(pts, dd0, "mapcolor_tmp0.mrc")
#			 drawmap(pts, dd1, "mapcolor_tmp1.mrc")
			continue
		pv=pval[did]
		minp=np.min(pv)-.01
		thr=np.arange(minp, np.max(pv), (np.max(pv)-minp)/ngray)[::-1]

		dd0=np.zeros(len(pts))-1
		dd1=np.zeros(len(pts))-1

		dd0[j0]=0
		dd1[j1]=0
#		 mxl=len(pts)
		connect=0
		xl=0
		for tr in thr:
			if connect:
				break
			for l in range(int((mxl+2)*3)):
				nbs0=np.any(dst[dd0==l]<=gdsz, axis=0) * (dd0<0) * (pval>=tr)
				nbs1=np.any(dst[dd1==l]<=gdsz, axis=0) * (dd1<0) * (pval>=tr)
				if len(grp)>0:
					nbs0*=ingroup
					nbs1*=ingroup
#				 print tr,np.sum(nbs0), np.sum(nbs1), dd0[j1], dd1[j0]

				if dd0[j1]>=0: nbs0*=False
				if dd1[j0]>=0: nbs1*=False

				dd0[nbs0]=l+1
				dd1[nbs1]=l+1

				if l>xl and np.sum(nbs0)<1 and np.sum(nbs1)<1: break

				if dd0[j1]>=0 and dd1[j0]>=0:
					connect=1
					break
			xl=l+1

		# ngb=np.sort(dst1, axis=1)[:,1]
		# si= np.array([r<=np.min(ri1[dst1[j]==ngb[j]]) for j,r in enumerate(ri1)]).astype(int)

		dd=dd0+dd1
		did=np.logical_or(dd0<0, dd1<0)
		dd[did]=-1
		
		if np.sum(dd>0)==0:
			print "No connected density at {}".format(i0)
			continue

		dd-=np.min(dd[dd>0])
		dd=1-dd/(np.max(dd)+1)
		dd=np.max(dd)-dd+1
		dd[did]=-1
		
		a=(dd<=np.maximum(dd0[j1], dd1[j0])) * (dd>0)
		
		mult=(float(gap)/realgap)
		shrtpath[a]+=mult/gap
		
		dd0[dd0>dd0[j1]]=-1
		a=a*(dd0>0)
		dd0=dd0/dd0[j1]*realgap+i0
		order[a]+=dd0[a]
		ordercount[a]+=1
#		 print i0, i1,np.sum(a),np.sum(shrtpath>0), l, mxl, connect


		if connect==0:
			print "Breaking point near residue {}".format(i0)
			#break

	#print np.sum(shrtpath>0)
	
	ordercount[ordercount==0]=1
	order/=ordercount
	
	return shrtpath, order


def short_path_iter(path, pts, dst, d0, pval, gap=5, gdsz=1.5, ngray=5., grp=[], shrtpath=[]):
	if len(shrtpath)==0:
		ss,oo=calc_shortest_path(path, pts, dst, d0, pval, gap, gdsz, ngray, grp)
	else:
		bb=shrtpath>0
		bkbone=pts[bb]
		dbb=dst[bb][:, bb]
		db0=d0[bb]
		cnt=np.argmin(db0, axis=0)
		newpath=bkbone[cnt]
		db0=scidist.cdist(bkbone, newpath)
#		 pvnew=pval[bb]
		pvnew=shrtpath[bb]
		if len(grp)>0: ngrp=grp[bb]
		else: ngrp=[]
		sp1,o=calc_shortest_path(newpath, bkbone, dbb, db0, pvnew,  gap, gdsz, ngray, ngrp)
		ss=np.zeros_like(shrtpath)
		ss[shrtpath>0]=sp1
		oo=np.zeros_like(shrtpath)
		oo[shrtpath>0]=o
	return ss,oo


def diffuse_val(gp0, dst0, sdst=2., gdst=3.):

	gp1=np.zeros_like(gp0)-1
	for i in range(len(dst0)):
		rid=np.logical_and(dst0[i]<2,abs(gp0-gp0[i])<3)
		g=gp0[rid]
		wt=np.exp(-.5*dst0[i][rid]**2)
		wt/=np.sum(wt)
		gp1[i]=np.sum(wt*g)
	return gp1



def build_path(pts, gi, oripath, pval=[],mingap=3., apix=1.):
	
	pathnew=oripath.copy()
	ptsapix=pts*apix
	for i in range(3,len(oripath)-3):
		ii=(abs(gi-i-.5)<.5)
		ii= np.where(ii)[0]
		ig=ii
		if len(ii)==0:

			print "Skipping residue {}".format(i)
			continue
	#	 ig=gs[ii]
		if i>0:
			nrm=np.linalg.norm(ptsapix[ig]-pathnew[i-1], axis=1)

			ii=ii[nrm>=min(mingap, np.max(nrm))]
			ig=ii
		if len(pval)>0:
			mi=np.argmax(pval[ig])
			pathnew[i]= ptsapix[ig[mi]]#*apix
		else:
			pathnew[i]= np.mean(ptsapix[ig], axis=0)#*apix

	
	return pathnew


def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--folder", type=str,help="folder to store the output", default="pathrefine_00")
	parser.add_argument("--mapin", type=str,help="input map file", default=None)
	parser.add_argument("--pathin", type=str,help="input path file", default=None)
	parser.add_argument("--sequence", type=str,help="sequence file input", default=None)
	parser.add_argument("--mapthr", type=float,help="threshold for the input map", default=-1)
	parser.add_argument("--niter", type=int,help="number of iterations", default=5)
	parser.add_argument("--nres", type=int,help="number of residue", default=-1)
	parser.add_argument("--lowres", action="store_true", default=False, help="low resolution mode. more robust but may result in poor performance.")
	parser.add_argument("--skiprefine", action="store_true", default=False, help="skip path refinement and only perform sequence assignment. Require existing refinement folder")

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	try: os.mkdir(options.folder)
	except: pass
	
	
	e=EMData(options.mapin)
	img=e.numpy().copy()
	apix=e["apix_x"]
	shape=img.shape
	print "Shape of the density map: {}. Apix is {:.2f}".format(shape, apix)
	if abs(apix-1.)>.2:
		print "Apix maybe too big or too small. Program may not work.. ~1 A/pixel is recommended.."
	
	if options.mapthr<0:
		thr=np.mean(img)+np.std(img)
		print "Set threshold to {:.2f}.".format(thr)
	else:
		thr=options.mapthr
	
	pts=np.where(img>thr)
	pts=np.array(pts).T
	pts=pts[:,::-1]
	pval=img[pts[:,2], pts[:,1], pts[:,0]]
	pts_input=pts.copy()
	pval_input=pval.copy()
	print "{} voxels above the threshold in the density map. Average intensity value is {:.2f}".format(len(pts), np.mean(pval))
	
	
	path=pdb2numpy(options.pathin)
	print "Read {} point in the input path".format(len(path))
	
	if options.nres<0:
		options.nres=len(path)
	
	if options.lowres:
		print "Interpolate initial path.."
		path=interp_points(path,options.nres)
	
		
	
	pathint=(path/apix).astype(int)
	ival=np.mean(img[pathint[:,2], pathint[:,1], pathint[:,0]])
	print "Average intensity value on the path is {:.2f}".format(ival)
	if ival<thr:
		print "Intensity value on path lower than the map threshold. Program will not work. Make sure the path is docked in the density map correctly. Exit."
		exit()
		
	numpy2pdb(path, "{}/path_init.pdb".format(options.folder))
	
	
	dst=scidist.cdist(pts, pts)
	dst[np.eye(len(dst),dtype=bool)]=np.inf
	d0=scidist.cdist(pts, pathint)
	grp=[]
	print "Starting refinement..."
	niter=options.niter
	for itr in range(niter):
		if options.skiprefine:
			break
		
		if itr>0:
			path=pp.copy()
			pathint=(path/apix).astype(int)
			pts=spts.copy()
			pval=img[pts[:,2], pts[:,1], pts[:,0]]

			dst=scidist.cdist(pts, pts)
			dst[np.eye(len(dst),dtype=bool)]=np.inf
			d0=scidist.cdist(pts, pathint)
		
		if options.lowres:
			grp=group_pts(pts, dst, d0)
		
		shrtpath,sp=short_path_iter(pathint, pts, dst, d0, pval, 3, 1., 5, grp)
		shrtpath,sp=short_path_iter(pathint, pts, dst, d0, pval, 5, 1., 5, grp,shrtpath)
		drawmap(pts, shrtpath,shape, apix, "{}/map_bkbone_{:02d}.hdf".format(options.folder, itr))

		for i in range(5):  sp=diffuse_val(sp, dst,5,3)
			
			
		drawmap(pts, sp,shape, apix,"{}/map_seq_{:02d}.hdf".format(options.folder, itr))
		ss=sp[sp>0]
		spts=pts[sp>0]
		gi=np.zeros_like(ss)
		gi[np.argsort(ss)]=np.arange(len(ss), dtype=float)/len(ss)*len(path)
		pp=build_path(spts, gi, path, mingap=1.0, apix=apix)
		#     pp=build_path(spts, gi, pval[sp>0], mingap=1.0)
		pp=interp_points(pp, len(pp))
		numpy2pdb(pp, "{}/path_refine_{:02d}.pdb".format(options.folder, itr))
		print "iteration {:d}, {:d} voxels on the backbone, Average C-alpha distance {:.2f}".format(itr, int(np.sum(shrtpath>0)), np.mean(np.linalg.norm(np.diff(pp, axis=0), axis=1)))
	
	
		
	pts=pts_input
	pval=pval_input

	e=EMData("{}/map_bkbone_{:02d}.hdf".format(options.folder, niter-1))
	img=e.numpy().copy()
	thr=.1
	bkbone=np.where(img>thr)
	bkbone=np.array(bkbone).T
	bkbone=bkbone[:,::-1]


	e=EMData("{}/map_seq_{:02d}.hdf".format(options.folder, niter-1))
	img=e.numpy().copy()
	sval=img[bkbone[:,2], bkbone[:,1], bkbone[:,0]]

	dst=scidist.cdist(pts, pts)
	d0=scidist.cdist(pts, bkbone)
	dst[np.eye(len(dst),dtype=bool)]=np.inf
	
		
	brh=np.zeros(len(pts))-1
	brh[np.argmin(d0, axis=0)]=0

	grp=np.zeros(len(pts))-1
	grp[np.argmin(d0, axis=0)]=sval


	for i in range(100):
		gi=np.logical_and(brh<0, np.any(dst[brh==i]<=1., axis=0))
		#     print np.sum(grp<0),np.sum(gi)
		if np.sum(gi)==0: break
		brh[gi]=i+1
		grp[gi]=grp[brh==i][np.argmin(dst[brh==i][:,gi], axis=0)]
		#     grp[gi]=
	
	
	drawmap(pts, brh, shape,apix,"{}/map_branch.hdf".format(options.folder))
	drawmap(pts, grp, shape,apix,"{}/map_branch_seq.hdf".format(options.folder))
	pp=build_path(bkbone, sval+1, path, mingap=1.0, apix=apix)
		
	
	scsz=np.zeros(len(path))
	sidechains=[]
	for i in range(len(path)):
		gid=np.logical_and(abs(grp-(i-.5))<.5, brh>1)
		dg=dst[gid][:,gid]
		ng=np.sum(gid)
		dgp=np.zeros(ng)
		scsz[i]=np.sum(pval[gid])
		
		if i<2 or i>=len(path)-2: continue
		
		for k in range(10):
			did=np.where(dgp==0)[0]
			if len(did)==0: break
			dgp[did[0]]=k+1
			for j in range(ng):
				di=np.logical_and(dgp<=0, np.any((dg[dgp==k+1])<2., axis=0))
				if np.sum(di)==0: break
				dgp[di]=k+1
			p=pts[gid][dgp==k+1]
			if len(p)>3:
				sc=[i, k, len(p)]
				sc.extend(np.mean(p, axis=0)*apix)
				sidechains.append(sc)
		
	sc=np.array(sidechains)
	sdst=scidist.cdist(sc[:,-3:],sc[:,-3:])
	sdst[np.eye(len(sc), dtype=bool)]=np.inf


	sa=np.array(np.where(sdst<2.)).T
	for s in sa:
		if s[1]>s[0]:
			nm= np.mean([sc[s[0]][-3:],sc[s[1]][-3:]], axis=0)
			sidechains[s[0]][3:]=nm
			sidechains[s[1]]=[]
	sidechains=[s for s in sidechains if len(s)==6]
	print "Found {:d} sidechain residues. Average C-alpha distance {:.2f}".format(len(sidechains), np.mean(np.linalg.norm(np.diff(pp, axis=0), axis=1)))
	numpy2pdb_sc(pp, "{}/path_refine_final.pdb".format(options.folder), bfac=scsz, sidechain=sidechains)
	

	if options.sequence==None:
		print "Done."
		E2end(logid)
		return
	
	
	
	seq=np.array(read_fasta(options.sequence))
	seqwt=molwt[seq]

	s0=scsz.copy()#[::-1]
	s1=seqwt.copy()

	# s1[s1<np.mean(s1)]=0
	s1-=40
	s1[s1<0]=0
	s1/=np.max(s1)
	s1=s1**2
	s0/=np.max(s0)

	#print np.max(s0), np.max(s1), np.std(s0), np.std(s1)

	s0=np.round(s0*19).astype(int)
	s1=np.round(s1*19).astype(int)
	
	
	lkup=np.array([chr(i+65) for i in range(20)])

	pw=[[],[]]
	for inv,is0 in enumerate([s0, s0[::-1]]):

		ss0=lkup[is0].tostring()
		ss1=lkup[s1].tostring()

		cmpdic={}
		for i in range(20):
			for j in range(20):
				cmpdic[(lkup[i], lkup[j])]=((19.-abs(i-j))/20.)




		pw[inv]=pairwise2.align.globaldx(ss0, ss1, cmpdic)
		#print pw[inv][0]


	inv=np.argmax([pw[i][0][2] for i in [0,1]])
	ppw=pw[inv]
	
	if inv>0: print "Chain is reversed.",
	print "Match score is {}".format(ppw[0][2]/float(len(path)))

	l=[[],[]]
	for k in [0,1]:
		
		for a in ppw[0][k]:
			if a=='-':
				l[k].append(0)
			else:
				l[k].append(ord(a)-65)
	
	
	

	ar=np.arange(len(s1))
	ind0=np.array([i!='-' for i in ppw[0][0]])
	ii0=np.zeros(len(ind0))-1
	ii0[ind0>0]=ar

	ind1=np.array([i!='-' for i in ppw[0][1]])
	ii1=np.zeros(len(ind1))-1
	# if inv==1: ar=ar[::-1]
	ii1[ind1>0]=ar

	seq_asign=ii1[ind0>0]
	if inv>0: seq_asign=seq_asign[::-1]
	#     seq_asign=len(path)-seq_asign-1
	#     seq_asign[seq_asign>=len(path)]=-1

	ss=s1.copy()
	# if inv==1: ss=ss[::-1]
	ss=np.append(ss,0)
	# ss=ss/np.max(ss)*20
	ss=np.array([ss[int(i)] for i in seq_asign])
	
		
	sv=sval+1
	print np.min(sv), np.max(sv)
	sq=np.arange(len(path))
	sid=np.where(seq_asign>=0)[0]
	# print [(sq[i], seq_asign[i]) for i in sid]
	# print sq[sid]
	# print seq_asign[sid]

	fun_x=interp1d(sq[sid], seq_asign[sid], fill_value='extrapolate')

	fx= fun_x(sv)
	# print np.min(fx), np.max(fx)
	# plt.plot(fx-(274-sv))
	# plt.plot(274-sv)
	if inv==0:
		pp0=pp.copy()
	else:
		pp0=pp[::-1]

	pp1=build_path(bkbone, fx+1, pp0, mingap=1.0, apix=apix)
	sc=np.array(sidechains)
	sc[:,0]= np.round(fun_x(sc[:,0]))
	
	
	scsz_srt=np.zeros_like(scsz)
	for i in range(len(path)):
		#print i, scsz[i], seq_asign[i], seqwt[int(seq_asign[i])], fun_x(i)
		try: scsz_srt[int(fun_x(i))]=scsz[i]
		except: pass
		
		
	numpy2pdb_sc(pp1, "{}/path_refine_seq.pdb".format(options.folder), bfac=scsz_srt, sidechain=sc, residue=seq)
	
	print "Done."
	
	E2end(logid)
	
def run(cmd):
	print cmd
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	