#!/usr/bin/env python
# Muyuan Chen 2024-06
from EMAN2 import *
import numpy as np
from sklearn.decomposition import PCA
import scipy.spatial.distance as scipydist
from sklearn import cluster,mixture
import umap

def main():
	
	usage=""" New single particle refinement routine. Still under construction. For simple usage,
	e2spa_refine2d.py --ptcl <particle list file> 
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ptcls", type=str,help="particle input", default=None)
	parser.add_argument("--setsf", type=str,help="structure factor for sharpening", default=None)
	parser.add_argument("--path", type=str,help="path", default=None)
	parser.add_argument("--init", type=str,help="initial averages", default=None)
	parser.add_argument("--mask", type=str,help="mask", default=None)
	parser.add_argument("--ncls", type=int,help="number of classes", default=64)
	parser.add_argument("--startiter", type=int,help="starting iteration", default=0)
	parser.add_argument("--aliref", type=int,help="number of alignment references", default=10)
	parser.add_argument("--maxres", type=float,help="max resolution", default=5)
	parser.add_argument("--niter", type=float,help="number of iterations", default=6)
	parser.add_argument("--parallel", type=str,help="Thread/mpi parallelism to use. Default is thread:12", default="thread:32")
	parser.add_argument("--tomo", action="store_true", default=False ,help="tomo averager.")

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	#### structure factor for sharpening
	s=np.loadtxt(options.setsf)
	x=np.arange(len(s))
	y=np.tanh(x*.1)
	s[:,1]-=np.min(s[:,1])
	s[:,1]*=y

	s[:,1]/=np.max(s[:,1])
	sfname="strucfac_hp.txt"
	np.savetxt(sfname, s)

	sf=XYData()
	sf.read_file(sfname)

	####
	path=options.path
	ptcls=options.ptcls
	
	run(f"e2proc2d.py {options.init} {path}/classes_00.hdf --inplace")
	cutoff=[ 0.1, 0.2, 0.3, 0.4]+[0.4]*options.niter
	cutoff=np.array(cutoff)
	
	e=EMData(f"{path}/classes_00.hdf",0,True)
	apix=e["apix_x"]
	maxabs=apix/options.maxres
	cutoff=np.minimum(cutoff, maxabs)
	print(cutoff)
	maxalires=max(options.maxres, 10)
	if options.mask:
		mask=EMData(options.mask)
	
	for itr in range(options.startiter, options.niter):
		
		fname=f"{path}/classes_{itr:02d}.hdf"
		n=EMUtil.get_image_count(fname)
		imgs=[]

		for i in range(n):
			data=EMData(fname,i)
			data.process_inplace("normalize.circlemean",{"radius":-10})
			# data.process_inplace("normalize.edgemean")
			data.process_inplace("mask.soft",{"outer_radius":-8})
			# data.process_inplace("filter.matchto",{"to":ref,"interpolate":1,"keephires":1})
			data.process_inplace("filter.setstrucfac",{"apix":data["apix_x"],"strucfac":sf})
			data.process_inplace("filter.lowpass.gauss",{"cutoff_abs":cutoff[itr]})
			
			if options.mask:
				data.process_inplace("normalize.edgemean")
				data.mult(mask)
			# data.process_inplace("mask.soft",{"outer_radius":-16})
			# data.process_inplace("normalize.circlemean",{"radius":-10})
			data.write_image(f"{path}/classes_{itr:02d}_sf.hdf", i)
			imgs.append(data)
			
		####
		cvs=[]
		for e in imgs:
			ny=e["ny"]
			ef=e.process("xform.centeracf")
			ef.process_inplace("mask.soft",{"inner_radius":ny//3})
			ef = ef.do_fft()
			cv = ef.calc_radial_dist(ny//2, 0, 1.0,True)
			cv=np.array(cv)/(ny**3)
			cvs.append(cv)
			
		cvs=np.array(cvs)
		sig=np.mean(cvs[:,8:], axis=1)
		
		####
		lst=[]
		for i in np.argsort(sig):
			l={"src":f"{path}/classes_{itr:02d}_sf.hdf", "idx":i}
			lst.append(l)
			
		save_lst_params(lst, f"{path}/classes_{itr:02d}_sort.lst")
		
		imgs=[imgs[i] for i in np.argsort(sig)[:20]]
		imgs_ali=[e.align("rotate_translate_tree", imgs[0], {"maxres":maxalires}) for e in imgs]
		imgsnp=np.array([e.numpy().copy() for e in imgs_ali])
		
		####
		pca=PCA(1)
		imgsflat=imgsnp.reshape((len(imgsnp), -1))
		p2=pca.fit_transform(imgsflat).flatten()
		i2=np.argsort(p2)
		
		tmp=[]
		for i in [0,-1]:
			imgs_ali=[e.align("rotate_translate_tree", imgs[i], {"maxres":maxalires}) for e in imgs]
			imgsnp=np.array([e.numpy().copy() for e in imgs_ali])
			tmp.append(imgsnp.reshape((len(imgsnp),-1)))
			
		imgsflat=np.hstack(tmp)
		pca=PCA(1)
		p2=pca.fit_transform(imgsflat).flatten()
		i2=np.argsort(p2)
		
		####
		imgs2=[]
		for i,k in enumerate(i2):
			if i==0:
				a=imgs[k]
			else:
				a=imgs[k].align("rotate_translate_tree", a, {"maxres":maxalires})
				
			a.write_image(f"{path}/classes_{itr:02d}_top20.hdf", i)
			imgs2.append(a)
			
		if itr+1==options.niter:
			break
			
		####
		nbasis=min(len(imgs2), 12)
		pca=PCA(nbasis)
		# imgs3=[m.process("filter.lowpass.gauss",{"cutoff_abs":.25}) for m in imgs2]
		imgs3=imgs2
		imgsnp=np.array([m.numpy().copy() for m in imgs3])
		imgsflat=imgsnp.reshape((len(imgsnp), -1))
		p2=pca.fit_transform(imgsflat)
		tosave=pca.components_
		tosave=np.append([np.mean(imgsflat, axis=0)], tosave, axis=0)

		for i,p in enumerate(tosave):
			m=p.reshape(imgsnp[0].shape)
			e=from_numpy(m).copy()
			e.set_attr_dict(imgs[0].get_attr_dict())
			e.process_inplace("normalize.edgemean")
			e.write_image(f"{path}/basis_{itr:02d}.hdf", i)
			
		####
		p=p2#[:,:3]
		dst=scipydist.cdist(p, p)
		keep=[np.argmin(i2)]
		nref=options.aliref
		for i in range(nref-1):
			d=dst[keep]
			d=np.min(d, axis=0)
			keep.append(np.argmax(d))
			# print(d)
			
		for i,k in enumerate(keep):
			e=EMData(f"{path}/classes_{itr:02d}_top20.hdf", k)
			e.write_image(f"{path}/aliref_{itr:02d}.hdf", i)

		####
		e2simmxcmd = f"e2simmx.py {path}/aliref_{itr:02d}.hdf {ptcls} {path}/simmx_{itr:02d}.hdf --saveali --cmp=frc --align=rotate_translate_tree:maxres={options.maxres} --aligncmp=frc:maxres={options.maxres} --parallel {options.parallel}"
		run(e2simmxcmd)
		
		run(f"e2basis.py projectrot {path}/basis_{itr:02d}.hdf {ptcls} {path}/simmx_{itr:02d}.hdf {path}/input_{itr+1:02d}_proj.hdf --oneout --mean1 --normproj --parallel {options.parallel}")
		
		####
		e=EMData(f"{path}/input_{itr+1:02d}_proj.hdf")
		e=e.numpy().copy()
		eg=e[:,4:]
		exf=e[:,:4]
		
		pca=umap.UMAP(n_components=2)
		p2=pca.fit_transform(eg)
		p2-=np.mean(p2, axis=0)
		p2/=np.std(p2, axis=0)
		np.savetxt(f"{path}/umap_{itr+1:02d}.txt", p2)
		
		ncls=options.ncls
		clust=cluster.KMeans(ncls)
		lbs=clust.fit_predict(p2)
		lbunq=np.unique(lbs)
		clsmx=f"{path}/classmx_{itr+1:02d}.hdf"
		
		####
		tosv=np.zeros((len(eg), 6), dtype=np.float32)
		tosv[:,0]=lbs
		tosv[:,1]=1
		tosv[:,2:]=exf
		for i in range(6):
			a=from_numpy(tosv[:,[i]])
			a.write_image(clsmx, i)
			
		avgr="mean"
		if options.tomo:
			avgr="mean.tomo"
			
		run(f"e2classaverage.py --input={ptcls} --classmx={clsmx} --output={path}/classes_{itr+1:02d}.hdf --center xform.centerofmass --iter=5 --align=rotate_translate_tree:maxres={options.maxres} --averager={avgr}  --keep=.8 --cmp=frc --aligncmp=frc:maxres={options.maxres} --parallel {options.parallel} --resultmx {path}/resultmx_{itr+1:02d}.hdf")
		
		
	
	E2end(logid)
	
	
if __name__ == '__main__':
	main()
	
