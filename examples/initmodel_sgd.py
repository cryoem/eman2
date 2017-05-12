#!/usr/bin/env python
# Muyuan Chen 2017-03
import numpy as np
from e2simmx import cmponetomany
from EMAN2 import *
import time
import threading
import Queue


def make3d(aptcls, sym="c1"):
	boxsize=aptcls[0]["nx"]
	pad=good_size(boxsize*3/2)
	recon=Reconstructors.get("fourier", {"sym":sym,"size":[pad,pad,pad]})

	# insert slices into initial volume
	recon.setup()
	for p in aptcls:
		p2=p.get_clip(Region(-(pad-boxsize)/2,-(pad-boxsize)/2,pad,pad))
		p3=recon.preprocess_slice(p2,p["xform.projection"])
		recon.insert_slice(p3,p["xform.projection"],p.get_attr_default("ptcl_repr",1.0))

	threed=recon.finish(True)
	threed.process_inplace("xform.applysym", {"sym":sym})
	threed.process_inplace("xform.centerofmass")
	threed=threed.get_clip(Region((pad-boxsize)/2,(pad-boxsize)/2,(pad-boxsize)/2,boxsize,boxsize,boxsize))
	threed.process_inplace("normalize")
	threed.process_inplace("mask.gaussian",{"inner_radius":boxsize/3.0,"outer_radius":boxsize/12.0})

	
	
	threed["apix_x"]=threed["apix_y"]=threed["apix_z"]=aptcls[0]["apix_x"]
	return threed

def makeprojs(threed, sym):
	origen=parsesym(sym)
	oris=origen.gen_orientations('eman',{"delta":9., "perturb":True})
	projs=[threed.project("standard",ort) for ort in oris]
	return projs


def do_ali(ptcls, projs):
	boxsize=ptcls[0]["nx"]
	bslst=[]
	quals=[]
	bss=.0
	pjs=[(p, None) for p in projs]
	for i in range(len(ptcls)):
		#sim=cmponetomany(pjs,ptcls[i],align=("rotate_translate_flip",{"maxshift":boxsize/5}),alicmp=("ccc",{}),ralign=("refine",{}),cmp=("frc",{"minres":80,"maxres":20}))
		sim=cmponetomany(pjs,ptcls[i],align=("rotate_translate_tree",{"maxshift":boxsize/5, "maxres":15}),alicmp=("ccc",{}),ralign=("refine",{}),cmp=("frc",{"minres":80,"maxres":10}))
		bs=min(sim)
		
		bss+=bs[0]
		bslst.append((bs[0],i))
		n=sim.index(bs)
#		 print n
		ptcls[i]["match_n"]=n
		ptcls[i]["match_qual"]=-bs[0]
		ptcls[i]["xform.projection"]=projs[n]["xform.projection"]# best orientation set in the original particle

	bslst.sort()# sorted list of all particle qualities
	bslst.reverse()
	aptcls=[]
	#for i in range(len(ptcls)*3/4):# We used to include 3/4 of the particles
	for i in range(len(ptcls)):
		n=ptcls[bslst[i][1]]["match_n"]
		quals.append(ptcls[bslst[i][1]]["match_qual"])
		aptcls.append(ptcls[bslst[i][1]].align("rotate_translate_tree",projs[n],{"maxshift":boxsize/5, "maxres":10},"frc",{}))
		aptcls[-1].process_inplace("normalize.toimage",{"to":projs[n]})
		#aptcls[-1].process_inplace("normalize")
		aptcls[-1].add(-aptcls[-1]["mean"])
	return aptcls

def random_sample(n,k):
	a=np.arange(n)
	np.random.shuffle(a)
	return a[:int(k)]

def make_model(jsd,myid, options):
	ptclname=options.ptcls
	num=EMUtil.get_image_count(ptclname)
	sym=options.sym
	nsample=options.batchsize
	learnrate=options.learnrate
	lrmult=options.lrdecay
	path=options.path
	
	aptcls=[EMData(ptclname, i) for i in random_sample(num, nsample)]
	origen=parsesym(sym)
	oris=origen.gen_orientations('rand',{"n":nsample+1, "phitoo":1})
	for i in range(nsample):
		aptcls[i]["xform.projection"]=oris[i]
		
	map0=make3d(aptcls,sym)
	map0.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1})
	map0.process_inplace("filter.lowpass.randomphase",{"cutoff_abs":.05})
	#map0.write_image("{}/model_00.hdf".format(path))
	boxsize=map0["nx"]
	niter=options.niter

	ncopy=(niter+1)*nsample/num+1
	samples=np.tile(np.arange(num),(ncopy,1))
	#print "ncopy:", ncopy
	for i in range(len(samples)): np.random.shuffle(samples[i])
	samples=samples.flatten()
	kk=0
	scr=[]
	for it in range(1,niter):
		
		pjs=None
		if options.writetmp: map0.write_image("{}/model_tmp_{:02d}.hdf".format(path,myid), it-1)
		pjs=makeprojs(map0, sym)
			
		if it==niter-1:
			if num<32:
				ptcls=[EMData(ptclname, i) for i in range(num)]
			else:
				smp=random_sample(num, 32)
				ptcls=[EMData(ptclname, i) for i in smp]
		else:
			ptcls=[EMData(ptclname, i) for i in samples[kk:kk+nsample]]
			
		for p in ptcls:
			if options.addnoise>0:
				p.process_inplace("math.addsignoise",{"noise":options.addnoise})
			p.process_inplace("filter.lowpass.gauss",{"cutoff_freq":.1})
			p.process_inplace("normalize.edgemean")
			
			
		kk+=nsample
	#	 print len(ptcls)
		aptcls=do_ali(ptcls, pjs)
		
		
			
		for i, p in enumerate(aptcls):
			scr.append(p.cmp("frc", pjs[p["match_n"]]))
		map1=make3d(aptcls,sym)	
			
		map0ft=map0.do_fft()
		map1ft=map1.do_fft()
		map1ft.process_inplace("mask.wedgefill",{"fillsource":map0ft, "thresh_sigma":.9})
		dmap=map1ft-map0ft
		
		
		#ddm=dmap*dmap
		#print "{:d}\t{:.3f}\t{:.3f}".format(it, ddm["mean_nonzero"], np.mean(scr))
		mapnewft=map0ft+learnrate*dmap
		mapnew=mapnewft.do_ift()
		#if it<niter/2:
			
			#mapnew.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.15})
		#else:
			#mapnew.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.3})
		learnrate*=lrmult
		
		
		mapnew.process_inplace("filter.lowpass.gauss",{"cutoff_freq":.1})
		#mapnew.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.3})
		mapnew.process_inplace("mask.auto3d",{"radius":boxsize/6,"threshold":map0["sigma_nonzero"]*.85,"nmaxseed":30,"nshells":boxsize/20,"nshellsgauss":boxsize/20})
		ddmap=mapnew-map0
		ddmap.mult(ddmap)
		if options.verbose>0:print myid, it, ddmap["mean"]
		
		
		map0=mapnew
		#if ddmap["mean"]<.005:
			#break
	
	compares=[]
	for i, p in enumerate(aptcls):
		compares.append(p)
		compares.append(pjs[p["match_n"]])
	scr=np.array(scr)
	jsd.put((map0, ddmap["mean"], np.mean(scr[-num:]), compares))
	pjs=None
	return 
	#launch_childprocess("e2proclst.py {}/model_*.hdf --create {}/models.lst".format(path, path));

def main():
	
	usage="""
	This program makes initial models using a (kind of) stochastic gradient descent approach.
	[prog] --ptcls <particle stack> 
	
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="Path to write initial model output. Default is initmodel_XX", default=None)
	parser.add_argument("--ptcls", type=str,help="Class average or particles input.", default=None)
	parser.add_argument("--sym", type=str,help="Symmetry", default='c1')
	parser.add_argument("--batchsize", type=int,help="Batch size of stochastic gradient desent. N particles are randomly selected to generate an initial model at each step.", default=5)
	parser.add_argument("--learnrate", type=float, help="Learning rate. i.e. how much the initial model changes toward the gradient direction in each iteration. Range from 0.0~1.0. Default is 0.3", default=.3)
	parser.add_argument("--lrdecay", type=float, help="Learning rate multiplier after each iteration.", default=1.)
	parser.add_argument("--addnoise", type=float, help="Add noise on particles at each iteration. Stablize convergence for some reason.", default=3.)
	parser.add_argument("--ntry", type=int,help="Number of tries.", default=10)
	parser.add_argument("--writetmp", action="store_true", default=False ,help="Write output for each iteration")
	parser.add_argument("--threads", type=int,help="threads", default=10)
	parser.add_argument("--niter", type=int,help="Number of iterations", default=20)
	parser.add_argument("--verbose", type=int,help="Verbose", default=0)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)


	tt0=time.time()

	path=options.path
	if path==None:
		for i in range(100):
			try:
				path="initmodel_{:02d}".format(i)
				os.mkdir(path)
				options.path=path
				print "Writing initial model in {}..".format(path)
				break
			except:
				continue
		else:
			print "Too many initmodel folders in the project, or something odd happened....Exit."
			exit()
		
	
	
	jsd=Queue.Queue(0)	
	NTHREADS=max(options.threads+1,2)
	thrds=[threading.Thread(target=make_model,args=(jsd, i, options)) for i in range(options.ntry)]	
	thrtolaunch=0
	models=[]
	scrs=[]
	grads=[]
	cmps=[]
	
	while thrtolaunch<len(thrds) or threading.active_count()>1:
		# If we haven't launched all threads yet, then we wait for an empty slot, and launch another
		# note that it's ok that we wait here forever, since there can't be new results if an existing
		# thread hasn't finished.
		if thrtolaunch<len(thrds) :
			while (threading.active_count()==NTHREADS ) : time.sleep(.1)
			thrds[thrtolaunch].start()
			thrtolaunch+=1
		else: time.sleep(1)
	
		while not jsd.empty():
			m, grad, scr,c=jsd.get()
			models.append(m)
			scrs.append(scr)
			grads.append(grad)
			cmps.append(c)

	for t in thrds:
		t.join()

	
	sid=np.argsort(scrs)
	for i,si in enumerate(sid):
		print i,si,grads[si],scrs[si]
		models[si].write_image("{}/model_{:02d}.hdf".format(path, i))
		for c in cmps[si]:
			c.write_image("{}/model_{:02d}_aptcl.hdf".format(path, i), -1)
		
		
	
	#m, grad, scr=make_model(options)
	print "total time:", time.time()-tt0


	E2end(logid)

	
def run(cmd):
	print cmd
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	