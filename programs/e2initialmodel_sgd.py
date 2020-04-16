#!/usr/bin/env python
# Muyuan Chen 2017-03
from past.utils import old_div
from future import standard_library
standard_library.install_aliases()
from builtins import range
import numpy as np
from EMAN2 import *
from EMAN2_utils import cmponetomany
import time
import threading
import queue


def make3d(aptcls, sym="c1"):
	boxsize=aptcls[0]["nx"]
	pad=good_size(old_div(boxsize*3,2))
	recon=Reconstructors.get("fourier", {"sym":sym,"size":[pad,pad,pad]})

	# insert slices into initial volume
	recon.setup()
	for p in aptcls:
		p2=p.get_clip(Region(old_div(-(pad-boxsize),2),old_div(-(pad-boxsize),2),pad,pad))
		p3=recon.preprocess_slice(p2,p["xform.projection"])
		recon.insert_slice(p3,p["xform.projection"],p.get_attr_default("ptcl_repr",1.0))

	threed=recon.finish(True)
	threed.process_inplace("xform.applysym", {"sym":sym})
	threed.process_inplace("xform.centerofmass")
	threed=threed.get_clip(Region(old_div((pad-boxsize),2),old_div((pad-boxsize),2),old_div((pad-boxsize),2),boxsize,boxsize,boxsize))
	threed.process_inplace("normalize")
	threed.process_inplace("mask.gaussian",{"inner_radius":old_div(boxsize,3.0),"outer_radius":old_div(boxsize,12.0)})

	
	
	threed["apix_x"]=threed["apix_y"]=threed["apix_z"]=aptcls[0]["apix_x"]
	return threed

def makeprojs(threed, options):
	origen=parsesym(options.sym)
	oris=origen.gen_orientations('eman',{"delta":9., "perturb":True})
	if options.fullcov:
		np.random.shuffle(oris)
		oris=oris[:options.batchsize]
	projs=[threed.project("standard",ort) for ort in oris]
	return projs

def do_ali_fullcov(ptcls, projs):
	boxsize=ptcls[0]["nx"]
	quals=[]
	bss=.0
	pts=[(p, None) for p in ptcls]
	for i in range(len(projs)):
		#sim=cmponetomany(pjs,ptcls[i],align=("rotate_translate_flip",{"maxshift":boxsize/5}),alicmp=("ccc",{}),ralign=("refine",{}),cmp=("frc",{"minres":80,"maxres":20}))
		sim=cmponetomany(pts,projs[i],align=("rotate_translate_tree",{"maxshift":old_div(boxsize,5), "maxres":20}),alicmp=("ccc",{}),ralign=("refine",{}),cmp=("frc",{"minres":80,"maxres":20}))
		bs=min(sim)
		
		bss+=bs[0]
		n=sim.index(bs)
#		 print n
		projs[i]["match_n"]=n
		projs[i]["match_qual"]=-bs[0]

	aptcls=[]
	#for i in range(len(ptcls)*3/4):# We used to include 3/4 of the particles
	for i in range(len(projs)):
		n=projs[i]["match_n"]
		quals.append(projs[i]["match_qual"])
		aptcls.append(ptcls[n].align("rotate_translate_tree",projs[i],{"maxshift":old_div(boxsize,5), "maxres":20},"frc",{}))
		aptcls[-1].process_inplace("normalize.toimage",{"to":projs[i]})
		#aptcls[-1].process_inplace("normalize")
		aptcls[-1]["match_n"]=i
		aptcls[-1].add(-aptcls[-1]["mean"])
		aptcls[-1]["xform.projection"]=projs[i]["xform.projection"]# best orientation set in the original
	return aptcls

def do_ali(ptcls, projs):
	boxsize=ptcls[0]["nx"]
	bslst=[]
	quals=[]
	bss=.0
	pjs=[(p, None) for p in projs]
	for i in range(len(ptcls)):
		#sim=cmponetomany(pjs,ptcls[i],align=("rotate_translate_flip",{"maxshift":boxsize/5}),alicmp=("ccc",{}),ralign=("refine",{}),cmp=("frc",{"minres":80,"maxres":20}))
		sim=cmponetomany(pjs,ptcls[i],align=("rotate_translate_tree",{"maxshift":old_div(boxsize,5), "maxres":20}),alicmp=("ccc",{}),ralign=("refine",{}),cmp=("frc",{"minres":160,"maxres":20}))
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
		aptcls.append(ptcls[bslst[i][1]].align("rotate_translate_tree",projs[n],{"maxshift":old_div(boxsize,5), "maxres":20},"frc",{}))
		aptcls[-1].process_inplace("normalize.toimage",{"to":projs[n]})
		#aptcls[-1].process_inplace("normalize")
		aptcls[-1].add(-aptcls[-1]["mean"])
	return aptcls

def random_sample(n,k):
	a=np.arange(n)
	np.random.shuffle(a)
	return a[:int(k)].tolist()

def make_model(jsd,myid, options):
	ptclname=options.ptcls
	num=EMUtil.get_image_count(ptclname)
	sym=options.sym
	nsample=options.batchsize
	learnrate=options.learnrate
	lrmult=options.lrdecay
	path=options.path
	
	ninit=min(32, num)
	aptcls=[EMData(ptclname, i) for i in random_sample(num, ninit)]
	origen=parsesym(sym)
	oris=origen.gen_orientations('rand',{"n":ninit+1, "phitoo":1})
	for i in range(ninit):
		aptcls[i]["xform.projection"]=oris[i]
		if options.shrink!=1:
			aptcls[i].process_inplace("math.fft.resample",{"n":options.shrink})
		aptcls[i].process_inplace("normalize")
	map0=make3d(aptcls,sym)
	map0.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1})
	map0.process_inplace("filter.lowpass.randomphase",{"cutoff_abs":.05})
	map0.process_inplace("normalize")
	#map0.write_image("{}/model_00.hdf".format(path))
	boxsize=map0["nx"]
	niter=options.niter

	ncopy=old_div((niter+1)*nsample,num)+10
	samples=np.tile(np.arange(num),(ncopy,1))
	#print "ncopy:", ncopy
	for i in range(len(samples)): np.random.shuffle(samples[i])
	samples=samples.flatten().tolist()
	kk=0
	scr=[]
	for it in range(1,niter):
		if it==niter-1:
			#### bigger batch size in the last iteration so we can sort the initial model stably 
			options.batchsize=max(64, options.batchsize)
		pjs=None
		if options.writetmp: map0.write_image("{}/model_tmp_{:02d}.hdf".format(path,myid), it-1)
		pjs=makeprojs(map0, options)
		
		
		if options.fullcov or it==niter-1:
			if num<32:
				ptcls=[EMData(ptclname, i) for i in range(num)]
			else:
				smp=random_sample(num, 32)
				ptcls=[EMData(ptclname, i) for i in smp]
		else:
			ptcls=[EMData(ptclname, i) for i in samples[kk:kk+options.batchsize]]
			
		for p in ptcls:
			if options.shrink!=1:
				p.process_inplace("math.fft.resample",{"n":options.shrink})
			
			p.process_inplace("normalize")
			if options.addnoise>0:
				p.process_inplace("math.addsignoise",{"noise":options.addnoise})
			p.process_inplace("filter.lowpass.gauss",{"cutoff_freq":old_div(1.,options.targetres)})
			#p.process_inplace("normalize.edgemean")
			p.process_inplace("normalize")
			
			
		kk+=options.batchsize
	#	 print len(ptcls)
		if options.fullcov:
			aptcls=do_ali_fullcov(ptcls, pjs)
		else:
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
		
		
		mapnew.process_inplace("filter.lowpass.gauss",{"cutoff_freq":old_div(1.,options.targetres)})
		#mapnew.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.3})
		#mapnew.process_inplace("mask.auto3d",{"radius":boxsize/6,"threshold":map0["sigma_nonzero"]*.85,"nmaxseed":30,"nshells":boxsize/20,"nshellsgauss":boxsize/20})
		ddmap=mapnew-map0
		ddmap.mult(ddmap)
		if options.verbose>0:print(myid, it, ddmap["mean"])
		
		
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
	This program makes initial models using a (kind of) stochastic gradient descent approach. It is recommended that the box size of particles is around 100. 
	[prog] --ptcls <particle stack> 
	
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	#parser.add_header(name="initialmodelheader", help='Options below this label are specific to e2initialmodel', title="### e2initialmodel options ###", row=1, col=0, rowspan=1, colspan=3)
	parser.add_argument("--path", type=str,help="Path to write initial model output. Default is initmodel_XX", default=None)
	parser.add_argument("--ptcls", type=str,help="Class average or particles input.", default=None, browser='EMBrowserWidget(withmodal=True,multiselect=False)', guitype='filebox', row=0, col=0, rowspan=1, colspan=3)
	parser.add_argument("--sym", type=str, default='c1', help = "Specify symmetry - choices are: c<n>, d<n>, h<n>, tet, oct, icos", guitype='symbox', row=1, col=0, rowspan=1, colspan=1)
	parser.add_argument("--batchsize", type=int,help="Batch size of stochastic gradient desent. N particles are randomly selected to generate an initial model at each step.", default=10, guitype='intbox', row=1, col=1, rowspan=1, colspan=1)
	parser.add_argument("--niter", type=int,help="Number of iterations", default=20, guitype='intbox', row=2, col=0, rowspan=1, colspan=1)
	parser.add_argument("--ntry", type=int,help="The number of different initial models to generate in search of a good one", default=10, guitype='intbox', row=2, col=1, rowspan=1, colspan=1)
	parser.add_argument("--learnrate", type=float, help="Learning rate. i.e. how much the initial model changes toward the gradient direction in each iteration. Range from 0.0~1.0. Default is 0.3", default=.3, guitype='floatbox', row=3, col=0, rowspan=1, colspan=1)
	parser.add_argument("--lrdecay", type=float, help="Learning rate multiplier after each iteration.", default=1., guitype='floatbox', row=3, col=1, rowspan=1, colspan=1)
	parser.add_argument("--addnoise", type=float, help="Add noise on particles at each iteration. Stablize convergence for some reason.", default=3., guitype='floatbox', row=4, col=0, rowspan=1, colspan=1)
	parser.add_argument("--shrink", type=int,help="shrinking factor", default=1, guitype='intbox', row=4, col=1, rowspan=1, colspan=1)
	parser.add_argument("--writetmp", action="store_true", default=False ,help="Write output for each iteration",guitype='boolbox', row=5, col=0, rowspan=1, colspan=1)
	parser.add_argument("--fullcov", action="store_true", default=False ,help="Assume the input particles covers most of the orientation of the model. This gives better performance when the model is relatively feature-less, but is more likely to fail when there are incorrect particles in the input.",guitype='boolbox', row=5, col=1, rowspan=1, colspan=1)
	parser.add_argument("--threads", type=int,help="threads", default=10, guitype='intbox', row=6, col=0, rowspan=1, colspan=1)
	parser.add_argument("--verbose", "-v", type=int,help="Verbose", default=0, guitype='intbox', row=6, col=1, rowspan=1, colspan=1)
	parser.add_argument("--targetres", type=float, help="Target resolution", default=20., guitype='floatbox', row=7, col=0, rowspan=1, colspan=1)
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
				print("Writing initial model in {}..".format(path))
				break
			except:
				continue
		else:
			print("Too many initmodel folders in the project, or something odd happened....Exit.")
			exit()
		
	
	
	jsd=queue.Queue(0)	
	NTHREADS=max(options.threads+1,2)
	thrds=[threading.Thread(target=make_model,args=(jsd, i, options)) for i in range(options.ntry)]	
	thrtolaunch=0
	models=[]
	scrs=[]
	grads=[]
	cmps=[]
	
	while thrtolaunch<len(thrds) or threading.active_count()>1 or not jsd.empty():
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
		print(i,si,grads[si],scrs[si])
		models[si].write_image("{}/model_{:02d}.hdf".format(path, i))
		for c in cmps[si]:
			c.write_image("{}/model_{:02d}_aptcl.hdf".format(path, i), -1)
		
		
	
	#m, grad, scr=make_model(options)
	print("total time:", time.time()-tt0)


	E2end(logid)

	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	
