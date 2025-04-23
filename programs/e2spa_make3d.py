#!/usr/bin/env python
from EMAN2 import *
import numpy as np
from EMAN2jsondb import JSTask
from scipy.signal import argrelextrema

def main():
	parser = EMArgumentParser(usage="")

	parser.add_argument("--output", default="threed.hdf", help="Output reconstructed volume file name.")
	parser.add_argument("--input", default=None, help="The input projections. Project should usually have the xform.projection header attribute, which is used for slice insertion")
	parser.add_argument("--sym", dest="sym", default="c1", help="Set the symmetry; if no value is given then the model is assumed to have no symmetry.\nChoices are: i, c, d, tet, icos, or oct.")

	parser.add_argument("--pad", default=-1,type=int, help="Will zero-pad images to the specifed size. ")
	
	parser.add_argument("--outsize", default=-1, type=int, help="Defines the dimensions of the final volume written to disk")
	
	parser.add_argument("--keep", type=str, dest="keep", help="The fraction of slices to keep, based on quality scores (1.0 = use all slices).",default=".9")
	
	parser.add_argument("--no_wt", action="store_true", dest="no_wt", default=False, help="This argument turns automatic weighting off causing all images to be weighted by 1.")
	
	parser.add_argument("--mode", type=str, default="trilinear", help="Fourier reconstruction 'mode' to use. The default should not normally be changed.")
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")

	parser.add_argument("--apix",metavar="A/pix",type=float,help="A/pix value for output, overrides automatic values",default=None)
	parser.add_argument("--tidrange", type=str,help="Range of tilt id to include for particles from tilt series. Specify two integers separated by ','.", default="-1,-1")
	
	parser.add_argument("--ref", type=str,help="Weight each particle using a specified reference map.", default=None)
	parser.add_argument("--minres", type=float,help="minimum resolution to compare when weighting by a reference map.", default=50)
	parser.add_argument("--maxres", type=float,help="maximum resolution to compare when weighting by a reference map.", default=-1)
	parser.add_argument("--parallel", type=str,help="Thread/mpi parallelism to use without shared memory. Each worker will reconstruct a map with a subset of particles and the results from workers will be averaged together with the corresponding Fourier weighting. Along with --threads, this allows having one worker per node using multiple threads.", default="thread:1")
	parser.add_argument("--threads", type=int,help="Number of threads using shared memory.", default=1)
	parser.add_argument("--setsf", type=str,help="Set structure factor from text file", default=None)
	
	parser.add_argument("--debug", action="store_true", default=False, help="Turn on debug mode. This will only process a small subset of the data.")
	parser.add_argument("--clsid", default=None, type=str, help="Only reconstruct a class of particles. Also take even/odd to reconstruct subsets of particles.")
	parser.add_argument("--listsel", default=None, type=str, help="only reconstruct particles of indices from the given list in a text file.")

	parser.add_argument("--p3did", type=int, help="only reconstruct images for one 3d particle of the given ID.",default=-1)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	logger=E2init(sys.argv,options.ppid)
	time0=time.time()

	# so it recognize even/odd or 0/1
	if options.clsid:
		options.clsid=options.clsid.replace("even","0").replace("odd","1")
		try: options.clsid=int(options.clsid)
		except:options.clsid=-1
	else:
		options.clsid=-1
		
	if options.sym.lower().startswith("h"):
		options.sym="c1"
		
	options.keep=[float(k) for k in options.keep.split(',')]
	if len(options.keep)<3: options.keep=[options.keep[0]]*3
	
	# get basic image parameters
	tmp=EMData(options.input,0,True)
	boxsz=tmp["nx"]
	if options.apix!=None : apix=options.apix
	else : apix=tmp["apix_x"]

	if options.outsize<0:
		options.outsize=boxsz
		
	if options.pad<0:
		options.pad=good_size(options.outsize*1.5)
		
		
	options.tidrange=[int(i) for i in options.tidrange.split(',')]
	if options.tidrange[0]>=0:
		print("including tilt ids from {} to {}".format(options.tidrange[0], options.tidrange[1]))
	
	data=initialize_data(options.input, options)
	padvol=options.pad
	prog=0.
	from EMAN2PAR import EMTaskCustomer
	if options.ref:
		print("weighting by reference...")
		ref=EMData(options.ref)
		ny=options.pad
		by=ref["ny"]
		ref=ref.get_clip(Region((by-ny)/2, (by-ny)/2,(by-ny)/2, ny, ny,ny))
		ref=ref.do_fft()
		ref.process_inplace("xform.phaseorigin.tocenter")
		ref.process_inplace("xform.fourierorigin.tocenter")
		
		etc=EMTaskCustomer(options.parallel, module="e2spa_make3d.WeightptclTask")
		num_cpus = etc.cpu_est()
		print("{} total CPUs available".format(num_cpus))
		tasks=[data[i::num_cpus] for i in range(num_cpus)]
		print("{} jobs".format(len(tasks)))
		
		tids=[]
		for t in tasks:
			task = WeightptclTask(t, ref, options)
			if options.debug:
				task.execute(print)
				return
			tid=etc.send_task(task)
			tids.append(tid)
		
		while 1:
			st_vals = etc.check_task(tids)
			E2progress(logger, .5*np.mean(st_vals)/100.)
			if np.min(st_vals) == 100: break
			time.sleep(5)
		prog=0.5
		
		wts=[0]*len(data)
		for i in tids:
			wt=etc.get_results(i)[1]
			for w in wt:
				try:
					wts[w[0]]=w[1]
				except:
					print(w)
				
		
		wts=np.array(wts)
		
		r0=0; r1=-1
		if options.minres>0:
			r0=int(apix*ny/options.minres)
			wts[:,:r0]=np.mean(wts[:,:r0], axis=0)
# 		
		if options.maxres>0:
			r1=int(apix*ny/options.maxres)
			wts[:,r1:]=np.mean(wts[:,r1:], axis=0)

		scr=np.mean(wts[:,r0:r1], axis=1)
		
		print(wts.shape, scr.shape)
		print(np.min(scr), np.mean(scr), np.max(scr))
		del etc
		
		#print(r0,r1)
		scrs=np.mean(wts[:,r0:r1], axis=1)
		if options.keep[1]<1:
			thr=np.sort(scrs)[int(len(scrs)*(1-options.keep[1]))-1]
			scrs[scrs<thr]=-1
		
		for i,d in enumerate(data):
			d["curve"]=wts[i]
			d["weight"]=float(scrs[i])
			
			
	etc=EMTaskCustomer(options.parallel, module="e2spa_make3d.Make3dTask")
	num_cpus = etc.cpu_est()
	print("{} total CPUs available".format(num_cpus))
	tasks=[data[i::num_cpus] for i in range(num_cpus)]
	print("{} jobs".format(len(tasks)))

	tids=[]
	for t in tasks:
		task = Make3dTask(t, options)
		if options.debug:
			task.execute(print)
			return
		tid=etc.send_task(task)
		tids.append(tid)

	while 1:
		st_vals = etc.check_task(tids)
		E2progress(logger, prog+(1-prog)*np.mean(st_vals)/100.)
		if np.min(st_vals) == 100: break
		time.sleep(5)
		
	output=EMData(padvol, padvol, padvol)
	normvol=EMData(padvol//2+1, padvol, padvol)
	output.to_zero()
	output.do_fft_inplace()
	normvol.to_zero()
	
	for i in tids:
		threed, norm=etc.get_results(i)[1]
		threed.process_inplace("math.multamplitude", {"amp":norm})
		output.add(threed)
		normvol.add(norm)
		
	normvol.process_inplace("math.reciprocal")
	output.process_inplace("math.multamplitude", {"amp":normvol})
	
	output.do_ift_inplace()
	output.depad()
	output.process_inplace("xform.phaseorigin.tocenter")
	
	del etc

	if options.verbose>0 : print("Finished Reconstruction")
	output["apix_x"]=output["apix_y"]=output["apix_z"]=apix
	sz=options.outsize
	output.clip_inplace(Region((padvol-sz)//2,(padvol-sz)//2,(padvol-sz)//2,sz,sz,sz))
	
	if os.path.isfile(options.output):
		os.remove(options.output)
		
	output.write_image(options.output,0)
	if options.setsf:
		launch_childprocess("e2proc3d.py {} {} --setsf {}".format(options.output,options.output,options.setsf))
		
	if options.verbose>0:
			print("Output File: "+options.output)

	E2end(logger)

	print("Reconstruction finishend ({:.1f} s)".format(time.time()-time0))

def initialize_data(inputfile, options):
	#returned list will contain dictionaries containing metadata about each image
	if options.listsel:
		sel=np.loadtxt(options.listsel)
		sel=sel.astype(int).tolist()
	else:
		sel=[]
		
	if inputfile.endswith(".lst"):
		rawdata=load_lst_params(inputfile, sel)
	else:
		n=EMUtil.get_image_count(inputfile)
		rawdata=[{"src":inputfile, "idx":i} for i in range(n)]
	trg=options.tidrange
	
	data=[]
	for dc in rawdata:
		if (options.clsid>=0) and ("class" in dc): 
			if dc["class"]!=options.clsid:
				continue
			
		if ("tilt_id" in dc) and (trg[0]>0):
			if (dc["tilt_id"]<trg[0]) or (dc["tilt_id"]>trg[1]):
				continue
			
		if ("ptcl3d_id" in dc) and (options.p3did>=0):
			if dc["ptcl3d_id"]!=options.p3did:
				continue
		
		data.append(dc)
		
	tokeep=np.ones(len(data), dtype=bool)
	print("{} particles total".format(np.sum(tokeep)))
	
	scrs=np.array([d["score"] if "score" in d else 2e5 for d in data])
	if np.max(scrs)>1e5:
		scrs=np.zeros(len(data))-1
	
	if np.std(scrs)>0 and ("ptcl3d_id" in data[0]):
		idx3d=np.array([d["ptcl3d_id"] for d in data])
		idx3d, invid=np.unique(idx3d, return_inverse=True)
		scr3d=np.array([np.mean(scrs[invid==i]) for i in range(len(idx3d))])
		thr=np.sort(scr3d)[int(len(scr3d)*options.keep[0])-1]
		for i,s in enumerate(scr3d):
			if s>thr:
				tokeep[invid==i]=False
		
		print("  Excluding particles on 3D particle score. Now {:.1f}%  left".format(100*np.mean(tokeep)))
	
	if np.std(scrs)>0:
		s=scrs[tokeep]
		thr=np.sort(s)[int(len(s)*options.keep[1])-1]
		tokeep[scrs>thr]=False
		print("  Excluding particles on 2D score. Now {:.1f}%  left".format(100*np.mean(tokeep)))
		
	if  "dxf" in data[0]:
		dt=np.array([d["dxf"].get_trans() for d in data])
		dt=np.linalg.norm(dt, axis=1)
		s=dt[tokeep]
		thr=np.sort(s)[int(len(s)*options.keep[2])-1]
		tokeep[dt>thr]=False
		print("  Excluding particles on translation. Now {:.1f}%  left".format(100*np.mean(tokeep)))
		
	data=[d for i,d in enumerate(data) if tokeep[i]]
	scrs=scrs[tokeep]
	if np.all(scrs==0):
		scrs+=1
	if np.min(scrs)>=0:
		print("positive score. assume this is weight...")
	else:
		scrs=-scrs
		scrs[scrs<0]=0
		scrs=scrs/np.max(scrs)
	
	if options.no_wt:
		scrs=scrs*0+1
	print("{} particle, score max {:.2f}, min {:.2f}".format(len(data),np.max(scrs), np.min(scrs)))
	
	for i in range(len(data)):
		data[i]["weight"]=float(scrs[i])
		data[i]["ii"]=i
		
	return data


def reconstruct(data,recon,pad,ref=None):

	for i,elem in enumerate(data):
		wt=elem["weight"]
		if wt<=0:
			continue
		
		try: img=EMData(elem["src"],elem["idx"])
		except:
			print(f"Couldn't load {elem['src']} {elem['idx']}. Skipping")
			continue
		if img["sigma"]==0 : continue
	
		img-=img.get_edge_mean()
		img.clip_inplace(Region((img["nx"]-pad)//2, (img["ny"]-pad)//2, pad, pad))
	
		if ("defocus" in elem) and abs(elem["defocus"])>1e-4 and img.has_attr("ctf"):
			ctf=img["ctf"]
			fft1=img.do_fft()
			flipim=fft1.copy()
			ctf.compute_2d_complex(flipim,Ctf.CtfType.CTF_SIGN)
			fft1.mult(flipim)
			ctf1=EMAN2Ctf(ctf)
			ctf1.defocus=ctf1.defocus+elem["defocus"]
			ctf1.compute_2d_complex(flipim,Ctf.CtfType.CTF_SIGN)
			fft1.mult(flipim)
			img=fft1.do_ift()
	
		img-=img.get_edge_mean()
		
		if "curve" in elem:
			c=elem["curve"]
			img.process_inplace("filter.radialtable", {"table":c.tolist()})
		
		if "xform.projection" in elem:
			pjxf=elem["xform.projection"]
		else:
			pjxf=img["xform.projection"]
			
		ts=Transform(pjxf)
		ts.set_rotation({"type":"eman"})
		ts.invert()   #### inverse the translation so make3d matches projection 
		img=recon.preprocess_slice(img,ts)
		

		recon.insert_slice(img,pjxf,wt)

	return



class Make3dTask(JSTask):


	def __init__(self, inp, options):

		data={"data":inp}
		JSTask.__init__(self,"Make3d",data,{},"")
		self.options=options


	def execute(self, callback):

		callback(0)
		data=self.data["data"]
		options=self.options
		Util.init_threads(options.threads)
		padvol=options.pad
		normvol=EMData(padvol//2+1, padvol, padvol)
		parms = {
			"size":[padvol]*3,
			"sym":options.sym,
			"mode":options.mode,
			"verbose":0,
			"quiet":1,
			"normout":normvol
			}

		recon=Reconstructors.get("fourier", parms)
		recon.setup()
		
		threads=[threading.Thread(target=reconstruct,
			    args=(data[i::options.threads],recon, padvol)) for i in range(options.threads)]

		for t in threads: t.start()
		for t in threads: t.join()

		#reconstruct(
			#data,
			#recon,
			#padvol,
			#)

		output = recon.finish(False)

		return (output, normvol)



class WeightptclTask(JSTask):
	
	
	def __init__(self, inp, ref, options):
		
		data={"data":inp, "ref":ref}
		JSTask.__init__(self,"Weightptcl",data,{},"")
		self.options=options
	
	
	def execute(self, callback):
		
		callback(0)
		data=self.data["data"]
		ref=self.data["ref"]
		options=self.options
		pad=options.pad
		wts=[]
		for i,elem in enumerate(data):
			img=EMData(elem["src"],elem["idx"])
			img-=img.get_edge_mean()
			img.clip_inplace(Region((img["nx"]-pad)//2, (img["ny"]-pad)//2, pad, pad))
			
			img=img.do_fft()
			img.process_inplace("xform.phaseorigin.tocenter")
			
			xf=Transform(elem["xform.projection"])
			pj=ref.project('gauss_fft',{"transform":xf, "returnfft":1})

			fsc=img.calc_fourier_shell_correlation(pj)
			w=np.array(fsc).reshape((3,-1))[1]
			
			#x=np.arange(len(w)*2)
			#p=argrelextrema(w, np.greater)[0]
			#p=p[w[p]>0]
			#p=p[p>4]
			#cfit = np.polyfit(x[p], np.log(w[p]), 1)
			#w=np.exp(x*cfit[0])*np.exp(cfit[1])
			#wts.append([elem["filenum"], c])
			
			wts.append([elem["ii"], w])
			
			if options.debug:
				print(elem["ii"], w, np.mean(w[2:70]))
		
		return wts


if __name__=="__main__":
	main()
