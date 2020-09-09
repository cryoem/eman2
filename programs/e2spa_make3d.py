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
	
	parser.add_argument("--keep", type=float, dest="keep", help="The fraction of slices to keep, based on quality scores (1.0 = use all slices). See keepsig.",default=.9)
	
	parser.add_argument("--no_wt", action="store_true", dest="no_wt", default=False, help="This argument turns automatic weighting off causing all images to be weighted by 1.")
	
	parser.add_argument("--mode", type=str, default="trilinear", help="Fourier reconstruction 'mode' to use. The default should not normally be changed.")
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")

	parser.add_argument("--apix",metavar="A/pix",type=float,help="A/pix value for output, overrides automatic values",default=None)
	
	parser.add_argument("--ref", type=str,help="ref", default=None)
	parser.add_argument("--minres", type=float,help="", default=200)
	parser.add_argument("--maxres", type=float,help="", default=5)
	parser.add_argument("--parallel", type=str,help="Thread/mpi parallelism to use", default=None)
	
	parser.add_argument("--debug", action="store_true", default=False, help="")

	# Database Metadata storage
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	logger=E2init(sys.argv,options.ppid)
	
	# get basic image parameters
	tmp=EMData(options.input,0,True)
	boxsz=tmp["nx"]
	if options.apix!=None : apix=options.apix
	else : apix=tmp["apix_x"]

	if options.pad<0:
		options.pad=good_size(boxsz*1.5)
		
	if options.outsize<0:
		options.outsize=boxsz
		
	data=initialize_data(options.input, options)
	padvol=options.pad

	from EMAN2PAR import EMTaskCustomer
	if options.ref:
		print("weighting by reference...")
		ref=EMData(options.ref)
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
			
		wts=[0]*len(data)
		for i in tids:
			wt=etc.get_results(i)[1]
			for w in wt:
				wts[w[0]]=w[1]
				
		wts=np.array(wts)
		del etc
		
		r0=int(apix*boxsz/options.minres)
		r1=int(apix*boxsz/options.maxres)
		print(r0,r1)
		scrs=np.mean(wts[:,r0:r1], axis=1)
		if options.keep<1:
			thr=np.sort(scrs)[int(len(scrs)*(1-options.keep))-1]
			scrs[scrs<thr]=-1
		
		for i,d in enumerate(data):
			#d["curve"]=wts[i]
			d["weight"]=float(scrs[i])
			
			
		#print(wts.shape)
		#w=from_numpy(wts).copy()
		#w.write_image(options.output.replace("threed","wtall"))
		#return
	
	
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
		E2progress(logger, .5+.5*np.mean(st_vals)/100.)
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
	if options.verbose>0:
			print("Output File: "+options.output)

	E2end(logger)

	print("Exiting")

def initialize_data(inputfile, options):
	#returned list will contain dictionaries containing metadata about each image

	n_input=EMUtil.get_image_count(inputfile)
	lst=LSXFile(inputfile)
	print(n_input," input images")

	data=[]
	xfkey=["type","alt","az","phi","tx","ty","tz","alpha","scale"]

	for i in range(n_input):
		
		lstinfo=lst.read(i)
		for d in lstinfo[2].split(';'):
			dc=eval(d)
			if "score" in dc:
				score=dc["score"]
			else:
				score=-1
				
			dcxf={k:dc[k] for k in dc.keys() if k in xfkey}
				
			elem={
				"xform":Transform(dcxf),
				"score":score,
				"filename":inputfile,
				"filenum":i
			}
			data.append(elem)
			
	scrs=np.array([d["score"] for d in data])
	if np.min(scrs)>=0:
		print("positive score. assume this is weight...")
	else:
		scrs=-scrs
		scrs[scrs<0]=0
		scrs=scrs/np.max(scrs)
	
	print("score max {:.2f}, min {:.2f}".format(np.max(scrs), np.min(scrs)))
	if options.keep<1:
		thr=np.sort(scrs)[int(len(scrs)*(1-options.keep))-1]
		scrs[scrs<thr]=-1
		
	
	for i in range(n_input):
		data[i]["weight"]=float(scrs[i])
		
	return data


def reconstruct(data,recon,pad,ref=None):

	for i,elem in enumerate(data):
		wt=elem["weight"]
		if wt<=0:
			continue
		
		img=EMData(elem["filename"],elem["filenum"])
		if img["sigma"]==0 : continue
	
		img-=img.get_edge_mean()
		img.clip_inplace(Region((img["nx"]-pad)//2, (img["ny"]-pad)//2, pad, pad))
		
		if "curve" in elem:
			c=elem["curve"]#*0+1
			#ctf=img["ctf"]
			#ctf.bfactor=0
			#sz=len(c)*2
			#apix=img["apix_x"]
			#ctf=abs(np.array(ctf.compute_1d(sz,1./(apix*sz),Ctf.CtfType.CTF_AMP)))
			##ci=np.where(np.diff(ctf)<0)[0][0]
			##ctf[:ci]=1
			#c=ctf*c
			img.process_inplace("filter.radialtable", {"table":c.tolist()})
			#wt=1
		
		ts=Transform(elem["xform"])
		ts.set_rotation({"type":"eman"})
		ts.invert()   #### inverse the translation so make3d matches projection 
		img=recon.preprocess_slice(img,ts)
		

		recon.insert_slice(img,elem["xform"],wt)

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

		padvol=options.pad
		normvol=EMData(padvol//2+1, padvol, padvol)
		parms = {
			"size":[padvol]*3,
			"sym":options.sym,
			"mode":options.mode,
			"verbose":options.verbose-1,
			"normout":normvol
			}

		recon=Reconstructors.get("fourier", parms)
		recon.setup()

		reconstruct(
			data,
			recon,
			padvol,
			)

		output = recon.finish(False)

		#callback(100)

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
		wts=[]
		for i,elem in enumerate(data):
			img=EMData(elem["filename"],elem["filenum"])
			img=img.do_fft()
			img.process_inplace("xform.phaseorigin.tocenter")
			
			xf=Transform(elem["xform"])
			pj=ref.project('gauss_fft',{"transform":xf, "returnfft":1})

			fsc=img.calc_fourier_shell_correlation(pj)
			w=np.array(fsc).reshape((3,-1))[1]
			wts.append([elem["filenum"], w])
			
			#x=np.arange(len(w)*2)
			#p=argrelextrema(w, np.greater)[0]
			#p=p[w[p]>0]
			#p=p[p>4]
			#cfit = np.polyfit(x[p], np.log(w[p]), 1)
			#c=np.exp(x*cfit[0])*np.exp(cfit[1])
			#wts.append([elem["filenum"], c])
			
			
		
		return wts


if __name__=="__main__":
	main()
