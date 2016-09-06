#!/usr/bin/env python
# Muyuan Chen 2016-08
from EMAN2 import *
import numpy as np
from Simplex import Simplex
from multiprocessing import Process, Queue
import time
import gc

def main():
	
	usage="align sub-tomogram by projections. called by refinetomo_easy.py "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--pjfile", type=str,help="projection file", default=None)
	parser.add_argument("--ptcl", type=str,help="particle file", default=None)
	parser.add_argument("--simmx", type=str,help="similarity matrix file", default=None)
	parser.add_argument("--lstin", type=str,help="input list file", default=None)
	parser.add_argument("--mapfile", type=str,help="initial map file", default=None)
	parser.add_argument("--lstout", type=str,help="output list file", default=None)
	parser.add_argument("--threads", type=int,help="number of threads", default=12)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	tstart=time.time()
	
	#### read orientation from projections
	print "Reading projections..."
	if options.pjfile:
		if not options.simmx:
			print "No simmx input..."
			exit()
			
		projfile=options.pjfile
		npj=EMUtil.get_image_count(projfile)
		proj_oris=[]
		for i in range(npj):
			e=EMData(projfile,i, True)
			tr=e["xform.projection"]
			proj_oris.append(tr)
		proj_oris=np.asarray(proj_oris)
	
	#### read info from particles
	print "Reading partiles..."
	ptclfile=options.ptcl
	npt=EMUtil.get_image_count(ptclfile)
	#npt=1000
	modelid=[]
	ptcl_ori=[]
	for i in range(npt):
		e=EMData(ptclfile,i, True)
		modelid.append(e["model_id"])
		ptcl_ori.append(e["xform.projection"])
	modelid=np.array(modelid)
	ptcl_ori=np.asarray(ptcl_ori)
	
	#### read simmx
	if options.simmx:
		print "Parsing simmilarity matrix..."
		simxfile=options.simmx
		simmxall=[]
		for i in range(6):
			ss=EMData(simxfile,i)
			simmxall.append(ss.numpy().copy())
		simmxall=np.asarray(simmxall)
		print "simmx shape:",simmxall.shape
		
		simnp=simmxall[0]
	else:
		if not options.lstin:
			print "Input should contain either simmx or lstin"
			exit()
		simmxall=None
		simnp=None
		lstin=LSXFile(options.lstin, True)
	
	e3d=EMData(options.mapfile)
	np3d=e3d.numpy().copy()
	mapft=get_fft(np3d)
	sz=mapft.shape[0]
	sli=np.indices((sz,sz))-sz/2
	sli=np.stack([sli[0], sli[1], np.zeros_like(sli[0])])
	
	x,y= np.indices((sz,sz))
	rr=np.sqrt((x - sz/2)**2 + (y - sz/2)**2).astype(int)
	rings=np.zeros((sz,sz,sz/2))
	for i in range(sz/2):
		rings[:,:,i]=(rr==i)

	
	#### refine orientation
	print "Refining orientation..."
	allmid=np.unique(modelid)
	jobs=[]
	inputs=[]
	rets=[]
	for mid in allmid:
		curidx=np.where(modelid==mid)[0]
		job={}
		if options.simmx:
			armin=np.argmin(simnp[curidx],axis=1)
			idx=np.argmin(simnp[curidx, armin])
			pjori=proj_oris[armin]
			cor=simnp[curidx, armin]
			idx=np.argmin(cor) #### index of the best slice
			job["usesimmx"]=True
			job["idx"]=idx 
			job["simmx"]=simmxall
			job["amin"]=armin[idx]
			job["proj_oris"]=proj_oris
		
		else:
			job["usesimmx"]=False
			tr=lstin.read(curidx[0])[2]
			job["trans0"]=eval(tr)
			
			
		job["modelid"]=mid
		job["rings"]=rings
		job["curidx"]=curidx
		job["ptcl_ori"]=ptcl_ori
		job["ptclfile"]=ptclfile
		job["sli"]=sli
		job["mapft"]=mapft
		job["3dmap"]=e3d
		que=Queue()
		prs=Process(target=refine_align, args=(job,que))
		jobs.append(prs)
		inputs.append(job)
		rets.append(que)
		#print mid
	
	#### hand built queue to avoid the memory leaking in from_numpy()
	numleft=len(jobs)
	running=[-1]*options.threads
	kj=0
	while(numleft>0):
		for i in range(options.threads):
			if running[i]>=0:
				if not jobs[running[i]].is_alive():
					#### job finished
					numleft-=1
					running[i]=-1
			if running[i]<0 and kj<len(jobs):
				#### start a job
				running[i]=kj
				kj+=1
				jobs[running[i]].start()
				#print "start: ", running[i], " on worker", i
		print  "Waiting for", numleft, "tasks to complete..."
		time.sleep(2)
	
	#### write output
	print "Writing output to ", options.lstout
	try: os.remove(options.lstout)
	except: pass
	
	lstfile=LSXFile(options.lstout, False)
	for oi in range(len(jobs)):
		olst=rets[oi].get()
		#print olst
		job=inputs[oi]
		curidx=job["curidx"]
		otr=Transform({"alt":olst[3], "az":olst[4], "phi":olst[5], "type":"eman"})
		vv=[olst[0], olst[1], olst[2]]
		for i in range(len(curidx)):
			tr= otr*ptcl_ori[curidx[i]]
			v=tr.transform(vv)
			tr.set_trans(v[0], v[1], v[2])
			lstfile.write(-1, curidx[i], ptclfile, tr.get_params("eman"))
		
		
	lstfile=None
	
	print time.time()-tstart
	E2end(logid)
	
def run(cmd):
	print cmd
	launch_childprocess(cmd)
	
def get_fft(img):
	return np.fft.fftshift(np.fft.fftn(np.fft.fftshift(img)))

def get_img(fft):
	return np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(fft)).real)	


def get_transform_simx(pti,pji, simmxall, proj_oris):
	pjo=proj_oris[pji].get_params('eman')
	dx=float(simmxall[1,pti,pji])
	dy=float(simmxall[2,pti,pji])
	dalpha=float(simmxall[3,pti,pji])
	dmir=simmxall[4,pti,pji]
	tr=Transform({"type":"2d", "alpha":dalpha, "tx": dx, "ty":dy, "mirror":int(dmir>0)})
	tr=tr.inverse().get_params("eman")
	pjo["phi"]=tr["phi"]
	pjo["tx"]=tr["tx"]
	pjo["ty"]=tr["ty"]
	pjo["mirror"]=tr["mirror"]
	return Transform(pjo)

def refine_align(job,ret):
	clskeep=.7
	#print job
	curidx=job["curidx"]
	ptcl_ori=job["ptcl_ori"]
	rings=job["rings"]
	ptclfile=job["ptclfile"]
	sli=job["sli"]
	mapft=job["mapft"]
	e3d=job["3dmap"]
	sz=mapft.shape[0]
	
	if job["usesimmx"]:
		idx=job["idx"]
		proj_oris=job["proj_oris"]
		oo=get_transform_simx(curidx[idx],job["amin"], job["simmx"],proj_oris)
		oi=oo*ptcl_ori[curidx[idx]].inverse()
		
	else:
		tr0=Transform(job["trans0"])
		oi=tr0*ptcl_ori[curidx[0]].inverse()
	
	oi=oi.get_params("eman")
	oilst=[oi["tx"], oi["ty"], oi["tz"], oi["alt"], oi["az"], oi["phi"]]
	
	imgs=[]
	for i in curidx:
		e=EMData(ptclfile,i)
		e.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/30.})
		e.process_inplace("normalize")
		imgs.append(e)
	# score=0
	
	def get_score(olst,data):
		score=[]
		oi=Transform({"alt":olst[3], "az":olst[4], "phi":olst[5], "type":"eman"})
		vv=[olst[0], olst[1], olst[2]]
		
		for i in range(len(curidx)):
			
			tr= oi*ptcl_ori[curidx[i]]
			
			#pp= e3d.project("standard", tr)
			
			surf= np.tensordot(np.asarray(tr.get_matrix_4x4()).reshape(4,4)[:3,:3], sli, axes=(0,0))
			ind=(surf+sz/2).astype(int)
			ind=np.clip(ind,0, sz-1)
			imft=mapft[ind[2], ind[1], ind[0]]
			v=tr.transform(vv)
			img=get_img(imft).T
			pp=from_numpy(img)
			pp.translate(v[0],v[1],0)
			#score.append(imgs[i].cmp("frc", pp,{'maxres':30}))
			score.append(imgs[i].cmp("ccc", pp))
			
			#e=imgs[i].copy()
			#e.translate(-v[0],-v[1],0)
			
			#tt=get_fft(e.numpy())
			#xreal=imft.real
			#ximag=imft.imag
			#treal=tt.real
			#timag=tt.imag
			
			

			#it1=xreal**2+ximag**2
			#it2=treal**2+timag**2
			#it1ring=np.sqrt(np.tensordot(it1,rings))
			#it2ring=np.sqrt(np.tensordot(it2,rings))
			#nrm=it1ring*it2ring


			#loss= - (np.tensordot((xreal*treal) + (ximag*timag),rings)/nrm)
			#score.append(np.mean(loss[:10]))
		
		#score/=float(len(curidx))
		score=np.sort(score)[:int(len(curidx)*clskeep)]
		score=np.mean(score)
		#print score
		
		return score
	
	incr=[20,20,20,90,180,180]
	
	#for k in range(10000): 
		#if k%100==0: print k/100
		#get_score(oilst,0)
	simp=Simplex(get_score,np.array(oilst), incr)
	
	locs=simp.minimize(maxiters=100,epsilon=.001,monitor=0)

	olst=locs[0]
	#gc.collect()
	print job["modelid"], locs[1], locs[2]
	ret.put(olst)
	return olst

if __name__ == '__main__':
	main()
	
