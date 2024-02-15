#!/usr/bin/env python
# Muyuan Chen 2022-07
from EMAN2 import *
import numpy as np

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--batch", type=int, help="batch size", default=50000)
	parser.add_argument("--niter", type=int, help="number of iterations", default=2)
	parser.add_argument("--load", action="store_true", default=False ,help="load.")
	parser.add_argument("--ptcl3d", action="store_true", default=False ,help="separate by 3d particles. batch will correspond to 3d patricles")
	parser.add_argument("--subtilt", action="store_true", default=False ,help="separate particles for subtilt refinement")
	parser.add_argument("--ppid", type=int,help="ppid...", default=-1)

	(options, args) = parser.parse_args()
	
	
	logid=E2init(sys.argv)
	
	cmd=args[0].split()
	pname=find_option(cmd, "--ptclsin")
	oname=find_option(cmd, "--ptclsout")
	dec=find_option(cmd, "--decoderout")
	enc=find_option(cmd, "--encoderout")
	midin=find_option(cmd, "--midin")
	midout=find_option(cmd, "--midout")
	citer=find_option(cmd, "--niter")
	diter=find_option(cmd, "--retraindec")
	learnrate=find_option(cmd, "--learnrate")
	path=pname[:pname.rfind('/')+1]
	batch=options.batch
	print("Writing tmp files in ",path)
	# print(int(diter))
	
	rawcmd=args[0]
	for rawiter in range(options.niter+1):
		if midout and rawiter==options.niter:
			rawcmd=rawcmd.replace(f"--niter {citer}", "--niter 0")
			rawcmd=rawcmd.replace(f"--learnrate {learnrate}", "--learnrate 0")
			if diter!=None and int(diter)>0:
				rawcmd+=" --skipdec"
			
		lst=load_lst_params(pname)
		print(f"List input {pname}, with {len(lst)} particles")
		
		if options.ptcl3d and "ptcl3d_id" in lst[0]:
			pids=np.array([a["ptcl3d_id"] for a in lst])
			uid=np.unique(pids)
			print(f"{len(uid)} 3d particles")
			p3did=[np.where(pids==u)[0] for u in uid]
			nbatch=ceil(len(uid)/batch)
		
		elif options.subtilt:
			srcs=np.unique([i["src"] for i in lst])
			p3d=np.unique([i["ptcl3d_id"] for i in lst])
			print(f"{len(srcs)} tomograms, {len(p3d)} 3d particles, {len(lst)} 2d particles")
			
			batches=[]
			current=[]
			for src in srcs:
				l2=[i for i in lst if i['src']==src]
				tids=np.unique([l["tilt_id"] for l in l2])
				for t in tids:
					l3=[l for l in l2 if l["tilt_id"]==t]
					if len(l3)+len(current)>options.batch:
						batches.append(current)
						current=[]
					current.extend(l3)
			
			batches.append(current)
			nbatch=len(batches)
			for i,b in enumerate(batches):print(i,len(b))
					
		else:
			options.ptcl3d=False		
			nbatch=ceil(len(lst)/batch)
			
		tmpfiles=[]
		if midin!=None:
			mid=np.loadtxt(midin)
			print(f"load --midin from {midin}, shape {mid.shape}")
		
		for it in range(nbatch):
			print("#######################")
			print(f"batch {it+1}/{nbatch}...")
			if options.ptcl3d:
				ul=p3did[it*batch:(it+1)*batch]
				ll=np.concatenate(ul).tolist()
				ll=[lst[i] for i in ll]
			elif options.subtilt:
				ll=batches[it]
			else:
				ll=lst[it*batch:(it+1)*batch]				
				
			print(f"{it}/{nbatch} batches, {len(ll)} particles")
			tmplst=path+f'tmp_input_{it:03d}.lst'
			tmpout=path+f'tmp_output_{it:03d}.lst'
			tmpmid=path+f'tmp_mid_{it:03d}.txt'
			save_lst_params(ll, tmplst)
			tmps=[tmplst,tmpout, tmpmid]
			cc=rawcmd.replace(pname, tmplst)
			if oname: cc=cc.replace(oname, tmpout)
			if midout: cc=cc.replace(midout, tmpmid)
			
			if midin!=None:
				tmpmidin=path+f'tmp_midin_{it:03d}.txt'
				mm=mid[it*batch:(it+1)*batch,:]
				np.savetxt(tmpmidin, mm)
				tmps.append(tmpmidin)
				cc=cc.replace(midin, tmpmidin)

			tmpfiles.append(tmps)
			
			if options.load or rawiter>0: 
				if enc!=None:
					cc=cc+f" --encoderin {enc}"
				if dec!=None:
					cc=cc+f" --decoderin {dec}"
					
			run(cc)
			
		if midout:
			midall=[]
			for tmp in tmpfiles:
				midall.append(np.loadtxt(tmp[2]))
			
			midall=np.concatenate(midall)
			midall[:,0]=np.arange(len(midall))
			print(f"final latent space saved to {midout}, shape {midall.shape}")
			np.savetxt(midout, midall)
				
		if oname:
			if "e2gmm_rigidbody.py" in rawcmd:
				#### the output lst will be broken into multiple patches
				for ip in range(100):
					om=f'{path}/tmp_output_p{ip:02d}_000.lst'
					if not os.path.isfile(om): break
				npatch=ip
				
				for ip in range(npatch):
					ptclsall=[]
					tmps=[f'{path}/tmp_output_p{ip:02d}_{it:03d}.lst' for it in range(nbatch)]
					for tmp in tmps:
						ptclsall.extend(load_lst_params(tmp))
						os.remove(tmp)
					
					o=oname.split('_')
					o.insert(-1,f'p{ip:02d}')
					om='_'.join(o)
					save_lst_params(ptclsall, om)
				
			else:
				ptclsall=[]
				for tmp in tmpfiles:
					ptclsall.extend(load_lst_params(tmp[1]))
					
				save_lst_params(ptclsall, oname)
			
		options.load=True
		
		for tmp in tmpfiles:
			for f in tmp:
				try: os.remove(f)
				except: pass
		
	print("batch processing finished.")
	if midout:
		print(f"latent space output saved to {midout}")
	if oname:
		print(f"aligned particles saved to {oname}")
	E2end(logid)
	
def find_option(cmd, op):
	
	ii=[i for i,c in enumerate(cmd) if c==op]
	if len(ii)!=1:
		print(f"e2gmm_batch: cannot find {op}")
		return None
	
	ii=ii[0]
	return cmd[ii+1]
	
if __name__ == '__main__':
	main()
	
