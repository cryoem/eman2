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
	parser.add_argument("--ppid", type=int,help="ppid...", default=-1)

	(options, args) = parser.parse_args()
	
	
	logid=E2init(sys.argv)
	
	cmd=args[0].split()
	pname=find_option(cmd, "--ptclsin")
	oname=find_option(cmd, "--ptclsout")
	dec=find_option(cmd, "--decoderout")
	enc=find_option(cmd, "--encoderout")
	midout=find_option(cmd, "--midout")
	citer=find_option(cmd, "--niter")
	path=pname[:pname.rfind('/')+1]
	print("Writing tmp files in ",path)
	
	rawcmd=args[0]
	for rawiter in range(options.niter+1):
		if rawiter==options.niter:
			rawcmd=rawcmd.replace(f"--niter {citer}", "--niter 0")
			
		lst=load_lst_params(pname)
		n=len(lst)
		print(f"List input {pname}, with {n} particles")
		
		batch=options.batch
		niter=ceil(n/batch)
		tmpfiles=[]
		
		for it in range(niter):
			ll=lst[it*batch:(it+1)*batch]
			print(it, len(ll))
			tmplst=path+f'tmp_input_{it:03d}.lst'
			tmpout=path+f'tmp_output_{it:03d}.lst'
			tmpmid=path+f'tmp_mid_{it:03d}.txt'
			save_lst_params(ll, tmplst)
			tmpfiles.append([tmplst,tmpout, tmpmid])
			
			cc=rawcmd.replace(pname, tmplst)
			if oname: cc=cc.replace(oname, tmpout)
			if midout: cc=cc.replace(midout, tmpmid)
			if options.load or rawiter>0: 
				cc=cc+f" --encoderin {enc} --decoderin {dec}"
					
			run(cc)
			
		if midout:
			midall=[]
			for tmp in tmpfiles:
				midall.append(np.loadtxt(tmp[2]))
			
			midall=np.concatenate(midall)
			midall[:,0]=np.arange(len(midall))
			print(midall.shape)
			np.savetxt(midout, midall)
				
		if oname:
			ptclsall=[]
			for tmp in tmpfiles:
				ptclsall.extend(load_lst_params(tmp[1]))
				
			save_lst_params(ptclsall, oname)
			
		options.load=True
		
		for tmp in tmpfiles:
			for f in tmp:
				try: os.remove(f)
				except: pass
		
	E2end(logid)
	
def find_option(cmd, op):
	
	ii=[i for i,c in enumerate(cmd) if c==op]
	if len(ii)!=1:
		print(f"error: cannot find {op}")
		return None
	
	ii=ii[0]
	return cmd[ii+1]
	
if __name__ == '__main__':
	main()
	
