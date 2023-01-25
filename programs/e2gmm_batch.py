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
	(options, args) = parser.parse_args()
	
	
	logid=E2init(sys.argv)
	
	cmd=args[0].split()
	pname=find_option(cmd, "--ptclsin")
	dec=find_option(cmd, "--decoderout")
	enc=find_option(cmd, "--encoderout")
	midout=find_option(cmd, "--midout")
	citer=find_option(cmd, "--niter")
	path=midout[:midout.rfind('/')+1]
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
		midall=[]
		for it in range(niter):
			ll=lst[it*batch:(it+1)*batch]
			print(it, len(ll))
			tmplst=path+f'tmp_input_{it:03d}.lst'
			tmpmid=path+f'tmp_mid_{it:03d}.txt'
			save_lst_params(ll, tmplst)
			tmpfiles.extend([tmplst, tmpmid])
			
			cc=rawcmd.replace(pname, tmplst)
			cc=cc.replace(midout, tmpmid)
			if options.load or it>0: 
				cc=cc+f" --encoderin {enc} --decoderin {dec}"
				
				
			run(cc)
			midall.append(np.loadtxt(tmpmid))
			
		midall=np.concatenate(midall)
		midall[:,0]=np.arange(len(midall))
		print(midall.shape)
		np.savetxt(midout, midall)
				
		options.load=True
		
		for f in tmpfiles:
			try: os.remove(f)
			except: pass
		
	E2end(logid)
	
def find_option(cmd, op):
	
	ii=[i for i,c in enumerate(cmd) if c==op]
	if len(ii)!=1:
		print(f"error: cannot find {op}")
	
	ii=ii[0]
	return cmd[ii+1]
	
if __name__ == '__main__':
	main()
	
