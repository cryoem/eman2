#!/usr/bin/env python
# Muyuan Chen 2023-03

import numpy as np
import matplotlib.pyplot as plt
from EMAN2 import *

if "CUDA_VISIBLE_DEVICES" not in os.environ:
	# so we can decide which gpu to use with environmental variable
	os.environ["CUDA_VISIBLE_DEVICES"]='0' 

#### do not occupy the entire GPU memory at once
##   seems necessary to avoid some errors...
os.environ["TF_FORCE_GPU_ALLOW_GROWTH"]='true' 

#### finally initialize tensorflow
import tensorflow as tf
#### we will import some functions from e2gmm_refine later
emdir=e2getinstalldir()
sys.path.insert(0,os.path.join(emdir,'bin'))

#### need to unify the float type across tenforflow and numpy
##   in theory float16 also works but it can be unsafe especially when the network is deeper...
floattype=np.float32

#### load the Gaussian model produced by e2gmm_refine.py
from e2gmm_refine_new import *


def xf2pts_3d(args):
	pts, ang=args

	azp=-ang[0]
	altp=ang[1]
	phip=-ang[2]

	#### rotate Gaussian positions
	matrix=make_matrix(azp, altp, phip)
	matrix=tf.reshape(matrix, shape=[3,3]) 
	matrix=tf.transpose(matrix)

	pts_rot=tf.matmul(pts, matrix)

	#### finally do the translation
	pts_rot_trans=tf.stack([pts_rot[:,0]+ang[3], pts_rot[:,1]+ang[4], pts_rot[:,2]+ang[5]], 1)

	#pts_rot_trans=pts_rot_trans*sz+sz/2
	return pts_rot_trans


@tf.function()
def xf2pts_3d_mult(args):
	return tf.vectorized_map(xf2pts_3d, args)
	
def main():
	
	usage="."
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--info3d", type=str,help=" particle 3d info input ", default=None)
	parser.add_argument("--ptclsin", type=str,help="2d particle input ", default=None)
	parser.add_argument("--ptclsout", type=str,help="aligned particle output", default=None)
	parser.add_argument("--model", type=str,help="model file as reference", default=None)
	parser.add_argument("--maxres", type=float,help="max resolution", default=None)
	parser.add_argument("--minres", type=float,help="min resolution", default=100)
	parser.add_argument("--clip", type=int,help="size of 3d particle", default=None)
	parser.add_argument("--niter", type=int,help="number of iteration", default=20)
	parser.add_argument("--batchsz", type=int,help="batch size", default=8)
	parser.add_argument("--localwt", type=float,help="weight of nearby particles. larger means central particle has higher weight. default is 2.", default=2)
	parser.add_argument("--nnb", type=int,help="number of neighbors to consider", default=5)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	
	logid=E2init(sys.argv)
	
	pts=np.loadtxt(options.model).astype(floattype)
	print("Gaussian model shape: ", pts.shape)

	##   turn model to tensorflow format
	pts=tf.constant(pts)

	#### load the alignment from spt refinement
	alipm=load_lst_params(options.ptclsin)
	info3d=load_lst_params(options.info3d)
	print("Load {} 2D particles".format(len(alipm)))
	
	e=EMData(options.ptclsin,0, True)
	raw_apix, raw_boxsz = e["apix_x"], e["ny"]
	options.maxboxsz=ceil(options.clip*raw_apix*2/options.maxres)//2*2
	
	options.apix=raw_apix*options.clip/options.maxboxsz
	options.trainmodel=False
	options.minpx=ceil(options.clip*raw_apix*2/options.minres)//2*2
	options.minpx=max(1, options.minpx)
	options.maxpx=options.maxboxsz//2
	data_cpx, xfsnp=load_particles(options)
	print(f"particle size {raw_boxsz}, clip to {options.maxboxsz}, shrink to {options.maxboxsz}")
	print(f"compare {options.minpx} to {options.maxpx} Fourier pixels")

	clipid=set_indices_boxsz(data_cpx[0].shape[1], options.apix, True)
	params=set_indices_boxsz(options.maxboxsz)
	
	print("prepare nearest neighbor distance")
	alldst=[]
	allidx=[]
	for ap in alipm:
		si=[i for i,a in enumerate(alipm) if a["src"]==ap["src"] and a["tilt_id"]==ap["tilt_id"]]
		si=np.array(si)
		sel=[alipm[i] for i in si]
		crd=np.array([info3d[a["ptcl3d_id"]]["coord"] for a in sel])
		c0=info3d[ap["ptcl3d_id"]]["coord"]
		d=np.linalg.norm(crd-c0, axis=1)
		di=np.argsort(d)[:options.nnb]
		d=(crd[di]-c0)/options.clip
		if len(di)<options.nnb:
			### need to fill in empty columns to keep the matrix shape
			nn=options.nnb-len(di)
			di=np.append(di, [di[0]]*nn)
			d=np.vstack([d, 100+np.zeros((nn, 3))])
			
		alldst.append(d)
		allidx.append(si[di])
		
	alldst=np.array(alldst)
	allidx=np.array(allidx)
	print(alldst.shape)
	
	#### rewrite xfsnp since we need tz
	xfs=[p["xform.projection"].get_params("eman") for p in alipm]
	for x in xfs: x["tz"]=0
	xfsnp=np.array([[x["az"],x["alt"],x["phi"], x["tx"], x["ty"], x["tz"]] for x in xfs], dtype=floattype)
	xfsnp[:,:3]=xfsnp[:,:3]*np.pi/180.
	xfsnp[:,3:]/=float(options.clip)
	xfsnp[:,4]*=-1 ## inconsistence between eman-numpy coordinates
	
	#### prepare batches
	drg=np.arange(len(allidx))
	bsz=options.batchsz
	nbatch=len(drg)//bsz+1
	drg=[drg[i*bsz:(i+1)*bsz] for i in range(nbatch)]
	
	print("Starting alignment...")
	niter=options.niter
	allxfs=[]
	for iis in drg:
		pid=allidx[iis]
		dst=alldst[iis]
		wt=np.exp(-np.linalg.norm(dst, axis=2)*options.localwt)
		wt/=np.sum(wt, axis=1)[:,None]
		
		opt=tf.keras.optimizers.Adam(learning_rate=1e-3) 
		ptr=tf.gather(data_cpx[0], pid.flatten())
		ptj=tf.gather(data_cpx[1], pid.flatten())
		ptcl_cpx=(ptr, ptj)
		# print(ptr.shape)
		
		xf=xfsnp[pid].reshape((-1,6))
		
		# print(xf.shape)
		
		xf0=xf*0
		xfvar=tf.Variable(np.zeros((len(iis),6), dtype=floattype))
		p0=tf.constant(tf.zeros((xf.shape[0],pts.shape[0], 5))+pts)
		p1=xf2pts_3d_mult((p0[:,:,:3], xf))
		p1=p1*[-1,1,1]+dst.reshape((-1,1,3))
		p1shp=p1.shape
		# print(p1.shape)
		
		cost=[]
		for it in range(niter):
			with tf.GradientTape() as gt:
				p2=tf.reshape(p1, (len(iis), -1,3))
				p2=xf2pts_3d_mult((p2, xfvar))
				
				p3=tf.reshape(p2, p1shp)
				
				p3=p3-dst.reshape((-1,1,3))
				p3=p3*[-1,1,1]
				p3=tf.concat([p3, p0[:,:,3:]], axis=2)
				
				proj_cpx=pts2img(p3, xf0)
				# print(proj_cpx[0].shape)
				fval=calc_frc(proj_cpx, ptcl_cpx, params["rings"], minpx=options.minpx, maxpx=options.maxpx)
				fval=tf.reshape(fval, (len(iis), -1))
				# print(fval.shape, wt.shape)
				loss=-tf.reduce_sum(fval*wt, axis=1)
				loss=tf.reduce_mean(loss)
				# print(loss)
				
			grad=gt.gradient(loss, xfvar)
			opt.apply_gradients([(grad, xfvar)])
			cost.append(loss)
			sys.stdout.write(f"\r batch {len(allxfs)}/{len(allidx)}, iter {it}/{niter}, loss {cost[0]:.4f} -> {loss:.4f} ")
			sys.stdout.flush()
			
		allxfs.extend(xfvar.numpy().tolist())
		
	allxfs=np.array(allxfs)
	xnp=allxfs.copy()
	xnp[:,:3]=xnp[:,:3]*180./np.pi
	xnp[:,3:]*=options.clip
	xfs=[Transform({"type":"eman", "az":x[0], "alt":x[1], 
				"phi":x[2], "tx":x[3], "ty":x[4], "tz":x[5]*0}) for x in xnp.tolist()]
	xfs=[x.inverse() for x in xfs]
	xfs=[x.get_params("eman") for x in xfs]
	for x in xfs:x["alt"]*=-1
	xfs=[Transform(x) for x in xfs]

	aliout=[a.copy() for a in alipm]
	for a, xf in zip(aliout, xfs):
		a["xform.projection"]=xf*a["xform.projection"]

			
	print("output written to", options.ptclsout)
	save_lst_params(aliout, options.ptclsout)
	
	return
			
			
			
if __name__ == '__main__':
	main()
	
