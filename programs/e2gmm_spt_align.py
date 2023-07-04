#!/usr/bin/env python
# Muyuan Chen 2023-03

import numpy as np
import matplotlib.pyplot as plt
from EMAN2 import *

if "CUDA_VISIBLE_DEVICES" not in os.environ:
	# so we can decide which gpu to use with environmental variable
	os.environ["CUDA_VISIBLE_DEVICES"]='1' 
	
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


def xf2pts_3d(pts, ang):
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

def main():
	
	usage="."
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ptclsin", type=str,help="2d particle input ", default=None)
	parser.add_argument("--ptclsout", type=str,help="aligned particle output", default=None)
	parser.add_argument("--model", type=str,help="model file as reference", default=None)
	parser.add_argument("--maxres", type=float,help="max resolution", default=None)
	parser.add_argument("--minres", type=float,help="min resolution", default=100)
	parser.add_argument("--clip", type=int,help="size of 3d particle", default=None)
	parser.add_argument("--niter", type=int,help="number of iteration", default=40)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	
	logid=E2init(sys.argv)
	
	pts=np.loadtxt(options.model).astype(floattype)
	print("Gaussian model shape: ", pts.shape)

	##   turn model to tensorflow format
	pts=tf.constant(pts)

	#### load the alignment from spt refinement
	##   make sure we have only one conformation for multiple 2D particles that correspond to one 3D particle
	alipm=load_lst_params(options.ptclsin)
	pids=np.array([a["ptcl3d_id"] for a in alipm])
	uid=np.unique(pids)
	p3did=[np.where(pids==u)[0] for u in uid]

	##   here we only use the near center tilts (-20 to +20 degrees)
	#tids=np.array([a["tilt_id"] for a in alipm])
	# tmid=np.sort(np.unique(tids))
	# tmid=tmid[len(tmid)//2]
	# trg=abs(tids-tmid)<90

	##   keep track of 3D particles from the 2D ones
	idx=np.unique(pids)
	idx3d=[]
	for i in idx:
		idx3d.append(np.where(pids==i)[0])

	print("{} 3D particles, each contain {} 2D particles".format(len(idx3d), len(idx3d[0])))
	
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
	dcpx=get_clip(data_cpx, params["sz"], clipid)

	
	niter=options.niter
	allxfs=[]
	for pid in p3did:
		opt=tf.keras.optimizers.Adam(learning_rate=1e-4) 
		ptr=tf.gather(data_cpx[0], pid)
		ptj=tf.gather(data_cpx[1], pid)
		ptcl_cpx=(ptr, ptj)
		xf=xfsnp[pid]
		xfvar=tf.Variable(np.zeros(6, dtype=floattype))
		
		cost=[]
		for it in range(niter):
			with tf.GradientTape() as gt:
				p1=xf2pts_3d(pts[:,:3], xfvar)
				p1=tf.concat((p1, pts[:,3:]), axis=1)
				proj_cpx=pts2img(p1[None,:,:], xf)
				fval=calc_frc(proj_cpx, ptcl_cpx, params["rings"], minpx=options.minpx, maxpx=options.maxpx)
				loss=-tf.reduce_mean(fval)
				
			grad=gt.gradient(loss, xfvar)
			opt.apply_gradients([(grad, xfvar)])
			cost.append(loss)
			sys.stdout.write(f"\r batch {len(allxfs)}/{len(p3did)}, iter {it}/{niter}, loss {cost[0]:.4f} -> {loss:.4f} ")
			sys.stdout.flush()
			
		allxfs.append(xfvar.numpy())

	allxfs=np.array(allxfs)
	print('\n', allxfs.shape)
	xnp=allxfs.copy()
	
	xnp[:,:3]=xnp[:,:3]*180./np.pi
	xnp[:,3:]*=options.clip
	xfs=[Transform({"type":"eman", "az":x[0], "alt":x[1], 
				"phi":x[2], "tx":x[3], "ty":x[4], "tz":x[5]}) for x in xnp.tolist()]
	aliout=[]
	for ip, pid in enumerate(p3did):
		ali=[alipm[i].copy() for i in pid]
		for a in ali:
			a["xform.projection"]=a["xform.projection"]*xfs[ip].inverse()
			
		aliout.extend(ali)
			
	print("output written to", options.ptclsout)
	save_lst_params(aliout, options.ptclsout)
	
	return
			
			
			
if __name__ == '__main__':
	main()
	
