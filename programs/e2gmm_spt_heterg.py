#!/usr/bin/env python
# Muyuan Chen 2024-02

import numpy as np
from EMAN2 import *

if "CUDA_VISIBLE_DEVICES" not in os.environ:
	# so we can decide which gpu to use with environmental variable
	os.environ["CUDA_VISIBLE_DEVICES"]='0'

#### do not occupy the entire GPU memory at once
##   seems necessary to avoid some errors...
os.environ["TF_FORCE_GPU_ALLOW_GROWTH"]='true'

#### finally initialize tensorflow
import tensorflow as tf
floattype=np.float32

#### we will import some functions from e2gmm_refine later
emdir=e2getinstalldir()
sys.path.insert(0,os.path.join(emdir,'bin'))

from e2gmm_refine_new import *
from e2gmm_spt_align import *
		

def main():
	
	usage="."
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ptclsin", type=str,help="particles input", default=None)
	# parser.add_argument("--ptclsout", type=str,help="particles output", default=None)
	parser.add_argument("--model", type=str,help="model input", default=None)
	parser.add_argument("--mask", type=str,help="masks for refinement. multiple files separated by ','", default=None)
	parser.add_argument("--niter", type=int, help="number of iteration",default=20)
	parser.add_argument("--clip", type=int, help="clip image to size",default=-1)
	
	parser.add_argument("--midin", type=str,help="tranform input", default=None)
	parser.add_argument("--midout", type=str,help="tranform output", default=None)

	parser.add_argument("--anchor", type=str,help="anchor points. will generate from model by default", default=None)
	# parser.add_argument("--xfin_starti", type=int,help="starting index for tranform input", default=0)

	parser.add_argument("--maxres", type=float, help="resolution",default=10.)
	parser.add_argument("--minres", type=float, help="resolution",default=200.)
	parser.add_argument("--learnrate", type=float, help="learning rate",default=1e-5)
	# parser.add_argument("--angle_range", type=str,help="search angle range", default=None)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	
	logid=E2init(sys.argv)
	
	alipm=load_lst_params(options.ptclsin)
	pids=np.array([a["ptcl3d_id"] for a in alipm])
	uid=np.unique(pids)
	p3did=[np.where(pids==u)[0] for u in uid]

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
	print(data_cpx[0].shape, xfsnp.shape)
	print(f"particle size {raw_boxsz}, clip to {options.clip}, shrink to {options.maxboxsz}")
	print(f"compare {options.minpx} to {options.maxpx} Fourier pixels")

	clipid=set_indices_boxsz(data_cpx[0].shape[1], options.apix, True)
	params=set_indices_boxsz(options.maxboxsz)	
	
	#### load model
	pts=np.loadtxt(options.model).astype(floattype)
	print("model shape",pts.shape)
	
	imsk=make_mask_gmm(options.mask, pts)
	print("masking {} out of {} points".format(np.sum(imsk), len(imsk)))
	
	#### anchor points
	if options.anchor==None:
		path=os.path.dirname(options.midout)
		pn=16
		km=KMeans(pn,max_iter=30)
		km.fit(pts[:,:3])
		pc=km.cluster_centers_
		pm=pts[imsk>.1]
		pn=16
		km=KMeans(pn,max_iter=30)
		km.fit(pm[:,:3])
		pc2=km.cluster_centers_

		pcx=np.vstack([pc, pc2])
		pp=np.hstack([pcx, np.zeros((len(pcx),1))+np.mean(pts[:,3]), np.zeros((len(pcx),1))+np.mean(pts[:,4])])
		np.savetxt(f"{path}/model_00_anchor.txt", pp)
		anchor=pp[:,:3].astype(floattype).copy()

	########
	pts=tf.constant(pts[None,...])
	print(pts.shape, anchor.shape)
	decode_model=build_decoder_anchor(pts, anchor, ninp=4)
	encode_model=build_encoder(nout=4, conv=False,ninp=len(pts[0]))

	batchsz=4
	print(len(p3did))
	print(np.array([len(i) for i in p3did]))

	allgrds=[]
	allscr=[]

	for pii in range(0, len(p3did), batchsz):
		pids=p3did[pii:pii+batchsz]
		pid=np.concatenate(pids)
		pn=np.array([len(i) for i in pids])

		ptr=tf.gather(data_cpx[0], pid)
		ptj=tf.gather(data_cpx[1], pid)
		xf=xfsnp[pid]

		ptcl_cpx=(ptr, ptj)

		with tf.GradientTape() as gt:
			pt=tf.Variable(tf.repeat(pts, len(pn), axis=0))
			pt2=tf.repeat(pt, pn, axis=0)

			imgs_cpx=pts2img(pt2, xf)
			fval=calc_frc(ptcl_cpx, imgs_cpx, params["rings"], minpx=options.minpx, maxpx=options.maxpx)

			loss=-tf.reduce_mean(fval)

		grad=tf.convert_to_tensor(gt.gradient(loss, pt))
		scr=tf.reshape(fval, (len(pids), -1))

		scr=tf.reduce_mean(scr, axis=1)

		allgrds.append(grad)
		allscr.append(scr)
		print(pii, len(p3did), np.mean(scr), end='\r')

	allgrds=np.concatenate(allgrds, axis=0)
	allscr=np.concatenate(allscr, axis=0)
	allgrds=allgrds/np.std(allgrds)
	print(allgrds.shape)
	np.savetxt(f"{path}/allgrds_00.txt", allgrds.reshape((len(allgrds), -1)))

	return


	#### refinement iterations
	for iip, pid in enumerate(p3did):
		ptr=tf.gather(data_cpx[0], pid)
		ptj=tf.gather(data_cpx[1], pid)
		ptcl_cpx=(ptr, ptj)
		xf=xfsnp[pid]
		
		ang_xfv=[]
		ang_loss=[]
		for ia, starta in enumerate(angrng):
			opt=tf.keras.optimizers.Adam(learning_rate=options.learnrate) 
			# xv=np.zeros((1,6), dtype=floattype)
			xv=xfin[iip][None, :].copy()
			# print(xv, xv.shape)
			xv[0,0]+=starta
			xfvar=tf.Variable(xv)
			
			for div in res_rng:
				cost=[]
				for it in range(options.niter):
					with tf.GradientTape() as gt:
						
						p1=rotpts_mult(pts[None, :,:3], xfvar[:,None,:], [imsk])
						p1=tf.concat((p1, pts[None,:,3:]), axis=2)
						proj_cpx=pts2img(p1, xf)
						
						fval=calc_frc(proj_cpx, ptcl_cpx, params["rings"], minpx=options.minpx, maxpx=options.maxpx//div)
						loss=-tf.reduce_mean(fval)

					if it>5 and loss>cost[-1]: break
					grad=gt.gradient(loss, xfvar)
					opt.apply_gradients([(grad, xfvar)])
					cost.append(loss)
					
				sys.stdout.write(f"\r batch {len(allxfs)}/{len(p3did)}, angle {ia}/{len(angrng)}, {xv[0,0]:.3f} -> {xfvar[0,0]:.3f} loss {cost[0]:.4f} -> {loss:.4f} ")
				sys.stdout.flush()
			
			ang_xfv.append(xfvar.numpy().copy())
			ang_loss.append(loss)
			
		xfv=ang_xfv[np.argmin(ang_loss)]
		allxfs.append(xfv)

	xfrot=np.vstack(allxfs).copy()
	print(xfrot.shape)
	xfrot=np.hstack([np.arange(len(xfrot))[:,None], xfrot])
	np.savetxt(options.xfout, xfrot)

	xfrot2=np.zeros((len(xfsnp), 6), dtype=floattype)
	for i,ip in enumerate(p3did):
		xfrot2[ip]=xfrot[i,1:]
		
	print(np.std(xfrot2))
	
	
	
	angles=tf.data.Dataset.from_tensor_slices((xfrot2, xfsnp)).batch(128)

	#### now convert rigid body movement to orientation change
	print("Converting back to orientation......")
	p00=pts
	cnt=tf.reduce_sum(p00[:,:3]*imsk[:,None], axis=0)/tf.reduce_sum(imsk)
	xfnew_all=[]

	for ang, xf in angles:
		azp=-ang[:,0]*np.pi
		altp=ang[:,1]*np.pi
		phip=-ang[:,2]*np.pi
		trans=ang[:,3:][:,None,:]*.2
		m=imsk[None,:,None]

		matrix=make_matrix(azp, altp, phip)
		matrix=tf.transpose(matrix)
		matrix=tf.reshape(matrix, shape=[-1, 3,3])
		matrix=tf.transpose(matrix, (0,2,1))

		cs=tf.matmul(-cnt[None,:], matrix)+cnt
		matrix4_rot=tf.concat([matrix, tf.zeros((matrix.shape[0], 3,1))], axis=2)
		t=tf.concat([cs+trans, tf.ones((matrix.shape[0],1,1))], axis=2)
		matrix4_rot=tf.concat([matrix4_rot, t], axis=1)

		azp=-xf[:,0]
		altp=xf[:,1]
		phip=-xf[:,2]
		trans=xf[:,3:][:,None,:]*[1,-1]

		matrix=make_matrix(azp, altp, phip)
		matrix=tf.transpose(matrix)
		matrix=tf.reshape(matrix, shape=[-1, 3,3])
		matrix=tf.transpose(matrix, (0,2,1))

		matrix4_proj=tf.concat([matrix, tf.zeros((matrix.shape[0], 3,1))], axis=2)
		t=tf.concat([trans, tf.zeros((matrix.shape[0],1,1)), tf.ones((matrix.shape[0],1,1))], axis=2)
		matrix4_proj=tf.concat([matrix4_proj, t], axis=1)

		matrix4_full=tf.matmul(matrix4_rot, matrix4_proj)

		t=tf.transpose(matrix4_full[:,:3,:3], (0,2,1))
		cos_altp=t[:,2,2]
		sin_altp=np.sqrt(1-cos_altp**2)
		cos_azp=t[:,2,1]/(-sin_altp)
		sin_azp=t[:,2,0]/sin_altp
		cos_phip=t[:,1,2]/sin_altp
		sin_phip=t[:,0,2]/sin_altp

		altp=tf.math.atan2(sin_altp, cos_altp)
		azp=tf.math.atan2(sin_azp, cos_azp)
		phip=tf.math.atan2(sin_phip, cos_phip)

		ts=matrix4_full[:,3]
		xfnew=tf.stack([-azp, altp, -phip, ts[:,0], -ts[:,1]])
		xfnew=tf.transpose(xfnew)

		xfnew_all.append(xfnew)


	xfnew_all=tf.concat(xfnew_all, axis=0)


	xnp=xfnew_all.numpy().copy()
	xnp[:,:3]=xnp[:,:3]*180./np.pi
	xnp[:,3:]*=options.clip
	xfs=[Transform({"type":"eman", "az":x[0], "alt":x[1], 
				"phi":x[2], "tx":x[3], "ty":x[4]}) for x in xnp.tolist()]

	lstin=load_lst_params(options.ptclsin)
	for i,xf in enumerate(xfs):
		lstin[i]["xform.projection"]=xf


	oname=options.ptclsout
	if os.path.isfile(oname): os.remove(oname)
	save_lst_params(lstin, oname)

	E2end(logid)
	
	
if __name__ == '__main__':
	main()
	
