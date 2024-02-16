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

#### we will import some functions from e2gmm_refine later
emdir=e2getinstalldir()
sys.path.insert(0,os.path.join(emdir,'bin'))

#### need to unify the float type across tenforflow and numpy
floattype=np.float32
from e2gmm_refine_new import *

@tf.function()
def rotpts(px, ang, msk):
	azp=-ang[:,0]*np.pi
	altp=ang[:,1]*np.pi
	phip=-ang[:,2]*np.pi
	trans=ang[:,3:][:,None,:]*.2
	m=msk[None,:,None]
	
	matrix=make_matrix(azp, altp, phip)
	matrix=tf.transpose(matrix)
	matrix=tf.reshape(matrix, shape=[-1, 3,3])
	matrix=tf.transpose(matrix, (0,2,1))

	cnt=tf.reduce_sum(px[0]*msk[:,None], axis=0)/tf.reduce_sum(msk)
	cs=tf.matmul(-cnt[None,:], matrix)+cnt
	t=tf.concat([cs+trans, tf.ones((matrix.shape[0],1,1))], axis=2)
	

	mt4=tf.concat([matrix, tf.zeros((matrix.shape[0], 3,1))], axis=2)
	mt4=tf.concat([mt4, t], axis=1)

	px4=tf.concat([px, tf.ones((px.shape[0], px.shape[1], 1))], axis=2)
	pr=tf.matmul(px4, mt4)[:,:,:3]
	
	return pr

@tf.function()
def rotpts_mult(px, ang, msk):
	m0=tf.ones_like(msk[0])
	p2=tf.zeros_like(px)
	for i,m in enumerate(msk):
		p1=rotpts(px, ang[:,i,:], m)
		m0-=m
		p2+=p1*m[None,:,None]
		
	p2+=px*m0[None,:,None]
	return p2
		

def main():
	
	usage="."
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ptclsin", type=str,help="particles input", default=None)
	parser.add_argument("--ptclsout", type=str,help="particles output", default=None)
	parser.add_argument("--model", type=str,help="model input", default=None)
	parser.add_argument("--mask", type=str,help="masks for refinement. multiple files separated by ','", default=None)
	parser.add_argument("--niter", type=int, help="number of iteration",default=20)
	parser.add_argument("--maxres", type=float, help="resolution",default=5.)
	parser.add_argument("--learnrate", type=float, help="learning rate",default=2e-3)
	parser.add_argument("--l2_weight", type=float, help="weighting factor for L2. default is 1",default=1.)
	parser.add_argument("--xf_file", type=str,help="file for the transform input/output. hdf format. will overwrite.", default=None)
	parser.add_argument("--xf_starti", type=int,help="first line ID for the xf_file. ", default=0)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	
	logid=E2init(sys.argv)
	maxres=options.maxres
	
	
	#### load particles and prepare some parameters
	fname=options.ptclsin
	e=EMData(fname, 0, True)
	raw_apix, raw_boxsz = e["apix_x"], e["ny"]

	maxboxsz=good_size(raw_apix/maxres*raw_boxsz*2)

	options.sym="c1"
	options.maxboxsz=maxboxsz
	options.trainmodel=False
	options.clip=-1

	####   load metadata first
	e=EMData(fname, 0, True)
	raw_apix, raw_boxsz = e["apix_x"], e["ny"]
	data_cpx, xfsnp = load_particles(options)

	print("Image size: ", data_cpx[0].shape)

	##   set up fourier indices for image generation/clipping later
	##   params is a dictionay that saves matrices for Fourier indexing
	apix=raw_apix*raw_boxsz/maxboxsz
	clipid=set_indices_boxsz(data_cpx[0].shape[1], apix, True)
	params=set_indices_boxsz(maxboxsz)
	
	#### load model
	pts=np.loadtxt(options.model).astype(floattype)
	pts=tf.constant(pts[None,:,:])
	
	#### read mask files and re-weight them
	p=pts[0].numpy()
	maskfile=options.mask.split(',')
	imsk=[ make_mask_gmm(m, p) for m in maskfile]
	npatch=len(imsk)
	
	m=np.sum(imsk, axis=0)
	m0=1-m
	m0[m0<0]=0
	m+=m0
	m[m==0]=1.
	for i in range(npatch):
		imsk[i]/=m
		
		
	#### initialize dataset
	options.minpx=2
	options.maxpx=ceil(raw_apix/maxres*raw_boxsz)
	options.batchsz=32
	
	last_xf=np.zeros((len(xfsnp), npatch, 6), dtype=floattype)
	if options.xf_file!=None:
		print(f"Loading previous transform from {options.xf_file}...")
		last_xf_raw=EMData(options.xf_file)
		last_xf_raw=last_xf_raw.numpy().copy()
		print("  reading lines from {} to {}".format(options.xf_starti, options.xf_starti+len(xfsnp)))
		last_xf=last_xf_raw[options.xf_starti:options.xf_starti+len(xfsnp)]
		last_xf=last_xf.reshape((len(xfsnp), npatch, 6)).astype(floattype)
		
	trainset=tf.data.Dataset.from_tensor_slices((data_cpx[0], data_cpx[1], xfsnp, last_xf))
	nsample=data_cpx[0].shape[0]
	trainset=trainset.batch(options.batchsz)
	nbatch=nsample//options.batchsz

	
	#### do one batch first to determine l2
	for ptr,ptj,xf,lx in trainset:
		opt=tf.keras.optimizers.Adam(learning_rate=options.learnrate) 
		ptcl_cpx=(ptr, ptj)
		
		xfvar=tf.Variable(np.zeros((xf.shape[0],npatch, 6), dtype=floattype))
		p0=tf.constant(tf.zeros((xf.shape[0],pts.shape[0], 5))+pts)
		
		cost=[]
		for itr in range(options.niter):
			with tf.GradientTape() as gt:
				p1=rotpts_mult(p0[:,:,:3], xfvar+lx, imsk)
				p2=tf.concat([p1, p0[:,:,3:]], axis=2)
				proj_cpx=pts2img(p2, xf)
				#print(proj_cpx[0].shape, ptcl_cpx[0].shape, params["rings"].shape)
				fval=calc_frc(proj_cpx, ptcl_cpx, params["rings"], minpx=options.minpx, maxpx=options.maxpx)
				loss=-tf.reduce_mean(fval)
				
				l2=tf.reduce_sum(xfvar**2)
				loss2=loss+l2*.0

			grad=gt.gradient(loss2, xfvar)
			opt.apply_gradients([(grad, xfvar)])
			cost.append(loss)
		
		break
	
	l2wd=float((cost[0]-cost[-1])/l2)*options.l2_weight
	print("Refine without L2:")
	print("  loss change {:.3f}, L2 change {:.3f}, use weight decay factor {:.3f}".format(float((cost[0]-cost[-1])), float(l2), l2wd))
		
	
	#### full refinement with l2 loss
	print("Full refinement......")
	xfvs=[]

	for ptr,ptj,xf,lx in trainset:
		opt=tf.keras.optimizers.Adam(learning_rate=options.learnrate) 
		ptcl_cpx=(ptr, ptj)
		xfvar=tf.Variable(np.zeros((xf.shape[0],npatch, 6), dtype=floattype))
		p0=tf.constant(tf.zeros((xf.shape[0],pts.shape[0], 5))+pts)
		
		cost=[]
		for itr in range(options.niter):
			with tf.GradientTape() as gt:
				p1=rotpts_mult(p0[:,:,:3], xfvar+lx, imsk)
				p2=tf.concat([p1, p0[:,:,3:]], axis=2)
				proj_cpx=pts2img(p2, xf)
				#print(proj_cpx[0].shape, ptcl_cpx[0].shape, params["rings"].shape)
				fval=calc_frc(proj_cpx, ptcl_cpx, params["rings"], minpx=options.minpx, maxpx=options.maxpx)
				loss=-tf.reduce_mean(fval)
				
				l2=tf.reduce_sum(xfvar**2)
				loss2=loss+l2*l2wd

			grad=gt.gradient(loss2, xfvar)
			opt.apply_gradients([(grad, xfvar)])
			cost.append(loss)

			sys.stdout.write("\r batch {}/{}, iter {}/{}, loss {:.3f} -> {:.3f}  l2: {:.3f} ".format(len(xfvs), nbatch, itr+1, options.niter, cost[0], loss, l2))
			sys.stdout.flush()

		xfvs.append(xfvar.numpy())

	xfrot=np.vstack(xfvs).copy()
	xfrot+=last_xf
	print()
	
	if options.xf_file!=None:
		print(f"Writing transform to {options.xf_file}...")
		lx=xfrot.reshape((len(xfsnp), -1)).copy()
		last_xf_raw[options.xf_starti:options.xf_starti+len(xfsnp)]=lx
		lx=from_numpy(last_xf_raw)
		lx.write_image(options.xf_file)
	
	angles=tf.data.Dataset.from_tensor_slices((xfrot, xfsnp)).batch(128)
	
	#### now convert rigid body movement to orientation change
	print("Converting back to orientation......")
	p00=pts[0]
	for pid in range(npatch):
		mm=imsk[pid]
		cnt=tf.reduce_sum(p00[:,:3]*mm[:,None], axis=0)/tf.reduce_sum(mm)
		xfnew_all=[]

		for ang0, xf in angles:
			ang=ang0[:,pid,:]
			azp=-ang[:,0]*np.pi
			altp=ang[:,1]*np.pi
			phip=-ang[:,2]*np.pi
			trans=ang[:,3:][:,None,:]*.2
			m=mm[None,:,None]

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
		xnp[:,3:]*=raw_boxsz
		xfs=[Transform({"type":"eman", "az":x[0], "alt":x[1], 
					"phi":x[2], "tx":x[3], "ty":x[4]}) for x in xnp.tolist()]

		lstin=load_lst_params(options.ptclsin)
		for i,xf in enumerate(xfs):
			lstin[i]["xform.projection"]=xf
		
		o=options.ptclsout.split('_')
		o.insert(-1,f'p{pid:02d}')
		oname='_'.join(o)
		if os.path.isfile(oname): os.remove(oname)
		save_lst_params(lstin, oname)
		sys.stdout.write("\r patch {}/{} ".format(pid+1, npatch))
		sys.stdout.flush()
	
	print()
		
	E2end(logid)
	
	
if __name__ == '__main__':
	main()
	
