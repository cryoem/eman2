#!/usr/bin/env python
# Muyuan Chen 2025-12
## jax reimplementation of gmm pipeline
## with new stuff..
from EMAN2 import *
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import scipy.spatial.distance as scipydist
from flax.serialization import to_state_dict, from_state_dict
import pickle
from EMAN2_utils import pdb2numpy

if "CUDA_VISIBLE_DEVICES" not in os.environ:
	# so we can decide which gpu to use with environmental variable
	os.environ["CUDA_VISIBLE_DEVICES"]='0' 
import jax
import jax.numpy as jnp
	
from flax import nnx
import flax.linen as nn
import optax
from jax import jit


def make_matrix_3d(angles):
	"""Generate a rotation matrix from Euler angles."""
	alpha, beta, gamma = angles[:3]

	R_x = jnp.array([[1, 0, 0],
					 [0, jnp.cos(alpha), -jnp.sin(alpha)],
					 [0, jnp.sin(alpha), jnp.cos(alpha)]])
	
	R_y = jnp.array([[jnp.cos(beta), 0, jnp.sin(beta)],
					 [0, 1, 0],
					 [-jnp.sin(beta), 0, jnp.cos(beta)]])
	
	R_z = jnp.array([[jnp.cos(gamma), -jnp.sin(gamma), 0],
					 [jnp.sin(gamma), jnp.cos(gamma), 0],
					 [0, 0, 1]])
	
	R = R_x @ R_y @ R_z  # Combined rotation
	return R


#### Symmetrize the Gaussian coordinates. Only works for c/d sym right now
##   note this only keep the (x,y,z) columns even with c1 symmetry
def get_sym_pts(sym, pts):
	if sym=="c1":
		return [tf.transpose(pts[:,:3])]
	nsym=Transform.get_nsym(sym)
	asym=[]
	if sym.startswith('d'):
		for i in range(nsym//2):
			a=np.pi*4/nsym*i
			asym.append(
				[pts[:,0]*np.cos(a)-pts[:,1]*np.sin(a),
				 pts[:,0]*np.sin(a)+pts[:,1]*np.cos(a),
				 pts[:,2]])
			asym.append(
				[pts[:,0]*np.cos(a)-pts[:,1]*np.sin(a),
				 -(pts[:,0]*np.sin(a)+pts[:,1]*np.cos(a)),
				 -pts[:,2]])
			

	if sym.startswith('c'):
		asym=[]
		for i in range(nsym):
			a=np.pi*2/nsym*i
			asym.append(
				[pts[:,0]*np.cos(a)-pts[:,1]*np.sin(a),
				 pts[:,0]*np.sin(a)+pts[:,1]*np.cos(a),
				 pts[:,2]])
			
	return asym


""" Project 3d Gaussian coordinates based on transforms to make projection
	input:  pts - ( number of Gaussian, 3 (x,y,z) )
			ang - ( 5 (az, alt, phi, tx, ty) )
"""
def xf2pts(pts, ang):
	azp=-ang[0]
	altp=ang[1]
	phip=-ang[2]

	matrix=jnp.stack([(jnp.cos(phip)*jnp.cos(azp) - jnp.cos(altp)*jnp.sin(azp)*jnp.sin(phip)),
	(jnp.cos(phip)*jnp.sin(azp) + jnp.cos(altp)*jnp.cos(azp)*jnp.sin(phip)),
	(jnp.sin(altp)*jnp.sin(phip)),

	(-jnp.sin(phip)*jnp.cos(azp) - jnp.cos(altp)*jnp.sin(azp)*jnp.cos(phip)),
	(-jnp.sin(phip)*jnp.sin(azp) + jnp.cos(altp)*jnp.cos(azp)*jnp.cos(phip)),
	(jnp.sin(altp)*jnp.cos(phip)),

	(jnp.sin(altp)*jnp.sin(azp)),
	(-jnp.sin(altp)*jnp.cos(azp)),
	jnp.cos(altp)], 0)

	#### rotate Gaussian positions
	matrix=jnp.reshape(matrix, shape=[3,3]) 
	matrix=jnp.transpose(matrix)

	pts_rot=jnp.matmul(pts[:,:3], matrix)

	#### finally do the translation
	pts_rot_trans=jnp.stack([pts_rot[:,0]+ang[3], 
							-pts_rot[:,1]+ang[4], 
							pts_rot[:,2]-ang[5], 
							pts[:,3], pts[:,4]
							], 1)
	#pts_rot_trans=jnp.stack([(pts_rot[:,0]+ang[3]), (-pts_rot[:,1])+ang[4]], 1)

	return pts_rot_trans

xf2pts_multi=jax.vmap(xf2pts, in_axes=(0,0))

def pts2img_one_sig(pts, ang):
	
	sz, idxft, rrft = params["sz"], params["idxft"], params["rrft"]
	xfo=params["xforigin"]
	
	bamp=jax.nn.relu(pts[:,3]).astype(np.complex64)[:,None]
	bsigma=jax.nn.relu(pts[:,4]).astype(np.complex64)*.001

	bpos=xf2pts(pts[:,:3],ang)
	bpos=bpos*params["sz_scale"]+sz/2
	bposft=bpos*np.pi*2

	cpxang_x=bposft[:,0,None] * idxft[0,[0],:]
	cpxang_y=bposft[:,1,None] * idxft[1,:,[0]]

	sig_x = jnp.exp(-rrft[0][None, :]*bsigma[:,None])
	sig_y = jnp.exp(-rrft[1][None, :]*bsigma[:,None])

	pgauss_x = jnp.exp(-1j* cpxang_x.astype(jnp.complex64))*bamp*sig_x
	pgauss_y = jnp.exp(-1j* cpxang_y.astype(jnp.complex64))*bamp*sig_y


	pgauss = jnp.matmul(jnp.transpose(pgauss_x), pgauss_y)
	pgauss = jnp.transpose(pgauss)
	return pgauss*xfo.astype(np.complex64)

pts2img=jax.vmap(pts2img_one_sig, in_axes=(0,0))

#### compute particle-projection FRC 
##   data_cpx, imgs_cpx - complex images in (real, imag) form
##   rings - indices of Fourier rings
##   return_curve - return the curve instead of average FRC score
##   minpx - skip the X initial low freq pixels
def calc_frc(data_cpx, imgs_cpx, rings, return_curve=False, minpx=1, maxpx=-1):
	mreal, mimag=imgs_cpx.real, imgs_cpx.imag
	dreal, dimag=data_cpx.real, data_cpx.imag
	#### normalization per ring
	nrm_img=mreal**2+mimag**2
	nrm_data=dreal**2+dimag**2

	nrm0=jnp.tensordot(nrm_img, rings, [[1,2],[0,1]])
	nrm1=jnp.tensordot(nrm_data, rings, [[1,2],[0,1]])

	nrm=jnp.sqrt(nrm0)*jnp.sqrt(nrm1)
	nrm=jnp.maximum(nrm, 1e-4) #### so we do not divide by 0

	#### average FRC per batch
	ccc=mreal*dreal+mimag*dimag
	frc=jnp.tensordot(ccc, rings, [[1,2],[0,1]])/nrm

	if return_curve:
		return frc
	else:
		frcval=jnp.mean(frc[:, minpx:maxpx], axis=1)
		return frcval

@jit
def get_slice(volume, ang):
	#### make projections of 3D volume
	
    R=make_matrix_3d(ang)
    rotated_coords = R @ ind_xy.T
    rotated_coords+=rawbox//2
    rotated_coords=jnp.concatenate([rotated_coords[:,:2], abs(rotated_coords[:,2:])], axis=1)
    
    ix=jnp.floor(rotated_coords).astype(int)
    ir=rotated_coords-jnp.floor(rotated_coords)
    ir=[ir, 1-ir]
    
    pj=jnp.zeros((rawbox, rawbox), dtype=np.complex64).flatten()
    
    for i0 in ind_interp[...,None]:
        w=ir[0]*i0+ir[1]*(1-i0)
        w=jnp.prod(w, axis=0)
        i=ix+i0
        pj=pj+volume[i[0], i[1], i[2]]*w

    pj=pj.reshape((rawbox, rawbox))
    return pj

get_slice_multi=jax.vmap(get_slice, in_axes=(None, 0))



#### load particles from file and fourier transform them
#### particles need to have their transform in file header or comment of list file
##   will also shrink particles and do FT
##   return Fourier images and transform matrix
def load_particles(options):
	fname=options.ptclsin
	boxsz=options.maxboxsz
	
	nptcl=EMUtil.get_image_count(fname)
	if options.id_keep: nptcl=len(options.id_keep)
	projs=[]
	hdrs=[]
	e=EMData(fname, 0, True)
	rawbox=e["nx"]
	if options.clip>0 and rawbox!=options.clip:
		cl=f" Clip to {options.clip}."
	else:
		cl=""
		
	print("Loading {} particles of box size {}.{} Shrink to {}".format(nptcl, rawbox, cl, boxsz))
	for j in range(0,nptcl,1000):
		print(f"\r {j}/{nptcl} ",end="")
		sys.stdout.flush()
		rg=range(j,min(j+1000,nptcl))
		if options.id_keep: rg=options.id_keep[j:j+1000]
		el=EMData.read_images(fname,rg)
		print(f"R     ",end="")
		sys.stdout.flush()
		for e in el:
			if options.clip>0 and rawbox!=options.clip:
				e.clip_inplace(Region((rawbox-options.clip)//2,(rawbox-options.clip)//2, options.clip,options.clip))
			if e["nx"]!=boxsz:
				e.process_inplace("math.fft.resample",{"n":e["nx"]/boxsz})
			
			hdrs.append(e.get_attr_dict())
			projs.append(e.numpy().copy())
			
	print(f"{nptcl}/{nptcl}")
	projs=np.array(projs, dtype=np.float32)
	
	data_cpx=np.fft.rfft2(np.fft.fftshift(projs*1e-3, axes=(1,2))).astype(np.complex64)

	if fname.endswith(".lst"):
		info=load_lst_params(fname, options.id_keep)
		xfs=[p["xform.projection"] for p in info]
			
	elif "xform.projection" in hdrs[0]:
		xfs=[p["xform.projection"] for p in hdrs]
		info=[{"src":fname, "idx":i, "xform.projection":x} for i,x in enumerate(xfs)]
		
	else:
		print("Error. No existing transform from particles...")
		exit()
		
	xfpm=[x.get_params("eman") for x in xfs]
	xfsnp=np.array([[x["az"],x["alt"],x["phi"], x["tx"], x["ty"], 0] for x in xfpm])
	xfsnp[:,:3]=xfsnp[:,:3]*np.pi/180.
	xfsnp[:,3:]/=float(options.clip)
	
	print("Data read complete")
	print("Data shape {}, xform shape {}".format(data_cpx.shape, xfsnp.shape))
	return data_cpx, xfsnp, info
	

#### compute fourier indices for image generation, clipping, and frc
#### pass indices in a dictionary
def set_indices_boxsz(boxsz, options, apix=0, return_freq=False):
	idx=np.indices((boxsz,boxsz))-boxsz//2
	idx=np.fft.ifftshift(idx)
	idx=idx[:,:,:boxsz//2+1]
	
	sz=sz_scale=boxsz
	sz_scale=boxsz*options.rawbox/options.clip
		
	idxft=(idx/sz).astype(np.float32)
	rrft=np.sqrt(np.sum(idx**2, axis=0)).astype(np.float32)## batch, npts, x-y
	rrft_x = idx[0,0,:]**2
	rrft_y = idx[1,:,0]**2
	rrft=[rrft_x, rrft_y]
	

	rr=np.round(np.sqrt(np.sum(idx**2, axis=0))).astype(int)
	rings=np.zeros((sz,sz//2+1,sz//2), dtype=np.float32) #### Fourier rings
	for i in range(sz//2):
		rings[:,:,i]=(rr==i)
		
	xvec=np.fromfunction(lambda i,j:1.0-2*((i+j)%2),(sz,sz//2+1),dtype=np.float32)
	global params
	params={"sz":sz, "idxft":idxft, "rrft":rrft, "rings":rings, "xforigin":xvec, "sz_scale":sz_scale}
	return params


#### training decoder on projections
def train_decoder(gen_model, trainset, params, options, pts=None):
	return

def save_ptcls_xform(xfsnp, boxsz, options, scr):
	xnp=xfsnp.copy()
	xnp[:,:3]=xnp[:,:3]*180./np.pi
	xnp[:,3:]*=boxsz
	xfs=[Transform({"type":"eman", "az":x[0], "alt":x[1], 
				"phi":x[2], "tx":x[3], "ty":x[4]}) for x in xnp.tolist()]

	oname=options.ptclsout
	print("saving aligned particles to {}".format(oname))
	lstin=load_lst_params(options.ptclsin)
	for i,xf in enumerate(xfs):
		lstin[i]["xform.projection"]=xf
		lstin[i]["score"]=-scr[i]
		
	if os.path.isfile(oname): os.remove(oname)
	save_lst_params(lstin, options.ptclsout)
	
	
def calc_gradient(data_cpx, data_xf, pts, params, options):
	print("Calculating gradient input...")
	allgrds=[]
	minpx,maxpx=options.minpx, options.maxpx_grad
	bsz=options.batchsz
	@jit
	def pt_to_frc(pt0, xf, pj_cpx):
		if options.spt:
			pt=jnp.repeat(pt0, xf.shape[1], axis=0)
			xf=xf.reshape((xf.shape[0]*xf.shape[1], xf.shape[2]))
			pj_cpx=pj_cpx.reshape((pj_cpx.shape[0]*pj_cpx.shape[1], pj_cpx.shape[2], pj_cpx.shape[3]))
		else:
			pt=pt0
			
		imgs_cpx=pts2img(pt, xf)
		fval=calc_frc(pj_cpx, imgs_cpx, params["rings"], minpx=minpx, maxpx=maxpx)
		return -jnp.mean(fval)
		
	grad_fn_00=jax.value_and_grad(pt_to_frc)
	
	if options.spt:
		nn=len(options.idx3d)
	else:
		nn=len(data_cpx)
	
	for ib in range(0, nn, bsz):
		if options.spt:
			i3=options.idx3d[ib:ib+bsz]
			pj_cpx=data_cpx[i3]
			xf=data_xf[i3]
			pt=np.repeat(pts[None,...], len(i3), axis=0)
			_, grd = grad_fn_00(pt, xf, pj_cpx)
		else:
			pj_cpx=data_cpx[ib:ib+bsz]
			xf=data_xf[ib:ib+bsz]
			pt=np.repeat(pts[None,...], xf.shape[0], axis=0)
			_, grd = grad_fn_00(pt, xf, pj_cpx)
		
		allgrds.append(grd)
	
	allgrds=np.concatenate(allgrds, axis=0)
	allgrds=allgrds/np.std(allgrds)
	allgrds=allgrds[:,:,:4]/np.std(allgrds, axis=(0,1))[:4]
	
	print("Gradient input shape:", allgrds.shape)
	return allgrds
	
	
def calc_conf(encode_model, enc_var, allgrds, rng, bsz=1000):
	
	## conformation output
	mid=[]
	for i in range(0, len(allgrds), bsz):
		a=allgrds[i:i+bsz]
		m = encode_model.apply(enc_var, a, training=False, rngs={"dropout": rng})
		mid.append(m)
		
	mid=np.concatenate(mid, axis=0)
	print("Encoder output shape:", mid.shape)
	return mid
	
def build_encoder(pts, allgrds, rng, options):
	
	print("Building Encoder...")
	nmid=options.nmid
	nhidden=options.ndense
	npt=pts.shape[0]

	class Encoder(nn.Module):
		@nn.compact
		def __call__(self, x, training=False):
			x = jnp.reshape(x, (x.shape[0],-1))
			for i in range(3):
				x = nn.Dense(nhidden)(x)
				x = nn.relu(x)
				#x = nn.Dropout(0.1, deterministic=not training)(x)
				
			x = nn.Dense(nmid, kernel_init=nn.initializers.normal(stddev=1e-3))(x)
			return x

	encode_model = Encoder()
	ag=allgrds[:options.batchsz]
	enc_var = encode_model.init(rng, ag)
	mid = encode_model.apply(enc_var, ag, training=True, rngs={"dropout": rng})
	
	print("Output shape:", mid.shape)	
	return encode_model, enc_var


def build_decoder(pts,  rng, options):
	
	print("Building Decoder...")
	nmid=options.nmid
	nhidden=options.ndense
	nmask=pts.shape[0]
	nfac=pts.shape[1]
	if options.fullrot>=0: 
		nfac+=1
		xfmask=np.zeros((1,nmask,6))
		xfmask[:,:,options.fullrot]=1
		
	class Decoder(nn.Module):
		@nn.compact
		def __call__(self, x, training=False):
			for i in range(3):
				x = nn.Dense(nhidden)(x)
				x = nn.relu(x)
				#x = nn.Dropout(0.1, deterministic=not training)(x)
				
			x = nn.Dense(nmask*nfac, kernel_init=nn.initializers.normal(stddev=1e-3))(x)
			# x = nn.Dense(nmask*nfac)(x)
			x = jnp.reshape(x, (x.shape[0],nmask, nfac))
			if options.fullrot>=0: 
				xa=jnp.arctan2(x[:,:,-1], x[:,:,0])
				x=x[:,:,:6]*(1-xfmask)
				x+=xa[:,:,None]*xfmask      
				
			return x

	mid=jnp.zeros((options.batchsz, nmid))
	decode_model = Decoder()
	dec_var = decode_model.init(rng, mid)
	out = decode_model.apply(dec_var, mid, training=True, rngs={"dropout": rng})

	print("Output shape:", out.shape)
	print("Deviation from neutral model: ", np.mean(abs(out)))	
	return decode_model, dec_var


def build_decoder_point_transformer(pts, rng, options):
	
	print("Building Point Transformer Decoder...")
	nmid=options.nmid
	nhidden=options.ndense
	npt=pts.shape[0]
	
	################## 
	pt0=pts[:,:3].copy()
	pns=[ 1024, 256, 64]
	nlayer=len(pns)
	ptpool=[pt0.copy()]
	ptidx=[np.arange(len(pt0), dtype=int)]

	for i in range(nlayer):
		pn=pns[i]
		km=KMeans(pn,max_iter=30, random_state=123)
		km.fit(ptpool[-1])
		pt1=km.cluster_centers_
		
		d=scipydist.cdist(pt1,pt0)
		d=np.argmin(d, axis=1)
		pt1=pt0[d]
		ptidx.append(d)

		print(pt1.shape, np.mean(pt1, axis=0))
		ptpool.append(pt1.copy())
		
	nlayer=len(ptpool)-1
	ptp=np.vstack(ptpool)
	
	################## 
	nnb=options.pt_neighbor
	nup=nnb//2
	nbs=[]
	pls=[]
	ups=[]
	upmat=[]

	for i,pt in enumerate(ptpool):
		
		dst01=scipydist.cdist(pt,pt)
		dstid0=np.argsort(dst01, axis=0)[:nnb].T
		nbs.append(dstid0)
		print(dstid0.shape)
		
		if i<len(ptpool)-1:
			dst01=scipydist.cdist(pt,ptpool[i+1])
			pool0=np.argsort(dst01, axis=0)[:nnb].T
			pls.append(pool0)
			# print(pool0.shape)
			
			pool1=np.argsort(dst01, axis=1)[:,:nup]
			
			ups.append(pool1)
			dwt=np.sort(dst01, axis=1)[:,:nup]
			
			dwt=np.exp(-(dwt**2)*50).astype(np.float32)
			# dwt=dwt/np.linalg.norm(dwt, axis=1)[:,None]
			dwt=dwt/np.sum(dwt, axis=1)[:,None]
			upmat.append(dwt.astype(np.float32))			
	
	
	################## 
	upmat1=[]
	ups1=[]
	for ip in range(2, len(ptpool)):
		dst01=scipydist.cdist(ptpool[0],ptpool[ip])
		pool1=np.argsort(dst01, axis=1)[:,:nup]
		
		ups1.append(pool1.copy())
		dwt=np.sort(dst01, axis=1)[:,:nup]
		dwt=np.exp(-(dwt**2)*50).astype(np.float32)
		dwt=dwt/np.sum(dwt, axis=1)[:,None]
		upmat1.append(dwt.astype(np.float32).copy())
		
		print(ip, len(ptpool[ip]), ups1[-1].shape, upmat1[-1].shape)
		
	################## 
	verbose=0

	class Decoder(nn.Module):
		@nn.compact
		def __call__(self, x, training=False):
			if verbose: print("###############")
			
			x1=x
			x2=jnp.zeros((len(x), pts.shape[0], pts.shape[1]), dtype=np.float32) + pts
			if verbose: print("x1,x2:",x1.shape, x2.shape)
			
			nout=256
			nk=512
			dout=5
			
			midshape=(ptpool[-1].shape[0],nout)
			x1 = nn.Dense(128)(x1)
			x1 = nn.relu(x1)
			x1 = nn.Dropout(0.1, deterministic=not training)(x1)
			## missing BatchNormalization
			x1 = nn.Dense(np.prod(midshape))(x1)
			x1=x1.reshape((-1, midshape[0], midshape[1]))

			y0=x1
			if verbose: print("x1: ", x1.shape)

			y3s=[]
			for li in range(1,len(nbs)):
				cur=-li
				
				if verbose: print("####", li, cur)
		
				y_q = nn.Dense(nk)(y0)
				y_k = nn.Dense(nk)(y0)
				y_k = y_k[:,nbs[cur]]

				y_qk=y_q[:,:,None,:]-y_k
				if verbose: print("y_q,k,qk:", y_q.shape, y_k.shape, y_qk.shape)

				p_r=x2[:, ptidx[cur]][:,:,:3]
				p_r=p_r[:, nbs[cur]]
				dp_r=ptpool[cur][:,None,:]-p_r
				if verbose: print("x2,p_r,dp_r: ", x2.shape, p_r.shape, dp_r.shape)
				
				dp_r = nn.Dense(3)(dp_r)
				## missing BatchNormalization
				dp_r = nn.relu(dp_r)
				dp_r = nn.Dense(nk)(dp_r)
				
				yqk_dp=y_qk+dp_r
				if verbose: print("dp_r,yqk_dp: ", dp_r.shape, yqk_dp.shape)

				w=yqk_dp
				## missing BatchNormalization
				w = nn.relu(w)
				w = nn.Dense(nk)(w)
				## missing BatchNormalization
				w = nn.relu(w) 
				w = nn.Dense(nk)(w)
				w = nn.softmax(w, axis=2)
				if verbose: print("w: ", w.shape)

				y_ak = nn.Dense(nk)(y_k) + dp_r
				if verbose: print("y_ak: ", y_ak.shape)

				y1 = jnp.sum(w*y_ak, axis=2)
				## missing BatchNormalization
				y1 = nn.Dropout(0.2, deterministic=not training)(y1)
				if verbose: print("y0, y1: ", y0.shape, y1.shape)
				
				y1 = jnp.concat([y0,y1], axis=2)

				yk=y1
				yk = nn.Dense(nk)(yk)
				yk = nn.Dropout(0.2, deterministic=not training)(yk)
				## missing BatchNormalization
				yk = nn.Dense(nout)(yk)
				if verbose: print("y0, yk: ", y0.shape, yk.shape)
					
				y0k=y0+yk
				
				y0=y0k[:, ups[cur]]
				y0=y0*upmat[cur][:,:,None]
				y0=jnp.sum(y0, axis=2)
				if verbose: print("y0: ", y0.shape)

				if li<len(nbs)-1:
					y3=y0k[:, ups1[cur]]
					y3=y3*upmat1[cur][:,:,None]
					y3=jnp.sum(y3, axis=2)
					y3s.append(y3)
					if verbose: print("y3: ", y3.shape)

			y0=jnp.concat([y0]+y3s, axis=2)
			if verbose: print("y0: ", y0.shape)
			
			y2 = nn.Dense(dout, use_bias=False, kernel_init=nn.initializers.normal(stddev=1e-4))(y0)
			y2=y2*0.5
			if verbose: print("y2: ", y2.shape)
			# y2=y0
			return y2

	mid=jnp.zeros((options.batchsz, nmid))
	decode_model = Decoder()
	dec_var = decode_model.init(rng, mid)
	out = decode_model.apply(dec_var, mid, training=True, rngs={"dropout": rng})

	print("Output shape:", out.shape)
	print("Deviation from neutral model: ", np.mean(abs(out)))	
	return decode_model, dec_var
	
	
#### train the conformation manifold from particles
def train_heterg(data_cpx, data_xf, allgrds, variables, encode_model, decode_model, params, pts, rng, options):
	print("Compiling training function...")
	bsz=options.batchsz
	pas=[int(i) for i in options.pas]
	pas=np.array(pas[:1]*3+pas[1:], dtype=np.float32)
	
	#######################
	def calc_loss_raw(variables, grd, dcpx, xf, rng, return_all=False):
		
		enc_var, dec_var=variables
		
		conf = encode_model.apply(enc_var, grd, training=True, rngs={"dropout": rng})
		cl=jnp.sqrt(jnp.sum(conf**2, axis=1))
		cl=jnp.mean(jnp.maximum(cl-1,0))

		conf=0.01*jax.random.normal(rng, conf.shape)+conf
		
		pout = decode_model.apply(dec_var, conf, training=True, rngs={"dropout": rng})
		p0=jnp.zeros((xf.shape[0], len(pts), 5))+pts

		pout=pout*pas
		pout=pout*options.mask[None,:,None]
		pout+=p0
				
		imgs_cpx=pts2img(pout, xf)
		fval=calc_frc(dcpx, imgs_cpx, params["rings"], minpx=options.minpx, maxpx=options.maxpx)
				
		loss=-jnp.mean(fval)+cl*.1
			
		return loss

	ib=0
	grd=allgrds[ib:ib+bsz]
	dcpx=data_cpx[ib:ib+bsz]
	xf=data_xf[ib:ib+bsz]
	
	calc_loss=jit(calc_loss_raw)
	grad_fn=jax.value_and_grad(calc_loss)
	loss, grads=grad_fn(variables, grd, dcpx, xf, rng)
	print("Initial loss:", loss)	
	
	#######################
	print("Now training...")
	allcost=[]
	optimizer = optax.adam(options.learnrate)
	opt_state = optimizer.init(variables)
	for itr in range(options.niter):
		cost=[]
		
		for ib in range(0, len(data_cpx), bsz):

			grd=allgrds[ib:ib+bsz]
			dcpx=data_cpx[ib:ib+bsz]
			xf=data_xf[ib:ib+bsz]
			
			rng=jax.random.fold_in(rng,123)    
			
			loss, grads=grad_fn(variables, grd, dcpx, xf, rng)
			updates, opt_state = optimizer.update(grads, opt_state)
			variables = optax.apply_updates(variables, updates)
			cost.append(loss)            

			if (ib//bsz)%10==0: 
				sys.stdout.write("\r {}/{}\t{:.4f}      ".format(ib, len(data_cpx), loss))
				sys.stdout.flush()

		sys.stdout.write("\r")
		ac=np.mean(cost)
		allcost.append(ac)
		print("iter {},  loss : {:.4f}".format(itr, ac))    
	
	return variables

def align_particles_00(data_cpx, data_xf, params, pts, lst_raw, options):
	#### simple optimization of projection angle
	#def calc_loss_3d_raw(xf_rot, xf_proj, dcpx, return_all=False):
		
		#xr=xf_rot[:,0]*np.array([1,1,1,1,1,0])
		
		#ptsx=jnp.repeat(pts[None,...], dcpx.shape[0], axis=0)
		#pts1=xf2pts_multi(ptsx, xr)*np.array([1,-1,1,1,1])
		
		#imgs_cpx=pts2img(pts1, xf_proj)
		#fval=calc_frc(dcpx, imgs_cpx, params["rings"], minpx=options.minpx, maxpx=options.maxpx)
		#if return_all: return fval
		#loss=-jnp.mean(fval)
		#return loss
		
	#### rigid body version
	def calc_loss_3d_raw(xf_rot, xf_proj, dcpx, l2wd=0, return_all=False):
		if options.focus:
			xf_rot=xf_rot*np.array([0,1])[None,:,None]
		xr=xf_rot.reshape((-1,6))
		ptsx=jnp.repeat(pts[None,...], dcpx.shape[0]*nmask, axis=0)
		
		pts1=xf2pts_multi(ptsx, xr)*np.array([1,-1,1,1,1])
		pts1=pts1.reshape((dcpx.shape[0], nmask, len(pts), 5))
		pts1=pts1*imasks[None,:,:,None]
		pts1=jnp.sum(pts1, axis=1)
		
		imgs_cpx=pts2img(pts1, xf_proj)
		fval=calc_frc(dcpx, imgs_cpx, params["rings"], minpx=options.minpx, maxpx=options.maxpx)
		
		if return_all: return fval
	
		l2=jnp.sum((xf_rot-xf_rot00)**2)
		loss=-jnp.mean(fval)+l2*l2wd
		return loss
	
	#### dealing with masks
	if options.niter<0: options.niter=50
	if options.learnrate<0: options.learnrate=1e-3
	
	options.focus=False
	if options.mask==None:
		nmask=1
		imasks=jnp.ones((nmask, pts.shape[0]))
		
	else:
		masks=EMData.read_images(options.mask)		
		imasks=np.array([make_mask_gmm(m, pts) for m in masks])
		
		if len(imasks)==1:
			nmask=2
			print(f"Focused refinement with mask size", int(np.sum(imasks[0])))
			imx=imasks[0]>.5
			imasks=np.zeros((2, len(pts)), dtype=float)
			imasks[0]=1-imx
			imasks[1]=imx
			options.focus=True
			
		else:
			nmask=len(imasks)
			imx=np.argmax(imasks, axis=0)
			imasks*=0
			imasks[imx, np.arange(imasks.shape[1])]=1.
		
			print(f"Using {nmask} masks. Sizes:", *np.sum(imasks, axis=1).astype(int))
		
	calc_loss_3d=jit(calc_loss_3d_raw)
	grad_fn_3d=jax.value_and_grad(calc_loss_3d)
	bsz=options.batchsz
	
	readxf=False
	if options.mask:
		n=nmask
		if options.focus: n-=1
		
		try:
			xfs=[Transform(lst_raw[0][f"xform.projection_{i:02d}"]) for i in range(n)]
			readxf=True
		except: 
			pass
		
		if readxf:
			print("Read transform per mask from lst file. xform for the first particle:")
			for x in xfs:
				print(x)
				
			xf_last=np.zeros((len(data_cpx), nmask, 6))
			for il,l in enumerate(lst_raw):
				xf0=l["xform.projection"]
				for i in range(n):
					xf1=l[f"xform.projection_{i:02d}"]
					dxf=xf0.inverse()*xf1
					
					x=dxf.get_params("eman")
					xnp=np.array([x["az"],x["alt"],x["phi"], x["tx"], x["ty"], x["tz"]])
					xnp[:3]=xnp[:3]*np.pi/180.
					xnp[3:]/=float(options.rawbox)
					if options.focus:
						im=i+1
					else:
						im=i
					xf_last[il, im]=xnp
			xf_last=jnp.array(xf_last)		
	
	if readxf==False:
		xf_last=jnp.zeros((len(data_cpx), nmask, 6))
		
	xf_new=[]
	score=[]
	
	###########	
	dcpx=data_cpx[:bsz]
	xf=data_xf[:bsz]
	xf_rot=xf_last[:bsz]*0.
	
	l2wd=0
	if options.mask: 
		optimizer = optax.adam(options.learnrate)
		opt_state = optimizer.init(xf_rot)
		xf_rot00=xf_rot.copy()
		cost=[]
		for itr in range(options.niter):
			
			loss, grads=grad_fn_3d(xf_rot, xf, dcpx, l2wd)
			updates, opt_state = optimizer.update(grads, opt_state)
			xf_rot=optax.apply_updates(xf_rot, updates)
			
			cost.append(loss)
		l2=jnp.sum((xf_rot-xf_rot00)**2)
		l2wd=float((cost[0]-cost[-1])/l2)*options.l2_weight
		print("Refine without L2:")
		print("  loss change {:.3f}, L2 change {:.3f}, use weight decay factor {:.3f}".format(float((cost[0]-cost[-1])), float(l2), l2wd))
		
	###########
	for ib in range(0, len(data_cpx), bsz):
		dcpx=data_cpx[ib:ib+bsz]
		xf=data_xf[ib:ib+bsz]
		xf_rot=xf_last[ib:ib+bsz]
		xf_rot00=xf_rot.copy()*0
		
		optimizer = optax.adam(options.learnrate)
		opt_state = optimizer.init(xf_rot)
		cost=[]
		for itr in range(options.niter):
			
			loss, grads=grad_fn_3d(xf_rot, xf, dcpx, l2wd)
			updates, opt_state = optimizer.update(grads, opt_state)			
			xf_rot=optax.apply_updates(xf_rot, updates)
			
			cost.append(loss)
		
		xf_new.append(xf_rot)
		
		scr=calc_loss_3d_raw(xf_rot, xf, dcpx, 0, True)
		score.append(scr)
		
		sys.stdout.write("\r {:>8}/{}\t{:.3f} -> {:.3f}         ".format(ib, len(data_cpx), cost[0], cost[-1]))
		sys.stdout.flush()

	xf_new=jnp.vstack(xf_new)
	score=np.concatenate(score)
	return xf_new, score

def align_particles(data_cpx, data_xf, params, pts, lst_raw, options):
	#### DNN version
	
	def calc_loss_3d_raw(var, grd, dcpx, xf_proj, rng, return_all=False):
		#### Main function for optimization
		#### Both DNN rigid body and orientation search.
		#### Both SPA and SPT. A bit over complicated...
		train=not return_all
		if options.align_mlp:
			enc_var, dec_var=var
			conf = encode_model.apply(enc_var, grd, training=train, rngs={"dropout": rng})
			conf=0.01*jax.random.normal(rng, conf.shape)+conf
			xf_rot = decode_model.apply(dec_var, conf, training=train, rngs={"dropout": rng})
			if options.rand_ang>0 and train:
				xf_rot=np.deg2rad(options.rand_ang)*jax.random.normal(rng, xf_rot.shape)+xf_rot
		else:
			xf_rot=var
			conf=0
		
		#### one focus mask. main body does not move
		xf_rot=xf_rot*focus_mask[None,:,None]
		if options.fullrot>=0:
			xf_rot=xf_rot*mask_param[None,None,:]
			#print(xf_rot)		
			
		#### SPT mode. One rotation per 3D particle
		if options.spt:
			ntilt=dcpx.shape[1]
			dcpx=dcpx.reshape((dcpx.shape[0]*dcpx.shape[1],dcpx.shape[2], dcpx.shape[3]))
			xf_proj=xf_proj.reshape((xf_proj.shape[0]*xf_proj.shape[1],xf_proj.shape[2]))
			
		xr=xf_rot.reshape((-1,6))
		ptsx=jnp.repeat(pts[None,...], xr.shape[0], axis=0)
		
		pts1=xf2pts_multi(ptsx, xr)*np.array([1,-1,1,1,1])
		pts1=pts1.reshape((-1, nmask, len(pts), 5))
		pts1=pts1*imasks[None,:,:,None]
		pts1=jnp.sum(pts1, axis=1)
		if options.spt:
			pts1=jnp.repeat(pts1, ntilt, axis=0)
		
		imgs_cpx=pts2img(pts1, xf_proj)
		fval=calc_frc(dcpx, imgs_cpx, params["rings"], minpx=options.minpx, maxpx=options.maxpx)
		
		if return_all: 
			return xf_rot, fval, conf
	
		loss=-jnp.mean(fval)
			
		return loss
	
	@jit
	def calc_loss_init(var, grd, xf, rng):
	#### initial training to match given xf
		enc_var, dec_var=var
		conf = encode_model.apply(enc_var, grd, training=False, rngs={"dropout": rng})
		xf_rot = decode_model.apply(dec_var, conf, training=False, rngs={"dropout": rng})		
		loss=jnp.mean(abs(xf_rot-xf))
		return loss

	print("#########################")
	if options.spt:
		pids=np.array([a["ptcl3d_id"] for a in lst_raw])
		uid=np.unique(pids)
		idx3d=[np.where(pids==u)[0] for u in uid]
		nt=np.mean([len(i) for i in idx3d])
		options.batchsz=int(options.batchsz//nt)
		options.idx3d=np.array(idx3d)
		print("SPT mode. Average {} 2D image per subtomogram. Use batch size {}".format(nt, options.batchsz))
		nn=len(options.idx3d)
	else:
		nn=len(data_cpx)
		
	rng=jax.random.key(0)
	print(f"Align {nn} particles...")
	mask_param=np.ones(6)
	if options.fullrot>=0:
		mask_param[:3]=0
		mask_param[options.fullrot]=1
	mask_param=jnp.array(mask_param, dtype=float)
		
	#####################
	#### dealing with masks
	if options.mask==None:
		nmask=1
		imasks=jnp.ones((nmask, pts.shape[0]))
		focus_mask=np.array([1],dtype=float)
	else:
		masks=EMData.read_images(options.mask)		
		masks=np.array([make_mask_gmm(m, pts) for m in masks])
		masks=np.concatenate([np.zeros((1, masks.shape[1])), masks])
		masks[0]=0.5  #### threshold for regions outside any masks
		
		nmask=len(masks)
		imasks=np.zeros((nmask, len(pts)), dtype=float)
		imx=np.argmax(masks, axis=0)
		imasks[imx, np.arange(imasks.shape[1])]=1.
		print(f"Focused refinement with given masks")
		print(f"Gaussian coverage of unmasked regions and {nmask-1} masks. Sizes:", *np.sum(imasks, axis=1).astype(int))
		focus_mask=np.ones(nmask,dtype=float)
		focus_mask[0]=0
		
	#####################
	#### load previous alignment 
	readxf=False
	xf_last=np.zeros((nn, nmask, 6))
	if "xf_all" in lst_raw[0] and options.no_load==False:
		print("Loading existing local transforms...")
		xfs=np.array([l["xf_all"] for l in lst_raw])
		xfs=xfs.reshape((len(xfs), -1, 6))
		
		if options.spt:
			xfs=xfs[options.idx3d[:,0]]
			print("Taking one xform per 3D particle, shape:", xfs.shape)
		else:
			print("Existing xform shape:", xfs.shape)
		xf_last[:,1:xfs.shape[1]+1]=xfs
		readxf=True
		
	if options.align_mlp:
		if options.niter<0: options.niter=20
		if options.learnrate<0: options.learnrate=1e-5
		
		allgrds=calc_gradient(data_cpx, data_xf, pts, params, options)
		
		if options.gradout:
			print("Save gradient to file:", options.gradout)
			np.save(options.gradout, allgrds)
			#np.save(options.gradout[:-4]+"_xf.npy", xf_last)
			exit()
			
		encode_model, enc_var=build_encoder(pts, allgrds, rng, options)
		decode_model, dec_var=build_decoder(xf_last[0], rng, options)
		if options.decoderin:
			f=open(options.decoderin, "rb")
			var=pickle.load(f)
			f.close()
			var = from_state_dict([enc_var, dec_var], var)
			enc_var, dec_var=var
			
		var=enc_var, dec_var
		
		if readxf:
			bsz=options.batchsz*8
			print("#### Matching DNN to existing xform...")
			optimizer = optax.adam(1e-5)
			opt_state = optimizer.init(var)
			grad_fn_init=jax.value_and_grad(calc_loss_init)	
			for itr in range(100):
				cost=[]
				for ib in range(0, nn, bsz):
					grd=allgrds[ib:ib+bsz]
					xf=xf_last[ib:ib+bsz]
						
					rng=jax.random.fold_in(rng,123)
					loss, grads=grad_fn_init(var, grd, xf, rng)
					updates, opt_state = optimizer.update(grads, opt_state)
					var = optax.apply_updates(var, updates)
					cost.append(loss)
				print(itr, np.mean(cost), end='\r')
			print()
			
	else:
		if options.niter<0: options.niter=50
		if options.learnrate<0: options.learnrate=1e-3
		
	bsz=options.batchsz
	calc_loss_3d=jit(calc_loss_3d_raw)
	grad_fn_3d=jax.value_and_grad(calc_loss_3d)
	
	print("#### Training to match images...")
	####################
	#### DNN for rigid body movement 
	if options.align_mlp:
		allcost=[]
		optimizer = optax.adam(options.learnrate)
		opt_state = optimizer.init(var)
		for itr in range(options.niter):
			cost=[]
			
			for ib in range(0, nn, bsz):
				rng=jax.random.fold_in(rng,123)    
				grd=allgrds[ib:ib+bsz]
				
				if options.spt:
					i3=options.idx3d[ib:ib+bsz]
					dcpx=data_cpx[i3]
					xf=data_xf[i3]
				else:
					dcpx=data_cpx[ib:ib+bsz]
					xf=data_xf[ib:ib+bsz]
					
				loss, grads=grad_fn_3d(var, grd, dcpx, xf, rng)
				updates, opt_state = optimizer.update(grads, opt_state)
				var = optax.apply_updates(var, updates)
				cost.append(loss)            

				if (ib//bsz)%10==0: 
					sys.stdout.write("\r {}/{}\t{:.4f}      ".format(ib, len(data_cpx), loss))
					sys.stdout.flush()

			sys.stdout.write("\r")
			ac=np.mean(cost)
			allcost.append(ac)
			print("iter {},  loss : {:.4f}".format(itr, ac))		
		
		####################
		print("Calculating final xforms...")
		xf_new=[]
		score=[]
		conf=[]
		for ib in range(0, nn, bsz):
			grd=allgrds[ib:ib+bsz]
			if options.spt:
				i3=options.idx3d[ib:ib+bsz]
				dcpx=data_cpx[i3]
				xf=data_xf[i3]			
			else:
				dcpx=data_cpx[ib:ib+bsz]
				xf=data_xf[ib:ib+bsz]
			xf_rot, scr, cf=calc_loss_3d_raw(var, grd, dcpx, xf, rng, True)
			xf_new.append(xf_rot)
			score.append(scr)
			conf.append(cf)
			
			sys.stdout.write("\r {:>8}/{}\t{:.3f}      ".format(ib, len(data_cpx), np.mean(scr)))
			sys.stdout.flush()
			
		xf_new=jnp.vstack(xf_new)
		score=np.concatenate(score)
		conf=np.vstack(conf)
		print('\nfinal xform and score shape:',xf_new.shape, score.shape)
		if options.decoderout: save_network(var, options)
		if options.midout: 
			np.savetxt(options.midout, np.hstack([np.arange(len(conf))[:,None],conf]))
		
	####################
	#### Simple orientation search here 
	else:
		xf_new=[]
		score=[]
		for ib in range(0, nn, bsz):
			xf_rot=xf_last[ib:ib+bsz]
			if options.spt:
				i3=options.idx3d[ib:ib+bsz]
				dcpx=data_cpx[i3]
				xf=data_xf[i3]
			else:
				dcpx=data_cpx[ib:ib+bsz]
				xf=data_xf[ib:ib+bsz]
			
			optimizer = optax.adam(options.learnrate)
			opt_state = optimizer.init(xf_rot)
			cost=[]
			
			for itr in range(options.niter):
				loss, grads=grad_fn_3d(xf_rot, 0, dcpx, xf, rng)
				updates, opt_state = optimizer.update(grads, opt_state)			
				xf_rot=optax.apply_updates(xf_rot, updates)
				cost.append(loss)
			
			xf_new.append(xf_rot)
			xf_rot, scr, cf=calc_loss_3d_raw(xf_rot, 0, dcpx, xf, rng, True)
			score.append(scr)
			
			sys.stdout.write("\r {:>8}/{}\t{:.3f} -> {:.3f}         ".format(ib, len(data_cpx), cost[0], cost[-1]))
			sys.stdout.flush()

		xf_new=jnp.vstack(xf_new)
		score=np.concatenate(score)
	print(xf_new.shape)
	#print(xf_new[:,-1])
	return xf_new, score

def train_model(params, pts, rng, options):
	#### build GMM from projections
	#######################
	def calc_loss_raw(variables, dcpx, xf, rng):
		
		p0=jnp.zeros((xf.shape[0], len(pts), 5))+pts
		x=jnp.ones((xf.shape[0], options.nmid))
		
		pout = model.apply(variables, x, training=False, rngs={"dropout": rng})
		pout+=p0
				
		imgs_cpx=pts2img(pout, xf)
		fval=calc_frc(dcpx, imgs_cpx, params["rings"], minpx=options.minpx, maxpx=options.maxpx)
				
		loss=-jnp.mean(fval)
			
		return loss

	model, variables=build_decoder(pts, rng, options)
	
	######################
	#### making projections from 3d volume
	sym=Symmetries.get("c1")
	xfs00=sym.gen_orientations("eman",{"delta":4})

	xfs=[x.get_params("xyz") for x in xfs00]
	xfs=np.array([[x["ztilt"], x["ytilt"], x["xtilt"]] for x in xfs])
	xfs=np.deg2rad(xfs[:,::-1])
	
	e=EMData(options.ptclsin)
	sz0=e["nx"]
	e.process_inplace("math.fft.resample",{"n":sz0/options.maxboxsz})
	e.process_inplace("normalize.edgemean")
	
	vol=e.numpy().T.copy()
	sz=vol.shape[0]
	pad=sz//4
	vol=np.pad(vol, pad)
	vol_ft=jnp.fft.fftn(jnp.fft.fftshift(vol))
	vol_ft=jnp.fft.fftshift(vol_ft)
	print(f"Input volume size: {sz0}, downsampled to {sz}, padded to {vol.shape[0]}")
	
	global ind_xy, ind_interp, ind_trans, rawbox
	#################
	# index for rotation
	rawbox=vol.shape[0]
	ind_xy=np.indices((rawbox, rawbox, 1))
	ind_xy=ind_xy.T.reshape((-1,3))
	ind_xy[:,:2]-=rawbox//2

	#################
	# index for interpolation
	ind_interp=np.indices((2,2,2))
	ind_interp=ind_interp.T.reshape(-1,3)
	
	bsz=32
	projs=[]
	for ib in range(0, len(xfs), bsz):
		x=xfs[ib:ib+bsz]
		projs.append(get_slice_multi(vol_ft, x))

	projs=jnp.concatenate(projs, axis=0)
	print("projections shape", projs.shape)
	projs_real=jnp.fft.ifftn(jnp.fft.ifftshift(projs, axes=(1,2)), axes=(1,2)).real
	projs_real=jnp.fft.ifftshift(projs_real, axes=(1,2))
	projs_real=projs_real[:, pad:-pad, pad:-pad]
	print("projections shape", projs_real.shape)
	try:os.remove("tmp_projs.hdf")
	except:pass
	for i,p in enumerate(projs_real):
		e=from_numpy(np.array(p))
		e.write_image("tmp_projs.hdf",i)
	
	data_cpx=jnp.fft.rfft2(jnp.fft.fftshift(projs_real*1e-3, axes=(1,2))).astype(np.complex64)
	
	xfpm=[x.get_params("eman") for x in xfs00]
	data_xf=np.array([[x["az"],x["alt"],x["phi"], x["tx"], x["ty"], 0] for x in xfpm])
	data_xf[:,:3]=data_xf[:,:3]*np.pi/180.
	data_xf[:,3:]/=float(rawbox)
	
	print(data_cpx.shape, data_xf.shape, options.maxboxsz)
	print(params["rings"].shape)
	#######################
	bsz=options.batchsz
	dcpx=data_cpx[:bsz]
	xf=data_xf[:bsz]
	calc_loss=jit(calc_loss_raw)
	grad_fn=jax.value_and_grad(calc_loss)
	loss, grads=grad_fn(variables, dcpx, xf, rng)
	print("Initial loss:", loss)	
	
	#######################
	print("Now training...")
	allcost=[]
	optimizer = optax.adam(options.learnrate)
	opt_state = optimizer.init(variables)
	for itr in range(options.niter):
		cost=[]
		
		for ib in range(0, len(data_cpx), bsz):

			dcpx=data_cpx[ib:ib+bsz]
			xf=data_xf[ib:ib+bsz]			
			rng=jax.random.fold_in(rng,123)    
			
			loss, grads=grad_fn(variables, dcpx, xf, rng)
			updates, opt_state = optimizer.update(grads, opt_state)
			variables = optax.apply_updates(variables, updates)
			cost.append(loss)            

			if (ib//bsz)%10==0: 
				sys.stdout.write("\r {}/{}       {:.4f}      ".format(ib, len(data_cpx), loss))
				sys.stdout.flush()

		sys.stdout.write("\r")
		ac=np.mean(cost)
		allcost.append(ac)
		print("iter {},  loss : {:.4f}".format(itr, ac))
	
	
	x=jnp.ones((xf.shape[0], options.nmid))	
	pout = model.apply(variables, x, training=False, rngs={"dropout": rng})
	pout = pout[0]+pts
	return pout

def save_network(variables, options):
	f=open(options.decoderout, "wb")
	variables00 = to_state_dict(variables)
	pickle.dump(variables00, f)
	f.close()
	print("Trained networks saved to:",options.decoderout)

def make_mask_gmm(mask, pts):
	
	m=mask.numpy().copy()
	p=pts[:,:3].copy()
	p=p[:,::-1]
	p[:,:2]*=-1
	p=(p+.5)*mask["nx"]
	
	o=np.clip(p, 0, mask["nx"]-1)
	o=np.round(o).astype(int)
	
	v=m[o[:,0], o[:,1], o[:,2]]
	imsk=jnp.array(v)
	
	return imsk

###############################
###############################
def main():
	
	usage="""
	Single particle alignment, reconstruction and heterogeneity analysis with Gaussian model and neural networks. 
	
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--model", type=str,help="load from an existing model file", default="")
	parser.add_argument("--modelout", type=str,help="output trained model file. only used when --projs is provided", default="")
	
	parser.add_argument("--decoderin", type=str,help="Rather than initializing neural network model, read an existing trained one", default=None)
	parser.add_argument("--decoderout", type=str,help="Save the trained model", default=None)
	
	parser.add_argument("--encoderin", type=str,help="does nothing.", default=None)
	parser.add_argument("--encoderout", type=str,help="does nothing..", default=None)
		
	parser.add_argument("--ptclsin", type=str,help="particles input for alignment", default="")
	parser.add_argument("--ptclsout", type=str,help="aligned particle output", default="")
	
	parser.add_argument("--learnrate", type=float,help="learning rate for model training only. ", default=-1)	
	parser.add_argument("--niter", type=int,help="number of iterations", default=-1)
	parser.add_argument("--batchsz", type=int,help="batch size", default=32)

	parser.add_argument("--maxres", type=float,help="maximum resolution. will overwrite maxboxsz. ", default=-1)
	parser.add_argument("--minres", type=float,help="minimum resolution. ", default=500)
	parser.add_argument("--maxgradres", type=float,help="maximum resolution for gradient. ", default=-1)

	parser.add_argument("--trainmodel", action="store_true", default=False ,help="align particles.")
	parser.add_argument("--align", action="store_true", default=False ,help="align particles.")
	parser.add_argument("--align_mlp", action="store_true", default=False ,help="align particles with a MLP.")
	parser.add_argument("--heter", action="store_true", default=False ,help="heterogeneity analysis.")
	parser.add_argument("--mask", type=str,help="Single mask for heterogeneity analysis or stack of masks for multi-body alignment. ", default=None)

	parser.add_argument("--midout", type=str,help="middle layer output", default="")
	parser.add_argument("--gradout", type=str,help="save gradient input to file. for testing only", default=None)

	parser.add_argument("--pas", type=str,help="choose whether to adjust position, amplitude, sigma. use 3 digit 0/1 input. default is 111", default="111")
	
	parser.add_argument("--nmid", type=int,help="size of the middle layer", default=3)
	parser.add_argument("--ndense", type=int,help="size of the layers between the middle and in/out, variable if -1. Default 512", default=512)
	parser.add_argument("--l2_weight", type=float, help="weighting factor for L2. default is 1",default=1.)

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	parser.add_argument("--pointtransformer", action="store_true",help="use a point transformer based network. ", default=False)
	parser.add_argument("--pt_neighbor", type=int, help="number of neighbor for point transformer model",default=16)
	parser.add_argument("--rand_ang", type=float, help="add small random number to angle output",default=-1)
	
	parser.add_argument("--spt", action="store_true",help="subtomogram mode. ", default=False)
	parser.add_argument("--spt_ntilt", type=int, help="number of tilt images to keep for the --spt mode",default=-1)
	parser.add_argument("--clip", type=int, help="clip input images to this size before resampling",default=-1)
	parser.add_argument("--fullrot", type=int, help="full rotation for one of the channel. testing only",default=-1)
	parser.add_argument("--no_load", action="store_true",help="no loading previous local transforms from the xf_all key. ", default=False)


	(options, args) = parser.parse_args()
	logid=E2init(sys.argv,options.ppid)
	
	gen_model=None
	
	#### load GMM from text file
	if options.model:
		if options.model.endswith(".pdb"):
			print("Generating Gaussian model from pdb...")
			p=pdb2numpy(options.model, allatom=True)			
			pts=np.zeros((len(p),5), dtype=np.float32)
			e=EMData(options.ptclsin, 0, True)
			p=p/e["ny"]/e["apix_x"]-0.5
			p[:,1:]*=-1
			pts[:,:3]=p
			pts[:,3]=.5
			pts[:,4]=1
		else:
			pts=np.loadtxt(options.model).astype(np.float32)
			
		options.npt=npt=len(pts)
		print("{} Gaussian loaded from {}".format(len(pts), options.model))
	
	#### box size and resolution
	e=EMData(options.ptclsin, 0, True)
	raw_apix, raw_boxsz = e["apix_x"], e["ny"]
	options.raw_apix=raw_apix
	if options.clip>0: 
		options.rawbox=raw_boxsz=options.clip
	else:
		options.clip=options.rawbox=raw_boxsz
		
	if options.maxres>0:
		maxboxsz=options.maxboxsz=ceil(raw_boxsz*raw_apix*2/options.maxres)//2*2
		print("using box size {}, max resolution {:.1f}".format(maxboxsz, options.maxres))
	else:
		maxboxsz=options.maxboxsz=raw_boxsz
		
	options.maxpx=options.maxboxsz//2
	options.minpx=ceil(raw_boxsz*raw_apix*2/options.minres)//2
	options.minpx=max(1, options.minpx)
	print("FRC compares from {} to {} Fourier pixels".format(options.minpx, options.maxpx))
	
	if options.maxgradres>0:
		options.maxpx_grad=ceil(raw_boxsz*raw_apix*2/options.maxgradres)//2
	else:
		options.maxpx_grad=options.maxpx
		
	options.apix=apix=raw_apix*raw_boxsz/maxboxsz
	params=set_indices_boxsz(maxboxsz, options)
	
	
	#### for SPT
	options.id_keep=[]
	if options.spt:
		#### parse SPT particles
		lst_raw=load_lst_params(options.ptclsin)
		pids=np.array([a["ptcl3d_id"] for a in lst_raw])
		uid=np.unique(pids)
		idx3d=[np.where(pids==u)[0] for u in uid]

		ln=np.max([len(i) for i in idx3d])
		print("{} 3D particles, each contain {:.1f} 2D particles".format(len(idx3d), ln))
		
		options.id_keep=[]
		torm=0
		if options.spt_ntilt>0:
			print(f"Keeping {options.spt_ntilt} 2D tilt per 3D particle")
			id_keep=[]
			for idx in idx3d:
				tid=np.array([lst_raw[i]["tilt_id"] for i in idx])
				tid=np.argsort(abs(tid-np.mean(tid)))[:options.spt_ntilt]
				if len(tid)!=options.spt_ntilt: 
					torm+=1
					continue
				options.id_keep.extend(idx[np.sort(tid)].tolist())
		else:
			
			idx3d_keep=[i.tolist() for i in idx3d if len(i)==ln]
			torm=len(idx3d)-len(idx3d_keep)
			options.id_keep=sum(idx3d_keep,[])
			options.spt_ntilt=ln
			
		if torm>0:
			print(f"Skipping {torm} 3D particles that does not have enouth tilts")
					
		lst_raw=[lst_raw[i] for i in options.id_keep]
		
	
	########
	rng=jax.random.key(0)	

	#### Decoder training from generated projections of a 3-D map
	#### Note that train_decoder takes options.decoderentropy into account internally
	if options.trainmodel:
		print("training GMM...")
		if options.niter<0: options.niter=40
		if options.learnrate<0: options.learnrate=1e-5
		
		pts_new=train_model(params, pts, rng, options)
		np.savetxt(options.modelout, pts_new)
		
		E2end(logid)
		return
		
	data_cpx, data_xf, lst_raw = load_particles(options)
	
	#### Align particles using GMM
	if options.align or options.align_mlp:
		
		xf_new, score=align_particles(data_cpx, data_xf, params, pts, lst_raw, options)
		if options.mask: xf_new=xf_new[:,1:]
			
		if options.spt:
			pids=np.array([a["ptcl3d_id"] for a in lst_raw])
			uid=np.unique(pids)
			uid={u:i for i,u in enumerate(uid)}
			
			#### reload all particles
			lst_raw=load_lst_params(options.ptclsin)
			xf_new2=np.zeros((len(lst_raw), xf_new.shape[1], xf_new.shape[2]))
			score2=np.zeros(len(lst_raw))+np.mean(score)*.5
			score2[options.id_keep]=score
			skipid=[]
			for i,l in enumerate(lst_raw):
				ip=l["ptcl3d_id"]
				if ip in uid:
					ip=uid[ip]
				else:
					if ip not in skipid:
						skipid.append(ip)
						print("skip 3D particle",ip)
					continue
				
				xf_new2[i]=xf_new[ip]
				
			#xf_new=np.repeat(xf_new, options.spt_ntilt, axis=0)
			xf_new=xf_new2
			score=score2
			
		xf_raw=[Transform(x["xform.projection"]) for x in lst_raw]
			
		for i,l in enumerate(lst_raw):
			l["score"]=-score[i]
			if xf_new.shape[1]>0:
				l["xf_all"]=xf_new[i].flatten().tolist()
			
		for ix in range(xf_new.shape[1]):
			xx=np.array(xf_new[:,ix].copy())
			xx[:,:3]=xx[:,:3]*180./np.pi
			xx[:,3:]*=options.clip
			
			xf2=[Transform({"type":"eman", "az":x[0], "alt":x[1], 
					"phi":x[2], "tx":x[3], "ty":x[4], "tz":x[5]}) for x in xx.tolist()]

			xf_combine=[y*x for x,y in zip(xf2, xf_raw)]
			for i,l in enumerate(lst_raw):
				if options.mask==None:
					l[f"xform.projection"]=xf_combine[i]
				else:
					l[f"xform.projection_{ix:02d}"]=xf_combine[i]
				
		save_lst_params(lst_raw, options.ptclsout)
		print("\nDone. Output written to", options.ptclsout)

		E2end(logid)
		return
	
	
	#### Heterogeneity analysis from particles	
	if options.heter:
		if options.niter<0: options.niter=20
		if options.learnrate<0: options.learnrate=1e-5
		
		allgrds=calc_gradient(data_cpx, data_xf, pts, params, options)
		encode_model, enc_var=build_encoder(pts, allgrds, rng, options)
		if options.pointtransformer:
			decode_model, dec_var=build_decoder_point_transformer(pts, rng, options)
		else:
			decode_model, dec_var=build_decoder(pts, rng, options)
		
		if options.decoderin:
			f=open(options.decoderin, "rb")
			variables=pickle.load(f)
			f.close()
			variables = from_state_dict([enc_var, dec_var], variables)
			enc_var, dec_var=variables

		variables=[enc_var, dec_var]
		nvar=sum(x.size for x in jax.tree.leaves(enc_var["params"]))
		print("number of encoder weights:", nvar)
		nvar=sum(x.size for x in jax.tree.leaves(dec_var["params"]))
		print("number of encoder weights:", nvar)
		
		if options.mask:
		
			msk=EMData(options.mask)
			options.mask=make_mask_gmm(msk, pts)
			print("masking {} out of {} points".format(np.sum(options.mask), len(options.mask)))
		
		variables=train_heterg(data_cpx, data_xf, allgrds, variables, encode_model, decode_model, params, pts, rng, options)
		
		enc_var, dec_var=variables
		if options.decoderout: save_network(variables, options)		
		mid=calc_conf(encode_model, enc_var, allgrds, rng)
		if options.midout: np.savetxt(options.midout, mid)
		
		E2end(logid)
		return
	
	
if __name__ == '__main__':
	main()
	
