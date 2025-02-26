#!/usr/bin/env python

from EMAN2 import *
import numpy as np
from scipy.optimize import curve_fit

if "CUDA_VISIBLE_DEVICES" not in os.environ:
	os.environ["CUDA_VISIBLE_DEVICES"]='0' 
	
import jax
import jax.numpy as jnp

from flax import nnx
import flax.linen as nn
import optax
from jax import jit

def main():
	
	usage="e2tomogram_refine.py <tilt series> "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--res", type=float,help="target resolution. default is 50", default=50)
	parser.add_argument("--niter", type=int,help="number of iteration", default=100)
	parser.add_argument("--refinerot", action="store_true",help="refine tilt axis rotation", default=False)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	init=False
	ntilt=[EMUtil.get_image_count(f) for f in args]
	maxtilt=int(np.max(ntilt))
	print(f"{len(args)} tilt series. max {maxtilt} tilt images")
	for fname in args:
		print(base_name(fname))
		refine_tomo_align(fname, options, maxtilt, init)
		init=True
	
	E2end(logid)
	return 

def refine_tomo_align(fname, options, maxtilt, init=False):
	
	info=dict(js_open_dict(info_name(fname)))
	imgs_raw=EMData.read_images(fname)
	loss=np.array(info["ali_loss"])
	tltpms=np.array(info["tlt_params"])
	ikeep=np.where(loss<500)[0]
	print(f"using {len(ikeep)} out of {len(imgs_raw)} tilts")
	imgs_raw=[imgs_raw[i] for i in ikeep]
	tltpm00=tltpms.copy()
	tltpms=tltpms[ikeep]
	
	sz=min(imgs_raw[0]["nx"],imgs_raw[0]["ny"])
	apix=imgs_raw[0]["apix_x"]
	res=options.res
	xx=apix/res*sz*2
	xx=int((xx//40+1)*40)
	scale=sz/xx
	print(f"tilt series pixel size {apix:.2f}, target resolution {res}")
	print(f"Scale input by {scale:.2f}, image size {xx}")
	
	for m in imgs_raw:
		m.process_inplace("math.fft.resample",{"n":scale})
		m.process_inplace("filter.highpass.gauss", {"cutoff_pixels":4})
		m.process_inplace("filter.lowpass.gauss", {"cutoff_freq":1./res})
		m.process_inplace("normalize")
		m.clip_inplace(Region((m["nx"]-xx)//2, (m["ny"]-xx)//2, xx,xx))
		
	imgs=jnp.array([m.numpy().copy() for m in imgs_raw])
	#print(imgs.shape)
	nimg=len(imgs)
	global realbox, pad, rawbox
	realbox=128
	pad=16
	rawbox=realbox+pad*2

	ntile=3
	b=64
	xy=np.indices((ntile,ntile)).reshape((2,-1)).T.astype(float)
	xy=xy/(ntile-1)*(imgs.shape[-1]-b*2)+b
	xy-=imgs.shape[-1]/2
	
	global ind_xy, ind_interp, ind_trans
	##################
	## index for rotation
	ind_xy=np.indices((rawbox, rawbox, 1))
	ind_xy=ind_xy.T.reshape((-1,3))
	ind_xy[:,:2]-=rawbox//2

	##################
	## index for interpolation
	ind_interp=np.indices((2,2,2))
	ind_interp=ind_interp.T.reshape(-1,3)

	##################
	## index for translation
	ind_trans=np.indices((rawbox, rawbox))
	ind_trans=ind_trans-rawbox/2
	
	
	##################
	#print(tltpms[:,2])
	if options.refinerot:
		print(f"initial tilt axis rotation {tltpms[1,2]:.2f}")
		tltrot=refine_tlt_rot(imgs_raw, tltpms, scale)
		tltpms[:,2]=tltrot
	
	
	xf0=tltpms.copy()
	xf1=np.zeros_like(xf0)
	xf1[:,:3]=np.deg2rad(xf0[:,[4,3,2]])
	xf1[:,3:]=xf0[:,[1,0]]
	allxf=xf1.copy()
	allxf=np.pad(allxf, [[0,maxtilt-nimg],[0,0]])
	#print(xf1.shape, tltpms.shape, allxf.shape)
	trans1=allxf[:,3:]/scale

	angoffset=np.zeros(3)
	
	##################
	mask_tilt=np.zeros(maxtilt)
	mask_tilt[:nimg]=1
	xf=allxf[:,:3]+angoffset
	data_cpx, imgs_clip, trans_offset=get_clips(imgs, xf, xy)
	
	b=int(realbox*.25)
	wid0=realbox//2-b//2
	wid1=realbox//2+b//2
	width_mask=np.zeros(realbox)
	width_mask[wid0:wid1]=1
	
	trans2=trans1.copy()
	global grad_fn_multi
	if init==False:    
		loss=calc_loss_multi(trans2, allxf[:,:3]+angoffset, data_cpx, trans_offset, mask_tilt, width_mask)
		#print(loss)
		grad_fn_multi=jax.value_and_grad(calc_loss_multi)
		loss, grads=grad_fn_multi(trans2, allxf[:,:3]+angoffset, data_cpx, trans_offset, mask_tilt, width_mask)
		
	##################
	optimizer = optax.adam(3e-2)
	opt_state = optimizer.init(trans2)
	cost=[]
	loss=calc_loss_multi(trans2, allxf[:,:3]+angoffset, data_cpx, trans_offset, mask_tilt, width_mask)
	print(f"initial loss: {loss:.4f}")
	
	for i in range(options.niter):
		loss, grads=grad_fn_multi(trans2, allxf[:,:3]+angoffset, data_cpx, trans_offset, mask_tilt, width_mask)
		updates, opt_state = optimizer.update(grads, opt_state)
		trans2 = optax.apply_updates(trans2, updates)
		
		if i%5==0:print(i, loss, end='\r')
		cost.append(loss)

	print(f"\nloss : {cost[0]:.3f} -> {cost[-1]:.3f}")
	##################
	stds2=[]
	for j in [0,1]:
		stds=[]
		if j==0:
			tt=trans1
		else:
			tt=trans2
		for ii,data in enumerate(data_cpx):
			volume=jnp.zeros((rawbox, rawbox, rawbox), dtype=np.complex64)
			weight=jnp.zeros((rawbox, rawbox, rawbox), dtype=np.float32)
		
			xf=jnp.concatenate([allxf[:,:3]+angoffset, tt+trans_offset[ii]], axis=1)
			volume, weight=insert_slice_multi(volume, weight, data, xf, mask_tilt)
			vol_nrm=get_volume(volume, weight)
			
			
			e=from_numpy(np.array(vol_nrm).T.copy())
			e.write_image(f"tmp_tomo_{j}.hdf", ii)

			
			ss=jnp.std(vol_nrm, axis=(0,1))
			stds.append(ss)
		
		stds=np.array(stds)
		stds2.append(stds)


	sdiff=[np.mean(s[:, wid0:wid1], axis=1) for s in stds2]
	for i in range(len(data_cpx)):
		s0=sdiff[0][i]
		s1=sdiff[1][i]
		print(f"tile {i} std: {s0:.3f} -> {s1:.3f}, diff {s1-s0:+.3f}")
		
	s0=np.mean(sdiff[0])
	s1=np.mean(sdiff[1])
	print(f"mean std: {s0:.3f} -> {s1:.3f}, diff {s1-s0:+.3f}")

	###################

	xf1=tltpms.copy()
	xf1[:,[1,0]]=trans2[:nimg]*scale
	xf1[:,[4,3,2]]+=np.rad2deg(angoffset)
	
	xf2=tltpm00.copy()
	xf2[ikeep]=xf1
	
	xf2=np.hstack([np.arange(len(xf2))[:,None], xf2])
	
	bm=base_name(fname)
	pmname=f"tmp_tltparams_{bm}.txt"
	np.savetxt(pmname, xf2)
	
	run(f"e2tomogram.py {fname} --bytile --niter 0 --noali --loadfile {pmname} --filterres {res} --notmp --normslice")
	return

def refine_tlt_rot(imgs, tltpms, scale):
	
	trans=tltpms[:,[0,1]]/scale
	imgs2=[]
	for i,im in enumerate(imgs):
		t=np.round(-trans[i]).astype(int).tolist()
		m=im.copy()    
		m.process_inplace("normalize.edgemean")
		m.process_inplace("xform", {"tx":t[0], "ty":t[1]})
		m.process_inplace("mask.soft",{"outer_radius":-16, "width":8})
		imgs2.append(m)
		
	data2=np.array([m.numpy().copy() for m in imgs2])
	sz=data2.shape[-1]
	ind=np.indices((sz,sz))
	ind-=sz//2
	r=np.linalg.norm(ind, axis=0)
	rmask=(r>sz//20).astype(float)
	
	ang=np.arctan2(ind[1], ind[0])
	ang=np.rad2deg(ang)%360
	
	data_fft=np.fft.fftshift(data2, axes=(1,2))
	data_fft=np.fft.fft2(data_fft, axes=(1,2))
	data_fft=np.fft.fftshift(data_fft, axes=(1,2))
	data_fft_abs=abs(data_fft)
	data_fft=data_fft/data_fft_abs

	sm2=np.mean(data_fft, axis=0)
	sm2=sm2*rmask
	
	nimg=len(tltpms)
	amax=np.round(tltpms[nimg//2,2]%360)
	
	w=20
	angrng=np.arange(amax-w, amax+w, .1)

	sm=abs(sm2).copy()
	cm2=[]
	for a in angrng:
		d=(ang-a)%360
		d=np.minimum(d, 360-d)
		c=np.mean(sm[d<1])
		cm2.append(c)
	cm2=np.array(cm2)
	cm2-=np.min(cm2)

	amax=angrng[np.argmax(cm2)]
	print("fine tilt axis search:",amax)
	
	p0=[np.max(cm2), angrng[np.argmax(cm2)], 1]
	coeff, var_matrix = curve_fit(fit_func, angrng, cm2, p0=p0)
	gfit = fit_func(angrng, *coeff)
	amax=coeff[1]
	print("curve fit:",amax)
	return amax

def fit_func(x, *p):
	A, mu, sigma = p
	y=sigma/(abs(x-mu)+sigma)
	y=y-np.min(y)
	y=y*A
	return y

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
	
	# R = R_z @ R_y @ R_x  # Combined rotation
	R = R_x @ R_y @ R_z  # Combined rotation
	# R = R_z @ R_x  # Combined rotation
	return R
	

@jit
def insert_slice(volume, weight, data, ang, mask=1):
	R=make_matrix_3d(ang)
	rotated_coords = R @ ind_xy.T
	rotated_coords+=rawbox//2
	rotated_coords=jnp.concatenate([rotated_coords[:,:2], abs(rotated_coords[:,2:])], axis=1)
	
	
	ix=jnp.floor(rotated_coords).astype(int)
	ir=rotated_coords-jnp.floor(rotated_coords)
	ir=[ir, 1-ir]
	
	x=ind_trans[0]*ang[3]+ind_trans[1]*ang[4]
	x=x.astype(np.complex64)
	x=jnp.exp(1j*np.pi*2.*x/rawbox)
	ds=(data*x).reshape((-1,))
		
	for i0 in ind_interp[...,None]:
		w=ir[0]*i0+ir[1]*(1-i0)
		w=jnp.prod(w, axis=0)
		w=w*mask
		i=ix+i0
		volume=volume.at[i[0], i[1], i[2]].add(ds*w)
		weight=weight.at[i[0], i[1], i[2]].add(w)

	return volume, weight

@jit
def insert_slice_multi(volume, weight, data, ang, mask=np.ones(100)):
	for d,a,m in zip(data, ang, mask):
		volume, weight=insert_slice(volume, weight, d, a, m)
	return volume, weight

@jit
def get_volume(volume, weight):
	# wt=weight.at[weight==0].set(1)
	wt=jnp.where(weight==0, 1, weight)
	vol_nrm=volume/wt
	vol_nrm=jnp.fft.ifftshift(vol_nrm)
	vol_nrm=jnp.fft.ifftshift(jnp.fft.ifftn(vol_nrm)).real
	vol_nrm=vol_nrm[pad:-pad, pad:-pad, pad:-pad]
	vol_nrm-=jnp.mean(vol_nrm[...,[0,-1]])
	vol_nrm=-nn.leaky_relu(-vol_nrm, 0.5)
	vol_nrm/=jnp.std(vol_nrm[...,[0,-1]])
	
	std_noise=2-(np.arange(realbox)/(realbox-1)*2-1)**2
	vol_nrm/=std_noise

	return vol_nrm

def get_clips(imgs, xf, xy):
	
	m=EMData(realbox, realbox)
	m.to_one()
	m.process_inplace("mask.decayedge2d", {"width":8})
	mask_edge=m.numpy().copy()
	mask_edge=np.pad(mask_edge, pad)
	
	imgs_clip=[]
	trans_offset=[]
	for pp in xy:
		pos=np.append(pp,0)
		pos2=np.array([(make_matrix_3d(x).T)@pos for x in xf])
		pos2=pos2[:,[1,0]]
		pos2=pos2+imgs.shape[-1]//2
		pos3=pos2[:,:2]
		pos3=np.round(pos3).astype(int)
		dx=pos2[:,:2]-pos3
		trans_offset.append(dx)
		pos3=np.clip(pos3, realbox//2, imgs.shape[-1]-realbox//2)
		
		m=np.array([imgs[i,p[0]-realbox//2:p[0]+realbox//2, p[1]-realbox//2:p[1]+realbox//2] for i,p in enumerate(pos3)])
		m=np.pad(m, [[0,0],[pad,pad],[pad,pad]] )
		m*=mask_edge
		imgs_clip.append(m)
		
	
	imgs_clip=np.array(imgs_clip)
	trans_offset=np.array(trans_offset)
	
	data_cpx=jnp.fft.fftshift(imgs_clip, axes=(2,3))
	data_cpx=jnp.fft.fft2(data_cpx)
	data_cpx=jnp.fft.fftshift(data_cpx, axes=(2,3))
	# data_cpx=np.pad(data_cpx, [[0,0],[0,maxtilt-nimg],[0,0],[0,0]])
	return data_cpx, imgs_clip, trans_offset

@jit
def calc_loss_multi(trans, xfrot, data_all, trans_offset, mask_tilt, width):
	score=[]
	for ii,data in enumerate(data_all):
		volume=jnp.zeros((rawbox, rawbox, rawbox), dtype=np.complex64)
		weight=jnp.zeros((rawbox, rawbox, rawbox), dtype=np.float32)
	
		xf=jnp.concatenate([xfrot, trans+trans_offset[ii]], axis=1)
		volume, weight=insert_slice_multi(volume, weight, data, xf, mask_tilt)
		vol_nrm=get_volume(volume, weight)
		
		std=jnp.std(vol_nrm, axis=(0,1))
		std=jnp.sum(std*width)/jnp.sum(width)
		score.append(std)

	score=jnp.array(score)
	loss=jnp.sum(score)
	loss-=jnp.min(score)
	loss=-loss/(len(data_all)-1)
	
	return loss


if __name__ == '__main__':
	main()
