#!/usr/bin/env python
# Muyuan Chen 2025-12
## jax reimplementation of gmm pipeline
## with new stuff..
from EMAN2 import *
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import scipy.spatial.distance as scipydist
from flax.serialization import to_state_dict, from_state_dict
import pickle

if "CUDA_VISIBLE_DEVICES" not in os.environ:
	# so we can decide which gpu to use with environmental variable
	os.environ["CUDA_VISIBLE_DEVICES"]='0' 
import jax
import jax.numpy as jnp
	
from flax import nnx
import flax.linen as nn
import optax
from jax import jit

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
	pts_rot_trans=jnp.stack([(pts_rot[:,0]+ang[3]), (-pts_rot[:,1])+ang[4]], 1)

	return pts_rot_trans


def pts2img_one_sig(pts, ang):
	# pts, ang = args
	sz, idxft, rrft = params["sz"], params["idxft"], params["rrft"]
	xfo=params["xforigin"]
	
	bamp=jax.nn.relu(pts[:,3]).astype(np.complex64)[:,None]
	bsigma=jax.nn.relu(pts[:,4]).astype(np.complex64)*.001

	bpos=xf2pts(pts[:,:3],ang)
	bpos=bpos*sz+sz/2
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

#### load particles from file and fourier transform them
#### particles need to have their transform in file header or comment of list file
##   will also shrink particles and do FT
##   return Fourier images and transform matrix
def load_particles(options):
	fname=options.ptclsin
	boxsz=options.maxboxsz
	
	nptcl=EMUtil.get_image_count(fname)
	projs=[]
	hdrs=[]
	e=EMData(fname, 0, True)
	rawbox=e["nx"]
	print("Loading {} particles of box size {}. shrink to {}".format(nptcl, rawbox, boxsz))
	for j in range(0,nptcl,1000):
		print(f"\r {j}/{nptcl} ",end="")
		sys.stdout.flush()
		el=EMData.read_images(fname,range(j,min(j+1000,nptcl)))
		print(f"R     ",end="")
		sys.stdout.flush()
		for e in el:
			if rawbox!=boxsz:
				e.process_inplace("math.fft.resample",{"n":e["nx"]/boxsz})
			
			hdrs.append(e.get_attr_dict())
			projs.append(e.numpy().copy())
	print(f"{nptcl}/{nptcl}")
	projs=np.array(projs, dtype=np.float32)
	
	data_cpx=np.fft.rfft2(np.fft.fftshift(projs*1e-3, axes=(1,2))).astype(np.complex64)

	xflst=False
	if fname.endswith(".lst"):
		info=load_lst_params(fname)
		if "xform.projection" in info[0]:
			xflst=True
			xfs=[p["xform.projection"].get_params("eman") for p in info]
			
	if xflst==False and ("xform.projection" in hdrs[0]):
		xflst=True
		xfs=[p["xform.projection"].get_params("eman") for p in hdrs]
		
	if xflst==False:
		xfs=[Transform().get_params("eman") for p in hdrs]
		print("No existing transform from particles...")
		
	xfsnp=np.array([[x["az"],x["alt"],x["phi"], x["tx"], x["ty"]] for x in xfs])
	xfsnp[:,:3]=xfsnp[:,:3]*np.pi/180.
	xfsnp[:,3:]/=float(rawbox)
	
	print("Data read complete")
	return data_cpx,xfsnp
	

#### compute fourier indices for image generation, clipping, and frc
#### pass indices in a dictionary
def set_indices_boxsz(boxsz, apix=0, return_freq=False):
	idx=np.indices((boxsz,boxsz))-boxsz//2
	idx=np.fft.ifftshift(idx)
	idx=idx[:,:,:boxsz//2+1]
	
	sz=boxsz
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
	params={"sz":sz, "idxft":idxft, "rrft":rrft, "rings":rings, "xforigin":xvec}
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
	minpx,maxpx=options.minpx, options.maxpx
	bsz=options.batchsz

	def pt_to_frc(pt, xf, pj_cpx):
		imgs_cpx=pts2img(pt, xf)
		fval=calc_frc(pj_cpx, imgs_cpx, params["rings"], minpx=minpx, maxpx=maxpx)
		return -jnp.mean(fval)
		
	grad_fn_00=jax.value_and_grad(pt_to_frc)

	for ib in range(0, len(data_cpx), bsz):
		
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
				x = nn.Dropout(0.1, deterministic=not training)(x)
				
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
	npt=pts.shape[0]

	class Decoder(nn.Module):
		@nn.compact
		def __call__(self, x, training=False):
			for i in range(3):
				x = nn.Dense(nhidden)(x)
				x = nn.relu(x)
				x = nn.Dropout(0.1, deterministic=not training)(x)
				
			x = nn.Dense(npt*5, kernel_init=nn.initializers.normal(stddev=1e-3))(x)
			x = jnp.reshape(x, (x.shape[0],npt, 5))
			x*=0.5
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

def save_network(variables, options):
	f=open(options.nnet_out, "wb")
	variables00 = to_state_dict(variables)
	pickle.dump(variables00, f)
	f.close()
	print("Trained networks saved to:",options.nnet_out)

###############################
###############################
def main():
	
	usage="""Single particle alignment, reconstruction and heterogeneity analysis with Gaussian model and neural networks. 
	
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--model", type=str,help="load from an existing model file", default="")
	parser.add_argument("--modelout", type=str,help="output trained model file. only used when --projs is provided", default="")
	
	parser.add_argument("--nnet_in", type=str,help="Rather than initializing neural network model, read an existing trained one", default="")
	parser.add_argument("--nnet_out", type=str,help="Save the trained model", default=None)
	
	parser.add_argument("--evalmodel", type=str,help="generate model projection images to the given file name", default="")
	
	parser.add_argument("--ptclsin", type=str,help="particles input for alignment", default="")
	parser.add_argument("--ptclsout", type=str,help="aligned particle output", default="")
	
	parser.add_argument("--learnrate", type=float,help="learning rate for model training only. Default is 1e-5. ", default=1e-5)	
	parser.add_argument("--niter", type=int,help="number of iterations", default=20)
	parser.add_argument("--batchsz", type=int,help="batch size", default=32)
	parser.add_argument("--maxboxsz", type=int,help="maximum fourier box size to use. 2 x target Fourier radius. ", default=64)
	parser.add_argument("--maxres", type=float,help="maximum resolution. will overwrite maxboxsz. ", default=-1)
	parser.add_argument("--maxgradres", type=float,help="maximum resolution for gradient. ", default=-1)
	parser.add_argument("--minres", type=float,help="minimum resolution. ", default=500)
	
	parser.add_argument("--trainmodel", action="store_true", default=False ,help="align particles.")
	parser.add_argument("--align", action="store_true", default=False ,help="align particles.")
	parser.add_argument("--heter", action="store_true", default=False ,help="heterogeneity analysis.")

	parser.add_argument("--midout", type=str,help="middle layer output", default="")
	parser.add_argument("--pas", type=str,help="choose whether to adjust position, amplitude, sigma. use 3 digit 0/1 input. default is 111", default="111")
	
	parser.add_argument("--nmid", type=int,help="size of the middle layer", default=3)
	parser.add_argument("--ndense", type=int,help="size of the layers between the middle and in/out, variable if -1. Default 512", default=512)

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	parser.add_argument("--pointtransformer", action="store_true",help="use a point transformer based network. ", default=False)
	parser.add_argument("--pt_neighbor", type=int, help="number of neighbor for point transformer model",default=16)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv,options.ppid)
	
	gen_model=None
	maxboxsz=options.maxboxsz
	
	## load GMM from text file
	if options.model:
		pts=np.loadtxt(options.model).astype(np.float32)
		options.npt=npt=len(pts)
		print("{} Gaussian loaded from {}".format(len(pts), options.model))
		
	
	#### load particles
	e=EMData(options.ptclsin, 0, True)
	raw_apix, raw_boxsz = e["apix_x"], e["ny"]
	options.raw_apix=raw_apix
	if options.maxres>0:
		maxboxsz=options.maxboxsz=ceil(raw_boxsz*raw_apix*2/options.maxres)//2*2
		print("using box size {}, max resolution {:.1f}".format(maxboxsz, options.maxres))

	options.maxpx=options.maxboxsz//2
	options.minpx=ceil(raw_boxsz*raw_apix*2/options.minres)//2
	options.minpx=max(1, options.minpx)
	print("FRC compares from {} to {} Fourier pixels".format(options.minpx, options.maxpx))
	
	data_cpx, data_xf = load_particles(options)
	options.apix=apix=raw_apix*raw_boxsz/maxboxsz
	params=set_indices_boxsz(maxboxsz)
	
	########
	rng=jax.random.key(0)	

	#### Decoder training from generated projections of a 3-D map or particles
	#### Note that train_decoder takes options.decoderentropy into account internally
	if options.trainmodel:
		pass
		
		
	#### Align particles using GMM
	if options.align:
		pass
		
	#### Heterogeneity analysis from particles	
	if options.heter:
		allgrds=calc_gradient(data_cpx, data_xf, pts, params, options)
		encode_model, enc_var=build_encoder(pts, allgrds, rng, options)
		if options.pointtransformer:
			decode_model, dec_var=build_decoder_point_transformer(pts, rng, options)
		else:
			decode_model, dec_var=build_decoder(pts, rng, options)
		
		variables=[enc_var, dec_var]
		nvar=sum(x.size for x in jax.tree.leaves(enc_var["params"]))
		print("number of encoder weights:", nvar)
		nvar=sum(x.size for x in jax.tree.leaves(dec_var["params"]))
		print("number of encoder weights:", nvar)
		
		variables=train_heterg(data_cpx, data_xf, allgrds, variables, encode_model, decode_model, params, pts, rng, options)
		
		enc_var, dec_var=variables
		if options.nnet_out: save_network(variables, options)		
		mid=calc_conf(encode_model, enc_var, allgrds, rng)
		if options.midout: np.savetxt(options.midout, mid)
		
	E2end(logid)
	
	
if __name__ == '__main__':
	main()
	
