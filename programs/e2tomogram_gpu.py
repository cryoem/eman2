#!/usr/bin/env python
# Muyuan Chen 2024-04
from EMAN2 import *
import numpy as np

if "CUDA_VISIBLE_DEVICES" not in os.environ:
    # so we can decide which gpu to use with environmental variable
    os.environ["CUDA_VISIBLE_DEVICES"]='0' 

#### do not occupy the entire GPU memory at once
##   seems necessary to avoid some errors...
os.environ["TF_FORCE_GPU_ALLOW_GROWTH"]='true' 
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' #### reduce log output

#### finally initialize tensorflow
import tensorflow as tf

floattype=np.float32
def main():
	
	usage="""
	This program takes unaligned tilt series, performs alignment, and generate a tomogram.
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--tltstep", type=float,help="Step between tilts. Ignored when rawtlt/mdoc is provided. Set to 0 if tilt present in header.", default=3.0)
	parser.add_argument("--tltax", type=float,help="Angle of the tilt axis.", default=None)
	(options, args) = parser.parse_args()

	logid=E2init(sys.argv)
	time0=time.time()
	path=os.path.dirname(args[0])
	##############################
	imgs=EMData.read_images(args[0])
	imgs=np.array([i.numpy().T.copy() for i in imgs], dtype=floattype)
	print("Input image shape:",imgs.shape)
	
	tltpm00=np.loadtxt(args[1])[:,1:]
	tltpm00[:,:2]/=16
	
	##############################
	global ind_xy, ind_interp, ind_trans,ind_vol, sz

	sz=256
	##################
	## index for rotation
	ind_xy=np.indices((sz, sz, 1))
	ind_xy=ind_xy.transpose((1,2,3,0)).reshape((-1,3))
	ind_xy=ind_xy.astype(floattype)
	ind_xy[:,:2]-=sz//2

	##################
	## index for interpolation
	ind_interp=np.indices((2,2,2))
	ind_interp=ind_interp.transpose((1,2,3,0)).reshape(-1,3)

	##################
	## index for translation
	ind_trans=np.indices((sz, sz))
	ind_trans=ind_trans.astype(np.float32)-sz/2
	
	##################
	## index for volume insertion
	ind_vol=np.indices((sz,sz,sz))
	ind_vol=ind_vol.transpose((1,2,3,0)).reshape((-1,3))
	
	
	##############################
	## Coarse alignment
	imgs_bin4=tf.keras.layers.AveragePooling2D(pool_size=(4,4))(imgs[...,None])[...,0]

	imsk=np.indices((256,256))/256-.5
	imsk=np.sum(imsk**2, axis=0)
	imsk=np.exp(-imsk*20)

	imgs_bin4=imgs_bin4-tf.reduce_mean(imgs_bin4, axis=(1,2))[:,None,None]
	imgs_bin4_ali=imgs_bin4*imsk[None,...]
	
	alltrans=tf.zeros((len(imgs),2), dtype=floattype)
	for itr in range(3):
		ts=coarse_align(imgs_bin4_ali)
		alltrans-=ts
		imgs_bin4_ali=tf.image.extract_glimpse(
			imgs_bin4[...,None],
			size=(256, 256),
			offsets=alltrans,
			centered=False, normalized=False, noise='zero')[...,0]
		imgs_bin4_ali*=imsk
		
		print("!!!",np.mean(np.linalg.norm(alltrans-tltpm00[:,:2], axis=1)))
		
	alltrans_coarse=alltrans.numpy().copy()
	dtime=time.time()-time0
	print("\nCoarse alignment done - {:.1f}s".format(dtime))
	
	##############################
	tltpm=np.zeros((len(imgs), 5), dtype=floattype)
	tltpm[:,2]=options.tltax
	tltpm[:,3]=np.arange(len(imgs))*options.tltstep
	tltpm[:,3]-=np.median(tltpm[:,3])
	tltpm[:,:2]=alltrans_coarse
	
	##############################
	## Patch alignment 2D
	imgs_bin2=tf.keras.layers.AveragePooling2D(pool_size=(2,2))(imgs[...,None])[...,0]
	tltpm_bin2=tltpm.copy()
	tltpm_bin2[:,:2]*=2
	
	### position of tiles
	ntile=7
	tid=np.indices((ntile, ntile)).T.reshape(-1,2)
	tid=tid/(ntile-1)-.5
	tid=tid*(imgs_bin2.shape[1]*1.0-sz)
	tpos=np.hstack([tid, np.zeros((len(tid),1))])
	tpos=tf.convert_to_tensor(tpos, dtype=floattype)
	
	for itr in range(3):
		all_ts=patch_track(imgs_bin2, tltpm_bin2, tpos, imsk)
		print(itr, np.mean(np.linalg.norm(all_ts, axis=2)))
		tltpm_bin2[:,:2]-=np.mean(all_ts, axis=0)
		print("!!!",np.mean(np.linalg.norm(tltpm_bin2[:,:2]/2-tltpm00[:,:2], axis=1)))
	
	##############################
	## Refine rotation
	for full in [0]:
		tpm=refine_rotation(tltpm_bin2, all_ts, tpos, fullrefine=full)
		tltpm_bin2=tpm.numpy().copy()
		for j in range(2):
			all_ts=patch_track(imgs_bin2, tltpm_bin2, tpos, imsk)
			tltpm_bin2[:,:2]-=np.mean(all_ts, axis=0)
		print(np.mean(np.linalg.norm(all_ts, axis=2)))
		print("!!!",np.mean(np.linalg.norm(tltpm_bin2[:,:2]/2-tltpm00[:,:2], axis=1)))
		
	tltpm=tltpm_bin2.copy()
	tltpm[:,:2]*=2
	
	dtime=time.time()-time0
	print("\nPatch track done - {:.1f}s".format(dtime))
	
	##############################
	for itr in range(3):
		all_ts=patch_track_3d(imgs_bin2, tltpm_bin2, tpos, imsk)
		tltpm_bin2[:,:2]-=all_ts
		print("!!!",np.mean(np.linalg.norm(tltpm_bin2[:,:2]/2-tltpm00[:,:2], axis=1)))
	
	dtime=time.time()-time0
	print("\n3D Patch track done - {:.1f}s".format(dtime))
	##############################
	
	# print(tltpm[:,3])
	tpm=np.hstack([np.arange(len(tltpm))[:,None], tltpm])
	tpm[:,1:3]*=4
	np.savetxt(f"{path}/tltparams_05.txt", tpm)
	
	imgs_ali=tf.image.extract_glimpse(imgs_bin2[...,None], size=(512, 512), offsets=tltpm[:,:2]/2., centered=False, normalized=False, noise='zero')[...,0]
	# print(tltpm[:,2])
	for i,m in enumerate(imgs_ali):
		e=from_numpy(m.numpy().T.copy())
		e.process_inplace("xform",{ "alpha":-float(tltpm[i,2])})
		e.write_image(f"{path}/twod_05.hdf", i)
		
	# return
	##############################
	## Final reconstruction...
	outx=outy=1024
	step=sz//2
	nstepx=int(outx/step/2)
	nstepy=int(outy/step/2)
	moretile=False

	allpk=[]
	tilepos=[]
	for stepx in range(-nstepx,nstepx+1):
		#### shift y by half a tile
		if moretile:
			yrange=range(-nstepy,nstepy+1)
		else: 
			yrange=range(-nstepy+stepx%2,nstepy+1,2)
		for stepy in yrange:
			pk=np.array([stepx*step,stepy*step,0])
			tilepos.append(pk)
			pk2=np.array(get_xf_pos(tltpm, pk.astype(floattype)), dtype=floattype)
			pk2=pk2+imgs.shape[1]//2-sz//2
			allpk.append(pk2)
			
	allpk=np.array(allpk)

	pkint=tf.round(allpk)
	pkres=allpk-tf.round(allpk)
	
	##############################
	tpos=np.array(tilepos)
	tpos=tpos+[outx//2, outy//2, sz//2]-sz//2
	fullvol=tf.zeros((outx, outy, sz), dtype=np.float32)
	
	for ii,pk in enumerate(pkint):
		fullvol=make_tile(imgs, pk, tpos[ii], tltpm, fullvol)
		print(ii, tpos[ii], end='\r')
		
		
	dtime=time.time()-time0
	print("\nReconstruction done - {:.1f}s".format(dtime))
	
	r=fullvol.numpy()
	b=from_numpy(r.transpose(2,1,0).copy())
	b.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.25})
	b.process_inplace("normalize.edgemean")
	b.write_compressed(args[2], 0, 8, nooutliers=True)
	
	
	dtime=time.time()-time0
	print("Finished. Total time: {:.1f}s".format(dtime))

	E2end(logid)
	return
	

def patch_track_3d(imgs_now, tltpm_start, tpos_now, imsk):
	

	trange=np.argsort(abs(tltpm_start[:,3]))

	vols=[tf.zeros((sz, sz, sz), dtype=np.complex64) for t in tpos_now]
	wts=[tf.zeros((sz, sz, sz), dtype=np.float32) for t in tpos_now]

	pk2=tf.convert_to_tensor([get_xf_pos(tltpm_start, tp)+imgs_now.shape[1]/2.-sz/2. for tp in tpos_now])
	alltrans=np.zeros((len(imgs_now), 2), dtype=floattype)
	for tid in trange:
		imgs_patch=[tf.image.extract_glimpse(imgs_now[tid:tid+1,...,None], size=(sz, sz), offsets=p[tid:tid+1], centered=False, normalized=False, noise='zero')[0,:,:,0] for p in pk2]
		imgs_patch=tf.convert_to_tensor(imgs_patch)
		imgs_patch=imgs_patch*imsk
		# print(imgs_patch.shape)
		
		data_cpx=tf.cast(imgs_patch, tf.complex64)*1e-3
		data_cpx=tf.signal.fftshift(data_cpx, axes=(1,2))
		data_cpx=tf.signal.fft2d(data_cpx)
			
		tpm_rot=tltpm_start[tid].copy()
		tpm_rot[:2]*=0
		
		trans=tf.zeros((2,), dtype=floattype)
		if tid!=trange[0]:
			pj_cpx=[]
			for iv in range(len(vols)):
				v=tf.math.divide_no_nan(vols[iv], tf.cast(wts[iv], tf.complex64))
				p=do_project(tpm_rot, v)
				pj_cpx.append(p)

			pj_cpx=tf.convert_to_tensor(pj_cpx)
			pj_cpx=tf.signal.ifftshift(pj_cpx, axes=(1,2))
			
			ccf=data_cpx*tf.math.conj(pj_cpx)
			ccf=tf.signal.ifft2d(ccf)
			ccf=tf.signal.ifftshift(ccf, axes=(1,2))
			ccf=tf.math.real(ccf)

			cc=tf.reduce_mean(ccf, axis=0)
			tx=tf.argmax(tf.math.reduce_max(cc, axis=1))
			ty=tf.argmax(tf.math.reduce_max(cc, axis=0))
			trans=tf.convert_to_tensor([tx, ty], dtype=floattype)-ccf.shape[1]/2
			
			tpm_rot[:2]+=trans

		data_cpx=tf.signal.fftshift(data_cpx, axes=(1,2))
		for iv in range(len(vols)):
			vols[iv], wts[iv]=do_insert(tpm_rot, data_cpx[iv], vols[iv], wts[iv]);
			
		alltrans[tid]=trans.numpy()
	return -alltrans

# @tf.function()
def refine_rotation(tltpm_now, all_ts, tpos, niter=50, fullrefine=0):
	pk0=[get_xf_pos(tltpm_now, tp) for tp in tpos]
	pk0=np.array(pk0, dtype=floattype)
	pk0=pk0-all_ts

	dps=[]
	xtlt=np.arange(len(tltpm_now), dtype=floattype)-len(tltpm_now)//2
	xtlt/=np.max(abs(xtlt))

	params=tf.Variable(np.zeros(10), dtype=floattype)
	tpm0=tf.constant(tltpm_now, dtype=floattype)

	opt=tf.keras.optimizers.Adam(learning_rate=0.1) 
	if fullrefine==0:
		pmask=[10,10,0,0,0,  10,10,1,0,0]
	elif fullrefine==1:
		pmask=[10,10,0,0,0,  10,10,1,1,1]
	else:
		pmask=[10,10,1,0,1,  10,10,1,1,1]
	
	for itr in range(niter):
		with tf.GradientTape() as gt:
			param3=xtlt[:,None]*params[None,:5]*pmask[:5]
			param3=param3+params[5:]*pmask[5:]
			# tpm_bias=tf.concat([tf.zeros((len(tltpm_now),2)), param3], axis=1)
			tpm=tpm0+param3

			pk1=[get_xf_pos(tpm, tp) for tp in tpos]
			pk1=tf.convert_to_tensor(pk1, dtype=floattype)
			dp=pk0-pk1

			loss=tf.reduce_mean(dp**2)

			dp=tf.reduce_mean(tf.linalg.norm(dp, axis=2))

		grad=gt.gradient(loss, [params])
		opt.apply_gradients(zip(grad, [params]))
		print(itr, dp, end='\r')
		
	
	print()
	return tpm
    

# @tf.function()
def patch_track(imgs_now, tltpm_now, tpos, imsk):
	
	
	all_ts=[]
	for tp in tpos:
		pk2=get_xf_pos(tltpm_now, tp)+imgs_now.shape[1]/2.-sz/2.
		imgs_tile=tf.image.extract_glimpse(imgs_now[:,:,:,None], size=(sz, sz), offsets=pk2, centered=False, normalized=False, noise='zero')[...,0]
		# print(imgs_tile.shape, imsk.shape)
		imgs_tile*=imsk
		data_cpx=tf.cast(imgs_tile, tf.complex64)*1e-3
		data_cpx=tf.signal.fftshift(data_cpx, axes=(1,2))
		data_cpx=tf.signal.fft2d(data_cpx)
		
		ccf=data_cpx[:-1]*tf.math.conj(data_cpx[1:])
		ccf=tf.signal.ifft2d(ccf)
		ccf=tf.signal.ifftshift(ccf, axes=(1,2))
		ccf=tf.math.real(ccf)

		tx=tf.argmax(tf.math.reduce_max(ccf, axis=2), axis=1)
		ty=tf.argmax(tf.math.reduce_max(ccf, axis=1), axis=1)
		trans=tf.cast(tf.stack([tx, ty], axis=1), floattype)-ccf.shape[1]/2
		ts=tf.concat([[[0,0]], trans], axis=0)
		ts=tf.cumsum(ts, axis=0)
		ts-=ts[len(ts)//2]
		all_ts.append(ts)
		
	all_ts=np.array(all_ts)
	return all_ts

# @tf.function()
def coarse_align(imgs_now):
	data_cpx=tf.cast(imgs_now, tf.complex64)*1e-3
	data_cpx=tf.signal.fftshift(data_cpx, axes=(1,2))
	data_cpx=tf.signal.fft2d(data_cpx)
	ccf=data_cpx[:-1]*tf.math.conj(data_cpx[1:])
	ccf=tf.signal.ifft2d(ccf)
	ccf=tf.signal.ifftshift(ccf, axes=(1,2))
	ccf=tf.math.real(ccf)

	tx=tf.argmax(tf.math.reduce_max(ccf, axis=2), axis=1)
	ty=tf.argmax(tf.math.reduce_max(ccf, axis=1), axis=1)
	trans=tf.cast(tf.stack([tx, ty], axis=1), floattype)-ccf.shape[1]/2
	ts=np.append([[0,0]], trans, axis=0)
	ts=np.cumsum(ts, axis=0)
	ts-=ts[len(ts)//2]
	return ts


#### get 2D position on a tilt given 3D location
@tf.function()
def get_xf_pos(tpm, pk):
	a=tpm[:,2]*np.pi/180.
	b=-tpm[:,3]*np.pi/180.
	c=tpm[:,4]*np.pi/180.

	p1=pk[:,None]

	############
	matrix=tf.stack([
	[tf.ones_like(c), tf.zeros_like(c), tf.zeros_like(c)],
	[tf.zeros_like(c), tf.cos(c), -tf.sin(c)],
	[tf.zeros_like(c), tf.sin(c),  tf.cos(c)],
	], 0)

	matrix=tf.reshape(matrix, shape=[3,3,-1]) 
	matrix=tf.transpose(matrix)
	p1=tf.matmul(matrix, p1)

	############
	matrix=tf.stack([
	[tf.cos(b), tf.zeros_like(b), -tf.sin(b)],
	[tf.zeros_like(b), tf.ones_like(b), tf.zeros_like(b)],
	[tf.sin(b), tf.zeros_like(b),  tf.cos(b)],
	], 0)

	matrix=tf.reshape(matrix, shape=[3,3,-1]) 
	matrix=tf.transpose(matrix)
	p1=tf.matmul(matrix, p1)

	############
	matrix=tf.stack([
	[tf.cos(a), -tf.sin(a), tf.zeros_like(a)],
	[tf.sin(a),  tf.cos(a), tf.zeros_like(a)],
	[tf.zeros_like(a), tf.zeros_like(a), tf.ones_like(a)],
	], 0)

	matrix=tf.reshape(matrix, shape=[3,3,-1]) 
	matrix=tf.transpose(matrix)
	p1=tf.matmul(matrix, p1)
	
	p1=p1[...,0]
	p1=p1[:,:2]+tpm[:,:2]
	return p1


@tf.function()
def make_tile(imgs, pk, tp, tltpm, fullvol):
	
	ts=tf.image.extract_glimpse(imgs[:,:,:,None], size=(sz, sz), offsets=pk, centered=False, normalized=False, noise='zero')[...,0]
	data_cpx=tf.cast(ts, tf.complex64)*1e-3
	data_cpx=tf.signal.fftshift(data_cpx, axes=(1,2))
	data_cpx=tf.signal.fft2d(data_cpx)
	data_cpx=tf.signal.fftshift(data_cpx, axes=(1,2))
	
	# xftile=tf.concat([xfsnp[:,:3], pkres[ii]], axis=1)
	xftile=tltpm*[0,0,1,1,1]
	
	vol=tf.zeros((sz, sz, sz), dtype=np.complex64)
	wt=tf.zeros((sz, sz, sz), dtype=np.float32)
	for i in range(len(tltpm)):
		vol, wt=do_insert(xftile[i], data_cpx[i], vol, wt);
	
	v=tf.math.divide_no_nan(vol, tf.cast(wt, tf.complex64))
	vr=tf.signal.ifftshift(v)
	vr=tf.signal.ifft3d(vr)
	vr=tf.signal.ifftshift(vr)
	vr=tf.math.real(vr)
	
	fullvol=tf.tensor_scatter_nd_add(fullvol, ind_vol+tp, tf.reshape(vr,(-1,)))
	return fullvol
	
@tf.function()
def do_insert(tpm, data, vol, wt):

	a=tpm[2]*np.pi/180.
	b=-tpm[3]*np.pi/180.
	c=tpm[4]*np.pi/180.

	############
	matrix=tf.stack([
	[1, 0, 0],
	[0, tf.cos(c), -tf.sin(c)],
	[0, tf.sin(c),  tf.cos(c)],
	], 0)
	
	matrix=tf.transpose(matrix)
	rotmat=matrix

	############
	matrix=tf.stack([
	[tf.cos(b), 0, -tf.sin(b)],
	[0, 1, 0],
	[tf.sin(b), 0,  tf.cos(b)],
	], 0)

	matrix=tf.transpose(matrix)
	rotmat=tf.matmul( matrix, rotmat)

	############
	matrix=tf.stack([
	[tf.cos(a), -tf.sin(a), 0],
	[tf.sin(a),  tf.cos(a), 0],
	[0,0,1],
	], 0)

	matrix=tf.transpose(matrix)
	rotmat=tf.matmul(matrix,rotmat)
	rotmat=tf.cast(rotmat, floattype)
	# return rotmat
	
	ind_rot=tf.matmul(ind_xy, rotmat)
	ind_rot+=sz//2
	ind_rot=tf.concat([ind_rot[:,:2], abs(ind_rot[:,2:])], axis=1)


	ix=tf.cast(tf.floor(ind_rot), np.int32)
	ir=ind_rot-tf.floor(ind_rot)
	ir=[ir, 1-ir]

	d=data
	x=tf.cast(ind_trans[0]*tpm[0]+ind_trans[1]*tpm[1], np.complex64)
	x=tf.exp(1j*np.pi*2.*x/sz)
	d=d*tf.cast(x, np.complex64)

	ds=tf.reshape(d,(-1,))

	for i0 in ind_interp:
		w=ir[0]*i0+ir[1]*(1-i0)
		w=tf.reduce_prod(w, axis=1)
		w2=tf.cast(w, np.complex64)
		vol=tf.tensor_scatter_nd_add(vol, ix+i0, ds*w2)
		wt=tf.tensor_scatter_nd_add(wt, ix+i0, w)

	return vol, wt


@tf.function()
def do_project(tpm, vol):

    a=tpm[2]*np.pi/180.
    b=-tpm[3]*np.pi/180.
    c=tpm[4]*np.pi/180.

    ############
    matrix=tf.stack([
    [1, 0, 0],
    [0, tf.cos(c), -tf.sin(c)],
    [0, tf.sin(c),  tf.cos(c)],
    ], 0)
    
    matrix=tf.transpose(matrix)
    rotmat=matrix

    ############
    matrix=tf.stack([
    [tf.cos(b), 0, -tf.sin(b)],
    [0, 1, 0],
    [tf.sin(b), 0,  tf.cos(b)],
    ], 0)

    matrix=tf.transpose(matrix)
    rotmat=tf.matmul(matrix, rotmat)

    ############
    matrix=tf.stack([
    [tf.cos(a), -tf.sin(a), 0],
    [tf.sin(a),  tf.cos(a), 0],
    [0,0,1],
    ], 0)

    matrix=tf.transpose(matrix)
    rotmat=tf.matmul(matrix, rotmat)
    rotmat=tf.cast(rotmat, floattype)
    # print(rotmat)
    
    ind_rot=tf.matmul(ind_xy, rotmat)
    ind_rot+=sz//2
    ind_rot=tf.concat([ind_rot[:,:2], abs(ind_rot[:,2:])], axis=1)


    ix=tf.cast(tf.floor(ind_rot), np.int32)
    ir=ind_rot-tf.floor(ind_rot)
    ir=[ir, 1-ir]

    x=tf.cast(ind_trans[0]*tpm[0]+ind_trans[1]*tpm[1], np.complex64)
    x=tf.exp(1j*np.pi*2.*x/sz)
    # d=d*tf.cast(x, np.complex64)

    # ds=tf.reshape(d,(-1,))
    pjs=[]
    for i0 in ind_interp:
        w=ir[0]*i0+ir[1]*(1-i0)
        w=tf.reduce_prod(w, axis=1)
        w2=tf.cast(w, np.complex64)
        pj=tf.gather_nd(vol, ix+i0)*w2
        pjs.append(pj)
        
    pj=tf.reduce_sum(pjs, axis=0)
    pj=tf.reshape(pj,(sz,sz))
    return pj


if __name__ == '__main__':
	main()
