#!/usr/bin/env python
from EMAN2 import *
import numpy as np

if "CUDA_VISIBLE_DEVICES" not in os.environ:
    # so we can decide which gpu to use with environmental variable
    os.environ["CUDA_VISIBLE_DEVICES"]='1' 

#### do not occupy the entire GPU memory at once
##   seems necessary to avoid some errors...
os.environ["TF_FORCE_GPU_ALLOW_GROWTH"]='true' 
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' #### reduce log output

#### finally initialize tensorflow
import tensorflow as tf
floattype=np.float32
def main():
	parser = EMArgumentParser(usage="")
	parser.add_argument("--input", default=None, help="The input projections. Project should usually have the xform.projection header attribute, which is used for slice insertion")
	parser.add_argument("--output", default=None, help="Output reconstructed volume file name.")
	parser.add_argument("--sym", default="c1", help="Symmetry.")	
	parser.add_argument("--batch", type=int, help="batch size",default=2000)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	logger=E2init(sys.argv,options.ppid)
	time0=time.time()


	global ind_xy, ind_interp, ind_trans, rawbox
	
	##################
	e=EMData(options.input, 0, True)
	rawbox=e["ny"]
	
	##################
	## index for rotation
	ind_xy=np.indices((rawbox, rawbox, 1))
	ind_xy=ind_xy.transpose((1,2,3,0)).reshape((-1,3))
	ind_xy=ind_xy.astype(floattype)
	ind_xy[:,:2]-=rawbox//2
	
	##################
	## index for interpolation
	ind_interp=np.indices((2,2,2))
	ind_interp=ind_interp.transpose((1,2,3,0)).reshape(-1,3)

	##################
	## index for translation
	ind_trans=np.indices((rawbox, rawbox))
	ind_trans=ind_trans.astype(np.float32)-rawbox/2
		
	##################
	vol=tf.zeros((rawbox, rawbox, rawbox), dtype=np.complex64)
	wt=tf.zeros((rawbox, rawbox, rawbox), dtype=np.float32)
	
	
	nptcl=EMUtil.get_image_count(options.input)
	for ip in range(0, nptcl, options.batch):
		projs=[]
		hdrs=[]
		el=EMData.read_images(options.input,range(ip,min(ip+options.batch,nptcl)))
		for e in el:
			hdrs.append(e.get_attr_dict())
			p=e.numpy().T.copy()*1e-3
			projs.append(p)

		projs=np.array(projs, dtype=np.complex64)
		projs=tf.signal.fftshift(projs, axes=(1,2))
		data_cpx=tf.signal.fft2d(projs)
		data_cpx=tf.signal.fftshift(data_cpx, axes=(1,2))
		# print(data_cpx.shape)
		
		##################
		xfsnp=[]
		nsym=Transform.get_nsym(options.sym)
			
		for h in hdrs:
			x=h["xform.projection"]
			xs=[x.get_sym(options.sym, s) for s in range(nsym)]
			xs=[p.get_params("eman") for p in xs]		
			xnp=np.array([[x["az"],x["alt"],x["phi"], x["tx"], x["ty"]] for x in xs], dtype=floattype)
			xfsnp.append(xnp)
			
		xfsnp=np.array(xfsnp, dtype=floattype)
		xfsnp[:,:,:3]=xfsnp[:,:,:3]*np.pi/180.
		xfsnp=tf.constant(xfsnp, dtype=floattype)

		
		##################
		
		for ii in range(len(xfsnp)):
			vol, wt=do_insert(xfsnp[ii], data_cpx[ii], vol, wt)
			if ii%200==0: print(ip+ii, time.time()-t0, end='\r')

	v=vol.numpy()
	w=wt.numpy()
	w[w==0]=1
	v/=w

	v=np.fft.ifftshift(v)
	r=np.fft.ifftshift(np.fft.ifftn(v)).real

	r=r.transpose(2,1,0).copy()

	b=from_numpy(r.copy())
	b.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.5})
	b.process_inplace("normalize.edgemean")
	b.write_image(options.output)


	E2end(logger)

	print("\nReconstruction finishend ({:.1f} s)".format(time.time()-time0))

@tf.function()
def do_insert(angs, data, vol, wt):
	
	for ang in angs:
		azp=-ang[2]
		altp=-ang[1]
		phip=-ang[0]

		matrix=tf.stack([(tf.cos(phip)*tf.cos(azp) - tf.cos(altp)*tf.sin(azp)*tf.sin(phip)),
		(tf.cos(phip)*tf.sin(azp) + tf.cos(altp)*tf.cos(azp)*tf.sin(phip)),
		(tf.sin(altp)*tf.sin(phip)),

		(-tf.sin(phip)*tf.cos(azp) - tf.cos(altp)*tf.sin(azp)*tf.cos(phip)),
		(-tf.sin(phip)*tf.sin(azp) + tf.cos(altp)*tf.cos(azp)*tf.cos(phip)),
		(tf.sin(altp)*tf.cos(phip)),

		(tf.sin(altp)*tf.sin(azp)),
		(-tf.sin(altp)*tf.cos(azp)),
		tf.cos(altp)], 0)

		matrix=tf.reshape(matrix, shape=[3,3]) 
		matrix=tf.transpose(matrix)

		ind_rot=tf.matmul(ind_xy, matrix)
		ind_rot+=rawbox//2
		ind_rot=tf.concat([ind_rot[:,:2], abs(ind_rot[:,2:])], axis=1)


		ix=tf.cast(tf.floor(ind_rot), np.int32)
		ir=ind_rot-tf.floor(ind_rot)
		ir=[ir, 1-ir]
		
		d=data
		x=tf.cast(ind_trans[0]*ang[3]+ind_trans[1]*ang[4], np.complex64)
		x=tf.exp(1j*np.pi*2.*x/rawbox)
		d=d*tf.cast(x, np.complex64)

		ds=tf.reshape(d,(-1,))

		for i0 in ind_interp:
			w=ir[0]*i0+ir[1]*(1-i0)
			w=tf.reduce_prod(w, axis=1)
			w2=tf.cast(w, np.complex64)
			vol=tf.tensor_scatter_nd_add(vol, ix+i0, ds*w2)
			wt=tf.tensor_scatter_nd_add(wt, ix+i0, w)

	return vol, wt
	
	
if __name__=="__main__":
	main()
