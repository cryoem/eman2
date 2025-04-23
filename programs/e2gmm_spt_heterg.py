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
	parser.add_argument("--mask", type=str,help="masks for refinement.", default=None)
	parser.add_argument("--niter", type=int, help="number of iteration",default=20)
	parser.add_argument("--clip", type=int, help="clip image to size",default=-1)
	
	# parser.add_argument("--midin", type=str,help="tranform input", default=None)
	parser.add_argument("--midout", type=str,help="tranform output", default=None)

	parser.add_argument("--encoderin", type=str,help="encoder input", default=None)
	parser.add_argument("--encoderout", type=str,help="encoder output", default=None)
	parser.add_argument("--decoderin", type=str,help="decoder input", default=None)
	parser.add_argument("--decoderout", type=str,help="decoder output", default=None)

	parser.add_argument("--anchor", type=str,help="anchor points save file. will generate from model by default", default=None)
	parser.add_argument("--n_anchor", type=int,help="number of anchor points. default 32", default=32)
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

	ln=np.mean([len(i) for i in idx3d])
	print("{} 3D particles, each contain {:.1f} 2D particles".format(len(idx3d), ln))
		
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
	
	path=os.path.dirname(options.midout)
	#### anchor points
	if options.anchor==None:
		options.anchor=f"{path}/model_00_anchor.txt"

	if os.path.isfile(options.anchor):
		anchor=np.loadtxt(options.anchor)
		anchor=anchor.astype(floattype).copy()
		print(f"load {len(anchor)} anchor points from {options.anchor}")

	else:
		
		pn=options.n_anchor//2
		km=KMeans(pn,max_iter=30)
		km.fit(pts[:,:3])
		pc=km.cluster_centers_

		pm=pts[imsk>.1]
		pn=options.n_anchor//2
		km=KMeans(pn,max_iter=30)
		km.fit(pm[:,:3])
		pc2=km.cluster_centers_

		pcx=np.vstack([pc, pc2])
		pp=np.hstack([pcx, np.zeros((len(pcx),1))+np.mean(pts[:,3]), np.zeros((len(pcx),1))+np.mean(pts[:,4])])

		anchor=pp[:,:3].astype(floattype).copy()
		np.savetxt(options.anchor, anchor)
		print(f"generated {len(anchor)} anchor points, saved to {options.anchor}")

	######## build encoder and decoder
	pts=tf.constant(pts[None,...])
	print(pts.shape, anchor.shape)
	if options.decoderin==None:
		decode_model=build_decoder_anchor(pts, anchor, ninp=4)
	else:
		decode_model=tf.keras.models.load_model(options.decoderin,compile=False)

	if options.encoderin==None:
		encode_model=build_encoder(nout=4, conv=False,ninp=len(pts[0]))
	else:
		encode_model=tf.keras.models.load_model(options.encoderin,compile=False)

	batchsz=4

	####### pre-compute gradients
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
		#print(pn, fval.shape, grad.shape)
		#scr=tf.reshape(fval, (len(pids), -1))
		#scr=tf.reduce_mean(scr, axis=1)
		scr=np.array([tf.reduce_sum(fval[x1-x0:x1]) for x0,x1 in zip(pn,np.cumsum(pn))])
		

		allgrds.append(grad)
		allscr.append(scr)
		print(pii, len(p3did), np.mean(scr), end='\r')

	allgrds=np.concatenate(allgrds, axis=0)
	allscr=np.concatenate(allscr, axis=0)
	allgrds=allgrds/np.std(allgrds)
	allgrds=allgrds.reshape((len(allgrds), -1)).astype(floattype)
	print(allgrds.shape)

	#### this initialize the weights in model
	conf=encode_model(allgrds[:2], training=True)
	pout=decode_model(conf, training=True)

	###############
	##### now train

	pas=[1,0,0]
	pas=tf.constant(np.array([pas[0],pas[0],pas[0],pas[1],pas[2]], dtype=floattype))

	wts=encode_model.trainable_variables + decode_model.trainable_variables[:-1]
	opt=tf.keras.optimizers.Adam(learning_rate=options.learnrate)

	for itr in range(options.niter):
		cost=[]
		for pii in range(0, len(p3did), batchsz):
			pids=p3did[pii:pii+batchsz]
			pid=np.concatenate(pids)
			pn=np.array([len(i) for i in pids])

			ptr=tf.gather(data_cpx[0], pid)
			ptj=tf.gather(data_cpx[1], pid)
			ptcl_cpx=(ptr, ptj)

			xf=xfsnp[pid]
			grd=allgrds[pii:pii+batchsz]

			with tf.GradientTape() as gt:

				conf=encode_model(grd, training=True)

				cl=tf.math.sqrt(tf.nn.relu(tf.reduce_sum(conf**2, axis=1)))
				cl=tf.reduce_mean(tf.maximum(cl-1,0))

				conf=0.01*tf.random.normal(conf.shape)+conf
				pout=decode_model(conf, training=True)

				pout=pout*imsk[None,:,None]
				pout=pout*pas
				pout+=pts

				pt2=tf.repeat(pout, pn, axis=0)

				imgs_cpx=pts2img(pt2, xf)
				fval=calc_frc(ptcl_cpx, imgs_cpx, params["rings"], minpx=options.minpx, maxpx=options.maxpx)

				loss=-tf.reduce_mean(fval)+cl*10

			grad=gt.gradient(loss, wts)
			opt.apply_gradients(zip(grad, wts))
			cost.append(loss)
			print(itr, pii, len(p3did), float(loss), end='\r')

		print(itr, np.mean(cost), "                       ")

	if options.niter>0:

		if options.decoderout!=None:
			if os.path.isfile(options.decoderout):
				os.remove(options.decoderout)
			decode_model.save(options.decoderout)
			print("Decoder saved as ",options.decoderout)
		if options.encoderout!=None:
			if os.path.isfile(options.encoderout):
				os.remove(options.encoderout)
			encode_model.save(options.encoderout)
			print("Encoder saved as ",options.encoderout)

	###############
	##### get latent output
	allconf=[]
	batchsz*=8
	for pii in range(0, len(p3did), batchsz):
		grd=allgrds[pii:pii+batchsz]
		conf=encode_model(grd, training=False)
		allconf.append(conf)

	allconf=np.concatenate(allconf, axis=0)
	allconf=np.hstack([np.arange(len(allconf))[:,None], allconf])

	print("conformation output", allconf.shape)
	np.savetxt(options.midout, allconf)

	E2end(logid)
	return
	
	
if __name__ == '__main__':
	main()
	
