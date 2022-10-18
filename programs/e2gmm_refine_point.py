#!/usr/bin/env python
# Muyuan Chen 2020-06
from EMAN2 import *
import numpy as np
from sklearn.decomposition import PCA
from EMAN2_utils import pdb2numpy

#### need to unify the float type across tenforflow and numpy
##   in theory float16 also works but it can be unsafe especially when the network is deeper...
floattype=np.float32

#### here we import tensorflow at the global scale so @tf.function works
##   although how much performance gain we get from @tf.function is questionable...
##   so here we need to make a fake tf module for --help so the CI works properly.
os.environ["CUDA_VISIBLE_DEVICES"]='0' 
os.environ["TF_FORCE_GPU_ALLOW_GROWTH"]='true' 
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' #### reduce log output
if ('-h' in sys.argv) or ('--help' in sys.argv):
	tf=type('empty', (object,), {})()
	tf.function=lambda f: f
	print("Printing help. Skip tensorflow import")
else:
	import tensorflow as tf

def main():

	usage="""Single particle alignment, reconstruction and heterogeneity analysis with Gaussian model and neural networks. There are a few modes.

	Reconstruction from projections or particles with transforms:
	e2gmm_refine.py --projs <projection file> --modelout <model output text file> --npts <number of Gaussian in model>

	Align particles using the model (local refinement only, still under testing):
	e2gmm_refine.py --model <model input> --ptclsin <particle input list> --ptclout <particle list output> --align

	Heterogeneity analysis:
	e2gmm_refine.py --model <model input> --ptclsin <particle input list>  --heter --midout <conformation output text file>

	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--chunk", type=str,help="In case of memory exhaustion, particles can be trained in chunks (similar to a batch). eg - 0,10 will process the first of 10 chunks. Be sure to store and restore the encoder and decoder if using this mode. Must be run sequentially 0-(n-1) without parallelism.", default=None)
	parser.add_argument("--sym", type=str,help="symmetry. currently only support c and d", default="c1")
	parser.add_argument("--model", type=str,help="load from an existing model file", default="")
	parser.add_argument("--segments", type=str,help="Divide the model into sequential domains. Comma separated list of integers. Each integer is the first sequence number of a new region, starting with 0",default=None)
	parser.add_argument("--modelout", type=str,help="output trained model file. only used when --projs is provided", default="")
	parser.add_argument("--decoderin", type=str,help="Rather than initializing the decoder from a model, read an existing trained decoder", default="")
	parser.add_argument("--decoderout", type=str,help="Save the trained decoder model. Filename should be .h5", default=None)
	parser.add_argument("--encoderin", type=str,help="Rather than initializing the encoder from scratch, read an existing trained encoder", default="")
	parser.add_argument("--encoderout", type=str,help="Save the trained encoder model. Filename should be .h5", default=None)
	parser.add_argument("--projs", type=str,help="projections with orientations (in hdf header or comment column of lst file) to train model", default="")
	parser.add_argument("--evalmodel", type=str,help="generate model projection images to the given file name", default="")
	parser.add_argument("--evalsize", type=int,help="Box size for the projections for evaluation.", default=-1)
	parser.add_argument("--ptclsin", type=str,help="particles input for alignment", default="")
	parser.add_argument("--ptclsout", type=str,help="aligned particle output", default="")
	parser.add_argument("--ptclsclip",type=int,help="clip particles to specified box size before use",default=-1)
	parser.add_argument("--learnrate", type=float,help="learning rate for model training only. Default is 1e-4. ", default=1e-4)
#	parser.add_argument("--sigmareg", type=float,help="regularizer for the sigma of gaussian width. Larger value means all Gaussian functions will have essentially the same width. Smaller value may help compensating local resolution difference.", default=.5)
	parser.add_argument("--modelreg", type=float,help="regularizer for for Gaussian positions based on the starting model, ie the result will be biased towards the starting model when training the decoder (0-1 typ). Default 0", default=0)
	parser.add_argument("--ampreg", type=float,help="regularizer for the Gaussian amplitudes in the first 1/2 of the iterations. Large values will encourage all Gaussians to have similar amplitudes. default = 0", default=0)
	parser.add_argument("--niter", type=int,help="number of iterations", default=10)
	parser.add_argument("--npts", type=int,help="number of points to initialize. ", default=-1)
	parser.add_argument("--batchsz", type=int,help="batch size", default=256)
	parser.add_argument("--minressz", type=int,help="Fourier diameter associated with minimum resolution to consider. ", default=4)
	parser.add_argument("--maxboxsz", type=int,help="maximum fourier box size to use. 2 x target Fourier radius. ", default=64)
	parser.add_argument("--maxres", type=float,help="maximum resolution. will overwrite maxboxsz. ", default=-1)
	parser.add_argument("--align", action="store_true", default=False ,help="align particles.")
	parser.add_argument("--heter", action="store_true", default=False ,help="heterogeneity analysis.")
	parser.add_argument("--decoderentropy", action="store_true", default=False ,help="This will train some entropy into the decoder using particles to reduce vanishing gradient problems")
	parser.add_argument("--perturb", type=float, default=0.1 ,help="Relative perturbation level to apply in each iteration during --heter training. Default = 0.1, decrease if models are too disordered")
	parser.add_argument("--conv", action="store_true", default=False ,help="Use a convolutional network for heterogeneity analysis.")
	parser.add_argument("--fromscratch", action="store_true", default=False ,help="start from coarse alignment. otherwise will only do refinement from last round")
	parser.add_argument("--ptclrepout", type=str,help="Save the per-particle representation input to the network to a file for later use", default="")
	parser.add_argument("--ptclrepin", type=str,help="Load the per-particle representation from a file rather than recomputing", default="")
	parser.add_argument("--midout", type=str,help="middle layer output", default="")
	parser.add_argument("--pas", type=str,help="choose whether to adjust position, amplitude, sigma. sigma is not supported in this version of the program. use 3 digit 0/1 input. default is 110, i.e. only adjusting position and amplitude", default="110")
	parser.add_argument("--nmid", type=int,help="size of the middle layer", default=4)
	parser.add_argument("--mask", type=str,help="remove points outside mask", default="")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	if options.chunk is not None:
		try: options.chunk=(int(options.chunk.split(",")[0]),int(options.chunk.split(",")[1]))
		except:
			print("Error: invalid chunk specification")
			sys.exit(1)
		if options.chunk[0]<0 or options.chunk[0]>options.chunk[1]-1 :
			print("Error: chunk specification out of range")
			sys.exit(1)

	logid=E2init(sys.argv,options.ppid)

	gen_model=None
	maxboxsz=good_size(options.maxboxsz)

	if options.pas[2]=="1" :
		print("ERROR: Sigma is not supported in this version of the program")
		sys.exit(1)

	## load GMM from text file
	if options.model:
		if options.model.endswith(".pdb"):

			p=pdb2numpy(options.model)
			pts=np.zeros((len(p),4))
			e=EMData(options.projs, 0, True)
			#sz=e["ny"]
			p=p/e["ny"]/e["apix_x"]-0.5
			pts[:,:3]=p
			pts[:,3]=.5
			#pts[:,4]=1

			print(pts)
		else:

			pts=np.loadtxt(options.model).astype(floattype)
			if pts.shape[1]==5: pts=pts[:,:4]
		npt=len(pts)
		print("{} gaussian in the model".format(len(pts)))

		if options.segments!=None:
			segments=[int(s) for s in options.segments.split(",")]
			if segments[0]!=0 : segments.insert(0,0)
			for i in range(len(segments-1)):
				if segments[i-1]>=segments[i] or segments[i]>npt:
					print("ERROR: segments must be a sequential list of integer model point numbers")
					sys.exit(1)
			options.nmid=len(segments)*2	# 2-D space for each feature
		else: segments=None
	else:
		pts=None

	#### This initializes the decoder directly from a set of coordinates
	#    This method (may be) used rather than saving the decoder model itself
	if (not options.decoderin) and options.model and options.projs:
		print("Recompute model from coordinates...")

		## duplicates points if we ask for more points than exist in text file
		if options.npts>len(pts):
			p=pts.copy()
			np.random.shuffle(p)
			p=np.repeat(p, 8, axis=0)
			p=p[:(options.npts-len(pts))]
			pts=np.concatenate([pts, p], axis=0)

		## randomize it a bit so we dont have all zero weights
		rnd=np.random.randn(pts.shape[0], pts.shape[1])*1e-3
		gen_model=build_decoder(options.nmid,pts+rnd,segments)
		print("{} gaussian in the model".format(len(pts)))

		## train the model from coordinates first
		conf=tf.zeros((1,options.nmid), dtype=floattype)
		#opt=tf.keras.optimizers.Adam(learning_rate=1e-4)
		opt=tf.keras.optimizers.Adamax()
		gen_model.compile(optimizer=opt, loss=tf.losses.MeanAbsoluteError())
		loss=[]
		for i in range(500):
			loss.append(gen_model.train_on_batch(conf, pts.reshape(1,pts.shape[0],pts.shape[1])))

		print("Abs loss from loaded model : {:.05f}".format(loss[-1]))

	# Read the complete decoder rather than reinitializing from the model
	if options.decoderin:
		gen_model=tf.keras.models.load_model(f"{options.decoderin}",compile=False)

	#### Decoder training from generated projections of a 3-D map or particles
	#### Note that train_decoder takes options.decoderentropy into account internally
	if options.projs:
		# The shape of the decoder is defined by the number of Gaussians (npts) and the number of latent variables (nmid)
		if gen_model==None:
			gen_model=build_decoder(options.nmid,options.npts)
		print("Train model from ptcl-xfrom pairs...")
		e=EMData(options.projs, 0, True)
		raw_apix, raw_boxsz = e["apix_x"], e["ny"]
		options.raw_apix=raw_apix
		if options.maxres>0:
#			maxboxsz=options.maxboxsz=ceil(raw_boxsz*raw_apix*2/options.maxres)//2*2
			maxboxsz=options.maxboxsz=good_size(ceil(raw_boxsz*raw_apix*2/options.maxres)//2*2)
			print("using box size {}, max resolution {:.1f}".format(maxboxsz, options.maxres))

		data_cpx, xfsnp, grpdct = load_particles(options.projs, maxboxsz, preclip=options.ptclsclip,chunk=options.chunk)
		apix=raw_apix*raw_boxsz/maxboxsz
		clipid=set_indices_boxsz(data_cpx[0].shape[1], apix, True)

		## training
		if options.niter>0:
			params=set_indices_boxsz(maxboxsz)
			dcpx=get_clip(data_cpx, params["sz"], clipid)
			trainset=tf.data.Dataset.from_tensor_slices((dcpx[0], dcpx[1], xfsnp))
			trainset=trainset.batch(options.batchsz)

			train_decoder(gen_model, trainset, params, options, pts)
			pout=gen_model(tf.zeros((1,options.nmid), dtype=floattype)).numpy()[0]

			pout[:,3]/=np.max(pout[:,3])

			if options.mask:
				## remove Gaussian points that falls outside the mask

				msk=EMData(options.mask)
				m=msk.numpy().copy()

				p=pout[:,:3].copy()
				p=p[:,::-1]
				p[:,:2]*=-1
				p=(p+.5)*msk["nx"]

				o=np.round(p).astype(int)
				v=m[o[:,0], o[:,1], o[:,2]]
				pout=pout[v>.9]

			#### save decoder if requested
			if options.decoderout!=None:
				gen_model.save(options.decoderout)
				print("Decoder saved as ",options.decoderout)

			#### save model to text file
			np.savetxt(options.modelout, pout)
			gen_model=build_decoder(options.nmid, pout,segments)


		#### make projection images from GMM
		if options.evalmodel:
			if options.evalsize<0:
				options.evalsize=raw_boxsz
			eval_model(gen_model, options)


	#### Load particles with xforms in header
	if options.ptclsin:
		e=EMData(options.ptclsin, 0)	# used to be header only, but inconvenient with clip possibility
		if options.ptclsclip>0 : e=e.get_clip(Region((e["nx"]-options.ptclsclip)//2,(e["ny"]-options.ptclsclip)//2,options.ptclsclip,options.ptclsclip))
		raw_apix, raw_boxsz = e["apix_x"], e["ny"]
		options.raw_apix=raw_apix
		if options.maxres>0:
			maxboxsz=options.maxboxsz=good_size(ceil(raw_boxsz*raw_apix*2/options.maxres)//2*2)
			print("using box size {}, max resolution {:.1f}".format(maxboxsz, options.maxres))

		data_cpx, xfsnp, grpdct = load_particles(options.ptclsin,maxboxsz,preclip=options.ptclsclip,chunk=options.chunk)
		apix=raw_apix*raw_boxsz/maxboxsz
		clipid=set_indices_boxsz(data_cpx[0].shape[1], apix, True)

	#### Align particles using GMM
	##   have not upgraded this part yet. probably still bugs left
	if options.ptclsin and options.align:
		pts=tf.constant(pts)
		print("Align particles...")
		if options.fromscratch:
			set_indices_boxsz(32)
			dcpx=get_clip(data_cpx, sz)
			xfs=coarse_align(dcpx, pts, options)

			xfsnp=np.zeros((data_cpx[0].shape[0], 4), dtype=floattype)
			frcs=np.zeros(data_cpx[0].shape[0], dtype=floattype)
			for xf in xfs:
				set_indices_boxsz(32)
				dcpx=get_clip(data_cpx, sz)
				xo, fc=refine_align(dcpx, xf, pts, options, lr=1e-2)
				fid=fc>frcs
				xfsnp[fid]=xo[fid]
				frcs[fid]=fc[fid]

			set_indices_boxsz(64)
			dcpx=get_clip(data_cpx, sz)
			xfsnp, frcs=refine_align(dcpx, xfsnp, pts, options)

		set_indices_boxsz(maxboxsz)
		dcpx=get_clip(data_cpx, sz)
		xfsnp, frcs=refine_align(dcpx, xfsnp, pts, options, lr=1e-4)

		save_ptcls_xform(xfsnp, raw_boxsz, options, frcs)

	#### Heterogeneity analysis from particles
	bsz=options.batchsz
	if options.ptclsin:
		pts=tf.constant(pts[None,:,:])
		params=set_indices_boxsz(maxboxsz)
		dcpx=get_clip(data_cpx, params["sz"], clipid)

		#### calculate d(FRC)/d(GMM) for each particle
		##   this will be the input for the deep network in place of the particle images
		if options.ptclrepin:
			## optionally load from saved files
			hdr=EMData(options.ptclrepin,0,1)
			nptcl=hdr["ny"]
			nvec=hdr["nx"]
			if options.chunk is not None:
				chunkn=nptcl//options.chunk[1]
				chunk0=chunkn*options.chunk[0]
				chunkn=min(nptcl-chunk0,chunkn)
				ag=EMData(options.ptclrepin,0,0,Region(0,chunk0,nvec,chunkn))
			else: ag=EMData(options.ptclrepin)

			allgrds=ag.numpy().copy()
			del ag
			allscr=allgrds[:,0]
			allgrds=allgrds[:,1:].reshape((len(allgrds), npt))
			print("Ptcl rep shape: ", allgrds.shape)

		else:
			trainset=tf.data.Dataset.from_tensor_slices((dcpx[0], dcpx[1], xfsnp))
			trainset=trainset.batch(bsz)
			allscr, allgrds=calc_gradient(trainset, pts, params, options )
#			allscr, allgrds=calc_gqual(trainset, pts, params, options )


			#### For tomographic data we sum the gradients over a 3-D particle, so all 2-D tilts have the same gradient
			if not grpdct is None:
				for k in grpdct.keys():
					# replaces the gradient for each 2D particle with the average gradient over all tilts in a 3D particle
					allgrds[grpdct[k]]=np.add.reduce(allgrds[grpdct[k]])*(10.0/len(grpdct[k]))	# 10 aribtrary value for test

			#### New idea. Keep only the amplitude gradient, gradient shape becomes (nptcl,ngauss) rather than (nptcl,ngauss,4)
			allgrds=allgrds[:,:,3]		# not necessary if using calc_gqual

			## save to hdf file
			if options.ptclrepout:
				#allgrds=allgrds.reshape((len(allgrds),-1))
				print("Ptcl rep shape: ", allgrds.shape)
				ag=from_numpy(np.hstack([allscr[:,None], allgrds]))
				ag.write_image(options.ptclrepout)
				del ag
				#allgrds=allgrds.reshape((len(allgrds), npt))

	if options.ptclsin and options.heter:
		#### build deep networks and make sure they work
		if options.encoderin:
			encode_model=tf.keras.models.load_model(f"{options.encoderin}",compile=False)
		else:
			encode_model=build_encoder(len(pts[0]),options.nmid)

		if options.decoderin:
			decode_model=tf.keras.models.load_model(f"{options.decoderin}",compile=False)
		else:
			decode_model=build_decoder(options.nmid, len(pts[0]), segments)

		mid=encode_model(allgrds[:bsz])
		print("Latent space shape: ", mid.shape)
		out=decode_model(mid)
		print("Output shape: ",out.shape)
		print("Deviation from neutral model: ", np.mean(abs(out-pts)))

		#### actual training
		ptclidx=allscr>-1
		trainset=tf.data.Dataset.from_tensor_slices((allgrds[ptclidx], dcpx[0][ptclidx], dcpx[1][ptclidx], xfsnp[ptclidx]))
		trainset=trainset.batch(bsz)

		train_heterg(trainset, pts, encode_model, decode_model, params, options)

		if options.decoderout!=None:
			decode_model.save(options.decoderout)
			print("Decoder saved as ",options.decoderout)

		if options.encoderout!=None:
			encode_model.save(options.encoderout)
			print("Encoder saved as ",options.encoderout)

		## conformation output
		mid=calc_conf(encode_model, allgrds[ptclidx], 1000)

		if options.midout:
			if options.chunk is not None:
				if options.chunk[0]==0: out=open(options.midout,"w")	# first batch erases the file
				else: out=open(options.midout,"a")	# later batches append

				for i in range(mid.shape[0]):
					out.write(f"{mid[i][0]:1.6f}\t{mid[i][1]:1.6f}\t{mid[i][2]:1.6f}\t{mid[i][3]:1.6f}\n")
				out.close()
			#sv=np.hstack([np.where(ptclidx)[0][:,None], mid])
			#print(mid.shape, sv.shape)
			#else: np.savetxt(options.midout, sv)

			print("Conformation output saved to {}".format(options.midout))

	if options.ptclrepin and options.encoderin and options.midout:
		hdr=EMData(options.ptclrepin,0,1)
		nptcl=hdr["ny"]
		nvec=hdr["nx"]
		if options.chunk is not None:
			chunkn=nptcl//options.chunk[1]
			chunk0=chunkn*options.chunk[0]
			chunkn=min(nptcl-chunk0,chunkn)
			ag=EMData(options.ptclrepin,0,0,Region(0,chunk0,nvec,chunkn))
		else: ag=EMData(options.ptclrepin)

		allgrds=ag.numpy().copy()
		del ag
		allscr=allgrds[:,0]
		print(allgrds.shape,allgrds[:,1:].shape,npt)
		allgrds=allgrds[:,1:].reshape((len(allgrds), npt))
		print("Ptcl rep shape: ", allgrds.shape)

		encode_model=tf.keras.models.load_model(f"{options.encoderin}",compile=False)

		mid=calc_conf(encode_model, allgrds, 1000)

		out=open(options.midout,"w")	# later batches append

		for i in range(mid.shape[0]):
			rep="\t".join([f"{v:1.6f}" for v in mid[i]])
			out.write(f"{i:d}\t{rep}\n")
		out.close()

	E2end(logid)


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

#### Project 3d Gaussian coordinates based on transforms to make projection
##   input:  pts - ( batch size, number of Gaussian, 3 (x,y,z) )
##                 ( number of Gaussian, 3) should also work
##           ang - ( batch size, 5 (az, alt, phi, tx, ty) )
#@tf.function
def xf2pts(pts, ang):

	#### input EMAN style euler angle (az, alt, phi) and make projection matrix
	##   note we need to be able to deal with a batch of particles at once
	##   so everything is in matrix form
	azp=-ang[:,0]
	altp=ang[:,1]
	phip=-ang[:,2]

	matrix=tf.stack([(tf.cos(phip)*tf.cos(azp) - tf.cos(altp)*tf.sin(azp)*tf.sin(phip)),
	(tf.cos(phip)*tf.sin(azp) + tf.cos(altp)*tf.cos(azp)*tf.sin(phip)),
	(tf.sin(altp)*tf.sin(phip)),

	(-tf.sin(phip)*tf.cos(azp) - tf.cos(altp)*tf.sin(azp)*tf.cos(phip)),
	(-tf.sin(phip)*tf.sin(azp) + tf.cos(altp)*tf.cos(azp)*tf.cos(phip)),
	(tf.sin(altp)*tf.cos(phip)),

	(tf.sin(altp)*tf.sin(azp)),
	(-tf.sin(altp)*tf.cos(azp)),
	tf.cos(altp)], 0)

	matrix=tf.transpose(matrix)
	matrix=tf.reshape(matrix, shape=[-1, 3,3]) #### Here we get a batch_size x 3 x 3 matrix

	#### rotate Gaussian positions
	##   here we try to make it also work when pts contains only the neutral model
	if len(pts.shape)>2:
		pts_rot=tf.tensordot(pts, matrix, [[2],[2]])
		pts_rot=tf.transpose(pts_rot, (0,2,1,3))

		#### the eye matrix here is mathematically unnecessary
		##   but somehow tensorflow 2.0 does not track gradient properly without it...
		##   shouldn't do much damage on the performance anyway
		e=tf.eye(pts.shape[0], dtype=bool)#.flatten()
		pts_rot=pts_rot[e]
		
	else:
		pts_rot=tf.tensordot(pts, matrix, [[1],[2]])
		pts_rot=tf.transpose(pts_rot, [1,0,2])

	#### finally do the translation
	tx=ang[:,3][:,None]
	ty=ang[:,4][:,None]
#	pts_rot_trans=tf.stack([(pts_rot[:,:,0]+tx), (-pts_rot[:,:,1])+ty], 2)
	pts_rot_trans=tf.stack([(-pts_rot[:,:,1])+ty,(pts_rot[:,:,0]+tx)], 2)

	#pts_rot_trans=pts_rot_trans*sz+sz/2
	return pts_rot_trans


#### make 2D projections from Gaussian coordinates in Fourier space
##   input:  pts - ( batch size, number of Gaussian, 5 (x,y,z,amp,sigma) )
##                 ( number of Gaussian, 3) should also work
##           ang - ( batch size, 5 (az, alt, phi, tx, ty) )
##        params - a dictionary of some Fourier indices for slicing
##                 sz - Fourier box size
##                 idxft - Fourier indices
##                 rrft - radial Fourier indices
##            lp - lowpass filter applied to the images
##                 this should not be necessary since we use FRC for loss
##                 but the dynamic range of values in Fourier space can sometimes be too high...
##           sym - symmetry string
#@tf.function
def pts2img(pts, ang, params, lp=.1, sym="c1"):
	bsz=ang.shape[0]
	sz, idxft, rrft=params["sz"], params["idxft"], params["rrft"]
	xfo=params["xforigin"]
	
	### initialize output and parse input
	imgs=tf.zeros((bsz, sz,sz), dtype=floattype)
	if len(pts.shape)>2 and pts.shape[0]>1:
		ni=pts.shape[1]
		pts=tf.reshape(pts, (-1, pts.shape[-1]))
		bamp=tf.reshape(pts[:, 3], (bsz,-1))
		multmodel=True
		
	else:
		bamp=pts[:, 3][None, :]
		multmodel=False
		
	### when a non c1 symmetry is provided, this will return a list of points
	##  one for each asymmetrical unit so we loop through them and sum the images
	p0=get_sym_pts(sym, pts)
	for p in p0:
		p=tf.transpose(p)
		if multmodel:
			p=tf.reshape(p, (bsz, ni, -1))

		## need to change from (-0.5, 0.5) to actual image coordinates
		bpos=xf2pts(p,ang)
		bpos=bpos*sz+sz/2

		bposf=tf.floor(bpos)
		bposi=tf.cast(bposf,tf.int32)	# integer index
		bposf=bpos-bposf				# remainder used for bilinear interpolation

		# messy tensor math here to implement bilinear interpolation
		bamp0=bamp*(1.0-bposf[:,:,0])*(1.0-bposf[:,:,1])	#0,0
		bamp1=bamp*(bposf[:,:,0])*(1.0-bposf[:,:,1])	#1,0
		bamp2=bamp*(bposf[:,:,0])*(bposf[:,:,1])		#1,1
		bamp3=bamp*(1.0-bposf[:,:,0])*(bposf[:,:,1])	#0,1
		bampall=tf.concat([bamp0,bamp1,bamp2,bamp3],1)
		bposall=tf.concat([bposi,bposi+(1,0),bposi+(1,1),bposi+(0,1)],1)
		imgs=tf.stack([tf.tensor_scatter_nd_add(imgs[i],bposall[i],bampall[i]) for i in range(imgs.shape[0])])

		#try: imgs=tf.tensor_scatter_nd_add(imgs,bposi,bamp)
		#except:
			#print(imgs.shape,bposi.shape,bamp.shape)
			#raise Exception

	fimgs=tf.signal.rfft2d(imgs)

	return (tf.math.real(fimgs)*xfo,tf.math.imag(fimgs)*xfo)

#### compute particle-projection FRC 
##   data_cpx, imgs_cpx - complex images in (real, imag) form
##   rings - indices of Fourier rings
##   return_curve - return the curve instead of average FRC score
##   minpx - skip the X initial low freq pixels
#@tf.function
def calc_frc(data_cpx, imgs_cpx, rings, return_curve=False,minpx=4):
	mreal, mimag=imgs_cpx
	dreal, dimag=data_cpx

	# debugging geometry
	#if random.randint(0,19)==0 :
		#cimg=tf.dtypes.complex(mreal,mimag)
		#img=tf.signal.irfft2d(cimg)
		#num=img.numpy()
		#emd=from_numpy(num)
		#emd.process_inplace("xform.phaseorigin.tocenter")
		#emd.write_image("tests.hdf",-1)

		#cimg=tf.dtypes.complex(dreal,dimag)
		#img=tf.signal.irfft2d(cimg)
		#num=img.numpy()
		#emd=from_numpy(num)
		#emd.process_inplace("xform.phaseorigin.tocenter")
		#emd.write_image("tests.hdf",-1)

	#### normalization per ring
	nrm_img=mreal**2+mimag**2
	nrm_data=dreal**2+dimag**2

	nrm0=tf.tensordot(nrm_img, rings, [[1,2],[0,1]])
	nrm1=tf.tensordot(nrm_data, rings, [[1,2],[0,1]])
	
	nrm=tf.sqrt(nrm0)*tf.sqrt(nrm1)
	nrm=tf.maximum(nrm, 1e-4) #### so we do not divide by 0
		
	#### average FRC per batch
	ccc=mreal*dreal+mimag*dimag
	frc=tf.tensordot(ccc, rings, [[1,2],[0,1]])/nrm
	
	if return_curve:
		return frc
	else:
		frcval=tf.reduce_mean(frc[:, minpx:], axis=1)
		return frcval
	
#### load particles from file and fourier transform them
#### particles need to have their transform in file header or comment of list file
##   will also shrink particles and do FT
##   return Fourier images and transform matrix
##   Two possible types of size reduction, shrinking to achieve boxsz and pre-clipping
##   produce a normal scaled particle in a specified box. This is used for subtomogram
##   data where there is generally padding on the particles.
def load_particles(fname, boxsz, preclip=-1,chunk=None):
	
	nptcl=EMUtil.get_image_count(fname)
	if chunk is not None:
		chunkn=nptcl//chunk[1]
		chunk0=chunkn*chunk[0]
		chunkn=min(nptcl-chunk0,chunkn)
	else: chunk0,chunkn=0,nptcl
	projs=[]
	hdrs=[]
	xfs=[]
	e=EMData(fname, 0, True)
	rawbox=e["nx"] if preclip<=0 else preclip
	print("Loading {} particles of box size {}. shrink to {}".format(nptcl, rawbox, boxsz))
	for j in range(chunk0,chunk0+chunkn,1000):
		print(f"\r {j-chunk0}/{chunkn} ",end="")
		sys.stdout.flush()
		el=EMData.read_images(fname,range(j,min(j+1000,chunk0+chunkn)))
		if el[0]["nx"]!=rawbox:
			for i in range(len(el)): el[i]=el[i].get_clip(Region((el[i]["nx"]-rawbox)//2,(el[i]["ny"]-rawbox)//2,rawbox,rawbox))

		print(f"R     ",end="")
		sys.stdout.flush()
		for e in el:
			xfs.append(e["xform.projection"].get_params("eman"))
			if rawbox!=boxsz:
				#### there is some fourier artifact with xform.scale. maybe worth switching to fft.resample?
				e.process_inplace("math.fft.resample",{"n":rawbox/boxsz})
				#nx=e["nx"]
				#e.clip_inplace(Region((nx-boxsz)//2,(nx-boxsz)//2, boxsz,boxsz))
				#e.process_inplace("xform.scale", {"scale":boxsz/rawbox,"clip":boxsz})
			e.process_inplace("xform.phaseorigin.tocorner")
			hdrs.append(e.get_attr_dict())
			projs.append(e.numpy().copy())
	print(f"  done       ")
	projs=np.array(projs)/1e3
	
	# for tomographic data we need to collect info on "grouping" of 2D particles into 3D particles
	# grpdct is keyed by 3D particle number and has a list of 2D partcile numbers (even after shuffling)
	if "ptcl3d_id" in hdrs[0]:
		grpdct={}
		for i,h in enumerate(hdrs):
			try: grpdct[h["ptcl3d_id"]].append(i)
			except: grpdct[h["ptcl3d_id"]]=[i]
		print(f"identified {len(grpdct)} 3-D particles, merging gradients over each 3-D particle")
	else: grpdct=None
		
	data_cpx=np.fft.rfft2(projs)
	data_cpx=(data_cpx.real.astype(floattype), data_cpx.imag.astype(floattype))

	#xflst=False
	#if fname.endswith(".lst"):		# "old style" LST files?
		#info=load_lst_params(fname)
		#if "xform.projection" in info[0]:
			#xflst=True
			#xfs=[p["xform.projection"].get_params("eman") for p in info]
			#if shuffle:
				#xfs=[xfs[i] for i in rndidx]

	#if xflst==False and ("xform.projection" in hdrs[0]):
		#xflst=True
		#xfs=[p["xform.projection"].get_params("eman") for p in hdrs]
		
	#if xflst==False:
		#xfs=[Transform().get_params("eman") for p in hdrs]
		#print("No existing transform from particles...")
		
	xfsnp=np.array([[x["az"],x["alt"],x["phi"], x["tx"], x["ty"]] for x in xfs], dtype=floattype)
	xfsnp[:,:3]=xfsnp[:,:3]*np.pi/180.
	xfsnp[:,3:]/=float(rawbox)
	
	print("Data read complete")
	return data_cpx,xfsnp,grpdct
	
#### do fourier clipping in numpy then export to tf to save memory...
def get_clip(datacpx, newsz, clipid):
	if (datacpx[0].shape[1])<=newsz:
		dc=datacpx
	else:
		s=newsz//2
		dc=[d[:,clipid[1]<s,:] for d in datacpx]
		dc=[d[:,:,clipid[0]<s+1] for d in dc]
	return [tf.constant(d) for d in dc]

#### compute fourier indices for image generation, clipping, and frc
#### pass indices in a dictionary
def set_indices_boxsz(boxsz, apix=0, return_freq=False):
	idx=np.indices((boxsz,boxsz))-boxsz//2
	idx=np.fft.ifftshift(idx)
	idx=idx[:,:,:boxsz//2+1]
	
	if return_freq:
		#global freq, clipid
		
		freq=np.fft.fftfreq(boxsz, apix)[:boxsz//2]
		ci=idx[0,0]
		cj=idx[1,:,0]
		cj[cj<0]+=1
		clipid=[abs(ci), abs(cj)]
		return clipid
	
	else:
		#global sz, idxft, rrft, rings
		sz=boxsz
		idxft=(idx/sz).astype(floattype)[:, None, :,:]
		rrft=np.sqrt(np.sum(idx**2, axis=0)).astype(floattype)## batch, npts, x-y

		rr=np.round(np.sqrt(np.sum(idx**2, axis=0))).astype(int)
		rings=np.zeros((sz,sz//2+1,sz//2), dtype=floattype) #### Fourier rings
		for i in range(sz//2):
			rings[:,:,i]=(rr==i)
		
		xvec=tf.constant(np.fromfunction(lambda i,j:1.0-2*((i+j)%2),(sz,sz//2+1),dtype=np.float32))
		phaseorg=tf.zeros((boxsz,boxsz//2+1))

		params={"sz":sz, "idxft":idxft, "rrft":rrft, "rings":rings, "xforigin":xvec}
		return params

def build_encoder(ninp,nmid):
	l2=tf.keras.regularizers.l2(1e-3)
	l1=tf.keras.regularizers.l1(1e-3)
	binit=tf.keras.initializers.RandomNormal(0,1e-2)
	kinit=tf.keras.initializers.HeNormal()

	print(f"Encoder {max(ninp//2,nmid)} {max(ninp//8,nmid)} {max(ninp//32,nmid)}")
	layers=[
	tf.keras.layers.Flatten(),
	tf.keras.layers.Dense(max(ninp//2,nmid*8), activation="relu", kernel_initializer=kinit, kernel_regularizer=l2,use_bias=True,bias_initializer=binit),
	tf.keras.layers.Dropout(.3),
	tf.keras.layers.Dense(max(ninp//8,nmid*4), activation="relu", kernel_initializer=kinit, kernel_regularizer=l2,use_bias=True),
	tf.keras.layers.Dense(max(ninp//32,nmid*2), activation="relu", kernel_initializer=kinit, kernel_regularizer=l2,use_bias=True),
	#tf.keras.layers.Dense(max(ninp//2,nmid*2), activation="tanh", kernel_regularizer=l1,use_bias=True,bias_initializer=kinit),
	#tf.keras.layers.Dropout(.3),
	#tf.keras.layers.Dense(max(ninp//8,nmid*2), activation="tanh", kernel_regularizer=l1,use_bias=True),
	#tf.keras.layers.Dense(max(ninp//32,nmid*2), activation="tanh", kernel_regularizer=l1,use_bias=True),
#	tf.keras.layers.BatchNormalization(),
	tf.keras.layers.Dense(nmid, kernel_regularizer=l2, kernel_initializer=kinit,use_bias=True),
	]
		
	encode_model=tf.keras.Sequential(layers)
	return encode_model

def init_segs(segs):
	nmid=len(segs)
	m=np.fromfunction(lambda i,j: ((j//4)==(i//2))*0.1+0.9,(nmid,nmid*2),dtype=float32)
	Localize1

#class Localize1(tf.keras.constraints.Constraint):
    #def __call__(self, w):

        #return w



#### build decoder network. 
# nmid is the number of elements in the middle layer, should match the number of middle layer elements in encoder, must also match number of segments if provided
# pt is either the integer number of gaussians or a numpy array of gaussian coordinates
# segs is integer gaussian numbers beginning each defined region or an empty list (or None)
def build_decoder(nmid, pt,segs ):
	if isinstance(pt, int):
		nout=pt
		initpts=False
		
	else:
		nout=len(pt)
		initpts=True
	
	if segs is None or len(segs)<2 : segs=None
	else: init_segs(segs)

	x0=tf.keras.Input(shape=(nmid))
	
	binit=tf.keras.initializers.RandomNormal(0,1e-2)
	kinit=tf.keras.initializers.HeNormal()
	l2=tf.keras.regularizers.l2(1e-3)
	l1=tf.keras.regularizers.l1(1e-3)
#	layer_output=tf.keras.layers.Dense(nout*4, kernel_initializer=kinit, activation="sigmoid",use_bias=True,kernel_constraint=Localize4())
	layer_output=tf.keras.layers.Dense(nout*4, kernel_initializer=binit, activation="sigmoid",use_bias=True)

	print(f"Decoder {max(nout//32,nmid)} {max(nout//8,nmid)} {max(nout//2,nmid)}")
	layers=[
		#tf.keras.layers.Dense(nmid*2,activation="relu",use_bias=True,bias_initializer=kinit,kernel_constraint=Localize1()),
		#tf.keras.layers.Dense(nmid*4,activation="relu",use_bias=True,kernel_constraint=Localize2()),
		#tf.keras.layers.Dense(nmid*8,activation="relu",use_bias=True,kernel_constraint=Localize3()),
		tf.keras.layers.Dense(nmid*2,activation="relu",kernel_initializer=kinit,use_bias=True,bias_initializer=binit),
		tf.keras.layers.Dense(nmid*4,activation="relu",kernel_initializer=kinit,use_bias=True),
		tf.keras.layers.Dense(nmid*8,activation="relu",kernel_initializer=kinit,use_bias=True),
		tf.keras.layers.Dropout(.3),
#		tf.keras.layers.BatchNormalization(),
		layer_output,
		tf.keras.layers.Reshape((nout,4))
	]

	
	## the five columns are for x,y,z,amp,sigma
	## the range for x,y,z is [-.5, .5]
	## range for amp is [0,1], sigma is [.5, 1.5]
	bshift=np.array([-.5,-.5,-.5,0]).astype(floattype)
	
	y0=x0
	for l in layers:
		y0=l(y0)
		
	y0=y0+tf.constant(bshift)
	gen_model=tf.keras.Model(x0, y0)
	
	## match the bias of the final layer to input
	## need to undo the sigmoid activation
	if initpts:
		bs=pt.copy()-bshift
		bs=np.clip(bs, 1e-6, 1-1e-6)
		bs=-np.log(1./bs-1)
		layer_output.bias.assign(bs.flatten())
	
	return gen_model


#### training decoder on projections
def train_decoder(gen_model, trainset, params, options, pts=None):
	"""pts input can optionally be used as a regularizer if they are known to be good"""
#	opt=tf.keras.optimizers.Adam(learning_rate=options.learnrate)
	opt=tf.keras.optimizers.Adamax()
	wts=gen_model.trainable_variables
	
	nbatch=0
	for t in trainset: nbatch+=1
	
	for itr in range(options.niter):
		cost=[]
		truecost=[]
		for pjr,pji,xf in trainset:
			if xf.shape[0]==1: continue
			pj_cpx=(pjr,pji)
			with tf.GradientTape() as gt:
				# training entropy into the decoder by training individual particles towards random points in latent space
#				if options.decoderentropy: conf=tf.random.normal((xf.shape[0],options.nmid),stddev=0.2)
				if options.decoderentropy: conf=tf.random.uniform((xf.shape[0],options.nmid),minval=-0.02, maxval=0.02)
				# normal behavior, training the neutral map to a latent vector of 0
				else: conf=tf.zeros((xf.shape[0],options.nmid), dtype=floattype)

				if options.decoderentropy: pout=gen_model(conf+tf.random.uniform((conf.shape[0],conf.shape[1]),minval=-0.02, maxval=0.02))
				else: pout=gen_model(conf)
				std=tf.reduce_mean(tf.math.reduce_std(pout, axis=1), axis=0)
				imgs_cpx=pts2img(pout, xf, params, sym=options.sym)
				fval=calc_frc(pj_cpx, imgs_cpx, params["rings"],minpx=options.minressz)
				loss=-tf.reduce_mean(fval)
#				l=loss+std[4]*options.sigmareg+std[3]*5*(options.niter-itr)/options.niter
				l=loss
				if options.modelreg>0: 
					#print(tf.reduce_sum(pout[0,:,:3]*pts[:,:3]),tf.reduce_sum((pout[0,:,:3]-pts[:,:3])**2),len(pts))
					l+=tf.reduce_sum((pout[0,:,:3]-pts[:,:3])**2)/len(pts)*options.modelreg*20.0		# factor of 20 is a rough calibration relative to the dynamic training
				if itr<options.niter//2: l+=std[3]*options.ampreg*options.ampreg
#				print(std)
			
			cost.append(loss)  
			truecost.append(l)
			grad=gt.gradient(l, wts)
			opt.apply_gradients(zip(grad, wts))
			
			sys.stdout.write("\r {}/{}\t{:.3f} ({:.3f})         ".format(len(cost), nbatch, loss,l))
			sys.stdout.flush()
		sys.stdout.write("\r")
		
		print("iter {}, loss : {:.4f} ({:.4f})         ".format(itr, np.mean(cost), np.mean(truecost)))

def eval_model(gen_model, options):
	
	imgs=[]
	xfs=[]
	symmetry=Symmetries.get(options.sym)
	xfs=symmetry.gen_orientations("eman", {"delta":5})
	xfs=[x.get_params("eman") for x in xfs]
	nxf=len(xfs)
	n=7-(nxf-1)%8
	xfs=xfs+[Transform().get_params("eman") for i in range(n)]
	print(nxf, len(xfs))
	
	xfsnp=np.array([[x["az"],x["alt"],x["phi"], x["tx"], x["ty"]] for x in xfs], dtype=floattype)
	xfsnp[:,:3]=xfsnp[:,:3]*np.pi/180.
	xfs=[]
	params=set_indices_boxsz(options.evalsize)
	trainset=tf.data.Dataset.from_tensor_slices((xfsnp))
	trainset=trainset.batch(8)
	for xf in trainset:
		conf=tf.zeros((xf.shape[0],options.nmid), dtype=floattype)
		pout=gen_model(conf)
		imgs_real, imgs_imag=pts2img(pout, xf, params, sym=options.sym)
		imgs_cpx=tf.complex(imgs_real, imgs_imag)
		imgs_out=tf.signal.irfft2d(imgs_cpx)
		imgs.append(imgs_out)
		xfs.append(xf.numpy())
		
	imgs=tf.concat(imgs, axis=0).numpy()
	xfs=np.concatenate(xfs, axis=0)
	imgs=imgs[:nxf]
	xfs=xfs[:nxf]
	
	xfs[:,:3]=xfs[:,:3]*180./np.pi
	xfs[:,3:]*=imgs.shape[-1]
	if os.path.isfile(options.evalmodel):
		os.remove(options.evalmodel)
	for m,x in zip(imgs, xfs):
		e=from_numpy(m)
		x=x.tolist()
		xf=Transform({"type":"eman", "az":x[0], "alt":x[1], "phi":x[2], "tx":x[3], "ty":x[4]})
		e["xform.projection"]=xf
		e["apix_x"]=e["apix_y"]=options.raw_apix
		e.write_image(options.evalmodel,-1)
		
	
def ccf_trans(ref,img):
	ref_real, ref_imag = ref
	img_real, img_imag = img
	ccfr=(img_real*ref_real+img_imag*ref_imag)
	ccfi=(img_real*ref_imag-img_imag*ref_real)
	ccf=tf.complex(ccfr, ccfi)
	ccf=tf.signal.ifftshift(tf.signal.irfft2d(ccf))
	s=ccf.shape[1]//4
	ccf=ccf[:, s:-s, s:-s]
	tx=tf.argmax(tf.math.reduce_max(ccf, axis=1), axis=1)
	ty=tf.argmax(tf.math.reduce_max(ccf, axis=2), axis=1)
	trans=tf.cast(tf.stack([tx, ty], axis=1), floattype)-ccf.shape[1]/2
	return -trans

def translate_image(img, trans):
	imgs_real, imgs_imag=img
	s=idxft[0]*trans[:,0,None,None]+idxft[1]*trans[:,1,None,None]
	s*=-2*np.pi
	imgs_real1=imgs_real*np.cos(s)-imgs_imag*np.sin(s)
	imgs_imag1=imgs_imag*np.cos(s)+imgs_real*np.sin(s)
	return(imgs_real1, imgs_imag1)

def coarse_align(dcpx, pts, options):
	print("coarse aligning particles")
	astep=7.5
	npt=len(dcpx[0])
	symmetry=Symmetries.get(options.sym)
	#xfs=symmetry.gen_orientations("rand", {"n":npt,"phitoo":1,"inc_mirror":1})
	
	
	xfs=symmetry.gen_orientations("eman", {"delta":astep,"phitoo":astep,"inc_mirror":1})
	xfs=[x.get_params("eman") for x in xfs]
	xfsnp=np.array([[x["az"],x["alt"],x["phi"],x["tx"], x["ty"]] for x in xfs], dtype=floattype)
	xfsnp[:,:3]=xfsnp[:,:3]*np.pi/180.
	#return [xfsnp]
	
	bsz=800
	allfrcs=np.zeros((npt, len(xfsnp)))
	alltrans=np.zeros((npt, len(xfsnp), 2))
	niter=len(xfsnp)//bsz
	for ii in range(niter):
		xx=xfsnp[ii*bsz:(ii+1)*bsz]
		projs_cpx=pts2img(pts, xx, sym=options.sym)
		
		for idt in range(npt):
			dc=list(tf.repeat(d[idt:idt+1], len(xx), axis=0) for d in dcpx)
			ts=ccf_trans(dc, projs_cpx)
			dtrans=translate_image(dc,ts)
			frcs=calc_frc(dtrans, projs_cpx,minpx=options.minressz)
			#mxxfs.append(np.argmax(frcs))
			allfrcs[idt,ii*bsz:(ii+1)*bsz]=frcs
			alltrans[idt,ii*bsz:(ii+1)*bsz]=ts

			sys.stdout.write("\r projs: {}/{}, data: {}/{}	".format(ii+1, niter, idt+1, npt))
			sys.stdout.flush()
	xfs=[]
	for ii in range(1):
		fid=np.argmax(allfrcs, axis=1)
		xf=xfsnp[fid]
		ts=alltrans[np.arange(npt),fid]
		xf[:,3:]-=ts/sz
		xfs.append(xf)
		allfrcs[np.arange(npt),fid]=-1
	
	print()
	return xfs

def refine_align(dcpx, xfsnp, pts, options, lr=1e-3):
	nsample=dcpx[0].shape[0]
	trainset=tf.data.Dataset.from_tensor_slices((dcpx[0], dcpx[1], xfsnp))
	trainset=trainset.batch(options.batchsz)
	nbatch=nsample//options.batchsz
	
	opt=tf.keras.optimizers.Adam(learning_rate=lr) 
	xfvs=[]
	frcs=[]
	niter=options.niter
	scr=[]
	
	for ptr,ptj,xf in trainset:
		ptcl_cpx=(ptr, ptj)
		xfvar=tf.Variable(xf)
		p=tf.constant(tf.zeros((xf.shape[0],pts.shape[0], 4))+pts)
		cost=[]
		for itr in range(niter):
			with tf.GradientTape() as gt:
				#xv=tf.concat([xf[:,:3], xfvar[:,3:]], axis=1)
				xv=xfvar
				proj_cpx=pts2img(p, xv, sym=options.sym)
				fval=calc_frc(proj_cpx, ptcl_cpx,minpx=options.minressz)
				loss=-tf.reduce_mean(fval)

			grad=gt.gradient(loss, xfvar)
			opt.apply_gradients([(grad, xfvar)])
			cost.append(loss)

			sys.stdout.write("\r batch {}/{}, iter {}/{}, loss {:.3f}   ".format(len(xfvs), nbatch, itr+1, niter, loss))
			sys.stdout.flush()
			
		xfvs.append(xfvar.numpy())
		frcs.extend(fval.numpy().tolist())
		scr.append([cost[0],cost[-1]])
		
	print()
	xfsnp1=np.vstack(xfvs)
	print(np.mean(abs(xfsnp1-xfsnp), axis=0))
	#xfsnp1[:,:3]=xfsnp[:,:3]
	xfsnp=xfsnp1.copy()
	
	frcs=np.array(frcs)
	
	scr=np.mean(scr, axis=0)
	print("Done. average loss from {:.3f} to {:.3f}".format(scr[0], scr[1]))
	
	return xfsnp, frcs

def save_ptcls_xform(xfsnp, boxsz, options, scr):
	xnp=xfsnp.copy()
	xnp[:,:3]=xnp[:,:3]*180./np.pi
	xnp[:,3:]*=boxsz
	xfs=[Transform({"type":"eman", "az":x[0], "alt":x[1], 
				"phi":x[2], "tx":x[3], "ty":x[4]}) for x in xnp.tolist()]

	oname=options.ptclsout
	print("saving aligned particles to {}".format(oname))
	if os.path.isfile(oname): os.remove(oname)
	lst=LSXFile(oname, False)
	lst0=LSXFile(options.ptclsin, True)
	print(scr)
	for i,xf in enumerate(xfs):
		l0=lst0.read(i)
		d=xf.get_params("eman")
		d["score"]=-scr[i]
		lst.write(-1, l0[0], l0[1], str(d))

	lst=None
	lst0=None
	

### Wrote this as an alterative to using the amplitude-only gradient, in the thought that a finite difference would provide
### better results for this purpose. In the end, it is actually almost indistinguishable, and takes longer to compute
def calc_gqual(trainset, pts, params, options):
	"""Instead of per gaussian variable gradient, try computing the per gaussian quality in neutral conformation"""
	allptqual=[]
	allscr=[]
	
	for pjr,pji,xf in trainset:
#		pt=tf.Variable(pts)
		pt=tf.Variable(tf.repeat(pts, xf.shape[0], axis=0))

		fval=np.zeros((xf.shape[0],pt.shape[1])) # batch size, # points
		pj_cpx=(pjr, pji)
#		print(pt.shape,xf.shape,pts.shape)
		imgs_cpx=pts2img(pt, xf, params, sym=options.sym)
		# FSC for full point projections
		base=calc_frc(pj_cpx, imgs_cpx, params["rings"],minpx=options.minressz)
		allscr.append(base)

		# FSC for omit-1 point projections
		for k in range(pt.shape[1]):
			#ptt=tf.Variable(np.delete(pts,k,0))
			ptt=tf.Variable(tf.repeat(np.delete(pts,k,1), xf.shape[0], axis=0))

#			print(ptt.shape,xf.shape,pts.shape,np.delete(pts,k,1).shape)
			imgs_cpx=pts2img(ptt, xf, params, sym=options.sym)
			fval[:,k]=calc_frc(pj_cpx, imgs_cpx, params["rings"],minpx=options.minressz)-base

		allptqual.append(fval)
		
	allptqual=np.concatenate(allptqual, axis=0)
	allscr=np.concatenate(allscr, axis=0)
	allptqual=allptqual/np.std(allptqual)
	print(" mean score: {:.3f}".format(np.mean(allscr)))
	return allscr, allptqual

### Particles are individually identified by the gradient of the energy with respect to the individual gaussians in the
### neutral configuration. For this reason, the automatic gaussian configuration generation will normally produce
### some excess gaussians in the very low density regions just outside the particle.
def calc_gradient(trainset, pts, params, options):
	allgrds=[]
	allscr=[]
	nbatch=0
	for t in trainset: nbatch+=1

	for pjr,pji,xf in trainset:
		pj_cpx=(pjr, pji)
		with tf.GradientTape() as gt:
			pt=tf.Variable(tf.repeat(pts, xf.shape[0], axis=0))

			imgs_cpx=pts2img(pt, xf, params, sym=options.sym)
			fval=calc_frc(pj_cpx, imgs_cpx, params["rings"],minpx=options.minressz)
			loss=-tf.reduce_mean(fval)

		grad=gt.gradient(loss, pt)*xf.shape[0]
		allgrds.append(grad.numpy().copy())
		allscr.append(fval.numpy().copy())
		sys.stdout.write("\r {}/{} : {:.4f}        ".format(len(allscr), nbatch, np.mean(fval)))
		sys.stdout.flush()


	allgrds=np.concatenate(allgrds, axis=0)
	allscr=np.concatenate(allscr, axis=0)
	allgrds=allgrds/np.std(allgrds)
	print(" mean score: {:.3f}".format(np.mean(allscr)))
	return allscr, allgrds

#### train the conformation manifold from particles
def train_heterg(trainset, pts, encode_model, decode_model, params, options):
	npt=pts.shape[1]
	pas=[int(i) for i in options.pas]
	pas=tf.constant(np.array([pas[0],pas[0],pas[0],pas[1]], dtype=floattype))
	
	## initialize optimizer
#	opt=tf.keras.optimizers.Adam(learning_rate=options.learnrate)
	opt=tf.keras.optimizers.Adamax()
	wts=encode_model.trainable_variables + decode_model.trainable_variables
	nbatch=0
	for t in trainset: nbatch+=1
	
	## Training
	allcost=[]
	for itr in range(options.niter):
		
		i=0
		cost=[]
		for grd,pjr,pji,xf in trainset:
			pj_cpx=(pjr, pji)
			with tf.GradientTape() as gt:
				## from gradient input to the latent space
				conf=encode_model(grd, training=True)
				
							
				## regularization of the latent layer range
				## ideally the output is within a 1-radius circle
				## but we want to make the contraint more soft so it won't affect convergence
				cl=tf.math.sqrt(tf.reduce_sum(conf**2, axis=1))
				cl=tf.reduce_mean(tf.maximum(cl-1,0))
				
				
				## perturb the conformation by a random value
				## similar to the variational autoencoder,
				## but we do not train the sigma of the random value here
				## since we control the radius of latent space already, this seems enough
				conf=(options.perturb*2.0/(itr+2.0))*tf.random.normal(conf.shape)+conf		# perturbation is reduced each iteration
#				conf=.1*tf.random.normal(conf.shape)+conf
				
				## mask out the target columns based on --pas
				pout=decode_model(conf, training=True)
				p0=tf.zeros((xf.shape[0],npt, 4))+pts
				pout=pout*pas+p0*(1-pas)
				
				## finally generate images and calculate frc
				imgs_cpx=pts2img(pout, xf, params, sym=options.sym)
				fval=calc_frc(pj_cpx, imgs_cpx, params["rings"],minpx=options.minressz)
				loss=-tf.reduce_mean(fval)+cl*1e-2
				
				if options.modelreg>0: 
					loss+=tf.reduce_sum((pout[:,:,:3]-pts[:,:,:3])**2)/len(pts)/xf.shape[0]*options.modelreg
			
			cost.append(loss)
			grad=gt.gradient(loss, wts)
			opt.apply_gradients(zip(grad, wts))
			
			i+=1
			if i%10==0: 
				sys.stdout.write("\r {}/{}\t{:.3f}         ".format(len(cost), nbatch, loss))
				sys.stdout.flush()
			
		sys.stdout.write("\r")
		
		print("iter {}, loss : {:.4f}".format(itr, np.mean(cost)))
		allcost.append(np.mean(cost))
		
	return allcost
	
def calc_conf(encode_model, allgrds, bsz=1000):
	
	## conformation output
	mid=[]
	for i in range(len(allgrds)//bsz+1):
		a=allgrds[i*bsz:(i+1)*bsz]
		m=encode_model(a)
		mid.append(m.numpy().copy())
		
	mid=np.concatenate(mid, axis=0)
	return mid
	

	
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	
