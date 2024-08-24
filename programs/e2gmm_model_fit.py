#!/usr/bin/env python
# Muyuan Chen 2023-09
from EMAN2 import *
import numpy as np
from Bio.PDB import *
import protein_constant as e2pc
import scipy.sparse as scipysparse
from sklearn.cluster import KMeans

floattype=np.float32
if "CUDA_VISIBLE_DEVICES" not in os.environ:
	# so we can decide which gpu to use with environmental variable
	os.environ["CUDA_VISIBLE_DEVICES"]='0' 
	
os.environ["TF_FORCE_GPU_ALLOW_GROWTH"]='true' 
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' #### reduce log output
import tensorflow as tf

from e2gmm_refine_new import set_indices_boxsz, load_particles, pts2img, calc_frc, get_clip
from e2gmm_model_refine import calc_clashid, calc_bond, calc_angle, find_clash, compile_chi_matrix, get_rotamer_angle, rotate_sidechain, get_info, calc_dihedral_tf, eval_rama,get_rama_types, calc_rotamer


class Layer_gather(tf.keras.Layer):
	def call(self, x, ii, axis):
		return tf.gather(x, ii, axis=1)

def build_decoder_CA(pts, icls, ninp=4, meanzero=False, freeamp=False):
	nc=int(np.max(icls))+1
	print("building decoder with {} Gaussian, using {} anchor points".format(len(pts[0]), nc))
	
	kinit=tf.keras.initializers.RandomNormal(0,1e-7)
	
	layers=[
		tf.keras.layers.Dense(128, activation="relu"),
		tf.keras.layers.Dense(256, activation="relu"),
		tf.keras.layers.Dense(512, activation="relu"),
		tf.keras.layers.Dropout(.2),
		tf.keras.layers.Dense((nc*5), activation="linear", kernel_initializer=kinit),
		tf.keras.layers.Reshape((nc,5))
		]

	x0=tf.keras.Input(shape=(ninp,))

	y0=x0
	for l in layers:
		y0=l(y0)
	
	y0=Layer_gather()(y0, icls, axis=1)
	y0=tf.keras.layers.Multiply()([y0, tf.constant([1,1,1,0,0])[None,None,:]])
 # 
	# y0=tf.gather(y0, icls, axis=1)
	# y0*=[1,1,1,0,0]
	# if freeamp:
	# 	lz0=tf.keras.layers.Dense((len(pts[0])*5), activation="linear", kernel_initializer=kinit)
	# 	y1=lz0(x0)
	# 	y1=tf.reshape(y1, (-1, len(pts[0]), 5))
	# 	y1=y1*[0,0,0,1,1]
	# 	y0+=y1
		
	if meanzero==False:
		y0=tf.keras.layers.Add()([y0, tf.constant(pts)])
	
	model=tf.keras.Model(x0, y0)
	return model
	
def build_decoder_BB(pts, msk_bb, ninp=4):
	kinit=tf.keras.initializers.RandomNormal(0,1e-7)
	npt=len(pts)
	layers=[
		tf.keras.layers.Dense(128, activation="relu"),
		tf.keras.layers.Dense(256, activation="relu"),
		tf.keras.layers.Dense(512, activation="relu"),
		tf.keras.layers.Dropout(.2),
		tf.keras.layers.Dense((npt*5), activation="linear", kernel_initializer=kinit),
		tf.keras.layers.Reshape((npt,5))
		]

	x0=tf.keras.Input(shape=(ninp,))

	y0=x0
	for l in layers:
		y0=l(y0)

	y0=tf.keras.layers.Multiply()([y0, tf.constant([1,1,1,0,0])[None,None,:]])
	y0=tf.keras.layers.Multiply()([y0, tf.constant(msk_bb[None,:,None])])
	
	model=tf.keras.Model(x0, y0)
	return model
	
def build_decoder_anchor(pts, capos, caidx, pcnt, ninp=4):
	# nc=int(np.max(icls))+1
	nc=len(pcnt)
	print("building decoder with {} Gaussian, using {} anchor points".format(len(pts[0]), nc))
	
	dwt=tf.reduce_sum((capos[None,:,:3]-pcnt[:,None,:3])**2, axis=2)
	dwt=tf.exp(-dwt*20)#*tf.constant(msk, dtype=float)
	dwt=tf.transpose(dwt)
	dwt=tf.cast(dwt, np.float32)
	
	layer_output=tf.keras.layers.Dense(len(capos), activation="linear", use_bias=False)
	
	kinit=tf.keras.initializers.RandomNormal(0,1e-7)
	layers=[
		tf.keras.layers.Dense(128, activation="relu"),
		tf.keras.layers.Dense(256, activation="relu"),
		tf.keras.layers.Dense(512, activation="relu"),
		tf.keras.layers.Dropout(.2),
		tf.keras.layers.Dense((nc*5), activation="linear", kernel_initializer=kinit),
		# tf.keras.layers.Reshape((nc,5))
		tf.keras.layers.Reshape((5, nc)),
		layer_output,
		tf.keras.layers.Permute((2,1))
		]

	x0=tf.keras.Input(shape=(ninp,))

	y0=x0
	for l in layers:
		y0=l(y0)

	layer_output.weights[0].assign(tf.transpose(dwt))
	y0=tf.gather(y0, caidx, axis=1)
	y0*=[1,1,1,0,0]
	y0+=pts
	
	
	model=tf.keras.Model(x0, y0)
	return model
	
def save_model_pdb(gen_model, gen_model_ca, options, fname, thetas=[]):
	
	conf=tf.zeros((1,options.nmid), dtype=floattype)+1.
	pout=gen_model(conf, training=False)
	pout=pout+gen_model_ca(conf, training=False)
	
	if len(thetas)>0:
		pout_rot=pout[:,:,:3]
		for ii in range(4):
			theta=thetas[ii][None,:]
			pout_rot=rotate_sidechain(pout_rot, options.rot_axis[ii], theta, options.rotmat_idx[ii])

		pout=tf.concat([pout_rot, pout[:,:,3:]], axis=-1)

	if options.writetxt:
		np.savetxt(fname+".txt", pout[0].numpy())
		
	atom_pos=(pout[:,:,:3]*[1,-1,-1]+0.5)*options.apix*options.maxboxsz
	atom_pos=atom_pos[0].numpy()

	if options.model.endswith(".cif"):
		pdbpar = MMCIFParser( QUIET = True) 
	else:
		pdbpar = PDBParser( QUIET = True) 
	pdbh = pdbpar.get_structure("model",options.model)

	atoms0=list(pdbh.get_atoms())
	for i,a in enumerate(atoms0):
		a.set_coord(atom_pos[i])
			  
	if options.model.endswith(".cif"):
		io=MMCIFIO()
		fname+=".cif"
	else:
		io=PDBIO()
		fname+=".pdb"
		
	io.set_structure(pdbh)
	io.save(fname)
	print(f"model saved to {fname}")
	
def save_weights(model, fname):
	if os.path.isfile(fname):
		os.remove(fname)
	model.save_weights(fname)
    
def main():
	
	usage="""
	GMM based model flexible fitting.
	
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="folder generated by e2gmm_compile_model", default=None)
	parser.add_argument("--model", type=str,help="alternative pdb model. default is model_input in path.", default=None)
	parser.add_argument("--map", type=str,help="map file for model fitting.", default=None)
	parser.add_argument("--resolution", type=float,help="target resolution.", default=4.)
	parser.add_argument("--learnrate", type=float,help="learning rate.", default=1e-5)
	parser.add_argument("--modelweight", type=float,help="weight for model. default 1", default=1)
	parser.add_argument("--npatch", type=int,help="number of patch for large scale flexible fitting. default is 64", default=64)
	parser.add_argument("--batchsz", type=int,help="batch size. default is 16", default=16)
	parser.add_argument("--niter", type=str,help="number of iterations. default is 20", default="20,20,20")
	parser.add_argument("--rebuild_rotamer", action="store_true", default=False ,help="rebuild all rotamers. slow. require high resolution map.")
	parser.add_argument("--writetxt", action="store_true", default=False ,help="write txt for model in addtion to pdb.")
	parser.add_argument("--load", action="store_true", default=False ,help="load existing data from gmm folder.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv,options.ppid)
	
	options.niter=[int(i) for i in options.niter.split(',')]
	
	path=options.path
	if options.model==None:
		if os.path.isfile(f"{path}/model_input.pdb"):
			options.model=f"{path}/model_input.pdb"
		else:
			options.model=f"{path}/model_input.cif"
	
	##########################
	options.cmd=' '.join(sys.argv)
	fm=f"{options.path}/0_fit_params.json"
	js=js_open_dict(fm)
	js.update(vars(options))
	js.close()	
	
	if options.model.endswith(".cif"):
		pdbpar = MMCIFParser( QUIET = True) 
	else:
		pdbpar = PDBParser( QUIET = True) 
	pdb = pdbpar.get_structure("model",options.model)
	atoms=options.atoms=list(pdb.get_atoms())
	atom_pos=np.array([a.get_coord() for a in atoms])
	residue=list(pdb.get_residues())
	
	print(f"Input model: {options.model}")
	print(f"  {len(residue)} residues, {len(atoms)} atoms.")
	
	options.ptclsin=f"{path}/map_projections.hdf"
	if not os.path.isfile(options.ptclsin):
		cmd=f"e2project3d.py {options.map} --outfile {options.ptclsin} --orientgen eman:delta=4 --parallel thread:12"
		run(cmd)
	else:
		print(f"Using existing projection file {options.ptclsin}...")	
		
	##########################################
	##   load metadata first
	e=EMData(options.ptclsin, 0, True)
	raw_apix, raw_boxsz = e["apix_x"], e["ny"]
	maxboxsz=ceil(raw_boxsz*raw_apix/options.resolution*2)//2*2
	options.maxboxsz=maxboxsz
	options.trainmodel=True
	options.clip=-1
	data_cpx, xfsnp = load_particles(options)
	print("Image size: ", data_cpx[0].shape)

	##   set up fourier indices for image generation/clipping later
	##   params is a dictionay that saves matrices for Fourier indexing
	options.apix=apix=raw_apix*raw_boxsz/maxboxsz
	clipid=set_indices_boxsz(data_cpx[0].shape[1], apix, True)
	params=set_indices_boxsz(maxboxsz)
	
	##########################################
	print("Initializing...")
	p=atom_pos.copy()
	pts=np.zeros((len(p),5), dtype=floattype)
	p=p/options.maxboxsz/options.apix-0.5
	p[:,1:]*=-1; pts[:,:3]=p[:,:3]
	pts[:,3]=.5; pts[:,4]=1
	
	resid=[get_info(a, True) for a in atoms]
	resid, caidx=np.unique(resid, return_inverse=True)
	capos=np.array([np.mean(pts[caidx==i], axis=0) for i in range(np.max(caidx)+1)])
	res_atom_dict={r:np.where(caidx==i)[0] for i,r in enumerate(resid)}
	print("Shape of CA model: ",capos.shape)
	msk_bb=np.array([a.id in ['C','CA','N','O'] for a in atoms], dtype=floattype)
	print(int(np.sum(msk_bb)), "atoms in backbone")
	
	reschain=[r.split('_')[0] for r in resid]
	reschain, reschainid=np.unique(reschain, return_inverse=True)
	ri=np.array([int(r.split('_')[1]) for r in resid])
	ri=ri+reschainid*np.max(ri)
	ri=ri-np.min(ri)
	ri=ri/np.max(ri)
	cp=np.hstack([capos, 10*reschainid[:,None], 10*ri[:,None]])
	# cp=capos#np.hstack([capos, 10*reschainid[:,None], 10*ri[:,None]])
	
	##########################################
	
	if options.niter[0]==0 and options.load:
		icls=np.loadtxt(f"{path}/patch_seg.txt").astype(int)
	else:
		km=KMeans(options.npatch,max_iter=30)
		km.fit(cp)
		pcnt=km.cluster_centers_
		icls=np.zeros(len(pts), dtype=int)
		for i in range(options.npatch):
			ii=np.where(km.labels_==i)[0]
			icls[np.isin(caidx, ii)]=i
		np.savetxt(f"{path}/patch_seg.txt", icls)
	
	# print(np.unique(icls, return_counts=True))
	gen_model=build_decoder_CA(pts[None,...], icls, meanzero=False, freeamp=False)
	# gen_model=build_decoder_anchor(pts[None,...], capos, caidx, pcnt)
	conf=tf.zeros((2,4), dtype=floattype)+1.
	pout=gen_model(conf)
	
	##########################################
	gen_model_ca=build_decoder_CA(pts[None,...], caidx, meanzero=True, freeamp=False)
	conf=tf.zeros((2,4), dtype=floattype)+1.
	d=gen_model_ca(conf)[0]
	
	dcpx=get_clip(data_cpx, params["sz"], clipid)
	trainset=tf.data.Dataset.from_tensor_slices((dcpx[0], dcpx[1], xfsnp))
	trainset=trainset.batch(options.batchsz)

	options.nmid=4
	options.minpx=4
	options.maxpx=maxboxsz//2
	nbatch=int(trainset.cardinality())
	opt=tf.keras.optimizers.Adam(learning_rate=options.learnrate) 
	
	##########################################
	maxpx_morph=options.maxpx//2
	etc=""
	print("Large scale model morphing...")
	wts=gen_model.trainable_variables
	for itr in range(options.niter[0]):
		cost=[]
		costetc=[]
		for pjr,pji,xf in trainset:
			if xf.shape[0]==1: continue
			pj_cpx=(pjr,pji)
			with tf.GradientTape() as gt:

				conf=tf.zeros((xf.shape[0],options.nmid), dtype=floattype)+1.
				pout=gen_model(conf, training=True)
				pout=pout+gen_model_ca(conf, training=True)

				imgs_cpx=pts2img(pout, xf)
				fval=calc_frc(pj_cpx, imgs_cpx, params["rings"], minpx=options.minpx, maxpx=maxpx_morph)
				loss=-tf.reduce_mean(fval)

			cost.append(loss)  

			assert np.isnan(loss.numpy())==False

			grad=gt.gradient(loss, wts)
			opt.apply_gradients(zip(grad, wts))	

			print("{}/{}\t{:.3f}{}         ".format(len(cost), nbatch, loss, etc), end='\r')

		print("iter {}, loss : {:.4f}  {} ".format(itr, np.mean(cost),etc))
		
	if options.niter[0]==0 and options.load:
		gen_model.load_weights(f"{path}/weights_morph.weights.h5")
	else:
		save_model_pdb(gen_model, gen_model_ca, options, f"{path}/fit_00")
		save_weights(gen_model,f"{path}/weights_morph.weights.h5")
	
	conf=tf.zeros((1,options.nmid), dtype=floattype)+1.
	pout=gen_model(conf, training=False)
	pout=pout+gen_model_ca(conf, training=False)
	atom_pos=(pout[:,:,:3]*[1,-1,-1]+0.5)*options.apix*options.maxboxsz
	
	##########################################
	options.bonds=bonds=np.loadtxt(f"{path}/model_bond.txt").astype(floattype)
	options.angle=angle=np.loadtxt(f"{path}/model_angle.txt").astype(floattype)
	options.vdwr_h=np.loadtxt(f"{path}/model_vdwr.txt").astype(floattype)
	options.idx_dih_chi=np.loadtxt(f"{path}/model_dih_chi.txt").astype(int)
	options.idx_dih_rama=np.loadtxt(f"{path}/model_rama_angle.txt").astype(int)
	options.idx_dih_plane=np.loadtxt(f"{path}/model_dih_plane.txt").astype(int)
	options.idx_dih_piptide=np.loadtxt(f"{path}/model_dih_piptide.txt").astype(int)
	if len(options.idx_dih_rama)==0:
		print("No protein")
		options.has_protein=False
	else:
		options.has_protein=True
		options.dih_type=get_rama_types(atoms, options.idx_dih_rama)
		
	ramalevel=[0.0005,0.02]
	options.rama_thr0=1-ramalevel[1]*1.1
	options.rama_thr1=1-ramalevel[0]*2	
	options.thr_plane=np.sin(10*np.pi/180)
	options.thr_piptide=np.sin(30*np.pi/180)
	
	options.clash_nb=128
	options.vdroverlap=.5	
	
	options.rota_thr=0.005
	options.small=small=0.001
	options.rota_thr=np.log(small+options.rota_thr)-np.log(small)
	
	bd=bonds[:,:2].astype(int)
	npt=len(atoms)
	maxlen=3
	bonds_matrix=scipysparse.csr_matrix((np.ones(len(bd)), (bd[:,0], bd[:,1])), shape=(npt,npt))

	d = scipysparse.csgraph.dijkstra(bonds_matrix, directed=False, limit=maxlen)
	options.connect_all=[np.where(i<=maxlen)[0] for i in d]
	options.clashid=calc_clashid(atom_pos[0].numpy(), options, pad=60)
	# print([len(i) for i in options.clashid])
	options.vdw_radius=tf.gather(options.vdwr_h, options.clashid)+options.vdwr_h[:,None]
	options.clash_omask=np.zeros_like(options.clashid)
		
	clash=find_clash(atom_pos, options)
	print(np.sum(clash.numpy()>0)//2, "clashes")	
		
		
	##########################################
	print("C-alpha model refinement...")
	opt=tf.keras.optimizers.Adam(learning_rate=options.learnrate) 
	wts=[]
	wts+=gen_model.trainable_variables
	wts+=gen_model_ca.trainable_variables
	save_weights(gen_model_ca, f"{path}/weights_ca_tmp.weights.h5")
	# weight_model=1e-4
	weight_model=0
	nstd=6.
	for itr in range(options.niter[1]):
		cost=[]
		costetc=[]
		for pjr,pji,xf in trainset:
			if xf.shape[0]==1: continue
			pj_cpx=(pjr,pji)
			with tf.GradientTape() as gt:

				conf=tf.zeros((xf.shape[0],options.nmid), dtype=floattype)+1.
				pout=gen_model(conf, training=True)
				pout=pout+gen_model_ca(conf, training=True)

				imgs_cpx=pts2img(pout, xf)
				fval=calc_frc(pj_cpx, imgs_cpx, params["rings"], minpx=options.minpx, maxpx=options.maxpx)
				loss=-tf.reduce_mean(fval)
				
				
				atom_pos=(pout[:,:,:3]*[1,-1,-1]+0.5)*options.apix*options.maxboxsz

				bond_len=calc_bond(atom_pos, bonds[:,:2].astype(int))
				bond_df=(bond_len-bonds[:,2])/bonds[:,3]
				bond_outlier=tf.maximum(abs(bond_df)-nstd, 0)
				bond_outlier=tf.reduce_mean(bond_outlier)*1000*20


				ang_val=calc_angle(atom_pos, angle[:,:3].astype(int))
				ang_df=(ang_val-angle[:,3])/angle[:,4]
				ang_outlier=tf.maximum(abs(ang_df)-nstd, 0)
				ang_outlier=tf.reduce_mean(ang_outlier)*1000*20

				clash=find_clash(atom_pos, options)
				
				nclash=np.sum(clash.numpy()>0)//2
				clash_score=tf.reduce_sum(clash)/conf.shape[0]/2.

				lossetc=0.
				lossetc+=bond_outlier
				lossetc+=ang_outlier
				lossetc+=clash_score
				

				l=loss*1+lossetc*weight_model

			cost.append(loss)
			costetc.append(lossetc)

			assert np.isnan(l.numpy())==False

			grad=gt.gradient(l, wts)
			opt.apply_gradients(zip(grad, wts))
			etc=""
			etc+=f", {bond_outlier:.3f}, {ang_outlier:.3f}, {clash_score:.3f}, {nclash}"


			print("{}/{}\t{:.3f}{}         ".format(len(cost), nbatch, loss, etc), end='\r')

		if itr==0:
			c0=np.mean(cost)
			c1=np.mean(costetc)
			
		elif itr==1:
			d0=c0-np.mean(cost)
			d1=np.mean(costetc)-c1
			weight_model=abs(d0/d1)*options.modelweight
			print(f"iter 0: FRC loss decrease {d0:.3e}, geometry loss increase {d1:.3e}")
			print(f"        using geometry weight {weight_model:.3e}")
			gen_model.load_weights(f"{path}/weights_morph.weights.h5")
			gen_model_ca.load_weights(f"{path}/weights_ca_tmp.weights.h5")
			print("iter,   loss,   bond outlier,  angle outlier, clash_score,  number of clash")
		else:
			print("iter {}, loss : {:.4f},  {} ".format(itr, np.mean(cost),etc))
		
	if options.niter[1]==0 and options.load:
		gen_model_ca.load_weights(f"{path}/weights_ca.weights.h5")
	else:
		save_model_pdb(gen_model, gen_model_ca, options, f"{path}/fit_01")
		save_weights(gen_model_ca, f"{path}/weights_ca.weights.h5")
		save_weights(gen_model, f"{path}/weights_morph.weights.h5")

	##########################################
	print("Compiling rotamers...")
	if len(options.idx_dih_chi)==0:
		print("No rotamers.")
		theta_all_out=[np.zeros(0, dtype=floattype) for ii in range(4)]
		
	else:
		options.chi_idx, options.chi_mat=compile_chi_matrix(options.idx_dih_chi)

		rot_axis=[]
		rotmat_idx=[]
		for i in range(4):
			ra, ri=get_rotamer_angle(atoms, i)    
			rot_axis.append(ra)
			rotmat_idx.append(ri)

		options.rot_axis=rot_axis
		options.rotmat_idx=rotmat_idx
		for i,m in enumerate(options.chi_mat):
			print(f"  rotamer chi = {i}, matrix shape: {m.shape}, total {len(rot_axis[i])} angles")

		### additional rotation for each chi angle
		theta_all_out=[np.zeros(len(rot_axis[ii]), dtype=floattype) for ii in range(4)]
		
		if options.rebuild_rotamer:
			print("Re-building rotamers...")
			chinid=[3,2,1,0]
			trainset1=tf.data.Dataset.from_tensor_slices((dcpx[0], dcpx[1], xfsnp))
			trainset1=trainset1.batch(32)
		else:
			chinid=[]
			
		##########################################
		for chin in chinid:
			print(f"chin = {chin}")
			conf=tf.zeros((1,options.nmid), dtype=floattype)+1.
			pout=gen_model(conf, training=False)
			pout=pout+gen_model_ca(conf, training=False)

			pout_rot=pout[:,:,:3]
			for ii in range(4):
				theta=theta_all_out[ii][None,:]
				pout_rot=rotate_sidechain(pout_rot, options.rot_axis[ii], theta, options.rotmat_idx[ii])

			pout=tf.concat([pout_rot, pout[:,:,3:]], axis=-1)


			atom_pos=(pout[:,:,:3]*[1,-1,-1]+0.5)*options.apix*options.maxboxsz

			clash=find_clash(atom_pos, options)
			nclash00=tf.reduce_sum(tf.sign(clash))

			for pjr,pji,xf in trainset1:
				pj_cpx=(pjr,pji)
				imgs_cpx=pts2img(pout, xf)
				fval=calc_frc(pj_cpx, imgs_cpx, params["rings"], minpx=options.minpx, maxpx=options.maxpx)
				loss0=-tf.reduce_mean(fval)
				break

			print(f"clash {nclash00:.0f}, loss {loss0:.4f}")

			ii_mat=options.chi_mat[chin].copy()
			ii=options.chi_idx[chin].copy()
			iipro=ii_mat[:,-1,-1]==0
			ii_mat=ii_mat[~iipro]
			ii=ii[~iipro]
			
			ii=ii.T.flatten()
			ii=options.idx_dih_chi[ii][:,:4]
			rota_id=ii[:,1].reshape(chin+1, -1)[0]
			at=[get_info(atoms[i],True) for i in rota_id]
			rid_atoms=[res_atom_dict[a] for a in at]

			pt=tf.gather(atom_pos, ii, axis=1)
			dih=calc_dihedral_tf(pt)%360
			dih=tf.reshape(dih, (atom_pos.shape[0], chin+1, -1))
			dih=tf.transpose(dih, (0,2,1))[0]
			
			theta_idx=[]
			keys=[get_info(atoms[i],True) for i in rota_id]
			for ii in range(chin+1): 
				ir=options.rot_axis[ii][:,0]
				dic={get_info(atoms[r],True):i for i,r in enumerate(ir)}
				ir=np.array([dic[k] for k in keys])
				theta_idx.append(ir)

			theta_idx=np.array(theta_idx).T

			conf=tf.zeros((1,options.nmid), dtype=floattype)+1.
			pout=gen_model(conf, training=False)
			pout=pout+gen_model_ca(conf, training=False)

			cost_rota=[]
			for rotai in range(ii_mat.shape[1]):

				rota_target=ii_mat[:,rotai,:chin+1]
				rota_cur=dih[:,:chin+1]
				ddr=rota_cur-rota_target
				ddr=ddr*np.pi/180.

				cost=[]
				for ir in range(len(ddr)):
					theta_all=[t.copy() for t in theta_all_out]

					for ii in range(chin+1):
						theta_all[ii][theta_idx[ir,ii]]=ddr[ir,ii]

					pout_rot=pout[:,:,:3]
					for ii in range(4):
						theta=theta_all[ii][None,:]
						pout_rot=rotate_sidechain(pout_rot, options.rot_axis[ii], theta, options.rotmat_idx[ii])

					atom_pos_rot=(pout_rot*[1,-1,-1]+0.5)*options.apix*options.maxboxsz
					clash=find_clash(atom_pos_rot, options)
					nclash=tf.reduce_sum(tf.sign(clash))

					pout1=tf.concat([pout_rot, pout[:,:,3:]], axis=-1)
					imgs_cpx=pts2img(pout1, xf)
					fval=calc_frc(pj_cpx, imgs_cpx, params["rings"], minpx=options.minpx, maxpx=options.maxpx)
					loss1=-tf.reduce_mean(fval)
					dloss=(loss1-loss0)*1e5
					print(f"rota {rotai}/{ii_mat.shape[1]}, {ir}/{len(ddr)}: clash {nclash:.0f}, loss diff {dloss:.4f}", end='\r')

					cost.append([nclash, dloss])

				cost=np.array(cost)
				c=cost[:,1]
				d=cost[:,0]-float(nclash00)
				d+=np.maximum(0, d-3)*10
				c+=d*.1
				
				cost_rota.append(c.copy())
			print()
			
			cr=np.array(cost_rota)
			ci=np.argmin(cr, axis=0)
			rota_target=ii_mat[np.arange(len(ii_mat)),ci,:chin+1]
			rota_cur=dih[:,:chin+1]
			ddr=rota_cur-rota_target
			ddr=ddr*np.pi/180.
			ddr=ddr.numpy()
			
			cc=np.min(cr, axis=0)
			# print("accepting {}/{} rotamers...".format(np.sum(cc>0), len(cc)))
			# ddr[cc>0,:]*=0

			theta_all=[t.copy() for t in theta_all_out]
			for ii in range(chin+1):
				theta_all[ii][theta_idx[:,ii]]=ddr[:,ii]

			pout_rot=pout[:,:,:3]
			for ii in range(4):
				theta=theta_all[ii][None,:]
				pout_rot=rotate_sidechain(pout_rot, options.rot_axis[ii], theta, options.rotmat_idx[ii])

			pout1=tf.concat([pout_rot, pout[:,:,3:]], axis=-1)
			atom_pos_rot=(pout_rot*[1,-1,-1]+0.5)*options.apix*options.maxboxsz
			clash=find_clash(atom_pos_rot, options)

			cl=clash[0].numpy()
			rid_clash=np.array([np.sum(cl[i]>0) for i in rid_atoms])    
			
			for ii in range(chin+1):
				it=theta_all[ii]!=0
				theta_all_out[ii][it]=theta_all[ii][it]

			conf=tf.zeros((1,options.nmid), dtype=floattype)+1.
			pout=gen_model(conf, training=False)
			pout=pout+gen_model_ca(conf, training=False)

			pout_rot=pout[:,:,:3]
			for ii in range(4):
				theta=theta_all_out[ii][None,:]
				pout_rot=rotate_sidechain(pout_rot, options.rot_axis[ii], theta, options.rotmat_idx[ii])

			pout1=tf.concat([pout_rot, pout[:,:,3:]], axis=-1)
		
			
			cost=[]
			for pjr,pji,xf in trainset:
				pj_cpx=(pjr,pji)
				imgs_cpx=pts2img(pout, xf)
				fval=calc_frc(pj_cpx, imgs_cpx, params["rings"], minpx=options.minpx, maxpx=options.maxpx)
				loss0=-tf.reduce_mean(fval)

				imgs_cpx=pts2img(pout1, xf)
				fval=calc_frc(pj_cpx, imgs_cpx, params["rings"], minpx=options.minpx, maxpx=options.maxpx)
				loss1=-tf.reduce_mean(fval)

				cost.append([loss0, loss1])

			cost=np.array(cost)
			d=cost[:,0]-cost[:,1]
			print(" loss from {:.3f} to {:.3f}".format(np.mean(cost[:,0]), np.mean(cost[:,1])))
		
			save_model_pdb(gen_model, gen_model_ca, options, f"{path}/fit_02", theta_all_out)
			
			np.savetxt(f"{path}/chi_theta.txt", np.concatenate(theta_all_out))
		
	if options.rebuild_rotamer==False and options.load==True:
		print("load rotamer from file")
		cc=np.loadtxt(f"{path}/chi_theta.txt").astype(floattype)
		k=0
		for i, t in enumerate(theta_all_out):
			theta_all_out[i]=cc[k:k+len(t)]
			k+=len(t)
		
	##########################################
	print("Full atom refinement...")	
	tf_theta=[]
	gen_model_bb=build_decoder_BB(pts[None,...],  msk_bb)
	if options.load:
		gen_model_bb.load_weights(f"{path}/weights_bb.weights.h5")
	
	save_weights(gen_model_bb, f"{path}/weights_bb.weights.h5")
	conf=tf.zeros((2,4), dtype=floattype)+1.
	d=gen_model_bb(conf)[0]
		
	opt=tf.keras.optimizers.Adam(learning_rate=options.learnrate) 
	wts=[]
	if options.has_protein:
		tf_theta=[tf.Variable(t) for t in theta_all_out]
		wts+=tf_theta
	wts+=gen_model.trainable_variables
	wts+=gen_model_ca.trainable_variables
	wts+=gen_model_bb.trainable_variables
	
	weight_model=0
	rota_score=rota_outlier=plane_score=rama_score=rama_outlier=0
	
	for itr in range(options.niter[2]):
			
		cost=[]
		costetc=[]
		for pjr,pji,xf in trainset:
			if xf.shape[0]==1: continue
			pj_cpx=(pjr,pji)
			with tf.GradientTape() as gt:

				conf=tf.zeros((xf.shape[0],options.nmid), dtype=floattype)+1.
				pout=gen_model(conf, training=True)
				pout=pout+gen_model_ca(conf, training=True)
				pout=pout+gen_model_bb(conf, training=True)
				
				if options.has_protein:
					pout_rot=pout[:,:,:3]
					for ii in range(4):
						theta=tf_theta[ii][None,:]
						theta=tf.repeat(theta, conf.shape[0], axis=0)
						pout_rot=rotate_sidechain(pout_rot, options.rot_axis[ii], theta, options.rotmat_idx[ii])

					pout=tf.concat([pout_rot, pout[:,:,3:]], axis=-1)
				
				imgs_cpx=pts2img(pout, xf)
				fval=calc_frc(pj_cpx, imgs_cpx, params["rings"], minpx=options.minpx, maxpx=options.maxpx)
				loss=-tf.reduce_mean(fval)
				
				####################
				atom_pos=(pout[:,:,:3]*[1,-1,-1]+0.5)*options.apix*options.maxboxsz
				lossetc=0.
				
				bond_len=calc_bond(atom_pos, bonds[:,:2].astype(int))
				bond_df=(bond_len-bonds[:,2])/bonds[:,3]
				bond_outlier=tf.maximum(abs(bond_df)-4.5, 0)
				bond_outlier=tf.reduce_mean(bond_outlier)*1000*20
				lossetc+=bond_outlier

				##############
				ang_val=calc_angle(atom_pos, angle[:,:3].astype(int))
				ang_df=(ang_val-angle[:,3])/angle[:,4]
				ang_outlier=tf.maximum(abs(ang_df)-4.5, 0)
				ang_outlier=tf.reduce_mean(ang_outlier)*1000*20
				lossetc+=ang_outlier
				
				##############
				clash=find_clash(atom_pos, options)
				clash_score=tf.reduce_sum(clash)/conf.shape[0]/2.
				lossetc+=clash_score
				
				if options.has_protein:
					##############
					pt=tf.gather(atom_pos, options.idx_dih_rama, axis=1)
					phi=calc_dihedral_tf(pt[:,:,:4])
					psi=calc_dihedral_tf(pt[:,:,4:])
					rama=eval_rama(phi, psi, options)
					rama_score=tf.reduce_mean(rama)*.1
					rama_outlier=tf.reduce_mean(tf.maximum(0,rama-options.rama_thr0))*500*5
					rama_outlier+=tf.reduce_mean(tf.maximum(0,rama-options.rama_thr1))*1000*1000*5

					lossetc+=rama_score 
					lossetc+=rama_outlier 
				
				##############
				pt=tf.gather(atom_pos, options.idx_dih_plane, axis=1)
				rot=calc_dihedral_tf(pt)
				rot=tf.sin(rot*np.pi/180.)
				rot=tf.maximum(0, abs(rot)-options.thr_plane)
				plane_score=tf.reduce_mean(rot)*1000
				
				if options.has_protein:
					##############
					pt=tf.gather(atom_pos, options.idx_dih_piptide, axis=1)
					rot=calc_dihedral_tf(pt)
					rot=tf.sin(rot*np.pi/180.)
					rot=tf.maximum(0, abs(rot)-options.thr_piptide)
					plane_score+=tf.reduce_mean(rot)*1000
				
				lossetc+=plane_score
				
				
				if options.has_protein:
					##############
					rota_out=[]
					for chin in range(4):
						d=calc_rotamer(atom_pos, chin, options)
						rota_out.append(d)

					big=np.log(options.small+1)-np.log(options.small)
					rota_out=tf.concat(rota_out, axis=1)            
					rota_score=tf.reduce_mean(big-rota_out)

					r1=tf.maximum(0,options.rota_thr-rota_out)
					rota_outlier=tf.reduce_mean(r1)*10000
					
					lossetc+=rota_score
					lossetc+=rota_outlier

				l=loss*1+lossetc*weight_model
				
			cost.append(loss)  
			costetc.append(lossetc)  

			assert np.isnan(l.numpy())==False

			grad=gt.gradient(l, wts)
			opt.apply_gradients(zip(grad, wts))
			etc=""
			etc+=f", {bond_outlier:.1f}, {ang_outlier:.1f}, {clash_score:.1f}"
			etc+=f", {rota_score:.1f}, {rota_outlier:.2f}, {plane_score:.1f}"
			etc+=f", {rama_score:.2f}, {rama_outlier:.1f}"


			print("{}/{}\t{:.3f}{}         ".format(len(cost), nbatch, loss, etc), end='\r')

		if itr==0:
			c0=np.mean(cost)
			c1=np.mean(costetc)
			# print(c0, c1)
			
		elif itr==1:
			# d0=cost[0]-cost[-1]
			# d1=costetc[0]-costetc[-1]
			d0=c0-np.mean(cost)
			d1=np.mean(costetc)-c1
			weight_model=abs(d0/d1)*options.modelweight
			print(f"\niter 1: FRC loss {c0:.4f}, decrease {d0:.3e}, geometry loss {c1:.2f}, increase {d1:.3e}")
			print(f"        using geometry weight {weight_model:.3e}")
			
			gen_model.load_weights(f"{path}/weights_morph.weights.h5")
			gen_model_ca.load_weights(f"{path}/weights_ca.weights.h5")
			gen_model_bb.load_weights(f"{path}/weights_bb.weights.h5")
			opt=tf.keras.optimizers.Adam(learning_rate=options.learnrate) 
			wts=[]
			if options.has_protein: 
				tf_theta=[tf.Variable(t) for t in theta_all_out]
				wts+=tf_theta
				
			wts+=gen_model.trainable_variables
			wts+=gen_model_ca.trainable_variables
			wts+=gen_model_bb.trainable_variables
			print("iter,   loss,   bond, angle, clash, rotamer score, outlier, plane score, rama score, outlier")
		else:
			print("iter {}, loss : {:.4f} {} ".format(itr, np.mean(cost),etc))
		
	save_model_pdb(gen_model, gen_model_ca, options, f"{path}/fit_03", tf_theta)
	if options.niter[2]>0:
		save_weights(gen_model, f"{path}/weights_morph.weights.h5")
		save_weights(gen_model_ca, f"{path}/weights_ca.weights.h5")
		save_weights(gen_model_bb,f"{path}/weights_bb.weights.h5")
		if options.has_protein: 
			np.savetxt(f"{path}/chi_theta.txt", np.concatenate([t.numpy() for t in tf_theta]))
		
	E2end(logid)
	
if __name__ == '__main__':
	main()
	
