#!/usr/bin/env python
# Muyuan Chen 2024-02

import numpy as np
from EMAN2 import *
import protein_constant as e2pc
from Bio.PDB import *

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
from e2gmm_model_refine import get_info

def main():
	
	usage="."
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--pdb", type=str,help="input model. pdb or cif", default=None)
	parser.add_argument("--mrc", type=str,help="input map for later comparison. any map format", default=None)

	parser.add_argument("--output", type=str,help="simulated map file output", default=None)
	parser.add_argument("--batchsz", type=int,help="", default=16)
	parser.add_argument("--res", type=float,help="max resolution for comparison", default=-1)
	parser.add_argument("--skiph",action="store_true",default=False,help="skip H atoms")

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	
	logid=E2init(sys.argv)

	pdbpar = PDBParser( QUIET = True)
	pdb = pdbpar.get_structure("model",options.pdb)

	residue=list(pdb.get_residues())

	if options.skiph:
		for r in residue:
			d=list(r.child_dict)
			for a in d:
				if a[0]=='H':# and len(a)>1:
					r.detach_child(a)

	atoms=list(pdb.get_atoms())
	atom_pos=np.array([a.get_coord() for a in atoms])
	print(f"{len(atoms)} atoms from {len(residue)} residues")

	map_proj=f"{options.output[:-4]}_projection.hdf"
	if os.path.isfile(map_proj):
		print(f"using existing map projection from {map_proj}")
	else:
		run(f"e2project3d.py {options.mrc} --outfile {map_proj} --orientgen=eman:delta=4 --parallel=thread:16")


	options.ptclsin=map_proj
	e=EMData(options.ptclsin, 0, True)
	raw_apix, raw_boxsz = e["apix_x"], e["ny"]
	print(f"pixel size {raw_apix}, box size {raw_boxsz}")

	if options.res<0:
		maxboxsz=options.maxboxsz=e["nx"]
		apix=options.apix=e["apix_x"]

	else:
		maxboxsz=options.maxboxsz=ceil(raw_boxsz*raw_apix/options.res*2)//2*2
		options.apix=apix=raw_apix*raw_boxsz/maxboxsz
		print(f"downsample to pixel size {apix}, box size {maxboxsz}")

	options.sym="c1"
	options.trainmodel=True
	options.clip=-1
	data_cpx, xfsnp = load_particles(options)

	clipid=set_indices_boxsz(data_cpx[0].shape[1], apix, True)
	params=set_indices_boxsz(maxboxsz)

	atom_type=[a.element for a in atoms]
	print("model include atom types",np.unique(atom_type))


	p0=atom_pos.copy()
	p0=p0/options.maxboxsz/options.apix-0.5
	p0[:,1:]*=-1

	a_fac=0.11130926
	s_fac=0.7090085
	pp=[]
	for i,k in enumerate(atom_type):
		sc=e2pc.form_fac[k]
		sc=np.array(sc).reshape((2,-1)).copy()/10
		p=p0[i]
		for c in sc.T:
			s=1./c[1]/s_fac
			a=np.sqrt(c[0]/c[1]/s_fac/a_fac)
			pp.append([p[0],p[1],p[2],a,s])
	pts=np.array(pp, dtype=floattype)
	print("model shape", pts.shape)


	trainset=tf.data.Dataset.from_tensor_slices((data_cpx[0], data_cpx[1], xfsnp))
	trainset=trainset.batch(16)
	nbatch=int(trainset.cardinality())

	#### simple network for amplitude optimization
	shp=pts.shape
	layers=[
		tf.keras.layers.Dense(64,activation="relu"),
		tf.keras.layers.Dense(128,activation="relu"),
		tf.keras.layers.Dense(256,activation="relu"),
		tf.keras.layers.Dropout(.1),
		tf.keras.layers.Dense(np.prod(shp)//5,activation="linear"),
		tf.keras.layers.Reshape((shp[0]//5, shp[1]))
	]
	model=tf.keras.Sequential(layers)
	o=model(tf.zeros((1,4), dtype=floattype))

	wts=model.trainable_variables
	opt=tf.keras.optimizers.Adam(learning_rate=1e-4)
	maxpx=ceil(maxboxsz/(1.0/(apix)))

	########################
	for itr in range(10):
		cost=[]
		for pjr,pji,xf in trainset:
			# if xf.shape[0]==1: continue
			pj_cpx=(pjr,pji)
			with tf.GradientTape() as gt:

				conf=tf.zeros((xf.shape[0],4), dtype=floattype)+1.
				pout=model(conf, training=True)*0.01
				pout=pout*[0,0,0,1,0]+1

				pout=tf.repeat(pout, 5, axis=1)
				pout=pout*pts

				imgs_cpx=pts2img(pout, xf)
				fval=calc_frc(pj_cpx, imgs_cpx, params["rings"], minpx=2, maxpx=maxpx)
				loss=-tf.reduce_mean(fval)

			cost.append(loss)

			assert np.isnan(loss.numpy())==False

			grad=gt.gradient(loss, wts)
			opt.apply_gradients(zip(grad, wts))

			print("{}/{}\t{:.3f}         ".format(len(cost), nbatch, loss), end='\r')

		print("iter {}, loss : {:.4f}  ".format(itr, np.mean(cost)))

	##################

	conf=tf.zeros((1,4), dtype=floattype)+1.
	pout=model(conf, training=False)*0.01
	pout=pout*[0,0,0,1,0]+1

	pout=tf.repeat(pout, 5, axis=1)
	pout=pout*pts

	p1=pout.numpy()[0]

	options.gmm_out=f"{options.output[:-4]}_model.txt"

	np.savetxt(options.gmm_out, p1)

	options.raw_apix=raw_apix
	options.evalsize=raw_boxsz
	options.evalmodel=f"{options.output[:-4]}_model_proj.hdf"
	eval_model(p1[None,...], options, usepts=True)


	run(f"e2spa_make3d.py --input {options.evalmodel} --output {options.output} --thread 32")
	run(f"e2proc3d.py {options.output} {options.output} --matchto {options.mrc}")
	run(f"e2proc3d.py {options.output} {options.output[:-4]}_fsc.txt --calcfsc {options.mrc}")
	run(f"e2proc3d.py {options.output} {options.output[:-4]}_sub.hdf --mult -1 --addfile {options.mrc} --process filter.lowpass.gauss:cutoff_abs=.4 --process normalize.edgemean")

	E2end(logid)
	
	
if __name__ == '__main__':
	main()
	
