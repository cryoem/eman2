#!/usr/bin/env python
# Muyuan Chen 2023-01

import numpy as np
#import matplotlib.pyplot as plt
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
##   in theory float16 also works but it can be unsafe especially when the network is deeper...
floattype=np.float32
from e2gmm_refine_new import *

def xf2pts_one(args):
	return xf2pts(args[0], args[1])

@tf.function
def xf2pts_mult(p0, xf):
	return tf.vectorized_map(xf2pts_one,(p0, xf))

def main():
	
	usage="Convert particle conformations from e2gmm_refine_new to orientations at the masked region."
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="path", default=None)
	parser.add_argument("--mask", type=str,help="mask", default=None)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	
	logid=E2init(sys.argv)
	path=options.path
	itr=0
	itnew=1
	
	for eo in ["even", "odd"]:
		print(f"Realigning {eo}...")
		pts=np.loadtxt(f"{path}/model_{itr:02d}_{eo}.txt")

		imsk=make_mask_gmm(options.mask, pts).numpy()
		#plt.scatter(p[:,1], p[:,2], c=imsk,  alpha=.2, cmap="RdBu_r")
		
		decode_model=tf.keras.models.load_model(f"{path}/dec_{eo}.h5",compile=False,custom_objects={"ResidueConv2D":ResidueConv2D})
		midraw=np.loadtxt(f"{path}/midall_{itr:02d}_{eo}.txt")[:,1:]

		pts=np.loadtxt(f"{path}/model_{itr:02d}_{eo}.txt")
		
		pfile=f"{path}/ptcls_{itr:02d}_{eo}.lst"
		lst=load_lst_params(pfile)
		xfs=[p["xform.projection"].get_params("eman") for p in lst]
		
		e=EMData(pfile, 0, True)
		rawbox=e["nx"]
		xfsnp=np.array([[x["az"],x["alt"],x["phi"], x["tx"], x["ty"]] for x in xfs], dtype=floattype)
		xfsnp[:,:3]=xfsnp[:,:3]*np.pi/180.
		xfsnp[:,3:]/=float(rawbox)

		bsz=128
		trainset=tf.data.Dataset.from_tensor_slices((midraw, xfsnp))
		trainset=trainset.batch(bsz)
		opt=tf.keras.optimizers.Adam(learning_rate=1e-3) 
		nbatch=len(xfsnp)//bsz
		niter=30
		xfvs=[]
		scores=[]
		for conf,xf in trainset:
			xfvar=tf.Variable(xf)
			p0=tf.constant(tf.zeros((xf.shape[0],pts.shape[0], 3))+pts[:,:3])

			pout=decode_model(conf)[:,:,:3]
			pj0=xf2pts_mult(pout, xf)


			cost=[]
			for it in range(niter):
				with tf.GradientTape() as gt:
					pj1=xf2pts_mult(p0, xfvar)
					dpj=pj1-pj0
					dpj=tf.reduce_sum(dpj**2, axis=2)*imsk
					# print(dpj.shape)
					loss1=tf.reduce_sum(dpj, axis=1)
					loss=tf.reduce_mean(loss1)*100
					if it==0: loss0=loss


				grad=gt.gradient(loss, xfvar)
				opt.apply_gradients([(grad, xfvar)])
				cost.append(loss)

			sys.stdout.write("\r batch {}/{}, loss {:.3f} -> {:.3f}   ".format(len(xfvs), nbatch, loss0, loss))
			sys.stdout.flush()

			xfvs.append(xfvar.numpy())
			scores.append(loss1.numpy())

		xfsnp1=np.vstack(xfvs)
		score=np.concatenate(scores)
		score-=np.max(score)
		score/=-np.min(score)
		
		
		xnp=xfsnp1.copy()
		xnp[:,:3]=xnp[:,:3]*180./np.pi
		xnp[:,3:]*=rawbox
		xfs=[Transform({"type":"eman", "az":x[0], "alt":x[1], 
					"phi":x[2], "tx":x[3], "ty":x[4]}) for x in xnp.tolist()]

		oname=f"{path}/ptcls_{itnew:02d}_{eo}.lst"
		print("saving aligned particles to {}".format(oname))
		lstin=load_lst_params(f"{path}/ptcls_{itr:02d}_{eo}.lst")
		#frcs=[s["score"] for s in lstin]
		for i,xf in enumerate(xfs):
			lstin[i]["xform.projection"]=xf
			lstin[i]["score"]=score[i]#frcs[i]

		if os.path.isfile(oname): os.remove(oname)
		save_lst_params(lstin, oname)


		run(f"e2spa_make3d.py --input {path}/ptcls_{itnew:02d}_{eo}.lst --output {path}/threed_{itnew:02d}_{eo}.hdf --keep .9 --parallel thread:32 --sym c1")
		run(f"e2proc3d.py {path}/threed_{itnew:02d}_{eo}.hdf {path}/threed_raw_{eo}.hdf")
		
		
	E2end(logid)
	
	
if __name__ == '__main__':
	main()
	
