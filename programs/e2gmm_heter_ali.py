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

def xf2pts(pts, ang):
	""" Project 3d Gaussian coordinates based on transforms to make projection
		input:  pts - ( number of Gaussian, 3 (x,y,z) )
				ang - ( 5 (az, alt, phi, tx, ty) )
	"""

	#### input EMAN style euler angle (az, alt, phi) and make projection matrix
	##   note we need to be able to deal with a batch of particles at once
	##   so everything is in matrix form	
	azp=-ang[0]
	altp=ang[1]
	phip=-ang[2]

	#### rotate Gaussian positions
	matrix=make_matrix(azp, altp, phip)
	matrix=tf.reshape(matrix, shape=[3,3]) 
	matrix=tf.transpose(matrix)

	pts_rot=tf.matmul(pts, matrix)
#     pts_rot=tf.transpose(pts_rot)

	#### finally do the translation
	pts_rot_trans=tf.stack([(pts_rot[:,0]+ang[3]), (-pts_rot[:,1])+ang[4]], 1)
	
	#pts_rot_trans=pts_rot_trans*sz+sz/2
	return pts_rot_trans


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
	parser.add_argument("--iterid", type=int, help="iteration id for midall_xx.txt",default=0)
	parser.add_argument("--newiter", type=int, help="iteration number to write. default is iterid+1",default=-1)
	parser.add_argument("--niter", type=int, help="number of iteration",default=30)
	parser.add_argument("--learnrate", type=float, help="learning rate",default=1e-3)
	parser.add_argument("--rigidbody", action="store_true", default=False ,help="for rigid body motion decoder")
	parser.add_argument("--nogoldstandard", action="store_true", default=False ,help="merge even/odd for heterogeneity analysis")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	
	logid=E2init(sys.argv)
	path=options.path
	itr=options.iterid
	if options.newiter<0:
		itnew=itr+1
	else:
		itnew=options.newiter
	
	if options.nogoldstandard:
		eoiter=[""]
	else:
		eoiter=["_even", "_odd"]
	for eo in eoiter:
		print(f"Realigning {eo}...")
		pts=np.loadtxt(f"{path}/model_{itr:02d}{eo}.txt")

		imsk=make_mask_gmm(options.mask, pts).numpy()
		
		if itr==0:
			decname=f"{path}/dec{eo}.h5"
		else:
			decname=f"{path}/dec_{itr:02d}{eo}.h5"
			
		decode_model=tf.keras.models.load_model(decname,compile=False,custom_objects={"ResidueConv2D":ResidueConv2D})
		midraw=np.loadtxt(f"{path}/midall_{itr:02d}{eo}.txt")[:,1:]
		
		pfile=f"{path}/ptcls_{itr:02d}{eo}.lst"
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
		
		nbatch=len(xfsnp)//bsz
		niter=options.niter
		xfvs=[]
		scores=[]
		for conf,xf in trainset:
			xfvar=tf.Variable(xf)
			opt=tf.keras.optimizers.Adam(learning_rate=options.learnrate) 
			p0=tf.constant(tf.zeros((xf.shape[0],pts.shape[0], 3))+pts[:,:3])
			
			if options.rigidbody:
				pout=decode_model(conf)
				pout=rotpts(p0[:,:,:3], pout, imsk) 
				
			else:
				pout=decode_model(conf)[:,:,:3]+p0
				
			pj0=xf2pts_mult(pout, xf)

			cost=[]
			for it in range(niter):
				with tf.GradientTape() as gt:
					pj1=xf2pts_mult(p0, xfvar)
					dpj=pj1-pj0
					dpj=tf.reduce_sum(dpj**2, axis=2)*imsk
					# print(dpj.shape)
					# loss1=tf.reduce_sum(dpj, axis=1)
					# loss=tf.reduce_mean(loss1)
					loss1=tf.reduce_mean(dpj, axis=1)
					loss=tf.reduce_mean(tf.math.sqrt(loss1))
					loss*=rawbox*e["apix_x"]
					if it==0: loss0=loss


				grad=gt.gradient(loss, xfvar)
				opt.apply_gradients([(grad, xfvar)])
				cost.append(loss)

			sys.stdout.write("\r batch {}/{}, RMSD loss {:.3f} -> {:.3f} A  ".format(len(xfvs), nbatch, loss0, loss))
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

		oname=f"{path}/ptcls_{itnew:02d}{eo}.lst"
		print("saving aligned particles to {}".format(oname))
		lstin=load_lst_params(f"{path}/ptcls_{itr:02d}{eo}.lst")
		#frcs=[s["score"] for s in lstin]
		for i,xf in enumerate(xfs):
			lstin[i]["xform.projection"]=xf
			lstin[i]["score"]=score[i]#frcs[i]

		if os.path.isfile(oname): os.remove(oname)
		save_lst_params(lstin, oname)

		if eo=="":
			for e in ["even", "odd"]:
				run(f"e2spa_make3d.py --input {path}/ptcls_{itnew:02d}.lst --output {path}/threed_{itnew:02d}_{e}.hdf --clsid {e} --keep .9 --parallel thread:32 --sym c1")
				run(f"e2proc3d.py {path}/threed_{itnew:02d}_{e}.hdf {path}/threed_raw_{e}.hdf")
		else:
			run(f"e2spa_make3d.py --input {path}/ptcls_{itnew:02d}{eo}.lst --output {path}/threed_{itnew:02d}{eo}.hdf --keep .9 --parallel thread:32 --sym c1")
			run(f"e2proc3d.py {path}/threed_{itnew:02d}{eo}.hdf {path}/threed_raw{eo}.hdf")
		
		
	E2end(logid)
	
	
if __name__ == '__main__':
	main()
	
