#!/usr/bin/env python
# Steve Ludtke 2026
#
from EMAN3 import *
from EMAN3jax import *
import jax
import jax.numpy as jnp
import jax.random as random
from flax import nnx
from flax.training import train_state
import optax
import os
import traceback
import time

#from sklearn.decomposition import PCA

try: os.mkdir(".jaxcache")
except: pass

# We cache the JIT compilation results to speed up future runs
jax.config.update("jax_compilation_cache_dir", "./.jaxcache")
jax.config.update("jax_persistent_cache_min_entry_size_bytes", -1)
jax.config.update("jax_persistent_cache_min_compile_time_secs", 2)
jax.config.update("jax_persistent_cache_enable_xla_caches", "xla_gpu_per_fusion_autotune_cache_dir")

jax.config.update("jax_default_matmul_precision", "float32")


def main():

	usage="""Starting with a delta function reconstruction and a stack of oriented particle data, this will attempt to learn the dynamics of the particle population
	by extending the concept of e3make3d_delta.py.

	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--sym", type=str,help="symmetry. currently only support c and d", default="c1")
	parser.add_argument("--model", type=str,help="Required: delta reconstruction (X,Y,Z,A) .txt file. This model will be used for representation of individual particles as well as the 3-D volume.", default=None)
	parser.add_argument("--ptcls", type=str,help="Required: particle data for training. Must be a .lst file with orientations.", default=None)
#	parser.add_argument("--segments", type=str,help="Divide the model into sequential domains. Comma separated list of integers. Each integer is the first sequence number of a new region, starting with 0",default=None)
	parser.add_argument("--decoderin", type=str,help="Rather than initializing the decoder from a model, read an existing trained decoder", default="")
	parser.add_argument("--decoderout", type=str,help="Save the trained decoder model. Filename should be .h5", default=None)
	parser.add_argument("--encoderin", type=str,help="Rather than initializing the encoder from scratch, read an existing trained encoder", default=None)
	parser.add_argument("--encoderout", type=str,help="Save the trained encoder model. Filename should be .h5", default=None)
	parser.add_argument("--ctf", type=int,help="0=no ctf, 1=single ctf, 2=layered ctf",default=0)
	parser.add_argument("--netstyle",type=str,help="Multiple network designs are available: leaky_5, relu_3, linear",default="leaky_5")
#	parser.add_argument("--learnrate", type=float,help="learning rate for model training only. Default is 1e-4. ", default=1e-4)
#	parser.add_argument("--sigmareg", type=float,help="regularizer for the sigma of gaussian width. Larger value means all Gaussian functions will have essentially the same width. Smaller value may help compensating local resolution difference.", default=.5)
#	parser.add_argument("--modelreg", type=float,help="regularizer for for Gaussian positions based on the starting model, ie the result will be biased towards the starting model when training the decoder (0-1 typ). Default 0", default=0)
#	parser.add_argument("--ampreg", type=float,help="regularizer for  Gaussian amplitudes. Large values will encourage all Gaussians towards 1.0 or -0.2. default = 0", default=0)
	parser.add_argument("--niter", type=int,help="number of iterations", default=32)
#	parser.add_argument("--npts", type=int,help="number of points to initialize. ", default=-1)
	parser.add_argument("--batchsz", type=int,help="particle batch size", default=1024)
#	parser.add_argument("--minressz", type=int,help="Fourier diameter associated with minimum resolution to consider. ", default=4)
#	parser.add_argument("--maxboxsz", type=int,help="maximum fourier box size to use. 2 x target Fourier radius. ", default=64)
#	parser.add_argument("--maxres", type=float,help="maximum resolution. will overwrite maxboxsz. ", default=-1)
#	parser.add_argument("--align", action="store_true", default=False ,help="align particles.")
#	parser.add_argument("--heter", action="store_true", default=False ,help="heterogeneity analysis.")
#	parser.add_argument("--decoderentropy", action="store_true", default=False ,help="This will train some entropy into the decoder using particles to reduce vanishing gradient problems")
#	parser.add_argument("--perturb", type=float, default=0.1 ,help="Relative perturbation level to apply in each iteration during --heter training. Default = 0.1, decrease if models are too disordered")
#	parser.add_argument("--conv", action="store_true", default=False ,help="Use a convolutional network for heterogeneity analysis.")
#	parser.add_argument("--fromscratch", action="store_true", default=False ,help="start from coarse alignment. otherwise will only do refinement from last round")
	parser.add_argument("--ptclrep", type=str,help="Save the per-particle representation to a file", default="ptclrep.bin")
	parser.add_argument("--latentout", type=str,help="middle layer output", default="")
#	parser.add_argument("--pas", type=str,help="choose whether to adjust position, amplitude, sigma. sigma is not supported in this version of the program. use 3 digit 0/1 input. default is 110, i.e. only adjusting position and amplitude", default="110")
	parser.add_argument("--nlatent", type=int,help="size of the middle layer. If model is grouped must be divisible by ngroup", default=4)
#	parser.add_argument("--mask", type=str,help="remove points outside mask", default="")
	parser.add_argument("--gpudev",type=int,help="GPU Device, default 0", default=0)
	parser.add_argument("--gpuram",type=int,help="Maximum GPU ram to allocate in MB, default=7168", default=7168)		# default should run on 8G cards, but probably not a good idea
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verbosity")

	(options, args) = parser.parse_args()

	logid=E3init(sys.argv,options.ppid)

	# This gives us access to all of the particle data and metadata
	cache=StackCache(options.ptcls)
	Nptcl=len(cache)			# number of particles, but note that only N*batchsize will be used in training

	# the starting model, which provides the origin in latent space
	points=Points(options.model)

	# critical for later in the program, this initializes the radius images for all of the samplings we will use
	for s in cache.sizes:
		rad_img_int(s)
		if options.ctf>0:
			rad2_img(s)

	# we determine the FSC threshold and weight for each lengthscale
	weights={}
	threshs={}
	for s in cache.sizes:
		ptcls=cache.read(s,range(0,len(cache),len(cache)//1000))	# We use 1000 distributed particles to assess the weight
		frcs=prj_frcs(points.jax,ptcls.jax,jnp.array(ptcls.metadata))
		try:
			thresh=1.25*np.std(frcs,0)/sqrt(len(ptcls))
			weight=1.0/np.array(thresh)		# this should make all of the standard deviations the same
			weight[0:2]=0			# low frequency cutoff
			weight[ptcls.shape[1]//2:]=0
			weight/=np.sum(weight)	# normalize to 1

			weight=jnp.array(weight*len(weight))	# the *len(weight) is dumb, but due to mean() being returned
		except:
			print(f"Weighting failed {s}")
			traceback.print_exc()
			weight=np.ones((len(frcs.shape[1])))
		threshs[s]=thresh
		weights[s]=weight
#		print(f"{s:3d}: {thresh}\n     {weight}")

#	print(weights[128],threshs[128])

	nans=set()
	batchsize=options.batchsz
	# Determine per-particle input representations -> ptclrep
	# grads=np.zeros((batchsize,len(points),4))
	# frcs=np.zeros(batchsize)
	# we recreate the ptclrep every time right now
	frcout=open("frcs.txt","w")
	ptclrep=np.memmap(options.ptclrep,mode="w+",dtype=np.float32,shape=(len(cache),len(points)))
	for bn in range(0,len(cache),batchsize):
		bnend=min(bn+batchsize,len(cache))
		ptcls=cache.read(128,range(bn,bnend))
		meta=ptcls.metadata
		mx2d=Orientations(meta[:,2:5]).to_mx2d(swapxy=True)
		tytx=jnp.array(meta[:,0:2])
		print(points.jax.shape,mx2d.shape,tytx.shape,ptcls.jax.shape)

		frcs,grads=prj_frc_loss_vmap(points.jax,mx2d,tytx,ptcls.jax,weights[128],threshs[128])
		grads*=points.jax.shape[0]
		ptclrep[bn:bnend]=grads[:,:,3]
		for f in frcs: frcout.write(f"{f:1.6f}\n")
		
# 		for n in range(bn,min(bn+1024,len(cache))):
# 			ns=n-bn
# #			print(bn,n,points.jax.shape,mx2d[:,:,ns:ns+1].shape,tytx[ns:ns+1].shape,ptcls.jax[ns:ns+1].shape,ptcls.shape)
# 			frc,grad=prj_frc_loss(points.jax,mx2d[:,:,ns:ns+1],tytx[ns:ns+1],ptcls.jax[ns:ns+1],weights[64],threshs[64])
# 			grads[n]=grad
# 			frcs[n]=frc

#	grads[np.isnan(grads)] = 0
#	frcs[np.isnan(frcs)] = 0
	# np.savetxt("gradsx.txt",grads[:,:,0])
	# np.savetxt("gradsy.txt",grads[:,:,1])
	# np.savetxt("gradsz.txt",grads[:,:,2])
	# np.savetxt("gradsa.txt",grads[:,:,3])
	# np.savetxt("frcs.txt",frcs)

	

	Npnt=ptclrep.shape[1]		# number of points
	netstyles={
		"leaky_5":([Npnt,Npnt//4,Npnt//8,max(Npnt//32,options.nlatent),max(Npnt//64,options.nlatent)],nnx.leaky_relu), 
		"relu_3":([Npnt,Npnt//4,Npnt//16],nnx.relu), 
		"linear":([Npnt],nnx.identity)
	}
	try: hidden,activation=netstyles[options.netstyle]
	except: error_exit(f"ERROR: available network styles are: {",".join(netstyles.keys())}")
	
	rngs = nnx.Rngs(int(time.time()))
	model = Autoencoder(Nin=ptclrep.shape[1], Nout=points.jax.size, Nlat=options.nlatent, hidden_dims=hidden, activation=activation, rngs=rngs)
	
	optimizer = nnx.Optimizer(model, optax.adam(learning_rate=1e-3),wrt=nnx.All)
	
	for epoch in range(options.niter):
		running_loss = 0.0
		
		for batch in range(Nptcl//batchsize):
			x = jnp.array(ptclrep[batch*batchsize:(batch+1)*batchsize])
			aux = cache.read(128,range(batch*batchsize,(batch+1)*batchsize))
			loss_val = train_step(model, optimizer, x, aux.jax, jnp.array(aux.metadata) )
			
			running_loss += float(loss_val)
		
			print(f"Epoch {epoch+1}/{num_epochs} | Loss: {running_loss / num_batches:.4f}")

	print(epoch,running_loss,options.niter,batchsize,Nptcl,Npnt)

	E3end(logid)

@nnx.jit
def train_step(model, optimizer, x, ptcl, meta):
	# 1. Define a local function that computes loss based ONLY on model state.
	#    nnx.value_and_grad differentiates with respect to the first argument's variables.
	def loss_fn(mdl):
		pointary = mdl(x)  # Forward pass
		return(sym_prj_frc_loss_ctf(pointary,meta[:,2:5],))
		return loss_function(predictions, aux)

	# 2. Compute loss value and gradients relative to model state
	loss, grads = nnx.value_and_grad(loss_fn)(model)

	# 3. Update the optimizer state and mutate the model's parameters in place
	optimizer.update(grads)

	return loss


class Encoder(nnx.Module):
	def __init__(self, input_dim: int, latent_dim: int, hidden_dims: list[int], activation, rngs: nnx.Rngs):
		layers = nnx.List()
		current_dim = input_dim

		for size in hidden_dims:
			layers.append(nnx.Linear(in_features=current_dim, out_features=size, rngs=rngs))
			layers.append(activation)
			current_dim = size

		layers.append(nnx.Linear(in_features=current_dim, out_features=latent_dim, rngs=rngs))
		self.layers = nnx.Sequential(*layers)

	def __call__(self, x):
		return self.layers(x)


class Decoder(nnx.Module):
	def __init__(self, latent_dim: int, output_dim: int, hidden_dims: list[int], activation, rngs: nnx.Rngs):
		layers = nnx.List()
		current_dim = latent_dim

		# Reverse the hidden dimensions for symmetric decoding
		for size in reversed(hidden_dims):
			layers.append(nnx.Linear(in_features=current_dim, out_features=size, rngs=rngs))
			layers.append(activation)
			current_dim = size

		# Map back to output dimensionality (Nout)
		layers.append(nnx.Linear(in_features=current_dim, out_features=output_dim, rngs=rngs))
		self.layers = nnx.Sequential(*layers)

	def __call__(self, z):
		return self.layers(z)

class Autoencoder(nnx.Module):
	def __init__(self, Nin: int, Nout: int, Nlat: int, hidden_dims: list[int], activation, rngs: nnx.Rngs):
		# Pass the same rngs container; NNX handles splitting under the hood
		self.encoder = Encoder(input_dim=Nin, latent_dim=Nlat, hidden_dims=hidden_dims, activation=activation,rngs=rngs)
		self.decoder = Decoder(latent_dim=Nlat, output_dim=Nout, hidden_dims=hidden_dims*4, activation=activation,rngs=rngs)

	def __call__(self, x):
		return self.decoder(self.encoder(x))

@jax.value_and_grad
def prj_frc_loss(points:jnp.array,mx2d,tytx,ptcls,weight,thresh):
	"""Aggregates the functions we need to calculate the gradient through. Computes the frc array resulting from the
	comparison of the Gaussians in gaus to one particle in a known orientation. Returns -frc since optax wants to minimize, not maximize"""

	ny=ptcls.shape[0]
	prj=point_project_single_fn(points,mx2d,ny,tytx)
	return -jax_frc_single(jax_fft2d(prj),ptcls,weight,thresh)

prj_frc_loss_vmap = jax.jit(jax.vmap(prj_frc_loss,in_axes=(None, 2, 0, 0, None, None),out_axes=0))

def prj_frcs(points:jnp.array,ptcls:jnp.array,meta:jnp.array):
	"""Computes the FRC between a 3-D model and a stack of projections. Instead of integrating to produce
	a loss function, this returns the individual FRC curves for statistical analysis"""
	mx2d=Orientations(meta[:,2:5]).to_mx2d(swapxy=True)
	ny=ptcls.shape[1]
	prjf=jax_fft2d_jit(point_project_simple_fn(points,mx2d,ny,meta[:,0:2]))
#	print(prjf.shape,ptcls.jax.shape)

	return jax_unweighted_frcs_jit(prjf,ptcls)





# Entry point
if __name__ == "__main__":
	main()

