#!/usr/bin/env python
# Steve Ludtke 2026
#
from EMAN3 import *
from EMAN3jax import *
import jax
import jax.numpy as jnp
import jax.random as random
from flax import nnx
from flax import nnx as nnx_flax
from flax.training import train_state
import optax
import os
import traceback
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
#	parser.add_argument("--net_style",type=str,help="Multiple network designs are available: leaky_5, relu_3, linear",default="leaky_5")
#	parser.add_argument("--learnrate", type=float,help="learning rate for model training only. Default is 1e-4. ", default=1e-4)
#	parser.add_argument("--sigmareg", type=float,help="regularizer for the sigma of gaussian width. Larger value means all Gaussian functions will have essentially the same width. Smaller value may help compensating local resolution difference.", default=.5)
#	parser.add_argument("--modelreg", type=float,help="regularizer for for Gaussian positions based on the starting model, ie the result will be biased towards the starting model when training the decoder (0-1 typ). Default 0", default=0)
#	parser.add_argument("--ampreg", type=float,help="regularizer for  Gaussian amplitudes. Large values will encourage all Gaussians towards 1.0 or -0.2. default = 0", default=0)
#	parser.add_argument("--niter", type=int,help="number of iterations", default=32)
#	parser.add_argument("--npts", type=int,help="number of points to initialize. ", default=-1)
#	parser.add_argument("--batchsz", type=int,help="batch size", default=192)
#	parser.add_argument("--minressz", type=int,help="Fourier diameter associated with minimum resolution to consider. ", default=4)
#	parser.add_argument("--maxboxsz", type=int,help="maximum fourier box size to use. 2 x target Fourier radius. ", default=64)
#	parser.add_argument("--maxres", type=float,help="maximum resolution. will overwrite maxboxsz. ", default=-1)
#	parser.add_argument("--align", action="store_true", default=False ,help="align particles.")
#	parser.add_argument("--heter", action="store_true", default=False ,help="heterogeneity analysis.")
#	parser.add_argument("--decoderentropy", action="store_true", default=False ,help="This will train some entropy into the decoder using particles to reduce vanishing gradient problems")
#	parser.add_argument("--perturb", type=float, default=0.1 ,help="Relative perturbation level to apply in each iteration during --heter training. Default = 0.1, decrease if models are too disordered")
#	parser.add_argument("--conv", action="store_true", default=False ,help="Use a convolutional network for heterogeneity analysis.")
#	parser.add_argument("--fromscratch", action="store_true", default=False ,help="start from coarse alignment. otherwise will only do refinement from last round")
#	parser.add_argument("--ptclrepout", type=str,help="Save the per-particle representation input to the network to a file for later use", default="")
#	parser.add_argument("--ptclrepin", type=str,help="Load the per-particle representation from a file rather than recomputing", default="")
#	parser.add_argument("--midout", type=str,help="middle layer output", default="")
#	parser.add_argument("--pas", type=str,help="choose whether to adjust position, amplitude, sigma. sigma is not supported in this version of the program. use 3 digit 0/1 input. default is 110, i.e. only adjusting position and amplitude", default="110")
#	parser.add_argument("--nmid", type=int,help="size of the middle layer. If model is grouped must be divisible by ngroup", default=4)
#	parser.add_argument("--mask", type=str,help="remove points outside mask", default="")
	parser.add_argument("--gpudev",type=int,help="GPU Device, default 0", default=0)
	parser.add_argument("--gpuram",type=int,help="Maximum GPU ram to allocate in MB, default=7168", default=7168)		# default should run on 8G cards, but probably not a good idea
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verbosity")

	(options, args) = parser.parse_args()

	logid=E3init(sys.argv,options.ppid)

	# This gives us access to all of the particle data and metadata
	cache=StackCache(options.ptcls)

	# the starting model, which provides the origin in latent space
	points=Gaussians(options.model)

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

	# Determine per-particle input representations
	for bn in range(0,len(cache),4096):
		ptcls=cache.read(64,range(bn,min(bn+4096,len(cache))))
		meta=ptcls.metadata
		mx2d=Orientations(meta[:,2:5]).to_mx2d(swapxy=True)
		tytx=jnp.array(meta[:,0:2])
		frc,grad=prj_frc_loss(points.jax,mx2d,tytx,ptcls.jax,weights[64],threshs[64])
		print(frc.shape,grad.shape,frc.min(),frc.max())

 #    # Initialize encoder
 #    encoder = Encoder(input_dim=input_dim, latent_dim=latent_dim, hidden_sizes=hidden_sizes)
 #    encoder_params = encoder.init(random.PRNGKey(0), random.normal(random.PRNGKey(1), (batch_size, input_dim)))
 #
 #    encoder_optimizer = optax.adam(learning_rate=1e-3)
 #    encoder_train_state = TrainState.create(
 #        apply_fn=encoder.apply,
 #        params=encoder_params,
 #        tx=encoder_optimizer
 #    )
 #
 #    # Initialize decoder
 #    decoder = Decoder(latent_dim=latent_dim, output_dim=output_dim, hidden_sizes=hidden_sizes)
 #    decoder_params = decoder.init(random.PRNGKey(0), random.normal(random.PRNGKey(1), (batch_size, latent_dim)))
 #
 #    decoder_optimizer = optax.adam(learning_rate=1e-3)
 #    decoder_train_state = TrainState.create(
 #        apply_fn=decoder.apply,
 #        params=decoder_params,
 #        tx=decoder_optimizer
 #    )
 #
 #
 #
	# latent=encoder_state.apply_fn(encoder_state.params, input_vector)
#	output=decoder_state.apply_fn(decoder_state.params, latent_vector)

	E3end(logid)

@jit
@jax.value_and_grad
def prj_frc_loss(points:jnp.array,mx2d,tytx,ptcls,weight,thresh):
	"""Aggregates the functions we need to calculate the gradient through. Computes the frc array resulting from the
	comparison of the Gaussians in gaus to particles in known orientations. Returns -frc since optax wants to minimize, not maximize"""

	ny=ptcls.shape[1]
	prj=gauss_project_simple_fn(points,mx2d,ny,tytx)
	return -jax_frc_jit_new(jax_fft2d(prj),ptcls,weight,thresh)

def prj_frcs(points:jnp.array,ptcls:jnp.array,meta:jnp.array):
	"""Computes the FRC between a 3-D model and a stack of projections. Instead of integrating to produce
	a loss function, this returns the individual FRC curves for statistical analysis"""
	mx2d=Orientations(meta[:,2:5]).to_mx2d(swapxy=True)
	ny=ptcls.shape[1]
	prjf=jax_fft2d_jit(gauss_project_simple_fn(points,mx2d,ny,meta[:,0:2]))
#	print(prjf.shape,ptcls.jax.shape)

	return jax_frcs_jit(prjf,ptcls)

class Encoder(nnx.Module):
    def __init__(self, input_dim, latent_dim, hidden_sizes, rngs):
        """
        Initialize the encoder.

        Args:
            input_dim: Dimension of input vectors
            latent_dim: Dimension of latent space
            hidden_sizes: Tuple of hidden layer sizes (one layer per size)
            rngs: NNX random keys (dropout, etc.)
        """
        layers = []
        current_dim = input_dim

        for size in hidden_sizes:
            layers.append(nnx_flax.Dense(current_dim, size))
            layers.append(nnx.relu)  # Activation
            current_dim = size

        # Final layer to latent dimension
        layers.append(nnx_flax.Dense(current_dim, latent_dim))

        self.layers = nnx.Sequential(layers)

    def __call__(self, x):
        return self.layers(x)


class Decoder(nnx.Module):
    def __init__(self, latent_dim, output_dim, hidden_sizes, rngs):
        """
        Initialize the decoder.

        Args:
            latent_dim: Dimension of latent space
            output_dim: Dimension of output vectors
            hidden_sizes: Tuple of hidden layer sizes (reversed from encoder)
            rngs: NNX random keys
        """
        layers = []
        current_dim = latent_dim

        for size in hidden_sizes:
            layers.append(nnx_flax.Dense(current_dim, size))
            layers.append(nnx.relu)  # Activation
            current_dim = size

        # Final layer to output dimension
        layers.append(nnx_flax.Dense(current_dim, output_dim))

        self.layers = nnx.Sequential(layers)

    def __call__(self, z):
        return self.layers(z)


"""
JAX Autoencoder with NNX (NumPy Neural Networks) API.
Features:
- Separate encoder and decoder with NNX Modules
- Two-phase training: decoder-only, then encoder-decoder
- Model save/load with NNX serialization
"""


# ============================================
# USER-PROVIDED FUNCTIONS (YOU WILL SUPPLY THESE)
# ============================================

def generate_input_vector(key):
    """
    Generate a single input vector of shape (N,)
    Replace with your actual implementation.

    Args:
        key: JAX random key

    Returns:
        Input vector of shape (N,)
    """
    N = 100  # Set your input dimension
    return random.normal(key, (N,))


def generate_output_vector(key, input_vec):
    """
    Generate a target output vector of shape (N, 4) given input vector.
    Replace with your actual implementation.

    Args:
        key: JAX random key
        input_vec: Input vector of shape (N,)

    Returns:
        Output vector of shape (N, 4)
    """
    N = 100
    output_dim = 4
    return random.normal(key, (N, output_dim))


def compute_loss(decoder_output, target_output):
    """
    Compute loss between decoder output and target output.
    Replace with your actual implementation.

    Args:
        decoder_output: Output from decoder
        target_output: Target output

    Returns:
        Scalar loss value
    """
    # Example: Mean squared error
    return jnp.mean((decoder_output - target_output) ** 2)




# ============================================
# TRAINING FUNCTIONS
# ============================================

# def train_decoder(
#     state: dict,
#     decoder: nnx.Module,
#     key,
#     latent_vectors,
#     output_vectors,
#     num_epochs: int = 1000,
#     batch_size: int = 32,
#     learning_rate: float = 1e-3
# ):
#     """
#     Train decoder alone with known latent vectors and outputs.
#
#     Args:
#         state: Training state dictionary
#         decoder: Decoder module instance
#         key: JAX random key
#         latent_vectors: Array of shape (num_samples, latent_dim)
#         output_vectors: Array of shape (num_samples, output_dim)
#         num_epochs: Number of training epochs
#         batch_size: Batch size for training
#         learning_rate: Learning rate for optimizer
#
#     Returns:
#         Updated training state
#     """
#
#     # JIT-compiled training step
#     @nnx.transformed(jit)
#     def train_step(decoder, state, batch_latent, batch_output):
#         def loss_fn(params):
#             preds = decoder(params, batch_latent)
#             loss = compute_loss(preds, batch_output)
#             return loss
#
#         loss, grads = nnx.grad(loss_fn, has_aux=False)(decoder, state)
#
#         # Update optimizer
#         nnx.optimizer_step(state, grads)
#
#         return loss, decoder
#
#     # JIT-compiled evaluation
#     @nnx.transformed(jit)
#     def evaluate(decoder, params, latent_batch, output_batch):
#         preds = decoder(params, latent_batch)
#         mse = compute_loss(preds, output_batch)
#         return mse
#
#     # Training loop
#     num_samples = latent_vectors.shape[0]
#     indices = jnp.arange(num_samples)
#
#     for epoch in range(num_epochs):
#         # Shuffle data
#         perm = random.shuffle(key, indices)
#         latent_shuffled = latent_vectors[perm]
#         output_shuffled = output_vectors[perm]
#
#         epoch_losses = []
#
#         # Mini-batch training
#         for i in range(0, num_samples, batch_size):
#             start = i
#             end = min(i + batch_size, num_samples)
#             batch_latent = latent_shuffled[start:end]
#             batch_output = output_shuffled[start:end]
#
#             # Get current parameters
#             params = nnx.get_model_state(decoder)
#             loss, decoder = train_step(decoder, state, batch_latent, batch_output)
#             epoch_losses.append(loss)
#
#         if (epoch + 1) % 100 == 0:
#             avg_loss = jnp.mean(jnp.array(epoch_losses))
#             print(f"Decoder Training - Epoch {epoch+1}/{num_epochs}, Loss: {avg_loss:.6f}")
#
#     return state
#
#
# def train_encoder_decoder(
#     encoder_state: dict,
#     decoder_state: dict,
#     encoder: nnx.Module,
#     decoder: nnx.Module,
#     key,
#     input_vectors,
#     num_epochs: int = 1000,
#     batch_size: int = 32,
#     learning_rate: float = 1e-3,
#     latent_dim: int = 10
# ):
#     """
#     Train complete encoder-decoder system end-to-end.
#
#     Args:
#         encoder_state: Training state for encoder
#         decoder_state: Training state for decoder
#         encoder: Encoder module instance
#         decoder: Decoder module instance
#         key: JAX random key
#         input_vectors: Array of shape (num_samples, input_dim)
#         num_epochs: Number of training epochs
#         batch_size: Batch size for training
#         learning_rate: Learning rate for optimizers
#         latent_dim: Expected latent vector dimension
#
#     Returns:
#         Tuple of (updated encoder_state, updated decoder_state)
#     """
#
#     @nnx.transformed(jit)
#     def end_to_end_train_step(
#         encoder,
#         encoder_state,
#         decoder,
#         decoder_state,
#         batch_input,
#         key
#     ):
#         def total_loss_fn(enc_params, dec_params):
#             # Encode
#             z = encoder(enc_params, batch_input)
#
#             # Decode
#             decoder_output = decoder(dec_params, z)
#
#             # Generate targets
#             target_outputs = vmap(generate_output_vector)(
#                 random.split(key, batch_input.shape[0]),
#                 batch_input
#             )
#
#             # Compute loss
#             loss = compute_loss(decoder_output, target_outputs)
#             return loss, z
#
#         (loss, z), encoder_grads, decoder_grads = nnx.value_and_grad(
#             total_loss_fn, has_aux=True
#         )(encoder, encoder_state, decoder, decoder_state, batch_input, key)
#
#         # Update both optimizers
#         nnx.optimizer_step(encoder_state, encoder_grads)
#         nnx.optimizer_step(decoder_state, decoder_grads)
#
#         return loss, z
#
#     # Training loop
#     num_samples = input_vectors.shape[0]
#     indices = jnp.arange(num_samples)
#
#     for epoch in range(num_epochs):
#         # Shuffle data
#         perm = random.shuffle(key, indices)
#         input_shuffled = input_vectors[perm]
#
#         epoch_losses = []
#
#         # Mini-batch training
#         for i in range(0, num_samples, batch_size):
#             start = i
#             end = min(i + batch_size, num_samples)
#             batch_input = input_shuffled[start:end]
#
#             loss, z = end_to_end_train_step(
#                 encoder, encoder_state,
#                 decoder, decoder_state,
#                 batch_input, key
#             )
#             epoch_losses.append(loss)
#
#         if (epoch + 1) % 100 == 0:
#             avg_loss = jnp.mean(jnp.array(epoch_losses))
#             print(f"End-to-End Training - Epoch {epoch+1}/{num_epochs}, Loss: {avg_loss:.6f}")
#
#     return encoder_state, decoder_state
#
#
# # ============================================
# # INFERENCE FUNCTIONS
# # ============================================
#
# def encode_input(encoder, encoder_state, input_vector):
#     """
#     Encode an input vector to latent space.
#
#     Args:
#         encoder: Encoder module instance
#         encoder_state: Encoder training state (for parameters)
#         input_vector: Input vector of shape (N,)
#
#     Returns:
#         Latent vector of shape (latent_dim,)
#     """
#     params = nnx.get_model_state(encoder)
#     return encoder(params, jnp.expand_dims(input_vector, 0))
#
#
# def decode_latent(decoder, decoder_state, latent_vector):
#     """
#     Decode a latent vector to output space.
#
#     Args:
#         decoder: Decoder module instance
#         decoder_state: Decoder training state (for parameters)
#         latent_vector: Latent vector of shape (latent_dim,)
#
#     Returns:
#         Output vector of shape (N, 4)
#     """
#     params = nnx.get_model_state(decoder)
#     return decoder(params, jnp.expand_dims(latent_vector, 0))
#
#
# # ============================================
# # MODEL SAVE/LOAD (NNX)
# # ============================================
#
# def save_model_to_disk(
#     encoder: nnx.Module,
#     encoder_state: dict,
#     decoder: nnx.Module,
#     decoder_state: dict,
#     file_path: str,
#     metadata: dict = None
# ):
#     """
#     Save the complete autoencoder to disk.
#
#     Args:
#         encoder: Encoder module
#         encoder_state: Encoder training state
#         decoder: Decoder module
#         decoder_state: Decoder training state
#         file_path: Path to save the model
#         metadata: Optional metadata dictionary
#     """
#     os.makedirs(os.path.dirname(file_path) or '.', exist_ok=True)
#
#     # Create save dict
#     save_dict = {
#         'encoder': encoder,
#         'encoder_state': encoder_state,
#         'decoder': decoder,
#         'decoder_state': decoder_state,
#     }
#
#     if metadata:
#         save_dict['metadata'] = metadata
#
#     # Serialize to msgpack
#     serialized = nnx.nxpack_save(save_dict)
#
#     with open(file_path, 'wb') as f:
#         f.write(serialized)
#
#     print(f"✓ Saved autoencoder to {file_path}")
#
#
# def load_model_from_disk(file_path: str):
#     """
#     Load the autoencoder from disk.
#
#     Args:
#         file_path: Path to model file
#
#     Returns:
#         Tuple of (encoder, encoder_state, decoder, decoder_state, metadata)
#     """
#     with open(file_path, 'rb') as f:
#         data = f.read()
#
#     save_dict = nnx.nxpack_restore(data)
#
#     encoder = save_dict['encoder']
#     encoder_state = save_dict['encoder_state']
#     decoder = save_dict['decoder']
#     decoder_state = save_dict['decoder_state']
#     metadata = save_dict.get('metadata', {})
#
#     return encoder, encoder_state, decoder, decoder_state, metadata
#
#
# # ============================================
# # MAIN EXAMPLE
# # ============================================
#
# def main():
#     """Example usage of the NNX autoencoder."""
#
#     # Configuration - ADJUST THESE VALUES
#     input_dim = 100  # Shape (N,)
#     output_dim = 4   # Shape (N, 4)
#     latent_dim = 10
#     hidden_sizes = (64, 32, 16)
#
#     # Number of samples for training
#     num_samples = 1000
#
#     # Initialize random key
#     key = random.PRNGKey(42)
#
#     # Generate sample data
#     print("Generating training data...")
#     keys = random.split(key, num_samples)
#
#     input_vectors = jnp.stack([generate_input_vector(k) for k in keys])
#     output_vectors = jnp.stack([
#         generate_output_vector(k, input_vectors[i])
#         for i, k in enumerate(keys)
#     ])
#
#     # ========================================
#     # PHASE 1: Decoder-only training
#     # ========================================
#     print(f"\n{'='*50}")
#     print("PHASE 1: Decoder-only training")
#     print(f"{'='*50}")
#
#     # Create random latent vectors for decoder training
#     latent_vectors = random.normal(random.PRNGKey(0), (num_samples, latent_dim))
#
#     # Create optimizer keys
#     keys = nnx.split(random.PRNGKey(0), 3)
#
#     # Initialize decoder
#     rngs = nnx.rngs(*keys[1:])
#     decoder = Decoder(latent_dim=latent_dim, output_dim=output_dim, hidden_sizes=hidden_sizes, rngs=rngs)
#
#     # Create decoder optimizer
#     decoder_optimizer = optax.adam(learning_rate=1e-3)
#     decoder_state = nnx.create_optimizer_state(decoder, decoder_optimizer)
#
#     # Train decoder
#     decoder_state = train_decoder(
#         decoder_state,
#         decoder,
#         random.PRNGKey(0),
#         latent_vectors,
#         output_vectors,
#         num_epochs=1000,
#         batch_size=32
#     )
#
#     # ========================================
#     # PHASE 2: Encoder-Decoder end-to-end training
#     # ========================================
#     print(f"\n{'='*50}")
#     print("PHASE 2: Encoder-Decoder end-to-end training")
#     print(f"{'='*50}")
#
#     # Initialize encoder
#     encoder = Encoder(input_dim=input_dim, latent_dim=latent_dim, hidden_sizes=hidden_sizes, rngs=rngs)
#     encoder_optimizer = optax.adam(learning_rate=1e-3)
#     encoder_state = nnx.create_optimizer_state(encoder, encoder_optimizer)
#
#     # Train encoder and decoder together
#     encoder_state, decoder_state = train_encoder_decoder(
#         encoder_state,
#         decoder_state,
#         encoder,
#         decoder,
#         random.PRNGKey(0),
#         input_vectors,
#         num_epochs=1000,
#         batch_size=32,
#         latent_dim=latent_dim
#     )
#
#     # ========================================
#     # SAVE THE TRAINED MODEL
#     # ========================================
#     print(f"\n{'='*50}")
#     print("SAVING MODEL")
#     print(f"{'='*50}")
#
#     save_model_to_disk(
#         encoder,
#         encoder_state,
#         decoder,
#         decoder_state,
#         save_path='autoencoder_complete.npz',
#         metadata={
#             'input_dim': input_dim,
#             'output_dim': output_dim,
#             'latent_dim': latent_dim,
#             'hidden_sizes': hidden_sizes,
#             'training_epochs': 1000,
#             'batch_size': 32
#         }
#     )
#
#     # ========================================
#     # INFERENCE EXAMPLE
#     # ========================================
#     print(f"\n{'='*50}")
#     print("INFERENCE EXAMPLE")
#     print(f"{'='*50}")
#
#     # Test on new sample
#     test_key = random.PRNGKey(42)
#     test_input = generate_input_vector(test_key)
#     print(f"Input vector shape: {test_input.shape}")
#
#     # Encode
#     latent = encode_input(encoder, encoder_state, test_input)
#     print(f"Latent vector shape: {latent.shape}")
#
#     # Decode
#     reconstructed = decode_latent(decoder, decoder_state, latent)
#     print(f"Reconstructed output shape: {reconstructed.shape}")


# Entry point
if __name__ == "__main__":
    main()

