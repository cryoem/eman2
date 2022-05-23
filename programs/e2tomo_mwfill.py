#!/usr/bin/env python
# Muyuan Chen 2020-03
from EMAN2 import *
from EMAN2_utils import *
import numpy as np

if "CUDA_VISIBLE_DEVICES" not in os.environ:
	# so we can decide which gpu to use with environmental variable
	os.environ["CUDA_VISIBLE_DEVICES"]='0' 
	
os.environ["TF_FORCE_GPU_ALLOW_GROWTH"]='true' 
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' #### reduce log output
if ('-h' in sys.argv) or ('--help' in sys.argv):
	tf=type('empty', (object,), {})()
	def empty():
		return lambda f: f
	tf.function=empty
	print(tf.function())
	print("Printing help. Skip tensorflow import")
else:
	import tensorflow as tf

def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--train", type=str,help="train on tomo", default=None)
	parser.add_argument("--boxsz", type=int,help="box size", default=128)
	#parser.add_argument("--gpuid", type=int,help="gpu id", default=None)
	parser.add_argument("--nsample", type=int,help="number of samples", default=2000)
	parser.add_argument("--learnrate", type=float,help="learning rate", default=2e-4)
	parser.add_argument("--load", type=str,help="load model", default=None)
	parser.add_argument("--applyto", type=str,help="apply to tomo", default=None)
	
	
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)

	if not os.path.isdir("neuralnets"):
		os.mkdir("neuralnets")
	
	#global tf
	#os.environ["TF_FORCE_GPU_ALLOW_GROWTH"]='true' 
	#tf=import_tensorflow(options.gpuid)
	
	if options.train!=None:
		gen_y2x=do_training(options)
		msave='neuralnets/mwfill2d_modelsave.h5'
		gen_y2x.save(msave)
		print("model saved to {}".format(msave))

	if options.load!=None:
		gen_y2x=tf.keras.models.load_model(options.load, compile=False, custom_objects={"FFTLayer": FFTLayer})
		
	if options.applyto!=None:
		fnames=options.applyto.split(',')
		for f in fnames:
			do_apply(f, gen_y2x, options)

	
		
	E2end(logid)

def do_apply(fname, model, options):
	t0=time.time()
	tomo=EMData(fname)
	sz=options.boxsz
	hdr=tomo.get_attr_dict()
	tomo=tomo.numpy().copy()
	tomo=np.clip(tomo/5., -1, 1)
	shp=np.array(tomo.shape)
	step=64
	shpout=[shp[0]//step, shp[1], shp[2]//step]

	coord=np.indices(shpout)
	coord=coord.reshape((3,-1)).T

	coord[:,2]=coord[:,2]*step+step//2
	coord[:,0]=(coord[:,0]+1)*sz//2
	coord+=sz//2


	testy=[]
	tomopad=np.pad(tomo, sz//2)
	print("Using {} tiles".format(len(coord)))
	for c in coord:
		img=tomopad[c[0]-sz//2:c[0]+sz//2, c[1], c[2]-sz//2:c[2]+sz//2].copy()
		testy.append(img)
		
	testy=np.array(testy).reshape((-1, sz,sz,1))

	y2xs=[]
	bsz=1000
	for i in range(len(coord)//bsz+1):
		y2x=func_iter(model, testy[i*bsz:(i+1)*bsz])
		y2xs.extend(y2x.numpy())
		
	y2x=np.array(y2xs)

	tomoout=np.zeros_like(tomopad)
	tomowt=np.zeros_like(tomopad)
	for i,c in enumerate(coord):
		tomoout[c[0]-sz//2:c[0]+sz//2, c[1], c[2]-sz//2:c[2]+sz//2]+=y2x[i,:,:,0]
		tomowt[c[0]-sz//2:c[0]+sz//2, c[1], c[2]-sz//2:c[2]+sz//2]+=1

	tomowt[tomowt==0]=1
	tomoout/=tomowt
	tomoout=tomoout[sz//2:-sz//2,sz//2:-sz//2,sz//2:-sz//2].copy()
	e=from_numpy(tomoout)
	e.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.5})
	
	outname=fname[:-4]+"_mw.hdf"
	e.set_attr_dict(hdr)
	e.write_image(outname)
	print("Output written to {}  ({:.1f}s)".format(outname, time.time()-t0))

### train generator
def do_training(options):
	
	tomo=EMData(options.train)
	sz=options.boxsz
		
	print("Selecting samples...")
	tomobin=tomo.process("math.meanshrink",{"n":4})
	tomobin.mult(-1)
	tomobin.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.45})
	tomobin.process_inplace("mask.onlypeaks")
	tomobin=tomobin.numpy().copy()
	tomo=tomo.numpy().copy()
	tomo=np.clip(tomo/5., -1, 1)
	
	shp=np.array(tomo.shape)
	
	tb=tomobin.flatten()
	tb=tb[tb>0]
	thr=-np.sort(-tb)[len(tb)//8]
	pts=np.array(np.where(tomobin>thr)).T
	pts*=4
	pts=pts[np.all(pts>sz//3, axis=1),:]
	pts=pts[np.all(pts<(shp-sz//3), axis=1),:]
	pts+=sz

	tomopad=np.pad(tomo, sz)

	nsample=options.nsample
	coord=np.random.randint(0,len(pts)-1, nsample)
	coord=pts[coord,:]
	
	### introduce some translational invarient
	coord0=np.vstack([coord+np.array(c)*sz//2 
		for c in [[0,0,0],[0,1,0],[0,-1,0],[0,0,-1],[0,0,1]]])
	coord1=np.vstack([coord+np.array(c)*sz//2 
		for c in [[0,0,0],[-1,0,0],[1,0,0],[0,0,-1],[0,0,1]]])
	
	trainx=[]
	for c in coord0:
		img=tomopad[c[0], c[1]-sz//2:c[1]+sz//2, c[2]-sz//2:c[2]+sz//2].copy()
		trainx.append(img)
		
	trainx=np.array(trainx).reshape((-1, sz,sz,1))

	trainy=[]
	for c in coord1:
		img=tomopad[c[0]-sz//2:c[0]+sz//2, c[1], c[2]-sz//2:c[2]+sz//2].copy()
		trainy.append(img)
		
	trainy=np.array(trainy).reshape((-1, sz,sz,1))

	### define missing wedge 
	idx=np.indices((sz,sz//2+1))
	idx=np.minimum(idx, sz-idx)
	rr=np.sqrt(np.sum(idx**2, axis=0))

	xft=tf.signal.rfft2d(trainy[:10,:,:,0])
	xft=np.mean(abs(xft), axis=0)
	rrint=np.round(rr).astype(int)
	nrm=np.bincount(rrint.flatten())

	tbin=np.bincount(rrint.flatten(), xft.flatten())
	sf=tbin/nrm
	sf[:5]=.01
	sf[sz//3:]=1e10
	xft=xft/sf[rrint]
	wedge=(xft>1.2).astype(np.float32)

	gen_x2y=model_addmw((sz,sz,1), wedge, noise=False)
	gen_y2x=model_genr((sz,sz,1))
	dis_x=model_dis((sz,sz,1))

	### tensorflow training functions
	lossabs=tf.losses.MeanAbsoluteError()
	lossbce=tf.keras.losses.BinaryCrossentropy()
	opt_dis=tf.keras.optimizers.Adam(learning_rate=1e-4) 
	opt_gen=tf.keras.optimizers.Adam(learning_rate=2e-4) 

	def train_gen(x,y):
		with tf.GradientTape() as gt:
			x2y=gen_x2y(x)
			x2y2x=func_iter(gen_y2x, x2y)
			l_gen=lossabs(x, x2y2x)
			
		var_gen=gen_y2x.trainable_variables
		gd_gen=gt.gradient(l_gen, var_gen)
		opt_gen.apply_gradients(zip(gd_gen, var_gen))
		
		return l_gen

	def train_dis(x,y):
		with tf.GradientTape() as gt:
			dx=dis_x(x)
			dxy=dis_x(y)
			l_dx=lossbce(tf.ones_like(dx), dx)
			l_dxy=lossbce(tf.zeros_like(dxy), dxy)
			l_dis=l_dx+l_dxy
				
		var_dis=dis_x.trainable_variables
		gd_dis=gt.gradient(l_dis, var_dis)
		opt_dis.apply_gradients(zip(gd_dis, var_dis))
		
		return l_dis
	
	def test_disx(y):
		y2x=func_iter(gen_y2x, y)
		dx=dis_x(y2x)
		return np.mean(dx)

	### setup training dataset
	batchsize=64
	trainset=tf.data.Dataset.from_tensor_slices((trainx, trainy))
	trainset=trainset.shuffle(500).batch(batchsize)
	nbatch=len(trainx)//batchsize
	
	print("Training discriminator...")
	for itr in range(3):
		cost=[]
		for x,y in trainset:
			l_dis=train_dis(x,y)
			cost.append(l_dis)
			sys.stdout.write("\r {}/{}    {:.3f}".format(len(cost), nbatch, l_dis))
			sys.stdout.flush()
			
		sys.stdout.write("\r")
		print("iter {}, loss : {:.3f}".format(itr, np.mean(cost)))
		

	print("Training generator...")
	vloss=[]
	for itr in range(50):
		cost=[]
		for x,y in trainset:
			l_gen=train_gen(x,y)
			cost.append(l_gen)
			sys.stdout.write("\r {}/{}    {:.3f}".format(len(cost), nbatch, l_gen))
			sys.stdout.flush()
			
		sys.stdout.write("\r")
		vloss.append(test_disx(trainy[:nsample//10]))
		print("iter {}, loss train: {:.3f}, valid : {:.3f}".format(itr, np.mean(cost), vloss[-1]))
		if itr>5 and (vloss[-1]>.95 or (vloss[-3]>vloss[-2] and vloss[-3]>vloss[-1])):
			#print("Valid curve flattened. Stop.")
			break
		
	print("Generating sample outputs...")
	outname="neuralnets/trainout_mwfill2d.hdf"
	if os.path.isfile(outname):
		os.remove(outname)
	
	y=trainy[:nsample//10]
	y2x=func_iter(gen_y2x, y)[:,:,:,0]
	for i,m in enumerate(y2x):
		
		e=from_numpy(y[i,:,:,0])
		e.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.5})
		e.write_image(outname,i*2)
		
		e=from_numpy(m.numpy())
		e.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.5})
		e.write_image(outname,i*2+1)

	return gen_y2x

def func_iter(func, x, itr=3):
	y=x
	for i in range(itr):
		y=func(y)
	return y

class FFTLayer(tf.keras.layers.Layer):

	def __init__(self, inv=0,**kwargs):
		super(FFTLayer, self).__init__()
		self.inv=inv

	def build(self,pm): 
		pass
	
	def get_config(self):
		config = super().get_config()
		config.update({
			"inv": self.inv,
		})
		return config

	def call(self, inputs):  
		if self.inv:
			return tf.signal.irfft2d(inputs)
		else:
			return tf.signal.rfft2d(inputs)

### return keras model that adds missing wedge to input
def model_addmw(input_shape, wedge, noise=False):
	
	x0=tf.keras.Input(shape=input_shape)
	sz=input_shape[0]
	if noise:
		ns=tf.zeros_like(x0)
		ns=tf.keras.layers.GaussianNoise(.05)(ns, training=True)
		x1=x0+ns
	else:
		x1=x0
	
	xft=FFTLayer(0)(x1[:,:,:,0])
	#xft=tf.signal.rfft2d(x1[:,:,:,0])
	xft=xft*wedge
	y0=FFTLayer(1)(xft)
	#y0=tf.signal.irfft2d(xft)
	y0=tf.reshape(y0, (-1,sz,sz,1))
	model=tf.keras.Model(x0, y0)
	return model
	
### return the generator model
def model_genr(input_shape):
	nbin=2
	x0=tf.keras.Input(shape=input_shape)
	sz=input_shape[0]
	
	xp0=tf.keras.layers.MaxPool2D((nbin,nbin))(x0)
	xft=FFTLayer(0)(xp0[:,:,:,0])
	#xft=tf.signal.rfft2d(xp0[:,:,:,0])
	xft=tf.keras.layers.Reshape((sz//nbin,sz//(nbin*2)+1,1))(xft)
	
	xft1=tf.concat([tf.math.real(xft), tf.math.imag(xft)], axis=-1)
	xft1=tf.keras.layers.Flatten()(xft1)
	

	xft2=tf.keras.layers.Dense(1024, activation='relu')(xft1)
	xft2=tf.keras.layers.Dense(1024, activation='relu')(xft1)
	nmid=(sz//nbin)*(sz//(nbin*2)+1)*2
	xft3=tf.keras.layers.Dense(4*nmid)(xft2)
	
	xft3=tf.keras.layers.Reshape((2,4,sz//nbin,sz//(nbin*2)+1))(xft3)
	xft3=tf.complex(xft3[:,0], xft3[:,1])
	
	#xft4=tf.signal.irfft2d(xft3)
	xft4=FFTLayer(1)(xft3)
	xft4=tf.transpose(xft4, (0,2,3,1))

	
	x1=tf.keras.layers.Conv2D(16, 5, activation="relu", strides=(2,2), padding="same", use_bias=False)(x0)
	x2=tf.keras.layers.Conv2D(16, 3, activation="relu", strides=(2,2), padding="same", use_bias=False)(x1)
	x3=tf.keras.layers.Conv2D(16, 3, activation="relu", strides=(2,2), padding="same", use_bias=False)(x2)

	x4=tf.keras.layers.Conv2D(16, 3, activation="relu", padding="same", use_bias=False)(x3)
	x4=tf.keras.layers.Concatenate()([x4, x3 ])
	
	y2=tf.keras.layers.Conv2DTranspose(16, 3, activation="relu", strides=(2,2), padding="same", use_bias=False)(x4)
	y2=tf.keras.layers.Concatenate()([x2, y2, ])

	y1=tf.keras.layers.Conv2DTranspose(16, 3, activation="relu", strides=(2,2), padding="same", use_bias=False)(y2)
	y1=tf.keras.layers.Concatenate()([x1, y1,xft4])
	y1=tf.keras.layers.Dropout(.2)(y1)
	
	y0=tf.keras.layers.Conv2DTranspose(16, 5, activation="relu", strides=(2,2), padding="same", use_bias=False)(y1)
	y0=tf.keras.layers.BatchNormalization()(y0)
	
	y0=tf.keras.layers.Conv2DTranspose(1, 5, activation="tanh",  padding="same")(y0)

	model=tf.keras.Model(x0, y0)
	return model

### return the discriminator model
def model_dis(input_shape):
	
	x0=tf.keras.Input(shape=input_shape)
	layers=[
		tf.keras.layers.Conv2D(16, 5, activation="relu", strides=(2,2), padding="valid"),
		tf.keras.layers.Conv2D(16, 3, activation="relu", strides=(2,2), padding="valid"),
		tf.keras.layers.Conv2D(16, 3, activation="relu", strides=(2,2), padding="valid"),
		tf.keras.layers.Flatten(),
		tf.keras.layers.Dense(128, activation="relu"),
		tf.keras.layers.BatchNormalization()
	]
	
	y0=x0
	for l in layers:
		y0=l(y0)
	
	x1=tf.keras.layers.MaxPool2D()(x0)
	xft=tf.abs(FFTLayer(0)(x1[:,:,:,0]))
	y1=tf.keras.layers.Flatten()(xft)
	y1=tf.keras.layers.Dense(128, activation="relu")(y1)
	y1=tf.keras.layers.BatchNormalization()(y1)
	
	y2=tf.keras.layers.Concatenate()([y0,y1])
	y3=tf.keras.layers.Dense(1, activation="sigmoid")(y2)

	model=tf.keras.Model(x0, y3)
	return model

if __name__ == '__main__':
	main()
	
