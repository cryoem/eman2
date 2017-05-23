#!/usr/bin/env python
# Muyuan Chen 2017-03
from EMAN2 import *
import numpy as np
import threading
import Queue
def main():
	
	usage=""" Same as the ConvNet boxer in e2boxer.py, but identifies two classes of particles. Only command line options are available. Note the hyper-parameters may not be optimized...
	
	To train the network: [prog] <ref0> <ref1> <ref_bad> 
	To box particles [prog] <micrographs> --autobox --boxsz # --threads #"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--autobox", action="store_true", default=False ,help="box particles")
	parser.add_argument("--boxsz", type=int,help="box size of input particles. Necessery for applying network", default=256)
	parser.add_argument("--threads", type=int,help="number of threads", default=10)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	
	boxerConvNet.do_import()
	if options.autobox:
		boxerConvNet.do_autobox_all(args, bxsz=options.boxsz, nthreads=options.threads)
		
	else:
		boxerConvNet.do_training(args)
	
	
	E2end(logid)
	
def run(cmd):
	print cmd
	launch_childprocess(cmd)
	


#####################
## Convolutional Neural Network boxer
##########
class boxerConvNet():
	
	
	#### import dependencies here. try not to break the whole program..
	@staticmethod
	def do_import():
		#try:
		
		global StackedConvNet, theano,T,conv,pool, save_model, load_model
		import theano
		import theano.tensor as T
		from theano.tensor.nnet import conv
		from theano.tensor.signal import pool
		from e2tomoseg_convnet import import_theano, StackedConvNet, save_model, load_model
		import_theano()
		boxerConvNet.import_done=True
		#return True
		#except: 
			#return False
	
	@staticmethod
	def do_training(args=None):
		
		refs0, refs1, badrefs=[EMData.read_images(a) for a in args]
		
			
		
		nnet_savename="nnet_pickptcls_2cls.hdf"
		bxsz=refs0[0]["nx"]
		sz=64
		shrinkfac=float(bxsz)/float(sz)
		
		print "Importing dependencies..."
		if not hasattr(boxerConvNet,'import_done'):
			if not boxerConvNet.do_import():
				print "Cannot import required dependencies..Stop."
				
		print "Setting up model ..."
		rng = np.random.RandomState(123)
		nkernel=[20,20,2]
		ksize=[15,15,15]
		poolsz=[2,1,1]
		batch_size=10
		
		image_shape=(batch_size, 1, sz, sz)
		convnet = StackedConvNet(
			rng,
			nkernel=nkernel,
			ksize=ksize,
			poolsz=poolsz,
			imageshape=image_shape
		)
		
		print "Pre-processing particles..."
		#### here we shrink the particles so they are 64x64
		#### and duplicate so there are more than 500 good and 500 bad particles
		
		nref_target=500
		
		data=[] ### particles in flattened numpy array
		lbs=[]  ### labels in flattened numpy array
		
		for label, refs in enumerate([badrefs,refs0, refs1]):
			nref=len(refs)
			if nref<5:
				print "Not enough references. Please box at least 5 good and 5 bad reference..."
				return []
			ncopy=nref_target/nref + 1
			
			for pp in refs:
				ptl=pp.process("math.fft.resample",{"n":shrinkfac})
				ptl.clip_inplace(Region(0,0, sz, sz))
				ptl.process_inplace("filter.highpass.gauss",{"cutoff_freq":0.005})
				ptl.process_inplace("filter.lowpass.gauss",{"cutoff_freq":0.05})
				for c in range(ncopy):
					
					tr=Transform()
					tr.set_rotation({"type":"2d","alpha":np.random.random()*360.0})
					img=ptl.process("xform",{"transform":tr})
					ar=img.numpy().copy()
					data.append(ar.flatten())
					lbs.append(label)
		
		
		rndid=range(len(data))
		np.random.shuffle(rndid)	
		data=[data[i] for i in rndid]
		lbs=[lbs[i] for i in rndid]
		data=np.asarray(data,dtype=theano.config.floatX)
		data/=np.std(data)
		data[data>2.]=2.
		data[data<-2.]=-2.
		lbs=np.asarray(lbs,dtype=int)
		train_set_x= theano.shared(data,borrow=True)
		
		#### make target output
		img=EMData(sz/2,sz/2)
		img.process_inplace("testimage.gaussian",{'sigma':5.})
		img.div(img["maximum"])
		gaus=img.numpy().copy().flatten()
		gaus=gaus.astype(theano.config.floatX)
		zero=np.zeros_like(gaus, dtype=theano.config.floatX)
		g0=np.hstack([gaus, zero])
		g1=np.hstack([zero, gaus])
		z0=np.hstack([zero, zero])
		lbarrs=np.array([z0, g0, g1])
		label_np=lbarrs[lbs]
		#print label_np.shape
		labels=theano.shared(label_np, borrow=True)
		
		
		print "Now Training..."
		classify=convnet.get_classify_func(train_set_x,labels,batch_size)
		learning_rate=0.002
		weightdecay=1e-5
		n_train_batches = len(data) / batch_size
		for epoch in xrange(20):
		# go through the training set
			c = []
			for batch_index in xrange(n_train_batches):
				err=classify(batch_index,
					lr=learning_rate,
					wd=weightdecay)
				c.append(err)

			learning_rate*=.96
			print 'Training epoch %d, cost ' % ( epoch),
			print np.mean(c),", learning rate",learning_rate

		
		save_model(convnet, nnet_savename)
		
		tstsz=100
		convnet.update_shape((tstsz, 1, sz, sz))
		test_cls = theano.function(
			inputs=[],
			outputs=convnet.clslayer.get_image(),
			givens={
				convnet.x: train_set_x[:tstsz]
			}
		)
		tstout=test_cls()
		trainoutfile="trainout_nnet.hdf"
		if os.path.isfile(trainoutfile):
			os.remove(trainoutfile)
		
		for i,mm in enumerate(tstout):
			t=train_set_x[i].eval().reshape(sz,sz)
			img=from_numpy(t)
			img.process_inplace("normalize")
			img.write_image(trainoutfile, -1)
			
			#print m.shape
			for m in mm:
				img=from_numpy(m)
				nx=img["nx"]
				img=img.get_clip(Region(-nx/2,-nx/2,nx*2,nx*2))
				img.scale(2.)
				#img.process_inplace("math.fft.resample",{"n":.5})
				img.mult(5)
				img.write_image(trainoutfile, -1)
			
		
	@staticmethod
	def do_autobox(micrograph,bxsz=256):

		
		nnet_savename="nnet_pickptcls_2cls.hdf"
		sz=64
		shrinkfac=float(bxsz)/float(sz)
		
		if os.path.isfile(nnet_savename)==False:
			print "Cannot find saved network, exit..."
			return 0
			
		#else:
		nx=int(micrograph["nx"]/shrinkfac)
		ny=int(micrograph["ny"]/shrinkfac)
		
			
		layers=boxerConvNet.load_network(nnet_savename, nx, ny)
		print "Applying neural net..."
		boxes=boxerConvNet.apply_network(micrograph, layers, shrinkfac, nx, ny)
		
		print "{} particles found..".format(len(boxes))
		
		
		return boxes
		
	
	@staticmethod
	def load_network(fname, nx, ny):
		print "Loading the Neural Net..."
			
		hdr=EMData(fname,0)
			
		ksize=hdr["ksize"]
		poolsz=hdr["poolsz"]
		labelshrink=np.prod(poolsz)
		k=1
		layers=[]
		for i in range(len(ksize)):
			layer={}
			b=EMData(fname,k)
			s=b["w_shape"]
			k+=1
			allw=[]
			layer["b"]=b
			layer["shp"]=s
			layer["pool"]=poolsz[i]
			for wi in range(s[0]*s[1]):
				w=EMData(fname,k)
				allw.append(w)
				k+=1
				
			allw=np.asarray(allw).reshape((s[0],s[1]))
			for wi in range(s[0]):
				
				for mi in range(s[1]):
					sw=allw[wi][mi]["nx"]
					allw[wi][mi]=allw[wi][mi].get_clip(Region((sw-nx)/2,(sw-ny)/2,nx,ny))
					
					allw[wi][mi].process_inplace("xform.phaseorigin.tocenter")
					#allw[wi][mi].do_fft_inplace()
					
			nx/=poolsz[i]
			ny/=poolsz[i]
			layer["allw"]=allw
			layers.append(layer)
		return layers
	
	@staticmethod
	def apply_network(micrograph, layers, shrinkfac, nx, ny):
		sz=64
		#### file name or EMData input
		if type(micrograph)==str:
			fm=EMData(micrograph)
		else:
			fm=micrograph.copy()
			
		#### preprocess..
		fm.process_inplace("math.fft.resample",{"n":shrinkfac})
		fm.clip_inplace(Region(0, 0, nx, ny))
		fm.process_inplace("filter.highpass.gauss",{"cutoff_freq":0.005})
		fm.process_inplace("filter.lowpass.gauss",{"cutoff_freq":0.05})
		fm.process_inplace("normalize")
		fm.process_inplace("threshold.clampminmax.nsigma", {"nsigma":2})
			
		#### apply network
		imgs=[fm]
		for layer in layers:
			
			s0=imgs[0]["nx"]
			s1=imgs[0]["ny"]
			
			imgout=[]
			allw=layer["allw"]
			s=layer["shp"]
			poolsz=layer["pool"]
			b=layer["b"]
			#print s,poolsz,s0,s1
			for wi in range(s[0]):
				
				cv=EMData(imgs[0])
				cv.to_zero()
				for mi in range(s[1]):
					ww=allw[wi][mi]
					#print ww.numpy().shape
					cv.add(imgs[mi].process("math.convolution",{"with":ww}))
				
				if poolsz>1:
					cv=cv.process("math.maxshrink",{"n":poolsz})
				cv.add(b[wi])
				cv.process_inplace("threshold.belowtozero")
				
				imgout.append(cv)
			
			imgs=imgout
			
		#### find boxes
		downsample=shrinkfac*2
		boxes=[]
		allpks=[]
		#imgs=[0,1]
		labels=["manual", "gauss"]
		for lb,img in enumerate(imgs):
			final=img.process("filter.lowpass.gauss",{"cutoff_abs":.2})
			#final.write_image("tmpout.hdf",lb)
			#final=EMData("tmpout.hdf", lb)
			
			
			try: thrn=params["threshold"]
			except:
				try: thrn=boxerConvNet.threshold.getValue()
				except:
					thrn=2.
			
			threshold=final["mean"]+final["sigma"]*thrn
			pks=final.peak_ccf(sz/4)
			pks=np.array(pks).reshape(-1,3)
			pks=np.hstack([pks, np.zeros((len(pks),1))+lb])
			#print pks.shape
			allpks.append(pks)
		
		allpks=np.vstack(allpks)
		sid=np.argsort(-allpks[:,0])
		pp=[[-sz, -sz]]
		for si in sid:
			p=allpks[si]
			if p[0]<threshold:
				break
			
			mindist=np.min(np.sum((np.array(pp)-np.array(p[1:3]))**2, axis=1))
			if mindist<sz*sz/4:
				#print np.sqrt(mindist), p
				continue
			pp.append(p[1:3])
			#boxes.append([int(p[1]*downsample),int(p[2]*downsample),"auto_convnet_p{:d}".format(int(p[3]))])
			boxes.append([int(p[1]*downsample),int(p[2]*downsample),labels[int(p[3])]])
		
		return boxes
	
	@staticmethod
	def do_autobox_all(filenames,bxsz=256, nthreads=10):
		jobs=[]
		prog=None
		#### get some parameters...
		nnet_savename="nnet_pickptcls_2cls.hdf"
		#bxsz=goodrefs[0]["nx"]
		sz=64
		shrinkfac=float(bxsz)/float(sz)
		
		## need the micrograph size to pad the kernels
		fsp=filenames[0]#.split()[1]
		hdr=EMData(fsp, 0, True)
		nx=int(hdr["nx"]/shrinkfac)
		ny=int(hdr["ny"]/shrinkfac)
		
		#### load network...
		layers=boxerConvNet.load_network(nnet_savename, nx, ny)
		#boxes=boxerConvNet.apply_network(micrograph, layers, shrinkfac, nx, ny)
		#### prepare the jobs..
		for i,fspl in enumerate(filenames):
			fsp=fspl#.split()[1]
			jobs.append((fsp, i, layers, shrinkfac, nx, ny))
		
		#### worker function
		def autobox_worker(que, job):
			fname, idx, layers, shrinkfac, nx, ny = job
			boxes=boxerConvNet.apply_network(fname, layers, shrinkfac, nx, ny)
			que.put((idx, fname,  boxes))
			return
		
		#### now start autoboxing...
		jsd=Queue.Queue(0)
		NTHREADS=max(nthreads+1,2)
		thrds=[threading.Thread(target=autobox_worker,args=(jsd,job)) for job in jobs]
		thrtolaunch=0
		while thrtolaunch<len(thrds) or threading.active_count()>1:
			#print thrtolaunch, len(thrds)
			if thrtolaunch<len(thrds) :
				
				while (threading.active_count()==NTHREADS ) : time.sleep(.1)
				print "Starting on img {}...".format(thrtolaunch)
				thrds[thrtolaunch].start()
				thrtolaunch+=1
			else: time.sleep(1)
		
			while not jsd.empty():
				idx, fsp, newboxes=jsd.get()
				print "{}) {} boxes -> {}".format(i,len(newboxes),fsp)
		
				# if we got nothing, we just leave the current results alone
				if len(newboxes)==0 : continue
			
				# read the existing box list and update
				db=js_open_dict(info_name(fsp))
				#try: 
					#boxes=db["boxes"]
					## Filter out all existing boxes for this picking mode
					#bname=newboxes[0][2]
					#boxes=[b for b in boxes if b[2]!=bname]
				#except:
					#boxes=[]
				boxes=[]	
				boxes.extend(newboxes)
				
				db["boxes"]=boxes
				db.close()
				if prog:
					prog.setValue(thrtolaunch-threading.active_count())
				
		return


if __name__ == '__main__':
	main()
	