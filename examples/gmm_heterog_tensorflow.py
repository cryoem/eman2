#!/usr/bin/env python
# Muyuan Chen 2018-01
from past.utils import old_div
from builtins import range
from builtins import object
from EMAN2 import *
import numpy as np
from EMAN2_utils import *
from sklearn.decomposition import PCA
import os
import json
os.environ["CUDA_VISIBLE_DEVICES"]="1"
import tensorflow as tf
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 

def main():
	
	usage="make 3d volume from 2d using gaussian balls"
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ptcls", type=str,help="Particle file.", default=None)
	parser.add_argument("--map3d", type=str,help="3D density map file. ", default=None)
	
	parser.add_argument("--initpos", type=str,help="file of initial ball position, from e2segment3d.py", default=None)
	parser.add_argument("--loadfrom", type=str, default=None,help="load json file from previous run")
	parser.add_argument("--moviepath", type=str, default="motion_00",help="path for the movie")
	
	parser.add_argument("--learnrate", type=float,help="learning rate", default=1e-1)
	parser.add_argument("--maxres", type=float,help="maximum resolution, default 10.", default=10.)
	parser.add_argument("--minres", type=float,help="minimum resolution, default 300.", default=300.)
	parser.add_argument("--batchsize", type=int,help="batch size. should be as large as possible, as long as they fit in the GPU.", default=500)
	
	parser.add_argument("--savefile", type=str,help="file to save the model", default="gmmsave.json")
	parser.add_argument("--sym", type=str,help="symmetry", default="c1")
	parser.add_argument("--path", type=str,help="path for output", default="gmm_00")
	parser.add_argument("--calcgrad", action="store_true", default=False ,help="calculate gradient for ball positions")
	parser.add_argument("--calcmotion", action="store_true", default=False ,help="calculate motion vector from gradients")
	parser.add_argument("--regress", action="store_true", default=False ,help="regress particles on the motion vector")
	parser.add_argument("--make3d", action="store_true", default=False ,help="make 3D from particles of different conformation")
	parser.add_argument("--nframe", type=int,help="number of frames for the 3D movie", default=3)
	parser.add_argument("--gausswidth", type=float,help="width of the gaussians", default=15.)
	parser.add_argument("--threads", type=int,help="number of threads for make3d", default=12)
	parser.add_argument("--listid", type=str, default=None ,help="list of gaussian indices to use while calculating motion. File containing int list separated by linebreak, or chimera list output.")

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	#### make sure input exist and load box size
	if options.map3d:
		e=EMData(options.map3d, 0, True)
	elif options.ptcls:
		e=EMData(options.ptcls, 0, True)
	else:
		print("No data input, exit...")
		return
		
	options.sz=sz=e["nx"]
	options.apix=apix=e["apix_x"]
	
	#### setup the model
	### load from existing json file
	if options.loadfrom: 	
		js=js_open_dict(options.loadfrom)
		gmm=GaussianModel(js.data, freq_bound=[options.minres, options.maxres])
	
	### load Gaussian position from text file and set width and amplitude by default
	elif options.initpos!=None:	
		cents=np.loadtxt(options.initpos)
		nball=len(cents)
		cents=cents[:,1:]-old_div(sz,2)
		## set width so model fills box
		initw=1.5*(old_div(np.cbrt(old_div(((sz*.6)**3),nball)),(4./3.*np.pi)))**2 

		wts=np.zeros(len(cents))+initw
		amp=np.ones(len(cents))
		
		model={
			"center":cents,
			"width":wts, 
			"boxsz":sz, 
			"sym":options.sym,
			"amp":amp,
			"apix":apix,
			"motvec":np.zeros_like(cents)
			}
		gmm=GaussianModel(model, freq_bound=[options.minres, options.maxres])
	
	else:
		print("No model input, exit...")
		return
	
	#### make working directory
	try: os.mkdir(options.path)
	except: pass
	
	jspath=os.path.join(options.path, "gmm_save.json")
	
	#### optimize model using an average map
	if options.map3d:
		fit_3d_map(gmm, options)
		gmm.save_model(jspath)
		
		print("Making density map for the model...")
		mdpath=os.path.join(options.path, "gmm_model_avg.hdf")
		e=gmm.get_map()
		e.write_image(mdpath)
		
		print("Fitting complete. Model saved to {}, density map saved to {}.".format(jspath, mdpath))
	#return
	#### real data input...
	if options.ptcls:
		#### load particles
		data,orient=load_particles(options.ptcls)
		
		#### calculate motion gradient
		print("Calculating per particle motion gradient...")
		grds=gmm.calc_grad_ptcls(data, orient, gmm.balls_pos, options.batchsize)
		grds=grds.reshape((-1, gmm.nball, 3))
		
		#### compute eigen-vector
		print("Calculation eigen-motion...")
		if options.listid==None:
			selid=[]
		else:
			#### parse chimera output file for selected atoms
			f=open(glst,'r')
			selid= np.array( [int(l[3:].strip("ALC .\n")) for l in f])
			f.close()
			print("Gaussians used for motion calculation: ", selid)
		
		vec=gmm.compute_eigen_vector(grds, selid)
		
		##### set the motion vector
		gmm.set_mot_vec(vec)
		gmm.save_model(jspath)
		mvpath=options.moviepath
		print("Making density map for the model...")
		try: os.mkdir(os.path.join(options.path,mvpath))
		except: pass
		mdpath=os.path.join(options.path,mvpath, "motion_model.hdf")
		try: os.remove(mdpath)
		except: pass
	
		for c in [-1.,-.5,0.,.5,1.]:
			gmm.set_conf(c)
			e=gmm.get_map()
			e.write_image(mdpath, -1)
		
		gmm.set_conf(0.)
		
		#### regress particles (only one iteration)
		print("Calculating particle conformations...")
		conf=gmm.calc_grad_ptcls(data, orient, gmm.conf_ballpos, options.batchsize)
		conf=conf.T[0]
		
		np.savetxt(os.path.join(options.path,"conf_ptcls.txt"),np.vstack([np.arange(len(conf)), conf]).T)
		
		make3d(options.ptcls, options, conf)
	
	E2end(logid)
	return

def make3d(ptcl_file, options, allconf):
	print("Making 3D movie...")
	mvpath=options.moviepath
	try: os.mkdir(os.path.join(options.path,mvpath))
	except: pass
	lstin=LSXFile(ptcl_file, True)
	mvlen=np.std(allconf)
	stepsz=old_div(mvlen,((options.nframe-1)/2.))
	framepos=np.arange(-mvlen,mvlen+1e-14, stepsz)+np.mean(allconf)
	print("Motion steps : Number of particles")
	winsz=stepsz*.6
	lstnames=[]
	lstlen=[]
	for kk,t in enumerate(framepos):
		lstlen.append(np.sum(np.abs(allconf-t)<winsz))
		print(t,"  :  ", lstlen[-1])
		idx= np.where(np.abs(allconf-t)<winsz)[0]
		fname=os.path.join(options.path, mvpath, "lst_{:02d}.lst".format(kk))
		try: os.remove(fname)
		except: pass
		lstnames.append(fname)
		lst=LSXFile(fname,False)
		for i in idx:
			ln=lstin.read(i)
			lst.write(-1,ln[0],ln[1],ln[2])
		lst=None
		
	lstin=None
	
	pad=good_size(options.sz*1.7)
	mapnames=[]
	for i,lname in enumerate(lstnames):
		outname=lname.replace("lst_", "map_")
		outname=outname[:-4]+".hdf"
		mapnames.append(outname)
		print("Making 3D for {}...".format(lname))
		cmd="make3dpar_rawptcls.py --input {} --output {} --sym {} --apix {} --pad {} --mode gauss_5 --threads {}".format(lname, outname, options.sym, options.apix, pad, options.threads)
		
		launch_childprocess(cmd)
		
	fname=os.path.join(options.path, mvpath, "maps_filt.hdf")
	for i,name in enumerate(mapnames):
		e=EMData(name)
		e.process_inplace("filter.lowpass.gauss",{"cutoff_pixels":15})
		e.process_inplace("normalize")
		e["conf_val"]=framepos[i]
		e.write_image(fname,i)
		
	print("Motion movie is saved to ", fname)



##### fit to average model
def fit_3d_map(gmm, options):
	
	print("Loading 3D map and making projections....")
	
	sym=parsesym(options.sym)
	xfs=sym.gen_orientations('eman', {"delta":12})
	# xfs=sym.gen_orientations('rand', {"n":300, "inc_mirror":1, "phitoo":1})
	imgs=[]
	imgs_r=[]
	imgs_i=[]
	oris=[]
	num=len(xfs)
	e=EMData(options.map3d)
	sz=e["nx"]
	print("Map size {}, symmetry {}, making {} projections..".format(sz, options.sym, len(xfs)))

	apix=e["apix_x"]
	for xf in xfs:
		pj=e.project('standard', xf)
		img=pj.numpy().copy()
		
		imgs.append(img)
		x=xf.get_params('eman')
		oris.append([x["az"],x["alt"],x["phi"], 1, x["tx"], x["ty"]])

	data=np.asarray(imgs,dtype=np.float32)
	data/=np.std(data)
	data=np.fft.fftshift(data, axes=[1,2])

	orient=np.asarray(oris,dtype=np.float32)
	orient[:,:3]=orient[:,:3]/180.0*3.14
	
	gmm.optimize(data, orient, options.learnrate)
		
### load particles (with orientations)
def load_particles(ptcl_file):
	
	
	print("loading particles...")
	ref=[]
	ori=[]
	imgs_r=[]
	imgs_i=[]
	num=EMUtil.get_image_count(ptcl_file)
	lst=LSXFile(ptcl_file, True)
	data_xfs=[]
	for i in range(num):
		e=EMData(ptcl_file,i)
		if e["sigma"]==0:
			print("Warning: empty particle (#{})".format(i))
			continue

		lstinfo=lst.read(i)
		pj=Transform(eval(lstinfo[2])).get_params("eman")
		data_xfs.append(pj)

		if pj["mirror"]>0: mr=-1
		else: mr=1
			
		ori.append([pj["az"],pj["alt"],pj["phi"],mr, pj["tx"], pj["ty"]])
		img=e.numpy().copy()
		
		
		ref.append(img)
		
	data=np.asarray(ref,dtype=np.float32)

	data/=np.std(data)
	data=np.fft.fftshift(data, axes=[1,2])
	orient= np.asarray(ori,dtype=np.float32)
	orient[:,:3]=orient[:,:3]/180.0*3.14
	
	return data, orient

	

	
class GaussianModel(object):
	
	def __init__(self, model, freq_bound=[300,10]):
		
		#############
		#### setup the model in tf
		#### model should contain center position, width and amplitude of the Gaussians
		#### freq_bound is the min and max resolution considered
		print("Setting up model...")
		
		#### load parameters from model
		self.sz=sz=model["boxsz"]
		self.apix=apix=model["apix"]
		self.sym=sym=str(model["sym"]).lower()
		self.balls_pos=balls_pos=tf.Variable(np.asarray(model["center"],dtype=np.float32))
		self.balls_wt=balls_wt=tf.Variable(np.asarray(model["width"],dtype=np.float32))
		self.balls_amp=balls_amp=tf.Variable(np.asarray(model["amp"],dtype=np.float32))
		self.vec_ballpos=vec_ballpos=tf.Variable(np.asarray(model["motvec"], dtype=np.float32))
		
		self.nball=nball=len(model["center"])
		## conformation value
		self.conf_ballpos=conf_ballpos=tf.Variable(0.)
		
		#### initialize tf session
		self.session=session=tf.Session()
		session.run(tf.global_variables_initializer())
		
		#### particle & orientation input
		self.data_input=data_input=tf.placeholder(tf.complex64, shape=[None, sz, sz], name='data_input')
		self.orient_input=orient_input=tf.placeholder(tf.float32, shape=[None, 6], name='orient_input')
		data_cpx=tf.fft2d(data_input) #### FFT the particles here

		#### now generate rotation matrix from eular angles
		ang=orient_input
		azp=ang[:,0]+np.pi
		altp=np.pi-ang[:,1]
		phip=np.pi-ang[:,2]-old_div(np.pi,2)

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
		
		#### transform Gaussian positions
		self.ballpos_conf=ballpos_conf=balls_pos+self.conf_ballpos*self.vec_ballpos
		ballpos_rot=tf.tensordot(ballpos_conf,matrix, [[1],[2]])
		ballpos_rot=tf.transpose(ballpos_rot, [1,0,2])

		tx=ang[:,5][:,None]
		ty=ang[:,4][:,None]
		mirror=ang[:,3][:,None]

		ballpos_rot_trans=tf.stack([ballpos_rot[:,:,0]+tx+old_div(sz,2), mirror*(ballpos_rot[:,:,1]+ty)+old_div(sz,2), ballpos_rot[:,:,2]], 0)
		ballpos_rot_trans=tf.transpose(ballpos_rot_trans, [1,2,0])
		bpos=ballpos_rot_trans[:,:, :2]
		
		#### make 2D projection images from the model
		grid_x, grid_y=np.indices((sz,sz))
		bpx=bpos[:,:,0][:,:, None, None]
		bpy=bpos[:,:,1][:,:, None, None]
		px=grid_x-bpx
		py=grid_y-bpy
		wt=balls_wt[:, None, None]
		amp=balls_amp[:, None, None]
		pgauss=0.5*np.sqrt(np.pi)*amp*tf.sqrt(wt)*tf.exp(old_div(-(px**2+py**2),wt))
		imgs=tf.reduce_sum(pgauss, axis=1)
		
		#### make 2D projections in Fourier space
		bposft=(bpos-sz/2.)*np.pi
		bpxft=bposft[:,:,0][:,:, None, None]
		bpyft=bposft[:,:,1][:,:, None, None]
		gridxft=(old_div((grid_x-sz/2.),sz)*2.).astype(np.float32)
		gridyft=(old_div((grid_x-sz/2.),sz)*2.).astype(np.float32)

		gridxft=np.fft.ifftshift(gridxft)
		gridyft=np.fft.ifftshift(gridyft)

		r2=gridxft**2+gridyft**2
		pgauss_amp=0.5*np.sqrt(np.pi)*amp*tf.sqrt(wt)*wt*tf.exp(-r2*wt/.4)
		pgauss_real=pgauss_amp*tf.cos(gridxft*bpxft+gridyft*bpyft)
		pgauss_imag=-pgauss_amp*tf.sin(gridxft*bpxft+gridyft*bpyft)

		imgs_real=tf.reduce_sum(pgauss_real, axis=1)
		imgs_imag=tf.reduce_sum(pgauss_imag, axis=1)
		imgs_cpx=tf.complex(imgs_real, imgs_imag)
		
		img_out=tf.real(tf.ifft2d(imgs_cpx))
		data_out=tf.real(tf.ifft2d(data_cpx))
		
		#### compute particle-projection FRC 
		x,y= grid_x, grid_y
		rr=np.sqrt((x - old_div(sz,2))**2 + (y - old_div(sz,2))**2).astype(int)
		rr=np.fft.fftshift(rr)
		
		rings=np.zeros((sz,sz,old_div(sz,2)), dtype=np.float32) #### Fourier rings
		for i in range(old_div(sz,2)):
			rings[:,:,i]=(rr==i)
		
		#### normalization per ring
		nrm_img=tf.real(imgs_cpx)*tf.real(imgs_cpx)+tf.imag(imgs_cpx)*tf.imag(imgs_cpx)
		nrm_data=tf.real(data_cpx)*tf.real(data_cpx)+tf.imag(data_cpx)*tf.imag(data_cpx)
		#### make sure normalization factors >= 0 
		## this really should not be necessary but x**2<0 does happen sometimes..
		## maybe some round up errors?
		nrm0=tf.nn.relu(tf.tensordot(nrm_img, rings, [[1,2],[0,1]])) 
		nrm1=tf.nn.relu(tf.tensordot(nrm_data, rings, [[1,2],[0,1]]))
		nrm=tf.sqrt(nrm0)*tf.sqrt(nrm1)
		nrm=tf.maximum(nrm, 1e-5) #### so we do not divide by 0
		#### min/max resolution considered
		freq=1./np.fft.fftfreq(sz, apix)[:old_div(sz,2)]
		fq= np.where(np.logical_and(freq<freq_bound[0], freq>freq_bound[1]))[0]
		#### average FRC per batch
		ccc=tf.real(imgs_cpx)*tf.real(data_cpx)+tf.imag(imgs_cpx)*tf.imag(data_cpx)
		frc=old_div(-tf.tensordot(ccc, rings, [[1,2],[0,1]]),nrm)
		self.score=score=tf.reduce_mean(frc[:,fq[0]:fq[-1]], axis=1)
		self.loss=loss=tf.reduce_mean(score)
		
	#### optimize model to fit data
	def optimize(self, data, orient, learnrate=1e-1):
		session=self.session
		loss=self.loss
		data_input=self.data_input
		orient_input=self.orient_input
		
		# train_step = tf.train.GradientDescentOptimizer(learnrate).minimize(loss, var_list=[balls_wt])
		optimizer=tf.train.AdamOptimizer(learnrate)
		train_wt = optimizer.minimize(loss, var_list=[self.balls_wt])
		train_amp = optimizer.minimize(loss, var_list=[self.balls_amp])
		train_pos = optimizer.minimize(loss, var_list=[self.balls_pos])

		session.run(tf.global_variables_initializer())
		
		print("iter\tscores:")
		for i in range(300):
    
			_, scr0=session.run((train_wt, loss), feed_dict={data_input:data, orient_input:orient})
			_, scr1=session.run((train_amp, loss), feed_dict={data_input:data, orient_input:orient})
			_, scr2=session.run((train_pos, loss), feed_dict={data_input:data, orient_input:orient})
			if i%100==0 or scr2<-.95:
				print(i, scr0, scr1, scr2)
			if scr2<-.95:
				break
				
		
	#### generate density map from model
	def get_map(self):
		
		session=self.session
		sz=self.sz
		bp=session.run(self.ballpos_conf)+old_div(sz,2)
		wt=session.run(self.balls_wt)
		amp=session.run(self.balls_amp)
		#print bp.shape, wt.shape
		mp=np.zeros((sz,sz,sz))
		ind_np=np.indices((sz,sz,sz)).astype(float)
		ind_np=np.transpose(ind_np,axes=(3,2,1,0))
		for i in range(len(bp)):
			p=bp[i]
			w=wt[i]
			a=amp[i]
			d=(ind_np-p)**2
			w=wt[i]
			mp+=a*np.exp(old_div(-np.sum(d,axis=3),(w)))
		e=from_numpy(mp)
		#     e.process_inplace("xform.applysym",{"sym":"d7"})
		return e

	#### save model to json file
	def save_model(self, fname="gmm_save.json"):
		
		
		session=self.session
		model={
			"center":session.run(self.balls_pos).tolist(),
			"width":session.run(self.balls_wt).tolist(), 
			"boxsz":self.sz, 
			"sym":self.sym,
			"amp":session.run(self.balls_amp).tolist(),
			"apix":self.apix,
			"motvec":session.run(self.vec_ballpos).tolist()
			}
		
		f=open(fname, 'w')
		json.dump(model, f)
		f.close()
	
	#### calculate motion gradient with respect of ''to'', from particles data
	def calc_grad_ptcls(self, data, orient, to, batchsize=500):
		
		n=batchsize
		#### batch size. since we compute jacobian of the batch, the value does not affect the result
		#### since we loop over each batch, having a larger batch size makes the process faster
		#### but make sure the whole batch fits in your GPU memory....
		
		data_input=self.data_input
		orient_input=self.orient_input
		score=self.score
		session=self.session
		
		grads=[]
		jac=calc_jacobian(score, to, n)
		for bi in range(len(data)//n):
			grads.append(session.run(jac, feed_dict={data_input:data[bi*n:(bi+1)*n], orient_input:orient[bi*n:(bi+1)*n]}))
			if bi%10==0:
				print("{}/{} batches finished".format(bi+1, len(data)//n+1))
		
		#### deal with the last batch
		rs=len(data)-len(grads)*n
		if rs>0:
			jac=calc_jacobian(score, to, rs)
			grads.append(session.run(jac, feed_dict={data_input:data[-rs:], orient_input:orient[-rs:]}))

		all_grd=np.vstack(grads)#.reshape((-1, nball, 3))
		
		return all_grd

	#### compute eign motion vector. 
	## selid: only consider a subset of gaussians
	def compute_eigen_vector(self, grds, selid=[]):
		
		session=self.session
		nball=self.nball
		
		## compensate for width and amp...
		aa=session.run(self.balls_amp)
		wt=session.run(self.balls_wt)
		wt/=np.max(wt)
		mlt=aa*wt*np.sqrt(wt)*np.exp(-wt)
		mlt/=np.max(mlt)

		if len(selid)==0:
			selid=np.arange(nball)
		
		ncopy=1 ### for symmetry
		idx=np.zeros((old_div(nball,ncopy), 3), dtype=bool)
		idx[selid]=1
		idx=idx.flatten()
		
		#### flatten and normalize the gradients
		grd_flatten=old_div(grds.reshape((-1, nball, 3)),mlt[None,:, None])
		grd_flatten=grd_flatten.reshape((-1, nball*3))
		
		#### now do pca
		#### only take the first eigen-vector for now...
		pca=PCA(1)
		ptmot=pca.fit_transform(grd_flatten[:,idx])
		pc=np.zeros(old_div(nball*3,ncopy))
		pc[idx]=pca.components_[0]
		pc=pc.reshape((old_div(nball,ncopy), -1))
		
		#### set the maximum amplitude of eigen-vecot to be the maximum width
		rr=np.sqrt(np.max(session.run(self.balls_wt)))
		pc=pc/np.sqrt(np.max(np.sum(pc**2, 1)))*rr/8.

		pc=pc.astype(np.float32)
		return pc
	
	#### set motion vector
	def set_mot_vec(self, vec):
		self.session.run(self.vec_ballpos.assign(np.asarray(vec, dtype=np.float32)))
	
	#### set conformation
	def set_conf(self, conf):
		self.session.run(self.conf_ballpos.assign(conf))
		#print self.session.run(self.conf_ballpos)
		#print self.session.run(self.ballpos_conf)
		#self.conf_ballpos=tf.Variable(conf)
		
		
	
#### jacobian calculation
#### copied from jeisses in https://github.com/tensorflow/tensorflow/issues/675
def calc_jacobian(y, x, n):

	loop_vars = [
		tf.constant(0, tf.int32),
		tf.TensorArray(tf.float32, size=n),
	]

	_, jacobian = tf.while_loop(
		lambda j, _: j < n,
		lambda j, result: (j+1, result.write(j, tf.gradients(y[j], x))),
		loop_vars)

	return jacobian.stack()
	
	
	
if __name__ == '__main__':
	main()
	
