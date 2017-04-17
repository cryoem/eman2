#!/usr/bin/env python
# Muyuan Chen 2015-08
from EMAN2 import *
import numpy as np
import random
import theano
import theano.tensor as T
from theano.tensor.shared_randomstreams import RandomStreams
from sklearn.decomposition import PCA


def main():
	
	usage="make 3d volume from 2d using gaussian balls"
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--learnrate", type=float,help="learning rate", default=None)
	parser.add_argument("--savefile", type=str,help="file to save the model", default="gmmsave.json")
	parser.add_argument("--sym", type=str,help="symmetry", default="c1")
	parser.add_argument("--map3d", type=str,help="3d density map file", default=None)
	parser.add_argument("--path", type=str,help="path for output", default="gmm_00")
	parser.add_argument("--initpos", type=str,help="file of initial ball position, from e2segment3d.py", default=None)
	parser.add_argument("--loadfrom", type=str, default=None,help="load json file from last run")
	parser.add_argument("--calcgrad", action="store_true", default=False ,help="calculate gradient for ball positions")
	parser.add_argument("--calcmotion", action="store_true", default=False ,help="calculate motion vector from gradients")
	parser.add_argument("--regress", action="store_true", default=False ,help="regress particles on the motion vector")
	parser.add_argument("--make3d", action="store_true", default=False ,help="make 3D from particles of different conformation")
	parser.add_argument("--nframe", type=int,help="number of frames for the 3D movie", default=6)
	parser.add_argument("--sz", type=int,help="size of box", default=-1)
	parser.add_argument("--gausswidth", type=float,help="width of the gaussians", default=15.)
	parser.add_argument("--threads", type=int,help="number of threads for make3d", default=12)
	parser.add_argument("--listid", type=str, default=None ,help="list of gaussian indices to use while calculating motion. File containing int list separated by linebreak, or chimera list output.")

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	#### Read basic info from particle header
	ptcl_file=args[0]
	e=EMData(ptcl_file,0, True)
	options.apix=e["apix_x"]
	if options.sz<0:
		options.sz=e["nx"]
		realsz=options.sz
	else:
		realsz=e["nx"]
	e=None
	
	#### make working directory
	try: os.mkdir(options.path)
	except: pass
	
	
	####
	print "Setting up model..."
	
	if options.loadfrom:
		js=js_open_dict(options.loadfrom)
		ballrec=BallsReconstruction(js.data)
		
	elif options.initpos!=None:
		cents=np.loadtxt(options.initpos)
		cents=cents[:,1:]-realsz/2
		wts=np.zeros(len(cents))
		model={"center":cents, "weight":wts, "boxsz":options.sz, "sym":options.sym,"width":options.gausswidth}
		ballrec=BallsReconstruction(model)
	
	else:
		print "No model input, exit..."
		return
	
	if options.map3d:
		if not options.learnrate:
			options.learnrate=.01
		train_3d(ballrec, options)
		write_3d(ballrec, options)
		save_model(ballrec,options)
	
	
	if options.calcgrad:
		data,orient=load_particles(ptcl_file)
		#data,orient=0,0
		grad=calc_grad(ballrec, options, data, orient)
	
	
	if options.calcmotion:
		if not options.calcgrad:
			gdsave=os.path.join(options.path,"grad_balls.txt")
			try:
				grad=np.loadtxt(gdsave)
			except:
				print "Cannot find gradients from {}".format(gdsave)
				return
		
		vec=calc_motion(ballrec, options, grad, options.listid)
	
	
	if options.regress:
		if not options.learnrate:
			options.learnrate=1e3
		if not options.calcgrad:
			data,orient=load_particles(ptcl_file)
		if not options.calcmotion:
			vsave=os.path.join(options.path,"motion_vec.txt")
			try:
				vec=np.loadtxt(vsave)
			except:
				print "Cannot find gradients from {}".format(vsave)
				return
			print "Load motion from ", vsave
		
		allconf=regress_particles(ballrec, options, data, orient, vec)
		
		
	if options.make3d:
		if not options.regress:
			cfsave=os.path.join(options.path,"conf_out.txt")
			try:
				allconf=np.loadtxt(cfsave)[:,2]
			except:
				print "Cannot load conformation from ", cfsave
			print "Loaded conformation for {} particles".format(len(allconf))
		
		make3d(ptcl_file, options, allconf)
	
	E2end(logid)
	return

def make3d(ptcl_file, options, allconf):
	print "Making 3D movie..."
	mvpath="motion_01"
	try: os.mkdir(os.path.join(options.path,mvpath))
	except: pass
	lstin=LSXFile(ptcl_file, True)
	mvlen=np.std(allconf)
	stepsz=mvlen/((options.nframe-1)/2.)
	framepos=np.arange(-mvlen,mvlen+.1, stepsz)+np.mean(allconf)
	print "Motion steps : Number of particles"
	winsz=stepsz*.6
	lstnames=[]
	lstlen=[]
	for kk,t in enumerate(framepos):
		lstlen.append(np.sum(np.abs(allconf-t)<winsz))
		print t,"  :  ", lstlen[-1]
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
		print "Making 3D for {}...".format(lname)
		cmd="make3dpar_rawptcls.py --input {} --output {} --sym {} --apix {} --pad {} --mode gauss_5 --threads {}".format(lname, outname, options.sym, options.apix, pad, options.threads)
		
		launch_childprocess(cmd)
		
	fname=os.path.join(options.path, mvpath, "maps_filt.hdf")
	for i,name in enumerate(mapnames):
		e=EMData(name)
		e.process_inplace("filter.lowpass.gauss",{"cutoff_pixels":15})
		e.process_inplace("normalize")
		e["conf_val"]=framepos[i]
		e.write_image(fname,i)
		
	print "Motion movie is saved to ", fname

def regress_particles(ballrec, options, data, orient, vec):
	
	print "Regress the particles..."
	num=data.shape[0].eval()
	all_conf=np.zeros(num, dtype=theano.config.floatX)
	ballrec.movvec.set_value(vec.astype(theano.config.floatX))
	index = T.lscalar()
	train_conf=theano.function(
		inputs=[ballrec.learnrate,index],
		outputs=ballrec.cost,
		updates=ballrec.update_conf,
		givens={
			ballrec.imgin: data[index],
			ballrec.orientin: orient[index],
			}
		)

	
	lstout=os.path.join(options.path,"conf_out.txt")
	try:os.remove(lstout)
	except: pass
	for i in range(num):
		cfval=all_conf[i]
		ballrec.conf.set_value(cfval)
		cc=[cfval]
		lr=options.learnrate
		for ti in range(2):
			loss=train_conf(lr,i)
			cc.append(ballrec.conf.get_value())
			lr*=.95
	#		 print i, conf.get_value(), cc[-1]-cc[-2]
			if abs(cc[-2]-cc[-1])<.001: break
	#	 print i, "  iter: ",ti, "conf: ", cc[-1]
		#print loss
		print i, ti, cc[-1],loss
		ff=open(lstout,'a')
		ff.write("{}\t{}\t{}\n".format(i, ti, cc[-1]))
		ff.close()
		all_conf[i]=cc[-1]
	
	print "Conformation of particles written to ", lstout
	return all_conf

def calc_motion(ballrec, options, grad, glst=None):

	npballs=ballrec.ballzero.get_value()
	if glst==None:
		glst=np.arange(len(npballs))
	else:
		f=open(glst,'r')
		glst= np.array( [int(l[3:].strip("ALC .\n")) for l in f])
		f.close()
		print "Gaussians used for motion calculation: ", glst
		
	gradsmall=grad.reshape((len(grad),len(npballs),3))[:,glst,:].reshape((len(grad),len(glst)*3))
	pca=PCA(3)
	ptmot=pca.fit_transform(gradsmall)
	mot=pca.components_[0].reshape((len(glst), 3))
	vec=np.zeros_like(npballs)
	vec[glst]+=mot
	mtsave=os.path.join(options.path,"motion_vec.txt")
	np.savetxt(mtsave,vec)
	print "Motion vector is saved to {}".format(mtsave)
	get_map=theano.function([], ballrec.map3d)
	ballrec.movvec.set_value(vec)
	mvlen=3.
	pp=np.arange(-mvlen,mvlen+.1, mvlen/((options.nframe-1)/2.))
	pp=pp.astype(theano.config.floatX)
	print "Motion steps are ",pp
	mpsave=os.path.join(options.path,"motion_model1.hdf")
	try: os.remove(mpsave)
	except: pass

	for i,p in enumerate(pp):
		ballrec.conf.set_value(p)
		mp=get_map()
		e=from_numpy(mp)
		e.write_image(mpsave, i)
	print "Motion maps are saved to {}".format(mpsave)
	return vec

def load_particles(ptcl_file):
	
	### load particles (with orientations)
	print "loading particles..."
	ref=[]
	ori=[]
	num=EMUtil.get_image_count(ptcl_file)
	#num=5000
	lst=LSXFile(ptcl_file, True)
	for i in range(num):
		e=EMData(ptcl_file,i)
		if e["sigma"]==0:
			print "Warning: empty particle (#{})".format(i)
			continue
		
		lstinfo=lst.read(i)
		pj=Transform(eval(lstinfo[2])).get_params("eman")
		
		if pj["mirror"]>0: mr=-1
		else: mr=1
		ori.append([pj["az"],pj["alt"],pj["phi"],mr, pj["tx"], pj["ty"]])
		a=e.numpy().copy()
		ref.append(a.flatten())
	data=np.asarray(ref,dtype=theano.config.floatX)
	orient= np.asarray(ori,dtype=theano.config.floatX)
	orient[:,:3]=orient[:,:3]/180.0*3.14
	#print np.std(data.get_value(),axis=1)
	data_t=theano.shared(value=data, name='imgs', borrow=True)
	orient_t=theano.shared(value=orient, name='orients', borrow=True)
	print "Finished reading {} particles with orientation".format(len(data))
	return data_t, orient_t

def numpy2pdb(data, fname):
	if data.shape[1]==3:
		data=np.hstack([data,np.zeros((len(data),1))])
	if data.shape[1]!=4:
		print "wrong dimension"
		return
	print data.shape
	data=np.asarray(data.astype(float).flatten())
	if len(data[0].shape)>0: data=data[0]
	lst=data.tolist()
#	 print lst
#	 print type(lst[0])
	p=PointArray()
	p.set_from(lst)
	p.save_to_pdb(fname)


def save_model(ballrec,options):
	if not (options.path in options.savefile):
		options.savefile=os.path.join(options.path,options.savefile)
	
	print "saving model to {}...".format(options.savefile)
	cents=ballrec.ballzero.get_value()
	wts=ballrec.weight.get_value()
	js=js_open_dict(options.savefile)
	js["center"]=cents.tolist()
	js["weight"]=wts.tolist()
	js["boxsz"]=options.sz
	js["sym"]=options.sym
	js["width"]=options.gausswidth
	js=None
	


#### calculate gradient of balls from ball positions
def calc_grad(ballrec, options, data, orient):
	npballs=ballrec.ballzero.get_value()
	npballs*=options.apix
	nballs=len(npballs)
	index = T.lscalar()
	num=data.shape[0].eval()
	print "Calculating gradient of {} particles...".format(num)
	get_grad_ballpos=theano.function(inputs=[index],
					outputs=[ballrec.cost,ballrec.grad_ballpos],
					givens={ 
						ballrec.imgin: data[index],
						ballrec.orientin: orient[index],
						}
					)
	
	get_2dimg_out=theano.function(inputs=[index],
					outputs=[ballrec.out],
					givens={ 
						#ballrec.imgin: data[index],
						ballrec.orientin: orient[index],
						}
					)
	gd_all=[]
	for b in range(num):
		if b<50:
			img=get_2dimg_out(b)[0]
			m0= data[b].eval().reshape((options.sz,options.sz))
			e=from_numpy(m0)
			e.write_image(options.path+"/tmpimgcmps.hdf",b*2)
			
			e=from_numpy(img.copy())
			e.write_image(options.path+"/tmpimgcmps.hdf",b*2+1)
			#exit()
		c,gd=get_grad_ballpos(b)
		gd=np.array(gd)
		print b,c,gd[0]
		gd_all.append(gd)
	print np.array(gd_all).shape
	gd_all=np.array(gd_all).reshape((int(num), nballs*3))	
	gdsave=os.path.join(options.path,"grad_balls.txt")
	np.savetxt(gdsave,gd_all)
	print "Gradient is saved to {}".format(gdsave)
	#gd_all=np.loadtxt(gdsave)
	
	#### variance of gradient
	nmm=np.mean(gd_all,axis=0)**2
	var=np.asmatrix(np.sum(nmm.reshape(npballs.shape),axis=1)).T
	var=var/np.max(var)*100
	tosave=np.hstack([npballs, var])
	varsave=os.path.join(options.path,"grad_amp.pdb")
	print "Gradient amplitude per Gaussian is saved to {}.".format(varsave)
	numpy2pdb(tosave, varsave)
	nmm=np.std(gd_all,axis=0)
	var=np.asmatrix(np.sum(nmm.reshape(npballs.shape),axis=1)).T
	var=var/np.max(var)*100
	tosave=np.hstack([npballs, var])
	varsave=os.path.join(options.path,"grad_std.pdb")
	numpy2pdb(tosave, varsave)
	print "Gradient variance per Gaussian is saved to {}.".format(varsave)
	
	return gd_all


#### train model using 3d map input
def train_3d(ballrec, options):
	print "Adjusting the model using 3D density map..."
		
	e=EMData(options.map3d)
	e.process_inplace("normalize")
	map3d=e.numpy().copy()
	map3d[map3d<0]=0
	map3d/=np.std(map3d)*10.
	train_map_p=theano.function(	inputs=[ballrec.learnrate],
					outputs=ballrec.map_err,
					updates=ballrec.map_update_p,
					givens={
						ballrec.target_map: map3d
						}
					)
	
	train_map_w=theano.function(	inputs=[ballrec.learnrate],
					outputs=ballrec.map_err,
					updates=ballrec.map_update_w,
					givens={
						ballrec.target_map: map3d
						}
					)
	lr=options.learnrate
	for j in range(5):
		print "iter {}, train weights...".format(j)
		for i in range(5):
			print "err: ",train_map_w(lr)
		print "		  train position..."
		for i in range(5):
			print "err: ",train_map_p(lr*.5)
		lr*=.5
	
	
#### write 3d output
def write_3d(ballrec,options):
	
	print "Writing density map output..."
	threedout=os.path.join(options.path,"gmm_3d.hdf")
	try: os.remove(threedout)
	except: pass
	get_map=theano.function([], ballrec.map3d)
		
	mapnp=get_map()
	print "Map shape:  ",np.shape(mapnp), ", apix: ", options.apix
	outmap=from_numpy(mapnp)
	outmap["apix_x"]=options.apix
	outmap["apix_y"]=options.apix
	outmap["apix_z"]=options.apix
	print "Map saved to ", threedout
	outmap.write_image(threedout)
	

	
class BallsReconstruction(object):
	
	def __init__(self, model):
		
		#### theano params 
		learnrate=T.scalar('lr')
		imgin=T.vector('imgin')
		ang=T.vector('orient')
		
		
		#### load model
		sz=self.boxsz=model["boxsz"]
		sym=str(model["sym"]).lower()
		balls_pos=np.asarray(model["center"],dtype=theano.config.floatX)
		balls_wt=np.asarray(model["weight"],dtype=theano.config.floatX)
		
		
		ballzero = theano.shared(value=balls_pos, name='balls', borrow=True)
		weight = theano.shared(value=balls_wt, name='wts', borrow=True)
		wts=T.nnet.sigmoid(weight)
		nballs=len(balls_pos)
		
		### motion vector of the balls
		mov_vec= np.zeros((nballs,3),dtype=theano.config.floatX)
		movvec=theano.shared(value=mov_vec, name="move_vec", borrow=True)
		
		
		cfval=np.array(0,dtype=theano.config.floatX)
		conf=theano.shared(value=cfval, name="conf", borrow=True)
		
		ball=ballzero+conf*movvec
		
		
		### deal with symmetry.. painfully...
		nsym=Transform.get_nsym(sym)
		if sym!="c1":
			
			if sym.startswith('d'):
				asym=[]
				for i in range(nsym/2):	
					a=6.28/(nsym/2)*i
					asym.append(T.stacklists(
						[ball[:,0]*T.cos(a)-ball[:,1]*T.sin(a),
						 ball[:,0]*T.sin(a)+ball[:,1]*T.cos(a),
						 ball[:,2]]))
					asym.append(T.stacklists(
						[ball[:,0]*T.cos(a)-ball[:,1]*T.sin(a),
						 -(ball[:,0]*T.sin(a)+ball[:,1]*T.cos(a)),
						 -ball[:,2]]))
				#print asym[0].T.eval()
				balls=T.concatenate(asym,axis=1).T
				nballs*=nsym
			
			if sym.startswith('c'):
				asym=[]
				for i in range(nsym):	
					a=6.28/(nsym)*i
					asym.append(T.stacklists(
						[ball[:,0]*T.cos(a)-ball[:,1]*T.sin(a),
						 ball[:,0]*T.sin(a)+ball[:,1]*T.cos(a),
						 ball[:,2]]))
				#print asym[0].T.eval()
				balls=T.concatenate(asym,axis=1).T
				nballs*=nsym
			ww=[wts for i in range(nsym)]
			wtss=T.concatenate(ww,axis=0)
		else:
			balls=ball
			wtss=wts
		print balls.shape.eval()
		print wtss.shape.eval()
		#numpy2pdb(balls.eval(), "tmp.pdb")
		#exit()

		### get the 3d density map for initial tunning
		ind_np=np.indices((sz,sz,sz)).astype(theano.config.floatX)
		ind_np=np.transpose(ind_np,axes=(3,2,1,0))
		ind=theano.shared(value=ind_np,borrow=True)
		def make_3d(p,w,den,ind):
			d=(ind-p)**2
			v=w*T.exp(-T.sum(d,axis=3)/(model["width"]))
			den+=v
			return den
		
		map_3d_all,update=theano.scan(fn=make_3d,
				outputs_info=T.zeros((sz,sz,sz)),
				sequences=[balls+sz/2,wtss],
				non_sequences=ind,
				)
		map_3d=map_3d_all[-1]
		#map_3d=map_3d.dimshuffle([2,1,0])

		self.target_map=T.tensor3('tar_map')
		self.map_err=T.sum((self.target_map-map_3d)**2)
		map_grad_w=T.grad(self.map_err, weight)
		self.map_update_w=[(weight,weight-map_grad_w*learnrate)]
		map_grad_p=T.grad(self.map_err, ballzero)
		self.map_update_p=[(ballzero,ballzero-map_grad_p*learnrate)]
		
		### make rotation matrix
		azp=ang[2]+3.14/2
		altp=3.14-ang[1]
		phip=6.28-ang[0]

		matrix=[(T.cos(phip)*T.cos(azp) - T.cos(altp)*T.sin(azp)*T.sin(phip)),
		(T.cos(phip)*T.sin(azp) + T.cos(altp)*T.cos(azp)*T.sin(phip)),
		(T.sin(altp)*T.sin(phip)),
		
		(-T.sin(phip)*T.cos(azp) - T.cos(altp)*T.sin(azp)*T.cos(phip)),
		(-T.sin(phip)*T.sin(azp) + T.cos(altp)*T.cos(azp)*T.cos(phip)),
		(T.sin(altp)*T.cos(phip)),
		
		(T.sin(altp)*T.sin(azp)),
		(-T.sin(altp)*T.cos(azp)),
		T.cos(altp)]

		mat=T.stacklists(matrix).T.reshape((3,3))
		newpos=T.dot(balls,mat)#+sz/2

		tx=ang[5]
		ty=ang[4]

		newpos=T.inc_subtensor(newpos[:,0],tx)
		newpos=T.inc_subtensor(newpos[:,1],ty)

		mirror=ang[3]
		newpos=T.set_subtensor(newpos[:,1], newpos[:,1]*mirror)
		newpos=newpos+sz/2
		#newpos[:,1]+=ty
		
		
		### gird for 2d images
		grid_x =T.arange(sz).astype(theano.config.floatX)
		grid_x=grid_x.repeat(sz,axis=0).reshape((sz,sz))
		grid_y=grid_x.copy().T
		
		### now make 2d projections
		
		iy=T.arange(nballs)

		def make_img(iy,r,x,y,pos, wt):
			ret=r+ wt[iy] *T.exp((-(x-pos[iy,0])**2 -(y-pos[iy,1])**2)/(15))
			return ret
			
		img,update=theano.scan(fn=make_img,
			outputs_info=T.zeros((sz,sz)),
			sequences=[iy],
			non_sequences=[grid_x,grid_y,newpos,wtss],
			)
		out=img[-1]#[:-1]
		
		
		
		xx=out.flatten()
		zz=imgin.flatten()
		L = T.sum(xx*zz)/T.sqrt(T.sum(xx**2))/T.sqrt(T.sum(zz**2))
		cost=-L#T.mean(L)
		grad_conf=T.grad(cost, conf)
		self.update_conf=[(conf, conf-learnrate*grad_conf)]

		self.grad_ballpos=T.grad(cost, ballzero)
		self.imgin=imgin
		self.orientin=ang
		self.balls=balls
		self.cost=cost
		self.learnrate=learnrate
		self.map3d=map_3d
		self.conf=conf
		#self.updates=updates
		#self.updates_ang=updates_ang
		self.out=out
		
		self.weight=weight
		self.ballzero=ballzero
		self.movvec=movvec
		
		
		
		#get_img=theano.function([],out)
		#get_pos=theano.function([],newpos)
		#get_cost=theano.function([imgin], cost)
		#get_gradconf=theano.function([imgin], grad_conf)
		#train_conf=theano.function([learnrate,imgin], cost, updates=update_conf)

		
		
if __name__ == '__main__':
	main()
	
