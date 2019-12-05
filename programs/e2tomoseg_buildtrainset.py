#!/usr/bin/env python
from __future__ import division
# Muyuan Chen 2015-03
from past.utils import old_div
from builtins import range
from EMAN2 import *
import numpy as np
import random
from eman2_gui.emapplication import EMApp
from eman2_gui.emimage import EMImageWidget
from eman2_gui.emboxerbase import *
import shutil


def main():
	
	usage="Generate training set for tomogram segmentation. Please run this program from the GUI in e2projectmanager.py."
	#print usage
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	#parser.add_header(name="tmpheader", help='temp label', title="### This program is NOT avaliable yet... ###", row=0, col=0, rowspan=1, colspan=2, mode="box,seg,set")
	#### boxing ####
	parser.add_argument("--boxing",action="store_true",help="Boxing particles.",default=False, guitype='boolbox', row=4, col=0, rowspan=1, colspan=1, mode='box[True]')
	parser.add_pos_argument(name="tomograms",help="List the file to process with e2spt_boxer here.", default="", guitype='filebox', browser="EMTomoTable(withmodal=True,multiselect=True,startpath=\"tomograms\")",  row=1, col=0,rowspan=1, colspan=3, mode="box")
	parser.add_argument("--boxsize","-B",type=int,help="Box size in pixels",default=-1, guitype='intbox', row=3, col=0, rowspan=1, colspan=3, mode="box")
	
	#### segment ####
	#parser.add_header(name="instruction", help='instruction', title="### Mark the target features white ###", row=0, col=0, rowspan=1, colspan=2, mode="seg")
	parser.add_argument("--segment",action="store_true",help="Segment particles.",default=False, guitype='boolbox', row=4, col=0, rowspan=1, colspan=1, mode='seg[True]')
	parser.add_pos_argument(name="particles",help="Particle file.", default="", guitype='filebox', browser="EMParticlesTable(withmodal=True)",  row=1, col=0,rowspan=1, colspan=3, mode="seg")
	parser.add_argument("--output", type=str,help="output file name. Default is the input particle file name plus _seg.hdf", default=None,guitype='strbox', row=3, col=0, rowspan=1, colspan=3, mode="seg")
	
	#### build set ####
	parser.add_argument("--buildset",action="store_true",help="Segment particles.",default=False, guitype='boolbox', row=7, col=0, rowspan=1, colspan=1, mode='set[True]')
	parser.add_argument("--particles_raw", type=str,help="Input raw particle file", default=None,guitype='filebox',browser="EMParticlesTable(withmodal=True,multiselect=True)", row=1, col=0, rowspan=1, colspan=3, mode="set")
	parser.add_argument("--particles_label", type=str,help="Input labels for particle file", default=None,guitype='filebox',browser="EMParticlesTable(withmodal=True,multiselect=True)", row=2, col=0, rowspan=1, colspan=3, mode="set")
	parser.add_argument("--boxes_negative", type=str,help="Input boxes of negative samples", default=None,guitype='filebox',browser="EMParticlesTable(withmodal=True,multiselect=True)", row=3, col=0, rowspan=1, colspan=3, mode="set")
	parser.add_argument("--ncopy",type=int,help="Number of copies for NEGATIVE samples. (number of copies of particles is calculated accordingly) ",default=10, guitype='intbox', row=5, col=0, rowspan=1, colspan=1, mode="set")
	parser.add_argument("--trainset_output", type=str,help="output file name of the training set.Default is the input particle file name plus _trainset.hdf", default=None,guitype='strbox', row=4, col=0, rowspan=1, colspan=3, mode="set")
	#parser.add_argument("--zthick",type=int,help="Thickness in z ",default=0, guitype='intbox', row=5, col=1, rowspan=1, colspan=1, mode="set")
	parser.add_argument("--validset",type=float,help="Propotion of particles in validation set. Default is 0.2 ",default=0.0, guitype='floatbox', row=7, col=1, rowspan=1, colspan=1, mode="set")

	##################
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	#### boxing ###
	if options.boxing:
		db = js_open_dict(EMBOXERBASE_DB)
		cache_box_size = True
		if options.boxsize == -1:
			cache_box_size = False
			options.boxsize = db.setdefault("box_size",64)

		application = EMApp()
		module = EMBoxerModule(args,options.boxsize)
		module.show_interfaces()

		application.execute(logid)
	
	
	#### segment ###
	if options.segment:
		filename=args[0]
		if options.output==None:
			options.output=filename[:-4]+"_seg.hdf"
		app = EMApp()
		img=EMData.read_images(filename)
		#print img[0]["mean"]
		w=EMImageWidget(data=img,old=None,app=app,force_2d=True)
		#w = EMWidgetFromFile(filename,application=app,force_2d=True)
		w.setWindowTitle(base_name(filename))
		w.set_scale(2.5)
		w.resize(400,400)
		w.show_inspector(1)
		ins=w.get_inspector()
		ins.mmtab.setCurrentIndex(5)
		ins.dtpenv.setText('100')
		w.set_mouse_mode(5)
		app.show_specific(w)
		app.exec_()
		try: os.remove(options.output)
		except:pass
		for e in img:
			e.process_inplace("threshold.belowtozero", {"minval":99})
			e.process_inplace("threshold.binary", {"value":1})
			e.write_image(options.output,-1)
		print("Segmented particles saved to {}.".format(options.output))
			
	#### build set ###
	if options.buildset:
		tomo_in=options.particles_raw
		seg_in=options.particles_label
		
		if tomo_in==None:
			print("particle_raw input required. exit.")
			return
		
		tomo_in=tomo_in.split(',')
		if seg_in==None:
			seg_in=[]
			print("Searching for particle label...")
			for t in tomo_in:
				s=t[:-4]+"_seg.hdf"
				if os.path.isfile(s):
					seg_in.append(s)
					print("\t{} -> {}".format(t, s))
				else:
					print("cannot find {}. exit".format(s))
					return
		
		else:
			seg_in=seg_in.split(',')
			if len(seg_in)!=len(tomo_in):
				print("number of particles does not match number of labels. exit")
				return
		
		
		neg_in=options.boxes_negative
		if neg_in==None:
			neg_in=[]
			p_copy=options.ncopy
			print("no negative particle input, making {} copies for particles".format(p_copy))
			
		else:
			neg_in=neg_in.split(',')
			
			n_ptcl=np.sum([EMUtil.get_image_count(t) for t in tomo_in])
			n_neg=np.sum([EMUtil.get_image_count(t) for t in neg_in])
			
			p_copy=options.ncopy*n_neg//n_ptcl
			print("total {} positive samples, {} negative samples".format(n_ptcl, n_neg)) 
			print("making {} copies for particles, and {} copies for negative samples".format(p_copy,options.ncopy))
			
		if options.trainset_output==None:
			options.trainset_output=tomo_in[0][:-4]+"_trainset.hdf"
		
		print("training set writting to {}".format(options.trainset_output))
		if os.path.isfile(options.trainset_output):
			os.remove(options.trainset_output)
		
		imgs=[]
		for ti, tomo in enumerate(tomo_in):
			
			n_ptcl=EMUtil.get_image_count(tomo)
			seg=seg_in[ti]
			
			for i in range(n_ptcl):
				t=EMData(tomo,i)
				#t=get_box(tomo_in,i,options.zthick)
				if t==None: continue
				s=EMData(seg, i)
				for c in range(p_copy):
					tr=Transform()
					rd=random.random()*360
					tr.set_rotation({"type":"2d","alpha":rd})
					e=t.process("xform",{"transform":tr})
					#e.process_inplace("normalize")
					imgs.append(e)
					#e.write_image(options.trainset_output,-1)
					e=s.process("xform",{"transform":tr})
					#e.write_image(options.trainset_output,-1)
					imgs.append(e)
		ngood=len(imgs)
		
		for neg in neg_in:
			s=EMData(neg,0)
			s.to_zero()
			n_neg=EMUtil.get_image_count(neg)
			for i in range(n_neg):
				t=EMData(neg, i)
				#t=get_box(neg_in,i,options.zthick)
				if t==None: continue
				for c in range(options.ncopy):
					tr=Transform()
					rd=random.random()*360
					tr.set_rotation({"type":"2d","alpha":rd})
					e=t.process("xform",{"transform":tr})
					#e.process_inplace("normalize")
					#e.write_image(options.trainset_output,-1)
					imgs.append(e)
					e=s.process("xform",{"transform":tr})
					#e.write_image(options.trainset_output,-1)
					imgs.append(e)

		print("Shuffling particles...")
		### randomize
		n=len(imgs)
		
		#n=EMUtil.get_image_count(options.trainset_output)
		idx=list(range(int((old_div(ngood,2))*(1-options.validset))))+list(range((old_div(ngood,2)), int(old_div(n,2)*(1-options.validset))))
		random.shuffle(idx)
		for i in idx:
			imgs[i*2]["valid_set"]=0
			imgs[i*2].write_image(options.trainset_output,-1)
			imgs[i*2+1]["valid_set"]=0
			imgs[i*2+1].write_image(options.trainset_output,-1)
			
		idx=list(range(int((old_div(ngood,2))*(1-options.validset)), old_div(ngood,2)))+list(range( int(old_div(n,2)*(1-options.validset)), old_div(n,2)))
		random.shuffle(idx)
		for i in idx:
			imgs[i*2]["valid_set"]=1
			imgs[i*2].write_image(options.trainset_output,-1)
			imgs[i*2+1]["valid_set"]=1
			imgs[i*2+1].write_image(options.trainset_output,-1)
			#e=EMData(options.trainset_output,i*2)
			#e.process_inplace("normalize")
			#e.write_image(tmpfile,-1)
			#e=EMData(options.trainset_output,i*2+1)
			#e.write_image(tmpfile,-1)
		#shutil.move(tmpfile,options.trainset_output)
		print("Generate a training set of {:d} samples.".format(old_div(n,2)))
		
	print("Done")
	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
def get_box(fname, idx, nz):
	if nz==0:
		return EMData(fname,idx)
	else:
		#print nz
		sep=2
		e=EMData(fname,idx,True)
		src=e["ptcl_source_image"]
		box=e["ptcl_source_coord"]
		boxz=e["ptcl_source_coord_z"]
		sz=e["nx"]
		out=EMData(sz,sz,nz*2+1)
		outnp=out.numpy()
		
		hdr=EMData(src,0,True)
		zmax=hdr["nz"]
		if boxz<nz*sep or boxz>zmax-nz*sep:
			print("skipping box ",idx)
			return None
		
		c=EMData(src,0,False,Region(box[0]-old_div(sz,2),box[1]-old_div(sz,2),boxz,sz,sz,1))
		outnp[nz]=c.numpy().copy()
		
		for i in range(nz):
			c=EMData(src,0,False,Region(box[0]-old_div(sz,2),box[1]-old_div(sz,2),boxz+(i+1)*sep,sz,sz,1))
			outnp[nz+i+1]=c.numpy().copy()
			c=EMData(src,0,False,Region(box[0]-old_div(sz,2),box[1]-old_div(sz,2),boxz-(i+1)*sep,sz,sz,1))
			outnp[nz-i-1]=c.numpy().copy()
			
		return out
	
	
	
	
	
	
if __name__ == '__main__':
	main()
	
