#!/usr/bin/env python

# Muyuan Chen 2016-01
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#
from EMAN2 import *
import numpy as np
#from multiprocessing import Pool
import threading
from Queue import Queue
#from EMAN2jsondb import JSTask
import time
import json
from shutil import copyfile

def main():
	progname = os.path.basename(sys.argv[0])
	usage="""prog --model model1.hdf,model2.hdf --oldpath refine_01
	Perform a 3d classification like e2refine_multi using the orientation of each particle in an e2refine_easy"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--newpath", type=str,help="Path to the classified results. Default = multinoali_XX", default=None)
	parser.add_argument("--oldpath", type=str,help="Path to the original refinement", default=None,guitype='filebox', filecheck=False,browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=2, col=0, rowspan=1, colspan=3)
	parser.add_argument("--models","--model", dest="model", type=str,help="Comma separated list of reference maps used for classification. If a single map is provided, data will be split into two groups based on similarity to the single map.", default=None,guitype='filebox', browser='EMModelsTable(withmodal=True,multiselect=True)', filecheck=False, row=7, col=0, rowspan=1, colspan=3)
	parser.add_argument("--simcmp",type=str,help="The name of a 'cmp' to be used in comparing the aligned images. eg- frc:minres=80:maxres=20. Default=ccc", default="ccc", guitype='strbox', row=10, col=0, rowspan=1, colspan=3)
	parser.add_argument("--threads", type=int,help="Number of threads.", default=4, guitype='intbox', row=12, col=0, rowspan=1, colspan=1)
	parser.add_argument("--parallel", type=str,help="Parallel option.", default="thread:4", guitype='strbox', row=13, col=0, rowspan=1, colspan=1)
	parser.add_argument("--iter", type=int,help="Number of iterations.", default=1, guitype='intbox', row=12, col=1, rowspan=1, colspan=1)
	parser.add_header(name="optheader", help='Optional parameters:', title="Optional:", row=14, col=0, rowspan=1, colspan=3)
	parser.add_argument("--mask",type=str,help="Name of an optional mask file. The mask is applied to the input models to focus the classification on a particular region of the map. Consider e2classifyligand.py instead.", default=None,guitype='filebox', browser='EMModelsTable(withmodal=True,multiselect=False)', filecheck=False, row=15, col=0, rowspan=1, colspan=3)
	parser.add_argument("--breaksym", action="store_true", default=False ,help="breaksym", guitype='boolbox', row=16, col=0, rowspan=1, colspan=1)
	parser.add_argument("--symcopy", action="store_true", default=False ,help="symcopy", guitype='boolbox', row=17, col=0, rowspan=1, colspan=1)
	parser.add_argument("--nomask", action="store_true", default=False ,help="no mask", guitype='boolbox', row=16, col=1, rowspan=1, colspan=1)

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	from EMAN2PAR import EMTaskCustomer
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)

	if not options.model:
		print "No model input. Exit."
		exit()

	inputmodel=options.model.split(',')
	modelstack=0
	if len(inputmodel)==1:
		num=EMUtil.get_image_count(inputmodel[0])
		if num>1:
			modelstack=num
			print "3D stack input. Perform multi-model refinement using existing alignment..."
			multimodel=True
			inputmodel=inputmodel*num
		else:

			multimodel=False
			print "One input model. Split the data by half accroding to the similarity to the input model..."
	else:
		multimodel=True
		print "Multiple input models. Perform multi-model refinement using existing alignment..."

	### make new folder
	if options.newpath == None:
		fls=[int(i[-2:]) for i in os.listdir(".") if i[:11]=="multinoali_" and len(i)==13 and str.isdigit(i[-2:])]
		if len(fls)==0 : fls=[0]
		options.newpath = "multinoali_{:02d}".format(max(fls)+1)

	print "Working directory: {}".format(options.newpath)
	try: os.mkdir(options.newpath)
	except:
		print "New path {} exist. Overwrite...".format(options.newpath)
		pass


	### read json file and parse some parameters
	with open(options.oldpath+"/0_refine_parms.json") as json_file:
		db = json.load(json_file)
	db=parse_json(db.copy())
	options.simcmp=parsemodopt(options.simcmp)
	if options.symcopy:
		options.breaksym=True
	sym=str(db["sym"])
	if db["breaksym"]:
		sym="c1"
	### copy the model to the new folder
	print "Preprocessing the input models..."
	#if options.mask:
		#options.mask="--multfile {}".format(options.mask)
	#else:
		#options.mask=""

	db_apix=db["apix"]
	if db_apix==0:
		e=EMData(inputmodel[0],0,True)
		db_apix=e["apix_x"]
	
	if options.nomask:
		automask="mask.soft:outer_radius=-1"
	else:
		automask=db["automask3d"]
	if multimodel:

		if modelstack>1:
			models=range(modelstack)
			for m in range(modelstack):
				outfile="{path}/model_input_{k}.hdf".format(path=options.newpath, k=m)
				run("e2proc3d.py {model} {out} --process=filter.lowpass.randomphase:cutoff_freq={freq} --apix={apix} --first {mi} --last {mi}".format(model=inputmodel[m],out=outfile,freq=1.0/(db["targetres"]*2),apix=db_apix, mi=m))
				inputmodel[m]=outfile
		else:

			models=range(len(inputmodel))
			for m in models:
				outfile="{path}/model_input_{k}.hdf".format(path=options.newpath, k=m)
				run("e2proc3d.py {model} {out} --process=filter.lowpass.randomphase:cutoff_freq={freq} --apix={apix} ".format(model=inputmodel[m],out=outfile,freq=1.0/(db["targetres"]*2),apix=db_apix))
				inputmodel[m]=outfile


	else:
		models=[0,1]
		outfile="{path}/model_input.hdf".format(path=options.newpath)
		run("e2proc3d.py {model} {out} --process=filter.lowpass.randomphase:cutoff_freq={freq} --apix={apix} ".format(model=inputmodel[0],out=outfile,freq=1.0/(db["targetres"]*2),apix=db_apix))
		inputmodel[0]=outfile

	output_3d=[]
	output_cls=[]
	input_eo_order={0:"even",1:"odd"}
	
	#### deal with breaksym
	if options.breaksym:
		if sym=="c1":
			print "Error: Breaksym option should be used based on a refinement with (non-c1) symmetry without breaksym. Exit."
			exit()
		if not multimodel:
			print "Breaksym does not support single model input yet...Exit."
			exit()
		print "Breaking symmetry..."
		#### original orientations
		origen=str(db["orientgen"])
		sym_object = parsesym(sym)
		[og_name,og_args] = parsemodopt(origen)
		oris0 = sym_object.gen_orientations(og_name, og_args)
		
		#### orientations with breaksym
		origen=origen+":breaksym=1:breaksym_real=1"
		sym_object = parsesym(sym)
		[og_name,og_args] = parsemodopt(origen)
		oris = sym_object.gen_orientations(og_name, og_args)
		
		nsym=Transform.get_nsym(sym)
		n=len(oris)
		n0=len(oris0)
		symdone=[-1]*n
		eulerlst=[]
		for i in range(n0):
			o=oris0[i]
			elst=[i]
			symdone[i]=i
			for s in range(1,nsym):
				ss=o.get_sym(sym, s)
				for j in range(n):
					if symdone[j]>=0:
						continue
					if oris[j]==ss:
						elst.append(j)
						symdone[j]=i
						break
			eulerlst.append(elst)
		print "Making {} projections in {} groups...".format(n, len(eulerlst))
		sym="c1"

		#for i in range(len(eulerlst)):
			#print i, eulerlst[i]
		#print symdone
		#exit()
	else:
		origen=db["orientgen"]
			
	if options.mask:
		mask3d="{path}/simmask.hdf".format(path=options.newpath)
		#maskfile="{path}/mask_projection.hdf".format(path=options.newpath)
		run("e2proc3d.py {msk0} {msk1} --apix {apix}".format(msk0=options.mask, msk1=mask3d, apix=db_apix))
		mskmodel=[]
		for m in inputmodel:
			mskmodel.append(m[:-4]+"_msk.hdf")
			run("e2proc3d.py {model} {mskmodel} --multfile {msk}".format(model=m, mskmodel=mskmodel[-1], msk=mask3d))
		run("e2project3d.py {model} --outfile  {mskfile} -f --orientgen {orient} --sym {sym} --parallel {para}".format(model=mask3d,mskfile=maskfile,orient=origen,sym=db["sym"],para=options.parallel))
		#simmask=EMData.read_images(maskfile)
	
	#### start iterations...
	for it in range(options.iter):
		print "Starting iteration {} ...".format(it)
		print "Making projections..."

		if it==0:
		#### first iteration. do one projection for even/odd
			if multimodel:
				projfile=[]
				
				for m in models:
					projfile.append("{path}/projections_{it:02d}_{k}.hdf".format(path=options.newpath, k=m, it=it))
					run("e2project3d.py {model}  --outfile {proj} -f --orientgen {orient} --sym {sym} --parallel {para}".format(		model=inputmodel[m],proj=projfile[-1],orient=origen,sym=db["sym"],para=options.parallel))
				#### make another set of projections for masked model
				if options.mask:
					projmskfile=[]
					for m in models:
						projmskfile.append("{path}/projections_msk_{it:02d}_{k}.hdf".format(path=options.newpath, k=m, it=it))
						run("e2project3d.py {model}  --outfile {proj} -f --orientgen {orient} --sym {sym} --parallel {para}".format(		model=mskmodel[m],proj=projmskfile[-1],orient=origen,sym=db["sym"],para=options.parallel))
						
			else:
				projfile=["{path}/projections_{it:02d}.hdf".format(path=options.newpath, it=it)]
				run("e2project3d.py {model}  --outfile {proj} -f --orientgen {orient} --sym {sym} --parallel {para}".format(		model=inputmodel[0],proj=projfile[0],orient=origen,sym=db["sym"],para=options.parallel))
				
				#### make another set of projections for masked model
				if options.mask:
					projmskfile=["{path}/projections_msk_{it:02d}.hdf".format(path=options.newpath, it=it)]
					run("e2project3d.py {model}  --outfile {proj} -f --orientgen {orient} --sym {sym} --parallel {para}".format(		model=mskmodel[0],proj=projmskfile[0],orient=origen,sym=db["sym"],para=options.parallel))

		output_3d.append({})
		output_cls.append({})
		### even/odd loop
		for eoid,eo in input_eo_order.items():

			if it>0:
				inputmodel=[output_3d[-2][eo][m] for m in models]
				print inputmodel
				multimodel=True
			#### make projections for even/odd
				projfile=["{path}/projections_{it:02d}_{k}_{eo}.hdf".format(path=options.newpath, k=m, it=it,eo=eo) for m in models]
				for m in models:
					run("e2project3d.py {model}  --outfile {proj} -f --orientgen {orient} --sym {sym} --parallel {para}".format(		model=inputmodel[m],proj=projfile[m],orient=origen,sym=db["sym"],para=options.parallel))
				if options.mask:
					mskmodel=[]
					projmskfile=[]
					for m in models:
						mskmodel.append(inputmodel[m][:-4]+"_msk.hdf")
						projmskfile.append("{path}/projections_msk_{it:02d}_{k}_{eo}.hdf".format(path=options.newpath, k=m, it=it,eo=eo))
						run("e2proc3d.py {model} {mskmodel} --multfile {msk}".format(model=inputmodel[m], mskmodel=mskmodel[-1], msk=mask3d))
						run("e2project3d.py {model}  --outfile {proj} -f --orientgen {orient} --sym {sym} --parallel {para}".format(		model=mskmodel[m],proj=projmskfile[-1],orient=origen,sym=db["sym"],para=options.parallel))

			oldmapfile=str(db["last_{}".format(eo)])
			ptclfile=str(db["input"][eoid])

			clsmx=oldmapfile.replace("threed","classmx")

			### old projection file is used for classaverage alignment
			oldprojfile=oldmapfile.replace("threed","projections")

			ncls=EMUtil.get_image_count(projfile[0])
			npt=EMUtil.get_image_count(ptclfile)
			newclsmx=["{path}/classmx_{it:02d}_{n}_{eo}.hdf".format(path=options.newpath,n=i,eo=eo,it=it) for i in models]
			classout=["{path}/classes_{it:02d}_{n}_{eo}.hdf".format(path=options.newpath,n=i,eo=eo,it=it) for i in models]
			threedout=["{path}/threed_{it:02d}_{n}_{eo}.hdf".format(path=options.newpath,n=i,eo=eo,it=it) for i in models]
			output_3d[-1][eo]=threedout
			output_cls[-1][eo]=classout

			#### get alignment from classmx file and calculate similarity
			print "Calculating similarity matrix..."
			cmxcls=EMData(clsmx,0)
			cmxtx=EMData(clsmx,2)
			cmxty=EMData(clsmx,3)
			cmxalpha=EMData(clsmx,4)
			cmxmirror=EMData(clsmx,5)

			projs=[]
			if options.mask:
				projs=[EMData.read_images(pj) for pj in projmskfile]
			else:
				projs=[EMData.read_images(pj) for pj in projfile]
			#for pj in projfile:
			#	projs.append(EMData.read_images(pj))
			#	if options.mask:
			#		for si in range(len(simmask)):
			#			projs[-1][si].mult(simmask[si])
			
			xforms=[]
			tasks=[]
			for i in range(npt):
				c=int(cmxcls[0,i])
				tr=Transform({"type":"2d","alpha":cmxalpha[0,i],"mirror":int(cmxmirror[0,i]),"tx":cmxtx[0,i],"ty":cmxty[0,i]})
				if options.breaksym:
					elst=eulerlst[symdone[c]]
					pjs=[]
					for k in range(len(projfile)):
						pjs.extend([projs[k][ec] for ec in elst])
				else:
					pjs=[projs[k][c] for k in range(len(projfile))]
				
				xforms.append({"ptclfile":ptclfile,"proj":pjs,"idx":i,"xform":tr,"cmp":options.simcmp})
				#tasks.append(MultiNoaliTask(xforms[-1]))
			t00=time.time()
			t01=t00
			
			#### the threading part is copied from e2spt_align
			jsd=Queue(0)
			NTHREADS=max(options.threads+1,2)
			thrds=[threading.Thread(target=do_compare,args=(jsd,xforms[i])) for i in xrange(npt)]
			thrtolaunch=0
			corr=[[] for i in range(npt)]
			while thrtolaunch<len(thrds) or threading.active_count()>1:
				# If we haven't launched all threads yet, then we wait for an empty slot, and launch another
				# note that it's ok that we wait here forever, since there can't be new results if an existing
				# thread hasn't finished.
				if thrtolaunch<len(thrds) :
					while (threading.active_count()==NTHREADS ) : time.sleep(.1)
					if time.time()-t01>3:
						print "Starting thread {}/{}".format(thrtolaunch,len(thrds))
						t01=time.time()
					thrds[thrtolaunch].start()
					thrtolaunch+=1
				else: time.sleep(1)
			
				while not jsd.empty():
					idx,ccc=jsd.get()
					#print idx,ccc
					corr[idx]=ccc
			
			print time.time()-t00
			np.savetxt("{path}/simmx_{it:02d}_{eo}.txt".format(path=options.newpath,eo=eo, it=it),corr)
			
			
			#corr=np.loadtxt("{path}/simmx_00_{eo}.txt".format(path=options.newpath,eo=eo))
			#corr=[[c] for c in corr]

			### classification
			print "Classifying particles..."
			if options.symcopy:
				cmxy=cmxcls["ny"]
				cmxtmp=EMData(nsym, cmxy)
			else:
				cmxtmp=cmxcls.copy()
			cmxtmp.to_zero()
			cmxtmp.sub(1)
			cmxout=[cmxtmp.copy() for s in models]



			if multimodel:
				if options.symcopy:
					for i in range(npt):
						c=np.asarray(corr[i])
						cls=np.argmin(c.reshape(len(corr[i])/len(models),len(models)),1)
						#### the class it belongs to
						v=cmxcls[0,i] 
						#### index of the class of each asym-unit
						el=np.asarray(eulerlst[symdone[int(v)]]) 
						for u in range(len(cls)):
							cmxout[int(cls[u])][u,i]=int(el[u])
						#exit()
					
				else:
					### simply classify
					clso=np.argmin(corr,1)
					cls=clso%len(models)
					clsm=clso//len(models)
					print eo,[float(sum(cls==k))/float(npt) for k in models]
					for i in range(npt):
						v=cmxcls[0,i]
						if options.breaksym:
							#print v,cls[i],clsm[i],symdone[int(v)]
							v=eulerlst[symdone[int(v)]][clsm[i]]
							
						for s in models:
							if s==cls[i]:
								cmxout[s][0,i]=int(v)
							else:
								cmxout[s][0,i]=-1
			else:
				### one model input, split the data to two halves
				for c in range(ncls):
					ss=[]
					ns=0
					for i in range(npt):
						v=cmxcls[0,i]
						cc=corr[i]
						#print cc
						if options.breaksym:
							v=eulerlst[symdone[int(v)]][np.argmin(cc)]
							cc=[np.min(cc)]
						if v==c:
							ss.extend(cc)
							ns+=1
						else:
							ss.extend([10]*len(cc))

					### split the data by half
					spt=int(ns*.5)
					for s in models:
						if s==0:
							toavg=np.argsort(ss)[:spt]
						else:
							toavg=np.argsort(ss)[spt:ns]
						#print toavg
						for i in toavg:
							cmxout[s][0,int(i)]=int(c)
					print c, ncls,cmxout[s]["mean_nonzero"]

			
			### write classmx
			for s in models:
				#print cmxout[s]["maximum"]
				cmxout[s].write_image(newclsmx[s])
				ns=EMUtil.get_image_count(clsmx)
				for i in range(1,ns):
					e=EMData(clsmx,i)
					if options.symcopy:
						enp=e.numpy()[:,0]
						anp=np.tile(enp,[nsym,1]).T.copy()
						a=from_numpy(anp)
						a.write_image(newclsmx[s],i)
					else:	
						e.write_image(newclsmx[s],i)
			#exit()
			print "Making class average and 3d map..."
			if len(projfile)==1: 
				tmpprojfile=projfile*len(models)
			else:
				tmpprojfile=projfile
			for s in models:
				### class average
				run("e2classaverage.py --input {inputfile} --classmx {clsmx} --decayedge --storebad --output {clsout} --ref {proj} --iter {classiter} -f --normproc {normproc} --averager {averager} {classrefsf} {classautomask} --keep {classkeep} {classkeepsig} --cmp {classcmp} --align {classalign} --aligncmp {classaligncmp} {classralign} {prefilt} --parallel {para}".format(
					inputfile=ptclfile, clsmx=newclsmx[s], clsout=classout[s], proj=tmpprojfile[s], classiter=db["classiter"], normproc=db["classnormproc"], averager=db["classaverager"], classrefsf=db["classrefsf"],
					classautomask=db["classautomask"],classkeep=db["classkeep"], classkeepsig=db["classkeepsig"], classcmp=db["classcmp"], classalign=db["classalign"], classaligncmp=db["classaligncmp"],
					classralign=db["classralign"], prefilt=db["prefilt"], para=options.parallel))

				### make 3d
				run("e2make3dpar.py --input {clsout} --sym {sym} --output {threed} {preprocess} --keep {m3dkeep} {keepsig} --apix {apix} --pad {m3dpad} --mode gauss_5 --threads {threads} ".format(
				clsout=classout[s],threed=threedout[s], sym=sym, recon=db["recon"], preprocess=db["m3dpreprocess"],  m3dkeep=db["m3dkeep"], keepsig=db["m3dkeepsig"],
				m3dpad=db["pad"],threads=options.threads, apix=db_apix))
		#exit()
		### post process
		print "Post processing..."
		if os.path.exists("strucfac.txt") :
			m3dsetsf="--setsf strucfac.txt"
		else:
			m3dsetsf=""

		for s in models:
			final3d="{path}/threed_{it:02d}_{n}.hdf".format(path=options.newpath,n=s, it=it)
			
			run("e2refine_postprocess.py --even {even3d} --odd {odd3d} --output {final3d} --automaskexpand {amaskxp} --align --mass {mass} --iter {it} {amask3d} {amask3d2} {m3dpostproc} {setsf} --sym={sym} --restarget={restarget} --underfilter".format(it=it,even3d=output_3d[-1]["even"][s], odd3d=output_3d[-1]["odd"][s], final3d=final3d, mass=db["mass"], amask3d=automask, sym=sym, amask3d2=db["automask3d2"], m3dpostproc=db["m3dpostprocess"], setsf=m3dsetsf,restarget=db["targetres"], amaskxp=db.setdefault("automaskexpand","0")))

			### copy the fsc files..
			fscs=["fsc_unmasked_{:02d}.txt".format(it),"fsc_masked_{:02d}.txt".format(it),"fsc_maskedtight_{:02d}.txt".format(it)]
			for fsc in fscs:
				fm=os.path.join(options.newpath,fsc)
				fmnew=os.path.join(options.newpath,fsc[:-4]+"_model_{:02d}.txt".format(s))
				try:
					copyfile(fm,fmnew)
					os.remove(fm)
				except: pass

			if it==options.iter-1:
				### make lists
				tmpcls=["tmpcls_even.lst","tmpcls_odd.lst"]
				tmpcls_m=[l.replace('.','_m1.') for l in tmpcls]
				run("e2classextract.py {clsfile} --refinemulti --setname {tmpcls}".format(clsfile=output_cls[-1]["even"][s],tmpcls=tmpcls[0]))
				run("e2classextract.py {clsfile} --refinemulti --setname {tmpcls}".format(clsfile=output_cls[-1]["odd"][s],tmpcls=tmpcls[1]))
				lstout="sets/{}_{}.lst".format(options.newpath,s)
				run("e2proclst.py {lst1} {lst2} --mergesort {lstout}".format(lst1=tmpcls_m[0], lst2=tmpcls_m[1], lstout=lstout))
				for l in tmpcls_m:
					try: os.remove(l)
					except: pass


	E2end(logid)

def parse_json(db):
	if db["classrefsf"] : db["classrefsf"]="--setsfref"
	else: db["classrefsf"]=""
	if db["classautomask"] : db["classautomask"]="--automask"
	else: db["classautomask"]=""
	if db["classkeepsig"] : db["classkeepsig"]="--keepsig"
	else: db["classkeepsig"]=""
	if db["prefilt"] : db["prefilt"]="--prefilt"
	else: db["prefilt"]=""
	if db["m3dpreprocess"]: db["m3dpreprocess"]="--preprocess "+ str(db["m3dpreprocess"])
	else: db["m3dpreprocess"]=""
	if db["m3dkeepsig"] : db["m3dkeepsig"]="--keepsig"
	else: db["m3dkeepsig"]=""
	if db["automask3d"]==None : db["automask3d"]=""
	else : db["automask3d"]="--automask3d "+str(db["automask3d"])
	if db["automask3d2"]==None : db["automask3d2"]=""
	else : db["automask3d2"]="--automask3d2 "+str(db["automask3d2"])
	if db["m3dpostprocess"]==None : db["m3dpostprocess"]=""
	else : db["m3dpostprocess"]="--m3dpostprocess "+ str(db["m3dpostprocess"])
	if db["classralign"]:
		db["classralign"]="--ralign {} --raligncmp {}".format(db["classralign"],db["classraligncmp"])
	else: db["classralign"]=""
	if db["classiter"]<0: db["classiter"]=0

	return db

def do_compare(jsd,data):

	e=EMData(data["ptclfile"],data["idx"])
	e.transform(data["xform"])
	ret=[]
	for pj in data["proj"]:
		ret.append(e.cmp(data["cmp"][0],pj,data["cmp"][1]))
	
	jsd.put((data["idx"],ret))
	return 

def run(cmd):
	print cmd
	launch_childprocess(cmd)


#class MultiNoaliTask(JSTask):
	
	#def __init__(self,data):

		#JSTask.__init__(self,"MultiNoaliTask",data,{},"")
	
	
	#def execute(self,callback=None):
		
		#data=self.data
		#e=EMData(data["ptclfile"],data["idx"])
		#e.transform(data["xform"])
		#ret=[]
		#for pj in data["proj"]:
			#ret.append(e.cmp(data["cmp"][0],pj,data["cmp"][1]))
		#return ret
		
		


if __name__ == '__main__':
	main()
