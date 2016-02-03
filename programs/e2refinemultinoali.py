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
from multiprocessing import Pool
import time
import json

def main():
	
	usage="""muticlass_noalign.py 
	Perform a 3d classification like e2refine_multi using the orientation of each particle in an e2refine_easy"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--newpath", type=str,help="path for the classificaton", default="newrefine_01")
	parser.add_argument("--oldpath", type=str,help="path for the original refinement", default="refine_01")
	parser.add_argument("--model", type=str,help="model file for classification", default=None)
	parser.add_argument("--threads", type=int,help="threads", default=12)
	
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	if not options.model:
		print "No model input. Exit."
		exit()
		
	inputmodel=options.model.split(',')
	
	if len(inputmodel)==1:
		multimodel=False
		print "One input model. Split the data by half accroding to the similarity to the input model..."
	elif len(inputmodel)==2:
		multimodel=True
		print "Two input models. Perform multi-model refinement using existing alignment..."
	else:
		print "Not implemented yet..."
	
	try: os.mkdir(options.newpath)
	except: 
		print "New path exist. Overwrite..."
		pass
	
	
	### read json file and parse some parameters
	with open(options.oldpath+"/0_refine_parms.json") as json_file:    
		db = json.load(json_file)
	db=parse_json(db.copy())
	
	
	### copy the model to the new folder
	input_eo_order={0:"even",1:"odd"}
	if multimodel:
		models=range(len(inputmodel))
		
		for m in models:
			outfile="{path}/model_input_{k}.hdf".format(path=options.newpath, k=m)
			run("e2proc3d.py {model} {out} --process=filter.lowpass.randomphase:cutoff_freq={freq} --apix={apix}".format(model=inputmodel[m],out=outfile,freq=1.0/(db["targetres"]*2),apix=db["apix"]))
		
	
	else:
		models=[0,1]
		outfile="{path}/model_input.hdf".format(path=options.newpath)
		run("e2proc3d.py {model} {out} --process=filter.lowpass.randomphase:cutoff_freq={freq} --apix={apix}".format(model=inputmodel[0],out=outfile,freq=1.0/(db["targetres"]*2),apix=db["apix"]))
		
		
	### make projections first, and use this for both even and odd
	
	print "Making projections..."
	if multimodel:
		projfile=[]
		for m in models:
			projfile.append("{path}/projections_00_{k}.hdf".format(path=options.newpath, k=m))
			run("e2project3d.py {model}  --outfile {proj} -f --orientgen {orient} --sym {sym} --parallel thread:{threads}".format(		model=inputmodel[m],proj=projfile[-1],orient=db["orientgen"],sym=db["sym"],threads=options.threads))
	else:
		projfile=["{path}/projections_00.hdf".format(path=options.newpath)]
		run("e2project3d.py {model}  --outfile {proj} -f --orientgen {orient} --sym {sym} --parallel thread:{threads}".format(		model=inputmodel[0],proj=projfile[0],orient=db["orientgen"],sym=db["sym"],threads=options.threads))
	
	output_3d={}
	output_cls={}
	for eoid,eo in input_eo_order.items():
		
		oldmapfile=str(db["last_{}".format(eo)])
		ptclfile=str(db["input"][eoid])
		
		clsmx=oldmapfile.replace("threed","classmx")
		
		### old projection file is used for classaverage alignment
		oldprojfile=oldmapfile.replace("threed","projections") 
		
		ncls=EMUtil.get_image_count(projfile[0])
		npt=EMUtil.get_image_count(ptclfile)
		newclsmx=["{path}/classmx_00_{n}_{eo}.hdf".format(path=options.newpath,n=i,eo=eo) for i in models]
		classout=["{path}/classes_00_{n}_{eo}.hdf".format(path=options.newpath,n=i,eo=eo) for i in models]
		threedout=["{path}/threed_00_{n}_{eo}.hdf".format(path=options.newpath,n=i,eo=eo) for i in models]
		output_3d[eo]=threedout
		output_cls[eo]=classout

		### get alignment from classmx file and calculate similarity
		print "calculating similarity..."
		cmxcls=EMData(clsmx,0)
		cmxtx=EMData(clsmx,2)
		cmxty=EMData(clsmx,3)
		cmxalpha=EMData(clsmx,4)
		cmxmirror=EMData(clsmx,5)
		
		projs=[]
		for pj in projfile:
			projs.append(EMData.read_images(pj))
		xforms=[]
		for i in range(npt):
			c=int(cmxcls[0,i])
			tr=Transform({"type":"2d","alpha":cmxalpha[0,i],"mirror":int(cmxmirror[0,i]),"tx":cmxtx[0,i],"ty":cmxty[0,i]})
			pjs=[projs[k][c] for k in range(len(projfile))]
			xforms.append({"ptclfile":ptclfile,"proj":pjs,"idx":i,"xform":tr})
		
		pool = Pool()
		corr=pool.map_async(do_compare, xforms)
		pool.close()
		while (True):
			if (corr.ready()): break
			remaining = corr._number_left
			print "Waiting for", remaining, "tasks to complete..."
			time.sleep(2)
		corr=corr.get()
		np.savetxt("{path}/simmx_00_{eo}.txt".format(path=options.newpath,eo=eo),corr)
		
		#corr=np.loadtxt("{path}/simmx_00_{eo}.txt".format(path=options.newpath,eo=eo))
		### classification 
		print "Classifying..."
		cmxtmp=cmxcls.copy()
		cmxtmp.to_zero()
		cmxtmp.sub(1)
		cmxout=[cmxtmp.copy(), cmxtmp.copy()]
		
		
		
		if multimodel:
			cls=np.argmin(corr,1)
			print eo,[float(sum(cls==k))/float(npt) for k in models]
			for i in range(npt):
				v=cmxcls[0,i]
				for s in models:
					if s==cls[i]:
						cmxout[s][0,i]=v
					else:
						cmxout[s][0,i]=-1
		else:
			
			for c in range(ncls):
				ss=[]
				ns=0
				for i in range(npt):
					v=cmxcls[0,i]
					if v==c:
						ss.append(corr[i])
						ns+=1
					else:
						ss.append([10]*len(corr[i]))
				
				### split the data by halv
				spt=int(ns*.5)
				for s in models:
					if s==0:
						toavg=np.argsort(ss)[:spt]
					else:
						toavg=np.argsort(ss)[spt:ns]
				
					for i in toavg:
						cmxout[s][0,i]=c
			
		
		### write classmx	
		for s in models:
			cmxout[s].write_image(newclsmx[s])
			ns=EMUtil.get_image_count(clsmx)
			for i in range(1,ns):
				e=EMData(clsmx,i)
				e.write_image(newclsmx[s],i)
		
		print "class averaging..."
		for s in models:	
			### class average
			run("e2classaverage.py --input {inputfile} --classmx {clsmx} --decayedge --storebad --output {clsout} --ref {proj} --iter {classiter} -f --normproc {normproc} --averager {averager} {classrefsf} {classautomask} --keep {classkeep} {classkeepsig} --cmp {classcmp} --align {classalign} --aligncmp {classaligncmp} {classralign} {prefilt} --parallel thread:{thrd}".format(
				inputfile=ptclfile, clsmx=newclsmx[s], clsout=classout[s], proj=oldprojfile, classiter=db["classiter"], normproc=db["classnormproc"], averager=db["classaverager"], classrefsf=db["classrefsf"],
				classautomask=db["classautomask"],classkeep=db["classkeep"], classkeepsig=db["classkeepsig"], classcmp=db["classcmp"], classalign=db["classalign"], classaligncmp=db["classaligncmp"],
				classralign=db["classralign"], prefilt=db["prefilt"], thrd=options.threads))
			
			### make 3d
			run("e2make3dpar.py --input {clsout} --sym {sym} --output {threed} {preprocess} --keep {m3dkeep} {keepsig} --apix {apix} --pad {m3dpad} --mode gauss_5 --threads {threads} ".format(
			clsout=classout[s],threed=threedout[s], sym=db["sym"], recon=db["recon"], preprocess=db["m3dpreprocess"],  m3dkeep=db["m3dkeep"], keepsig=db["m3dkeepsig"],
			m3dpad=db["pad"],fillangle=db["astep"] ,threads=options.threads, apix=db["apix"]))
 
	### post process
	print "Post processing..."
	if os.path.exists("strucfac.txt") :
		m3dsetsf="--setsf strucfac.txt"
	else:
		m3dsetsf=""
		
	for s in models:
		final3d="{path}/threed_00_{n}.hdf".format(path=options.newpath,n=s)
		run("e2refine_postprocess.py --even {even3d} --odd {odd3d} --output {final3d} --automaskexpand {amaskxp} --align --mass {mass} --iter 0 {amask3d} {amask3d2} {m3dpostproc} {setsf} --sym={sym} --restarget={restarget} --underfilter".format(even3d=output_3d["even"][s], odd3d=output_3d["odd"][s], final3d=final3d, mass=db["mass"], amask3d=db["automask3d"], sym=db["sym"], amask3d2=db["automask3d2"], m3dpostproc=db["m3dpostprocess"], setsf=m3dsetsf,restarget=db["targetres"], amaskxp=db["automaskexpand"]))
		
		### make lists
		tmpcls=["tmpcls_even.lst","tmpcls_odd.lst"]
		tmpcls_m=[l.replace('.','_m1.') for l in tmpcls]
		run("e2classextract.py {clsfile} --refinemulti --setname {tmpcls}".format(clsfile=output_cls["even"][s],tmpcls=tmpcls[0]))
		run("e2classextract.py {clsfile} --refinemulti --setname {tmpcls}".format(clsfile=output_cls["odd"][s],tmpcls=tmpcls[1]))
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

def do_compare(data):
	
	e=EMData(data["ptclfile"],data["idx"])
	e.transform(data["xform"])
	ret=[]
	for pj in data["proj"]:
		ret.append(e.cmp("ccc",pj))
	return ret
	

def run(cmd):
	print cmd
	launch_childprocess(cmd)
	
if __name__ == '__main__':
	main()
	