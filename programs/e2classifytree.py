#!/usr/bin/env python

#
# Author: Muyuan Chen, April 2015
# Copyright (c) 2000-2007 Baylor College of Medicine
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

# e2classifytree.py  April 2015 Muyuan Chen
# This program classify particles using a binary tree


from EMAN2 import *
import os
import numpy as np
from EMAN2jsondb import JSTask
	
def main():
	
	usage="""e2classifytree.py <projection> <particle> [options]
	
	Classify particles using a binary tree. Can be used as an alternative for e2simmx2stage.py + e2classify.py.
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--threads", type=int,help="", default=12)
	parser.add_argument("--nodes", type=str,help="", default="nodes.hdf")
	#parser.add_argument("--clsmx", type=str,help="", default="clsmx.hdf")
	parser.add_argument("--output", type=str,help="", default="clsmx.hdf")
	parser.add_argument("--align",type=str,help="The name of an 'aligner' to use prior to comparing the images", default=None)
	parser.add_argument("--aligncmp",type=str,help="Name of the aligner along with its construction arguments",default="dot")
	parser.add_argument("--ralign",type=str,help="The name and parameters of the second stage aligner which refines the results of the first alignment", default=None)
	parser.add_argument("--raligncmp",type=str,help="The name and parameters of the comparitor used by the second stage aligner. Default is dot.",default="dot")
	parser.add_argument("--cmp",type=str,help="The name of a 'cmp' to be used in comparing the aligned images", default="dot:normalize=1")
	parser.add_argument("--cmpdiff", action="store_true", default=False ,help="Compare using the difference of the two children")
	parser.add_argument("--incomplete", type=int,help="The degree of incomplete allowed in the tree on each level", default=0)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--parallel", default=None, help="parallelism argument")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()
	E2n=E2init(sys.argv,options.ppid)
	
	options.align=parsemodopt(options.align)
	options.aligncmp=parsemodopt(options.aligncmp)
	options.ralign=parsemodopt(options.ralign)
	options.raligncmp=parsemodopt(options.raligncmp)
	options.cmp=parsemodopt(options.cmp)
	
	projs=args[0]
	#projsimmx=args[1]
	ptcl=args[1]
	npj=EMUtil.get_image_count(projs)
	npt=EMUtil.get_image_count(ptcl)
	if options.parallel==None:
		par="thread:{:d}".format(options.threads)
	else:
		par=options.parallel
		
	### Build tree
	### always overwrite the tree here now
	#if not os.path.isfile(options.nodes):
	print "Building binary tree..."
	buildtree(projs,par,options.nodes,options.incomplete,options.verbose)
	#else:
		#print "Using existing tree..."
	
	## Generate children pairs for comparison
	print "Generating children pairs for comparison..."
	if options.cmpdiff:
		nodepath= os.path.dirname(options.nodes)
		masktmp='/'.join([nodepath,"tmp_msk.hdf"])
		if os.path.isfile(masktmp): os.remove(masktmp)
		cmptmp='/'.join([nodepath,"tmp_cmp.hdf"])
		if os.path.isfile(cmptmp):
			os.remove(cmptmp)
		makechildpair(options.nodes, cmptmp, masktmp)
	else:
		masktmp=None
		cmptmp=None
	
	E2progress(E2n,0.5)
	#exit()
	print "Starting classification..."
	### Classify particles
	
		
	clsmx=[EMData(1,npt) for i in range(7)]
	nnod=EMUtil.get_image_count(options.nodes)
	if options.parallel :
		from EMAN2PAR import EMTaskCustomer
		etc=EMTaskCustomer(options.parallel)
		tasks=[]
		step=50
		tt=[range(i,i+step) for i in range(0,npt-step,step)]
		tt.append(range(tt[-1][-1]+1,npt))
		
		for it in tt:
			tasks.append(TreeClassifyTask(ptcl, it, options.nodes, options.align, options.aligncmp, options.cmp, options.ralign, options.raligncmp, cmptmp, masktmp))
		
		taskids=etc.send_tasks(tasks)
		ptclpernode=[0 for i in range(nnod)]
		nfinished=0
		while len(taskids)>0 :
			haveprogress=False
			time.sleep(3)
			curstat=etc.check_task(taskids)
			for i,j in enumerate(curstat):
				if j==100 :
					haveprogress=True
					rslt=etc.get_results(taskids[i])
					rslt= rslt[1]
					for r in rslt:
						nfinished+=1
						if options.verbose>0: print "Particle:",r["id"],"\tnodes:",r["choice"]
						for c in r["choice"]:
							ptclpernode[c]+=1
						clsmx[0].set_value_at(0,r["id"],r["cls"])
						for nt in range(1,7):
							clsmx[nt].set_value_at(0,r["id"],r["simmx"][nt])
			
			taskids=[j for i,j in enumerate(taskids) if curstat[i]!=100]
			if haveprogress: print "{:d}/{:d} finished".format(nfinished,npt)
			E2progress(E2n, 0.5 + float(nfinished)/npt)
			
		for i in range(nnod):
			ndtmp=EMData(options.nodes,i,True)
			ndtmp["tree_nptls"]=ptclpernode[i]
			ndtmp.write_image(options.nodes,i)
	
	else:
		
		### To record the number of particles in each branch of the tree
		for i in range(nnod):
			ndtmp=EMData(options.nodes,i,True)
			ndtmp["tree_nptls"]=0
			ndtmp.write_image(options.nodes,i)
		t={}
		clsmx=[EMData(1,npt) for i in range(7)]
		for i in range(options.threads):
			ai=[x for x in range(npt) if x%options.threads==i]
			t[i]=threading.Thread(target=classify,args=(ptcl,ai,options.nodes,clsmx,options.align,options.aligncmp,options.cmp,options.ralign,options.raligncmp,cmptmp,masktmp))
			t[i].start()
		for i in range(options.threads):
			t[i].join()
		
	if os.path.isfile(options.output):
		os.remove(options.output)
	for  i in clsmx:
		i.write_image(options.output,-1)
	
	if options.cmpdiff:	
		os.remove(cmptmp)
		os.remove(masktmp)
	print "Finished~"
	E2progress(E2n,1.0)
	E2end(E2n)
	

	
def buildtree(projs,par,nodes,incomplete,verbose):
	simxorder={0:"tx",1:"ty",2:"alpha",3:"mirror",4:"scale"}
	nodepath= os.path.dirname(nodes)
	tmpsim='/'.join([nodepath,"tmp_simmix.hdf"])
	par="--parallel "+par
	### Building the similarity matrix for all projections
	#cmd="e2simmx.py {pj} {pj} {smx} --align=rotate_translate_flip --aligncmp=sqeuclidean:normto=1 --cmp=sqeuclidean --saveali -v {vb:d} --force {parallel}".format(pj=projs,smx=tmpsim, parallel=par, vb=verbose-1)
	cmd="e2simmx.py {pj} {pj} {smx} --align=rotate_translate_flip --aligncmp=sqeuclidean:normto=1 --cmp=frc:maxres=10.0 --ralign=refine --raligncmp=frc:maxres=10.0 --saveali -v {vb:d} --force {parallel}".format(pj=projs,smx=tmpsim, parallel=par, vb=verbose-1)
	print cmd
	launch_childprocess(cmd)
	
	### Initialize buttom level nodes
	if (os.path.isfile(nodes)):
		os.remove(nodes)
	
	tr=Transform()
	npj=EMUtil.get_image_count(projs)
	for i in range(npj):
		pj=EMData(projs,i)
		pj.process_inplace("normalize.edgemean")
		pj["tree_children"]=[-1,-1]
		pj["tree_transform"]=tr
		pj.write_image(nodes,i)
		
	simmx=EMData(tmpsim,0)
	dst=EMNumPy.em2numpy(simmx)
	epms=[EMData(tmpsim,i+1) for i in range(5)]
	pms=[EMNumPy.em2numpy(i) for i in epms]
	ai =range(dst[0].size)		# index of each node in "nodes.hdf"
	big=(dst.max()+1)
	tmplist="tmplist.lst"

	dst+=np.identity(dst[0].size)*big
	om=0
	
	for k in range(dst[0].size-1):
		#print len(ai),dst[0].size-om 

		### Finish one level, move to next
		if (dst[0].size-om < 2+incomplete):
			#omat.fill(0)
			#print ar
			om=0
			if os.path.isfile(tmplist): os.remove(tmplist)
			rr=LSXFile(tmplist)
			
			## FIXME I do not know why I need to write twice here, but e2simmx reads one less line if I don't
			for r,a in enumerate(ai):
				rr.write(r,a,nodes)
				rr.write(r,a,nodes)
			
			if len(ai)<10:
				par=""
			cmd="e2simmx.py {lst} {lst} {sim} --align=rotate_translate_flip --aligncmp=sqeuclidean:normto=1  --ralign=refine --raligncmp=frc:maxres=10.0 --cmp=frc:maxres=10.0 --saveali -v {vb:d} --force {parallel}".format(lst=tmplist, sim=tmpsim, parallel=par, vb=verbose-1)
			launch_childprocess(cmd)
			#launch_childprocess("e2simmx.py {lst} {lst} {sim} --align=rotate_translate_flip --aligncmp=sqeuclidean:normto=1 --cmp=sqeuclidean --saveali -v {vb:d} --force {parallel}".format(lst=tmplist, sim=tmpsim, parallel=par, vb=verbose-1))
			simmx=EMData(tmpsim,0)
			dst=EMNumPy.em2numpy(simmx)
			dst+=np.identity(dst[0].size)*big
			epms=[EMData(tmpsim,i+1) for i in range(5)]
			#print 'size',epms[0]['nx'],epms[0]['ny']
			pms=[EMNumPy.em2numpy(i) for i in epms]
		
		### Find the two children
		x,y = np.where(dst == np.min(dst))
		x=x[0]
		y=y[0]
		
		### Do averaging
		if verbose>0: print "Averaging ",ai[x],ai[y]," to ",npj+k
		
		alipm=[a[x,y] for a in pms]
		alidict={"type":"2d"}
		for i,a in simxorder.items():
			alidict[a]=float(alipm[i])
		alidict["mirror"]=int(alidict["mirror"])
		#print x,y,alipm,alidict
		tr=Transform(alidict)
		img1=EMData(nodes,ai[x])
		img2=EMData(nodes,ai[y])
		img1["tree_transform"]=tr
		img1.write_image(nodes,ai[x])
		img1.process_inplace("xform",{"transform":tr})
		img1.add(img2)
		img1.div(2)
		#img1.process_inplace("normalize.edgemean")
		img1["tree_children"]=[ai[x],ai[y]]
		img1.write_image(nodes,-1)
		om+=1
		
		
		### Update the distance matrix
		for nrw in range(5):
			pms[nrw][x,:]=0
			pms[nrw][:,x]=0
			pms[nrw]=np.delete(pms[nrw],y,0)
			pms[nrw]=np.delete(pms[nrw],y,1)
			
		dst[x,:]=big
		dst[:,x]=big
		dst=np.delete(dst,y,0)
		dst=np.delete(dst,y,1)
		
		ai[x]=npj+k
		del ai[y]
	
	os.remove(tmpsim)
	os.remove(tmplist)
	
	return 	

### Do ref-target comparison
def compare(ref,target,options):
	
	align=options["align"]
	alicmp=options["alicmp"]
	cmp=options["cmp"]
	ralign=options["ralign"]
	alircmp=options["alircmp"]
	
	if align[0] :
		ref.del_attr("xform.align2d")
		ta=ref.align(align[0],target,align[1],alicmp[0],alicmp[1])
		#if verbose>3: print ta.get_attr("xform.align2d")
		#ta.debug_print_params()

		if ralign and ralign[0]:
		
			ralign[1]["xform.align2d"] = ta.get_attr("xform.align2d")
			ref.del_attr("xform.align2d")
			ta = ref.align(ralign[0],target,ralign[1],alircmp[0],alircmp[1])

		t =  ta.get_attr("xform.align2d")
		t.invert()
		p = t.get_params("2d")

		scr=(target.cmp(cmp[0],ta,cmp[1]),1,p["tx"],p["ty"],p["alpha"],p["mirror"],p["scale"])
#			

	else :
		scr=(target.cmp(cmp[0],ref,cmp[1]),1,0,0,0,1,1)

	return scr
	
### Compare the two children of current node
def cmpchild(ref,nimg,limg,rimg,mask,options):
	
	align=options["align"]
	alicmp=options["alicmp"]
	ralign=options["ralign"]
	alircmp=options["alircmp"]
	
	### Align to current node first
	if align[0] :
		ref.del_attr("xform.align2d")
		ta=ref.align(align[0],nimg,align[1],alicmp[0],alicmp[1])
		#if verbose>3: print ta.get_attr("xform.align2d")
		#ta.debug_print_params()

		if ralign and ralign[0]:
		
			ralign[1]["xform.align2d"] = ta.get_attr("xform.align2d")
			ref.del_attr("xform.align2d")
			ta = ref.align(ralign[0],nimg,ralign[1],alircmp[0],alircmp[1])
	
	### Compare the aligned particle with the two children 
	

	#tmp="tmp.hdf"
	#ta.write_image(tmp,-1)
	ta.mult(mask)
	#ta.write_image(tmp,-1)
	#limg.write_image(tmp,-1)
	#rimg.write_image(tmp,-1)
	
	cmp=parsemodopt("dot:normalize=1")
	dl=ta.cmp(cmp[0],limg,cmp[1])
	dr=ta.cmp(cmp[0],rimg,cmp[1])
	
	return dl,dr
	
	
### Classify each particle using the tree
def classify(ptcl,ai,nodes,clsmx,align,alicmp,cmp,ralign,alircmp,cmptmp,masktmp):
	#tmp="tmp.hdf"
	#if os.path.isfile(tmp): os.remove(tmp)
	#ai=[1]
	options={"align":align, "alicmp":alicmp, "cmp":cmp, "ralign":ralign, "alircmp":alircmp}
	rt=EMUtil.get_image_count(nodes)-1
	rimg=EMData()
	rimg.read_image(nodes,rt,True)
	nimg=EMData()
	rl=rimg["tree_children"][0]
	rr=rimg["tree_children"][1]
	for pp in ai:
		ni=rt
		nl=rl
		nr=rr
		prob=EMData(ptcl,pp)
		#prob.process_inplace("normalize")
		choice=[]
		while(1):
			choice.append(ni)
			nimg.read_image(nodes,ni,True)
			nimg["tree_nptls"]+=1
			nimg.write_image(nodes,ni)
			if (nimg["tree_children"][0]<0):
				break
			if cmptmp==None:
				limg=EMData(nodes,nimg["tree_children"][0])
				rimg=EMData(nodes,nimg["tree_children"][1])
				dl=compare(prob,limg,options)
				dr=compare(prob,rimg,options)
				dl=dl[0]
				dr=dr[0]
			else:
				nimg.read_image(nodes,ni)
				limg=EMData(cmptmp,nimg["tree_children"][0])
				rimg=EMData(cmptmp,nimg["tree_children"][1])
				mask=EMData(masktmp,ni)
				dl,dr=cmpchild(prob,nimg,limg,rimg,mask,options)
				
			if dl<dr:
				ni=nimg["tree_children"][0]
			else:
				ni=nimg["tree_children"][1]
		print "Particle",pp,"nodes",choice
		nimg=EMData(nodes,ni)
		pm=compare(prob,nimg,options)
		clsmx[0].set_value_at(0,pp,ni)
		for nt in range(1,7):
			clsmx[nt].set_value_at(0,pp,pm[nt])
	
	
def makechildpair(nodes, cmptmp, masktmp):
	for ni in range(EMUtil.get_image_count(nodes)):
		nimg=EMData(nodes,ni)
		if (nimg["tree_children"][0]<0):
			nimg.to_zero()
			nimg.write_image(masktmp,ni)
			continue
		limg=EMData(nodes,nimg["tree_children"][0])
		rimg=EMData(nodes,nimg["tree_children"][1])
		#print ni,nimg["tree_children"][0],nimg["tree_children"][1]
		tr=limg["tree_transform"]
		limg.process_inplace("xform",{"transform":tr})
		msk=nimg.copy()
		msk.process_inplace("threshold.belowtominval",{"minval":msk["sigma"],"newval":-1})
		msk.process_inplace("threshold.binary")
		nimg.sub(limg)
		nimg.process_inplace("math.absvalue")
		nimg.process_inplace("threshold.belowtozero",{"minval":nimg["sigma"]})
		nimg.mult(msk)
		nimg.write_image(masktmp,ni)
		 
		limg.mult(nimg)
		rimg.mult(nimg)
		limg.write_image(cmptmp,nimg["tree_children"][0])
		rimg.write_image(cmptmp,nimg["tree_children"][1])
		
		
class TreeClassifyTask(JSTask):
	
	def __init__(self,ptcl,ptid,nodes,align=None,alicmp=("dot",{}),cmp=("dot",{}), ralign=None, alircmp=("dot",{}),cmptmp=None,masktmp=None):
		rt=EMUtil.get_image_count(nodes)
		if cmptmp==None or masktmp==None:
			### Compare to the two children seperately 
			data={"images":["cache",ptcl,ptid], "nodes":["cache",nodes,0,rt]}
			cmpdiff=False
		else:
			### Mask out the difference between the two children
			cn=EMUtil.get_image_count(cmptmp)
			mn=EMUtil.get_image_count(masktmp)
			data={"images":["cache",ptcl,ptid], "nodes":["cache",nodes,0,rt], "cmptmp":["cache",cmptmp,0,cn], "masktmp":["cache",masktmp,0,mn] }
			cmpdiff=True
			
		JSTask.__init__(self,"TreeClassify",data,{},"")
		self.options={"align":align, "alicmp":alicmp, "cmp":cmp, "ralign":ralign, "alircmp":alircmp,"cmpdiff":cmpdiff, "id":ptid}
	
	
	def execute(self,callback=None):
		options=self.options
		rst=[]
		for tt in self.data["images"][2]:
			
			ni=self.data["nodes"][3]-1
			nimg=EMData()
			prob=EMData(self.data["images"][1],tt)
			
			choice=[]
			while(1):
				### Compare through the tree
				choice.append(ni)
				nimg.read_image(self.data["nodes"][1],ni,True)
				#nimg["tree_nptls"]+=1
				#nimg.write_image(self.data["nodes"][1],ni)
				
				if (nimg["tree_children"][0]<0):
					### At the leaf of the tree, stop.
					break
				
				if options["cmpdiff"]:
					nimg.read_image(self.data["nodes"][1],ni)
					mask=EMData(self.data["masktmp"][1],ni)
					limg=EMData(self.data["cmptmp"][1],nimg["tree_children"][0])
					rimg=EMData(self.data["cmptmp"][1],nimg["tree_children"][1])
					dl,dr=cmpchild(prob,nimg,limg,rimg,mask,options)
				else:
					limg=EMData(self.data["nodes"][1],nimg["tree_children"][0])
					rimg=EMData(self.data["nodes"][1],nimg["tree_children"][1])
					dl=compare(prob,limg,options)
					dr=compare(prob,rimg,options)
					dl=dl[0]
					dr=dr[0]
					
				if dl<dr:
					ni=nimg["tree_children"][0]
				else:
					ni=nimg["tree_children"][1]
			
			#print pp,choice
			nimg=EMData(self.data["nodes"][1],ni)
			pm=compare(prob,nimg,options)
		
			rst.append({"cls":ni,"simmx":pm, "id":tt,"choice":choice})
		return rst
		#clsmx[0].set_value_at(0,pp,ni)
		#for nt in range(1,7):
			#clsmx[nt].set_value_at(0,pp,pm[nt])
		
		
		
		
		

if __name__ == '__main__':
	main()
	 
