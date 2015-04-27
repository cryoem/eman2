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

	
def main():
	
	########################
	### Classify by averaging
	
	usage="e2classifytree.py <projection> <particle> [options]"
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
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
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
	
	### Build tree
	if not os.path.isfile(options.nodes):
		buildtree(projs,options.threads,options.nodes)
		
	E2progress(E2n,0.6)
	print "Start classification..."
	## Classify particles
	t={}
	clsmx=[EMData(1,npt) for i in range(7)]
	for i in range(options.threads):
		ai=[x for x in range(npt) if x%options.threads==i]
		t[i]=threading.Thread(target=classify,args=(ptcl,ai,options.nodes,clsmx,options.align,options.aligncmp,options.cmp,options.ralign,options.raligncmp))
		t[i].start()
	for i in range(options.threads):
		t[i].join()
		
	if os.path.isfile(options.output):
		os.system('rm '+options.output)
	for  i in clsmx:
		i.write_image(options.output,-1)
		
	print "Finished~"
	E2progress(E2n,1.0)
	E2end(E2n)
	

	
def buildtree(projs,thread,nodes):
	simxorder={0:"tx",1:"ty",2:"alpha",3:"mirror",4:"scale"}
	tmpsim="tmp_simmix.hdf"
	
	### Building the similarity matrix for all projections
	cmd="e2simmx.py {pj} {pj} {smx} --align=rotate_translate_flip --aligncmp=sqeuclidean:normto=1 --cmp=sqeuclidean --saveali --force --parallel=thread:{thr:d}".format(pj=projs,smx=tmpsim,thr=thread)
	os.system(cmd)
	
	### Initialize buttom level nodes
	if (os.path.isfile(nodes)):
		os.system("rm {}".format(nodes))
		
	npj=EMUtil.get_image_count(projs)
	for i in range(npj):
		pj=EMData(projs,i)
		pj.process_inplace("normalize.edgemean")
		pj["tree"]=[-1,-1]
		pj.write_image(nodes,i)
		
	simmx=EMData(tmpsim,0)
	dst=EMNumPy.em2numpy(simmx)
	epms=[EMData(tmpsim,i+1) for i in range(5)]
	pms=[EMNumPy.em2numpy(i) for i in epms]
	ai =range(dst[0].size)		# index of each node in "nodes.hdf"
	big=(dst.max()+1)


	dst+=np.identity(dst[0].size)*big
	om=0
	
	for k in range(dst[0].size-1):

		### Finish one level, move to next
		if (dst[0].size-om < 2):
			#omat.fill(0)
			#print ar
			om=0
			if os.path.isfile("rr.lst"): os.system("rm rr.lst")
			rr=LSXFile("rr.lst")
			for r in range(len(ai)):
				for t in range(2): rr.write(r,ai[r],nodes)
				
			os.system("e2simmx.py rr.lst rr.lst {} --align=rotate_translate_flip --aligncmp=sqeuclidean:normto=1 --cmp=sqeuclidean --saveali --force --parallel=thread:{:d}".format(tmpsim,thread))
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
		print "Averaging ",ai[x],ai[y]," to ",npj+k
		
		alipm=[a[x,y] for a in pms]
		alidict={"type":"2d"}
		for i,a in simxorder.items():
			alidict[a]=float(alipm[i])
		alidict["mirror"]=int(alidict["mirror"])
		#print x,y,alipm,alidict
		tr=Transform(alidict)
		img1=EMData(nodes,ai[x])
		img2=EMData(nodes,ai[y])
		img1.process_inplace("xform",{"transform":tr})
		img1.add(img2)
		img1.div(2)
		img1.process_inplace("normalize.edgemean")
		img1["tree"]=[ai[x],ai[y]]
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
	
	os.system("rm {}".format(tmpsim))
	
	return 	


def compare(ref,target,align=None,alicmp=("dot",{}),cmp=("dot",{}), ralign=None, alircmp=("dot",{})):
	
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

		scr=(target.cmp(cmp[0],ta,cmp[1]),p["tx"],p["ty"],p["alpha"],p["mirror"],p["scale"])
#			ta.write_image("dbug.hdf",-1)

#				print ta["source_n"],target["source_n"]
			#psub=target.process("math.sub.optimal",{"ref":ta})
			#nout=ta["source_n"]*3
			#ta.write_image("dbug_%d.hdf"%target["source_n"],nout)
			#target.write_image("dbug_%d.hdf"%target["source_n"],nout+1)
			#psub.write_image("dbug_%d.hdf"%target["source_n"],nout+2)


	else :
		scr=(target.cmp(cmp[0],ref,cmp[1]),0,0,0,1.0,False)

	return scr, p
	
	
	##print pj,n
	#align=parsemodopt("rotate_translate_flip")
	#alicmp=parsemodopt(acmp)
	#ccmp=parsemodopt(acmp)
	#ref.del_attr("xform.align2d")
	#ta=ref.align(align[0],target,align[1],alicmp[0],alicmp[1])
	#ralign["xform.align2d"] = ta.get_attr("xform.align2d")
	#ref.del_attr("xform.align2d")
	#ta = ref.align(ralign,target,ralign[1],alircmp[0],alircmp[1])
	
	
	#t=apt.get_attr("xform.align2d")
	#pm=t.get_params('2d')
	##apt=apt.align("refinecg",b,{},alncmp[0],alncmp[1])
	##apt.process_inplace("normalize.edgemean")
	#scr=b.cmp(ccmp[0],apt,ccmp[1])
	#return scr,pm
	
	
### Classify each particle using the tree
def classify(ptcl,ai,nodes,clsmx,align,alicmp,cmp,ralign,alircmp):
	
	rt=EMUtil.get_image_count(nodes)-1
	rimg=EMData()
	rimg.read_image(nodes,rt,True)
	nimg=EMData()
	rl=rimg["tree"][0]
	rr=rimg["tree"][1]
	for pp in ai:
		ni=rt
		nl=rl
		nr=rr
		prob=EMData(ptcl,pp)
		while(1):
			nimg.read_image(nodes,ni,True)
			if (nimg["tree"][0]<0):
				break
			limg=EMData(nodes,nimg["tree"][0])
			rimg=EMData(nodes,nimg["tree"][1])
			dl=compare(prob,limg,align,alicmp,cmp,ralign,alircmp)
			dr=compare(prob,rimg,align,alicmp,cmp,ralign,alircmp)
			if dl<dr:
				ni=nimg["tree"][0]
			else:
				ni=nimg["tree"][1]
		#print 
		#print pp,ni
		nimg=EMData(nodes,ni)
		d,pm=compare(prob,nimg,align,alicmp,cmp,ralign,alircmp)
		#print pm
		clsmx[0].set_value_at(0,pp,ni)
		clsmx[1].set_value_at(0,pp,1)
		clsmx[2].set_value_at(0,pp,pm['tx'])
		clsmx[3].set_value_at(0,pp,pm['ty'])
		clsmx[4].set_value_at(0,pp,pm['alpha'])
		clsmx[5].set_value_at(0,pp,pm['mirror'])
		clsmx[6].set_value_at(0,pp,pm['scale'])
	
	
	
	
if __name__ == '__main__':
	main()
	 
