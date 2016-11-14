#!/usr/bin/env python
# Muyuan Chen 12/2014
# Minor Modification to the pathwalking structure


import EMAN2
from EMAN2 import *
import numpy as np
import math
import random
from sys import stdout

def distance(point1, point2):
    return math.sqrt(sum([(x[0]-x[1])**2 for x in zip(point1, point2)]))

def distance_lines(p1,p2,q1,q2):
	a=p2-p1
	b=q2-q1
	c=q1-p1
	crossab=np.cross(a,b)
	dist=abs(np.dot(c,crossab))/np.linalg.norm(crossab)
	return dist

def read_pdb(filename):
	
	points = {}
	atomnumber=np.array([])
	pdbfile = open(filename, "r")
	lines = pdbfile.readlines()
	pdbfile.close()

	chainnum=1
	for line in (i for i in lines if i.startswith("ATOM  ")):
		lchain=line[21].strip() or None
		break
	
	
	count = 0
	for line in (i for i in lines if i.startswith("ATOM  ")):
		#latomtype = line[12:15].strip()
		nchain = line[21].strip() or None			
		if nchain != lchain:
			chainnum+=1
			pos=[-1,-1,-1]
			points[count]=(pos)
			count+=1
			atomnumber=np.append(atomnumber,0)
		
		atomnumber=np.append(atomnumber,int(line[22:27]))
		pos = [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())]
		points[count] = (pos)
		count += 1
		lchain=nchain
	
	print "Read ",chainnum," chains, ",count-chainnum+1," atoms."
	return points,count,atomnumber
    
def write_pdb(filename,ncent,SX,SY,SZ,color,atomnumber,vx):
	shp=ncent.shape
	if color==None:
		color=np.zeros(shp[0],float)
	out = open(filename,"w")
	nchn=97
	chain = chr(nchn)
	count=0
	for atom in range(shp[0]):
		if (ncent[atom,0]<0):
			out.write("""TER  %6d      ALA %s%4d\n"""%(count, chain, atom))
			nchn+=1
			chain=chr(nchn)
			continue
		out.write(
			"ATOM %6d  CA  ALA %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f     S_00  0\n"
			%(atomnumber[atom], chain,atomnumber[atom] ,(ncent[atom,0]-SX/2)*vx[0], (ncent[atom,1]-SY/2)*vx[1], (ncent[atom,2]-SZ/2)*vx[2], 1, color[atom]) 
			)
		count+=1
	out.write("""TER  %6d      ALA %s%4d\n"""%(shp[0], chain, atom))
	out.write("END")
	out.close()

#def read_helix(filename):
    #helixfile = open(filename,"r")
    #lines = helixfile.readlines()
    #hlxid=set()
    #for ln in lines:
		#hlxid.add(int(ln[0:10].strip()))
		#hlxid.add(int(ln[10:15].strip()))
    #return hlxid

def find_neighbors(ncent):
	naa=ncent.shape[0]
	neighbors=np.empty((naa,3),int)
	for i in range(1,naa-1):
		dist=np.empty(naa,float)
		for j in range(naa):
			dist[j]=abs(distance(ncent[i,:],ncent[j,:]))
		dista=dist.argsort()
		neighbors[i,:]=[i,dista[1],dista[2]]
	return neighbors
	
def get_density(mrc,posi,posb,posf):
	midp=10
	point1=posi
	pointo=[posb,posf]
	mrcp=0
	for point2 in pointo:
		a=[(x[1]-x[0])/midp for x in zip(point1, point2)]
		mpt=0
		p=point1
		for i in range(midp/2):
			np=[(x[0]+i*x[1]) for x in zip(p,a)]
			mpt+=mrc[int(np[0]),int(np[1]),int(np[2])]
			p=np
		mrcp+=mpt
	if mrcp==0:
		mrcp=10000
		#print posi,posb,posf
	return mrcp/midp/2


def read_helix(filename):
	print "Reading helix atoms..."
	points = {}
	count=0
	pdbfile = open(filename, "r")
	lines = pdbfile.readlines()
	pdbfile.close()
	nhlx=0
	hlxid=set()
	for line in (i for i in lines if i.startswith("HELIX  ")):
		atomid=[int(line[21:27].strip()), int(line[33:38].strip())]
		for atline in (i for i in lines if i.startswith("ATOM  ")):
			cn=int(atline[22:30].strip())
			if cn>=atomid[0] and cn<=atomid[1]:
				if atline[13:15]=="CA":
					hlxid.add(cn)
					pos = [float(atline[30:38].strip()), float(atline[38:46].strip()), float(atline[46:54].strip()),nhlx]
					points[count] = pos
					count += 1
		nhlx+=1
	return hlxid

		
		
def calc_score(ncent,i,mrc,wd=1,wm=8,wa=2,ww=0,bl=0):
	
	posi=ncent[i]
	posb=ncent[i-1]
	posf=ncent[i+1]
	
	da=[(x[0]-x[1]) for x in zip(posi, posb)]
	db=[(x[0]-x[1]) for x in zip(posi, posf)]
	ang=np.dot(da,db)/(distance(posb,posi)*distance(posf,posi))	# angle
	mrcp=get_density(mrc,posi,posb,posf)# density at atom
	dista=abs(distance(posb,posi)-bl)
	distb=abs(distance(posi,posf)-bl)
	distc=abs(distance(posb,posf)-3.78)
	
	have_alter=1
	for p in [posb,posf]:
		nsum=0
		for c in ncent:
			if distance(c,p) < 5:
				nsum+=1
		if nsum<5:
			have_alter=0
	
	if have_alter==0: wd=wd/2
	dis=dista+distb
	score=wd*dis+(wm/mrcp)+wa*ang+ww*distc
	#print atomnumber[i],wd*dis,wm/mrcp,wa*ang,nsum
	#
	
	if (posi[0]<0) or (posb[0]<0) or (posf[0]<0):
		score=0
	return score

def main():

	usage = """ Usage...
	"""
	parser = EMAN2.EMArgumentParser(usage=usage,version=EMAN2.EMANVERSION)
	parser.add_argument("--output", type=str,help="Output pdb file")
	parser.add_argument("--mapin", type=str,help="mapin file for input",default=None)
	parser.add_argument("--pdbin", type=str,help="pdb file for input",default=None)
	parser.add_argument("--edgefile", type=str ,help="Fixed edges",default=None)
	parser.add_argument("--thr", type=float,help="Threshold for removing Ca atoms", default=-1)
	parser.add_argument("--dochange", type=int,help="Just mark the sidechain atoms or remove them", default=0)
	parser.add_argument("--mode", type=str,help="Mode: e/s/d/c")


	
	(options, args) = parser.parse_args()
	eval_atom=0
	shaking=0
	decision=0
	crossbond=0
	if options.mode.startswith("e"):
		eval_atom=1
	elif options.mode.startswith("s"):
		shaking=1
	elif options.mode.startswith("d"):
		decision=1
	elif options.mode.startswith("c"):
		crossbond=1
	
	
	points,count,atomnumber=read_pdb(options.pdbin)
	if options.edgefile<>None:
		hlxid=read_helix(options.edgefile)
	else:
		hlxid=set()
	bondl=0#3.78
	naa=count
	dochange=1
	
	
	mrc=EMData(options.mapin)
	SX=mrc.get_xsize()
	SY=mrc.get_ysize()
	SZ=mrc.get_zsize()
	apix_x=mrc["apix_x"]
	apix_y=mrc["apix_y"]
	apix_z=mrc["apix_z"]
	mrcdata=np.empty((SX,SY,SZ),float)
	for i in range(SX):
		for j in range(SY):
			for k in range(SZ):
				mrcdata[i,j,k]=mrc.get_value_at(i,j,k)
				
	ncent=np.empty((naa,3),float)
	
	for i in range(naa):
		ncent[i,:]=points[i]/np.array([apix_x,apix_y,apix_z])#+np.array([SX/2,SY/2,SZ/2])
	thresh=.3
	neark=5
	noiset=1000
	gaussig=1
	atomcolor=np.zeros(naa,float)
	#neighbors=find_neighbors(ncent)
	
	
	if crossbond==1:
		csbd=np.empty(naa,int)
		for i in range(1,naa-1):
			pos1=(ncent[i]+ncent[i+1])/2
			mind=10
			for j in range(i+2,naa-1):
				pos2=(ncent[j]+ncent[j+1])/2
				dp=distance(pos1,pos2)
				dl=distance_lines(ncent[i],ncent[i+1],ncent[j],ncent[j+1])
				if(dp<3):
					dd=dl
				else:
					dd=dp
				if(dd<mind):
					mind=dd
					csbd[i]=j
			atomcolor[i]=mind
			if (options.thr<0):
				options.thr=np.mean(atomcolor)+2*np.std(atomcolor)
				print options.thr
			if (atomcolor[i]<options.thr):
				print i, csbd[i], atomcolor[i]
				atomcolor[i]=options.thr
				atomcolor[i+1]=options.thr
				atomcolor[csbd[i]]=options.thr
				atomcolor[csbd[i]+1]=options.thr
				for j in range(i+1,(1+csbd[i]+i+1)/2):
					k=csbd[i]+i+1-j
					
					tmp=np.copy(ncent[j])
					ncent[j]=ncent[k]
					ncent[k]=tmp
					#print j,ncent[j],k,ncent[k]
				
				
				
	
	if decision==1:
		# delete side chains ( sharp angle, two neighbors are close to each other)
		dellist=[]
		orina=naa
		for i in range(1,naa-1):
			ang=calc_score(ncent,i,mrcdata,0,0,1)
			wid=calc_score(ncent,i,mrcdata,0,0,0,1)
			inhlx=(atomnumber[i+1] in hlxid) and (atomnumber[i-1] in hlxid) and (atomnumber[i-1] in hlxid)
			if ang>cos(options.thr) and wid<4 and inhlx==False:
				atomcolor[i]=10
				dellist.append(i)
		
		print np.hstack((atomnumber[dellist],atomcolor[dellist]))
		ncent=np.delete(ncent,dellist,axis=0)
		atomnumber=np.delete(atomnumber,dellist)
		atomcolor=np.delete(atomcolor,dellist)
		
		naa= len(ncent)
		# long bonds
		acent=np.copy(ncent)
		k=1
		nn=max(atomnumber)
		distb=np.empty(naa,float)
		for i in range(1,naa-2):
			posi=ncent[i]
			#print i+1
			posf=ncent[i+1]
			distb[i]=abs(distance(posi,posf))
		srtdst=np.sort(distb)
		thr=srtdst[naa*2-orina-1]
		for i in range(1,naa-2):
			posi=ncent[i]
			posf=ncent[i+1]
			if (distb[i]>thr and posi[0]>0 and posf[0]>0):
				pos=(posi+posf)/2
				print i,k,distb[i],get_density(mrc,pos,posi,posf)
				acent=np.insert(acent,i+k,pos,axis=0)
				atomnumber=np.insert(atomnumber,i+k,nn+1)
				atomcolor=np.insert(atomcolor,i+k,100)
				nn+=1
				k+=1
				
		ncent=np.copy(acent)
		print len(ncent)
		#for i in range(1,naa-1):
			#posi=ncent[i]
			#posb=ncent[i-1]
			#posf=ncent[i+1]
			#dista=abs(distance(posb,posi))
			#distb=abs(distance(posi,posf))
			#atomcolor[i]=max(dista,distb)
	
	
	if shaking==1:
		mrc=EMData(options.mrcin)
		atoms=PointArray()
		atoms.read_from_pdb(options.pdbin)
		atoms.sim_set_pot_parms(3.78,10, 5, 0, 0.0,10,mrc,2,0)


		for iters in range(500):
			print iters
			atoms.sim_printstat()
			atoms.sim_minstep_seq(.01)
			#if iters%100==0:
				#atoms.save_to_pdb("a"+str(iters)+".pdb") 
		atoms.save_to_pdb(options.output) 
		exit()
				
	
	
	if eval_atom==1:
		#if options.dochange:
			#dellist=[]
			#for i in range(1,naa-1):
				#ang=calc_score(ncent,i,mrcdata,0,0,1)
				#wid=calc_score(ncent,i,mrcdata,0,0,0,1)
				#inhlx=(atomnumber[i+1] in hlxid) and (atomnumber[i-1] in hlxid) and (atomnumber[i-1] in hlxid)
				#if ang>cos(.6) and wid<4 and inhlx==False:
					#atomcolor[i]=10
					#dellist.append(i)
			
			#print np.hstack((atomnumber[dellist],atomcolor[dellist]))
			#ncent=np.delete(ncent,dellist,axis=0)
			#atomnumber=np.delete(atomnumber,dellist)
			#atomcolor=np.delete(atomcolor,dellist)
			#naa= len(ncent)
			
		
		
		
		thresh=options.thr
		while 1:
			maxscore=0
			for i in range(1,naa-1):			
				score=calc_score(ncent,i,mrcdata,3,2,1,-3)
				if (atomnumber[i] in hlxid):
					score=0
				atomcolor[i]=score
				if score>maxscore:
					maxscore=score
					dellist=i
					
			
			if thresh==-1:
				thresh=np.mean(atomcolor)+np.std(atomcolor)*2
				
			if (options.dochange==0):
				break
			maxscore=max(atomcolor)
			if maxscore< thresh:
				break
			
			
			print np.hstack((atomnumber[dellist],atomcolor[dellist]))
			ncent=np.delete(ncent,dellist,axis=0)
			atomnumber=np.delete(atomnumber,dellist)
			atomcolor=np.delete(atomcolor,dellist)
			#print atomnumber[dochange],maxscore
			#ncent=np.delete(ncent,dochange,axis=0)
			naa-=1
			
	#atomcolor-=min(atomcolor)
	#atomcolor=atomcolor/max(atomcolor)*99.9
	print np.mean(atomcolor)+np.std(atomcolor)*2
	write_pdb(options.output,ncent,0,0,0,atomcolor,atomnumber,[apix_x,apix_y,apix_z])

    
    
if __name__ == '__main__':
    main()
	