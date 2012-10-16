#!/usr/bin/env python
###	procpdb.py	Steven Ludtke	2/2002

#N procpdb.py
#F This program is for performing simple operations on PDB files
#T LM
#1
#P <infile>	Input file
#P [<outfile>]	Output file
#P [rot=<alt,az,phi>]	Rotation in EMAN convention
#P [trans=<dx,dy,dz>]	Translate (before rotation)
#P [posttrans=<dx,dy,dz>]	Post-translation
#P [centeratoms]	Center based on the mean atom position
#P [centerelec]	Center based on the center of electron charge
#P [centermass]	Center based on the center of mass
#P [animorph=<n>,<apix>,<vecfile>]	This will use a morph vector file (segment3d) to morph atom positions
#P [apix=<A/pix>]	multiply trans by apix
#P [split]	Split file at TER records. Results imperfect.
#P [include=[helix,sheet,other]]
#P [chains=<chainltr>]	eg - "ABO", for extracting a portion of a complex
#D Simple program for manipulating PDB filess in the EMAN convention

from EMAN2 import *
import os
import sys
import random
import time
import string
from os import system
from os import unlink
from math import *
#from sys import argv


#		   1		 2		 3		 4		 5		 6		 7
# 01234567890123456789012345678901234567890123456789012345678901234567890123456789
# ATOM	  1  N   ASP L   1	  57.169  62.936   7.756  1.00 30.65	  1ACY 129
# HEADER	COMPLEX(ANTIBODY/HIV-1 FRAGMENT)		10-FEB-94   1ACY	  1ACY   2

# Note the unusual entries 'A' for ambiguous, treated as Nitrogen and 'w' for Water
atomdefs={'H':(1.0,1.00794),'C':(6.0,12.0107),'A':(7.0,14.00674),'N':(7.0,14.00674),'O':(8.0,15.9994),'P':(15.0,30.973761),
	'S':(16.0,32.066),'W':(18.0,1.00794*2.0+15.9994),'AU':(79.0,196.96655) }

#### MAIN############################

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage:\nprocpdb.py <input>\nprocpdb.py <input> <output> [rot=<alt,az,phi>] [trans=<dx,dy,dz>] [centeratoms] [centerelec] [centermass] [apix=<A/pixel>]\n."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	####################
	parser.add_argument("--animorph", "-AN", type=str, help="This will use a morph vector file (segment3d) to morph atom positions,#P [animorph=<n>,<apix>,<vecfile>]",default=None)
	parser.add_argument("--apix", "-A", type=float, help="apix", default=1.0)
	parser.add_argument("--scale", "-S", type=float, help="scale", default=1.0)

	parser.add_argument("--center", "-C", type=str, help="center of the rotation, (0,0,0)", default='0.0,0.0,0.0')
	
	parser.add_argument("--chains",type=str,help="String list of chain identifiers to include, eg 'ABEFG'", default=None)
	parser.add_argument("--posttrans", "-PR", type=str, help="posttransform, (0,0,0)",default='0,0,0')
	parser.add_argument("--qrot", type=str, help="spin,sigrot rotation, (angle,x,y,z)",default='0,0,0,0')
	parser.add_argument("--split", "-spl",type=float, help="split file at TER records. Results imperfect.", default=0)
	parser.add_argument("--trans", "-TR", type=str, help="transform, (0,0,0)",default='0,0,0')
	parser.add_argument("--transcenter", "-TC", type=int, help="transform center,1 atom, 2 electron dentisty 3 mass", default=0)
	parser.add_argument("--include", type=str,help="savetype", default=["helix","sheet","other"])
	parser.add_argument("--quiet",action="store_true",default=False,help="Verbose is the default")
	parser.add_argument("--mirror",type=str, help="mirror",default='False')
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--type", "-T", type=str, help="convention type", default='eman')
	parser.add_argument("--rot", "-R", type=str, help="rotation, (0,0,0)",default='0,0,0')
	#eman input, default setting
	parser.add_argument("--az", "-az", type=float, help="az in eman convention.", default=0)
	parser.add_argument("--alt", "-alt", type=float, help="alt in eman convention.", default=0)
	parser.add_argument("--phi", "-phi", type=float, help="phi.", default=0)
	#imagic
	parser.add_argument("--alpha", "-alpha", type=float, help="alpha in imagic convention.", default=0)
	parser.add_argument("--beta", "-beta", type=float, help="beta in imagic convention.", default=0)
	parser.add_argument("--gamma", "-gamma", type=float, help="gamma in imagic convention.", default=0)
#spider
	parser.add_argument("--theta", "-theta", type=float, help="theta.", default=0)
	parser.add_argument("--psi", "-psi", type=float, help="psi in spider convention.", default=0)
#xyz
	parser.add_argument("--xtilt", "-xtilt", type=float, help="xtilt in xyz convention.", default=0)
	parser.add_argument("--ytilt", "-ytilt", type=float, help="ytilt in xyz convention.", default=0)
	parser.add_argument("--ztilt", "-ztilt", type=float, help="ztilt in xyz convention.", default=0)
#mrc
	parser.add_argument("--omega", "-omega", type=float, help="omega.", default=0)
#quaternion
	parser.add_argument("--e0", "-e0", type=float, help="e0 in quaternion convention.", default=0)
	parser.add_argument("--e1", "-e1", type=float, help="e1 in quaternion convention.", default=0)
	parser.add_argument("--e2", "-e2", type=float, help="e2 in quaternion convention.", default=0)
	parser.add_argument("--e3", "-e3", type=float, help="e3 in quaternion convention.", default=0)
#spin
	parser.add_argument("--n1", "-n1", type=float, help="n1.", default=0)
	parser.add_argument("--n2", "-n2", type=float, help="n2.", default=0)
	parser.add_argument("--n3", "-n3", type=float, help="n3.", default=0)
#sigrot
	parser.add_argument("--q", "-q", type=float, help="q in sgirot convention.", default=0)
#matrix
	parser.add_argument("--matrix", "-matrix", type=str, help="transform matrix.", default='0,0,0,0,0,0,0,0,0,0,0,0')
	
	(options, args) = parser.parse_args()
	
	savetypes=options.include
	apix=options.apix
	chains=options.chains
	ptrans=options.posttrans
	qrot=options.qrot
	rot=options.rot
	spl=options.split
	trans=options.trans
	center=options.transcenter
	cen=options.center
	animorph=options.animorph
	chains=options.chains
	type=options.type
	mirror=options.mirror
	scale=options.scale
	
	az=options.az
	alt=options.alt
	phi=options.phi
	alpha=options.alpha
	beta=options.beta
	gamma=options.gamma
	theta=options.theta
	psi=options.psi
	omega=options.omega
	e0=options.e0
	e1=options.e1
	e2=options.e2
	e3=options.e3
	n1=options.n1
	n2=options.n2
	n3=options.n3
	q=options.q
	
	xtilt=options.xtilt
	ytilt=options.ytilt
	ztilt=options.ztilt
	
	matrix=options.matrix
	
	print "rot="
	print rot
	
	s=rot.split(',')
	rot=(float(s[0]),float(s[1]),float(s[2]))
	
	s=matrix.split(',')
	mat=(float(s[0]),float(s[1]),float(s[2]),float(s[3]),float(s[4]),float(s[5]),float(s[6]),float(s[7]),float(s[8]),float(s[9]),float(s[10]),float(s[11]))
	s=qrot.split(',')
	qrot=(float(s[0])/57.29578,float(s[1]),float(s[2]),float(s[3]))
	s=cen.split(',')
	cen=(float(s[0]),float(s[1]),float(s[2]))
	s=ptrans.split(',')
	ptrans=(float(s[0]),float(s[1]),float(s[2]))
	s=trans.split(',')
	trans=(float(s[0]),float(s[1]),float(s[2]))
	
	if (mirror=="True")or(mirror=="true"): mirror=1
	else:mirror=0 
	
	if (type=='eman'):
		t=Transform()
		t.set_params({"type":type,"az":az,"alt":alt,"phi":phi,"scale":scale,"mirror":mirror,"tx":trans[0],"ty":trans[1],"tz":trans[2]})
		print type
		v1=(1.6,-85.45,44.62)
		#v2=(1.65,-84.3,45.57)
		#v3=(0.33,-84.25,46.32)
		v1_transformed=t.transform(v1)
		#v2_transformed=t.transform(v2)
		#v3_transformed=t.transform(v3)
		print v1_transformed 
	
	elif (type=='imagic'):
		t=Transform()
		t.set_params({"type":type,"alpha":alpha,"beta":beta,"gamma":gamma,"scale":scale,"mirror":mirror,"tx":trans[0],"ty":trans[1],"tz":trans[2]})
		
	#spider input
	elif (type=='spider'):
		t=Transform()
		t.set_params({"type":type,"phi":phi,"theta":theta,"psi":psi,"scale":scale,"mirror":mirror,"tx":trans[0],"ty":trans[1],"tz":trans[2]})
	#xyz input
	elif (type=='xyz'):
		
		t=Transform()
		t.set_params({"type":type,"xtilt":xtilt,"ytilt":ytilt,"ztilt":ztilt,"scale":scale,"mirror":mirror,"tx":trans[0],"ty":trans[1],"tz":trans[2]})
	#mrc input
	elif (type=='mrc'):	
		t=Transform()
		t.set_params({"type":type,"phi":phi,"theta":theta,"omega":omega,"scale":scale,"mirror":mirror,"tx":trans[0],"ty":trans[1],"tz":trans[2]})
	#quaternion input
	elif (type=='quaternion'):
		t=Transform()
		t.set_params({"type":type,"e0":e0,"e1":e1,"e2":e2,"e3":e3,"scale":scale,"mirror":mirror,"tx":trans[0],"ty":trans[1],"tz":trans[2]})
	#spin
	elif (type=='spin'):
		t=Transform()
		t.set_params({"type":type,"Omega":omega,"n1":n1,"n2":n2,"n3":n3,"scale":scale,"mirror":mirror,"tx":trans[0],"ty":trans[1],"tz":trans[2]})
		
	#sgirot
	elif (type=='sgirot'):
		t=Transform()
		t.set_params({"type":type,"q":q,"n1":n1,"n2":n2,"n3":n3,"scale":scale,"mirror":mirror,"tx":trans[0],"ty":trans[1],"tz":trans[2]})
		#matrix
	elif (type=='matrix'):
		t=Transform()
		t=Transform([mat[0],mat[1],mat[2],mat[3],mat[4],mat[5],mat[6],mat[7],mat[8],mat[9],mat[10],mat[11]])
		
	else: 
		print"get error, please input the right convention, example eman, imagic, spider, mrc, xyz, spin, sgirot, quaternion"
		sys.exit(1)

######################################
	if len(args)<2 : 
		parser.error("Input and output files required") 
		sys.exit(1)

	inp = open(args[0], 'r')

	lines = inp.readlines()
	inp.close()

	helixlist = []
	sheetlist = []
	title = ""
	atoms = 0
	mass = 0
	electrons = 0
# center of atoms,center of electron density,center of mass,min,max
	xval = [0, 0, 0, None, None]
	yval = [0, 0, 0, None, None]
	zval = [0, 0, 0, None, None]

	if (spl) :
		for i in range(len(lines)):
			if (lines[i][:4] == "ATOM" or lines[i][:6] == "HETATM") :
				start = i
				break

		sn = 0
		out = None
		for i in lines[start:]:
			if i[:3] == "ter" or i[:3] == "TER" or sn == 0 :
				if (out) : out.close()
				sn += 1
				out = open("%s.%d.ent" % (args[1], sn), "w")
				for j in lines[:start]: out.write(j)
			else : out.write(i)

		sys.exit(1)

	for i in lines:
		if (i[:6] == 'HEADER') : title = i[10:72].strip()
		if (i[:5] == 'HELIX') :
			helixlist.extend(range(int(i[21:25]), int(i[33:37]) + 1))
		if (i[:5] == 'SHEET') :
			sheetlist.extend(range(int(i[22:26]), int(i[33:37]) + 1))
		if ((i[:4] == 'ATOM' or i[:6] == 'HETATM') and (not chains or (i[21] in chains))) :
			a = i[12:14].strip()
			res = int(i[22:26].strip())
			if (res in helixlist) :
				 if (not "helix" in savetypes) : continue
			elif (res in sheetlist) :
				if (not "sheet" in savetypes) : continue
			elif (not "other" in savetypes) : continue
	
			x = float(i[30:38])
			y = float(i[38:46])
			z = float(i[46:54])
			atm = (0, 0)
			try:
				atm = atomdefs[a]
				mass += atm[1]
				electrons += atm[0]
				atoms += 1.0
			except:
				print "Unknown Atom '%s' ignored." % a
	
			xval[0] += x
			yval[0] += y
			zval[0] += z
			xval[1] += x * atm[0]
			yval[1] += y * atm[0]
			zval[1] += z * atm[0]
			xval[2] += x * atm[1]
			yval[2] += y * atm[1]
			zval[2] += z * atm[1]
			if (xval[3] == None) :
				xval[3] = x
				yval[3] = y
				zval[3] = z
			xval[3] = min(xval[3], x)
			yval[3] = min(yval[3], y)
			zval[3] = min(zval[3], z)
			xval[4] = max(xval[4], x)
			yval[4] = max(yval[4], y)
			zval[4] = max(zval[4], z)
	
	xval[0] /= atoms
	yval[0] /= atoms
	zval[0] /= atoms
	xval[1] /= electrons
	yval[1] /= electrons
	zval[1] /= electrons
	xval[2] /= mass
	yval[2] /= mass
	zval[2] /= mass
	
	print title
	print "%1.0f atoms   %1.0f electrons   mass= %1.3f kDa" % (atoms, electrons, mass / 1000.0)
	print "atom center = (%1.2f,%1.2f,%1.2f)" % (xval[0], yval[0], zval[0])
	print "electron density center = (%1.2f,%1.2f,%1.2f)" % (xval[1], yval[1], zval[1])
	print "center of mass = (%1.2f,%1.2f,%1.2f)" % (xval[2], yval[2], zval[2])
	print "x range: %1.1f to %1.1f  (%1.1f)" % (xval[3], xval[4], xval[4] - xval[3])
	print "y range: %1.1f to %1.1f  (%1.1f)" % (yval[3], yval[4], yval[4] - yval[3])
	print "z range: %1.1f to %1.1f  (%1.1f)" % (zval[3], zval[4], zval[4] - zval[3])
	
	#		   1		 2		 3		 4		 5		 6		 7
	# 01234567890123456789012345678901234567890123456789012345678901234567890123456789
	# ATOM	  1  N   ASP L   1	  57.169  62.936   7.756  1.00 30.65	  1ACY 129
	# HEADER	COMPLEX(ANTIBODY/HIV-1 FRAGMENT)		10-FEB-94   1ACY	  1ACY   2
	
	if (center == 1) : trans = (-xval[0], -yval[0], -zval[0])
	elif (center == 2) : trans = (-xval[1], -yval[1], -zval[1])
	elif (center == 3) : trans = (-xval[2], -yval[2], -zval[2])
		
	#print "\n len(args)=..................."
	#print len(args)
	
	if (len(args) > 1) :	
		if (animorph) :
			v = open(animorph[2], "r")
			vecs = v.readlines()
			v.close()
			
			apix = float(animorph[1])
			vecs = [map(lambda y:float(y), x.split()) for x in vecs]
			for i, v in enumerate(vecs):
				vecs[i] = [apix * (v[3]), apix * (v[4]), apix * (v[5]), apix * (v[0] - v[3]), apix * (v[1] - v[4]), apix * (v[2] - v[5])]
	#			vecs[i]=[apix*(v[0]+v[3])/2,apix*(v[1]+v[4])/2,apix*(v[2]+v[5])/2,apix*(v[3]-v[0]),apix*(v[4]-v[1]),apix*(v[5]-v[2])]
			
			for m in range(int(animorph[0])):
				aw = m / (float(animorph[0]) - 1.0)
				out = open("%s.%02d%s" % (args[1][:-4], m, args[1][-4:]), "w")
				for i in lines:
					if (i[:6] == 'HEADER') : title = i[10:72].strip()
					if (i[:4] == 'ATOM' or i[:6] == 'HETATM') :
						if (not chains or (i[21] in chains)):
							res = int(i[22:26].strip())
							if (res in helixlist) :
								if (not "helix" in savetypes) : continue
							elif (res in sheetlist) :
								if (not "sheet" in savetypes) : continue
							elif (not "other" in savetypes) : continue
				
							x = float(i[30:38])
							y = float(i[38:46])
							z = float(i[46:54])
							
							dx, dy, dz, nm = 0, 0, 0, 0
	#						print "-------------------"
							for v in vecs:
								w = 1.0 / ((x - v[0]) ** 4 + (y - v[1]) ** 4 + (z - v[2]) ** 4)
	#							print w
								dx += w * v[3]
								dy += w * v[4]
								dz += w * v[5]
								nm += w
							dx /= nm
							dy /= nm
							dz /= nm
							out.write(i[:30] + " %7.2f %7.2f %7.2f" % (x + dx * aw, y + dy * aw, z + dz * aw) + i[54:])
					else : out.write(i)
				out.close()
				print "%d/%d complete" % (m, int(animorph[0]))
					
		else:
			out=open(args[1],"w")
			print args[1]
			
			for i in lines:
				if (i[:6]=='HEADER') : title=i[10:72].strip()
				if (i[:4]=='ATOM' or i[:6]=='HETATM') :
					if (not chains or (i[21] in chains)):
						res=int(i[22:26].strip())
						if (res in helixlist) :
							if (not "helix" in savetypes) : continue
						elif (res in sheetlist) :
							if (not "sheet" in savetypes) : continue
						elif (not "other" in savetypes) : continue
						x=float(i[30:38])+(apix*trans[0])
						y=float(i[38:46])+(apix*trans[1])
						z=float(i[46:54])+(apix*trans[2])
						v1=(x,y,z)
						v1_transformed=t.transform(v1)
						x4=v1_transformed[0]
						y4=v1_transformed[1]
						z4=v1_transformed[2]
						out.write(i[:30]+" %7.2f %7.2f %7.2f"%(x4,y4,z4)+i[54:])
						
				else : out.write(i)
	   
			out.close()
			
if __name__ == "__main__":
	main()