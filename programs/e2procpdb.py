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
from math import *
import os
import sys

atomdefs={'H':(1.0,1.00794),'C':(6.0,12.0107),'A':(7.0,14.00674),'N':(7.0,14.00674),'O':(8.0,15.9994),'P':(15.0,30.973761),
	'S':(16.0,32.066),'W':(18.0,1.00794*2.0+15.9994),'AU':(79.0,196.96655) }

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
	parser.add_argument("--trans", "-TR", type=str, help="transform, (0,0,0)",default='0,0,0')
	parser.add_argument("--include", type=str,help="savetype", default=["helix","sheet","other"])
	parser.add_argument("--mirror",type=str, help="mirror",default='False')
	parser.add_argument("--type", "-T", type=str, help="convention type", default='eman')
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
	
	matrix=options.matrix
	trans=options.trans
	mirror=options.mirror
	
	s=matrix.split(',')
	mat=(float(s[0]),float(s[1]),float(s[2]),float(s[3]),float(s[4]),float(s[5]),float(s[6]),float(s[7]),float(s[8]),float(s[9]),float(s[10]),float(s[11]))
	s=trans.split(',')
	trans=(float(s[0]),float(s[1]),float(s[2]))
	if (mirror=="True")or(mirror=="true"): mirror=1
	else:mirror=0 

	if len(args)<2 : 
		parser.error("Input and output files required")
		sys.eixt(1)

	print args[0]
	print args[0]
	if (options.type=='eman'):
		t=Transform()
		t.set_params({"type":options.type,"az":options.az,"alt":options.alt,"phi":options.phi,"scale":options.scale,"mirror":mirror,"tx":trans[0],"ty":trans[1],"tz":trans[2]})
		print type
		v1=(1.6,-85.45,44.62)
		v1_transformed=t.transform(v1)
		print v1_transformed 
	elif (options.type=='imagic'):
		t=Transform()
		t.set_params({"type":options.type,"alpha":options.alpha,"beta":options.beta,"gamma":options.gamma,"scale":options.scale,"mirror":mirror,"tx":trans[0],"ty":trans[1],"tz":trans[2]})
	elif (options.type=='spider'):
		t=Transform()
		t.set_params({"type":options.type,"phi":options.phi,"theta":options.theta,"psi":options.psi,"scale":options.scale,"mirror":mirror,"tx":trans[0],"ty":trans[1],"tz":trans[2]})
	elif (options.type=='xyz'):
		t=Transform()
		t.set_params({"type":options.type,"xtilt":options.xtilt,"ytilt":options.ytilt,"ztilt":options.ztilt,"scale":options.scale,"mirror":mirror,"tx":trans[0],"ty":trans[1],"tz":trans[2]})
	elif (options.type=='mrc'):	
		t=Transform()
		t.set_params({"type":options.type,"phi":options.phi,"theta":options.theta,"omega":options.omega,"scale":options.scale,"mirror":mirror,"tx":trans[0],"ty":trans[1],"tz":trans[2]})
	elif (options.type=='quaternion'):
		t=Transform()
		t.set_params({"type":options.type,"e0":options.e0,"e1":options.e1,"e2":options.e2,"e3":options.e3,"scale":options.scale,"mirror":mirror,"tx":trans[0],"ty":trans[1],"tz":trans[2]})
	elif (options.type=='spin'):
		t=Transform()
		t.set_params({"type":options.type,"Omega":options.omega,"n1":options.n1,"n2":options.n2,"n3":options.n3,"scale":options.scale,"mirror":mirror,"tx":trans[0],"ty":trans[1],"tz":trans[2]})
	elif (options.type=='sgirot'):
		t=Transform()
		t.set_params({"type":options.type,"q":options.q,"n1":options.n1,"n2":options.n2,"n3":options.n3,"scale":options.scale,"mirror":mirror,"tx":trans[0],"ty":trans[1],"tz":trans[2]})
	elif (options.type=='matrix'):
		t=Transform()
		t=Transform([mat[0],mat[1],mat[2],mat[3],mat[4],mat[5],mat[6],mat[7],mat[8],mat[9],mat[10],mat[11]])
	else: 
		print"get error, please input the right convention, example eman, imagic, spider, mrc, xyz, spin, sgirot, quaternion"
		sys.exit(1)

	inp = open(args[0], 'r')
	lines = inp.readlines()
	inp.close()
	
	outputlines=pdb_transform(t,lines,options.center,options.include,options.animorph,options.apix,options.chains,trans)
	print "output="
	print len(outputlines)

	out=open(args[1],"w")
	for i in outputlines: out.write(i)
	out.close()


def pdb_transform(t,lines,center=0,savetypes=["helix","sheet","other"],animorph=None, apix=1.0,chains=None,trans=(0,0,0)):
	
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
		
		
	if (center == 1) : trans = (-xval[0], -yval[0], -zval[0])
	elif (center == 2) : trans = (-xval[1], -yval[1], -zval[1])
	elif (center == 3) : trans = (-xval[2], -yval[2], -zval[2])
	
	if (animorph) :
		print "animorph"
		s=animorph.split(',')
		animorph=(float(s[0]),float(s[1]),str(s[2]))
		
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
			counter=0
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
						
						x4=x+dx*aw
						y4=y+dy*aw
						z4=z+dz*aw
						i=i[:30]+" %7.2f %7.2f %7.2f"%(x4,y4,z4)+i[54:]
						lines[counter]=i
				counter=counter+1
			
			print "%d/%d complete" % (m, int(animorph[0]))
		return lines

	else: 
		print "transform"
		print t
		print len(lines)
		counter=0
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
					v1_transformed=t*v1
					x4=v1_transformed[0]
					y4=v1_transformed[1]
					z4=v1_transformed[2]
					
					i=i[:30]+" %7.2f %7.2f %7.2f"%(x4,y4,z4)+i[54:]
					lines[counter]=i
			counter = counter+1
		return lines
						
if __name__ == "__main__":
	main()
