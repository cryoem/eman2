#!/usr/bin/env python

# Author: Ian Rees (ian.rees@bcm.edu), 08/01/2010
# Copyright (c) 2000-2011 Baylor College of Medicine

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


import EMAN2
import collections
import math
import random
import sys
import optparse
import os
import commands
import sys
import operator
import copy
import re
import numpy
import json




def cross_product(a,b):
    """Cross product of two 3-d vectors. from http://www-hep.colorado.edu/~fperez/python/python-c/weave_examples.html"""
    cross = [0]*3
    cross[0] = a[1]*b[2]-a[2]*b[1]
    cross[1] = a[2]*b[0]-a[0]*b[2]
    cross[2] = a[0]*b[1]-a[1]*b[0]
    return numpy.array(cross)



def norm_vector(a,b):
	v1=(a[0]-b[0],a[1]-b[1],a[2]-b[2])
	lengthV=math.sqrt(v1[0]**2+v1[1]**2+v1[2]**2)
	vec=(v1[0]/lengthV,v1[1]/lengthV,v1[2]/lengthV)
	return vec





class PathWalker(object):

	def __init__(self, filename=None, start=None, end=None, edgefile=None, edges=None, dmin=2.0, dmax=4.2, average=3.78, atomtype='CA', chain=None, noise=0, solver=False):

		# Run parameters
		self.dmin = dmin
		self.dmax = dmax
		self.average = average	
		self.filename = filename
		self.basename = self.get_basename()
		self.chain = chain
		self.atomtype = atomtype
		self.noise = noise
		self.solver = solver

		if self.atomtype in ['all', 'None', '']:
			self.atomtype = None
			
		# Read PDB file
		self.points = self.read_pdb(atomtype=self.atomtype, chain=self.chain)
		
		# Point graph
		self.itree = {}

		# Calculated distances
		self.distances = {}
		self.weighted = {}

		if self.average > self.dmax:
			self.dmax = self.average
		if self.average < self.dmin:
			self.dmin = self.average

		for i in self.points.keys():
			self.itree[i] = set()
				
		for count1, point1 in self.points.items():
			for count2, point2 in self.points.items():
				d, w = self.calcweight(point1, point2)
				self.distances[(count1, count2)] = d
				self.weighted[(count1, count2)] = w
				if self.cutoff(d):
					self.itree[count1].add(count2)
				# print count1, count2, self.distances[(count1, count2)], self.weighted[(count1, count2)]


		# Read an edge fragment file... 1 string of points per line, separated space
		self.fixededges = self.read_fixed(edgefile)
		
		# ... add any additional edges, and/or start+end
		if edges:
			self.fixededges.extend(edges)

		self.start = min(self.points)
		self.end = max(self.points)
		if start != None:
			self.start = start
		if end != None:
			self.end = end

		print "Note: linking start/end: ", self.start, self.end
		self.fixededges.append((self.start, self.end))
	
		# Process forced edges
		for link in self.fixededges:
			self.itree[link[0]].add(link[1])
			self.itree[link[1]].add(link[0])
			self.weighted[(link[0], link[1])] = 0
			self.weighted[(link[1], link[0])] = 0

		# Some useful statistics
		self.endpoints = self.find_endpoints()
		self.branches = self.find_branches()
		self.total_branches = len(self.branches)
		self.total_endpoints = len(self.endpoints)
		self.total_nodes = len(self.itree)
		self.total_edges = len(self.distances)




	def run(self):
		# Write the processed PDB file
		self.stats()
		
		self.write_pdb(self.points.keys(), tree=self.itree)

		if not self.solver:
			return

		ret = self.get_json()
		#print ret
		
		# Solve paths
		ret['paths'] = self.solve(solver=self.solver)

		# We need to rotate the path so the start atom is first..
		paths = ret.get('paths', {})
		for pathid in sorted(paths.keys()):
			print "--- Evaluating path %s of %s ---"%(pathid, len(paths))
			d = paths[pathid]

			# Rotate the path so the starting element is first in the list
			rot = d['path'].index(self.start or d['path'][0])
			path = collections.deque(d['path'])
			path.rotate(-rot)
			d['path'] = list(path)
		
			# Evaluate path and calculate CA ramachandran angles
			d.update(self.evaluate_path(d['path']))	
			d.update(self.ca_ram(d['path']))


		# Write output
		self.write_pdbs(paths, filename=self.basename+".out.pdb")			
		self.write_json(ret)
		
		

	
	def read_fixed(self, edgefile):
		# Edge file format:
		# 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37
		# 39 40 41 42 43
		# ...
		fixededges = []
		if not edgefile:
			return fixededges
		f = open(edgefile)
		fragments = [i.strip() for i in f.readlines()]
		f.close()
		for fragment in fragments:
			fragment = map(int, fragment.split())
			for i in range(len(fragment)-1):
				fixededges.append((fragment[i], fragment[i+1]))
		return fixededges
	
	
	
	def get_json(self):
		keys = [
			'total_nodes', 
			'total_edges', 
			'total_branches', 
			'total_endpoints', 
			'edges_fixed', 
			'noise', 
			'filename', 
			'basename', 
			'dmax', 
			'dmin', 
			'average',
			'chain',
			'atomtype',
			'solver',
			'start',
			'end',
			'fixededges',
			'endpoints',
			'branches']
			
		ret = {}
		for key in keys:
			i = getattr(self, key, None)
			if i != None:
				ret[key] = i

		return ret
		


	def stats(self):
		print "\n=== Statistics for %s ==="%self.filename
		print "Total Nodes: %s"%len(self.itree)
		print "Total Edges: %s"%len(self.distances)
		
		print "Endpoints: %s"%self.total_endpoints
		print "\t", self.endpoints
		
		print "Branches: %s"%self.total_branches
		print "\t", self.branches
		
		print "Edges:"
		for k in sorted(self.points.keys()):
			print "\t%-10s%s"%(k, self.itree[k])




	def solve(self, solver="exhaustive"):
		solver = solver.lower()
		method = getattr(self, '_solve_%s'%solver, None)
		if method:
			return method()
		# else:
		#	raise ValueError, "Unknown solver: %s"%solver
		return {}




	# def children(start, tree, minlength=1):
	def _solve_exhaustive(self, minlength=1):

		minlength = len(self.itree) - 1
		print "\n=== Exhaustive Search ==="

		endpaths = set()
		to_crawl = collections.deque()
		to_crawl.append((self.start,))
	
		while to_crawl:
			current = to_crawl.pop()
			children = self.itree.get(current[-1], ())
			# print endpaths

			if not children - set(current) and len(current) > minlength:
				print "Found path of len", len(current)
				endpaths.add(current)
				continue
		
			for child in children:
				if child not in current:
					to_crawl.append(current + (child,))
	
		ret = {}
		for count, path in enumerate(endpaths):
			ret[count] = {'path':path, 'solver':'exhaustive', 'score':0}
	
		return ret
		


	def _solve_concorde(self):
		print "\n=== Solving TSP with Concorde ==="
		
		tspfile = self.basename+".tsp"
		outfile = self.basename+".out"
		
		print "TSPLib file: %s"%tspfile
		print "Results file: %s"%outfile
		
		self.write_tsplib(filename=tspfile)
		os.system("concorde -x -m -o %s %s"%(outfile, tspfile))

		return self._readtour_concorde(outfile)




	def _solve_lkh(self):
		print "\n=== Solving TSP with LKH ==="
		
		tspfile = self.basename+".tsp"
		lkhfile = self.basename+".lkh"
		outfile = self.basename+".out"

		lkh = """PROBLEM_FILE = %s\nOUTPUT_TOUR_FILE = %s\nPRECISION = 100\n"""%(tspfile,outfile)
		f = open(lkhfile, "w")
		f.write(lkh)
		f.close()

		self.write_tsplib(filename=tspfile)
		os.system("LKH %s"%(lkhfile))
		
		return self._readtour_lkh(outfile)
		
		
		
		
	def _readtour_lkh(self, tourfile):		
		# TSP Tour file format..
		f = open(tourfile)
		r = [i.strip() for i in f.readlines()]
		f.close()

		score_tsp = filter(lambda i:i.startswith('COMMENT : Length'), r)
		if score_tsp:
			score_tsp = score_tsp[0].partition('=')[2].strip()
			score_tsp = float(score_tsp)

		path = map(int, r[r.index("TOUR_SECTION")+1:r.index("-1")])
		
		# We need to convert the TSP #'s back to atom #s
		keyorder = sorted(self.points.keys())				
		path = [keyorder[i-1] for i in path]
		
		# Save output
		ret = {
			'path': path,
			'score': score_tsp,
			'solver':'lkh'
		}
		return {0:ret}




	def _readtour_concorde(self, tourfile):
		# Concorde Tour file format:
		# numpoints
		# point1, point2, point3...
		f = open(tourfile)
		r = f.readlines()
		f.close()
		
		path = [map(int, i.split()) for i in r][1:]
		path = reduce(operator.concat, r)
		path = [i+1 for i in r]

		# Save output
		ret = {
			'path': path,
			'score': 0,
			'solver':'concorde'
		}
		return {0:ret}

		


	def evaluate_path(self, path=None):
		
		print "\n=== Evaluating Path ==="
		
		def fmt(*args):
			r = ['%-10s'%i for i in args]
			return "".join(r)
		
		print fmt("PDB Edge", "", "Distance", "Weight")
		breaks = []

		pmin = min(path)
		pmax = max(path)
		distances = []
		r = []
		for i in range(len(path)-1):
			p1 = path[i]
			p2 = path[i+1]
			d1 = '%0.2f'%self.distances.get((p1, p2))
			d2 = self.weighted.get((p1, p2))
			distances.append(self.distances.get((p1, p2)))
				
				
			if (path[i] in [pmin, pmax] or path[i] in [path[i+1]+1, path[i+1]-1]):
				r.append(True)
				b = ""
			else:
				r.append(False)
				breaks.append((d1, d2))
				b = "*"
				
			print fmt(p1, p2, d1, d2, b)

		r2 = []
		j = 0
		for i in range(len(r)-1):
			if r[i]==True and r[i+1]==True: j += 1
			else: j = 0
			r2.append(j)
		
		rmax = max(r2)
	 	score_fragment = float(rmax) / float(len(r2))

		correctbonds = len(path)-len(breaks)
		score_path = float(correctbonds) / len(path)

		print "\nPath quality statistics:"
		print "\tNoise: %s"%self.noise
		print "\tLongest continuous fragment score: %s"%score_fragment
		print "\tNumber of breaks: %s out of %s bonds"%(len(breaks), len(path))
		print "\tPath score: %0.5f"%(score_path)
		print "\tMinimum distance:", min(distances)
		print "\tMaximum distance:", max(distances)
		
		# Save output
		ret = {
			'path_length': len(path),
			'path_breaks': len(breaks),
			'score_fragment': score_fragment,
			'score_path': score_path,
			'dmin': min(distances),
			'dmax': max(distances)
		}
		return ret
		
		
			

	def cutoff(self, d):
		return self.dmin <= d <= self.dmax
	



	def calcweight(self, point1, point2):
		d = self.distance(point1, point2)
		if self.cutoff(d):
			w = int(math.fabs(d-self.average)*100)
			# w = 0
		else:
			w = int((math.fabs(d-self.average)*100)**2)
		if w > 100000:
			w = 100000
		return d, w




	def distance(self, point1, point2):
		return math.sqrt(sum([(x[0]-x[1])**2 for x in zip(point1, point2)]))




	def find_endpoints(self):
		"""Return list of points that have != 2 connections as possible termini; note: circular paths may have 2 links..."""
		return [k for k,v in self.itree.items() if len(self.itree.get(k)) == 1]




	def find_branches(self):
		"""Return points with > 2 connections"""
		return [k for k,v in self.itree.items() if len(self.itree.get(k)) > 2]
	
	
	

	def get_basename(self):
		basename = os.path.basename(self.filename).split(".")
		basename = ".".join(basename[:basename.index("pdb")])
		r = "%s-(\d+).tsp"%basename
		r = re.compile(r)
		runs = map(int, [(r.findall(i) or [0])[0] for i in os.listdir(".")])
		if not runs:
			runcount = 0
		else:
			runcount = max(runs)+1

		return "%s-%s"%(basename, runcount)
		


	def write_json(self, soln, filename=None):
		outfile = filename or self.basename + ".json"
		print "Writing solution dictionary to %s"%outfile
		f = open(outfile, 'w')
		json.dump(soln, f, indent=1)
		f.close()

					
					

	def write_pdb(self, path, filename=None, tree=None):
		d = {0:{'path':path}}
		return self.write_pdbs(d, filename=filename, tree=tree)
		
		
		

	def write_pdbs(self, paths, filename=None, tree=None):
		outfile = filename or self.basename + ".pdb"
		out = open(outfile,"w")
		chain = "A"

		print "\n=== Writing %s ==="%outfile

		for pathid in sorted(paths.keys()):
			d = paths[pathid]
			path = d['path']
			count = 1
			out.write(
				"MODEL     %4d\n"%pathid
			)
			for atom in path:
				point = self.points.get(atom)
				out.write(
					"ATOM %6d  C   ALA %s%4d    %8.3f%8.3f%8.3f %5.2f     0     S_00  0\n"
					%(atom, chain, atom, point[0], point[1], point[2], len(self.itree.get(atom,[])))
					)
				count += 1
			out.write("""TER  %6d      ALA %s%4d\n"""%(count, chain, atom))

			if tree:
				connected = []
				count = 0
				for k,v in tree.items():
					for v2 in v:
						if (k,v2) in connected or (v2,k) in connected:
							continue
						count+=1
						out.write('CONECT %4d %4d\n'%(k,v2))
						connected.append((k,v2))
				print "Wrote %s edges"%count
		
			# out.write('CONECT %4d %4d\n'%(atoms[0], atoms[-1]))
			out.write("ENDMDL\n")


		print "Wrote %s points"%len(path)


		out.write("END")
		out.close()




	def makenoise(self):
		return random.gauss(self.average, self.noise) - self.average



	def read_pdb(self, atomtype=None, chain=None):
		print "\n=== Reading %s (atom type %s, chain %s) ==="%(self.filename, atomtype, chain)
		points = {}

		pdbfile = open(self.filename, "r")
		lines = pdbfile.readlines()
		pdbfile.close()

		atomtypes = set()
		chains = set()
		for line in (i for i in lines if i.startswith("ATOM  ")):
			atomtypes.add(line[12:15].strip() or None)
			chains.add(line[21].strip() or None)
			
		print "Found chains %s and atomtypes %s"%(chains, atomtypes)
		if len(atomtypes) == 1:
			print "Only one atomtype present; using all atoms in PDB file"
			atomtype = None
			
		if len(chains) == 1:
			print "Only one chain present"
			chain = None	

		count = 1
		for line in (i for i in lines if i.startswith("ATOM  ")):
			latomtype = line[12:15].strip()
			lchain = line[21].strip() or None			

			if atomtype and atomtype != latomtype:
				continue
			if chain and chain != lchain:
				continue

			pos = int(line[23:27])
			x = [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())]
			points[pos] = tuple([i+self.makenoise() for i in x])
			count += 1

		if len(points) == 0:
			raise Exception, "No atoms found in PDB file! Are chain and atomtype correct?"
		
		print "Loaded %s C-a points"%len(points)
		
		return points





	def write_tsplib(self, filename=None):
		filename = filename or self.basename+".tsp"
		fout = open(filename, "w")
		
		# TSPLib format header
		header = [
			"NAME: %s"%self.filename,
			"TYPE: TSP",
			"COMMENT: %s"%self.filename,
			"DIMENSION: %s"%len(self.points),
			"EDGE_WEIGHT_TYPE: EXPLICIT",
			"EDGE_WEIGHT_FORMAT: FULL_MATRIX",
			""
		]
		
		fout.write("\n".join(header))

		# TSPlib expects edges to start at #1
		keyorder = sorted(self.points.keys())		
		
		if self.fixededges:
			fout.write("FIXED_EDGES_SECTION\n")
			for i in self.fixededges:
				fout.write("%s %s\n"%(keyorder.index(i[0])+1, keyorder.index(i[1])+1))
			fout.write("-1\n")

		fout.write("EDGE_WEIGHT_SECTION\n")

		for point1 in keyorder:
			row = []
			for point2 in keyorder:
				d = self.weighted.get((point1, point2))
				row.append(d)

			fout.write("%s\n"%" ".join(map(str, row)))


		fout.write("EOF")
		fout.close()
		
		return filename
	
	
	
	
	def ca_ram(self, path, filename=None):

		print "\n=== C-a Ramachandran ==="
		
		points = map(self.points.get, path)
		# filename = filename or self.basename+".out.angles"

		out=[]
		index=1
		while index < len(points)-2:
			#calculate Ca-Ca*-Ca angle
			v1=norm_vector(points[index-1],points[index])
			v2=norm_vector(points[index+1],points[index])
			dp=numpy.dot(v1,v2)
			if dp > 1:
				dp=1
			if dp<-1:
				dp=-1
			CaCaCa=math.acos(dp)*(180/math.pi)
			#print index+1, CaCaCa, 

		
			#calculate Ca-Ca*-Ca-Ca torsion 
			v1=norm_vector(points[index-1],points[index])
			v2=norm_vector(points[index+1],points[index])
			xp1=cross_product(v1,v2)
			lengthxp1=math.sqrt(xp1[0]**2+xp1[1]**2+xp1[2]**2)
			xp1n=(xp1[0]/lengthxp1, xp1[1]/lengthxp1, xp1[2]/lengthxp1)
	
			v3=norm_vector(points[index],points[index+1])
			v4=norm_vector(points[index+2],points[index+1])
			xp2=cross_product(v3,v4)
			lengthxp2=math.sqrt(xp2[0]**2+xp2[1]**2+xp2[2]**2)
			xp2n=(xp2[0]/lengthxp2, xp2[1]/lengthxp2, xp2[2]/lengthxp2) 
	
			dpxps=numpy.dot(xp2n,xp1n)
			if dpxps > 1:
				dp=1
			if dpxps<-1:
				dpxps=-1
			CaCaCaCa=math.acos(dpxps)*(180/math.pi)
			# print CaCaCa, CaCaCaCa
			out.append((CaCaCa, CaCaCaCa))
			index=index+1
		
		#f = open(filename, 'w')
		#for i in out:
		#	f.write("%s %s\n"%(i[0],i[1]))
		#f.close()
		return {'ca_angles': out}
		
		
	
# parse helices
# a = [i.strip().split() for i in helix.split("\n")]
# for i in b: print " ".join(map(str,range(int(i[5]), int(i[8])+1)))	
	

def main():
	usage = """%prog [options] <pdb file>
	Find paths between two atoms in a PDB model
	"""

	parser = optparse.OptionParser(usage=usage, version=EMAN2.EMANVERSION)
	parser.add_option("--start", type="int",help="Start ATOM")
	parser.add_option("--end", type="int",help="End ATOM")	
	parser.add_option("--average", type="float",help="Average Ca-Ca length", default=3.78)
	parser.add_option("--dmin", type="float",help="Mininum Ca-Ca length", default=2.0)
	parser.add_option("--dmax", type="float",help="Maximum Ca-Ca length", default=4.2)
	parser.add_option("--noise", type="float",help="Add Gaussian Noise", default=0.0)
	parser.add_option("--solver", type="str" ,help="Run TSP Solver: concorde or lkh")
	parser.add_option("--atomtype", type="str" ,help="Load Atom Type. Default: 'CA'. Options: 'C' or 'all'", default="CA")	
	parser.add_option("--chain", type="str" ,help="Load Chain. Default: load all chains")	
	parser.add_option("--edgefile", type="str" ,help="Load fixed fragment file; one sequence of forced connections per line, separated by space.")	
	parser.add_option("-e", "--edge", action="append", help="Forced edge: e.g. -e1,3")
	parser.add_option("--iterations", type="int",help="Iterations", default=1)	
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n",type="int", default=0, help='verbose level [0-9], higher number means higher level of verboseness')
	
	(options, args) = parser.parse_args()

	filename = args[0]
	for j in range(options.iterations):
		if options.iterations > 1:
			print "Iteration:", j
	
		pw = PathWalker(
			filename=filename, 
			start=options.start, 
			end=options.end, 
			edgefile=options.edgefile, 
			dmin=options.dmin, 
			dmax=options.dmax, 
			average=options.average, 
			noise=options.noise, 
			atomtype=options.atomtype, 
			chain=options.chain, 
			solver=options.solver
		)
		pw.run()



# If executed as a program
if __name__ == '__main__':
	main()








