#!/usr/bin/env python

import os, sys, commands
from sys import argv
from EMAN2 import *
import matplotlib.pyplot as plt

def main():
	name = argv[1]

	print "I am analyzing particle ", name
	particle = EMData(name,0)
	particle = particle.process("filter.lowpass.gauss",{"cutoff_freq":0.02})
	particle = particle.process('normalize')

	radius = particle.get_xsize()/2

	plot_name = name.replace('.hdf','_rdPLOT.png')

	values = particle.calc_radial_dist(radius, 0, 1, 1)

	print "The name for the plot is", plot_name
	
	txtfilename = plot_name.replace('.png','.txt')
	
	f = open(txtfilename,'w')
	lines = []
	for value in values:
		line = str(value) + '\n'
		lines.append(line)
	f.writelines(lines)
	f.close()
	
	print "The values calculated for the Radial Density of THIS particle are", values

	plt.plot(values, linewidth=1)
	plt.title("Radial density plot")
	plt.ylabel("Density (arbitrary units)")
	plt.xlabel("Radius [pixels]")
	a = plt.gca()
	a.set_xlim([0,radius])
	plt.savefig(plot_name)
	plt.clf()
	
	return()
	
if __name__ == '__main__':
	main()
