#!/usr/bin/env python
#
# Author: John Flanagan, 04/10/2012 (jfflanag@bcm.edu)
# Copyright (c) 2000-2012 Baylor College of Medicine
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
#
#

def db_compute_fsc(a, b, apix, pathname, dbname):
	fsc = a.calc_fourier_shell_correlation(b)
	third = len(fsc)/3
	xaxis = fsc[0:third]
	plot = fsc[third:2*third]
	error = fsc[2*third:]
		
	convergence_db_name = "bdb:"+pathname+"#convergence.results"
	db = db_open_dict(convergence_db_name)
	
	tmpaxis = [x/apix for x in xaxis]
	db[dbname] = [tmpaxis,plot]
	#db["error_"+s+"_fsc"] = [xaxis,error] #we're not plotting the errors
	db_close_dict(convergence_db_name)
	
	
def get_e2refine_even_odd_results_list(keys):
	'''
	Extract the names from the keys that match the e2resolution.py output naming convention
	(keys is a list of keys in the convergence.results dictionary, in a refinement directory)
	'''
	solns = []
	for k in keys:
		if k[0:13] == "conv_even_odd":
			solns.append(k)
	solns.sort()
	return solns
	
def get_e2eotest_results_list(keys):
	'''
	Extract the names from the keys that match the e2eotest.py output naming convention
	(keys is a list of keys in the convergence.results dictionary, in a refinement directory)
	'''
	solns = []
	for k in keys:
		if len(k) > 7 and k[0:8] == "even_odd":
			solns.append(k)
	solns.sort()
	return solns
	
def get_e2resolution_results_list(keys):
	'''
	Extract the names from the keys that match the e2resolution.py output naming convention
	(keys is a list of keys in the convergence.results dictionary, in a refinement directory)
	'''
	solns = []
	for k in keys:
		if len(k) > 6 and k[-7:] == "res_fsc":
			solns.append(k)
	solns.sort()
	return solns
	
def get_convergence_results_list(keys):
	'''
	Extract the names from the keys that match the e2refine.py convergence plot output naming convention
	(keys is a list of keys in the convergence.results dictionary, in a refinement directory)
	'''
	solns = []
	if "init_00_fsc" in keys:
		solns.append("init_00_fsc")
		
	i = 0
	while True:
		s1 = str(i)
		s2 = str(i+1)
		if len(s1) == 1: s1 = "0"+s1
		if len(s2) == 1: s2 = "0"+s2
		k = s1+"_"+s2+"_fsc"
		if k in keys:
			solns.append(k)
		else:
			break

		i += 1

	return solns