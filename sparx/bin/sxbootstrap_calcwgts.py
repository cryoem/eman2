#!/usr/bin/env python

from EMAN2 import *
from sparx import *
from time import time
from optparse import OptionParser
from math import pi,sin,cos

angle_to_rad = pi/180.0

def getvec( phi, tht ):

	if tht > 180.0:
		tht = tht - 180.0
		phi = phi + 180.0
	elif tht > 90.0:
		tht = 180.0 - tht
		phi = phi + 180.0

	assert tht <=90.0

	tht *= angle_to_rad
	phi *= angle_to_rad

	x = sin(tht)*cos(phi) 
	y = sin(tht)*sin(phi)
	z = cos(tht)

	return (x,y,z)

def nearest_ang( vecs, phi, tht ) :
	vec = getvec( phi, tht )

	best_s = -1.0
	best_i = -1

	for i in xrange( len(vecs) ):
		s = vecs[i][0]*vec[0] + vecs[i][1]*vec[1] + vecs[i][2]*vec[2]
		if s > best_s:
			best_s = s
			best_i = i

	return best_i

def cmpccc( a, b ):
	if( a[0] > b[0] ) :
		return 1

	if( a[0] < b[0] ) :
		return -1

	return 0


def exclude( prjccc, nexclude, bothside ):
	nprj = len(prjccc)
	tmp = [None] * nprj
	for i in xrange(nprj):
		tmp[i] = [ float(prjccc[i]), i ]

	tmp.sort( cmpccc )

	exclist = []
	for i in xrange( min(nprj,nexclude) ):
		print "ATTENTION: prj #%6d ccc %5.2f will be excluded" % ( tmp[i][1], tmp[i][0] )
		exclist.append( tmp[i][1] )
		

	if bothside:
		nexclude = min(nprj, nexclude)
		for i in xrange( nexclude ):
			id = nprj - nexclude + i

			if exclist.count(id) > 0 : 
				continue

			print "ATTENTION: prj #%6d ccc %5.2f will be excluded" % ( tmp[id][1], tmp[id][0] )
			exclist.append( tmp[id][1] )
	
	return exclist







def bootstrap_calcwgts( prjfile, wgtfile, voronoi, delta, refvol=None, fl=None, fh=None, CTF=False, nexclude=0, bothside=False, MPI=False, verbose=True ):
	from projection import prep_vol,prgs

	if MPI:
		from mpi import mpi_comm_rank, mpi_comm_size, mpi_reduce
		from mpi import MPI_INT, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD
		myid = mpi_comm_rank( MPI_COMM_WORLD )
		ncpu = mpi_comm_size( MPI_COMM_WORLD )
	else:
		myid = 0
		ncpu = 1

	if verbose:
		finf = open( "calcwgts%04d.txt" % myid, "w" )
	nprj = EMUtil.get_image_count( prjfile )

	beg, end = MPI_start_end( nprj, ncpu, myid )
	if verbose:
		finf.write( "begin, end: %d %d\n" %(beg, end) )
		finf.flush()

	if voronoi:
		angs = [0.0] * (2*nprj)
		for iprj in xrange(beg,end):
			prj = get_im( prjfile, iprj )
			phi,tht,psi,s2x,s2y = get_params_proj( prj )
			angs[2*iprj] = phi
			angs[2*iprj+1] = tht
		
		if MPI:
			angs = mpi_reduce( angs, 2*nprj, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD )
		if myid == 0:
			tmp = [0.0]*len(angs)
			for i in xrange( len(angs) ):
				tmp[i] = float(angs[i])
			wgts = Util.cml_weights( tmp )
	
	else:
		eve_angs = even_angles( delta, 0.0, 89.99, method='P' )
		eve_vecs = [None]*len(eve_angs)
		for i in xrange( len(eve_angs) ):
			eve_vecs[i] = getvec( eve_angs[i][0], eve_angs[i][1] )

		nang = len(eve_angs)
		occurs = [0]*nang
		oriids = [0]*nprj


	if not(refvol is None):
		volft,kb = prep_vol( get_image(refvol) )
		prjccc = [0.0]*nprj


	for iprj in xrange(beg, end):

		if not voronoi:
			prj = get_im( prjfile, iprj )
			phi,tht,psi,s2x,s2y = get_params_proj( prj )
			aid = nearest_ang( eve_vecs, phi, tht )

			if verbose:
				finf.write( "prj %6d: %10.3f %10.3f close to ori %6d: %10.3f %10.3f\n" % (iprj, phi, tht, aid, eve_angs[aid][0], eve_angs[aid][1]))
				finf.flush()

			oriids[iprj] = aid
			occurs[aid] += 1

		if not(refvol is None):
			refprj = prgs( volft, kb, [phi,tht,psi,-s2x,-s2y] )
			if( fl is None or fh is None ):
				fltprj = prj
			else:
				fltprj = filt_tanl( prj, fl, fh )

			if CTF:
				ctf_params = prj.get_attr( "ctf" )
				refprj = filt_ctf( refprj, ctf_params )

			prjccc[iprj] =  ccc( refprj, fltprj )
	if MPI:
		if not(voronoi):
			occurs = mpi_reduce( occurs, len(occurs), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD )
			oriids = mpi_reduce( oriids, len(oriids), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD )

		if not(refvol is None):
			prjccc = mpi_reduce( prjccc, len(prjccc), MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD )

	if myid==0:

		if not(refvol is None):
			exclist = exclude( prjccc, nexclude, bothside )

			for id in exclist:
				if voronoi:
					wgts[id] = -.0
				else:
					oid = oriids[id]
					occurs[oid] -= 1
					oriids[id] = -1


		if verbose and (not voronoi):
			for i in xrange( nang ):
				finf.write( "ori %6d: %d\n" % (i, occurs[i]) )
				finf.flush()
		'''
		for i in xrange( len(occurs) ):
			if int(occurs[i])==0:
				#print "Warning: at least one orientation has no projections assigned. increase delta"
				pass
		'''

		os = open( wgtfile, "w" )
		for i in xrange(nprj):

			if voronoi:
				wgt = wgts[i]
			else:
				aid = oriids[i]

				if aid == -1:
					wgt = 0.0
				else:
					wgt = 1.0/int(occurs[aid])
			print >> os, wgt
		



def main():

	import sys

        arglist = []
        for arg in sys.argv:
	    arglist.append( arg )

	progname = os.path.basename(arglist[0])
	usage = progname + " prjstack wgtfile [refvol] --voronoi --MPI --CTF --delta --fl --fh --nexclude --exclude_bothside"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--voronoi", action="store_true", default=False, help="use voronoi diagram to create weights")
	parser.add_option("--delta", type="float", default=1.0, help="for weights calculation")
	parser.add_option("--CTF", action="store_true", default=False, help="whether consider CTF" )
	parser.add_option("--MPI", action="store_true", default=False, help="whether running mpi version")
	parser.add_option("--fl",  type="float",   default=None, help="before calculate projection apply filtration to projection" )
	parser.add_option("--fh",  type="float",   default=None, help="before calculate projection apply filtration to projection" )
	parser.add_option("--nexclude", type="int", default=None, help="number of projections (with lowest CCC) to be excluded." )
	parser.add_option("--exclude_bothside", action="store_true", default=False, help="exclude also projection with highest CCC." )


	(options, args) = parser.parse_args( arglist[1:] )

	if( len(args) !=2 and len(args) != 3):
		print "usage: " + usage
		return None

	prjfile = args[0]
	wgtfile = args[1]
	if len(args)==3:
		refvol = args[2]
	else:
		refvol = None

	if options.MPI:
		from mpi import mpi_init, mpi_comm_rank, mpi_comm_size, MPI_COMM_WORLD
		sys.argv = mpi_init( len(sys.argv), sys.argv )
	
		from utilities import init_mpi_bdb
		init_mpi_bdb()

	bootstrap_calcwgts( prjfile, wgtfile, options.voronoi, options.delta,refvol, options.fl, options.fh, options.CTF, options.nexclude, options.exclude_bothside, options.MPI )


if __name__ == "__main__":
	main()

