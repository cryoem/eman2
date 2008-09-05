#!/usr/bin/env python

from EMAN2 import *
from sparx import *
from time import time
from optparse import OptionParser

def bootstrap_genbuf( prjfile, bufprefix, beg, end, finfo=None ):
	start_time = time()
	istore = newfile_store( bufprefix, 4 )
	for i in xrange( beg, end ):
		prj = get_im( prjfile, i ) 
		istore.add_image( prj )
		if not(info is None) and (i%100==99):
			finfo.write( "%6d buffered, time: %10.3f\n" % (i, time()-start_time) )
			finfo.flush()


def bootstrap_insert( bufprefix, fftvols, wgtvols, mults, info=None):
	ostore = newfile_store( bufprefix, 4 )
	blocksize = 250
	nvol = len(fftvols)
	nprj = len(mults[0])
	nblock = (nprj-1)/blocksize + 1

	overall_start = time()

	for iblock in xrange(nblock):
    		if iblock==nblock - 1:
        		pbeg = iblock*blocksize
        		pend = nprj
    		else:
        		pbeg = iblock*blocksize
        		pend = pbeg + blocksize

		start_time = time()		
		ostore.read( pend - pbeg )
 		if not(info is None):
				info.write("    block %d read.\t time: %10.3f\n" % (iblock, time() - start_time) )
				info.flush()




		start_time = time()
    		for ivol in xrange(nvol):
        		start_time = time()
        		ostore.add_tovol( fftvols[ivol], wgtvols[ivol], mults[ivol], pbeg, pend )
       		if not(info is None):
			info.write("    block %d inserted.\t time: %10.3f\n" % (iblock, time() - start_time) )
			info.flush()
        

	if not(info is None):
		info.write("    Project inserted.\t time: %10.3f\n" % (time() - overall_start) )
		info.flush()

def bootstrap_finish( rectors, fftvols, wgtvols, idingrp, grpsize, grpcomm, volfile, iter, info=None ):
	from time import time

	nvol = len(fftvols)
	nround = nvol/grpsize

	assert nvol%grpsize==0

	for iround in xrange(nround):
		for icpu in xrange(grpsize):
			ivol = grpsize*iround + icpu
			start_time = time()
			reduce_EMData_to_root( fftvols[ivol], idingrp, icpu, grpcomm )
			reduce_EMData_to_root( wgtvols[ivol], idingrp, icpu, grpcomm )

			if not(info is None):
				info.write( "    fftvol %d reduced.\t time: %10.3f\n" % (ivol, time()-start_time) )
				info.flush()

	
		start_time = time()
		iwrite = nround*iter + iround

		ivol = grpsize*iround + idingrp
		vol = rectors[ivol].finish()
		vol.write_image( volfile, nround*iter+iround )
		info.write( "    vol %d reconstructed.\t time: %10.3f\n" % (ivol, time()-start_time) )


def bootstrap_prepare( prjfile, nvol, snr ):
	nx = get_im( prjfile, 0 ).get_xsize()
	fftvols = [None]*nvol
	wgtvols = [None]*nvol
	rectors = [None]*nvol
	for i in xrange(nvol):
		fftvols[i] = EMData()
		wgtvols[i] = EMData()
		params = {"size":nx, "npad":4, "snr":snr, "weight":wgtvols[i], "fftvol":fftvols[i]}
		rectors[i] = Reconstructors.get( "nn4_ctf", params )
		rectors[i].setup()

	return rectors, fftvols, wgtvols

def select_id( acc_prbs ):
    r = random()
    i = 0

    while i < len(acc_prbs) and r >= acc_prbs[i]:
        i+=1

    assert r < acc_prbs[i]
    return i


def bootstrap_mults( accu_prbs, nvol, beg, end):
	from random import randint
	mults = []
	nprj = len(accu_prbs)
	for ivol in xrange(nvol):
		m = [0]*nprj
		for i in xrange(nprj):
			ip = select_id( accu_prbs )
			m[ip] += 1
	
		assert( sum(m)==nprj )
		mults.append( m[beg:end] )

	return mults

def prepare_wgts( wgts ):
	nprj = len(wgts)
	accu_prbs = [0.0] * nprj
	wsum = 0.0
	for i in xrange(nprj):
		wsum += wgts[i]
		accu_prbs[i] = wsum

	for i in xrange(nprj):
		accu_prbs[i] /= wsum

	return accu_prbs 		


def bootstrap( prjfile, wgts, volprefix, bufprefix, nprj, niter, nbufvol, snr, genbuf ) :
	from mpi import mpi_comm_rank, mpi_comm_size, mpi_comm_split, MPI_COMM_WORLD
	from random import seed

	if( nprj==-1) :
		nprj = EMUtil.get_image_count( prjfile )
	myid = mpi_comm_rank( MPI_COMM_WORLD )
	ncpu = mpi_comm_size( MPI_COMM_WORLD )
	accu_prbs = prepare_wgts( wgts )


	finfo=open( "progress%04d.txt" % myid, "w" )

	if genbuf:
		bootstrap_genbuf( prjfile, bufprefix, 0, nprj, finfo )
	
	volfile = "%s%04d.hdf" % (volprefix, myid)
	for iter in xrange(niter):
		finfo.write( "Iteration %d: \n" % iter )
		finfo.flush()

		iter_start = time()
		mults = bootstrap_mults(accu_prbs, nbufvol, 0, nprj)

		assert len(mults)==nbufvol
		rectors, fftvols, wgtvols = bootstrap_prepare( prjfile, nbufvol, snr )
		bootstrap_insert( bufprefix, fftvols, wgtvols, mults, finfo )
		bootstrap_finish( rectors, fftvols, wgtvols, myid, ncpu, MPI_COMM_WORLD, volfile, iter, finfo )
		
		finfo.write( "time for iteration: %10.3f\n" % (time() - iter_start) )
		finfo.flush()

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


def bootstrap_calcwgts( prjfile, wgtfile, delta ):
	from mpi import mpi_comm_rank, mpi_comm_size, mpi_reduce
	from mpi import MPI_INT, MPI_SUM, MPI_COMM_WORLD
	

	myid = mpi_comm_rank( MPI_COMM_WORLD )
	ncpu = mpi_comm_size( MPI_COMM_WORLD )

	finf = open( "calcwgts%04d.txt" % myid, "w" )
	nprj = EMUtil.get_image_count( prjfile )

	beg, end = MPI_start_end( nprj, ncpu, myid )
	finf.write( "begin, end: %d %d\n" %(beg, end) )
	finf.flush()

	eve_angs = even_angles( delta, 0.0, 90.0 )
	eve_vecs = [None]*len(eve_angs)
	for i in xrange( len(eve_angs) ):
		eve_vecs[i] = getvec( eve_angs[i][0], eve_angs[i][1] )

	nang = len(eve_angs)
	occurs = [0]*nang
	angids = [0]*nprj

	for iprj in xrange(beg, end):
		prj = get_im( prjfile, iprj )
		phi = prj.get_attr( "phi" )
		tht = prj.get_attr( "theta" )
		aid = nearest_ang( eve_vecs, phi, tht )

		finf.write( "prj %6d: %10.3f %10.3f close to ori %6d: %10.3f %10.3f\n" % (iprj, phi, tht, aid, eve_angs[aid][0], eve_angs[aid][1]))
		finf.flush()

		angids[iprj] = aid
		occurs[aid] += 1

	occurs = mpi_reduce( occurs, len(occurs), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD )
	angids = mpi_reduce( angids, len(angids), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD )

	if myid==0:
		sumoccur = 0
		for i in xrange( nang ):
			finf.write( "ori %6d: %d\n" % (i, occurs[i]) )
			finf.flush()
			sumoccur += int(occurs[i])
		assert sumoccur==nprj

		for i in xrange( len(occurs) ):
			if int(occurs[i])==0:
				print "Error: some orientation has not projections. increase delta"
				return None	


		os = open( wgtfile, "w" )
		for i in xrange(nprj):
			aid = angids[i]
			wgt = 1.0/int(occurs[aid])/nang
			print >> os, wgt
		



def main():
	from mpi import mpi_init
	import sys

	sys.argv = mpi_init( len(sys.argv), sys.argv )

        arglist = []
        for arg in sys.argv:
	    arglist.append( arg )

	progname = os.path.basename(arglist[0])
	usage = progname + " prjstack wgtfile [volprefix bufprefix] --niter --nbufvol --nprj --snr --genbuf --calcwgts --delta"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--nprj", type="int", default=-1, help="  number of projections used by default all projections will be used")
	parser.add_option("--snr",  type="float", default=1.0, help="signal-to-noise ratio" )
	parser.add_option("--delta", type="float", default=1.0, help="for weights calculation")
        parser.add_option("--genbuf", action="store_true", default=False, help="whether generating buffer")
	parser.add_option("--calcwgt", action="store_true", default=False, help="calculate weights or load weights from file" )
	parser.add_option("--bufprefix", type="string", help="the location of the scratch file" )
	parser.add_option("--nbufvol", type="int", help="number of fftvol in the memory" )
	parser.add_option("--niter", type="int", help="number of iteration. Each iteration, nbufvol of bootstrap volumes will be generated.")

	(options, args) = parser.parse_args( arglist[1:] )

	if( len(args) !=2 ):
		print "usage: " + usage
		return None

	prjfile = args[0]

	if( options.calcwgt ):
		wgtfile = args[1]
		bootstrap_calcwgts( prjfile, wgtfile, options.delta )
	else :
		wgts = read_text_file( argv[1], 0 )
		volprefix = args[2]
		bootstrap( prjfile, wgts, volprefix, bufprefix, options.nprj, options.niter, options.nbufvol, options.snr, options.genbuf )

if __name__ == "__main__":
	main()

