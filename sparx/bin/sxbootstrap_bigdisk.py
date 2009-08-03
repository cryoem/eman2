#!/usr/bin/env python

from EMAN2 import *
from sparx import *
from time import time
from optparse import OptionParser

def bootstrap_genbuf( prjfile, bufprefix, beg, end, CTF, npad, info=None ):
	start_time = time()
	istore = newfile_store( bufprefix, npad, CTF )
	for i in xrange( beg, end ):
		prj = get_im( prjfile, i )
		istore.add_image( prj, prj.get_attr("xform.projection") )
		if( not(info is None) and ((i%100==99 or i==end-1))):
			finfo.write( "%6d buffered, time: %10.3f\n" % (i+1, time()-start_time) )
			finfo.flush()


def bootstrap_insert( bufprefix, fftvols, wgtvols, mults, CTF, npad, info=None):
	ostore = newfile_store( bufprefix, npad, CTF )
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
			t = time()
			info.write("        block %d read.   \t time: %10.3f %10.3f\n" % (iblock, t-start_time, t-overall_start) )
			info.flush()

		start_time = time()
    		for ivol in xrange(nvol):
        		ostore.add_tovol( fftvols[ivol], wgtvols[ivol], mults[ivol], pbeg, pend )
       		if not(info is None):
			t = time()
			info.write("        block %d inserted.\t time: %10.3f %10.3f\n" % (iblock, t-start_time, t-overall_start) )
			info.flush()
        

	if not(info is None):
		info.write("    Project inserted.\t time: %10.3f\n" % (time() - overall_start) )
		info.flush()

def bootstrap_finish( rectors, fftvols, wgtvols, volfile, iter, info=None ):
	from time import time
	overall_start = time()
	nvol = len(fftvols)
	for ivol in xrange(nvol):
		start_time = time()
		iwrite = nvol*iter + ivol

		vol = rectors[ivol].finish()
		vol.write_image( volfile, iwrite )
		if not(info is None):
			t = time()
			info.write( "        vol %d reconstred.\t time: %10.3f %10.3f\n" % (ivol, t-start_time, t-overall_start) )
			info.flush()

	if not(info is None):
		info.write( "    Volume finished.\t time: %10.3f\n" % (time() - overall_start) )



def bootstrap_prepare( prjfile, nvol, snr, CTF, npad ):
	nx = get_im( prjfile, 0 ).get_xsize()
	fftvols = [None]*nvol
	wgtvols = [None]*nvol
	rectors = [None]*nvol
	for i in xrange(nvol):
		fftvols[i] = EMData()
		wgtvols[i] = EMData()
		if CTF:
			params = {"size":nx, "npad":npad, "snr":snr, "weight":wgtvols[i], "fftvol":fftvols[i]}
			rectors[i] = Reconstructors.get( "nn4_ctf", params )
		else:
			params = {"size":nx, "npad":npad, "snr":snr, "weight":wgtvols[i], "fftvol":fftvols[i]}
			rectors[i] = Reconstructors.get( "nn4", params )
		
		rectors[i].setup()

	return rectors, fftvols, wgtvols

def select_id( acc_prbs, r=None ):
	from random import random
	if r is None:
		r = random()

	if r < acc_prbs[0]:
		return 0

	assert r <= acc_prbs[-1]

	n = len(acc_prbs)

	curt = n/2
	size = n/4

	while True :
		assert curt>=0 and curt < n
		if r < acc_prbs[curt-1]:
			curt -= size
		elif r >= acc_prbs[curt]:
			curt += size
		else:
			return curt

		size = size/2
		if size==0: 
			size = 1

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

def write_mults( fmults, iter, mults ):
	fmults.write( "iter: %d\n" % iter )

	for ivol in xrange( len(mults) ):
		fmults.write( "vol: %d\n" % ivol )
		for i in xrange( len(mults[ivol]) ):
			fmults.write( "%4d " % mults[ivol][i] )
			if (i+1)%16==0:
				fmults.write( "\n" )
	
	fmults.write( "\n" )
	fmults.flush()


def bootstrap( prjfile, wgts, outdir, bufprefix, nbufvol, nvol, seedbase, snr, genbuf, ngroup, CTF, npad, MPI, verbose = 0 ) :
	from random import seed
	import os

	nprj = EMUtil.get_image_count( prjfile )

	if MPI:
		from mpi import mpi_barrier, mpi_comm_rank, mpi_comm_size, mpi_comm_split, MPI_COMM_WORLD
		myid = mpi_comm_rank( MPI_COMM_WORLD )
		ncpu = mpi_comm_size( MPI_COMM_WORLD )
	else:
		myid = 0
		ncpu = 1

	accu_prbs = prepare_wgts( wgts )


	if myid==0:
		if os.path.exists(outdir):
			os.system( "rm -rf " + outdir )
		os.system( "mkdir " + outdir )

	if MPI:
		mpi_barrier( MPI_COMM_WORLD )

	if(verbose == 1):  finfo=open( os.path.join(outdir, "progress%04d.txt" % myid), "w" )
	else:              finfo = None

	groupsize = ncpu/ngroup
	if genbuf and (myid%groupsize==0):
		bootstrap_genbuf( prjfile, bufprefix, 0, nprj, CTF, npad, finfo )
	
	if genbuf and MPI:
		mpi_barrier( MPI_COMM_WORLD )

	myseed = seedbase + 10*myid

	seed( seedbase + 10*myid )
	volfile = os.path.join(outdir, "bsvol%04d.hdf" % myid)

	niter = nvol/ncpu/nbufvol
	for iter in xrange(niter):
		if(verbose == 1):
			finfo.write( "Iteration %d: \n" % iter )
			finfo.flush()

		iter_start = time()
		mults = bootstrap_mults(accu_prbs, nbufvol, 0, nprj)

		assert len(mults)==nbufvol
		rectors, fftvols, wgtvols = bootstrap_prepare( prjfile, nbufvol, snr, CTF, npad )
		bootstrap_insert( bufprefix, fftvols, wgtvols, mults, CTF, npad, finfo )
		bootstrap_finish( rectors, fftvols, wgtvols, volfile, iter, finfo )
		rectors = None
		fftvols = None
		wgtvols = None
		if(verbose == 1):
			finfo.write( "time for iteration: %10.3f\n" % (time() - iter_start) )
			finfo.flush()

def main():

	import sys

        arglist = []
        for arg in sys.argv:
	    arglist.append( arg )

	progname = os.path.basename(arglist[0])
	usage = progname + " prjstack wgtfile outdir bufprefix --nvol --nbufvol --seedbase --snr --genbuf --ngroup --npad --CTF"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--nvol",     type="int",                         help="number of bootstrap volumes to be generated.")
	parser.add_option("--nbufvol",  type="int",                         help="number of fftvol in the memory" )
	parser.add_option("--genbuf",   action="store_true", default=False, help="whether generate the buffer")
	parser.add_option("--ngroup",   type="int",          default=1,     help="how many groups to use (each group will share a buffer)")
	parser.add_option("--CTF",      action="store_true", default=False, help="whether consider CTF" )
	parser.add_option("--snr",      type="float",        default=1.0,   help="signal-to-noise ratio" )
	parser.add_option("--npad",     type="int",          default=4,     help="times of padding" )
	parser.add_option("--seedbase", type="int",                         help="random seed base" )
	parser.add_option("--MPI",      action="store_true", default=False, help="use MPI")
	parser.add_option("--verbose",  type="int",          default=0,     help="verbose level: 0 no, 1 yes" )

	(options, args) = parser.parse_args( arglist[1:] )

	if( len(args) !=2 and len(args) != 4):
		print "usage: " + usage
		return None

	prjfile = args[0]

	if options.MPI and options.genbuf:   ERROR('Generation of buffer does not have MPI version, swith off MPI flag', "sxbootstrap_bigdisk", 1)
	if options.MPI:
		from mpi import mpi_init, mpi_comm_rank, mpi_comm_size, MPI_COMM_WORLD
		sys.argv = mpi_init( len(sys.argv), sys.argv )

		from utilities import init_mpi_bdb
		init_mpi_bdb()

	wgts = read_text_file( args[1], 0 )
	outdir = args[2]
	bufprefix = args[3]
	bootstrap( prjfile, wgts, outdir, bufprefix, options.nbufvol, options.nvol, options.seedbase, options.snr, options.genbuf, options.ngroup, options.CTF, options.npad, options.MPI, options.verbose )


if __name__ == "__main__":
	main()

