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
		if not(info is None):
			finfo.write( "%6d buffered, time: %10.3f\n" % (i, time()-start_time) )
			finfo.flush()


def bootstrap_insert( bufprefix, fftvols, wgtvols, mults, info=None):
	ostore = newfile_store( bufprefix, 4 )
	blocksize = 250
	nvol = len(fftvols)
	nprj = len(mults)
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


    		for ivol in xrange(nvol):
        		start_time = time()
        		ostore.add_tovol( fftvols[ivol], wgtvols[ivol], mults, pbeg, pend )
        		#if not(info is None):
			#	info.write("    Volume #%d, prj (%6d, %6d).\t time: %10.3f\n" % (ivol, pbeg, pend, time() - start_time) )
			#	info.flush()

	if not(info is None):
		info.write("    Project inserted.\t time: %10.3f\n" % (time() - overall_start) )
		info.flush()

def bootstrap_finish( rectors, fftvols, wgtvols, myid, comm, info=None ):
	from time import time

	nvol = len(fftvols)
	for ivol in xrange(nvol):
		start_time = time()
		reduce_EMData_to_root( fftvols[ivol], myid, ivol, comm )
		reduce_EMData_to_root( wgtvols[ivol], myid, ivol, comm )
		if not(info is None):
			info.write( "    fftvol %d reduced.\t time: %10.3f\n" % (ivol, time()-start_time) )
			info.flush()

	
	if myid < nvol:
		start_time = time()
		vol = rectors[myid].finish()
		info.write( "    vol %d reconstructed.\t time: %10.3f\n" % (myid, time()-start_time) )
		return vol
	else:
		return None


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



def bootstrap( prjfile, volprefix, nprj, snr, genbuf ) :
	from mpi import mpi_comm_rank, mpi_comm_size, mpi_comm_split, MPI_COMM_WORLD

	if( nprj==-1) :
		nprj = EMUtil.get_image_count( prjfile )
	myid = mpi_comm_rank( MPI_COMM_WORLD )
	ncpu = mpi_comm_size( MPI_COMM_WORLD )
	bufprefix="/scratch/tmp/tmpslice%04d" % myid
	grpsize = 4
	idingrp = myid%grpsize
	nvol = 4

	assert ncpu%grpsize==0

	finfo = open( "progress%04d.txt" % myid, "w" )


	beg, end = MPI_start_end( nprj, grpsize, idingrp )
	if genbuf:
		bootstrap_genbuf( prjfile, bufprefix, beg, end, finfo )

	finfo.write( "begin, end: %6d %6d\n" %(beg, end) )
	finfo.flush()

	# bootstrap calculation
	grpcomm = mpi_comm_split( MPI_COMM_WORLD, myid/grpsize, myid%grpsize )

	
	niter = 1
	volfile = "%s%04d.hdf" % (volprefix, myid)
	for iter in xrange(niter):
		finfo.write( "Iteration %d: \n" % iter )
		finfo.flush()

		iter_start = time()
		mults = [1]*nprj
		rectors, fftvols, wgtvols = bootstrap_prepare( prjfile, nvol, snr )
		bootstrap_insert( bufprefix, fftvols, wgtvols, mults[beg:end], finfo )
		vol = bootstrap_finish( rectors, fftvols, wgtvols, idingrp, grpcomm, finfo )
		if not(vol is None):
			vol.write_image( volfile, iter )
		
		finfo.write( "time for iteration: %10.3f\n" % (time() - iter_start) )
		finfo.flush()


def main():
	from mpi import mpi_init
	import sys

	sys.argv = mpi_init( len(sys.argv), sys.argv )

        arglist = []
        for arg in sys.argv:
	    arglist.append( arg )

	progname = os.path.basename(arglist[0])
	usage = progname + " prjstack volprefix --nprj --snr --genbuf"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--nprj", type="int", default=-1, help="  number of projections used by default all projections will be used")
	parser.add_option("--snr",  type="float", default=1.0, help="signal-to-noise ratio" )
        parser.add_option("--genbuf", action="store_true", default=False, help="whether generating buffer")

	(options, args) = parser.parse_args( arglist[1:] )

	if( len(args) !=2 ):
		print "usage: " + usage
		return None

	prjfile = args[0]
	volprefix = args[1]
	bootstrap( prjfile, volprefix, options.nprj, options.snr, options.genbuf )

if __name__ == "__main__":
	main()
 
'''
v1 = getImage( "vol_direct.hdf" )
for ivol in xrange(nvol):
    start_time = time()
    v = reconstructors[ivol].finish()
    print "Volume #", ivol, " reconstructed. time: ", time() - start_time 

    print "ccc(v1,v2): ", ccc( v1, v )
'''
