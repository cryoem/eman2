#!/usr/bin/env python

from string import atoi, atof, split
from sys import argv, exit



def histogram1d( data, nbin ) :
    fmax = max( data )
    fmin = min( data )

    binsize = (fmax - fmin)/(nbin-1)

    start = fmin-binsize/2.0

    region = [None]*nbin
    hist = [None]*nbin
    for i in xrange(nbin):
          region[i] = start + (i+0.5)*binsize
          hist[i] = 0

    for d in data:
        id = int( (d-start)/binsize )
        hist[id]+=1
    return region,hist



if len(argv) !=4:
    print "Usage: histogram file columnid nbin"
    exit(-1)

f = open( argv[1], 'r' )
cid = atoi( argv[2] )
nbin = atoi( argv[3] )

data = []
line = f.readline()
while len(line)>0:
    items = split( line )

    data.append( atof(items[cid-1]) )

    line = f.readline()

region,hist = histogram1d( data, nbin )

for i in xrange( len(hist) ):
    print '%6d %10.5f %6d' % (i, region[i], hist[i] )



