#!/usr/bin/env python

from string import atof, split, atoi
from math import sqrt

from sys import argv,exit

if len(argv)!=3:
    print "Usage: txtccc.py txt1:columnid txt2:columnid"
    exit(-1)


file1,sid1 = split( argv[1], ':' )
file2,sid2 = split( argv[2], ':' )

f1 = open( file1, 'r' )
f2 = open( file2, 'r' )

id1 = atoi( sid1 )
id2 = atoi( sid2 )


sumx2 = 0.0
sumy2 = 0.0
sumxy = 0.0
sumx  = 0.0
sumy  = 0.0
n = 0

line1 = f1.readline()
while len(line1)>0:
    line2 = f2.readline()
    if len(line2)==0:
        print "Error: two files different length"
	exit(-2)

    items1 = split( line1 )
    items2 = split( line2 )

    x = atof( items1[id1-1] )
    y = atof( items2[id2-1] )

    sumx += x
    sumy += y
    sumx2 += x*x
    sumy2 += y*y
    sumxy += x*y
    n += 1
    line1 = f1.readline()


avgx2 = sumx2/n
avgy2 = sumy2/n
avgxy = sumxy/n

avgx = sumx/n
avgy = sumy/n


sigmax = sqrt( avgx2 - avgx*avgx )
sigmay = sqrt( avgy2 - avgy*avgy )

ccc = ( avgxy - avgx*avgy ) / sigmax / sigmay

a = (avgxy - avgx*avgy) / (sigmax*sigmax)

b = avgy - a*avgx


print 'n, ccc, a, b: %4d %10.5f %10.5e %10.5e' % (n, ccc, a, b)



