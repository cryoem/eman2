#!/bin/env python

# $Id$

# todo: verify the filters who have the same names in proc3d
#       and proc2d have the same implementation
#
# todo: lp, hp, tlp vs apix

from EMAN2 import *
from optparse import OptionParser
import sys
import os.path

def main():
    progname = os.path.basename(sys.argv[0])
    usage = progname + " options inputfile outputfile"
    parser = OptionParser(usage)

    parser.add_option("--shrink", type="int", action="append", help=".")
    parser.add_option("--scale", type="float", action="append", help="")
    parser.add_option("--clip", type="float", nargs=6, action="append", help="")
    parser.add_option("--clipc", type="float", nargs=3, action="append", help="")
    parser.add_option("--fftclip", type="float", nargs=3, action="append", help="")
    parser.add_option("--filter", type="string", action="append",
                      help="apply a filter. FILTER=filtername:param1=val1:param2=val2")
    
    parser.add_option("--apix", type="float", help="")
    parser.add_option("--origin", type="float", nargs=3, help="")
    parser.add_option("--mult", type="float", help="")
    parser.add_option("--add", type="float", help="")
    
    parser.add_option("--calcsf", type="string", nargs=2, help="")
    parser.add_option("--tophalf", action="store_true", help="")
    parser.add_option("--icos5fhalfmap", action="store_true", help="")
    parser.add_option("--outtype", type="string", help="output image format, mrc, imagic, hdf, etc")
    
    append_options = ["clip", "filter", "clipc", "fftclip", "shrink", "scale"]
    optionlist = []
    
    for arg1 in sys.argv[1:]:
        if arg1[0] == "-":
            optionlist.append(arg1.lstrip("-"))

    (options, args) = parser.parse_args()
    
    if len(args) != 2:
        print "usage: " + usage
        print "Please run '" + progname + " -h' for detailed options"
        sys.exit(1)
    
    infile = args[0]
    outfile = args[1]

    index_d = {}
    for append_option in append_options:
        index_d[append_option] = 0

    data = EMData()
    data.read_image(infile)
    if options.apix > 0:
        data.set_xpixel(apix)
        data.set_ypixel(apix)
        data.set_zpixel(apix)
        
    if options.origin:
        data.set_xyz_origin(options.origin[0], options.origin[1], options.origin[2])

    if options.icos5fhalfmap:
        print "not implemnted yet"

    x = data.get_xsize()
    y = data.get_ysize()
    z = data.get_zsize()
    
    print "Original image : %dx%dx%d Mean=%1.3g Sigma=%1.3g Min=%1.3g Max=%1.3g\n" % \
    (x, y, z, data.get_mean(), data.get_sigma(), data.get_min(), data.get_max())
    
    if options.calcsf:
        dataf = data.do_fft()
        data.gimme_fft()
        curve = dataf.calc_radial_func(y, 0, 0.5)
        Util.save_data(0, 1.0/(apix*2.0*y), curve, options.calcsf);

    xc = -1
    yc = -1
    zc = -1

    nx = x
    ny = y
    nz = z
    
    if options.clip:
        (nx, ny, nz, xc, yc, zc) = options.clip
        
    if not (xc>=0 and yc>=0 and zc>=0 and xc<x and yc<y and zc<z):
        xc = x/2
        yc = y/2
        zc = z/2

    if options.scale > 1:
        data = data.clip(xc-nx/2, yc-ny/2,zc-nz/2,nx,ny,nz)

    if options.scale != 1:
        data.scale(options.scale)
        
    if options.scale < 1 || (scale == 1 && (x!=nx || y!=ny ||z!=nz)):
        data = data.clip(xc-nx/2,yc-ny/2,zc-nz/2,nx,ny,nz)

    if options.shrink > 1:
        data.median_shrink(options.shrink)
        nx = data.get_xsize()
        ny = data.get_ysize()
        nz = data.get_zsize()

    if options.mult:
        data.mult(options.mult)

    if options.add:
        data.add(options.add)

    if options.fftclip:
        fnx = options.fftclip[0]
        fny = options.fftclip[1]
        fnz = options.fftclip[2]
        
        fft = data.do_fft()
        padfft = fft.clip(0, 0, 0, fnx+2, fny, fnz)
        data = padfft->do_ift()
        
    # next tophalf
    if options.tophalf:
        half = data.get_top_half()
        data = half

    if options.outtype:
        data.write_image(outfile, 0, EMUtil.get_image_ext_type(options.outtype))
        


if __name__ == "__main__":
    main()
    

