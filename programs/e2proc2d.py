#!/bin/env python

# $Id$

from EMAN2 import *
from optparse import OptionParser
import sys
import os.path
import math
import random

# usage: e2proc2d.py --clip 200 200 --outtype mrc input.mrc output.mrc

def read_listfile(listfile, excludefile, nimg):
    imagelist = None
    infile = None
    exclude = 0
    
    if listfile:
        infile = listfile
    elif excludefile:
        infile = excludefile
        exclude = 1
        
    if infile:
        try:
            lfp = open(infile, "rb")
        except IOError:
            print "Error: couldn't read list file '%s'" % infile
            sys.exit(1)
        
        if exclude:
            imagelist = [1] * nimg
        else:
            imagelist = [0] * nimg

        imagelines = lfp.readlines()
        for line in imagelines:
            if line[0] != "#":
                n = int(line)
                if n >= 0 and n < nimg:
                    if options.exclude:
                        imagelist[n] = 0
                    else:
                        imagelist[n] = 1
    return imagelist



def main():
    progname = os.path.basename(sys.argv[0])
    usage = progname + " options inputfile outputfile"

    parser = OptionParser(usage,version=EMANVERSION)
  
    parser.add_option("--apix", type="float", help="the Angstrom/pixel for S scaling")
    parser.add_option("--average", action="store_true", help="Averages all input images (without alignment) and writes a single (normalized) output image")
    parser.add_option("--calcsf", type="string", nargs=2, help="calculate a radial structure factor for the image and write it to the output file, must specify apix. divide into <n> angular bins")    
    parser.add_option("--clip", type="float", nargs=2, action="append", help="Define the output image size. CLIP=xsize ysize")
    parser.add_option("--ctfsplit", action="store_true", help="Splits the input file into output files with the same CTF parameters")
    parser.add_option("--exclude", type="string", help="Excludes image numbers in EXCLUDE file")
    parser.add_option("--fftavg", type="string", help="Incoherent Fourier average of all images and write a single power spectrum image")
    parser.add_option("--filter", type="string", action="append", help="apply a filter. FILTER=filtername:param1=val1:param2=val2")
    parser.add_option("--first", type="int", help="the first image in the input to process [0 - n-1])")
    parser.add_option("--inplace", action="store_true", help="Output overwrites input, USE SAME FILENAME, DO NOT 'clip' images.")
    parser.add_option("--interlv", type="string", help="Specifies a 2nd input file. Output will be 2 files interleaved.")
    parser.add_option("--last", type="int", help="the last image in the input to process")
    parser.add_option("--list", type="string", help="Works only on the image numbers in LIST file")
    parser.add_option("--meanshrink", type="int", action="append", help="Reduce an image size by an integral scaling factor using average. Clip is not required.")
    parser.add_option("--mraprep",  action="store_true", help="this is an experimental option")
    parser.add_option("--norefs", action="store_true", help="Skip any input images which are marked as references (usually used with classes.*)")
    parser.add_option("--outtype", type="string", help="output image format, mrc, imagic, hdf, etc")
    parser.add_option("--plt", type="string", help="output the orientations in IMAGIC .plt file format")
    parser.add_option("--radon",  action="store_true", help="Do Radon transform")
    parser.add_option("--rfp",  action="store_true", help="this is an experimental option")
    parser.add_option("--scale", type="float", action="append", help="Scale by specified scaling factor. Clip must also be specified to change the dimensions of the output map.")
    parser.add_option("--selfcl", type="int", nargs=2, help="Output file will be a 180x180 self-common lines map for each image.")
    parser.add_option("--setsfpairs",  action="store_true", help="Applies the radial structure factor of the 1st image to the 2nd, the 3rd to the 4th, etc") 
    parser.add_option("--shrink", type="int", action="append", help="Reduce an image size by an integral scaling factor, uses median filter. Clip is not required.")
    parser.add_option("--split", type="int", help="Splits the input file into a set of n output files")
    parser.add_option("--verbose", type="int", help="verbose level [1-5]")

    append_options = ["clip", "filter", "meanshrink", "shrink", "scale"]

    optionlist = get_optionlist(sys.argv[1:])
    
    (options, args) = parser.parse_args()
    
    if len(args) != 2:
        print "usage: " + usage
        print "Please run '" + progname + " -h' for detailed options"
        sys.exit(1)
    
    infile = args[0]
    outfile = args[1]

    average = None
    fftavg = None
    
    n0 = 0
    n1 = -1

    if options.first:
        n0 = options.first
    if options.last:
        n1 = options.last

    sfout_n = 0
    sfout = None
    sf_amwid = 0

        
    MAXMICROCTF = 1000
    defocus_val = [0] * MAXMICROCTF
    bfactor_val = [0] * MAXMICROCTF

    pltfp = None
    if options.plt:
        pltfp = open(options.plt, "wb")
    
    if options.verbose:
        Log.set_log_level(options.verbose)
    
    d = EMData()
    nimg = EMUtil.get_image_count(infile)
    if nimg <= n1 or n1 < 0:
        n1 = nimg - 1

    ld = EMData()
    print "nimg = ", nimg

    imagelist = read_listfile(options.list, options.exclude, nimg)
    sfcurve1 = None

    index_d = {}
    for append_option in append_options:
        index_d[append_option] = 0
    
    for i in range(n0, n1+1):
        if imagelist and (not imagelist[i]):
            continue

        if options.split and options.split > 1:
            outfile = outfile[:-4] + ".%02d.img" % (i % split)

        d.read_image(infile, i)

        if pltfp:
            r = d.get_rotation()
            pi2d = 180/math.pi
            pltfp.write("%f,%f,%f \n" % (r.eman_phi() * pi2d,
                                         r.eman_alt() * pi2d,
                                         r.eman_az() * pi2d))
            continue
        
        nx = d.get_xsize()
        ny = d.get_ysize()
        
        sigma = d.get_attr("sigma").__float__()
        if sigma == 0:
            print "Warning: sigma = 0 for image " + i
            continue

        for option1 in optionlist:

            if option1 == "filter":
                fi = index_d[option1]
                (filtername, param_dict) = parse_filter_params(options.filter[fi])
                d.filter(filtername, param_dict)
                index_d[option1] += 1

            elif option1 == "norefs" and d.get_average_nimg() <= 0:
                continue
            
            elif option1 == "ctfsplit":

                if i == n0 or (not EMUtil.is_same_ctf(d, ld)):
                    ctf = d.get_ctf()
            
                    for j in range(1, options.ctfsplit):
                        if defocus_val[j] == ctf.get_defocus() and bfactor_val[j] == ctf.get_bfactor():
                            break
                    if options.ctfsplit <= j:
                        options.ctfsplit = j + 1
                        print "New CTF at " + i
                
                    defocus_val[j] = ctf.get_defocus()
                    bfactor_val[j] = ctf.get_bfactor()
                    outfile = outfile + ".%02d.img" % j
                    ld = d.copy(False, False)
            
                
            elif option1 == "setsfpairs":
                dataf = d.do_fft()
                d.gimme_fft()
                x0 = 0
                step = 0.5
                
                if i%2 == 0:
                    sfcurve1 = dataf.calc_radial_dist(nx, x0, step)
                else:
                    sfcurve2 = dataf.calc_radial_dist(nx, x0, step)
                    for j in range(nx):
                        if sfcurve1[j] > 0 and sfcurve2[j] > 0:
                            sfcurve2[j] = sqrt(sfcurve1[j] / sfcurve2[j])
                        else:
                            sfcurve2[j] = 0;

		    dataf.apply_radial_func(x0, step, sfcurve2);
                    d = dataf.do_ift();
                    dataf.gimme_fft();
                    
            elif option1 == "rfp":
                e = d.make_rotational_footprint()
                e.append_image("rfp.hed")

            elif option1 == "scale":
                scale_f = options.scale[index_d[option1]]
                if scale_f != 1.0:
                    d.scale(scale_f)
                index_d[option1] += 1
                 
            elif option1 == "clip":
                ci = index_d[option1]
                (clipx, clipy) = options.clip[ci]
				
                e = d.get_clip(Region((nx-clipx)/2, (ny-clipy)/2, clipx, clipy))
                e.set_attr("avgnimg", d.get_attr("avgnimg"))
                d = e
                index_d[option1] += 1
            
            elif option1 == "shrink":
                shrink_f = options.shrink[index_d[option1]]
                if shrink_f > 1:
                     d.median_shrink(shrink_f)
                index_d[option1] += 1

            elif option1 == "meanshrink":
                mshrink = options.meanshrink[index_d[option1]]
                if mshrink > 1:
                    d.mean_shrink(mshrink)
                index_d[option1] += 1
        
            elif option1 == "selfcl":
                scl = options.selfcl[0] / 2
                sclmd = options.selfcl[1]
                sc = EMData()
            
                if sclmd == 0:
                    sc.common_lines_real(d, d, scl, true)
                else:
                    e = d.copy()
                    e.filter("Phase180")
                
                    if sclmd == 1:
                        sc.common_lines(e, e, sclmd, scl, true)
                        sc.filter("LinearXform", Dict("shift", EMObject(-90.0), "scale", EMObject(-1.0)))		
                    elif sclmd == 2:
                        sc.common_lines(e, e, sclmd, scl, true)
                    else:
                        print "Error: invalid common-line mode '" + sclmd + "'"
                        sys.exit(1)
                
            elif option1 == "radon":
                r = d.do_radon()
                d = r
            
            elif option1 == "average":
                if not average:
                    average = d.copy(False, False)
                else:
                    average.add(d)
                continue

            elif option1 == "fftavg":
                if not fftavg:
                    fftavg = EMData()
                    fftavg.set_size(nx+2, ny)
                    fftavg.set_complex(1)
                    fftavg.to_zero()
                d.filter("EdgeMeanMask")
                d.filter("NormalizeStd")
                df = d.do_fft()
                df.mult(df.get_ysize())
                fftavg.add_incoherent(df)
                d.gimme_fft
                continue

            elif option1 == "calcsf":
                sfout_n = int(options.calcsf[0])
                sfout = options.calcsf[1]
                sf_amwid = 2 * math.pi / sfout_n
                    
                dataf = d.do_fft()
                d.gimme_fft()
                curve = dataf.calc_radial_dist(ny, 0, 0.5)
                outfile2 = sfout
                
                if n1 != 0:
                    outfile2 = sfout + ".%03d" % (i+100)
                    
                sf_dx = 1.0 / (apix * 2.0 * ny)
                Util.save_data(0, sf_dx, curve, outfile2)
	    
                if sfout_n > 0:
                    for j in range(0, sfout_n):
                        curve = dataf.calc_radial_dist(ny, 0, 0.5, j * sf_amwid, sf_amwid)                    
                        outfile2 = os.path.basename(sfout) + "-" + str(j) + "-" + str(sfout_n) + ".pwr"
                        if n1 != 0:
                            outfile2 = outfile2 + ".%03d" % (i+100)
                        Util.save_data(0, sf_dx, curve, outfile2)
                
                
            elif option1 == "interlv":
                d.read_image(options.interlv, i)
                d.append_image(outfile, IMAGIC)
            
            elif option1 == "outtype":
                if options.outtype in ["mrc", "pif", "png", "pgm"]:
                    if n1 != 0:
                        outfile = "%03d." % (i + 100) + outfile
                elif options.outtype == "spidersingle":
                    if n1 != 0:
                        spiderformat = "%s%%0%dd.spi" % (outfile, int(log10(i))+1)
                        outfile = spiderformat % i

                if options.inplace:
                        d.write_image(outfile, i)
                elif options.mraprep:
                        outfile = outfile + "%04d" % i + ".lst"
                        options.outtype = "lst"
                
                d.write_image(outfile, 0, EMUtil.get_image_ext_type(options.outtype))
                
    #end of image loop

    if average:
        average.set_average_nimg(n1-n0+1)
        average.filter("NormalizeStd");
	average.append_image(outfile);

    if options.fftavg:
        ffgavg.mult(1.0 / sqrt(n1 - n0 + 1))
        fftavg.write_image(options.fftavg, 0, MRC)
    
    n_outimg = EMUtil.get_image_count(args[1])
    print str(n_outimg) + " images"
    
    
if __name__ == "__main__":
    main()
    

