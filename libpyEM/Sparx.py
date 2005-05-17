# Helper functions from Pawel Penczek
# Please do not alter this file without permision from the author.

from EMAN2 import *
from math import *
from time import *
from random import *

def getImage(imagename, nx = 0, ny = 1, nz = 1):
    """Read an image from the disk or assign existing object to the output.

    Usage: myimage = readImage("path/to/image")
    or     myimage = readImage(name_of_existing_image)
    """
    if type(imagename) == type(""):
        e = EMData()
        e.read_image(imagename)
    elif not imagename:
        e = EMData()
        if (nx > 0):
            e.set_size(nx, ny, nz)
    else:
        e = imagename
    return e

def dropImage(imagename,destination):
    """Write an image to the disk or assign to an output object.

    Usage:  dropImage(name_of_existing_image,"path/to/image")
    or      dropImage(name_of_existing_image,myimage)
    """
    if type(destination) == type(""):
        imagename.write_image(destination,0,SINGLE_SPIDER)
    else:
        destination = EMData()
        destination = imagename

def descriptive_statistics(image):
    """Calculate the descriptive statistics of an image.

    Usage: [mean, sigma, xmin, xmax, nx, ny, nz =] descriptive_statistics(image object)
           or
           [mean, sigma, xmin, xmax, nx, ny, nz =] descriptive_statistics(image filename)

    Purpose: calculate basic statistical characteristics of an image.
    """
    e = getImage(image)
    mean = e.get_attr("mean")
    sigma = e.get_attr("sigma")
    imin = e.get_attr("minimum")
    imax = e.get_attr("maximum")
    nx = e.get_xsize()
    ny = e.get_ysize()
    nz = e.get_zsize()
    print "Image size: nx = %i, ny = %i, nz = %i" % (nx, ny, nz)
    print "avg = %g, std dev = %g, min = %g, max = %g" % (mean, sigma, imin, imax)
    return mean,sigma,imin,imax, nx, ny, nz
    

def printRow(input, ix=0, iz=0):
    """Print the data in slice iz, row ix of an image to standard out.

    Usage: printRow(image, ix, iz)
       or
           printRow("path/to/image", ix, iz)
    """
    image=getImage(input)
    nx = image.get_xsize()
    ny = image.get_ysize()
    nz = image.get_zsize()
    print "(z = %d slice, x = %d row)" % (iz, ix)
    line = []
    for iy in xrange(ny):
        line.append("%12.5g  " % (image.get_value_at(ix,iy,iz)))
        if ((iy + 1) % 5 == 0):
            line.append("\n   ")
    line.append("\n")
    print "".join(line)

def printCol(input, iy=0, iz=0):
    """Print the data in slice iz, column iy of an image to standard out.

       Usage: printCol(image, iy, iz)
          or
              printCol("path/to/image", iy, iz)
    """
    image=getImage(input)
    nx = image.get_xsize()
    ny = image.get_ysize()
    nz = image.get_zsize()
    print "(z = %d slice, y = %d col)" % (iz, iy)
    line = []
    for ix in xrange(ny):
        line.append("%12.5g  " % (image.get_value_at(ix,iy,iz)))
        if ((ix + 1) % 5 == 0):
            line.append("\n   ")
    line.append("\n")
    print "".join(line)

def printSlice(input, iz=0):
    """Print the data in slice iz of an image to standard out.

    Usage: printImage(image, int)
       or
           printImage("path/to/image", int)
    """
    image=getImage(input)
    nx = image.get_xsize()
    ny = image.get_ysize()
    nz = image.get_zsize()
    print "(z = %d slice)" % (iz)
    line = []
    for ix in xrange(nx):
        for iy in xrange(ny):
            line.append("%12.5g  " % (image.get_value_at(ix,iy,iz)))
            if ((iy + 1) % 5 == 0):
                line.append("\n   ")
        line.append("\n")
    print "".join(line)

def printImage(input):
    """Print the data in an image to standard out.

    Usage: printImage(image)
       or
           printImage("path/to/image")
    """
    image=getImage(input)
    nz = image.get_zsize()
    for iz in xrange(nz):
        printSlice(input, iz)


def add_series(file_pattern,i1,i2,average,variance):
    """ Calculate average and variance files for an image series
    
    Usage:  add_series("img****.ext",i1,i2,average,variance)
      i1 - first file in image series
      i2 - last file in image series
      average and variance are output objects, or, if written as "a", are output disk files
      
    """
    fname = Util.parse_spider_fname(file_pattern,[i1]) #f=file_pattern[i1]
    ave = getImage(fname)
    var = ave*ave  #pow(ave,2.0)
    descriptive_statistics(ave)

    # process the remaining files
    for index in range(i1+1,i2+1):
        fname = Util.parse_spider_fname(file_pattern,[index])
        e = getImage(fname)
        ave = ave + e
        var = var + e*e  #pow(e,2.0)
    
    print "sum"
    descriptive_statistics(ave)
    ii=i2-i1+1
    ave = ave/ii
    print "average"
    descriptive_statistics(ave)
    #var = (var - pow(ave,2)/ii)/(ii-1)
    var = (var - ave*ave*ii)/(ii-1)
    print "variance"
    descriptive_statistics(var)

    dropImage(ave,average)
    dropImage(var,variance)


def do_reconstruction(filepattern, start, end, anglelist, symangs=[0.,0.,0.]):
    """Perform a 3-D reconstruction using Pawel's FFT Back Projection algoritm.
       
       Input:
         filepattern -- string such as "foo{****}.ext" that will be
                        used to determine the filenames of the 
                        projections to be read in.
        start        -- initial integer value to put in the field
        end          -- final integer value to put in the field
        anglelist    -- flat list of euler angles (in degrees), with the
                        number of euler angles equal to the number of 
                        projections to be read in.
        symangs      -- flat list of euler angles (in degrees)
                        corresponding to symmetries.
    
       Return:  3d reconstructed volume image

       Usage:
         
         anglelist = getAngles("myangles.txt") # not yet written
         symangs = getSymmetries("mysyms.txt") # not yet written
         filepattern = "proj{****}.hrs"
         start = 0
         end = 5087
         vol = do_reconstruction(filepattern, start, end, anglelist, symangs)
    """
    from math import radians
    npad = 4
    # convert angles to transform (rotation) objects
    nangles = len(anglelist) / 3
    rotations = []
    for i in range(nangles):
        phi = radians(anglelist[3*i])
        theta = radians(anglelist[3*i+1])
        psi = radians(anglelist[3*i+2])
        Ttype = Transform3D.EulerType.SPIDER
        rotations.append(Transform3D(Ttype, phi, theta, psi))
        
    # read first image to determine the size to use
    projname = Util.parse_spider_fname(filepattern,[start]) 
    first = getImage(projname)
    size = first.get_xsize()
    # sanity check -- image must be square
    if first.get_xsize() != first.get_ysize():
        print "Image projections must be square!"
        # FIXME: throw exception instead
        return None
    del first # don't need it any longer
    # reconstructor
    r = Reconstructors.get("PawelBackProjection", {"size":size, "npad":npad})
    r.setup()
    for i in range(start, end+1):
        projname = Util.parse_spider_fname(filepattern,[i])
        projection = getImage(projname)
        r.insert_slice(projection, rotations[i])
    v = r.finish()
    return v

def create_write_projections(volume, filepattern, anglelist, radius):
    nangles = len(anglelist) / 3
    for i in range(nangles):
        myangles = anglelist[3*i:3*(i+1)] # just a single slice of phi, theta, psi
        myparams = {"angletype":"SPIDER",
                    "anglelist":myangles,
                    "radius":radius}
        proj = volume.project("Pawel",myparams)
        projname = Util.parse_spider_fname(filepattern, [i])
        proj.write_image(projname, 0, EMUtil.ImageType.IMAGE_SINGLE_SPIDER)

def do_alignment(exptpattern, start, end, refpattern, alipattern, anglelist):
    newangles = []
    for i in range(start, end+1):
        exptname = Util.parse_spider_fname(exptpattern, [i])
        aliname  = Util.parse_spider_fname(alipattern, [i])
        exptimage = getImage(exptname)
        nangles = len(anglelist) / 3
        for ref in range(nangles):
            refname = Util.parse_spider_fname(refpattern, [ref])
            refimage = getImage(refname)
            #  do something real here
        # this next bit is utter rubbish just so the code "works"
        aliimage = exptimage 
        aliimage.write_image(aliname, 0, EMUtil.ImageType.IMAGE_SINGLE_SPIDER)
        newangles.append(1.0)
        newangles.append(2.0)
        newangles.append(3.0)
    return newangles

