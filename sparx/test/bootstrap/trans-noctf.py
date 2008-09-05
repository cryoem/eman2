from sparx import *
from EMAN2 import *
from sparx import *

from sys import  *
from string import *

i_ctfs = open( "CTF/ctfs.ttG" )
ignore = i_ctfs.readline()

pixel = 1.22
Cs = 2.0
voltage = 300.0
amp_contrast = 0.1

total_proj = 0
proj_out = "efg-all-noctf-part2.hdf"

for ii in xrange(47,61) :

    proj_in = "DATA/data-122-%3d.ttG.spi" % ii
    proj_in = replace( proj_in, ' ', '0' )

    fn_angl = "FINAL/defgrp%3d/angles171.ttG" % ii
    fn_angl = replace( fn_angl, ' ', '0' )

    fn_tran = "FINAL/defgrp%3d/trans171.ttG" % ii
    fn_tran = replace( fn_tran, ' ', '0' )

    nimage  = EMUtil.get_image_count( proj_in )

    f_trans = open( fn_tran, "r" )
    ignore = f_trans.readline()

    f_angle = open( fn_angl, "r" )
    ignore = f_angle.readline()

    defocus = atof( split( i_ctfs.readline() )[2] )

    print "Converting ", proj_in, ", defocus is ", defocus 

    for i in xrange(nimage):
        if( i%60 == 0 ) :
            stdout.write( " %6d  " % i )
            stdout.flush()

        trans_line = f_trans.readline()
        angle_line = f_angle.readline()

        data = EMData()
        data.read_image( proj_in, i )


        theta = atof( angle_line[16:27] )
        phi   = atof( angle_line[27:43] )

        psi = atof( trans_line[ 7:19] )
        sx  = atof( trans_line[19:31] ) * 2.0
        sy  = atof( trans_line[31:43] ) * 2.0
        mir = atof( trans_line[43:55] )

        nx = data.get_xsize()
        kb = kbt(nx)

        data = rotshift2dg(data, psi, sx, sy,kb)

        if mir > 0.5 : 
            data = data.process("mirror",{"axis":'x'})
    
        btwl = filt_btwl( data, 0.12, 0.13 )
        deci = Util.decimate( btwl, 4, 4, 1 )

        deci.set_attr_dict({'phi':phi, 'theta':theta, 'psi':0.0, 'sx':0.0, 'sy':0.0, 'mirror':0.0})
        deci.set_attr_dict({'defocus':defocus, 'amp_contrast':amp_contrast, 'voltage':voltage, 'Cs':Cs, 'pixel':pixel*4})
        deci.set_attr_dict({'b_factor':0.0})
        deci.set_attr_dict({'active':1, 'ctf_applied':0.0})

        deci.write_image( proj_out, total_proj )
        total_proj = total_proj + 1

        stdout.write( "." )
        stdout.flush( )

        if( i%60 == 59 ) :
            stdout.write( "\n" )

    stdout.write( "\n" )


print "Totally ", total_proj, " wrote"

