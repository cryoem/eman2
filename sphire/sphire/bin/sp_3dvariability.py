#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
from past.utils import old_div
#
# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Tapu Shaikh 2019 (tapu.shaikh@mpi-dortmund.mpg.de)
# Author: Adnan Ali 2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Luca Lusnig 2019 (luca.lusnig@mpi-dortmund.mpg.de)
# Author: Toshio Moriya 2019 (toshio.moriya@kek.jp)
#
# Copyright (c) 2019 Max Planck Institute of Molecular Physiology
# Copyright (c) 2012 The University of Texas - Houston Medical School
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#
import EMAN2
import EMAN2_cppwrap
import EMAN2db
import mpi
import numpy
import optparse
from ..libpy import sp_applications
from ..libpy import sp_filter
from ..libpy import sp_fundamentals
from ..libpy import sp_global_def
from ..libpy import sp_logger
from ..libpy import sp_morphology
from ..libpy import sp_reconstruction
from ..libpy import sp_statistics
from ..libpy import sp_utilities
import sys
import time
from builtins import range
import os


mpi.mpi_init(0, [])

"""
Instruction:
nohup mpirun -np  64   --hostfile ./node4567.txt  sx3dvariability.py bdb:data/data \
  --output_dir=var3d  --window=300 --var3D=var.hdf --img_per_grp=100 \
      --CTF>var3d/printout &

 1. The order of applying the input parameters is:  1. window; 2. decimation. The low-pass 
    filter is with absolute frequency unit and applied to the decimated images. So, it is 
    equivalent to the low-pass filter of decimate*fl with respetive to the original image 
    size.

 2. Always use small decimation rate in the first run if there is no prior info about 3D \
   variability analysis of the data and the memory requirement. In addition, one can skip 
   the low-pass filtration when decimation rate is small, which significantly speeds up computation.

 3. The program will check the user provided parameters such that the final image size is a 
    product of small primes such as 2, 3, 5...

"""


def check_output_format(input_var, filename):
    out_data_3d = EMAN2.EMData(10, 10, 10)
    out_data_3d += 1
    try:
        out_data_3d.write_image("test_"+filename)
    except Exception as e:
        msg="The parameter '"+input_var+"' has a not valid extension or the file is not writable. Actual namefile is: '"+filename
        sp_global_def.ERROR(msg,action=1)
        return False

    # i cannot just remove it because in case of multiple cpu use it could crash
    try:
        os.remove("test_"+filename)
    except Exception as useless_e:
        pass
    return True



def run():
    def params_3D_2D_NEW(phi, theta, psi, s2x, s2y, mirror):
        # the final ali2d parameters already combine shifts operation first and rotation operation second for parameters converted from 3D
        if mirror:
            m = 1
            alpha, sx, sy, scalen = sp_utilities.compose_transform2(
                0, s2x, s2y, 1.0, 540.0 - psi, 0, 0, 1.0
            )
        else:
            m = 0
            alpha, sx, sy, scalen = sp_utilities.compose_transform2(
                0, s2x, s2y, 1.0, 360.0 - psi, 0, 0, 1.0
            )
        return alpha, sx, sy, m

    progname = optparse.os.path.basename(sys.argv[0])
    usage = (
        progname
        + " prj_stack  --ave2D= --var2D=  --ave3D= --var3D= --img_per_grp= --fl=  --aa=   --sym=symmetry --CTF"
    )
    parser = optparse.OptionParser(usage, version=sp_global_def.SPARXVERSION)

    parser.add_option(
        "--output_dir", type="string", default="./", help="Output directory"
    )
    parser.add_option(
        "--ave2D",
        type="string",
        default=False,
        help="Write to the disk a stack of 2D averages",
    )
    parser.add_option(
        "--var2D",
        type="string",
        default=False,
        help="Write to the disk a stack of 2D variances",
    )
    parser.add_option(
        "--ave3D",
        type="string",
        default=False,
        help="Write to the disk reconstructed 3D average",
    )
    parser.add_option(
        "--var3D",
        type="string",
        default=False,
        help="Compute 3D variability (time consuming!)",
    )
    parser.add_option(
        "--img_per_grp",
        type="int",
        default=100,
        help="Number of neighbouring projections.(Default is 100)",
    )
    parser.add_option(
        "--no_norm",
        action="store_true",
        default=False,
        help="Do not use normalization.(Default is to apply normalization)",
    )
    # parser.add_option("--radius", 	    type="int"         ,	default=-1   ,				help="radius for 3D variability" )
    parser.add_option(
        "--npad",
        type="int",
        default=2,
        help="Number of time to pad the original images.(Default is 2 times padding)",
    )
    parser.add_option(
        "--sym", type="string", default="c1", help="Symmetry. (Default is no symmetry)"
    )
    parser.add_option(
        "--fl",
        type="float",
        default=0.0,
        help="Low pass filter cutoff in absolute frequency (0.0 - 0.5) and is applied to decimated images. (Default - no filtration)",
    )
    parser.add_option(
        "--aa",
        type="float",
        default=0.02,
        help="Fall off of the filter. Use default value if user has no clue about falloff (Default value is 0.02)",
    )
    parser.add_option(
        "--CTF",
        action="store_true",
        default=False,
        help="Use CFT correction.(Default is no CTF correction)",
    )
    # parser.add_option("--MPI" , 		action="store_true",	default=False,				help="use MPI version")
    # parser.add_option("--radiuspca", 	type="int"         ,	default=-1   ,				help="radius for PCA" )
    # parser.add_option("--iter", 		type="int"         ,	default=40   ,				help="maximum number of iterations (stop criterion of reconstruction process)" )
    # parser.add_option("--abs", 		type="float"   ,        default=0.0  ,				help="minimum average absolute change of voxels' values (stop criterion of reconstruction process)" )
    # parser.add_option("--squ", 		type="float"   ,	    default=0.0  ,				help="minimum average squared change of voxels' values (stop criterion of reconstruction process)" )
    parser.add_option(
        "--VAR",
        action="store_true",
        default=False,
        help="Stack of input consists of 2D variances (Default False)",
    )
    parser.add_option(
        "--decimate",
        type="float",
        default=0.25,
        help="Image decimate rate, a number less than 1. (Default is 0.25)",
    )
    parser.add_option(
        "--window",
        type="int",
        default=0,
        help="Target image size relative to original image size. (Default value is zero.)",
    )
    # parser.add_option("--SND",			action="store_true",	default=False,				help="compute squared normalized differences (Default False)")
    # parser.add_option("--nvec",			type="int"         ,	default=0    ,				help="Number of eigenvectors, (Default = 0 meaning no PCA calculated)")
    parser.add_option(
        "--symmetrize",
        action="store_true",
        default=False,
        help="Prepare input stack for handling symmetry (Default False)",
    )
    parser.add_option(
        "--overhead", type="float", default=0.5, help="python overhead per CPU."
    )

    (options, args) = parser.parse_args()
    #####
    # from mpi import *

    #  This is code for handling symmetries by the above program.  To be incorporated. PAP 01/27/2015

    # Set up global variables related to bdb cache
    if sp_global_def.CACHE_DISABLE:
        sp_utilities.disable_bdb_cache()

    # Set up global variables related to ERROR function
    sp_global_def.BATCH = True

    # detect if program is running under MPI
    RUNNING_UNDER_MPI = "OMPI_COMM_WORLD_SIZE" in optparse.os.environ
    if RUNNING_UNDER_MPI:
        sp_global_def.MPI = True
    if options.output_dir == "./":
        current_output_dir = optparse.os.path.abspath(options.output_dir)
    else:
        current_output_dir = options.output_dir
    if options.symmetrize:

        if mpi.mpi_comm_size(mpi.MPI_COMM_WORLD) > 1:
            sp_global_def.ERROR("Cannot use more than one CPU for symmetry preparation")

        if not optparse.os.path.exists(current_output_dir):
            optparse.os.makedirs(current_output_dir)
            sp_global_def.write_command(current_output_dir)

        if optparse.os.path.exists(
            optparse.os.path.join(current_output_dir, "log.txt")
        ):
            optparse.os.remove(optparse.os.path.join(current_output_dir, "log.txt"))
        log_main = sp_logger.Logger(sp_logger.BaseLogger_Files())
        log_main.prefix = optparse.os.path.join(current_output_dir, "./")

        instack = args[0]
        sym = options.sym.lower()
        if sym == "c1":
            sp_global_def.ERROR("There is no need to symmetrize stack for C1 symmetry")

        line = ""
        for a in sys.argv:
            line += " " + a
        log_main.add(line)

        if instack[:4] != "bdb:":
            # if output_dir =="./": stack = "bdb:data"
            stack = "bdb:" + current_output_dir + "/data"
            sp_utilities.delete_bdb(stack)
            junk = sp_utilities.cmdexecute("sp_cpy.py  " + instack + "  " + stack)
        else:
            stack = instack

        qt = EMAN2_cppwrap.EMUtil.get_all_attributes(stack, "xform.projection")

        na = len(qt)
        ts = sp_utilities.get_symt(sym)
        ks = len(ts)
        angsa = [None] * na

        for k in range(ks):
            # Qfile = "Q%1d"%k
            # if options.output_dir!="./": Qfile = os.path.join(options.output_dir,"Q%1d"%k)
            Qfile = optparse.os.path.join(current_output_dir, "Q%1d" % k)
            # delete_bdb("bdb:Q%1d"%k)
            sp_utilities.delete_bdb("bdb:" + Qfile)
            # junk = cmdexecute("e2bdb.py  "+stack+"  --makevstack=bdb:Q%1d"%k)
            junk = sp_utilities.cmdexecute(
                "e2bdb.py  " + stack + "  --makevstack=bdb:" + Qfile
            )
            # DB = db_open_dict("bdb:Q%1d"%k)
            DB = EMAN2db.db_open_dict("bdb:" + Qfile)
            for i in range(na):
                ut = qt[i] * ts[k]
                DB.set_attr(i, "xform.projection", ut)
                # bt = ut.get_params("spider")
                # angsa[i] = [round(bt["phi"],3)%360.0, round(bt["theta"],3)%360.0, bt["psi"], -bt["tx"], -bt["ty"]]
            # write_text_row(angsa, 'ptsma%1d.txt'%k)
            # junk = cmdexecute("e2bdb.py  "+stack+"  --makevstack=bdb:Q%1d"%k)
            # junk = cmdexecute("sxheader.py  bdb:Q%1d  --params=xform.projection  --import=ptsma%1d.txt"%(k,k))
            DB.close()
        # if options.output_dir =="./": delete_bdb("bdb:sdata")
        sp_utilities.delete_bdb("bdb:" + current_output_dir + "/" + "sdata")
        # junk = cmdexecute("e2bdb.py . --makevstack=bdb:sdata --filt=Q")
        sdata = "bdb:" + current_output_dir + "/" + "sdata"
        sp_global_def.sxprint(sdata)
        junk = sp_utilities.cmdexecute(
            "e2bdb.py   " + current_output_dir + "  --makevstack=" + sdata + " --filt=Q"
        )
        # junk = cmdexecute("ls  EMAN2DB/sdata*")
        # a = get_im("bdb:sdata")
        a = sp_utilities.get_im(sdata)
        a.set_attr("variabilitysymmetry", sym)
        # a.write_image("bdb:sdata")
        a.write_image(sdata)

    else:

        myid = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
        number_of_proc = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
        main_node = 0
        shared_comm = mpi.mpi_comm_split_type(
            mpi.MPI_COMM_WORLD, mpi.MPI_COMM_TYPE_SHARED, 0, mpi.MPI_INFO_NULL
        )
        myid_on_node = mpi.mpi_comm_rank(shared_comm)
        no_of_processes_per_group = mpi.mpi_comm_size(shared_comm)
        masters_from_groups_vs_everything_else_comm = mpi.mpi_comm_split(
            mpi.MPI_COMM_WORLD, main_node == myid_on_node, myid_on_node
        )
        color, no_of_groups, balanced_processor_load_on_nodes = sp_utilities.get_colors_and_subsets(
            main_node,
            mpi.MPI_COMM_WORLD,
            myid,
            shared_comm,
            myid_on_node,
            masters_from_groups_vs_everything_else_comm,
        )
        overhead_loading = options.overhead * number_of_proc
        # memory_per_node  = options.memory_per_node
        # if memory_per_node == -1.: memory_per_node = 2.*no_of_processes_per_group
        keepgoing = 1

        current_window = options.window
        current_decimate = options.decimate

        if len(args) == 1:
            stack = args[0]
        else:
            sp_global_def.sxprint("Usage: " + usage)
            sp_global_def.sxprint(
                "Please run '" + progname + " -h' for detailed options"
            )
            sp_global_def.ERROR(
                "Invalid number of parameters used. Please see usage information above."
            )
            return

        t0 = time.time()
        # obsolete flags
        options.MPI = True
        # options.nvec = 0
        options.radiuspca = -1
        options.iter = 40
        options.abs = 0.0
        options.squ = 0.0

        if options.fl > 0.0 and options.aa == 0.0:
            sp_global_def.ERROR(
                "Fall off has to be given for the low-pass filter", myid=myid
            )

        # if options.VAR and options.SND:
        # 	ERROR( "Only one of var and SND can be set!",myid=myid )

        if options.VAR and (options.ave2D or options.ave3D or options.var2D):
            sp_global_def.ERROR(
                "When VAR is set, the program cannot output ave2D, ave3D or var2D",
                myid=myid,
            )

        # if options.SND and (options.ave2D or options.ave3D):
        # 	ERROR( "When SND is set, the program cannot output ave2D or ave3D", myid=myid )

        # if options.nvec > 0 :
        # 	ERROR( "PCA option not implemented", myid=myid )

        # if options.nvec > 0 and options.ave3D == None:
        # 	ERROR( "When doing PCA analysis, one must set ave3D", myid=myid )

        if current_decimate > 1.0 or current_decimate < 0.0:
            sp_global_def.ERROR(
                "Decimate rate should be a value between 0.0 and 1.0", myid=myid
            )

        if current_window < 0.0:
            sp_global_def.ERROR(
                "Target window size should be always larger than zero", myid=myid
            )

        if myid == main_node:
            img = sp_utilities.get_image(stack, 0)
            nx = img.get_xsize()
            ny = img.get_ysize()
            if min(nx, ny) < current_window:
                keepgoing = 0
        keepgoing = sp_utilities.bcast_number_to_all(
            keepgoing, main_node, mpi.MPI_COMM_WORLD
        )
        if keepgoing == 0:
            sp_global_def.ERROR(
                "The target window size cannot be larger than the size of decimated image",
                myid=myid,
            )

        options.sym = options.sym.lower()
        # if global_def.CACHE_DISABLE:
        # 	from utilities import disable_bdb_cache
        # 	disable_bdb_cache()
        # global_def.BATCH = True

        if myid == main_node:
            if not optparse.os.path.exists(current_output_dir):
                optparse.os.makedirs(
                    current_output_dir
                )  # Never delete output_dir in the program!

        img_per_grp = options.img_per_grp
        # nvec        = options.nvec
        radiuspca = options.radiuspca
        # if os.path.exists(os.path.join(options.output_dir, "log.txt")): os.remove(os.path.join(options.output_dir, "log.txt"))
        log_main = sp_logger.Logger(sp_logger.BaseLogger_Files())
        log_main.prefix = optparse.os.path.join(current_output_dir, "./")

        error = 0
        if myid == main_node:
            # check extension output files
            valid_output=True
            if options.var3D:
                valid_output = valid_output and check_output_format(input_var="--var3D", filename=options.var3D)

            if options.ave3D:
                valid_output = valid_output and check_output_format(input_var="--ave3D", filename=options.ave3D)

            if options.var2D:
                valid_output = valid_output and check_output_format(input_var="--var2D", filename=options.var2D)

            if options.ave2D:
                valid_output = valid_output and check_output_format(input_var="--ave2D", filename=options.ave2D)

            if valid_output is False:
                error = 1



            line = ""
            for a in sys.argv:
                line += " " + a
            log_main.add(line)
            log_main.add("-------->>>Settings given by all options<<<-------")
            log_main.add("Symmetry             : %s" % options.sym)
            log_main.add("Input stack          : %s" % stack)
            log_main.add("Output_dir           : %s" % current_output_dir)

            if options.ave3D:
                log_main.add("Ave3d                : %s" % options.ave3D)
            if options.var3D:
                log_main.add("Var3d                : %s" % options.var3D)
            if options.ave2D:
                log_main.add("Ave2D                : %s" % options.ave2D)
            if options.var2D:
                log_main.add("Var2D                : %s" % options.var2D)
            if options.VAR:
                log_main.add("VAR                  : True")
            else:
                log_main.add("VAR                  : False")
            if options.CTF:
                log_main.add("CTF correction       : True  ")
            else:
                log_main.add("CTF correction       : False ")

            log_main.add("Image per group      : %5d" % options.img_per_grp)
            log_main.add("Image decimate rate  : %4.3f" % current_decimate)
            log_main.add("Low pass filter      : %4.3f" % options.fl)
            current_fl = options.fl
            if current_fl == 0.0:
                current_fl = 0.5
            log_main.add(
                "Current low pass filter is equivalent to cutoff frequency %4.3f for original image size"
                % round((current_fl * current_decimate), 3)
            )
            log_main.add("Window size          : %5d " % current_window)
            log_main.add("sx3dvariability begins")

        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        error = sp_utilities.bcast_number_to_all(
            error, source_node=0, mpi_comm=mpi.MPI_COMM_WORLD
        )

        if error == 1:
            return

        symbaselen = 0
        if myid == main_node:
            nima = EMAN2_cppwrap.EMUtil.get_image_count(stack)
            img = sp_utilities.get_image(stack)
            nx = img.get_xsize()
            ny = img.get_ysize()
            nnxo = nx
            nnyo = ny
            if options.sym != "c1":
                imgdata = sp_utilities.get_im(stack)
                try:
                    i = imgdata.get_attr("variabilitysymmetry").lower()
                    if i != options.sym:
                        sp_global_def.ERROR(
                            "The symmetry provided does not agree with the symmetry of the input stack",
                            myid=myid,
                        )
                except:
                    sp_global_def.ERROR(
                        "Input stack is not prepared for symmetry, please follow instructions",
                        myid=myid,
                    )
                i = len(sp_utilities.get_symt(options.sym))
                if (old_div(nima, i)) * i != nima:
                    sp_global_def.ERROR(
                        "The length of the input stack is incorrect for symmetry processing",
                        myid=myid,
                    )
                symbaselen = old_div(nima, i)
            else:
                symbaselen = nima
        else:
            nima = 0
            nx = 0
            ny = 0
            nnxo = 0
            nnyo = 0
        nima = sp_utilities.bcast_number_to_all(nima)
        nx = sp_utilities.bcast_number_to_all(nx)
        ny = sp_utilities.bcast_number_to_all(ny)
        nnxo = sp_utilities.bcast_number_to_all(nnxo)
        nnyo = sp_utilities.bcast_number_to_all(nnyo)
        if current_window > max(nx, ny):
            sp_global_def.ERROR("Window size is larger than the original image size")

        if current_decimate == 1.0:
            if current_window != 0:
                nx = current_window
                ny = current_window
        else:
            if current_window == 0:
                nx = int(nx * current_decimate + 0.5)
                ny = int(ny * current_decimate + 0.5)
            else:
                nx = int(current_window * current_decimate + 0.5)
                ny = nx
        symbaselen = sp_utilities.bcast_number_to_all(symbaselen)

        # check FFT prime number
        is_fft_friendly = nx == sp_fundamentals.smallprime(nx)

        if not is_fft_friendly:
            if myid == main_node:
                log_main.add(
                    "The target image size is not a product of small prime numbers"
                )
                log_main.add("Program adjusts the input settings!")
            ### two cases
            if current_decimate == 1.0:
                nx = sp_fundamentals.smallprime(nx)
                ny = nx
                current_window = nx  # update
                if myid == main_node:
                    log_main.add("The window size is updated to %d." % current_window)
            else:
                if current_window == 0:
                    nx = sp_fundamentals.smallprime(int(nx * current_decimate + 0.5))
                    current_decimate = old_div(float(nx), nnxo)
                    ny = nx
                    if myid == main_node:
                        log_main.add(
                            "The decimate rate is updated to %f." % current_decimate
                        )
                else:
                    nx = sp_fundamentals.smallprime(
                        int(current_window * current_decimate + 0.5)
                    )
                    ny = nx
                    current_window = int(old_div(nx, current_decimate) + 0.5)
                    if myid == main_node:
                        log_main.add(
                            "The window size is updated to %d." % current_window
                        )

        if myid == main_node:
            log_main.add("The target image size is %d" % nx)

        if radiuspca == -1:
            radiuspca = old_div(nx, 2) - 2
        if myid == main_node:
            log_main.add("%-70s:  %d\n" % ("Number of projection", nima))
        img_begin, img_end = sp_applications.MPI_start_end(nima, number_of_proc, myid)

        """Multiline Comment0"""

        if options.VAR:  # 2D variance images have no shifts
            varList = []
            # varList   = EMData.read_images(stack, range(img_begin, img_end))
            for index_of_particle in range(img_begin, img_end):
                image = sp_utilities.get_im(stack, index_of_particle)
                if current_window > 0:
                    varList.append(
                        sp_fundamentals.fdecimate(
                            sp_fundamentals.window2d(
                                image, current_window, current_window
                            ),
                            nx,
                            ny,
                        )
                    )
                else:
                    varList.append(sp_fundamentals.fdecimate(image, nx, ny))

        else:
            if myid == main_node:
                t1 = time.time()
                proj_angles = []
                aveList = []
                tab = EMAN2_cppwrap.EMUtil.get_all_attributes(stack, "xform.projection")
                for i in range(nima):
                    t = tab[i].get_params("spider")
                    phi = t["phi"]
                    theta = t["theta"]
                    psi = t["psi"]
                    x = theta
                    if x > 90.0:
                        x = 180.0 - x
                    x = x * 10000 + psi
                    proj_angles.append([x, t["phi"], t["theta"], t["psi"], i])
                t2 = time.time()
                log_main.add(
                    "%-70s:  %d\n" % ("Number of neighboring projections", img_per_grp)
                )
                log_main.add("...... Finding neighboring projections\n")
                log_main.add("Number of images per group: %d" % img_per_grp)
                log_main.add("Now grouping projections")
                proj_angles.sort()
                proj_angles_list = numpy.full((nima, 4), 0.0, dtype=numpy.float32)
                for i in range(nima):
                    proj_angles_list[i][0] = proj_angles[i][1]
                    proj_angles_list[i][1] = proj_angles[i][2]
                    proj_angles_list[i][2] = proj_angles[i][3]
                    proj_angles_list[i][3] = proj_angles[i][4]
            else:
                proj_angles_list = 0
            proj_angles_list = sp_utilities.wrap_mpi_bcast(
                proj_angles_list, main_node, mpi.MPI_COMM_WORLD
            )
            proj_angles = []
            for i in range(nima):
                proj_angles.append(
                    [
                        proj_angles_list[i][0],
                        proj_angles_list[i][1],
                        proj_angles_list[i][2],
                        int(proj_angles_list[i][3]),
                    ]
                )
            del proj_angles_list
            proj_list, mirror_list = sp_utilities.nearest_proj(
                proj_angles, img_per_grp, range(img_begin, img_end)
            )
            all_proj = []
            for im in proj_list:
                for jm in im:
                    all_proj.append(proj_angles[jm][3])
            all_proj = list(set(all_proj))
            index = {}
            for i in range(len(all_proj)):
                index[all_proj[i]] = i
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            if myid == main_node:
                log_main.add(
                    "%-70s:  %.2f\n"
                    % ("Finding neighboring projections lasted [s]", time.time() - t2)
                )
                log_main.add(
                    "%-70s:  %d\n"
                    % ("Number of groups processed on the main node", len(proj_list))
                )
                log_main.add(
                    "Grouping projections took:  %12.1f [m]"
                    % (old_div((time.time() - t2), 60.0))
                )
                log_main.add("Number of groups on main node: ", len(proj_list))
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

            if myid == main_node:
                log_main.add("...... Calculating the stack of 2D variances \n")
            # Memory estimation. There are two memory consumption peaks
            # peak 1. Compute ave, var;
            # peak 2. Var volume reconstruction;
            # proj_params = [0.0]*(nima*5)
            aveList = []
            varList = []
            # if nvec > 0: eigList = [[] for i in range(nvec)]
            dnumber = len(all_proj)  # all neighborhood set for assigned to myid
            pnumber = len(proj_list) * 2.0 + img_per_grp  # aveList and varList
            tnumber = dnumber + pnumber
            vol_size2 = old_div(nx ** 3 * 4.0 * 8, 1.0e9)
            vol_size1 = old_div(2.0 * nnxo ** 3 * 4.0 * 8, 1.0e9)
            proj_size = old_div(
                nnxo * nnyo * len(proj_list) * 4.0 * 2.0, 1.0e9
            )  # both aveList and varList
            orig_data_size = old_div(nnxo * nnyo * 4.0 * tnumber, 1.0e9)
            reduced_data_size = old_div(nx * nx * 4.0 * tnumber, 1.0e9)
            full_data = numpy.full((number_of_proc, 2), -1.0, dtype=numpy.float16)
            full_data[myid] = orig_data_size, reduced_data_size
            if myid != main_node:
                sp_utilities.wrap_mpi_send(full_data, main_node, mpi.MPI_COMM_WORLD)
            if myid == main_node:
                dummy =None
                for iproc in range(number_of_proc):
                    if iproc != main_node:
                        dummy = sp_utilities.wrap_mpi_recv(iproc, mpi.MPI_COMM_WORLD)
                        full_data[numpy.where(dummy > -1)] = dummy[
                            numpy.where(dummy > -1)
                        ]
                del dummy
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            full_data = sp_utilities.wrap_mpi_bcast(
                full_data, main_node, mpi.MPI_COMM_WORLD
            )
            # find the CPU with heaviest load
            minindx = numpy.argsort(full_data, 0)
            heavy_load_myid = minindx[-1][1]
            total_mem = sum(full_data)
            if myid == main_node:
                if current_window == 0:
                    log_main.add(
                        "Nx:   current image size = %d. Decimated by %f from %d"
                        % (nx, current_decimate, nnxo)
                    )
                else:
                    log_main.add(
                        "Nx:   current image size = %d. Windowed to %d, and decimated by %f from %d"
                        % (nx, current_window, current_decimate, nnxo)
                    )
                log_main.add("Nproj:       number of particle images.")
                log_main.add("Navg:        number of 2D average images.")
                log_main.add("Nvar:        number of 2D variance images.")
                log_main.add(
                    "Img_per_grp: user defined image per group for averaging = %d"
                    % img_per_grp
                )
                log_main.add(
                    "Overhead:    total python overhead memory consumption   = %f"
                    % overhead_loading
                )
                log_main.add(
                    "Total memory) = 4.0*nx^2*(nproj + navg +nvar+ img_per_grp)/1.0e9 + overhead: %12.3f [GB]"
                    % (total_mem[1] + overhead_loading)
                )
            del full_data
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            if myid == heavy_load_myid:
                log_main.add(
                    "Begin reading and preprocessing images on processor. Wait... "
                )
                ttt = time.time()
            # imgdata = EMData.read_images(stack, all_proj)
            imgdata = [None for im in range(len(all_proj))]
            for index_of_proj in range(len(all_proj)):
                # image = get_im(stack, all_proj[index_of_proj])
                if current_window > 0:
                    imgdata[index_of_proj] = sp_fundamentals.fdecimate(
                        sp_fundamentals.window2d(
                            sp_utilities.get_im(stack, all_proj[index_of_proj]),
                            current_window,
                            current_window,
                        ),
                        nx,
                        ny,
                    )
                else:
                    imgdata[index_of_proj] = sp_fundamentals.fdecimate(
                        sp_utilities.get_im(stack, all_proj[index_of_proj]), nx, ny
                    )

                if current_decimate > 0.0 and options.CTF:
                    ctf = imgdata[index_of_proj].get_attr("ctf")
                    ctf.apix = old_div(ctf.apix, current_decimate)
                    imgdata[index_of_proj].set_attr("ctf", ctf)

                if myid == heavy_load_myid and index_of_proj % 100 == 0:
                    log_main.add(
                        " ...... %6.2f%% "
                        % (old_div(index_of_proj, float(len(all_proj))) * 100.0)
                    )
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            if myid == heavy_load_myid:
                log_main.add(
                    "All_proj preprocessing cost %7.2f m"
                    % (old_div((time.time() - ttt), 60.0))
                )
                log_main.add("Wait untill reading on all CPUs done...")
            """Multiline Comment1"""
            if not options.no_norm:
                mask = sp_utilities.model_circle(old_div(nx, 2) - 2, nx, nx)
            if myid == heavy_load_myid:
                log_main.add("Start computing 2D aveList and varList. Wait...")
                ttt = time.time()
            inner = old_div(nx, 2) - 4
            outer = inner + 2
            xform_proj_for_2D = [None for i in range(len(proj_list))]
            for i in range(len(proj_list)):
                ki = proj_angles[proj_list[i][0]][3]
                if ki >= symbaselen:
                    continue
                mi = index[ki]
                dpar = EMAN2_cppwrap.Util.get_transform_params(
                    imgdata[mi], "xform.projection", "spider"
                )
                phiM, thetaM, psiM, s2xM, s2yM = (
                    dpar["phi"],
                    dpar["theta"],
                    dpar["psi"],
                    -dpar["tx"] * current_decimate,
                    -dpar["ty"] * current_decimate,
                )
                grp_imgdata = []
                for j in range(img_per_grp):
                    mj = index[proj_angles[proj_list[i][j]][3]]
                    cpar = EMAN2_cppwrap.Util.get_transform_params(
                        imgdata[mj], "xform.projection", "spider"
                    )
                    alpha, sx, sy, mirror = params_3D_2D_NEW(
                        cpar["phi"],
                        cpar["theta"],
                        cpar["psi"],
                        -cpar["tx"] * current_decimate,
                        -cpar["ty"] * current_decimate,
                        mirror_list[i][j],
                        )
                    if thetaM <= 90:
                        if mirror == 0:
                            alpha, sx, sy, scale = sp_utilities.compose_transform2(
                                alpha, sx, sy, 1.0, phiM - cpar["phi"], 0.0, 0.0, 1.0
                            )
                        else:
                            alpha, sx, sy, scale = sp_utilities.compose_transform2(
                                alpha,
                                sx,
                                sy,
                                1.0,
                                180 - (phiM - cpar["phi"]),
                                0.0,
                                0.0,
                                1.0,
                                )
                    else:
                        if mirror == 0:
                            alpha, sx, sy, scale = sp_utilities.compose_transform2(
                                alpha, sx, sy, 1.0, -(phiM - cpar["phi"]), 0.0, 0.0, 1.0
                            )
                        else:
                            alpha, sx, sy, scale = sp_utilities.compose_transform2(
                                alpha,
                                sx,
                                sy,
                                1.0,
                                -(180 - (phiM - cpar["phi"])),
                                0.0,
                                0.0,
                                1.0,
                            )
                    imgdata[mj].set_attr(
                        "xform.align2d",
                        EMAN2_cppwrap.Transform(
                            {
                                "type": "2D",
                                "alpha": alpha,
                                "tx": sx,
                                "ty": sy,
                                "mirror": mirror,
                                "scale": 1.0,
                            }
                        ),
                    )
                    grp_imgdata.append(imgdata[mj])
                if not options.no_norm:
                    for k in range(img_per_grp):
                        ave, std, minn, maxx = EMAN2_cppwrap.Util.infomask(
                            grp_imgdata[k], mask, False
                        )
                        grp_imgdata[k] -= ave
                        grp_imgdata[k] = old_div(grp_imgdata[k], std)
                if options.fl > 0.0:
                    for k in range(img_per_grp):
                        grp_imgdata[k] = sp_filter.filt_tanl(
                            grp_imgdata[k], options.fl, options.aa
                        )

                #  Because of background issues, only linear option works.
                if options.CTF:
                    ave, var = sp_statistics.aves_wiener(
                        grp_imgdata, SNR=1.0e5, interpolation_method="linear"
                    )
                else:
                    ave, var = sp_statistics.ave_var(grp_imgdata)
                # Switch to std dev
                # threshold is not really needed,it is just in case due to numerical accuracy something turns out negative.
                var = sp_morphology.square_root(sp_morphology.threshold(var))

                sp_utilities.set_params_proj(ave, [phiM, thetaM, 0.0, 0.0, 0.0])
                sp_utilities.set_params_proj(var, [phiM, thetaM, 0.0, 0.0, 0.0])

                aveList.append(ave)
                varList.append(var)
                xform_proj_for_2D[i] = [phiM, thetaM, 0.0, 0.0, 0.0]

                """Multiline Comment2"""
                if (myid == heavy_load_myid) and (i % 100 == 0):
                    log_main.add(
                        " ......%6.2f%%  " % (old_div(i, float(len(proj_list))) * 100.0)
                    )
            del imgdata, grp_imgdata, cpar, dpar, all_proj, proj_angles, index
            if not options.no_norm:
                del mask
            if myid == main_node:
                del tab
            #  At this point, all averages and variances are computed
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

            if myid == heavy_load_myid:
                log_main.add(
                    "Computing aveList and varList took %12.1f [m]"
                    % (old_div((time.time() - ttt), 60.0))
                )

            xform_proj_for_2D = sp_utilities.wrap_mpi_gatherv(
                xform_proj_for_2D, main_node, mpi.MPI_COMM_WORLD
            )
            if myid == main_node:
                sp_utilities.write_text_row(
                    [str(entry) for entry in xform_proj_for_2D],
                    optparse.os.path.join(current_output_dir, "params.txt"),
                )
            del xform_proj_for_2D
            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            if options.ave2D:
                if myid == main_node:
                    log_main.add("Compute ave2D ... ")
                    km = 0
                    for i in range(number_of_proc):
                        if i == main_node:
                            for im in range(len(aveList)):
                                aveList[im].write_image(
                                    optparse.os.path.join(
                                        current_output_dir, options.ave2D
                                    ),
                                    km,
                                )
                                km += 1
                        else:
                            nl = mpi.mpi_recv(
                                1,
                                mpi.MPI_INT,
                                i,
                                sp_global_def.SPARX_MPI_TAG_UNIVERSAL,
                                mpi.MPI_COMM_WORLD,
                            )
                            nl = int(nl[0])
                            for im in range(nl):
                                ave = sp_utilities.recv_EMData(i, im + i + 70000)
                                """Multiline Comment3"""
                                tmpvol = sp_fundamentals.fpol(ave, nx, nx, 1)
                                tmpvol.write_image(
                                    optparse.os.path.join(
                                        current_output_dir, options.ave2D
                                    ),
                                    km,
                                )
                                km += 1
                else:
                    mpi.mpi_send(
                        len(aveList),
                        1,
                        mpi.MPI_INT,
                        main_node,
                        sp_global_def.SPARX_MPI_TAG_UNIVERSAL,
                        mpi.MPI_COMM_WORLD,
                    )
                    for im in range(len(aveList)):
                        sp_utilities.send_EMData(
                            aveList[im], main_node, im + myid + 70000
                        )
                        """Multiline Comment4"""
                if myid == main_node:
                    sp_applications.header(
                        optparse.os.path.join(current_output_dir, options.ave2D),
                        params="xform.projection",
                        fimport=optparse.os.path.join(current_output_dir, "params.txt"),
                    )
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            if options.ave3D:
                t5 = time.time()
                if myid == main_node:
                    log_main.add("Reconstruct ave3D ... ")
                ave3D = sp_reconstruction.recons3d_4nn_MPI(
                    myid, aveList, symmetry=options.sym, npad=options.npad
                )
                sp_utilities.bcast_EMData_to_all(ave3D, myid)
                if myid == main_node:
                    if current_decimate != 1.0:
                        ave3D = sp_fundamentals.resample(
                            ave3D, old_div(1.0, current_decimate)
                        )
                    ave3D = sp_fundamentals.fpol(
                        ave3D, nnxo, nnxo, nnxo
                    )  # always to the orignal image size
                    sp_utilities.set_pixel_size(ave3D, 1.0)
                    ave3D.write_image(
                        optparse.os.path.join(current_output_dir, options.ave3D)
                    )
                    log_main.add(
                        "Ave3D reconstruction took %12.1f [m]"
                        % (old_div((time.time() - t5), 60.0))
                    )
                    log_main.add(
                        "%-70s:  %s\n"
                        % ("The reconstructed ave3D is saved as ", options.ave3D)
                    )

            mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
            del ave, var, proj_list, stack, alpha, sx, sy, mirror, aveList
            """Multiline Comment5"""

            if options.ave3D:
                del ave3D
            if options.var2D:
                if myid == main_node:
                    log_main.add("Compute var2D...")
                    km = 0
                    for i in range(number_of_proc):
                        if i == main_node:
                            for im in range(len(varList)):
                                tmpvol = sp_fundamentals.fpol(varList[im], nx, nx, 1)
                                tmpvol.write_image(
                                    optparse.os.path.join(
                                        current_output_dir, options.var2D
                                    ),
                                    km,
                                )
                                km += 1
                        else:
                            nl = mpi.mpi_recv(
                                1,
                                mpi.MPI_INT,
                                i,
                                sp_global_def.SPARX_MPI_TAG_UNIVERSAL,
                                mpi.MPI_COMM_WORLD,
                            )
                            nl = int(nl[0])
                            for im in range(nl):
                                ave = sp_utilities.recv_EMData(i, im + i + 70000)
                                tmpvol = sp_fundamentals.fpol(ave, nx, nx, 1)
                                tmpvol.write_image(
                                    optparse.os.path.join(
                                        current_output_dir, options.var2D
                                    ),
                                    km,
                                )
                                km += 1
                else:
                    mpi.mpi_send(
                        len(varList),
                        1,
                        mpi.MPI_INT,
                        main_node,
                        sp_global_def.SPARX_MPI_TAG_UNIVERSAL,
                        mpi.MPI_COMM_WORLD,
                    )
                    for im in range(len(varList)):
                        sp_utilities.send_EMData(
                            varList[im], main_node, im + myid + 70000
                        )  # What with the attributes??
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
                if myid == main_node:
                    sp_applications.header(
                        optparse.os.path.join(current_output_dir, options.var2D),
                        params="xform.projection",
                        fimport=optparse.os.path.join(current_output_dir, "params.txt"),
                    )
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        if options.var3D:
            if myid == main_node:
                log_main.add("Reconstruct var3D ...")
            t6 = time.time()
            # radiusvar = options.radius
            # if( radiusvar < 0 ):  radiusvar = nx//2 -3
            res = sp_reconstruction.recons3d_4nn_MPI(
                myid, varList, symmetry=options.sym, npad=options.npad
            )
            # res = recons3d_em_MPI(varList, vol_stack, options.iter, radiusvar, options.abs, True, options.sym, options.squ)
            if myid == main_node:
                if current_decimate != 1.0:
                    res = sp_fundamentals.resample(res, old_div(1.0, current_decimate))
                res = sp_fundamentals.fpol(res, nnxo, nnxo, nnxo)
                sp_utilities.set_pixel_size(res, 1.0)
                res.write_image(
                    os.path.join(current_output_dir, options.var3D)
                )
                log_main.add(
                    "%-70s:  %s\n"
                    % ("The reconstructed var3D is saved as ", options.var3D)
                )
                log_main.add(
                    "Var3D reconstruction took %f12.1 [m]"
                    % (old_div((time.time() - t6), 60.0))
                )
                log_main.add(
                    "Total computation time %f12.1 [m]"
                    % (old_div((time.time() - t0), 60.0))
                )
                log_main.add("sx3dvariability finishes")

        if RUNNING_UNDER_MPI:
            sp_global_def.MPI = False

        sp_global_def.BATCH = False

def main():
    sp_global_def.print_timestamp("Start")
    run()
    sp_global_def.print_timestamp("Finish")
    mpi.mpi_finalize()

if __name__ == "__main__":
    main()
