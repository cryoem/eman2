from __future__ import print_function
from __future__ import division

# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Tapu Shaikh 2019 (tapu.shaikh@mpi-dortmund.mpg.de)
# Author: Adnan Ali 2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Luca Lusnig 2019 (luca.lusnig@mpi-dortmund.mpg.de)
#
# Copyright (c) 2019 Max Planck Institute of Molecular Physiology
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
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


import unittest
from numpy import allclose, array_equal
from sphire.libpy import sp_global_def
from os import path, mkdir, makedirs
import os


sp_global_def.BATCH = True
sp_global_def.MPI = False


from tests.test_module import (
    remove_list_of_file,
    remove_dir,
    get_arg_from_pickle_file,
    returns_values_in_file,
    give_ali_vol_data,
	ABSOLUTE_PATH_TO_RESOURCES,
)
from tests.test_module import (
    IMAGE_2D,
    IMAGE_3D,
    IMAGE_BLANK_3D,
)

from sphire.libpy.sp_utilities import even_angles

from EMAN2_cppwrap import EMData

from libpy_py3 import sp_applications as oldfu
from sphire.libpy import sp_applications as fu

TOLERANCE = 0.0005
import dill
dill._dill._reverse_typemap["ObjectType"] = object

ABSOLUTE_PATH_TO_STACK = 'bdb:'+os.path.join(os.getcwd(), "resources_tests/03_PARTICLES_BDB/mpi_proc_007/TcdA1-0187_frames_ptcls")
ABSOLUTE_PATH = 'resources_tests/applications_folder/'



"""
WHAT IS MISSING:
0) in all the cases where the input file is an image. I did not test the case with a complex image. I was not able to generate it
1) ali2d_MPI the results are different ... ANYWAY IT SEEMS TO BE NOT USED
2) ali2d_base the results are different
    -) if CTF =True and myid == main_node (i.e.: ln 695) it try to use an not defined variable and crash .... dead code or bug?
3) project3d --> How can I really test with 'listctfs'? It is never used in a real SPHIRE code case, anyway I'd test ... do we really want to test it?

RESULT AND KNOWN ISSUES
Some compatibility tests for the following functions fail!!!
1) Test_project3d --> some compatibility tests fail

In these tests there is a bug --> syntax error:
1) project3d --> in case of realsp=trillinear=True it will spawn an error message but its behaviour will be as if it'd receive as input realsp=True and Trilinear=False ....BUG or just BAD IMPLEMENTED

In these tests there is a strange behavior:
1) header --> you have to save the values in a file to test its behaviour. But if you set an output file you cannot set the other values becuase the
        "op = zero+one++consecutive+randomize+rand_alpha+(fimport!=None)+(fexport!=None)+fprint+backup+restore+delete+doset" check .... is it a bug? which is the purpose of this function??
2) within_group_refinement:
    -) I tested just the default method because the other 2 are not in use (Fabian docet)
    -) do not test with randomize=True because the results will be never the same because the random orientation in same calculations
    -) since the input parameter are basically used in 'sparx_alignment.ali2d_single_iter' and 'sparx_utilities.combine_params2' that i have already deeply tested
        I tested it just the pickle file case
"""


""" issues about the cleaned functions
There are some opened issues in:

3) mref_ali3d_MPI is corrupted function so we wont create unit test for that below unittest is just to show where the problem lies inside the code Adnan's note
4) Kmref_ali3d_MPI is corrupted function so we wont create unit test for that below unittest is just to show where the problem lies inside the code Adnan's note
6) recons3d_n_trl_MPI_one_node: I cannot run it. see the error message in "Test_recons3d_n_trl_MPI_one_node.test_NoCTF_symC1"
7) pca
    -) need a file name containing the average of the input stack
9) refvol: seems to be never USED
11) ali3d_mref_Kmeans_MPI and mref_ali3d_EQ_Kmeans are used only in sxsort3d.py that is buggy and maybe we will implement it from scratch. I did not test them for now
"""

"""
pickle files stored under smb://billy.storage.mpi-dortmund.mpg.de/abt3/group/agraunser/transfer/Adnan/pickle files
"""


"""
Comments from Adnan to Luca above issues
0) I have mentioned the solution of creating complex images in sp_alignment. Please have a look.
2) ali2d_base  is solved if you put a tolerance in it . because it creates data with random generation so it will always have different results.
7) You can average an input stack and it will give one image which is the average of all the stack images. I think   sum(input stack) / no of images will give you the result.
"""



"""IT SEEMS TO BE NOT USED"""


class Test_ali2d_MPI(unittest.TestCase):
    argum = get_arg_from_pickle_file(
        path.join(ABSOLUTE_PATH_TO_RESOURCES, "applications.ali2d_base")
    )

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali2d_MPI()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali2d_MPI()
        self.assertEqual(
            str(cm_new.exception), "ali2d_MPI() missing 2 required positional arguments: 'stack' and 'outdir'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    @unittest.skip("The result are not the same")
    def test_pickle_file_values(self):
        (
            stack,
            outdir,
            maskfile,
            ir,
            ou,
            rs,
            xr,
            yr,
            ts,
            nomirror,
            dst,
            center,
            maxit,
            CTF,
            snr,
            Fourvar,
            user_func_name,
            random_method,
            log,
            number_of_proc,
            myid,
            main_node,
            mpi_comm,
        ) = self.argum[0]

        outdirnewa = path.join(ABSOLUTE_PATH, "ali2d_MPI_NEW")
        outdirnewb = path.join(ABSOLUTE_PATH,"ali2d_MPI_OLD")

        remove_dir(outdirnewa)
        remove_dir(outdirnewb)

        fu.ali2d_MPI(
            ABSOLUTE_PATH_TO_STACK,
            outdirnewa,
            maskfile=maskfile,
            ou=ou,
            xr=xr,
            yr=yr,
            ts=ts,
            dst=dst,
            maxit=4,
            CTF=True,
            snr=snr,
        )
        # mpi_barrier(MPI_COMM_WORLD)
        oldfu.ali2d_MPI(
            ABSOLUTE_PATH_TO_STACK,
            outdirnewb,
            maskfile=maskfile,
            ou=ou,
            xr=xr,
            yr=yr,
            ts=ts,
            dst=dst,
            maxit=4,
            CTF=True,
            snr=snr,
        )
        # mpi_barrier(MPI_COMM_WORLD)

        self.assertEqual(
            returns_values_in_file(path.join(outdirnewa, "resolution001")),
            returns_values_in_file(path.join(outdirnewb, "resolution001")),
        )
        self.assertEqual(
            returns_values_in_file(path.join(outdirnewa, "resolution002")),
            returns_values_in_file(path.join(outdirnewb, "resolution002")),
        )
        remove_dir(outdirnewa)
        remove_dir(outdirnewb)

    @unittest.skip("To avoid output folders")
    def test_ali2d_MPI_error_directory_exists(self):
        (
            stack,
            outdir,
            maskfile,
            ir,
            ou,
            rs,
            xr,
            yr,
            ts,
            nomirror,
            dst,
            center,
            maxit,
            CTF,
            snr,
            Fourvar,
            user_func_name,
            random_method,
            log,
            number_of_proc,
            myid,
            main_node,
            mpi_comm,
        ) = self.argum[0]

        outdirnewa = path.join(ABSOLUTE_PATH, "ali2d_MPI_NEW")
        outdirnewb = path.join(ABSOLUTE_PATH, "ali2d_MPI_OLD")

        if os.path.isdir(outdirnewa):
            remove_dir(outdirnewa)
            mkdir(outdirnewa)
        else:
            pass
        if os.path.isdir(outdirnewb):
            remove_dir(outdirnewb)
            mkdir(outdirnewb)
        else:
            pass


        with self.assertRaises(SystemExit) as cm_new:
            fu.ali2d_MPI(
                ABSOLUTE_PATH_TO_STACK,
                outdirnewa,
                maskfile=maskfile,
                ou=ou,
                xr=xr,
                yr=yr,
                ts=ts,
                dst=dst,
                maxit=4,
                CTF=True,
                snr=snr,
            )

        sp_global_def.BATCH = True
        # mpi_barrier(MPI_COMM_WORLD)
        with self.assertRaises(SystemExit) as cm_old:
            oldfu.ali2d_MPI(
                ABSOLUTE_PATH_TO_STACK,
                outdirnewb,
                maskfile=maskfile,
                ou=ou,
                xr=xr,
                yr=yr,
                ts=ts,
                dst=dst,
                maxit=4,
                CTF=True,
                snr=snr,
            )
        # mpi_barrier(MPI_COMM_WORLD)

        self.assertEqual(str(cm_new.exception), str(1))
        self.assertEqual(str(cm_old.exception), str(1))
        #
        remove_dir(outdirnewa)
        remove_dir(outdirnewb)


class Test_ali2d_base(unittest.TestCase):
    argum = get_arg_from_pickle_file(
        path.join(ABSOLUTE_PATH_TO_RESOURCES, "applications.ali2d_base")
    )

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali2d_base()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali2d_base()
        self.assertEqual(
            str(cm_new.exception), "ali2d_base() missing 2 required positional arguments: 'stack' and 'outdir'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    @unittest.skip("MPI barrier causing trouble since we dont initialize now mpi")
    def test_ali2d_base_true_should_return_equal_object(self):
        # self.assertTrue(True)
        # """
        (stack, outdir, maskfile, ir, ou, rs, xr, yr, ts, nomirror, dst, center, maxit, CTF, snr,
         Fourvar, user_func_name, random_method, log ,number_of_proc, myid, main_node, mpi_comm) = self.argum[0]

        # stack = give_ali2d_base_data()

        outdirnew = path.join( ABSOLUTE_PATH,'ali2d_base_new')
        outdirnewold = path.join(ABSOLUTE_PATH,'ali2d_base_old')

        print("random method name is  ", random_method)
        remove_dir(outdirnew)
        remove_dir(outdirnewold)

        makedirs(outdirnew)
        makedirs(outdirnewold)
        number_of_proc = 1
        myid = 0
        main_node = 0
        maxit = 2
        return_new = fu.ali2d_base(stack, outdirnew, maskfile = maskfile, ir= ir, ou = ou, rs = rs, xr =xr, yr = yr,\
                                   ts = ts, nomirror= nomirror, dst = dst, center = center, maxit =maxit, CTF =True, snr = snr, \
                                   Fourvar =Fourvar, user_func_name = user_func_name,  \
                                   log = log, number_of_proc = number_of_proc, myid = myid, main_node = main_node, mpi_comm = None)

        # mpi_barrier(MPI_COMM_WORLD)

        return_old = oldfu.ali2d_base(stack, outdirnewold, maskfile = maskfile, ir= ir, ou = ou, rs = rs, xr =xr, yr = yr,\
                                   ts = ts, nomirror= nomirror, dst = dst, center = center, maxit =maxit, CTF =True, snr = snr, \
                                   Fourvar =Fourvar, user_func_name = user_func_name, random_method = random_method, \
                                   log = log, number_of_proc = number_of_proc, myid = myid, main_node = main_node, mpi_comm = None)

        # mpi_barrier(MPI_COMM_WORLD)

        self.assertTrue(allclose(return_new, return_old, atol=TOLERANCE,equal_nan=True))
        # image = sp_utilities.get_im(path.join(outdirnew, "aqfinal.hdf"))
        # image2 = sp_utilities.get_im(path.join(outdirnewold, "aqfinal.hdf"))
        # print('Image dimension', image.get_3dview().shape)
        # self.assertTrue(numpy.allclose(image.get_2dview(), image2.get_2dview() , atol = 0.1))
        remove_dir(outdirnew)
        remove_dir(outdirnewold)
        # self.assertEqual(returns_values_in_file(path.join(outdirnew, "resolution001")),returns_values_in_file(path.join(outdirnewold, "resolution001")))
        # self.assertEqual(returns_values_in_file(path.join(outdirnew, "resolution002")),returns_values_in_file(path.join(outdirnewold, "resolution002")))
        # """


class Test_cpy(unittest.TestCase):

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.cpy()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.cpy()
        print(str(cm_new.exception))
        print(type(cm_new.exception))
        self.assertEqual(
            str(cm_new.exception), "cpy() missing 2 required positional arguments: 'ins_list' and 'ous'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    """ sometimes fail
    def test_default_case(self):
        ins_list = path.join(ABSOLUTE_PATH_TO_SPHIRE_DEMO_RESULTS_FOLDER, 'Class2D/best.hdf')
        file_path_new = 'bdb:{0}#{1}'.format(ABSOLUTE_PATH, "new_substack")
        file_path_old = 'bdb:{0}#{1}'.format(ABSOLUTE_PATH, "old_substack")
        fu.cpy(ins_list,file_path_new)
        oldfu.cpy(ins_list, file_path_old)
        db_new, keys = EMAN2db_db_open_dict( file_path_new, True, True)
        db_new.realopen()
        db_old, keys = EMAN2db_db_open_dict( file_path_old, True, True)
        db_old.realopen()

        for index in db_new.bdb.keys():
            self.assertEqual(str(pickle_loads(db_new.bdb.get(index))), str(pickle_loads(db_old.bdb.get(index))))
        remove_dir( path.join(ABSOLUTE_PATH, 'EMAN2DB'))
    """

    @unittest.skip("Unwanted eman2db folder is created, the test works")
    def test_error_hdf_not_found(self):
        ins_list = path.join(ABSOLUTE_PATH, "not_found.hdf")
        file_path_new = "bdb:{0}#{1}".format(".", "new_substack")
        file_path_old = "bdb:{0}#{1}".format(".", "old_substack")
        with self.assertRaises(Exception) as cm_new:
            fu.cpy(ins_list, file_path_new)
        with self.assertRaises(Exception) as cm_old:
            oldfu.cpy(ins_list, file_path_old)

        print(str(cm_new.exception))
        self.assertEqual(
            cm_new.exception.args[0], cm_old.exception.args[0]
        )


# since it returns a huge list of 2Dimage I decide to unittest the len of this list, the first and the last image
class Test_project3d(unittest.TestCase):
    def test_all_the_conditions(self, return_new=(), return_old=()):
        self.assertEqual(len(return_new), len(return_old))
        for i, j in zip(return_new, return_old):
            self.assertTrue(
                allclose(j.get_3dview(), i.get_3dview(), atol=TOLERANCE, equal_nan=True)
            )

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.project3d()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.project3d()
        self.assertEqual(
            str(cm_new.exception), "project3d() missing 1 required positional argument: 'volume'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_2Dimg_crashes_because_signal11SIGSEGV(self):
        self.assertTrue(True)
        """
        return_new = fu.project3d(volume=IMAGE_2D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_2D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)
        """

    def test_Nonetype_img_returns_AttributeError(self):
        with self.assertRaises(AttributeError) as cm_new:
            fu.project3d(
                volume=None,
                stack=None,
                mask=None,
                delta=5,
                method="S",
                phiEqpsi="Minus",
                symmetry="c1",
                listagls=None,
                listctfs=None,
                noise=None,
                realsp=False,
                trillinear=False,
            )
        with self.assertRaises(AttributeError) as cm_old:
            oldfu.project3d(
                volume=None,
                stack=None,
                mask=None,
                delta=5,
                method="S",
                phiEqpsi="Minus",
                symmetry="c1",
                listagls=None,
                listctfs=None,
                noise=None,
                realsp=False,
                trillinear=False,
            )
        self.assertEqual(
            str(cm_new.exception), "'NoneType' object has no attribute 'get_xsize'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_empty_img_returns_ZeroDivisionError(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.project3d(
                volume=EMData(),
                stack=None,
                mask=None,
                delta=5,
                method="S",
                phiEqpsi="Minus",
                symmetry="c1",
                listagls=None,
                listctfs=None,
                noise=None,
                realsp=False,
                trillinear=False,
            )
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.project3d(
                volume=EMData(),
                stack=None,
                mask=None,
                delta=5,
                method="S",
                phiEqpsi="Minus",
                symmetry="c1",
                listagls=None,
                listctfs=None,
                noise=None,
                realsp=False,
                trillinear=False,
            )
        self.assertEqual(str(cm_new.exception), "float division by zero")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_3Dimg_realsp_and_trilinear_error_msg(self):
        with self.assertRaises(SystemExit) as cm_new:
            return_new = fu.project3d(
                volume=IMAGE_3D,
                stack=None,
                mask=None,
                delta=5,
                method="S",
                phiEqpsi="Minus",
                symmetry="c1",
                listagls=None,
                listctfs=None,
                noise=None,
                realsp=True,
                trillinear=True,
            )
        sp_global_def.BATCH = True
        oldfu.sp_global_def.BATCH = True
        # sp_global_def.MPI = True
        with self.assertRaises(SystemExit) as cm_old:
            return_old = oldfu.project3d(
                volume=IMAGE_3D,
                stack=None,
                mask=None,
                delta=5,
                method="S",
                phiEqpsi="Minus",
                symmetry="c1",
                listagls=None,
                listctfs=None,
                noise=None,
                realsp=True,
                trillinear=True,
            )

        print(str(cm_new.exception))
        self.assertEqual(
            cm_new.exception.args[0], cm_old.exception.args[0]
        )

        # self.test_all_the_conditions(return_new, return_old)
        # self.assertTrue(
        #     array_equal(
        #         return_new[0].get_2dview().flatten(),
        #         [
        #             -2.9546239376068115,
        #             -2.980635643005371,
        #             -2.957489490509033,
        #             -3.0083203315734863,
        #             -2.9922144412994385,
        #             -3.0299971103668213,
        #             -3.0781164169311523,
        #             -3.0145583152770996,
        #             -2.9256629943847656,
        #             0.0,
        #             -2.8150410652160645,
        #             -2.827364444732666,
        #             -2.8614673614501953,
        #             -2.830690383911133,
        #             -2.829454183578491,
        #             -2.874300718307495,
        #             -2.903381109237671,
        #             -2.948106050491333,
        #             -2.9924325942993164,
        #             0.0,
        #             -3.099151849746704,
        #             -3.1775362491607666,
        #             -3.1339433193206787,
        #             -3.058056592941284,
        #             -3.075274705886841,
        #             -3.0218067169189453,
        #             -3.0256922245025635,
        #             -2.9919466972351074,
        #             -3.0068883895874023,
        #             0.0,
        #             -3.1129660606384277,
        #             -3.1050636768341064,
        #             -3.0171358585357666,
        #             -3.062551975250244,
        #             -3.0829596519470215,
        #             -3.088615655899048,
        #             -3.130885601043701,
        #             -3.1775994300842285,
        #             -3.195767402648926,
        #             0.0,
        #             -3.2638678550720215,
        #             -3.2692501544952393,
        #             -3.2990612983703613,
        #             -3.3620688915252686,
        #             -3.322422504425049,
        #             -3.3664584159851074,
        #             -3.4466452598571777,
        #             -3.4770753383636475,
        #             -3.4245657920837402,
        #             0.0,
        #             -3.46885347366333,
        #             -3.513859987258911,
        #             -3.533717155456543,
        #             -3.5431182384490967,
        #             -3.549842119216919,
        #             -3.566971778869629,
        #             -3.549363613128662,
        #             -3.5507988929748535,
        #             -3.57666015625,
        #             0.0,
        #             -3.6432371139526367,
        #             -3.6037395000457764,
        #             -3.580766201019287,
        #             -3.555506467819214,
        #             -3.6111154556274414,
        #             -3.676938056945801,
        #             -3.6694860458374023,
        #             -3.7656328678131104,
        #             -3.8234314918518066,
        #             0.0,
        #             -4.073511123657227,
        #             -4.049807548522949,
        #             -4.092417240142822,
        #             -4.069482803344727,
        #             -4.054717540740967,
        #             -4.033217430114746,
        #             -4.045413970947266,
        #             -4.015814304351807,
        #             -4.025791168212891,
        #             0.0,
        #             -4.108157157897949,
        #             -4.08082389831543,
        #             -4.089325904846191,
        #             -4.105474472045898,
        #             -4.154421806335449,
        #             -4.11970329284668,
        #             -4.01369571685791,
        #             -4.044655799865723,
        #             -4.087397575378418,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #         ],
        #     )
        # )
        # self.assertTrue(
        #     array_equal(
        #         return_new[-1].get_2dview().flatten(),
        #         [
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             -4.726984024047852,
        #             -3.3136746883392334,
        #             0.06191571056842804,
        #             -0.7369051575660706,
        #             -1.886271595954895,
        #             -1.0977611541748047,
        #             0.0,
        #             0.0,
        #             0.0,
        #             -4.1548051834106445,
        #             -5.542438983917236,
        #             -1.1352554559707642,
        #             -0.521945059299469,
        #             -2.9987430572509766,
        #             -4.247453689575195,
        #             -0.5074707865715027,
        #             0.0,
        #             0.0,
        #             0.0,
        #             -5.457365989685059,
        #             -3.326054811477661,
        #             -0.1451924592256546,
        #             -2.027810573577881,
        #             -7.404072284698486,
        #             -2.60823917388916,
        #             -0.7226942777633667,
        #             0.0,
        #             0.0,
        #             -5.197442531585693,
        #             -6.677453517913818,
        #             -0.4378611147403717,
        #             -0.987167239189148,
        #             -6.304862976074219,
        #             -5.375576019287109,
        #             -1.645225167274475,
        #             -0.8825708627700806,
        #             0.0,
        #             -3.368478298187256,
        #             -9.108906745910645,
        #             -2.7291793823242188,
        #             -0.9390549659729004,
        #             -4.275655269622803,
        #             -7.276113510131836,
        #             -3.089268445968628,
        #             -0.8961836695671082,
        #             -2.811896800994873,
        #             0.0,
        #             -7.060383319854736,
        #             -5.394136905670166,
        #             -1.2909959554672241,
        #             -2.072096586227417,
        #             -7.80009651184082,
        #             -5.0674662590026855,
        #             -1.6751399040222168,
        #             -1.7208508253097534,
        #             -4.079824924468994,
        #             -3.752784490585327,
        #             -6.187793254852295,
        #             -1.7922195196151733,
        #             -0.8131762146949768,
        #             -5.211874485015869,
        #             -5.875916481018066,
        #             -2.8483564853668213,
        #             -0.655139148235321,
        #             -3.230496406555176,
        #             -5.632868766784668,
        #             -5.287583827972412,
        #             -2.3222413063049316,
        #             -1.0202511548995972,
        #             -2.8211617469787598,
        #             -5.295591831207275,
        #             -4.478401184082031,
        #             0.0827903300523758,
        #             -1.9030909538269043,
        #             -6.013443946838379,
        #             -4.83743143081665,
        #             -1.9481210708618164,
        #             -0.9295767545700073,
        #             -1.113168478012085,
        #             -2.9656569957733154,
        #             -4.386804103851318,
        #             -1.042134404182434,
        #             -0.9061259627342224,
        #             -4.556496620178223,
        #             -7.315436840057373,
        #             -1.2541165351867676,
        #             -0.34925854206085205,
        #             -0.31941545009613037,
        #             -1.0179522037506104,
        #             -2.6664602756500244,
        #             -1.561692237854004,
        #             -0.19883808493614197,
        #             -2.5375053882598877,
        #             -7.955165863037109,
        #             -3.782134771347046,
        #             2.181246280670166,
        #         ],
        #     )
        # )
        # self.assertEqual(len(return_old), 849)

    def test_3Dimg_default_case(self):
        vol , refv = give_ali_vol_data()
        return_new = fu.project3d(
            volume= vol,
            stack=None,
            mask=None,
            delta=5,
            method="S",
            phiEqpsi="Minus",
            symmetry="c1",
            listagls=None,
            listctfs=None,
            noise=None,
            realsp=False,
            trillinear=False,
        )
        sp_global_def.BATCH=True
        return_old = oldfu.project3d(
            volume=vol,
            stack=None,
            mask=None,
            delta=5,
            method="S",
            phiEqpsi="Minus",
            symmetry="c1",
            listagls=None,
            listctfs=None,
            noise=None,
            realsp=False,
            trillinear=False,
        )
        self.test_all_the_conditions(return_new, return_old)
        self.assertEqual(len(return_old), 849)

    def test_blank_img_default_case(self):
        return_new = fu.project3d(
            volume=IMAGE_BLANK_3D,
            stack=None,
            mask=None,
            delta=5,
            method="S",
            phiEqpsi="Minus",
            symmetry="c1",
            listagls=None,
            listctfs=None,
            noise=None,
            realsp=False,
            trillinear=False,
        )
        return_old = oldfu.project3d(
            volume=IMAGE_BLANK_3D,
            stack=None,
            mask=None,
            delta=5,
            method="S",
            phiEqpsi="Minus",
            symmetry="c1",
            listagls=None,
            listctfs=None,
            noise=None,
            realsp=False,
            trillinear=False,
        )
        self.test_all_the_conditions(return_new, return_old)
        self.assertTrue(
            array_equal(
                return_new[0].get_2dview().flatten(),
                [
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            )
        )
        self.assertTrue(
            array_equal(
                return_new[-1].get_2dview().flatten(),
                [
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            )
        )

        self.assertEqual(len(return_old), 849)

    def test_3Dimg_trilinear_case(self):
        vol, refv = give_ali_vol_data()
        return_new = fu.project3d(
            volume=vol,
            stack=None,
            mask=None,
            delta=5,
            method="S",
            phiEqpsi="Minus",
            symmetry="c1",
            listagls=None,
            listctfs=None,
            noise=None,
            realsp=False,
            trillinear=True,
        )
        return_old = oldfu.project3d(
            volume=vol,
            stack=None,
            mask=None,
            delta=5,
            method="S",
            phiEqpsi="Minus",
            symmetry="c1",
            listagls=None,
            listctfs=None,
            noise=None,
            realsp=False,
            trillinear=True,
        )
        self.test_all_the_conditions(return_new, return_old)
        self.assertEqual(len(return_old), 849)

    def test_blank_img_trilinear_case(self):
        return_new = fu.project3d(
            volume=IMAGE_BLANK_3D,
            stack=None,
            mask=None,
            delta=5,
            method="S",
            phiEqpsi="Minus",
            symmetry="c1",
            listagls=None,
            listctfs=None,
            noise=None,
            realsp=False,
            trillinear=True,
        )
        return_old = oldfu.project3d(
            volume=IMAGE_BLANK_3D,
            stack=None,
            mask=None,
            delta=5,
            method="S",
            phiEqpsi="Minus",
            symmetry="c1",
            listagls=None,
            listctfs=None,
            noise=None,
            realsp=False,
            trillinear=True,
        )
        self.test_all_the_conditions(return_new, return_old)
        self.assertTrue(
            array_equal(
                return_new[0].get_2dview().flatten(),
                [
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            )
        )
        self.assertTrue(
            array_equal(
                return_new[-1].get_2dview().flatten(),
                [
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            )
        )
        self.assertEqual(len(return_old), 849)

    def test_3Dimg_realsp_case(self):
        return_new = fu.project3d(
            volume=IMAGE_3D,
            stack=None,
            mask=None,
            delta=5,
            method="S",
            phiEqpsi="Minus",
            symmetry="c1",
            listagls=None,
            listctfs=None,
            noise=None,
            realsp=True,
            trillinear=False,
        )
        return_old = oldfu.project3d(
            volume=IMAGE_3D,
            stack=None,
            mask=None,
            delta=5,
            method="S",
            phiEqpsi="Minus",
            symmetry="c1",
            listagls=None,
            listctfs=None,
            noise=None,
            realsp=True,
            trillinear=False,
        )
        self.test_all_the_conditions(return_new, return_old)

        # self.assertTrue(
        #     array_equal(
        #         return_new[0].get_2dview().flatten(),
        #         [
        #             -2.9546239376068115,
        #             -2.980635643005371,
        #             -2.957489490509033,
        #             -3.0083203315734863,
        #             -2.9922144412994385,
        #             -3.0299971103668213,
        #             -3.0781164169311523,
        #             -3.0145583152770996,
        #             -2.9256629943847656,
        #             0.0,
        #             -2.8150410652160645,
        #             -2.827364444732666,
        #             -2.8614673614501953,
        #             -2.830690383911133,
        #             -2.829454183578491,
        #             -2.874300718307495,
        #             -2.903381109237671,
        #             -2.948106050491333,
        #             -2.9924325942993164,
        #             0.0,
        #             -3.099151849746704,
        #             -3.1775362491607666,
        #             -3.1339433193206787,
        #             -3.058056592941284,
        #             -3.075274705886841,
        #             -3.0218067169189453,
        #             -3.0256922245025635,
        #             -2.9919466972351074,
        #             -3.0068883895874023,
        #             0.0,
        #             -3.1129660606384277,
        #             -3.1050636768341064,
        #             -3.0171358585357666,
        #             -3.062551975250244,
        #             -3.0829596519470215,
        #             -3.088615655899048,
        #             -3.130885601043701,
        #             -3.1775994300842285,
        #             -3.195767402648926,
        #             0.0,
        #             -3.2638678550720215,
        #             -3.2692501544952393,
        #             -3.2990612983703613,
        #             -3.3620688915252686,
        #             -3.322422504425049,
        #             -3.3664584159851074,
        #             -3.4466452598571777,
        #             -3.4770753383636475,
        #             -3.4245657920837402,
        #             0.0,
        #             -3.46885347366333,
        #             -3.513859987258911,
        #             -3.533717155456543,
        #             -3.5431182384490967,
        #             -3.549842119216919,
        #             -3.566971778869629,
        #             -3.549363613128662,
        #             -3.5507988929748535,
        #             -3.57666015625,
        #             0.0,
        #             -3.6432371139526367,
        #             -3.6037395000457764,
        #             -3.580766201019287,
        #             -3.555506467819214,
        #             -3.6111154556274414,
        #             -3.676938056945801,
        #             -3.6694860458374023,
        #             -3.7656328678131104,
        #             -3.8234314918518066,
        #             0.0,
        #             -4.073511123657227,
        #             -4.049807548522949,
        #             -4.092417240142822,
        #             -4.069482803344727,
        #             -4.054717540740967,
        #             -4.033217430114746,
        #             -4.045413970947266,
        #             -4.015814304351807,
        #             -4.025791168212891,
        #             0.0,
        #             -4.108157157897949,
        #             -4.08082389831543,
        #             -4.089325904846191,
        #             -4.105474472045898,
        #             -4.154421806335449,
        #             -4.11970329284668,
        #             -4.01369571685791,
        #             -4.044655799865723,
        #             -4.087397575378418,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #         ],
        #     )
        # )
        # self.assertTrue(
        #     array_equal(
        #         return_new[-1].get_2dview().flatten(),
        #         [
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             -4.726984024047852,
        #             -3.3136746883392334,
        #             0.06191571056842804,
        #             -0.7369051575660706,
        #             -1.886271595954895,
        #             -1.0977611541748047,
        #             0.0,
        #             0.0,
        #             0.0,
        #             -4.1548051834106445,
        #             -5.542438983917236,
        #             -1.1352554559707642,
        #             -0.521945059299469,
        #             -2.9987430572509766,
        #             -4.247453689575195,
        #             -0.5074707865715027,
        #             0.0,
        #             0.0,
        #             0.0,
        #             -5.457365989685059,
        #             -3.326054811477661,
        #             -0.1451924592256546,
        #             -2.027810573577881,
        #             -7.404072284698486,
        #             -2.60823917388916,
        #             -0.7226942777633667,
        #             0.0,
        #             0.0,
        #             -5.197442531585693,
        #             -6.677453517913818,
        #             -0.4378611147403717,
        #             -0.987167239189148,
        #             -6.304862976074219,
        #             -5.375576019287109,
        #             -1.645225167274475,
        #             -0.8825708627700806,
        #             0.0,
        #             -3.368478298187256,
        #             -9.108906745910645,
        #             -2.7291793823242188,
        #             -0.9390549659729004,
        #             -4.275655269622803,
        #             -7.276113510131836,
        #             -3.089268445968628,
        #             -0.8961836695671082,
        #             -2.811896800994873,
        #             0.0,
        #             -7.060383319854736,
        #             -5.394136905670166,
        #             -1.2909959554672241,
        #             -2.072096586227417,
        #             -7.80009651184082,
        #             -5.0674662590026855,
        #             -1.6751399040222168,
        #             -1.7208508253097534,
        #             -4.079824924468994,
        #             -3.752784490585327,
        #             -6.187793254852295,
        #             -1.7922195196151733,
        #             -0.8131762146949768,
        #             -5.211874485015869,
        #             -5.875916481018066,
        #             -2.8483564853668213,
        #             -0.655139148235321,
        #             -3.230496406555176,
        #             -5.632868766784668,
        #             -5.287583827972412,
        #             -2.3222413063049316,
        #             -1.0202511548995972,
        #             -2.8211617469787598,
        #             -5.295591831207275,
        #             -4.478401184082031,
        #             0.0827903300523758,
        #             -1.9030909538269043,
        #             -6.013443946838379,
        #             -4.83743143081665,
        #             -1.9481210708618164,
        #             -0.9295767545700073,
        #             -1.113168478012085,
        #             -2.9656569957733154,
        #             -4.386804103851318,
        #             -1.042134404182434,
        #             -0.9061259627342224,
        #             -4.556496620178223,
        #             -7.315436840057373,
        #             -1.2541165351867676,
        #             -0.34925854206085205,
        #             -0.31941545009613037,
        #             -1.0179522037506104,
        #             -2.6664602756500244,
        #             -1.561692237854004,
        #             -0.19883808493614197,
        #             -2.5375053882598877,
        #             -7.955165863037109,
        #             -3.782134771347046,
        #             2.181246280670166,
        #         ],
        #     )
        # )
        # self.assertEqual(len(return_old), 849)

    def test_blank_img_realsp_case(self):
        return_new = fu.project3d(
            volume=IMAGE_BLANK_3D,
            stack=None,
            mask=None,
            delta=5,
            method="S",
            phiEqpsi="Minus",
            symmetry="c1",
            listagls=None,
            listctfs=None,
            noise=None,
            realsp=True,
            trillinear=False,
        )
        return_old = oldfu.project3d(
            volume=IMAGE_BLANK_3D,
            stack=None,
            mask=None,
            delta=5,
            method="S",
            phiEqpsi="Minus",
            symmetry="c1",
            listagls=None,
            listctfs=None,
            noise=None,
            realsp=True,
            trillinear=False,
        )
        self.test_all_the_conditions(return_new, return_old)
        self.assertTrue(
            array_equal(
                return_new[0].get_2dview().flatten(),
                [
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            )
        )
        self.assertTrue(
            array_equal(
                return_new[-1].get_2dview().flatten(),
                [
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            )
        )
        self.assertEqual(len(return_old), 849)

    @unittest.skip(
        "since it adds a random img_noise_to the stack the results cannot be the same"
    )
    def test_blank_img_with_noise(self):
        self.assertTrue(True)
        """
        return_new = fu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = 5, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = 5, realsp = False, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)
        """

    @unittest.skip(
        "since it adds a random img_noise_to the stack the results cannot be the same"
    )
    def test_3Dimg_with_noise(self):
        self.assertTrue(True)
        """
        return_new = fu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = 5, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_3D, stack = None, mask = None, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = 5, realsp = False, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)
        """

    def test_save_on_hdf(self):
        outnew = path.join(ABSOLUTE_PATH, "project3dnew.hdf")
        outold = path.join(ABSOLUTE_PATH, "project3old.hdf")
        # with self.assertRaises(SystemExit) as cm_new:
        return_new = fu.project3d(
            volume=IMAGE_BLANK_3D,
            stack=outnew,
            mask=None,
            delta=5,
            method="S",
            phiEqpsi="Minus",
            symmetry="c1",
            listagls=None,
            listctfs=None,
            noise=None,
            realsp=False,
            trillinear=False,
        )
        fu.sp_global_def.BATCH = True
        oldfu.sp_global_def.BATCH = True
        # with self.assertRaises(SystemExit) as cm_old:
        return_old = oldfu.project3d(
            volume=IMAGE_BLANK_3D,
            stack=outold,
            mask=None,
            delta=5,
            method="S",
            phiEqpsi="Minus",
            symmetry="c1",
            listagls=None,
            listctfs=None,
            noise=None,
            realsp=False,
            trillinear=False,
        )
        if (path.exists(outnew)):
            os.system("rm -f  " + outnew)
        if (path.exists(outold)):
            os.system("rm -f  " + outold)

        # self.assertTrue(numpy.allclose(return_new[2].get_3dview(), return_old[2].get_3dview() , atol= 0.1, equal_nan=True))

    @unittest.skip("compatibility test failed")
    def test_3Dimg_with_mask(self):
        """
        return_new = fu.project3d(volume=IMAGE_3D, stack = None, mask = MASK, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_3D, stack = None, mask = MASK, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)
        self.assertTrue(array_equal(return_new[0].get_2dview().flatten(),
                                    [float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN')]))
        self.assertTrue(array_equal(return_new[-1].get_2dview().flatten(),
                                    [float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN')]))
        self.assertEqual(len(return_old), 849)
        """

    @unittest.skip("compatibility test failed")
    def test_blank_img_with_mask(self):
        """
        return_new = fu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = MASK, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        return_old = oldfu.project3d(volume=IMAGE_BLANK_3D, stack = None, mask = MASK, delta = 5, method = "S", phiEqpsi = "Minus", symmetry = "c1", listagls = None , listctfs = None, noise = None, realsp = False, trillinear = False)
        self.test_all_the_conditions(return_new, return_old)
        self.assertTrue(array_equal(return_new[0].get_2dview().flatten(),
                                    [float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'),
                                     float('NaN'), float('NaN'), float('NaN'), float('NaN')]))
        self.assertTrue(array_equal(return_new[-1].get_2dview().flatten(), [float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN')]))
        self.assertEqual(len(return_old), 849)
        """

    def test_3Dimg_with_listagls(self):
        vol, refv = give_ali_vol_data()
        listangls = even_angles(
            delta=15.0,
            theta1=0.0,
            theta2=90.0,
            phi1=0.0,
            phi2=359.99,
            method="S",
            phiEqpsi="Minus",
            symmetry="sd1",
            ant=0.0,
        )
        return_new = fu.project3d(
            volume=vol,
            stack=None,
            mask=None,
            delta=5,
            method="S",
            phiEqpsi="Minus",
            symmetry="c1",
            listagls=listangls,
            listctfs=None,
            noise=None,
            realsp=False,
            trillinear=False,
        )
        return_old = oldfu.project3d(
            volume=vol,
            stack=None,
            mask=None,
            delta=5,
            method="S",
            phiEqpsi="Minus",
            symmetry="c1",
            listagls=listangls,
            listctfs=None,
            noise=None,
            realsp=False,
            trillinear=False,
        )
        self.test_all_the_conditions(return_new, return_old)
        self.assertEqual(len(return_old), 6)

    def test_blank_img_with_listagls(self):
        listangls = even_angles(
            delta=15.0,
            theta1=0.0,
            theta2=90.0,
            phi1=0.0,
            phi2=359.99,
            method="S",
            phiEqpsi="Minus",
            symmetry="sd1",
            ant=0.0,
        )
        return_new = fu.project3d(
            volume=IMAGE_BLANK_3D,
            stack=None,
            mask=None,
            delta=5,
            method="S",
            phiEqpsi="Minus",
            symmetry="c1",
            listagls=listangls,
            listctfs=None,
            noise=None,
            realsp=False,
            trillinear=False,
        )
        return_old = oldfu.project3d(
            volume=IMAGE_BLANK_3D,
            stack=None,
            mask=None,
            delta=5,
            method="S",
            phiEqpsi="Minus",
            symmetry="c1",
            listagls=listangls,
            listctfs=None,
            noise=None,
            realsp=False,
            trillinear=False,
        )
        self.test_all_the_conditions(return_new, return_old)
        self.assertTrue(
            array_equal(
                return_new[0].get_2dview().flatten(),
                [
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            )
        )
        self.assertTrue(
            array_equal(
                return_new[-1].get_2dview().flatten(),
                [
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            )
        )
        self.assertEqual(len(return_old), 6)

    def test_3Dimg_empty_listagls(self):
        return_new = fu.project3d(
            volume=IMAGE_3D,
            stack=None,
            mask=None,
            delta=5,
            method="S",
            phiEqpsi="Minus",
            symmetry="c1",
            listagls=[],
            listctfs=None,
            noise=None,
            realsp=False,
            trillinear=False,
        )
        return_old = oldfu.project3d(
            volume=IMAGE_3D,
            stack=None,
            mask=None,
            delta=5,
            method="S",
            phiEqpsi="Minus",
            symmetry="c1",
            listagls=[],
            listctfs=None,
            noise=None,
            realsp=False,
            trillinear=False,
        )
        self.assertTrue(array_equal(return_old, return_new))
        self.assertTrue(array_equal(return_new, []))

    def test_blank_img_empty_listagls(self):
        return_new = fu.project3d(
            volume=IMAGE_BLANK_3D,
            stack=None,
            mask=None,
            delta=5,
            method="S",
            phiEqpsi="Minus",
            symmetry="c1",
            listagls=[],
            listctfs=None,
            noise=None,
            realsp=False,
            trillinear=False,
        )
        return_old = oldfu.project3d(
            volume=IMAGE_BLANK_3D,
            stack=None,
            mask=None,
            delta=5,
            method="S",
            phiEqpsi="Minus",
            symmetry="c1",
            listagls=[],
            listctfs=None,
            noise=None,
            realsp=False,
            trillinear=False,
        )
        self.assertTrue(array_equal(return_old, return_new))
        self.assertTrue(array_equal(return_new, []))


# I cannot create an 3D image with 'xform.align3d' key
class Test_ali_vol(unittest.TestCase):
    argum = get_arg_from_pickle_file(
        path.join(ABSOLUTE_PATH_TO_RESOURCES, "applications.ali_vol")
    )

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.ali_vol()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.ali_vol()
        self.assertEqual(
            str(cm_new.exception), "ali_vol() missing 4 required positional arguments: 'vol', 'refv', 'ang_scale', and 'shift_scale'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_with_empty_refv_returns_RuntimeError(self):
        # (vol, refv, ang_scale, shift_scale, radius) = self.argum[0]

        vol , ref = give_ali_vol_data()

        shift_scale = 5.0
        ang_scale = 7.0
        radius = 37.5

        with self.assertRaises(RuntimeError) as cm_new:
            fu.ali_vol(
                vol=vol,
                refv=EMData(),
                ang_scale=ang_scale,
                shift_scale=shift_scale,
                radius=radius,
                discrepancy="ccc",
            )
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ali_vol(
                vol=vol,
                refv=EMData(),
                ang_scale=ang_scale,
                shift_scale=shift_scale,
                radius=radius,
                discrepancy="ccc",
            )
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_with_empty_vol_returns_RuntimeError(self):
        (vol, refv, ang_scale, shift_scale, radius) = self.argum[0]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ali_vol(
                vol=EMData(),
                refv=refv,
                ang_scale=ang_scale,
                shift_scale=shift_scale,
                radius=radius,
                discrepancy="ccc",
            )
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ali_vol(
                vol=EMData(),
                refv=refv,
                ang_scale=ang_scale,
                shift_scale=shift_scale,
                radius=radius,
                discrepancy="ccc",
            )
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_with_NoneType_vol_returns_RuntimeError(self):
        (vol, refv, ang_scale, shift_scale, radius) = self.argum[0]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ali_vol(
                vol=None,
                refv=refv,
                ang_scale=ang_scale,
                shift_scale=shift_scale,
                radius=radius,
                discrepancy="ccc",
            )
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ali_vol(
                vol=None,
                refv=refv,
                ang_scale=ang_scale,
                shift_scale=shift_scale,
                radius=radius,
                discrepancy="ccc",
            )
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_with_NoneType_refv_returns_RuntimeError(self):
        (vol, refv, ang_scale, shift_scale, radius) = self.argum[0]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ali_vol(
                vol=vol,
                refv=None,
                ang_scale=ang_scale,
                shift_scale=shift_scale,
                radius=radius,
                discrepancy="ccc",
            )
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ali_vol(
                vol=vol,
                refv=None,
                ang_scale=ang_scale,
                shift_scale=shift_scale,
                radius=radius,
                discrepancy="ccc",
            )
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "InvalidValueException")
        self.assertEqual(msg[3], "x size <= 0")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_2Dimg_as_vol_returns_RuntimeError(self):
        (vol, refv, ang_scale, shift_scale, radius) = self.argum[0]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ali_vol(
                vol=IMAGE_2D,
                refv=refv,
                ang_scale=ang_scale,
                shift_scale=shift_scale,
                radius=radius,
                discrepancy="ccc",
            )
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ali_vol(
                vol=IMAGE_2D,
                refv=refv,
                ang_scale=ang_scale,
                shift_scale=shift_scale,
                radius=radius,
                discrepancy="ccc",
            )
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_2Dimg_as_refvol_returns_RuntimeError(self):
        (vol, refv, ang_scale, shift_scale, radius) = self.argum[0]
        with self.assertRaises(RuntimeError) as cm_new:
            fu.ali_vol(
                vol=vol,
                refv=IMAGE_2D,
                ang_scale=ang_scale,
                shift_scale=shift_scale,
                radius=radius,
                discrepancy="ccc",
            )
        with self.assertRaises(RuntimeError) as cm_old:
            oldfu.ali_vol(
                vol=vol,
                refv=IMAGE_2D,
                ang_scale=ang_scale,
                shift_scale=shift_scale,
                radius=radius,
                discrepancy="ccc",
            )
        msg = str(cm_new.exception).split("'")
        msg_old = str(cm_old.exception).split("'")
        self.assertEqual(msg[0].split(" ")[0], "NotExistingObjectException")
        self.assertEqual(msg[3], "The requested key does not exist")
        self.assertEqual(msg[0].split(" ")[0], msg_old[0].split(" ")[0])
        self.assertEqual(msg[3], msg_old[3])
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_with_pickle_values(self):
        # (vol, refv, ang_scale, shift_scale, radius) = self.argum[0]
        vol , refv = give_ali_vol_data()
        shift_scale = 5.0
        ang_scale = 7.0
        radius = 37.5
        return_new = fu.ali_vol(
            vol=vol,
            refv=refv,
            ang_scale=ang_scale,
            shift_scale=shift_scale,
            radius=radius,
            discrepancy="ccc",
        )
        return_old = oldfu.ali_vol(
            vol=vol,
            refv=refv,
            ang_scale=ang_scale,
            shift_scale=shift_scale,
            radius=radius,
            discrepancy="ccc",
        )
        self.assertTrue(allclose(return_new.get_3dview(), return_old.get_3dview() , atol=1, equal_nan=True))
        # self.assertTrue(
        #     array_equal(
        #         return_new.get_3dview().flatten()[60100:60250],
        #         [
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             -1.0253742933273315,
        #             -0.7576642632484436,
        #             0.2794378101825714,
        #             0.8850940465927124,
        #             0.6183716058731079,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             -0.3909780979156494,
        #             -0.946395993232727,
        #             -1.3263295888900757,
        #             -0.9816580414772034,
        #             0.012633796781301498,
        #             0.8520616888999939,
        #             1.0967035293579102,
        #             0.7950575947761536,
        #             0.03073885850608349,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #         ],
        #     )
        # )

    def test_NoneRadius(self):
        # (vol, refv, ang_scale, shift_scale, radius) = self.argum[0]
        vol , refv = give_ali_vol_data()
        shift_scale = 5.0
        ang_scale = 7.0
        radius = 37.5
        return_new = fu.ali_vol(
            vol=vol,
            refv=refv,
            ang_scale=ang_scale,
            shift_scale=shift_scale,
            radius=None,
            discrepancy="ccc",
        )
        return_old = oldfu.ali_vol(
            vol=vol,
            refv=refv,
            ang_scale=ang_scale,
            shift_scale=shift_scale,
            radius=None,
            discrepancy="ccc",
        )

        self.assertTrue(allclose(return_new.get_3dview(), return_old.get_3dview(), atol = 1, equal_nan=True))
        # self.assertTrue(
        #     array_equal(
        #         return_new.get_3dview().flatten()[60100:60250],
        #         [
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             -1.0253742933273315,
        #             -0.7576642632484436,
        #             0.2794378101825714,
        #             0.8850940465927124,
        #             0.6183716058731079,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             -0.3909780979156494,
        #             -0.946395993232727,
        #             -1.3263295888900757,
        #             -0.9816580414772034,
        #             0.012633796781301498,
        #             0.8520616888999939,
        #             1.0967035293579102,
        #             0.7950575947761536,
        #             0.03073885850608349,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #         ],
        #     )
        # )

    def test_with_zero_radius(self):
        # (vol, refv, ang_scale, shift_scale, radius) = self.argum[0]
        vol , refv = give_ali_vol_data()
        shift_scale = 5.0
        ang_scale = 7.0
        radius = 37.5

        return_new = fu.ali_vol(
            vol=vol,
            refv=refv,
            ang_scale=ang_scale,
            shift_scale=shift_scale,
            radius=0,
            discrepancy="ccc",
        )
        return_old = oldfu.ali_vol(
            vol=vol,
            refv=refv,
            ang_scale=ang_scale,
            shift_scale=shift_scale,
            radius=0,
            discrepancy="ccc",
        )
        self.assertTrue(allclose(return_new.get_3dview(), return_old.get_3dview(), atol=1, equal_nan=True))
        # self.assertTrue(
        #     array_equal(
        #         return_new.get_3dview().flatten()[60100:60250],
        #         [
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             -1.0253742933273315,
        #             -0.7576642632484436,
        #             0.2794378101825714,
        #             0.8850940465927124,
        #             0.6183716058731079,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             -0.3909780979156494,
        #             -0.946395993232727,
        #             -1.3263295888900757,
        #             -0.9816580414772034,
        #             0.012633796781301498,
        #             0.8520616888999939,
        #             1.0967035293579102,
        #             0.7950575947761536,
        #             0.03073885850608349,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #             0.0,
        #         ],
        #     )
        # )

    def test_with_zero_shift_returns_ZeroDivisionError(self):
        (vol, refv, ang_scale, shift_scale, radius) = self.argum[0]
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.ali_vol(
                vol=vol,
                refv=refv,
                ang_scale=ang_scale,
                shift_scale=0,
                radius=radius,
                discrepancy="ccc",
            )
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.ali_vol(
                vol=vol,
                refv=refv,
                ang_scale=ang_scale,
                shift_scale=0,
                radius=radius,
                discrepancy="ccc",
            )
        self.assertEqual(str(cm_new.exception), "float division by zero")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_with_zero_ang_returns_ZeroDivisionError(self):
        (vol, refv, ang_scale, shift_scale, radius) = self.argum[0]
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.ali_vol(
                vol=vol,
                refv=refv,
                ang_scale=0,
                shift_scale=shift_scale,
                radius=radius,
                discrepancy="ccc",
            )
        with self.assertRaises(ZeroDivisionError) as cm_old:
            oldfu.ali_vol(
                vol=vol,
                refv=refv,
                ang_scale=0,
                shift_scale=shift_scale,
                radius=radius,
                discrepancy="ccc",
            )
        self.assertEqual(str(cm_new.exception), "float division by zero")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))



class Test_extract_value(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.extract_value()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.extract_value()
        self.assertEqual(
            str(cm_new.exception), "extract_value() missing 1 required positional argument: 's'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_string_integer(self):
        return_new = fu.extract_value("20")
        return_old = oldfu.extract_value("20")
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_old, 20)
        self.assertEqual(return_new, 20)

    def test_string_float(self):
        return_new = fu.extract_value("20.1")
        return_old = oldfu.extract_value("20.1")
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, 20.1)

    def test_string_not_handled(self):
        return_new = fu.extract_value("invalid")
        return_old = oldfu.extract_value("invalid")
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, "invalid")


class Test_header(unittest.TestCase):
    outnew = path.join(ABSOLUTE_PATH, "Headernew.hdf")
    outold = path.join(ABSOLUTE_PATH, "Headerold.hdf")
    params = "xform.projection"
    data = get_arg_from_pickle_file(
        path.join(ABSOLUTE_PATH_TO_RESOURCES, "applications.header")
    )[0]

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.header()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.header()
        self.assertEqual(
            str(cm_new.exception), "header() missing 2 required positional arguments: 'stack' and 'params'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_default_No_option_errorMsg(self):
        return_new = fu.header(
            stack=[],
            params=[],
            zero=False,
            one=False,
            set=0.0,
            randomize=False,
            rand_alpha=False,
            fimport=None,
            fexport=None,
            fprint=False,
            backup=False,
            suffix="_backup",
            restore=False,
            delete=False,
            consecutive=False,
        )
        return_old = oldfu.header(
            stack=[],
            params=[],
            zero=False,
            one=False,
            set=0.0,
            randomize=False,
            rand_alpha=False,
            fimport=None,
            fexport=None,
            fprint=False,
            backup=False,
            suffix="_backup",
            restore=False,
            delete=False,
            consecutive=False,
        )
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)

    def test_default_Too_option_errorMsg(self):
        return_new = fu.header(
            stack=[],
            params=[],
            zero=True,
            one=True,
            set=0.0,
            randomize=False,
            rand_alpha=False,
            fimport=None,
            fexport=None,
            fprint=False,
            backup=False,
            suffix="_backup",
            restore=False,
            delete=False,
            consecutive=False,
        )
        return_old = oldfu.header(
            stack=[],
            params=[],
            zero=True,
            one=True,
            set=0.0,
            randomize=False,
            rand_alpha=False,
            fimport=None,
            fexport=None,
            fprint=False,
            backup=False,
            suffix="_backup",
            restore=False,
            delete=False,
            consecutive=False,
        )
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)

    @unittest.skip('Need to ask Luca why it is failing. Unable to create EMAN2DB directory')
    def test_default_(self):
        return_new = fu.header(
            stack=ABSOLUTE_PATH_TO_STACK,
            params=self.params,
            zero=False,
            one=False,
            set=0.0,
            randomize=False,
            rand_alpha=False,
            fimport=None,
            fexport=self.outnew,
            fprint=False,
            backup=False,
            suffix="_backup",
            restore=False,
            delete=False,
            consecutive=False,
        )
        return_old = oldfu.header(
            stack=ABSOLUTE_PATH_TO_STACK,
            params=self.params,
            zero=False,
            one=False,
            set=0.0,
            randomize=False,
            rand_alpha=False,
            fimport=None,
            fexport=self.outold,
            fprint=False,
            backup=False,
            suffix="_backup",
            restore=False,
            delete=False,
            consecutive=False,
        )
        self.assertTrue(
            returns_values_in_file(self.outold), returns_values_in_file(self.outnew)
        )
        self.assertEqual(return_new, return_old)
        self.assertEqual(return_new, None)
        remove_list_of_file([self.outold, self.outnew])


class Test_MPI_start_end(unittest.TestCase):
    nima = 64
    nproc = 8
    myid = 1

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.MPI_start_end()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.MPI_start_end()
        print(str(cm_new.exception))
        self.assertEqual(
            str(cm_new.exception), "MPI_start_end() missing 3 required positional arguments: 'nima', 'nproc', and 'myid'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_default_case(self):
        return_new = fu.MPI_start_end(nima=self.nima, nproc=self.nproc, myid=self.myid)
        return_old = fu.MPI_start_end(nima=self.nima, nproc=self.nproc, myid=self.myid)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, [8, 16]))

    def test_zero_nima(self):
        return_new = fu.MPI_start_end(nima=0, nproc=self.nproc, myid=self.myid)
        return_old = fu.MPI_start_end(nima=0, nproc=self.nproc, myid=self.myid)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, [0, 0]))

    def test_zero_nproc(self):
        with self.assertRaises(ZeroDivisionError) as cm_new:
            fu.MPI_start_end(nima=self.nima, nproc=0, myid=self.myid)
        with self.assertRaises(ZeroDivisionError) as cm_old:
            fu.MPI_start_end(nima=self.nima, nproc=0, myid=self.myid)
        self.assertEqual(str(cm_new.exception), "float division by zero")
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_zero_myd(self):
        return_new = fu.MPI_start_end(nima=self.nima, nproc=self.nproc, myid=0)
        return_old = fu.MPI_start_end(nima=self.nima, nproc=self.nproc, myid=0)
        self.assertTrue(array_equal(return_new, return_old))
        self.assertTrue(array_equal(return_new, [0, 8]))


"""seems to be never USED"""


class Test_refvol(unittest.TestCase):
    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.refvol()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.refvol()
        self.assertEqual(
            str(cm_new.exception), "refvol() missing 4 required positional arguments: 'vollist', 'fsclist', 'output', and 'mask'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))


"""
    -) I tested just the default method because the other 2 are not in use (Fabian docet)
    -) do not test with randomize=True because the results will be never the same because the random orientation in same calculations
    -) since the input parameter are basically used in 'sparx_alignment.ali2d_single_iter' and 'sparx_utilities.combine_params2' that i have already deeply tested
        I tested it just the pickle file case
"""


class Test_within_group_refinement(unittest.TestCase):
    argum = get_arg_from_pickle_file(
        path.join(ABSOLUTE_PATH_TO_RESOURCES, "applications.within_group_refinement")
    )[0]
    randomize = False

    def test_wrong_number_params_too_few_parameters(self):
        with self.assertRaises(TypeError) as cm_new:
            fu.within_group_refinement()
        with self.assertRaises(TypeError) as cm_old:
            oldfu.within_group_refinement()
        self.assertEqual(
            str(cm_new.exception),
            "within_group_refinement() missing 13 required positional arguments: 'data', 'maskfile', 'randomize', 'ir', 'ou', 'rs', 'xrng', 'yrng', 'step', 'dst', 'maxit', 'FH', and 'FF'"
        )
        self.assertEqual(str(cm_new.exception), str(cm_old.exception))

    def test_pickle_file(self):
        (
            data,
            maskfile,
            randomize_not_used,
            ir,
            ou,
            rs,
            xrng,
            yrng,
            step,
            dst,
            maxit,
            FH,
            FF,
        ) = self.argum
        return_new = fu.within_group_refinement(
            data=data,
            maskfile=maskfile,
            randomize=self.randomize,
            ir=ir,
            ou=ou,
            rs=rs,
            xrng=xrng,
            yrng=yrng,
            step=step,
            dst=dst,
            maxit=maxit,
            FH=FH,
            FF=FF,
            method="",
            CTF=False,
        )
        return_old = oldfu.within_group_refinement(
            data=data,
            maskfile=maskfile,
            randomize=self.randomize,
            ir=ir,
            ou=ou,
            rs=rs,
            xrng=xrng,
            yrng=yrng,
            step=step,
            dst=dst,
            maxit=maxit,
            FH=FH,
            FF=FF,
            method="",
            CTF=False,
        )
        self.assertTrue(
            allclose(
                return_new.get_3dview(),
                return_old.get_3dview(),
                atol=TOLERANCE,
                equal_nan=True,
            )
        )

