

try:
  from pyStarDB import sp_pystardb as star
  STAR_AVAILABLE = True
except ImportError:
  STAR_AVAILABLE = False

import EMAN2
import EMAN2db
import EMAN2_cppwrap
from EMAN2 import Transform
from EMAN2 import EMUtil
from EMAN2 import EMAN2Ctf

import time
import unittest
import numpy  as np


inputfile = '/home/adnan/PycharmProjects/newrelion/Refine3D/job056/run_data.star'
outputfile = '/home/adnan/Desktop/star_file_project/testing.star'
bdbfile = 'bdb:/home/adnan/PycharmProjects/DoseWeighting/Relion_Stackv5/MotionCorr/job049/MOVIES_RELION/Particles/sphire_relion_stack'


class Test_sp_window_with_star(unittest.TestCase):

    def test_something(self):

        # In case of appending an existing star file
        local_bdb_stack = EMAN2db.db_open_dict(inputfile)
        test_dict = {}
        for i in range(29):
            test_dict[str(i)] = i

        newdataframe = local_bdb_stack.append(test_dict, ignore_index = True)

        # In case of creating a new one
        star_stack =  EMAN2db.db_open_dict(outputfile, write_new = True)
        star_stack['0'] = test_dict

        pass


    def test_new_clas(self):
        # aa = EMAN2db.db_open_dict(inputfile)
        # dd = star.StarFile(inputfile)
        # bb = EMAN2db.Pd_to_Db_conversion(dd)

        aa = EMAN2db.db_open_dict(inputfile)

        EMAN2db.db_close_dict(inputfile)

