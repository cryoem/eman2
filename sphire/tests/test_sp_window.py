

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
import pandas as pd
import os


inputfile = '/home/adnan/PycharmProjects/newrelion/Refine3D/job056/run_data.star'
outputfile = '/home/adnan/Desktop/star_file_project/testing.star'
bdbfile = 'bdb:/home/adnan/PycharmProjects/DoseWeighting/Relion_Stackv5/MotionCorr/job049/MOVIES_RELION/Particles/sphire_relion_stack'


class Test_sp_window_with_star(unittest.TestCase):

    def test_something(self):
        EMAN2db.db_get_all_attributes(inputfile, "ptcl_source_image")
        # In case of appending an existing star file
        local_bdb_stack = EMAN2db.db_open_dict(inputfile)
        test_dict = {}
        for i in range(29):
            test_dict[str(i)] = i

        # In case of creating a new one
        star_stack =  EMAN2db.db_open_dict(outputfile)
        star_stack[0] = test_dict

        pass


    def test_pandas_merging(self):
        fname = "name.star"
        fname1 = "name1.star"

        b = star.StarFile(fname)
        a = pd.DataFrame([[0, 1], [2, 3]], columns=['_c1', '_c2'])
        b.update('my_tag', a, True)

        c = star.StarFile(fname1)
        d = pd.DataFrame([[9, 3], [3, 6]], columns=['_c1', '_c2'])
        c.update('my_tag', d, True)

        newd =  b + c
        b += c
        print(newd)
        print(b)
        self.assertTrue(b , newd)


    def test_pandas_add_test(self):
        fname = "name.star"
        fname1 = "name1.star"

        b = star.StarFile(fname)
        a = pd.DataFrame([[0, 1], [2, 3]], columns=['_c1', '_c2'])
        b.update('my_tag', a, True)

        c = star.StarFile(fname1)
        d = pd.DataFrame([[9, 3], [3, 6]], columns=['_c1', '_c2'])
        c.update('my_tag', d, True)

        newd =  b + [c, c , c]
        b += [c, c , c]
        print(newd)
        print(b)
        self.assertTrue(b , newd)


    def test_markus_specification(self):

        fname = "name.star"
        fname1 = "name1.star"

        try:
            os.remove(fname)
            os.remove(fname1)
        except FileNotFoundError:
            pass
        x = pd.DataFrame([[0, 1], [1, 2]], columns=['_c1', '_c2'])
        y = pd.DataFrame([[3, 4], [5, 6]], columns=['_c1', '_c2'])
        a = star.StarFile(fname)
        a.update('a', x, True)
        b = star.StarFile(fname1)
        b.update('a', y, True)

        c = star.StarFile.add_star([a, b])

        c.line_dict['is_loop'] = True
        # c.write_star_file("testing.star")
        self.assertTrue(np.array_equal(c['a']['_c1'].values, [0, 1, 3, 5]) )
        self.assertTrue(np.array_equal(c['a']['_c2'].values, [1, 2, 4, 6]))

        data = {'_cc21': [3, 2, 1, 0], '_cc2': [2, 5, 7, 9]}

        d = star.StarFile(None)
        star.StarFile.add_star([a, b], d)
        self.assertTrue(np.array_equal(d['a']['_c1'].values, [0, 1, 3, 5]) )
        self.assertTrue(np.array_equal(d['a']['_c2'].values, [1, 2, 4, 6]))

        e = a + b
        self.assertTrue(np.array_equal(e['a']['_c1'].values, [0, 1, 3, 5]) )
        self.assertTrue(np.array_equal(e['a']['_c2'].values, [1, 2, 4, 6]))

        f = a
        f += b
        self.assertTrue(np.array_equal(f['a']['_c1'].values, [0, 1, 3, 5]))
        self.assertTrue(np.array_equal(f['a']['_c2'].values, [1, 2, 4, 6]))


    # def check_get_index(self):
    #     from sp_star import

























