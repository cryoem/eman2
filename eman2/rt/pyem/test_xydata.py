#!/usr/bin/env python

#
# Author: Liwei Peng, 01/30/2005 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

from EMAN2 import *
import unittest
import testlib
from optparse import OptionParser

IS_TEST_EXCEPTION = False

class TestXYData(unittest.TestCase):
    """this is the unit test for XYData class"""
    
    def test_pair(self):
        """test XYData.Pair class ..........................."""
        p = XYData.Pair(1.1, 2.2)
        self.assertAlmostEqual(1.1, p.x, 3)
        self.assertAlmostEqual(2.2, p.y, 3)
        
        p2 = XYData.Pair(1.0, 2.0)
        self.assertTrue(p2<p)
    
    def test_get_size(self):
        """test get_size() function ........................."""
        xy = XYData()
        xy.set_size(100)
        self.assertEqual(100, xy.get_size())
        
    def test_get_x(self):
        """test get_x() function ............................"""
        xy = XYData()
        lx = range(10)
        ly = (i+100 for i in range(10))
        xy.set_xy_list(lx, ly)
        
        xy.update();
        self.assertAlmostEqual(5.0, xy.get_x(5), 3)
        self.assertAlmostEqual(105.0, xy.get_y(5), 3)
        self.assertAlmostEqual(106.0, xy.get_yatx(6), 3)
        
    def test_get_xlist(self):
        """test get_xlist() function ........................"""
        xy = XYData()
        lx = range(10)
        ly = (i+100 for i in range(10))
        xy.set_xy_list(lx, ly)
        
        xlist = xy.get_xlist()
        ylist = xy.get_ylist()
        for i in xrange(len(xlist)):
            self.assertAlmostEqual(i, xlist[i], 3)
            self.assertAlmostEqual(100+i, ylist[i], 3)
        
    def test_read_write_file(self):
        """test read/write file ............................."""
        xy = XYData()
        lx = range(10)
        ly = (i+100 for i in range(10))
        xy.set_xy_list(lx, ly)
        
        file = 'xydata.dat'
        xy.write_file(file)
        xy2 = XYData()
        xy2.read_file(file)
        self.assertAlmostEqual(5.0, xy2.get_x(5), 3)
        self.assertAlmostEqual(105.0, xy2.get_y(5), 3)
        self.assertAlmostEqual(106.0, xy2.get_yatx(6), 3)
        testlib.safe_unlink(file)
        
def test_main():
    p = OptionParser()
    p.add_option('--t', action='store_true', help='test exception', default=False )
    global IS_TEST_EXCEPTION
    opt, args = p.parse_args()
    if opt.t:
        IS_TEST_EXCEPTION = True
    Log.logger().set_level(-1)  #perfect solution for quenching the Log error information, thank Liwei
    suite = unittest.TestLoader().loadTestsFromTestCase(TestXYData)
    unittest.TextTestRunner(verbosity=2).run(suite)

if __name__ == '__main__':
    test_main()
