#!/usr/bin/env python

#
# Author: Grant Tang, 09/03/2005 (gtang@bcm.edu)
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
import unittest,os,sys
import testlib
from pyemtbx.exceptions import *
from math import pi
from optparse import OptionParser

IS_TEST_EXCEPTION = False

class TestAverager(unittest.TestCase):
    """averager test"""
	
    #def test_ImageAverager(self):
    #    """test ImageAverager ..............................."""
    
    def test_AbsMaxMinAverager(self):
		"""test AbsMaxMinAverager ..........................."""
		
		e = EMData(3,3)
		for i in range(3):
			for j in range(3):
				e.set_value_at(i,j,2.0)
		e.set_value_at(0,0,100)
		
		f = EMData(3,3)
		for i in range(3):
			for j in range(3):
				f.set_value_at(i,j,10.0)
		
		#test default find maximum values
		avgr = Averagers.get("absmaxmin")
		avgr.add_image(e)
		avgr.add_image(f)
		avg = avgr.finish()
		data = avg.get_2dview()
		self.assertAlmostEqual(data[0,0], 100.0, 3)
		self.assertAlmostEqual(data[2,2], 10.0, 3)
		
		#test finding minimum values
		avgr2 = Averagers.get("absmaxmin", {"min":1})
		avgr2.add_image(e)
		avgr2.add_image(f)
		avg2 = avgr2.finish()
		data2 = avg2.get_2dview()
		self.assertAlmostEqual(data2[0,0], 10.0, 3)
		self.assertAlmostEqual(data2[2,2], 2.0, 3)

def test_main():
    p = OptionParser()
    p.add_option('--t', action='store_true', help='test exception', default=False )
    global IS_TEST_EXCEPTION
    opt, args = p.parse_args()
    if opt.t:
        IS_TEST_EXCEPTION = True
    Log.logger().set_level(-1)  #perfect solution for quenching the Log error information, thank Liwei
    suite = unittest.TestLoader().loadTestsFromTestCase(TestAverager)
    unittest.TextTestRunner(verbosity=2).run(suite)

if __name__ == '__main__':
    test_main()
