#!/usr/bin/env    python

from EMAN2 import *
import unittest,os,sys
from test import test_support
import testlib
from pyemtbx.exceptions import *

class TestAverager(unittest.TestCase):
    """averager test"""
    
    def test_ImageAverager(self):
        """test ImageAverager ..............................."""
        
   
def test_main():
    test_support.run_unittest(TestAverager)

if __name__ == '__main__':
    test_main()