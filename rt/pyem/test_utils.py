#!/bin/env python

from EMAN2 import *
import unittest
from test import test_support
import os


class TestUtils(unittest.TestCase):

    def test_is_file_exist(self):
        result1 = Util.is_file_exist(TestUtil.get_debug_image("search.dm3"))
        self.assertEqual(result1, True)
        
        result2 = Util.is_file_exist(TestUtil.get_debug_image("__nosuchafile__.dm3"))
        self.assertEqual(result2, False)
        
        result3 = Util.is_file_exist("")
        self.assertEqual(result3, False)


    def test_sstrcmp(self):
        e1 = Util.sstrncmp("helloworld", "hello");
        self.assertEqual(e1, True);

        e2 = Util.sstrncmp("foobar", "bar");
        self.assertEqual(e2, False);

        e3 = Util.sstrncmp("", "bar");
        self.assertEqual(e3, False);

        e4 = Util.sstrncmp("bar", "");
        self.assertEqual(e4, True);

        e5 = Util.sstrncmp("", "");
        self.assertEqual(e5, True);

    def test_int2str(self):
        s1 = Util.int2str(123)
        self.assertEqual(s1, "123")

        s2 = Util.int2str(-1)
        self.assertEqual(s2, "-1")

        
    def test_change_filename_ext(self):
        file1 = Util.change_filename_ext("hello.c", "cpp")
        file2 = Util.change_filename_ext("hello.", "mrc")
        file3 = Util.change_filename_ext("hello", "cpp")

        file4 = Util.change_filename_ext("hello.c", "")
        file5 = Util.change_filename_ext("hello.", "")
        file6 = Util.change_filename_ext("hello", "")

        self.assertEqual(file1, "hello.cpp")
        self.assertEqual(file2, "hello.mrc")
        self.assertEqual(file3, "hello.cpp")
        self.assertEqual(file4, "hello.c")
        self.assertEqual(file5, "hello.")
        self.assertEqual(file6, "hello")


    def test_remove_filename_ext(self):
        s1 = Util.remove_filename_ext("hello.cpp")
        self.assertEqual(s1, "hello")

        s1 = Util.remove_filename_ext("hello.")
        self.assertEqual(s1, "hello")

        s1 = Util.remove_filename_ext("hello")
        self.assertEqual(s1, "hello")

        s1 = Util.remove_filename_ext("")
        self.assertEqual(s1, "")
    

    def test_get_filename_ext(self):
        s1 = Util.get_filename_ext("hello.cpp")
        self.assertEqual(s1, "cpp")
        
        s1 = Util.get_filename_ext("hello.")
        self.assertEqual(s1, "")
        
        s1 = Util.get_filename_ext("hello")
        self.assertEqual(s1, "")
        
        s1 = Util.get_filename_ext("")
        self.assertEqual(s1, "")	
    

    def test_sbasename(self):
        b1 = Util.sbasename("hello.c")
        self.assertEqual(b1, "hello.c")

        b2 = Util.sbasename("/tmp/hello.mrc")
        self.assertEqual(b2, "hello.mrc")

        b3 = Util.sbasename("./test/hello")
        self.assertEqual(b3, "hello")


    
        

    

        
def test_main():
    test_support.run_unittest(TestUtils)

if __name__ == '__main__':
    test_main()


