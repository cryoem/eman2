#!/usr/bin/env python

from EMAN2 import *

def print_dict(d):
    print "    max = ", d["maximum"].get_float()
    print "    min = ", d["minimum"].get_float()
    print "average = ", d["average"].get_float()
    print
    
def main():
    e = EMData()
    e.read_image("/home/lpeng/images/tablet.mrc")
    d = e.get_attr_dict()
    print_dict(d)
    
    e += 10
    d1 = e.get_attr_dict()
    print_dict(d1)

    e += 100
    d2 = e.get_attr_dict()
    print_dict(d2)


if __name__ == "__main__":
    main()


