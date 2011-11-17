#!/usr/bin/env python

EMANVERSION="EMAN 2.0.4"
CVSDATESTAMP="$Date$"

def main():
    print EMANVERSION + ' (CVS' + CVSDATESTAMP[6:-2] +')' 

if __name__== "__main__":
    main()
