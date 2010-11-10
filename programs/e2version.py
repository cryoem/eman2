#!/usr/bin/env python

EMANVERSION="EMAN2 RC3"
CVSDATESTAMP="$Date$"

def main():
    print EMANVERSION + ' (CVS' + CVSDATESTAMP[6:-2] +')' 

if __name__== "__main__":
    main()
