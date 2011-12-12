#!/usr/bin/env python

import os
import sys

EMANVERSION="EMAN 2.0.4"
CVSDATESTAMP="$Date$"

def main():
    print EMANVERSION + ' (CVS' + CVSDATESTAMP[6:-2] +')' 
    
    if os.name=='posix':
    	cmd = 'cat /etc/system-release'
    	fin,fout = os.popen4(cmd)
    	result = fout.read().strip()
    	print 'Your EMAN2 is running on: ', result, os.uname()[2], os.uname()[-1]
    elif os.name=='nt':
		ver = sys.getwindowsversion()
		ver_format = ver[3], ver[0], ver[1]
		win_version = {
					(1, 4, 0): '95',
					(1, 4, 10): '98',
					(1, 4, 90): 'ME',
					(2, 4, 0): 'NT',
					(2, 5, 0): '2000',
					(2, 5, 1): 'XP',
					(2, 5, 2): '2003',
					(2, 6, 0): '2008',
					(2, 6, 1): '7'
				}
		
		if win_version.has_key(ver_format):
			winsysver = 'Windows' + ' ' + win_version[ver_format]
		else:
			winsysver = 'Windows'
    	
		if 'PROGRAMFILES(X86)' in os.environ:
			winsystype = '64bit'
		else:
			winsystype = '32bit'
    	
		print 'Your EMAN2 is running on: ', winsysver, winsystype
    	
    print 'Your Python version is: ', os.sys.version.split()[0]

if __name__== "__main__":
    main()
