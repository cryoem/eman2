from __future__ import print_function

import sys
import os
import clipboard

src = ' '.join(sys.argv[1:])

dest = src.split('/')
print("src: %s" % src)
print("dest: %s" % dest)
dest = "old_div({0}, {1})".format(dest[0], dest[1])

print(dest)
clipboard.copy(dest)
