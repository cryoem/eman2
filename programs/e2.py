#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

import os

e2real = os.path.join(os.path.abspath(os.path.dirname(__file__)), "e2_real.py")
os.execlp("ipython","ipython","-i","--gui=qt",e2real)
