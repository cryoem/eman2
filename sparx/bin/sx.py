#!/usr/bin/env python

import os

e2real = os.path.join(os.path.abspath(os.path.dirname(__file__)), "sx_real.py")
os.execlp("ipython","ipython","-i",e2real)
