#!/usr/bin/env python

import os

e2real=os.getenv("EMAN2DIR")+"/bin/sx_real.py"
os.execlp("ipython","ipython","-i",e2real)