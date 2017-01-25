#!/usr/bin/env python

import os

e2real=os.getenv("EMAN2DIR")+"/bin/e2_real.py"
os.execlp("ipython","ipython","-i","--gui=qt",e2real)