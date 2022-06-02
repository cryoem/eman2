#! /usr/bin/env python
from __future__ import print_function
from __future__ import division
#
# Do not run this script directly. You could run this as
# ipython --gui=qt -i e2_real.py
# otherwise it is normally run via the e2.py script
#
# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Tapu Shaikh 2019 (tapu.shaikh@mpi-dortmund.mpg.de)
# Author: Adnan Ali 2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Luca Lusnig 2019 (luca.lusnig@mpi-dortmund.mpg.de)
# Author: Toshio Moriya 2019 (toshio.moriya@kek.jp)
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
#
# Copyright (c) 2019 Max Planck Institute of Molecular Physiology
# Copyright (c) 2000-2006 Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

import EMAN2
from EMAN2 import *
import PyQt5
import eman2_gui.emapplication
import eman2_gui.emimage
import os
from sphire.libpy import sp_global_def
from sphire.libpy.sp_sparx import *

def main():
    GUIUSE = True

    try:
        if EMAN2.get_platform() == "Linux" and os.getenv("DISPLAY") == None:
            raise Exception

        # import IPython.lib.inputhook

        app = eman2_gui.emapplication.EMApp()
        # IPython.lib.inputhook.enable_qt4(app)

        def ipy_on_timer():
            eman2_gui.emimage.image_update()

        ipytimer = PyQt5.QtCore.QTimer()
        ipytimer.timeout.connect(ipy_on_timer)
        ipytimer.start(200)

        EMAN2.GUIMode = True
        EMAN2.app = app
    except:
        GUIUSE = False


    if GUIUSE:
        print("Welcome to the interactive SPARX-GUI Python interface, provided by ipython")
    else:
        print(
            "Welcome to the interactive SPARX-NoGUI Python interface, provided by ipython"
        )

    print("  ", sp_global_def.SPARXVERSION)

if __name__ == "__main__":
	main()
