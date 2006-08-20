#!/usr/bin/env python

import sys
import os
from optparse import OptionParser
from EMAN2 import *
from emimagemx import *
from emimage2d import *
import pyshed
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
	
def imgupd():
	for i in EMImage2D.allim.keys():
		if i.data.get_attr("changecount")!=i.changec :
			i.setData(i.data)

	for i in EMImageMX.allim.keys():
		try:
			if len(i.data)!=i.nimg : i.setData(i.data)
		except:
			pass
		for j in i.changec.keys():
			upd=0
			if j.get_attr("changecount")!=i.changec[j] :
				upd=1
				break
		if upd : i.setData(i.data)


if __name__ == "__main__":
	app = QtGui.QApplication([])
	window = pyshed.Shell()
	window.show()
	
	tmr=QtCore.QTimer()
	tmr.setInterval(300)
	tmr.connect(tmr,QtCore.SIGNAL("timeout()"), imgupd)
	tmr.start()
	
	sys.exit(app.exec_())
