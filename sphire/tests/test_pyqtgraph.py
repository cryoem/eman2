# -*- coding: utf-8 -*-
"""
Demonstrates GLVolumeItem for displaying volumetric data.
"""
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.opengl as gl
import sp_utilities
import numpy as np

app = QtGui.QApplication([])
w = gl.GLViewWidget()
w.orbit(256, 256)
w.setCameraPosition(0, 0, 0)
w.opts['distance'] = 400
w.show()
w.setWindowTitle('Plotting a Volume')


g = gl.GLGridItem()
g.scale(60, 60, 1)
w.addItem(g)


data = sp_utilities.get_im('/home/adnan/DemoResults/12_POSTREFINER/vol_combined.hdf')
da = data.get_3dview()
d2 = np.empty(da.shape + (4,), dtype=np.ubyte)
d2[..., 0] = da * (255./(da.max()/1))
d2[..., 1] = d2[..., 0]
d2[..., 2] = d2[..., 0]
d2[..., 3] = d2[..., 0]
d2[..., 3] = (d2[..., 3].astype(float) / 255.)**2 * 255
d2[d2>125] = 0

# RGB orientation lines (optional)
d2[:, 0, 0] = [255, 0, 0, 255]
d2[0, :, 0] = [0, 255, 0, 255]
d2[0, 0, :] = [0, 0, 255, 255]



v = gl.GLVolumeItem(d2, sliceDensity=1, smooth=True, glOptions='additive')
v.translate(-d2.shape[0]/2, -d2.shape[1]/2, -150)
w.addItem(v)

## Start Qt event loop unless running in interactive mode.
if __name__ == '__main__':
    import sys

    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
