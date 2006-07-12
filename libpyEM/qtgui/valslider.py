import sys
from PyQt4 import QtCore, QtGui

def clamp(x0,val,x1):
	return max(min(val,x1),x0)


class ValSlider(QtGui.QWidget):
	"""The valslider class represents a connected text widget and horizontal slider.
	setValue(float) - to programatically change the value
	emit valueChanged(float)
	"""
	def __init__(self, parent, range=None, label=None):
		if not parent: raise Exception,"ValSliders must have parents"
		QtGui.QWidget.__init__(self,parent)
		
		if range : self.range=list(range)
		else : self.range=[0,1.0]
		self.value=0.0
		self.ignore=0
		
		self.hboxlayout = QtGui.QHBoxLayout(self)
		self.hboxlayout.setMargin(0)
		self.hboxlayout.setSpacing(6)
		self.hboxlayout.setObjectName("hboxlayout")
		
		self.label = QtGui.QLabel(self)
		
		if label:
			sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Policy(5),QtGui.QSizePolicy.Policy(0))
#			sizePolicy.setHorizontalStretch(1)
#			sizePolicy.setVerticalStretch(0)
#			sizePolicy.setHeightForWidth(self.text.sizePolicy().hasHeightForWidth())
			self.label.setSizePolicy(sizePolicy)
			self.label.setMinimumSize(QtCore.QSize(30,0))
			self.label.setObjectName("label")
			self.hboxlayout.addWidget(self.label)
		
			self.label = QtGui.QLabel(self)
			self.label.setText(label)
		
		self.text = QtGui.QLineEdit(self)
		
		sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Policy(7),QtGui.QSizePolicy.Policy(0))
		sizePolicy.setHorizontalStretch(1)
		sizePolicy.setVerticalStretch(0)
		sizePolicy.setHeightForWidth(self.text.sizePolicy().hasHeightForWidth())
		self.text.setSizePolicy(sizePolicy)
		self.text.setMinimumSize(QtCore.QSize(80,0))
		self.text.setObjectName("text")
		self.hboxlayout.addWidget(self.text)
		
		self.slider = QtGui.QSlider(self)
		
		sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Policy(7),QtGui.QSizePolicy.Policy(0))
		sizePolicy.setHorizontalStretch(7)
		sizePolicy.setVerticalStretch(0)
		sizePolicy.setHeightForWidth(self.slider.sizePolicy().hasHeightForWidth())
		self.slider.setSizePolicy(sizePolicy)
		self.slider.setMinimumSize(QtCore.QSize(100,0))
		self.slider.setMaximum(4095)
		self.slider.setSingleStep(16)
		self.slider.setPageStep(256)
		self.slider.setOrientation(QtCore.Qt.Horizontal)
		self.slider.setObjectName("slider")
		self.hboxlayout.addWidget(self.slider)
		
		QtCore.QObject.connect(self.text, QtCore.SIGNAL("editingFinished()"), self.textChange)
		QtCore.QObject.connect(self.slider, QtCore.SIGNAL("valueChanged(int)"), self.sliderChange)

	def setRange(self,minv,maxv):
		self.range=[minv,maxv]
		self.updates()

	def setValue(self,val):
		self.value=val
		self.updateboth()
		self.emit(QtCore.SIGNAL("valueChanged"),self.value)
		
	def textChange(self):
		if self.ignore : return
		x=self.text.text()
		if len(x)==0 : return
		if x[0]=='<' : 
			try: self.range[1]=float(x[1:])
			except: pass
			self.updateboth()
		elif x[0]=='>' : 
			try: self.range[0]=float(x[1:])
			except: pass
			self.updateboth()
		else:
			try:
				self.value=float(x)
				self.updates()
				self.emit(QtCore.SIGNAL("valueChanged"),self.value)
			except:
				self.updateboth()
				
	def sliderChange(self,x):
		if self.ignore : return
		self.value=(self.slider.value()/4095.0)*(self.range[1]-self.range[0])+self.range[0]
		self.updatet()
		self.emit(QtCore.SIGNAL("valueChanged"),self.value)
		
		
	def updates(self):
		self.ignore=1
		self.slider.setValue(clamp(0,(self.value-self.range[0])/(self.range[1]-self.range[0])*4095.0,4095.0))
		self.ignore=0

	def updatet(self):
		self.text.setText(str(self.value)[:self.text.width()/10-1])
		
	def updateboth(self):
		self.updates()
		self.updatet()
		