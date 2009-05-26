#!/usr/bin/env python
#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
#
#

# emimage.py  12/07/2005   Steven Ludtke

# This is a WxPython widget for display of a single EMData object

import sys
import wx
from time import time
from weakref import WeakKeyDictionary
from imagecontrol import *
from libpyGLUtils2 import GLUtil

ts0,ts1,ts2=0,0,0

class EMImage(wx.Dialog):
	"""This class is an image display window for a single EMData object
	created using the WxPython toolkit"""
	allim=WeakKeyDictionary()
	def __init__(self,emdata=None):
		"""Creates a new EMImage instance and optionally initializes it"""
		if emdata: wsize=(emdata.get_xsize(),emdata.get_ysize())
		else: wsize=(200,200)
		sset=wx.SystemSettings
		dim=(sset.GetMetric(wx.SYS_SCREEN_X),sset.GetMetric(wx.SYS_SCREEN_Y))
		wsize=(min(wsize[0],dim[0]-50),min(wsize[1],dim[1]-50))
		
		wx.Dialog.__init__(self,None,size=wsize,style=wx.DEFAULT_DIALOG_STYLE+wx.RESIZE_BORDER)
		
		self.origin=[0,0]		# origin in the image being displayed
		self.scale=1.0			# display scale
		self.mingrayval=0.0		# value in EMData corresponding to min gray level
		self.maxgrayval=1.0		# value in EMData corresponding to max gray level
		self.data=None			# pointer to the displayed EMData object
		self.Show()
		self.bm=None
		self.control=None
		self.dofft=0
		
		if emdata: self.setdata(emdata)
		wx.EVT_PAINT(self, self.OnPaint)
		wx.EVT_SIZE(self, self.OnResize)
		wx.EVT_MOUSE_EVENTS(self, self.OnMouse)
		EMImage.allim[self]=0
		del self		# decrement the reference count so allim doesn't prevent deletion
		
#	def __del__(self):
#		EMImage.allim.remove(self)
	
	def changed(self):
		"""This should be called whenever the image needs to be redisplayed. If the actual image
		data changes, setdata should be called again"""
		
		if not self.data: return
#		print "  %s (%f %f) (%f %f)                \r"%(str(self.origin),self.mingrayval,self.maxgrayval,self.data.get_attr("minimum"),self.data.get_attr("maximum")),
#		sys.stdout.flush()
		if self.origin[0]<0 : self.origin[0]=0
		if self.origin[1]<0 : self.origin[1]=0
		d=GLUtil.render_amp8(self.data,self.origin[0],self.origin[1],self.GetSize()[0],self.GetSize()[1],self.GetSize()[0]*3,self.scale,1,254, self.mingrayval,self.maxgrayval,1)
		
		image=wx.EmptyImage(self.GetSize()[0],self.GetSize()[1])
		image.SetDataBuffer(d)
		self.bm=image.ConvertToBitmap()
		
		self.Refresh(0)
		self.Update()
		self.changec=self.data.get_attr("changecount")
		
	def setdata(self,emdata):
		"""called whenever the image to be displayed changes. A simple call to changed will not
		modify the display limits"""
		self.rdata=emdata
#		if self.control: self.control.set_target(self)
		if self.dofft:
			self.data=emdata.do_fft()
			self.data.set_value_at(0,0,0,0)
			self.data.set_value_at(1,0,0,0)
			self.data.process_inplace("xform.fourierorigin",{})
			self.data=self.data.get_fft_amplitude()
		else:
			self.data=self.rdata
		self.mingrayval=self.data.get_attr("minimum")
		self.maxgrayval=self.data.get_attr("maximum")
		self.changed()
		
	def setfft(self,dof):
		self.dofft=dof
		self.setdata(self.rdata)
	
	def setminmax(self,m0,m1):
		self.mingrayval=m0
		self.maxgrayval=m1
		self.changed()
	
	def OnResize(self,event):
		self.changed()
	
	def OnMouse(self,event):
#		print event.m_x,event.m_y,event.LeftIsDown(),event.RightIsDown(),event.RightDown()
		if event.RightDClick() :
			self.origin=[0,0]
			self.down=None
			self.changed()
		elif event.RightDown() :
			self.down=(event.m_x,event.m_y)
			self.begin=self.origin
		elif event.MiddleDown() :
			if not self.control :self.control=ImageControl(None)
			self.control.set_target(self)
			self.control.Show()
		elif event.Dragging() and event.RightIsDown():
			self.origin=[self.begin[0]-(event.m_x-self.down[0]),self.begin[1]+(event.m_y-self.down[1])]
			self.changed()
		elif event.RightUp():
			self.down=None
	
	def OnPaint(self,event):
		dc=wx.BufferedPaintDC(self)
		if not self.bm :
			dc.Clear()
			return
		dc.DrawBitmap(self.bm,0,0,0)
#		print "*",
#		sys.stdout.flush()
#		dc.SetFont(self.font)
#		dc.SetBackground(wx.WHITE_BRUSH)
#		dc.Clear()
#		dc.SetPen(wx.BLACK_PEN)
#		dc.SetBrush(wx.TRANSPARENT_BRUSH)
