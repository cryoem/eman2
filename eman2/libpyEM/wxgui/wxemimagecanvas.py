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

import wx
from math import *

DEBUG = 0
DEBUG2 = 0
BUFFERED = 1

class EMImageCanvas(wx.Window):
	def __init__(self, parent,id=-1,size=(100,100)):
		wx.Window.__init__(self,parent,id, (0,0), size=size, style=wx.SIMPLE_BORDER)

		self.target = 0

		# Give temporary default values
		self.scale = 1
		self.center = [-self.GetSize()[0]/2,self.GetSize()[1]/2]

		self.linecoords = ((0,0),(0,0),(0,0),(0,0),(0,0))

		self.lowerleft = self.linecoords[0]
		self.begin=[0,0]	
		self.down=[0,0]
		self.x_per = 0.5
		self.y_per = 0.5

		self.center_line1 = ((0,0),(0,0))
		self.center_line2 = ((0,0),(0,0))
		
		# Stuff for drawing.
		self.maxWidth=self.GetSize()[0]
		self.maxHeight=self.GetSize()[1]
		self.x = self.y = 0
		self.drawing = False
		self.SetBackgroundColour("WHITE")

		# Initial draw.
		if BUFFERED:
			self.buffer = wx.EmptyBitmap(self.maxWidth, self.maxHeight)
			dc = wx.BufferedDC(None,self.buffer)
			dc.SetBackground(wx.Brush(self.GetBackgroundColour()))
			dc.Clear()
			self.DoDrawing(dc)

		self.Bind(wx.EVT_LEFT_DOWN, self.OnLeftButtonEvent)
		self.Bind(wx.EVT_MOTION,    self.OnLeftButtonEvent)
		self.Bind(wx.EVT_PAINT,		self.OnPaint)

	def set_target(self,t):
		self.target=t
		self.scale = self.target.scale
		self.origin = self.target.origin
		self.ReDraw()

	def getWidth(self):
		return self.maxWidth

	def getHeight(self):
		return self.maxHeight

	def OnPaint(self,event):
		"""Paint event."""
		if BUFFERED:
			dc = wx.BufferedPaintDC(self,self.buffer, wx.BUFFER_VIRTUAL_AREA)
		else:
			dc = wx.PaintDC(self)
			self.PrepareDC(dc)
			self.DoDrawing(dc)


	def MakeBox(self):

		if self.target:

			width = (self.target.GetSize()[0] * 1/self.target.scale) / self.target.emdatasize[0] * self.getWidth() / 2
			height = (self.target.GetSize()[1] * 1/self.target.scale) / self.target.emdatasize[1] * self.getHeight() / 2

			#print "frame size: " + str(self.target.GetSize())

			self.center[0] = int(-(self.x_per * float(self.GetSize()[0])))
			self.center[1] = int(-((self.y_per - 1.0) * float(self.GetSize()[0])))

			self.linecoords = (((-self.center[0] - width),(self.center[1] - height)),
			 	  ((-self.center[0] + width ),(self.center[1] - height)),
				  ((-self.center[0] + width ),(self.center[1] + height )),
				  ((-self.center[0] - width),(self.center[1] + height )),
				  ((-self.center[0] - width),(self.center[1] - height)))

			self.center_line2 = ((0,0),(0,0))

			self.center_line1 = ((-self.center[0]-5,self.center[1]-5),(-self.center[0]+5,self.center[1]+5))
			self.center_line2 = ((-self.center[0]-5,self.center[1]+4),(-self.center[0]+5,self.center[1]-6))


	def DoDrawing(self,dc,printing=False):
		"""Draw all the boxes and lines for each screen refresh."""
		self.MakeBox()

		dc.SetPen(wx.Pen("BLACK"))

		dc.DrawLines(self.linecoords, 0)

		dc.SetPen(wx.Pen("RED"))

		dc.DrawLines(self.center_line1, 0)
		dc.DrawLines(self.center_line2, 0)



		#dc.DrawLines(self.linecoords2, 0)


	def SetXY(self,event):
		"""Set x and y coords."""
		self.x,self.y = self.ConvertEventCoords(event)

	def ConvertEventCoords(self,event):
		"""Convert to coords."""
		xView, yView = self.GetViewStart()
		xDelta, yDelta = self.GetScrollPixelsPerUnit()
		return (event.GetX() + (xView * xDelta),
						event.GetY() + (yView * yDelta))

	def ReDraw(self):
		"""Force screen redraw, e.g., when option is changed."""
		self.drawing = True
		if BUFFERED:
			cdc = wx.ClientDC(self)
			self.PrepareDC(cdc)
			dc = wx.BufferedDC(cdc, self.buffer)
			dc.Clear()
		else:
			dc = wx.ClientDC(self)
			self.PrepareDC(dc)

		dc.BeginDrawing()
		self.DoDrawing(dc)
		dc.EndDrawing()
		self.drawing = False
	


	def OnLeftButtonEvent(self,event):
		if event.RightDClick():
			self.center=[0,0]
			self.down=None
			self.UpdateCanvas()
			self.ReDraw()
		elif event.RightDown():
			self.down=(event.m_x,event.m_y)
			self.begin=self.center
		elif event.Dragging() and event.RightIsDown():
			self.center=[self.begin[0]-(event.m_x-self.down[0]),self.begin[1]+(event.m_y-self.down[1])]
			self.UpdateCanvas()
		elif event.RightUp():
			self.down=None

	def UpdateCanvas(self):

		if self.target:

			self.x_per = -self.center[0] / float(self.GetSize()[0])
			self.y_per = 1.0 - self.center[1] / float(self.GetSize()[1])

			self.target.x_per = self.x_per
			self.target.y_per = self.y_per

			self.target.origin = ((self.x_per * self.target.GetSize()[1]),(self.y_per * self.target.GetSize()[1]))

			self.ReDraw()
			self.target.changed()
