#!/usr/bin/env python
#
# Author: Steven Ludtke (sludtke@bcm.edu)
# Copyright (c) 2011- Baylor College of Medicine
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
from EMAN2 import *
from EMAN2db import e2gethome
from EMAN2jsondb import js_open_dict
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt, QTimer, QBuffer
from .emapplication import EMApp
from .emimage2d import *
from .emimagemx import *
from .empdbitem3d import *
from .emplot2d import *
from .emhist import *
from .emplot3d import *
from .emscene3d import *
from .emdataitem3d import *
from libpyUtils2 import EMUtil
from .matching import matches_pats
from .valslider import StringBox
import os
import threading
import weakref
import random
import struct
import socket
from PIL import Image
import io


class EMDisplayServerWidget(QtWidgets.QWidget):
	"""This widget listens for connections on port 31980, and opens appropriate display widgets
	of appropriate types for visualization of Jupyter Lab sessions"""
	widget_types={
		"image":EMImage2DWidget,
		"imagemx":EMImageMXWidget,
		"volume":EMScene3D,
		"plot2d":EMPlot2DWidget,
		"plot3d":EMPlot3DWidget,
		"histogram":EMHistogramWidget
		}
	
	def __init__(self,port=31980):
		self.port=port
		self.app=EMApp.instance()

		# although this looks dumb it is necessary to break Python's issue with circular imports(a major weakness of Python IMO)
		#global emscene3d, emdataitem3d
		#from . import emscene3d
		#from . import emdataitem3d

		QtWidgets.QWidget.__init__(self,None)
		self.gbl = QtWidgets.QGridLayout(self)
		self.gbl.setContentsMargins(8, 8, 8, 8)
		self.gbl.setSpacing(6)

		# widget list
		self.wdglist=QtWidgets.QListWidget(self)
		self.wdglist.setSizePolicy(QtWidgets.QSizePolicy.Preferred,QtWidgets.QSizePolicy.Expanding)
		#for i in imagenames:
			#self.setlist.addItem(i)
		self.gbl.addWidget(self.wdglist,0,0)
		
		self.wdglist.currentRowChanged[int].connect(self.wdgSel)
	
		self.widgets={wt:{} for wt in self.widget_types}		# dict of dicts for widgets organized by type
		#self.threads=[]
		#self.qin=queue.Queue(0)		# queue for incoming messages. All GUI activity handled in main thread
		self.quit=False				# set to trigger thread exit

		# Timer event triggers GUI updates in the main thread
		self.msg=None				# the receiving thread sets this when a packet comes in
		self.reply=None				# the main thread sets this with a response for the requesting program
		self.timer = QtCore.QTimer()
		self.timer.timeout.connect(self.time_out)
		self.timer.start(0.5)		# check for new data twice a second
		
		# Launch the server thread
		self.serverthread=threading.Thread(target=self.listenthread)
		self.serverthread.start()
		
	def closeEvent(self,event) :
		print("Server exiting. Waiting for network shutdown.")
		self.quit=True
		self.serverthread.join()
	
	def rawtopng(self,depth,nx,ny,data):
		if depth==1: 
			with io.BytesIO() as bio:
				Image.frombytes("L",(nx,ny),data).save(bio,"png")
				png=bio.getvalue()
				bio.close()
		else: 
			with io.BytesIO() as bio:
				Image.frombytes("RGB",(nx,ny),data).save(bio,"png")
				png=bio.getvalue()
				bio.close()
		return png
	
	def time_out(self):
		""" This is the routine which actually displays incoming data.
		This runs in the main thread, so it can talk to Qt safely."""
		if self.quit: self.app.quit()
		if self.msg==None: return
		try: 
			data,vtype,vname,dname,settings=self.msg
		except:
			print("Error :",self.msg)
			self.reply=["error","Bad Parameters"]
			self.msg=None
			return
		
		self.msg=None
		
		if vtype not in self.widget_types:
			self.reply=["error","bad visualization type"]
			return
		
		if vname==None: vname="default"
		if vname not in self.widgets[vtype]:
			self.widgets[vtype][vname]=self.widget_types[vtype]()		# create a new display widget
			self.wdglist.addItem(f"{vtype}_{vname}")

		widget=self.widgets[vtype][vname]
		widget.show()
		widget.raise_()
	
		if vtype=="image":
			widget.set_data(data,dname)
			#depth,nx,ny,raw=widget.render_bitmap()
			#png=self.rawtopng(depth,nx,ny,raw)
		elif vtype=="imagemx":
			widget.set_data(data, dname)
		elif vtype=="volume":
			pass
		elif vtype=="plot2d":
			widget.set_data(data,key=dname)
#			widget.set_data(data,key=dname,replace=False,quiet=False,color=-1,linewidth=1,linetype=-2,symtype=-2,symsize=10,comments=None)
		elif vtype=="plot3d":
			widget.set_data(data,key=dname)
		elif vtype=="histogram":
			pass
		
		pix=widget.renderPixmap()
		nx=pix.size().width()
		ny=pix.size().height()
		buffer=QBuffer();
		buffer.open(QBuffer.ReadWrite)
		pix.save(buffer,"PNG")
		png=buffer.data().data()	# this is correct. The first data gets the QByteArray, the second gets it's data

		self.reply=["png",(nx,ny,png)]
		
	# EMImage2DWidget - set_data(self,incoming_data,file_name="",retain_current_settings=True, keepcontrast=False, xyz=2)
	# EMImageMXWidget - def set_data(self, obj, filename = '', update_gl = True, soft_delete = False) :
	# EMScene3D - insertNewNode(self, name, node, parentnode=None, parentidx=None)
	# EMDataItem3D - setData(self, data, n=0)
	# EMPlot2DWidget - set_data(self,input_data,key="data",replace=False,quiet=False,color=-1,linewidth=1,linetype=-2,symtype=-2,symsize=10,comments=None)
	# EMPlot3DWidget - set_data(self,input_data,key="data",replace=False,quiet=False,color=-1,linewidth=1,linetype=-2,symtype=-2,symsize=10,comments=None)
	# EMHistogramWidget - set_data(self,input_data,key="data",replace=False,quiet=False,color=-1,alpha=0.8,rwidth=0.8)
	
	def wdgSel(self,row):
		"""When the user highlights a new widget"""
		curlabel=str(self.wdglist.currentItem().text())
		wtype,wname=curlabel.split("_",1)
		wdg=self.widgets[wtype][wname]
		wdg.show()
		wdg.raise_()
		
	def listenthread(self):
		"""This thread listens for new connections and completes a single transaction for each connection"""
		lsock=socket.socket() 
		lsock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
		lsock.settimeout(2.0)
		try: lsock.bind(("",self.port))
		except:
			print("Error: cannot bind port ",self.port)
			self.quit=True
			return
		lsock.listen(3)          # Wait for 'connect' requests
		print("listening on port ",self.port)

		# _very_ basic security
		magicpath=f"{e2gethome()}/.eman2/server_magic"
		if os.path.isfile(magicpath):
			with open(magicpath,"rb") as fin:
				self.magic=fin.read(8)
		else:
			self.magic=struct.pack("Q",random.getrandbits(64))
			out=open(magicpath,"wb")
			out.write(self.magic)
			fid=out.fileno()
			try: os.fchmod(fid,stat.S_IRUSR|stat.S_IWUSR)
			except: pass
			out.close()

		while True:
			try:
				sock2=lsock.accept()		# accept the connection (new socket)
			except:
				if self.quit: 
					lsock.close()
					break		# gives us the opportunity for the program to exit every 2 seconds
#				lsock.shutdown(socket.SHUT_RDWR)
#				print("timeout")
				continue
			msg=sock2[0].recv(8)	# receive 8 bytes as a security check
			if msg!=self.magic:
				time.sleep(3)
				sock2[0].close()			# no reply on failure
				print(f'Security mismatch {struct.unpack("Q",msg)} {struct.unpack("Q",self.magic)}')
				continue
				
			# if we got here we have a good connection
			#self.threads.append(threading.Thread(target=self.recvthread,args=[sock2]))
			
			try:
				self.msg=sock_recvobj(sock2[0])		# pass the data to the main thread
			except:
				import traceback
				traceback.print_exc()
				print("Error with:",self.msg)
				sock2[0].close()
				continue
			
			while True:
				time.sleep(0.1)
				if self.reply!=None: 
					sock_sendobj(sock2[0],self.reply)
					self.reply=None
					break
			
			
	#def recvthread(self,sock):
		#"""This thread receives new data to display, and responds with an image or an error"""
		#data=recvobj(sock)
		
		
