#!/usr/bin/env python
# Muyuan Chen 2016-09
from EMAN2 import *
import numpy as np
from emapplication import EMApp
from emimage2d import EMImage2DWidget
from emshape import EMShape

import PyQt4
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from PyQt4.QtCore import QTimer

def main():
	
	usage="[prog] <2D image file> "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)


	
	filename=args[0]
	
	app = EMApp()
	img=EMData(filename)
	
	
	#print img[0]["mean"]
	w=EMBreakBrick(img,app)
	#w.set_data(img,filename)
	
	app.show_specific(w)
	app.exec_()
	
	
	
	
	
	E2end(logid)
	






class EMBreakBrick(EMImage2DWidget):
	
	def __init__(self,img, app):
		EMImage2DWidget.__init__(self,img,application=app)
		
		self.sx=img["nx"]
		self.sy=img["ny"]
		self.set_scale(.3)
		
		minxy=min(self.sx, self.sy)
		self.bar_len=minxy*.1
		self.bar_ypos=-self.sy*.2
		self.bar_thick=20
		self.bar_xpos=self.sx/2
		self.bar=EMShape()
		self.barspeed=0.02*minxy
		
		self.score=0
		self.score_label=EMShape()
		
		
		self.ball=EMShape()
		self.ball_rad=20
		self.ball_pos=np.array([self.bar_xpos, self.bar_ypos+self.bar_thick+self.ball_rad])
		self.ball_speed=minxy*.01		
		self.set_shapes({0:self.bar, 1:self.ball,2:self.score_label})
		
		self.game_started=False
		
		self.data.mult(-1)
		#self.data.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.05})
		#self.data.process_inplace("normalize")
		#self.data.process_inplace("threshold.belowtozero",{"minval":1})
		#print self.data["mean_nonzero"],self.sx,self.sy
		self.data.div(self.data["mean_nonzero"]*self.sx*self.sy)
		self.data.mult(1000)
		self.auto_contrast()
		self.data_grad=self.data.process("math.gradient.direction")
		self.del_msk=self.data.copy()
		self.del_msk.to_one()
		
		self.update_bar()
		self.update_score()
		self.update_ball()
		print "break bricks.."
	
	def update_score(self):
		self.score_label.setShape(["label", 1,1,1, -self.sx*.1, self.sy*1.1, "{:.02f}".format(self.score), 100,5])
		self.update_shapes({2:self.score_label})
		self.updateGL()
	
	def update_bar(self):
		self.bar.setShape(["line",1,1,1,self.bar_xpos-self.bar_len,self.bar_ypos,self.bar_xpos+self.bar_len,self.bar_ypos,self.bar_thick])
		
		self.update_shapes({0:self.bar})
		self.updateGL()
		
		if not self.game_started:
			self.ball_pos=np.array([self.bar_xpos, self.bar_ypos+self.bar_thick+self.ball_rad])
			self.update_ball()
	
	
	def update_ball(self):
		self.ball.setShape(["circle",1,1,1,self.ball_pos[0], self.ball_pos[1], self.ball_rad,2])
		#self.set_shapes({1:self.ball})
		self.update_shapes({1:self.ball})
		self.updateGL()
	
	def start_game(self):
		print "Start~"
		self.ball_ori=45./180.*np.pi
		self.ball_vec=np.array([self.ball_speed*np.cos(self.ball_ori),self.ball_speed*np.sin(self.ball_ori)])
		
		self.game_timer=QTimer()
		self.game_timer.timeout.connect(self.time_update)
		self.game_timer.start(30)
		self.game_started=True
		
	def time_update(self):
		#print "timer~~"
		self.ball_move()
		self.update_ball()
	
	def ball_move(self):
		newpos=self.ball_pos+self.ball_vec
		if newpos[0]>self.sx or newpos[0]<0: self.ball_vec[0]*=-1
		elif newpos[1]>self.sy: self.ball_vec[1]*=-1
		elif abs(newpos[1]-self.bar_ypos)<self.bar_thick and abs(newpos[0]-self.bar_xpos)<self.bar_len:
				self.ball_vec[1]*=-1
		
		elif newpos[1]<self.bar_ypos-200:
			print "Lose..."
			self.game_timer.stop()
		
		else:
			p=newpos.astype(int)
			if p[1]<0 or p[0]<0 or p[0]>self.sx or p[1]>self.sy:
				val=0
			else:
				val=self.data.get_value_at(p[0],p[1])
			if val>0:
				grad=self.data_grad.get_value_at(p[0],p[1])
				#print val, grad
				newori=(grad*2-self.ball_ori)% (np.pi*2)
				oridiff=abs(newori-self.ball_ori)*180./np.pi
				if oridiff>180: oridiff=oridiff-180
				#print self.ball_ori*180./np.pi,newori*180./np.pi,oridiff
				if oridiff<30:
					#print "???"
					newori+=.5*np.pi*np.sign(newori-self.ball_ori)
				self.ball_ori=newori
				self.ball_vec=np.array([self.ball_speed*np.cos(self.ball_ori),self.ball_speed*np.sin(self.ball_ori)])
				
				self.del_msk.to_one()
				self.del_msk.process_inplace("mask.soft", {"dx":p[0]-self.sx/2, "dy": p[1]-self.sy/2, "outer_radius":30})
				delimg=self.data*self.del_msk
				delval=delimg["mean"]*self.sx*self.sy
				self.score+=delval
				self.update_score()
				
				self.data.sub(delimg)
				#self.data_grad=self.data.process("math.gradient.direction")
				self.force_display_update()
				
			
		
		
		self.ball_pos+=self.ball_vec
		
		
	
	def keyPressEvent(self,event):

		if event.key() == Qt.Key_Right:
			if self.bar_xpos+self.bar_len<self.sx:
				self.bar_xpos+=self.barspeed
				self.update_bar()
		elif event.key() == Qt.Key_Left:
			if self.bar_xpos-self.bar_len>0:
				self.bar_xpos-=self.barspeed
				self.update_bar()
		elif event.key() ==Qt.Key_Space:
			if not self.game_started:
				self.start_game()
			
			
	
	
if __name__ == '__main__':
	main()
	