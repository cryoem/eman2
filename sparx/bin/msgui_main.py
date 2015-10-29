#!/usr/bin/env python
#
# Author: Toshio Moriya, 10/23/2006
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#

import sys
import os
# import global_def
from global_def import *
from PyQt4.Qt import *
from PyQt4 import QtGui
from PyQt4 import QtCore
# import subprocess
from subprocess import *
from EMAN2 import *
from sparx import *
from EMAN2_cppwrap import *

# ========================================================================================
# Layout of the Pop Up window SXPopup_info; started by the function info of the main window
class SXPopup_info(QWidget):
	def __init__(self):
		QWidget.__init__(self)
		#Here we just set the window title and  3 different labels, with their positions in the window
		self.setWindowTitle('MPI-Sparx GUI Info Page')
		title1=QtGui.QLabel('<b>MPI-Sparx GUI Prototype</b>', self)
		title1.move(10,10)
		title2=QtGui.QLabel('<b>Authors: Toshio Moriya</b> ', self)
		title2.move(10,40)
		title3=QtGui.QLabel('For more information visit:\nhttp://sparx-em.org/sparxwiki ', self)
		title3.move(10,70)

# ========================================================================================
class SXcmd_token:
	def __init__(self):
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.name = ''			# Use name of argument or option in command line
		self.label = ''			# Use an user friendly name of argument or option
		self.help = ''		    # Help info
		self.type = ''		    # sxtype_arg or sxtype_option
		self.group = ''			# sxgroup_main or sxgroup_advance
		self.widget = ''		# sxwidget_line_edit or sxwidget_check_box
		self.default = None		# default value
		self.gui = None
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

# ========================================================================================
# sx_start
# Provides all necessary functionarity
# tabs only contains gui and knows how to layout them
class SXPopup(QWidget):
	def __init__(self, sxscript, short_info):
		QWidget.__init__(self)

		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.sxscript = sxscript
		self.short_info = short_info
		self.mpi_support = False
		self.mpi_add_flag = False
		self.cmd_token_list = []	
		self.cmd_token_dict = {} # For ease of access, use dictionary version of cmd_token_list using label for the key 
		
		self.sxcmd_line = ""
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

		# Set parser from sx*.py script
		self.extract_info_from_parser()
		
		self.setWindowTitle(self.sxscript)
		self.tab_main = SXTab_main(self)
		self.tab_advance = SXTab_advance(self)
		self.tab_main.w1 = self.tab_advance
		self.TabWidget = QtGui.QTabWidget()
		self.TabWidget.insertTab(0,self.tab_main,'Main')
		self.TabWidget.insertTab(1,self.tab_advance,'Advanced')
		self.TabWidget.resize(730,860)
		self.TabWidget.show()

	# NOTE: 20015/10/23 Toshio Moriya
	# For now, please overwrite this method in the child class
	# In the future, it should import the sx"*.py cript at run time
	def get_option_parser(self):
		assert False, 'Unreachable Code!!!'
		return None
		
	def set_group_cmd_token_list(self, parser_option_list, sxcmd_token_type):		
		# Loop through parser_option_list
		for parser_option in parser_option_list:
			# Do not add for the following case
			if parser_option.help == None or parser_option.help == "" or parser_option.help[0] != "<":
				continue
			
			# Create new SXcmd_token for this parser_option
			cmd_token = SXcmd_token()
			cmd_token.name = parser_option.dest

			# Extract label and help from parser_option help
			assert(parser_option.help[0] == "<")
			cmd_token.label = parser_option.help[1:].split(">")[0]
			cmd_token.help = parser_option.help
			
			# For command type and command prefix, use passed value 
			cmd_token.type = sxcmd_token_type
			
			# Check if this is advanced option, and set group
			if parser_option.help[-1 * len("(advanced)"):] != "(advanced)":
				cmd_token.group = 'sxgroup_main'
			else:
				cmd_token.group = 'sxgroup_advance'
			
			# Set gui type and default value for parser_option
			if parser_option.action != "store_true":
				cmd_token.widget = 'sxwidget_line_edit'
				cmd_token.default = str(parser_option.default)
			else:
				cmd_token.widget = 'sxwidget_check_box'
				cmd_token.default = parser_option.default
			
			self.cmd_token_list.append(cmd_token)
			self.cmd_token_dict[cmd_token.label] = cmd_token
				
	# Set self.mpi_support and self.sxcmd_token_list by extracting info from option parser
	def extract_info_from_parser(self):
		
		parser = self.get_option_parser()

		# Check if the script supports MPI
		mpi_group = parser.get_option_group("--sxgui_mpi_options")
		self.mpi_support = mpi_group.get_option("--MPI_support").default
		if self.mpi_support:
			self.mpi_add_flag =  mpi_group.get_option("--MPI_add_flag").default
		else:
			self.mpi_add_flag = False
			
		# Create GUI components for command arguments and options
		self.sxcmd_token_list = []
		self.set_group_cmd_token_list(parser.get_option_group("--sxgui_arguments").option_list, 'sxtype_arg')
		self.set_group_cmd_token_list(parser.option_list, 'sxtype_option')
	
	def map_gui_to_cmd_line(self):
		# Add program name to command line
		sxcmd_line = "%s.py" % self.sxscript
		
		# Loop through all command tokens
		for cmd_token in self.cmd_token_list:
			if cmd_token.type == 'sxtype_arg':
				assert cmd_token.widget == 'sxwidget_line_edit'
				assert type(cmd_token.gui) is QtGui.QLineEdit
				sxcmd_line += " %s" % (cmd_token.gui.text())
			elif cmd_token.type == 'sxtype_option':
				if cmd_token.widget == 'sxwidget_line_edit':
					sxcmd_line += " --%s=%s" % (cmd_token.name, cmd_token.gui.text())
				elif cmd_token.widget == 'sxwidget_check_box':
					if cmd_token.gui.checkState() == Qt.Checked:
						sxcmd_line += " --%s" % (cmd_token.name)
				else:
					assert False, 'Unreachable Code!!!'
			else:
				assert False, 'Unreachable Code!!!'
		
		# Add MPI flag only if mpi is supported and number of MPI processer (np) is larger than 1
		# NOTE: 2015/10/27 Toshio Moriya
		# This is not elegant but can be removed when --MPI flag is removed from all sx*.py scripts 
		if self.mpi_support and int(str(self.tab_main.mpi_nproc_edit.text())) > 1 and self.mpi_add_flag:
			sxcmd_line += ' --MPI'

		self.sxcmd_line = sxcmd_line
		
	def generate_cmd_line(self):
	
		self.map_gui_to_cmd_line()
		
		file_name_out = QtGui.QFileDialog.getSaveFileName(self, "Generate Command Line", options = QtGui.QFileDialog.DontUseNativeDialog)
		if file_name_out != '':
			file_out = open(file_name_out,'w')
			file_out.write(self.sxcmd_line + '\n')
			file_out.close()
			print 'Saved the following command to %s:' % file_name_out
			print self.sxcmd_line

	def execute_cmd(self):
		# Set self.sxcmd_line
		self.map_gui_to_cmd_line()
		
		# out_dir = str(self.cmd_token_dict['Output Directory'].gui.text())
		# if os.path.exists(out_dir):
		# 	print "Output directory " + out_dir + " already exists!"
		# 	return
		
		# If mpi is not supported set number of MPI processer (np) to 1
		np = 1
		if self.mpi_support:
			np = int(str(self.tab_main.mpi_nproc_edit.text()))
				
		# Case 1: queue submission is enabled
		if self.tab_main.qsub_enable_checkbox.checkState() == Qt.Checked:
				
			# If number of MPI processer (np) is 1, add MPI flag to sxcmd_line
			# NOTE: 2015/10/27 Toshio Moriya
			# This is not elegant but can be removed when --MPI flag is removed from all sx*.py scripts 
			if np == 1 and self.mpi_add_flag:
				self.sxcmd_line += ' --MPI'
			
			# Create script for queue submission from a give template
			assert(os.path.exists(self.tab_main.qsub_script_edit.text()))
			file_template = open(self.tab_main.qsub_script_edit.text(),'r')
			
			file_name_qsub_script = 'qsub_' + str(self.tab_main.qsub_job_name_edit.text()) + '.sh'
			file_qsub_script = open(file_name_qsub_script,'w')
			
			for line_io in file_template:
				if line_io.find('XXX_SXMPI_NPROC_XXX') != -1:
					line_io = line_io.replace('XXX_SXMPI_NPROC_XXX', str(np))
				if line_io.find('XXX_SXMPI_JOB_NAME_XXX') != -1:
					line_io = line_io.replace('XXX_SXMPI_JOB_NAME_XXX', str(self.tab_main.qsub_job_name_edit.text()))
				if line_io.find('XXX_SXCMD_LINE_XXX') != -1:
					line_io = line_io.replace('XXX_SXCMD_LINE_XXX', self.sxcmd_line)
				
				file_qsub_script.write(line_io)
				
			file_template.close()
			file_qsub_script.close()
			
			# Generate command line for queue submission
			cmd_line = str(self.tab_main.qsub_cmd_edit.text()) + ' ' + file_name_qsub_script

			print 'Wrote the following sparx command line in the queue submission script: '
			print self.sxcmd_line
			print 'Submitted a job by the following command: '
			print cmd_line
		
		# Case 2: queue submission is disabled, but MPI is enabled
		elif self.mpi_support:
			assert self.tab_main.qsub_enable_checkbox.checkState() != Qt.Checked

			cmd_line = self.sxcmd_line
			# Add MPI execution to command line only if number of MPI processer (np) is larger than 1
			if np > 1:
				cmd_line = str(self.tab_main.mpi_cmd_line_edit.text())
				# If empty string is entered, use a default template
				if cmd_line == '':
					cmd_line = 'mpirun -np XXX_SXMPI_NPROC_XXX XXX_SXCMD_LINE_XXX'
				if cmd_line.find('XXX_SXMPI_NPROC_XXX') != -1:
					cmd_line = cmd_line.replace('XXX_SXMPI_NPROC_XXX', str(np))
				if cmd_line.find('XXX_SXCMD_LINE_XXX') != -1:
					cmd_line = cmd_line.replace('XXX_SXCMD_LINE_XXX', self.sxcmd_line)					
			
			print 'Executed the following command: '
			print cmd_line
		
		process = subprocess.Popen(cmd_line, shell=True)
		self.emit(QtCore.SIGNAL("process_started"), process.pid)
		
	def save_params(self):		
		file_name_out = QtGui.QFileDialog.getSaveFileName(self, "Save Parameters", options = QtGui.QFileDialog.DontUseNativeDialog)
		if file_name_out != '':
			file_out = open(file_name_out,'w')
		
			# Write script name for consistency check upon loading
			file_out.write('%s \n' % (self.sxscript))	

			# Loop through all command tokens
			for cmd_token in self.cmd_token_list:
				val_str = ''
				if cmd_token.widget == 'sxwidget_line_edit':
					assert type(cmd_token.gui) is QtGui.QLineEdit
					val_str = str(cmd_token.gui.text())
				elif cmd_token.widget == 'sxwidget_check_box':
					assert type(cmd_token.gui) is QtGui.QCheckBox
					if cmd_token.gui.checkState() == Qt.Checked:
						val_str = 'YES'
					else:
						val_str = 'NO'
				else:
					assert False, 'Unreachable Code!!!'
				file_out.write('%s == %s \n' % (cmd_token.label, val_str))				

			# Write MPI parameters 
			file_out.write('%s == %s \n' % ('MPI processors', str(self.tab_main.mpi_nproc_edit.text())))		
			file_out.write('%s == %s \n' % ('MPI Command Line Template', str(self.tab_main.mpi_cmd_line_edit.text())))		

			# Write Qsub paramters 
			if self.tab_main.qsub_enable_checkbox.checkState() == Qt.Checked:
				val_str = 'YES'
			else:
				val_str = 'NO'
			file_out.write('%s == %s \n' % ('Submit Job to Queue', val_str))	
			file_out.write('%s == %s \n' % ('Job Name', str(self.tab_main.qsub_job_name_edit.text())))		
			file_out.write('%s == %s \n' % ('Submission Command', str(self.tab_main.qsub_cmd_edit.text())))		
			file_out.write('%s == %s \n' % ('Submission Script Template', str(self.tab_main.qsub_script_edit.text())))		

			file_out.close()

	def load_params(self):		
		file_name_in = QtGui.QFileDialog.getOpenFileName(self, "Load parameters", options = QtGui.QFileDialog.DontUseNativeDialog)
		if file_name_in != '':
			file_in = open(file_name_in,'r')
		
			# Check if this parameter file is for this sx script
			line_in = file_in.readline()
			if line_in.find(self.sxscript) != -1:
				# loop through the other lines
				for line_in in file_in:
					# Extract label (which should be left of '=='). Also strip the ending spaces
					label_in = line_in.split('==')[0].strip() 			
					# Extract value (which should be right of '=='). Also strip all spaces
					val_str_in = line_in.split('==')[1].strip() 
			
					if label_in == "MPI processors":
						self.tab_main.mpi_nproc_edit.setText(val_str_in)
					elif label_in == "MPI Command Line Template":
						self.tab_main.mpi_cmd_line_edit.setText(val_str_in)
					elif label_in == "Submit Job to Queue":
						if val_str_in == 'YES':
							self.tab_main.qsub_enable_checkbox.setChecked(True)
						else:
							assert val_str_in == 'NO'
							self.tab_main.qsub_enable_checkbox.setChecked(False)
					elif label_in == "Job Name":
						self.tab_main.qsub_job_name_edit.setText(val_str_in)
					elif label_in == "Submission Command":
						self.tab_main.qsub_cmd_edit.setText(val_str_in)
					elif label_in == "Submission Script Template":
						self.tab_main.qsub_script_edit.setText(val_str_in)
					else:
						# Get corresponding cmd_token
						cmd_token = self.cmd_token_dict[label_in]
						assert cmd_token.label == label_in
			
						# Set gui values
						if cmd_token.widget == 'sxwidget_line_edit':
							assert type(cmd_token.gui) is QtGui.QLineEdit
							cmd_token.gui.setText(val_str_in)
						elif cmd_token.widget == 'sxwidget_check_box':
							assert type(cmd_token.gui) is QtGui.QCheckBox
							if val_str_in == 'YES':
								cmd_token.gui.setChecked(True)
							else:
								assert val_str_in == 'NO'
								cmd_token.gui.setChecked(False)
						else:
							assert False, 'Unreachable Code!!!'
			else:
				QMessageBox.warning(self, 'Fail to load paramters', 'The specified file is not paramter file for %s.' % self.sxscript)
	
			file_in.close()

#	def choose_file(self):
#		#opens a file browser, showing files only in .hdf format
#		file_name = QtGui.QFileDialog.getOpenFileName(self, "Open Data File", "", "HDF files (*.hdf)")
#		#after the user selected a file, we obtain this filename as a Qstring
#		a=QtCore.QString(file_name)
#		print a
#		#we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
#		self.stacknameedit.setText(str(a))
		
#		#Function choose_file started when  the  open_file of the  Poptwodali window is clicked (same as above but for bdb files(maybe we can combine these two into one function)
#	def choose_file1(self):
#		file_name1 = QtGui.QFileDialog.getOpenFileName(self, "Open Data File", "EMAN2DB/", "BDB FILES (*.bdb)" )
#		a=QtCore.QString(file_name1)
#		b=os.path.basename(str(a))
#		c=os.path.splitext(b)[0]
#		d="bdb:"+c
#		print d
#		self.stacknameedit.setText(d)
	
#	def show_output_info(self):
#		QMessageBox.information(self, "sx* output",'outdir is the name of the output folder specified by the user. If it does not exist, the directory will be created. If it does exist, the program will crash and an error message will come up. Please change the name of directory and restart the program.')
		
class SXTab_main(QWidget):

	def __init__(self, parent = None):
		QWidget.__init__(self, parent)
		
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		self.sxpop = parent
		
		# layout parameters
		self.y1 = 10
		self.y2 = self.y1 + 98
		self.y4 = self.y2 + 450
		self.y5 = self.y4 + 95
		self.yspc = 4
		
		self.x1 = 10
		self.x2 = self.x1 + 200
		self.x3 = self.x2+145
		self.x4 = self.x3+100
		self.x5 = 230
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		
		# # Set the window title
		# self.setWindowTitle(parent.sxscript)
		# Set a label and its position in this tab
		temp_label = QtGui.QLabel('<b>%s</b> - %s' % (parent.sxscript, parent.short_info), self)
		temp_label.move(self.x1,self.y1)
		self.y1 += 50
		
		# Add load paramater button 
		self.load_params_btn = QPushButton("Load parameters", self)
		self.load_params_btn.move(self.x1-5,self.y1)
		self.load_params_btn.setToolTip('Load gui parameter settings to retrieve a previously-saved one')
		self.connect(self.load_params_btn, SIGNAL("clicked()"), parent.load_params)
		
		# Add gui components for editing command args and options
		for cmd_token in parent.cmd_token_list:
			if cmd_token.group == 'sxgroup_main':			
				# Create label 
				gui_label = QtGui.QLabel(cmd_token.label, self)
				gui_label.move(self.x1,self.y2)
				
				# Create gui component and associate it to this cmd_token 
				if cmd_token.widget == "sxwidget_line_edit":
					cmd_token.gui = QtGui.QLineEdit(self)
					cmd_token.gui.setText(cmd_token.default)
				elif cmd_token.widget == "sxwidget_check_box":
					cmd_token.gui = QtGui.QCheckBox("", self)
					cmd_token.gui.setCheckState(cmd_token.default)
				else:
					assert False, 'Unreachable Code!!!'
			
				cmd_token.gui.move(self.x2,self.y2 - 7)
				cmd_token.gui.setToolTip(cmd_token.help)		
	
				self.y2 = self.y2+25

		# Add space
		self.y2 = self.y2+25*1

		# Add gui components for MPI related paramaters if necessary
		temp_label = QtGui.QLabel('MPI processors', self)
		temp_label.move(self.x1,self.y2)
		self.mpi_nproc_edit = QtGui.QLineEdit(self)
		self.mpi_nproc_edit.setText('1')
		self.mpi_nproc_edit.move(self.x2,self.y2)
		self.mpi_nproc_edit.setToolTip('The number of processors to use. Default is single processor mode')
		
		self.y2 =self.y2+25

		temp_label = QtGui.QLabel('MPI Command Line Template', self)
		temp_label.move(self.x1,self.y2)
		self.mpi_cmd_line_edit = QtGui.QLineEdit(self)
		self.mpi_cmd_line_edit.setText('')
		self.mpi_cmd_line_edit.move(self.x2,self.y2)
		self.mpi_cmd_line_edit.setToolTip('The template of MPI command line (e.g. "mpirun -np XXX_SXMPI_NPROC_XXX --host n0,n1,n2 XXX_SXCMD_LINE_XXX"). If empty, use "mpirun -np XXX_SXMPI_NPROC_XXX XXX_SXCMD_LINE_XXX"')

		self.y2 =self.y2+25
		
		# If MPI is not supported, disable this widget
		self.set_widget_enable_state(self.mpi_nproc_edit, parent.mpi_support)
		self.set_widget_enable_state(self.mpi_cmd_line_edit, parent.mpi_support)

		# Add gui components for queue submission (qsub)
		is_qsub_enabled = False
		temp_label = QtGui.QLabel('Submit Job to Queue', self)
		temp_label.move(self.x1,self.y2)
		self.qsub_enable_checkbox = QtGui.QCheckBox("", self)
		self.qsub_enable_checkbox.setCheckState(is_qsub_enabled)
		self.qsub_enable_checkbox.stateChanged.connect(self.set_qsub_enable_state) # To control enable state of the following qsub related widgets
		self.qsub_enable_checkbox.move(self.x2,self.y2)
		self.qsub_enable_checkbox.setToolTip('Submit job to queue')
		
		self.y2 =self.y2+25
		
		temp_label = QtGui.QLabel('Job Name', self)
		temp_label.move(self.x1,self.y2)
		self.qsub_job_name_edit = QtGui.QLineEdit(self)
		self.qsub_job_name_edit.setText(parent.sxscript)
		self.qsub_job_name_edit.move(self.x2,self.y2)
		self.qsub_job_name_edit.setToolTip('Name of this job')

		self.y2 =self.y2+25

		temp_label = QtGui.QLabel('Submission Command', self)
		temp_label.move(self.x1,self.y2)
		self.qsub_cmd_edit = QtGui.QLineEdit(self)
		self.qsub_cmd_edit.setText('qsub')
		self.qsub_cmd_edit.move(self.x2,self.y2)
		self.qsub_cmd_edit.setToolTip('Name of submission command to queue job')

		self.y2 =self.y2+25

		temp_label = QtGui.QLabel('Submission Script Template', self)
		temp_label.move(self.x1,self.y2)
		self.qsub_script_edit = QtGui.QLineEdit(self)
		self.qsub_script_edit.setText('msgui_qsub.sh')
		self.qsub_script_edit.move(self.x2,self.y2)
		self.qsub_script_edit.setToolTip('File name of submission script template (e.g. $EMAN2DIR/bin/msgui_qsub.sh')

		self.y2 =self.y2+25
		
		# Initialize enable state of qsub related widgets
		self.set_qsub_enable_state()
		
		# Add space
		self.y2 = self.y2+25*1

		# Add save paramater button 
		self.save_params_btn = QPushButton("Save Parameters", self)
		# self.save_params_btn.move(self.x1-5,  self.y4)
		self.save_params_btn.move(self.x1-5,  self.y2)
		self.save_params_btn.setToolTip('Save gui parameter settings')
		self.connect(self.save_params_btn, SIGNAL("clicked()"), parent.save_params)
		
		# self.y4 = self.y4+30
		self.y2 = self.y2+30

		self.cmd_line_btn = QPushButton("Generate command line", self)
		# self.cmd_line_btn.move(self.x1-5,  self.y4)
		self.cmd_line_btn.move(self.x1-5,  self.y2)
		self.cmd_line_btn.setToolTip('Generate command line from gui parameter settings')
		self.connect(self.cmd_line_btn, SIGNAL("clicked()"), parent.generate_cmd_line)
		
		self.y2 = self.y2+30

		# Add a run button
		self.execute_btn = QtGui.QPushButton('Run %s' % parent.sxscript, self)
		# make 3D textured push button look
		s = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px}"
		self.execute_btn.setStyleSheet(s)
		# self.execute_btn.move(self.x5,  self.y5)
		self.execute_btn.move(self.x5,  self.y2)
		self.connect(self.execute_btn, SIGNAL("clicked()"), parent.execute_cmd)

	def set_widget_enable_state(self, widget, is_enabled):
		# Set enable state and background color of widget according to enable state
		bg_color = Qt.white
		if is_enabled == False:
			bg_color = Qt.gray
		
		widget.setEnabled(is_enabled)
		widget_palette = widget.palette()
		widget_palette.setColor(widget.backgroundRole(), bg_color)
		widget.setPalette(widget_palette)

	def set_qsub_enable_state(self):
		is_enabled = False
		if self.qsub_enable_checkbox.checkState() == Qt.Checked:
			is_enabled = True
		
		# Set enable state and background color of mpi related widgets
		if self.sxpop.mpi_support:
			self.set_widget_enable_state(self.mpi_cmd_line_edit, not is_enabled)
		
		# Set enable state and background color of qsub related widgets
		self.set_widget_enable_state(self.qsub_job_name_edit, is_enabled)
		self.set_widget_enable_state(self.qsub_cmd_edit, is_enabled)
		self.set_widget_enable_state(self.qsub_script_edit, is_enabled)
		
class SXTab_advance(QWidget):
	def __init__(self, parent = None):
		QWidget.__init__(self)
				
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# layout parameters
		self.y1=10
		self.yspc = 4
		
		self.x1 = 20
		self.x2 = self.x1+280
		self.x3 = self.x2+145
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		
		# Set the window title
		#self.setWindowTitle('%s advanced parameter selection' % parent.sxscript)
		# Set a label and its position in this tab
		title1=QtGui.QLabel('<b>%s</b> - Set advanced parameters' % parent.sxscript, self)
		title1.move(self.x1,self.y1)
		self.y1 = self.y1+25
		
		# Add gui components for editing command args and options
		for cmd_token in parent.cmd_token_list:
			if cmd_token.group == 'sxgroup_advance':
				# Create label 
				gui_label = QtGui.QLabel(cmd_token.label, self)
				gui_label.move(self.x1,self.y1)
		
				# Create gui component and associate it to this cmd_token 
				if cmd_token.widget == "sxwidget_line_edit":
					cmd_token.gui = QtGui.QLineEdit(self)
					cmd_token.gui.setText(cmd_token.default)
				elif cmd_token.widget == "sxwidget_check_box":
					cmd_token.gui = QtGui.QCheckBox("", self)
					cmd_token.gui.setCheckState(cmd_token.default)
				else:
					assert False, 'Unreachable Code!!!'

				cmd_token.gui.move(self.x2,self.y1)
				cmd_token.gui.setToolTip(cmd_token.help)		

				self.y1 = self.y1+25
# sx_end

# ========================================================================================
# sxcter_start
class SXPopup_sxcter(SXPopup):
	def __init__(self):
		SXPopup.__init__(self, "sxcter", 'Automated estimation of CTF parameters with error assessment')

	def get_option_parser(self):
		assert self.sxscript == "sxcter"
		
		# Set parser from sx*.py script
		from msgui_parser_sxcter import main as parser_main
		parser = parser_main(["sxcter.py"])
		
		return parser
# sxcter_end				

# ========================================================================================
# sxwindow_start
class SXPopup_sxwindow(SXPopup):
	def __init__(self):
		SXPopup.__init__(self, "sxwindow", 'Window out particles with known coordinates from a micrograph')

	def get_option_parser(self):
		assert self.sxscript == "sxwindow"
		# Set parser from sx*.py script
		from msgui_parser_sxwindow import main as parser_main
		parser = parser_main(["sxwindow.py"])
		
		return parser
# sxwindow_end				

# ========================================================================================
# sxisac_start
class SXPopup_sxisac(SXPopup):
	def __init__(self):
		SXPopup.__init__(self, "sxisac", 'Perform Iterative Stable Alignment and Clustering (ISAC) on a 2-D image stack')

	def get_option_parser(self):
		assert self.sxscript == "sxisac"
		# Set parser from sx*.py script
		from msgui_parser_sxisac import main as parser_main
		parser = parser_main(["sxisac.py"])
		
		return parser
# sxisac_end				

# ========================================================================================
# sxviper_start
class SXPopup_sxviper(SXPopup):
	def __init__(self):
		SXPopup.__init__(self, "sxviper", 'Use the "common line" algorithm to assign initial values of phi, theta, psi to 2D average projections.')
		
	def get_option_parser(self):
		assert self.sxscript == "sxviper"
		# Set parser from sx*.py script
		from msgui_parser_sxviper import main as parser_main
		parser = parser_main(["sxviper.py"])
		
		return parser
# sxviper_end	

# ========================================================================================
# sxmeridien_start
class SXPopup_sxmeridien(SXPopup):
	def __init__(self):
		SXPopup.__init__(self, "sxmeridien", 'Performs 3D structure refinement.')
		
	def get_option_parser(self):
		assert self.sxscript == "sxmeridien"
		# Set parser from sx*.py script
		from msgui_parser_sxmeridien import main as parser_main
		parser = parser_main(["sxmeridien.py"])
		
		return parser

	def map_gui_to_cmd_line(self):
		#Here we just read in all user inputs in the line edits of the Poptwodali window
		stack       = self.cmd_token_dict["Particle Stack"].gui.text()  # self.stacknameedit.text()
		output      = self.cmd_token_dict["Output Directory"].gui.text()  # self.foldernameedit.text()
		radius      = self.cmd_token_dict["Particle Radius"].gui.text()  # self.radiusedit.text()
		inires      = self.cmd_token_dict["Initial Resolution"].gui.text()  # self.iniresedit.text()
		inivol      = self.cmd_token_dict["Initial 3D Structure"].gui.text()  # self.inivoledit.text()
		sym         = self.cmd_token_dict["Point-Group Symmetry"].gui.text()  # self.symedit.text()
		pwreference = self.cmd_token_dict["Power Spectrum"].gui.text()  # self.pwreferenceedit.text()
		mask        = self.cmd_token_dict["3D Mask"].gui.text()  # self.maskedit.text()
		
		if mask == "None":
			mask = ""
		
		cmd1 = "sxmeridien.py "+str(stack) +" "+ str(output)+" "+ str(inivol)
		
		args = " --radius=" + str(radius) +  " --sym=" + str(sym) + " --smear"

		cmd1 += args
		if len(str(mask)) > 0:
				cmd1 += " --mask3D=" +str(mask)
		if len(str(pwreference)) > 0:
			cmd1 += " --pwreference=" + str(pwreference)
		if len(str(inires)) > 0:
			cmd1 += " --inires=" + str(inires)

		CTF        = self.cmd_token_dict["Use CTF Correction"].gui.checkState() # self.ctfchkbx.checkState()
		# --function=[/home/justus,foo,bar] will result in loading the function bar() from the module foo.py residing in the directory /home/justus.
		# --function=foo
		userf_string = str(self.cmd_token_dict["User Function Name"].gui.text())
		userf_tokens = userf_string.strip().replace("[", "").replace("]", "").split(',')
		if len(userf_tokens) == 1:
			userf    = userf_tokens[0] # self.w1.usrfuncedit.text()
			userfile = "" # self.w1.usrfuncfileedit.text()
		elif len(userf_tokens) == 3:
			userf    = userf_tokens[2] # self.w1.usrfuncedit.text()
			userfile = userf_tokens[0] + "/" + userf_tokens[1] # self.w1.usrfuncfileedit.text()
		else:
			ERROR("msgui_main", "Invalid string format for \"User Function Name\"", 1)
							
		if CTF == Qt.Checked:
			cmd1 = cmd1 + " --CTF"
		if len(userfile) < 1:
			if(len(userf) > 0):   cmd1 += " --function="+str(userf)
		else:
			userfile = str(userfile)
			# break it down into file name and directory path
			rind = userfile.rfind('/')
			
			if rind == -1:
				userfile = os.path.abspath(userfile)
				rind = userfile.rfind('/')

			fname = userfile[rind+1:]
			fname, ext = os.path.splitext(fname)
			fdir = userfile[0:rind]
			cmd1 += " --function=\"[" +fdir+","+fname+","+str(userf)+"]\""
				
		nproc = self.tab_main.mpi_nproc_edit.text() # self.nprocedit.text()
		
		# self.savedparmsdict = {'stackname':str(stack),'foldername':str(output),'inivol':str(inivol),\
		# 						'radius':str(radius),'sym':str(sym),'inires':str(inires),'nproc':str(nproc),'pwreference':str(pwreference),\
		# 						'maskname':str(mask),"ctf":CTF,"usrfunc":str(userf), "usrfuncfile":str(userfile)}
		# 
		# self.w1.savedparmsdict = self.savedparmsdict
		# 
		# if int(str(nproc)) > 0:
		# 		cmd1 = "mpirun -np " + str(nproc) + " " + cmd1
		# else:  ERROR("sxgui","numper of processors has to be specified",1)
		if int(str(nproc)) <= 0:
			ERROR("msgui_main", "number of processors has to be specified", 1)
		
		self.sxcmd_line = cmd1
# sxmeridien_end				

# ========================================================================================
# sx3dvariability_start
class SXPopup_sx3dvariability(SXPopup):
	def __init__(self):
		SXPopup.__init__(self, "sx3dvariability", '3D local variability in the real space.')
		
	def get_option_parser(self):
		assert self.sxscript == "sx3dvariability"
		# Set parser from sx*.py script
		from msgui_parser_sx3dvariability import main as parser_main
		parser = parser_main(["sx3dvariability.py"])
		
		return parser
# sx3dvariability_end	

# ========================================================================================
# sxlocres_start
class SXPopup_sxlocres(SXPopup):
	def __init__(self):
		# SXPopup.__init__(self, "sxlocres", 'Compute local FSC resolution in real space unsing half-maps, within area outlined by the maskfile and within regions wn x wn x wn')
		SXPopup.__init__(self, "sxlocres", 'Compute local FSC resolution in real space unsing half-maps.')
		
	def get_option_parser(self):
		assert self.sxscript == "sxlocres"
		# Set parser from sx*.py script
		from msgui_parser_sxlocres import main as parser_main
		parser = parser_main(["sxlocres.py"])
		
		return parser
		
	def map_gui_to_cmd_line(self):
		#Here we just read in all user inputs in the line edits of the Poptwodali window
		vol1     = self.cmd_token_dict["First Half-set Volume"].gui.text()  # self.vol1edit.text()
		vol2     = self.cmd_token_dict["Second Half-set Volume"].gui.text()  # self.vol2edit.text()
		mask     = self.cmd_token_dict["Mask Volume"].gui.text()  # self.maskedit.text()
		output   = self.cmd_token_dict["Local Resolution Volume"].gui.text()  # self.outputedit.text()
		wn       = self.cmd_token_dict["Window Size in Voxels"].gui.text()  # self.wnedit.text()
		step     = self.cmd_token_dict["Resolution Step in Voxels"].gui.text()  # self.stepedit.text()
		cutoff   = self.cmd_token_dict["Resolution Cutoff"].gui.text()  # self.cutoffedit.text()
		radius   = self.cmd_token_dict["Particule Radius"].gui.text()  # self.radiusedit.text()
		fsc      = self.cmd_token_dict["FSC File Name"].gui.text()  # self.fscedit.text()
		# nproc    = self.cmd_token_dict["Particle Stack"].gui.text()  # self.nprocedit.text()
		
		cmd1 = "sxlocres.py "+str(vol1) +" "+ str(vol2)
		if(len(str(mask)) > 0): cmd1 += " "+ str(mask)
		cmd1 += " "+ str(output)
		# if int(str(nproc)) >1:
		# 		cmd1 = "mpirun -np " + str(nproc) + " " + cmd1
		# cmd1 += " --wn="+ str(wn)+" --step="+ str(step) + " --cutoff="+str(cutoff)+ " --radius="+str(radius)+ " --fsc="+str(fsc)
		cmd1 += " --wn="+ str(wn)+" --step="+ str(step) + " --cutoff="+str(cutoff)+ " --radius="+str(radius)
		if fsc != "None":
			cmd1 += " --fsc="+str(fsc)
			
		# if int(str(nproc)) >1:
		# 	cmd1 += "  --MPI"
		
		# self.savedparmsdict = {'vol1':str(vol1),'vol2':str(vol2),'mask':str(mask),'output':str(output),'wn':str(wn),\
		# 				'step':str(step),'cutoff':cutoff,'radius':str(radius),'fsc':str(fsc),'nproc':nproc}
		
		# if writefile:		
		# 	(fname,stat)= QInputDialog.getText(self,"Generate Command Line","Enter name of file to save command line in",QLineEdit.Normal,"")
		# 	if stat:
		# 		f = open(fname,'a')
		# 		f.write(cmd1)
		# 		f.write('\n')
		# 		f.close()
		# 
		# print cmd1
		# self.cmd = cmd1
		
		self.sxcmd_line = cmd1
# sxlocres_end	

# ========================================================================================
# sxfilterlocal_start
class SXPopup_sxfilterlocal(SXPopup):
	def __init__(self):
		# SXPopup.__init__(self, "sxfilterlocal", 'Locally filter input volume based on values within the associated local resolution volume (sxlocres.py) within area outlined by the maskfile.')
		SXPopup.__init__(self, "sxfilterlocal", 'Locally filter input volume based on local resolution.')
		
	def get_option_parser(self):
		assert self.sxscript == "sxfilterlocal"
		# Set parser from sx*.py script
		from msgui_parser_sxfilterlocal import main as parser_main
		parser = parser_main(["sxfilterlocal.py"])
		
		return parser
# sxfilterlocal_end	

# ========================================================================================
# sxsort3d_start
class SXPopup_sxsort3d(SXPopup):
	def __init__(self):
		# SXPopup.__init__(self, "sxsort3d", 'Sorts out possible conformations from one heterogenous data set whose xform.projection parameters are already determined using K-means, and Equal K-means method.')
		SXPopup.__init__(self, "sxsort3d", 'Sorts out possible conformations from one heterogenous data set.')
		
	def get_option_parser(self):
		assert self.sxscript == "sxsort3d"
		# Set parser from sx*.py script
		from msgui_parser_sxsort3d import main as parser_main
		parser = parser_main(["sxsort3d.py"])
		
		return parser
# sxsort3d_end	

# ========================================================================================
###MAIN WINDOW		(started by class App)
#This class includes the layout of the main window; within each class, i name the main object self, to avoid confusion)			
class MainWindow(QtGui.QWidget):
	def __init__(self, parent=None):
		QtGui.QWidget.__init__(self, parent)
		# self.setStyleSheet('background-image: url("1.png")')
		# Set the title of the window
		self.setWindowTitle('MPI-SPARX GUI (Alpha Version)')
		self.setAutoFillBackground(True)				
		palette = QPalette(self)				
		palette.setBrush(QPalette.Background, QBrush(QPixmap(get_image_directory() + "sxgui.py_main_window_background_image.png")))				
		# palette.setBrush(QPalette.Background, QBrush(QPixmap("Fig6.png")))
		# palette.setBrush(QPalette.Background, QBrush(QPixmap("spaxgui02.png")))
		self.setPalette(palette)
		
		# --------------------------------------------------------------------------------
		# General 
		# --------------------------------------------------------------------------------
		# Add title label and set position and font style
		title=QtGui.QLabel("<span style='font-size:18pt; font-weight:600; color:#aa0000;'><b>PROGRAMS </b></span><span style='font-size:12pt; font-weight:60; color:#aa0000;'>(shift-click for wiki)</span>", self)
		title.move(17,47)
		QtGui.QToolTip.setFont(QtGui.QFont('OldEnglish', 8))

		# Add Push button to display popup window for info about the application
		self.btn_info = QPushButton(self)
		self.connect(self.btn_info, SIGNAL("clicked()"), self.info)
		icon = QIcon(get_image_directory() + "sparxicon.png") # Decorates the button with the sparx image
		self.btn_info.setIcon(icon)
		self.btn_info.move(5, 5)
		self.btn_info.setToolTip('Info Page')

		# Add Close button
		self.btn_quit = QPushButton("Close", self)
		self.btn_quit.setToolTip('Close SPARX GUI ')
		self.btn_quit.move(65, 5)
		self.connect(self.btn_quit, QtCore.SIGNAL('clicked()'),QtGui.qApp, QtCore.SLOT('quit()'))

		# --------------------------------------------------------------------------------
		# Add buttons for SPA related sx*.py processes
		# --------------------------------------------------------------------------------
		# self.y2 = 65
		self.y2 = 95
		
		temp_btn = QPushButton("sxcter", self)
		temp_btn.move(10, self.y2)
		temp_btn.setToolTip('Automated estimation of CTF parameters with error assessment')
		self.connect(temp_btn, SIGNAL("clicked()"), self.sxcter)

		self.y2 += 30

		temp_btn = QPushButton("sxwindow", self)
		temp_btn.move(10, self.y2)
		temp_btn.setToolTip('Window out particles with known coordinates from a micrograph')
		self.connect(temp_btn, SIGNAL("clicked()"), self.sxwindow)

		self.y2 += 30

		temp_btn = QPushButton("sxisac", self)
		temp_btn.move(10, self.y2)
		temp_btn.setToolTip('Perform Iterative Stable Alignment and Clustering (ISAC) on a 2-D image stack')
		self.connect(temp_btn, SIGNAL("clicked()"), self.sxisac)

		self.y2 += 30

		temp_btn = QPushButton("sxviper", self)
		temp_btn.move(10, self.y2)
		temp_btn.setToolTip('Use the "common line" algorithm to assign initial values of phi, theta, psi to 2D average projections.')
		self.connect(temp_btn, SIGNAL("clicked()"), self.sxviper)
		
		self.y2 += 30

		temp_btn = QPushButton("sxmeridien", self)
		temp_btn.move(10, self.y2)
		temp_btn.setToolTip('3D structure refinement')
		self.connect(temp_btn, SIGNAL("clicked()"), self.sxmeridien)

		self.y2 += 30

		temp_btn = QPushButton("sx3dvariability", self)
		temp_btn.move(10, self.y2)
		temp_btn.setToolTip('3D local variability in the real space.')
		self.connect(temp_btn, SIGNAL("clicked()"), self.sx3dvariability)
		
		self.y2 += 30

		temp_btn = QPushButton("sxlocres", self)
		temp_btn.move(10, self.y2)
		temp_btn.setToolTip('Local resolution (FSC)')
		self.connect(temp_btn, SIGNAL("clicked()"), self.sxlocres)

		self.y2 += 30

		temp_btn = QPushButton("sxfilterlocal", self)
		temp_btn.move(10, self.y2)
		temp_btn.setToolTip('3D local filter')
		self.connect(temp_btn, SIGNAL("clicked()"), self.sxfilterlocal)

		self.y2 += 30

		temp_btn = QPushButton("sxsort3d", self)
		temp_btn.move(10, self.y2)
		temp_btn.setToolTip('Sorts out possible conformations from one heterogenous data set whose xform.projection parameters are already determined using K-means, and Equal K-means method.')
		self.connect(temp_btn, SIGNAL("clicked()"), self.sxsort3d)
		
		# --------------------------------------------------------------------------------
		# Window settings
		# --------------------------------------------------------------------------------
		# Set the width and height of the main window
		self.resize(300,400)

	# ====================================================================================
	# Click actions: The following functions are associated with the click event of push buttons (btn##) on the main window. 
	def sxcter(self):
		modifiers = QtGui.QApplication.keyboardModifiers()
		if modifiers == QtCore.Qt.ShiftModifier:
			os.system("python -m webbrowser %ssxcter" % SPARX_DOCUMENTATION_WEBSITE)
			return
			
		self.w = SXPopup_sxcter()
	
	def sxwindow(self):
		modifiers = QtGui.QApplication.keyboardModifiers()
		if modifiers == QtCore.Qt.ShiftModifier:
			os.system("python -m webbrowser %ssxwindow" % SPARX_DOCUMENTATION_WEBSITE)
			return

		self.w = SXPopup_sxwindow()

	def sxisac(self):
		modifiers = QtGui.QApplication.keyboardModifiers()
		if modifiers == QtCore.Qt.ShiftModifier:
			os.system("python -m webbrowser %ssxisac" % SPARX_DOCUMENTATION_WEBSITE)
			return

		self.w = SXPopup_sxisac()

	def sxviper(self):
		modifiers = QtGui.QApplication.keyboardModifiers()
		if modifiers == QtCore.Qt.ShiftModifier:
			os.system("python -m webbrowser %ssxviper" % SPARX_DOCUMENTATION_WEBSITE)
			return
		
		self.w = SXPopup_sxviper()

	def sxmeridien(self):
		modifiers = QtGui.QApplication.keyboardModifiers()
		if modifiers == QtCore.Qt.ShiftModifier:
			os.system("python -m webbrowser %ssxmeridien" % SPARX_DOCUMENTATION_WEBSITE)
			return

		self.w = SXPopup_sxmeridien()

	def sx3dvariability(self):
		modifiers = QtGui.QApplication.keyboardModifiers()
		if modifiers == QtCore.Qt.ShiftModifier:
			os.system("python -m webbrowser %ssx3dvariability" % SPARX_DOCUMENTATION_WEBSITE)
			return

		self.w = SXPopup_sx3dvariability()
		
	def sxlocres(self):
		modifiers = QtGui.QApplication.keyboardModifiers()
		if modifiers == QtCore.Qt.ShiftModifier:
			os.system("python -m webbrowser %ssxlocres" % SPARX_DOCUMENTATION_WEBSITE)
			return
		
		self.w = SXPopup_sxlocres()
		
	def sxfilterlocal(self):
		modifiers = QtGui.QApplication.keyboardModifiers()
		if modifiers == QtCore.Qt.ShiftModifier:
			os.system("python -m webbrowser %ssxfilterlocal" % SPARX_DOCUMENTATION_WEBSITE)
			return
		
		self.w = SXPopup_sxfilterlocal()

	def sxsort3d(self):
		modifiers = QtGui.QApplication.keyboardModifiers()
		if modifiers == QtCore.Qt.ShiftModifier:
			os.system("python -m webbrowser %ssxsort3d" % SPARX_DOCUMENTATION_WEBSITE)
			return
		
		self.w = SXPopup_sxsort3d()
	
	#This is the function info, which is being started when the Pushbutton btn_info of the main window is being clicked
	def info(self):
		#print "Opening a new popup window..."
		#opens the window SXPopup_info, and defines its width and height
		#The layout of the SXPopup_info window is defined in class SXPopup_info(QWidget Window)
		self.w = SXPopup_info()
		self.w.resize(250,200)
		self.w.show()

#  this is the main class of the program
#  Here we provide the necessary imports. The basic GUI widgets are located in QtGui module.
class App(QApplication):
	def __init__(self, *args):
		QApplication.__init__(self, *args)
		#here we define the main window (class MainWindow)
		self.main = MainWindow()
		# self.main.resize(400,450)
		self.main.resize(1000, 755)
		#here we define that when all windows are closed, function byebye of class App will be started
		self.connect(self, SIGNAL("lastWindowClosed()"), self.byebye )
		#hshows main window
		self.main.show()
		#function byebye (just quit)  
	def byebye( self ):
		print' bye bye!'
		self.exit(0)

#  Necessary for execution of the program
def main(args):
	from optparse import OptionParser
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <image>
		
Automatic and manual particle selection. This version is specifically aimed at square boxes
for single particle analysis."""
	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--demo", type="string", default="",   help="Name of the demo whose input parameters should be used to initialize the GUI fields: --demo=mpibdb means the input parameters in demo/mpi_bdb will be used, --demo=mpibdbctf means the input parameters in demo/mpi_bdb_ctf will be used")
	global options
	(options, args) = parser.parse_args()
	global DEMO_mpibdbctf
	DEMO_mpibdbctf = 'mpibdbctf'
	global DEMO_mpibdb
	DEMO_mpibdb = 'mpibdb'
	
	global app
	app = App(args)
	app.setWindowIcon(QIcon(get_image_directory()+"sparxicon.png"))
	
	app.main.show()
	app.main.raise_()
	app.exec_()

if __name__ == "__main__":
	main(sys.argv)
