#!/usr/bin/env python
#
# Author: Toshio Moriya, 11/11/2015 (toshio.moriya@mpi-dortmund.mpg.de)
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPHIRE software packages have some GPL dependencies,
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
from subprocess import *
from functools import partial  # Use to connect event-source widget and event handler
from PyQt4.Qt import *
from PyQt4 import QtGui
from PyQt4 import QtCore
from EMAN2 import *
from EMAN2_cppwrap import *
from global_def import *
from sparx import *

# ========================================================================================
class SXcmd_token:
	def __init__(self):
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.key_base = ""           # key base name of command token (argument or option) in command line
		self.key_prefix = ""         # key prefix of of command token. None for argument, "--" or "-" for option
		self.label = ""              # User friendly name of argument or option
		self.help = ""               # Help info
		self.group = ""              # Tab group: main or advanced
		self.is_required = False     # Required argument or options. No default value are available 
		self.default = ""            # Default value
		self.type = ""               # Type of value
		self.restore = ""            # Restore value
		self.is_in_io = False        # <Used only in wikiparser.py> To check consistency between "usage in command line" and list in "== Input ==" and "== Output ==" sections
		self.restore_widget = None   # <Used only in sxgui.py> Restore widget instance associating with this command token
		self.widget = None           # <Used only in sxgui.py> Widget instance associating with this command token
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		
	def initialize_edit(self, key_base):
		self.key_base = key_base
		self.key_prefix = None
		self.label = None
		self.help = None
		self.group = None
		self.is_required = None
		self.default = None
		self.type = None

# ========================================================================================
class SXcmd:
	def __init__(self, category = "", is_submittable = True):
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.name = ""                        # Name of this command (i.e. name of sx*.py script but without .py extension), used for generating command line
		self.mode = ""                        # key base name of a command token, defining mode/subset of this command. For fullset command, use empty string
		self.label = ""                       # User friendly name of this command
		self.short_info = ""                  # Short description of this command
		self.mpi_support = False              # Flag to indicate if this command suppors MPI version
		self.mpi_add_flag = False             # DESIGN_NOTE: 2015/11/12 Toshio Moriya. This can be removed when --MPI flag is removed from all sx*.py scripts 
		self.category = category              # Category of this command; pipe (pipeline), util (utility)
		self.is_submittable = is_submittable  # External GUI Application (e.g. sxgui_cter.py) should not be submitted to job queue
		self.token_list = []                  # list of command tokens. Need this to keep the order of command tokens
		self.token_dict = {}                  # dictionary of command tokens, organised by key base name of command token. Easy to access a command token but looses their order
		self.btn = None                       # <Used only in sxgui.py> QPushButton button instance associating with this command
		self.widget = None                    # <Used only in sxgui.py> SXCmdWidget instance associating with this command
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		
	def get_mode_name_for(self, target_name):
		mode_name = self.name
		if self.mode != "":
			if target_name in ["file_path"]:
				mode_name = "%s_%s" % (self.name, self.mode)
			elif target_name in ["human"]:
				mode_name = "%s %s%s" % (self.name, self.token_dict[self.mode].key_prefix, self.mode)
				
		return mode_name

# ========================================================================================
def construct_sxcmd_list():
	sxcmd_list = []
	
	# Actual sx commands are inserted into the following section by wikiparser.py.
	# @@@@@ START_INSERTION @@@@@
	# @@@@@ END_INSERTION @@@@@
	
	# Create command token dictionary for each SXcmd instance
	for sxcmd in sxcmd_list:
		for sxcmd_token in sxcmd.token_list:
			# Handle very special cases
			if sxcmd_token.type == "function":
				n_widgets = 2 # function type has two line edit boxes
				sxcmd_token.label = [sxcmd_token.label, "enter name of external file with .py extension containing user function"]
				sxcmd_token.help = [sxcmd_token.help, "(leave blank if file is not external to sphire)"]
				sxcmd_token.default = [sxcmd_token.default, "None"]
				sxcmd_token.restore = sxcmd_token.default
			# else: Do nothing for the other types
			
			# Register this to command token dictionary
			sxcmd.token_dict[sxcmd_token.key_base] = sxcmd_token
		
	return sxcmd_list

# ========================================================================================
class SXconst:
	def __init__(self):
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.key = ""                # <Used only in sxgui.py> key of constant parameter
		self.label = ""              # <Used only in sxgui.py> User friendly name of constant parameter
		self.help = ""               # <Used only in sxgui.py> Help info
		self.register = ""           # <Used only in sxgui.py> Default value
		self.type = ""               # <Used only in sxgui.py> Type of value
		self.register_widget = None  # <Used only in sxgui.py> Restore widget instance associating with this command token
		self.widget = None           # <Used only in sxgui.py> Widget instance associating with this constant parameter
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

# ========================================================================================
class SXconst_set:
	def __init__(self):
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.name = ""               # <Used only in sxgui.py> Name of this constant parameter set
		self.label = ""              # <Used only in sxgui.py> User friendly name of this set
		self.short_info = ""         # <Used only in sxgui.py> Short description of this set
		self.list = []               # <Used only in sxgui.py> list of constant parameters. Need this to keep the order of constant parameters
		self.dict = {}               # <Used only in sxgui.py> dictionary of constant parameters, organised by keys of constant parameters. Easy to access each constant parameter but looses their order
		self.btn = None              # <Used only in sxgui.py> QPushButton button instance associating with this set
		self.window = None           # <Used only in sxgui.py> Widget instance associating with this set
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

# ========================================================================================
def construct_sxconst_set():
	
	# Create the project constant parameter set for project settings
	sxconst_set = SXconst_set(); sxconst_set.name = "Project Settings"; sxconst_set.label = "Project Settings"; sxconst_set.short_info = "Set constant parameter values for this project. These constants will be used as default values of associated arugments and options in command settings. However, the setting here is not required to run commands."
	sxconst = SXconst(); sxconst.key = "protein"; sxconst.label = "protein name"; sxconst.help = "a valid string for file names on your OS."; sxconst.register = "MY_PROTEIN"; sxconst.type = "string"; sxconst_set.list.append(sxconst); sxconst_set.dict[sxconst.key] = sxconst
	sxconst = SXconst(); sxconst.key = "apix"; sxconst.label = "micrograph pixel size [A/pixel]"; sxconst.help = "in Angstrom/pixel"; sxconst.register = "1.0"; sxconst.type = "float"; sxconst_set.list.append(sxconst); sxconst_set.dict[sxconst.key] = sxconst
	sxconst = SXconst(); sxconst.key = "ctfwin"; sxconst.label = "CTF window size [pixel]"; sxconst.help = "in pixel. it should be slightly larger than particle box size"; sxconst.register = "512"; sxconst.type = "int"; sxconst_set.list.append(sxconst); sxconst_set.dict[sxconst.key] = sxconst
	sxconst = SXconst(); sxconst.key = "box"; sxconst.label = "particle box size [pixel]" ; sxconst.help = " in pixel/voxel"; sxconst.register = "0"; sxconst.type = "int"; sxconst_set.list.append(sxconst); sxconst_set.dict[sxconst.key] = sxconst
	sxconst = SXconst(); sxconst.key = "radius"; sxconst.label = "protein particle radius [pixel]"; sxconst.help = "in pixel"; sxconst.register = "0"; sxconst.type = "int"; sxconst_set.list.append(sxconst); sxconst_set.dict[sxconst.key] = sxconst
	sxconst = SXconst(); sxconst.key = "sym"; sxconst.label = "point-group symmetry"; sxconst.help = "e.g. c1, c4, d5"; sxconst.register = "c1"; sxconst.type = "string"; sxconst_set.list.append(sxconst); sxconst_set.dict[sxconst.key] = sxconst
	sxconst = SXconst(); sxconst.key = "mass"; sxconst.label = "protein molecular mass [kDa]"; sxconst.help = "in kDa"; sxconst.register = "0.0"; sxconst.type = "float"; sxconst_set.list.append(sxconst); sxconst_set.dict[sxconst.key] = sxconst
	sxconst = SXconst(); sxconst.key = "config"; sxconst.label = "imaging configrations"; sxconst.help = "a free-style string for your record. please use it to describe the set of imaging configrations used in this project (e.g. types of microscope, detector, enegy filter, abbration corrector, phase plate, and etc."; sxconst.register = "MY_MICROSCOPE"; sxconst.type = "int"; sxconst_set.list.append(sxconst); sxconst_set.dict[sxconst.key] = sxconst
	
	return sxconst_set
	
# ========================================================================================
class SXLookFeelConst:
	# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	# static class variables
	# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	default_bg_color = QColor(195, 195, 230) # Blueish Non-Transparent
	cmd_settings_bg_color = QColor(195, 195, 230, 240) # cmd_settings_bg_color = QColor(195, 195, 230, 175) # Blueish Transparent
	
	grid_margin = 12 # grid_margin = 8
	grid_spacing = 6
	
	# sxcmd_area_min_width = 1080 # Best for MAC OSX
	# sxcmd_area_min_height = 1080 # Best for MAC OSX
	sxcmd_area_min_width = 1140 # Best for Linux
	sxcmd_area_min_height = 1080 # Best for Linux
	
	# menu_area_min_width = 240 # Best for Mac OSX
	menu_area_min_half_width = 120 # Best for Linux
	menu_area_min_width = menu_area_min_half_width * 2 + grid_spacing * 1 # Best for Linux
	menu_area_min_height = sxcmd_area_min_height
	
	# sxconst_set_area_min_width = 440 # Best for Mac OSX
	sxconst_set_area_min_width = 480 # Best for Linux
	sxconst_set_area_min_height = sxcmd_area_min_height
	
	# sxconst_set_window_width = sxconst_set_area_min_width + grid_margin * 5 # Best setting for MAC OSX
	sxconst_set_window_width = sxconst_set_area_min_width + grid_margin * 5 # Best setting for Linux
	sxconst_set_window_height = sxconst_set_area_min_height + grid_margin * 2
	
	# sxmain_window_width = sxcmd_area_min_width + menu_area_min_width + grid_margin * (7 + 1) # Best setting for MAC OSX
	sxmain_window_width = sxcmd_area_min_width + menu_area_min_width + grid_margin * (7 + 7) # Best setting for Linux
	sxmain_window_height = sxcmd_area_min_height + grid_margin * 2
	
	sxutil_window_width = sxmain_window_width
	sxutil_window_height = sxmain_window_height
	
	project_dir = "sxgui_settings"

# ========================================================================================
# Provides all necessary functionarity
# tabs only provides widgets and knows how to layout them
class SXCmdWidget(QWidget):
	def __init__(self, sxconst_set, sxcmd, parent = None):
		super(SXCmdWidget, self).__init__(parent)
		
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.sxconst_set = sxconst_set
		self.sxcmd = sxcmd
		self.sxcmd_tab_main = None
		self.sxcmd_tab_advance = None
		self.child_application_list = []
		
		self.gui_settings_file_path = "%s/gui_settings_%s.txt" % (SXLookFeelConst.project_dir, self.sxcmd.get_mode_name_for("file_path"))
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		
		# Set grid layout
		grid_layout = QGridLayout(self)
		# grid_layout.setMargin(SXLookFeelConst.grid_margin)
		# grid_layout.setSpacing(SXLookFeelConst.grid_spacing)
		
		self.setAutoFillBackground(True)
		palette = QPalette()
		palette.setBrush(QPalette.Background, QBrush(SXLookFeelConst.cmd_settings_bg_color))
		self.setPalette(palette)
		
		self.sxcmd_tab_main = SXCmdTab("Main", self)
		self.sxcmd_tab_advance = SXCmdTab("Advanced", self)
#		self.sxcmd_tab_main.w1 = self.sxcmd_tab_advance
# 		self.tab_widget = QTabWidget(self)
		self.tab_widget = QTabWidget()
		self.tab_widget.insertTab(0, self.sxcmd_tab_main, self.sxcmd_tab_main.name)
		self.tab_widget.insertTab(1, self.sxcmd_tab_advance, self.sxcmd_tab_advance.name)
		# self.tab_widget.setAutoFillBackground(True)
		# palette = self.tab_widget.palette()
		# palette.setBrush(QPalette.Background, QBrush(SXLookFeelConst.cmd_settings_bg_color))
		# self.tab_widget.setPalette(palette)
#		self.tab_widget.resize(880,860) # self.tab_widget.resize(900,1080)
#		self.tab_widget.show()
		grid_layout.addWidget(self.tab_widget, 0, 0)
		
#		# Load the previously saved parameter setting of this sx command
#		if os.path.exists(self.gui_settings_file_path):
#			self.read_params(self.gui_settings_file_path)
	
	def map_widgets_to_sxcmd_line(self):
		# Add program name to command line
		sxcmd_line = "%s.py" % self.sxcmd.name
		
		# Loop through all command tokens
		for sxcmd_token in self.sxcmd.token_list:
			# First, handle very special cases
			if sxcmd_token.type == "function":
				user_func_name_index = 0
				external_file_path_index = 1
				user_func_name = str(sxcmd_token.widget[user_func_name_index].text())
				external_file_path = str(sxcmd_token.widget[external_file_path_index].text())
				
				# This is not default value
				if external_file_path not in ["", sxcmd_token.default[external_file_path_index]]:
					# Case 1: User specified an exteranl function different from default or empty string
					if os.path.splitext(external_file_path)[1] != ".py": 
						QMessageBox.warning(self, "Invalid paramter value", "Exteranl File Path (%s) should include the python script extension (.py)." % (external_file_path))
						return ""
					dir_path, file_basename = os.path.split(external_file_path)
					file_basename = file_basename.replace(".py", "")
					sxcmd_line += " %s%s=[%s,%s,%s]" % (sxcmd_token.key_prefix, sxcmd_token.key_base, dir_path, file_basename, user_func_name)
				elif user_func_name != sxcmd_token.default[user_func_name_index]:
					# Case 2: User specified an internal function different from default
					sxcmd_line += " %s%s=%s" % (sxcmd_token.key_prefix, sxcmd_token.key_base, user_func_name)
				# else: User left default value. Do nothing
			# Then, handle the other cases//
			else:
				if sxcmd_token.type == "bool":
					if not ((sxcmd_token.widget.checkState() == Qt.Checked) == sxcmd_token.default and sxcmd_token.is_required == False): 
						### if (sxcmd_token.widget.checkState() == Qt.Checked) == sxcmd_token.default and sxcmd_token.is_required == True:  # Add this token to command line
						### if (sxcmd_token.widget.checkState() == Qt.Checked) != sxcmd_token.default and sxcmd_token.is_required == True:  # Add this token to command line
						### if (sxcmd_token.widget.checkState() == Qt.Checked) != sxcmd_token.default and sxcmd_token.is_required == False: # Add this token to command line
						sxcmd_line += " %s%s" % (sxcmd_token.key_prefix, sxcmd_token.key_base)
					#else: 
						### if (sxcmd_token.widget.checkState() == Qt.Checked) == sxcmd_token.default and sxcmd_token.is_required == False: # Do not add this token to command line
				else:
					if sxcmd_token.widget.text() == sxcmd_token.default:
						### if sxcmd_token.widget.text() == sxcmd_token.default and sxcmd_token.is_required == True:  # Error case
						if sxcmd_token.is_required == True: 
							QMessageBox.warning(self, "Invalid paramter value", "Token (%s) of command (%s) is required. Please set the value for this." % (sxcmd_token.label, self.sxcmd.get_mode_name_for("message_output")))
							return ""
						### if sxcmd_token.widget.text() == sxcmd_token.default and sxcmd_token.is_required == False: # Do not add this token to command line
						# else: # assert(sxcmd_token.is_required == False) # Do not add to this command line
					else: # sxcmd_token.widget.text() != sxcmd_token.default
						### if sxcmd_token.widget.text() != sxcmd_token.default and sxcmd_token.is_required == True:  # Add this token to command line
						### if sxcmd_token.widget.text() != sxcmd_token.default and sxcmd_token.is_required == False: # Add this token to command line
						
						# For now, using line edit box for the other type
						widget_text = str(sxcmd_token.widget.text())
						if sxcmd_token.type not in ["int", "float", "apix", "wn", "box", "radius"]:
							# Always enclose the string value with single quotes (')
							widget_text = widget_text.strip("\'")  # make sure the string is not enclosed by (')
							widget_text = widget_text.strip("\"")  # make sure the string is not enclosed by (")
							widget_text = "\'%s\'" % (widget_text) # then, enclose the string value with single quotes (')
						
						if sxcmd_token.key_prefix == "":
							sxcmd_line += " %s" % (widget_text)
						elif sxcmd_token.key_prefix == "--":
							sxcmd_line += " %s%s=%s" % (sxcmd_token.key_prefix, sxcmd_token.key_base, widget_text)
						else:
							ERROR("Logical Error: Encountered unexpected prefix for token (%s) of command (%s). Consult with the developer." % (sxcmd_token.key_base, self.sxcmd.get_mode_name_for("human")), "%s in %s" % (__name__, os.path.basename(__file__)))
						# else: # assert(sxcmd_token.widget.text() == sxcmd_token.default) # Do not add to this command line
		
		return sxcmd_line
	
	def generate_cmd_line(self):
		# Generate SX command line 
		sxcmd_line = self.map_widgets_to_sxcmd_line()
		
		if sxcmd_line:
			# SX command line is not empty
			# If mpi is not supported set number of MPI processer (np) to 1
			np = 1
			if self.sxcmd.mpi_support:
				# mpi is supported
				np = int(str(self.sxcmd_tab_main.mpi_nproc_edit.text()))
				# DESIGN_NOTE: 2015/10/27 Toshio Moriya
				# Since we now assume sx*.py exists in only MPI version, always add --MPI flag if necessary
				# This is not elegant but can be removed when --MPI flag is removed from all sx*.py scripts 
				if self.sxcmd.mpi_add_flag:
					sxcmd_line += " --MPI"
					
				# DESIGN_NOTE: 2016/02/11 Toshio Moriya
				# Ideally, the following exceptional cases should not handled in here 
				# because it will remove the generality from the software design
				required_key_base = None
				if self.sxcmd.name == "sxisac":
					required_key_base = "indep_run"
				elif self.sxcmd.name == "sxviper":
					required_key_base = "nruns"
				elif self.sxcmd.name == "sxrviper":
					required_key_base = "n_shc_runs"
				# else: # Do nothing
				
				if required_key_base != None:
					required_divisor = int(str(self.sxcmd.token_dict[required_key_base].widget.text()))
					required_label =  self.sxcmd.token_dict[required_key_base].label
					if required_divisor == 0:
						QMessageBox.warning(self, "Invalid paramter value", "\"%s\" must be larger than 0. Please check the setting" % (required_label))
						return "" 
					
					valid_np = np
					if valid_np % required_divisor != 0:
						if valid_np < required_divisor:
							valid_np = required_divisor
						else:
							valid_np = valid_np - (valid_np % required_divisor)
						QMessageBox.warning(self, "Invalid paramter value", "The number of \"MPI processes\" (%d) is invalid. It MUST BE multiplicity of \"%s\" (%d). Please check the setting. A close valid number is %d." % (np, required_label, required_divisor,valid_np))
						return "" 
							
			# else: assert(np == 1) # because the "MPI Processes" is disabled for sx*.py process which does not support mpi
				
			# Generate command line according to the case
			cmd_line = ""
			if self.sxcmd_tab_main.qsub_enable_checkbox.checkState() == Qt.Checked:
				# Case 1: queue submission is enabled (MPI can be supported or unsupported)
				# Create script for queue submission from a give template
				if os.path.exists(self.sxcmd_tab_main.qsub_script_edit.text()) != True: 
					QMessageBox.warning(self, "Invalid paramter value", "Invalid file path for qsub script template (%s)." % (self.sxcmd_tab_main.qsub_script_edit.text()))
					return "" 
					
				file_template = open(self.sxcmd_tab_main.qsub_script_edit.text(),"r")
				# Extract command line from qsub script template 
				for line in file_template:
					if line.find("XXX_SXCMD_LINE_XXX") != -1:
						cmd_line = line.replace("XXX_SXCMD_LINE_XXX", sxcmd_line)
						if cmd_line.find("XXX_SXMPI_NPROC_XXX") != -1:
							cmd_line = cmd_line.replace("XXX_SXMPI_NPROC_XXX", str(np))
						if cmd_line.find("XXX_SXMPI_JOB_NAME_XXX") != -1:
							cmd_line = cmd_line.replace("XXX_SXMPI_JOB_NAME_XXX", str(self.sxcmd_tab_main.qsub_job_name_edit.text()))
				file_template.close()
			elif self.sxcmd.mpi_support:
				# Case 2: queue submission is disabled, but MPI is supported
				if self.sxcmd_tab_main.qsub_enable_checkbox.checkState() == Qt.Checked: ERROR("Logical Error: Encountered unexpected condition for sxcmd_tab_main.qsub_enable_checkbox.checkState. Consult with the developer.", "%s in %s" % (__name__, os.path.basename(__file__)))
				# Add MPI execution to command line
				cmd_line = str(self.sxcmd_tab_main.mpi_cmd_line_edit.text())
				# If empty string is entered, use a default template
				if cmd_line == "":
					cmd_line = "mpirun -np XXX_SXMPI_NPROC_XXX XXX_SXCMD_LINE_XXX"
				if cmd_line.find("XXX_SXMPI_NPROC_XXX") != -1:
					cmd_line = cmd_line.replace("XXX_SXMPI_NPROC_XXX", str(np))
				if cmd_line.find("XXX_SXCMD_LINE_XXX") != -1:
					cmd_line = cmd_line.replace("XXX_SXCMD_LINE_XXX", sxcmd_line)
			else: 
				# Case 3: queue submission is disabled, and MPI is not supported
				if self.sxcmd_tab_main.qsub_enable_checkbox.checkState() == Qt.Checked: ERROR("Logical Error: Encountered unexpected condition for sxcmd_tab_main.qsub_enable_checkbox.checkState. Consult with the developer.", "%s in %s" % (__name__, os.path.basename(__file__)))
				# Use sx command as it is
				cmd_line = sxcmd_line
		else:
			# SX command line is be empty because an error happens in map_widgets_to_sxcmd_line
			cmd_line = ""
		
		return cmd_line
	
	def execute_cmd_line(self):
		# Generate command line 
		cmd_line = self.generate_cmd_line()
		
		if cmd_line:
			# Command line is not empty
			# First, check existence of outputs
			for sxcmd_token in self.sxcmd.token_list:
				if sxcmd_token.type == "output":
					if os.path.exists(sxcmd_token.widget.text()):
						# DESIGN_NOTE: 2015/11/24 Toshio Moriya
						# This special case needs to be handled with more general method...
						if self.sxcmd.name in ["sxisac", "sxviper", "sxrviper", "sxmeridien", "sxsort3d"]:
							reply = QMessageBox.question(self, "Output Directory/File", "Output Directory/File (%s) already exists. Do you really want to run the program with continue mode?" % (sxcmd_token.widget.text()), QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
							if reply == QMessageBox.No:
								return
							# else: # Do nothing
						else:
							QMessageBox.warning(self, "Output Directory/File", "Output Directory/File (%s) already exists. Please change the name and try it again. Aborting execution ..." % (sxcmd_token.widget.text()))
							return
			
			# If mpi is not supported set number of MPI processer (np) to 1
			np = 1
			if self.sxcmd.mpi_support:
				np = int(str(self.sxcmd_tab_main.mpi_nproc_edit.text()))
			
			if self.sxcmd_tab_main.qsub_enable_checkbox.checkState() == Qt.Checked:
				# Case 1: queue submission is enabled (MPI can be supported or unsupported)
				# Create script for queue submission from a give template
				template_file_path = self.sxcmd_tab_main.qsub_script_edit.text()
				if os.path.exists(template_file_path) == False: 
					QMessageBox.warning(self, "Invalid paramter value", "Invalid file path for qsub script template (%s). Aborting execution ..." % (template_file_path))
					return
				file_template = open(self.sxcmd_tab_main.qsub_script_edit.text(),"r")
				file_name_qsub_script = "qsub_" + str(self.sxcmd_tab_main.qsub_job_name_edit.text()) + ".sh"
				file_qsub_script = open(file_name_qsub_script,"w")
				for line_io in file_template:
					if line_io.find("XXX_SXCMD_LINE_XXX") != -1:
						line_io = cmd_line
					else:
						if line_io.find("XXX_SXMPI_NPROC_XXX") != -1:
							line_io = line_io.replace("XXX_SXMPI_NPROC_XXX", str(np))
						if line_io.find("XXX_SXMPI_JOB_NAME_XXX") != -1:
							line_io = line_io.replace("XXX_SXMPI_JOB_NAME_XXX", str(self.sxcmd_tab_main.qsub_job_name_edit.text()))
					file_qsub_script.write(line_io)
				file_template.close()
				file_qsub_script.close()
				# Generate command line for queue submission
				cmd_line_in_script = cmd_line
				cmd_line = str(self.sxcmd_tab_main.qsub_cmd_edit.text()) + " " + file_name_qsub_script
				print "Wrote the following command line in the queue submission script: "
				print cmd_line_in_script
				print "Submitted a job by the following command: "
				print cmd_line
			else:
				# Case 2: queue submission is disabled (MPI can be supported or unsupported)
				if self.sxcmd_tab_main.qsub_enable_checkbox.checkState() == Qt.Checked: ERROR("Logical Error: Encountered unexpected condition for sxcmd_tab_main.qsub_enable_checkbox.checkState. Consult with the developer.", "%s in %s" % (__name__, os.path.basename(__file__)))
				print "Executed the following command: "
				print cmd_line
		
			# Execute the generated command line
 			process = subprocess.Popen(cmd_line, shell=True)
			self.emit(SIGNAL("process_started"), process.pid)
			if self.sxcmd.is_submittable == False:
				assert(self.sxcmd.mpi_support == False)
				# Register to This is a GUI application
				self.child_application_list.append(process)
			
			# Save the current state of GUI settings
			if os.path.exists(SXLookFeelConst.project_dir) == False:
				os.mkdir(SXLookFeelConst.project_dir)
			self.write_params(self.gui_settings_file_path)
		# else: SX command line is be empty because an error happens in generate_cmd_line. Let's do nothing
	
	def save_cmd_line(self):
		# Generate command line 
		cmd_line = self.generate_cmd_line()
		if cmd_line:
			file_path_out = QFileDialog.getSaveFileName(self, "Generate Command Line", options = QFileDialog.DontUseNativeDialog)
			if file_path_out != "":
				file_out = open(file_path_out,"w")
				file_out.write(cmd_line + "\n")
				file_out.close()
				print "Saved the following command to %s:" % file_path_out
				print cmd_line
				
				# Save the current state of GUI settings
				if os.path.exists(SXLookFeelConst.project_dir) == False:
					os.mkdir(SXLookFeelConst.project_dir)
				self.write_params(self.gui_settings_file_path)
		# else: Do nothing
	
	def write_params(self, file_path_out):
		file_out = open(file_path_out,"w")
		
		# Write script name for consistency check upon loading
		file_out.write("@@@@@ %s gui settings - " % (self.sxcmd.get_mode_name_for("human")))
		# file_out.write(EMANVERSION + " (CVS" + CVSDATESTAMP[6:-2] +")")
		file_out.write(EMANVERSION + " (GITHUB: " + DATESTAMP +")" )
		file_out.write(" @@@@@ \n")
		
		# Define list of (tab) groups
		group_main = "main"
		group_advanced = "advanced"
		
		# Loop through all groups. First write out values of widgets in main tab, then ones in advanced
		for group in [group_main, group_advanced]:
			# Loop through all command tokens
			for cmd_token in self.sxcmd.token_list:
				if cmd_token.group == group:
					# First, handle very special cases
					if cmd_token.type == "function":
						# This type has two line edit boxes as a list of widget
						n_widgets = 2
						for widget_index in xrange(n_widgets):
							val_str = str(cmd_token.widget[widget_index].text()) 
							file_out.write("<%s> %s (default %s) == %s \n" % (cmd_token.key_base, cmd_token.label[widget_index], cmd_token.default[widget_index], val_str))
					# Then, handle the other cases
					else:
						val_str = ""
						if cmd_token.type == "bool":
							if cmd_token.widget.checkState() == Qt.Checked:
								val_str = "YES"
							else:
								val_str = "NO"
						else:
							# The other type has only one line edit box
							val_str = str(cmd_token.widget.text())
						
						if cmd_token.is_required == False:
							file_out.write("<%s> %s (default %s) == %s \n" % (cmd_token.key_base, cmd_token.label, cmd_token.default, val_str))
						else:
							file_out.write("<%s> %s (default required %s) == %s \n" % (cmd_token.key_base, cmd_token.label, cmd_token.type, val_str))
				# else: do nothig
			
		# At the end of parameter file...
		# Write MPI parameters 
		file_out.write("%s == %s \n" % ("MPI processors", str(self.sxcmd_tab_main.mpi_nproc_edit.text())))
		file_out.write("%s == %s \n" % ("MPI Command Line Template", str(self.sxcmd_tab_main.mpi_cmd_line_edit.text())))
		# Write Qsub paramters 
		if self.sxcmd_tab_main.qsub_enable_checkbox.checkState() == Qt.Checked:
			val_str = "YES"
		else:
			val_str = "NO"
		file_out.write("%s == %s \n" % ("Submit Job to Queue", val_str))	
		file_out.write("%s == %s \n" % ("Job Name", str(self.sxcmd_tab_main.qsub_job_name_edit.text())))
		file_out.write("%s == %s \n" % ("Submission Command", str(self.sxcmd_tab_main.qsub_cmd_edit.text())))
		file_out.write("%s == %s \n" % ("Submission Script Template", str(self.sxcmd_tab_main.qsub_script_edit.text())))
		
		file_out.close()
	
	def read_params(self, file_path_in):
		file_in = open(file_path_in,"r")
	
		# Check if this parameter file is for this sx script
		line_in = file_in.readline()
		if line_in.find("@@@@@ %s gui settings" % (self.sxcmd.get_mode_name_for("human"))) != -1:
			n_function_type_lines = 2
			function_type_line_counter = 0
			# loop through the rest of lines
			for line_in in file_in:
				# Extract label (which should be left of "=="). Also strip the ending spaces
				label_in = line_in.split("==")[0].strip()
				# Extract value (which should be right of "=="). Also strip all spaces
				val_str_in = line_in.split("==")[1].strip() 
				
				if label_in == "MPI processors":
					self.sxcmd_tab_main.mpi_nproc_edit.setText(val_str_in)
				elif label_in == "MPI Command Line Template":
					self.sxcmd_tab_main.mpi_cmd_line_edit.setText(val_str_in)
				elif label_in == "Submit Job to Queue":
					if val_str_in == "YES":
						self.sxcmd_tab_main.qsub_enable_checkbox.setChecked(True)
					else: # assert(val_str_in == "NO")
						self.sxcmd_tab_main.qsub_enable_checkbox.setChecked(False)
				elif label_in == "Job Name":
					self.sxcmd_tab_main.qsub_job_name_edit.setText(val_str_in)
				elif label_in == "Submission Command":
					self.sxcmd_tab_main.qsub_cmd_edit.setText(val_str_in)
				elif label_in == "Submission Script Template":
					self.sxcmd_tab_main.qsub_script_edit.setText(val_str_in)
				else:
					# Extract key_base of this command token
					target_operator = "<"
					item_tail = label_in.find(target_operator)
					if item_tail != 0: 
						QMessageBox.warning(self, "Invalid Paramter File Format", "Command token entry should start from \"%s\" for key base name in line (%s). The format of this file might be corrupted. Please save the paramater file again." % (target_operator, line_in))
					label_in = label_in[item_tail + len(target_operator):].strip() # Get the rest of line
					target_operator = ">"
					item_tail = label_in.find(target_operator)
					if item_tail == -1: 
						QMessageBox.warning(self, "Invalid Paramter File Format", "Command token entry should have \"%s\" closing key base name in line (%s) The format of this file might be corrupted. Please save the paramater file again." % (target_operator, line_in))
					key_base = label_in[0:item_tail]
					# Get corresponding cmd_token
					if key_base not in self.sxcmd.token_dict.keys(): 
						QMessageBox.warning(self, "Invalid Paramter File Format", "Invalid base name of command token \"%s\" is found in line (%s). This paramter file might be imcompatible with the current version. Please save the paramater file again." % (key_base, line_in))
					cmd_token = self.sxcmd.token_dict[key_base]
					# First, handle very special cases
					if cmd_token.type == "function":
						cmd_token.widget[function_type_line_counter].setText(val_str_in)
						function_type_line_counter += 1
						function_type_line_counter %= n_function_type_lines # function have two line edit boxes
					# Then, handle the other cases
					else:
						if cmd_token.type == "bool":
							# construct new widget(s) for this command token
							if val_str_in == "YES":
								cmd_token.widget.setChecked(Qt.Checked)
							else: # val_str_in == "NO"
								cmd_token.widget.setChecked(Qt.Unchecked)
						else:
							# For now, use line edit box for the other type
							cmd_token.widget.setText(val_str_in)
						
		else:
			QMessageBox.warning(self, "Fail to load paramters", "The specified file is not paramter file for %s." % self.sxcmd.get_mode_name_for("human"))
		
		file_in.close()
	
	def save_params(self):
		file_path_out = str(QFileDialog.getSaveFileName(self, "Save Parameters", options = QFileDialog.DontUseNativeDialog))
		if file_path_out != "":
			self.write_params(file_path_out)
	
	def load_params(self):
		file_path_in = str(QFileDialog.getOpenFileName(self, "Load parameters", options = QFileDialog.DontUseNativeDialog))
		if file_path_in != "":
			self.read_params(file_path_in)
	
	def select_file(self, target_widget, file_format = ""):
		file_path = ""
		if file_format == "bdb":
			file_path = str(QFileDialog.getOpenFileName(self, "Select BDB File", "", "BDB files (*.bdb)", options = QFileDialog.DontUseNativeDialog))
			# Use relative path. 
			if file_path:
				file_path = "bdb:./" + os.path.relpath(file_path).replace("EMAN2DB/", "#").replace(".bdb", "")
				file_path = file_path.replace("/#", "#")
				# If the input directory is the current directory, use the simplified DBD file path format
				if file_path.find(".#") != -1:
					file_path = file_path.replace(".#", "")
		elif file_format == "py":
			file_path = str(QFileDialog.getOpenFileName(self, "Select Python File", "", "PY files (*.py)", options = QFileDialog.DontUseNativeDialog))
			# Use full path
		elif file_format == "pdb":
			file_path = str(QFileDialog.getOpenFileName(self, "Select PDB File", "", "PDB files (*.pdb *.pdb1)", options = QFileDialog.DontUseNativeDialog))
			# Use relative path. 
			if file_path:
				file_path = os.path.relpath(file_path)
		else:
			if file_format:
				file_path = str(QFileDialog.getOpenFileName(self, "Select %s File" % (file_format.upper()), "", "%s files (*.%s)"  % (file_format.upper(), file_format), options = QFileDialog.DontUseNativeDialog))
			else:
				file_path = str(QFileDialog.getOpenFileName(self, "Select File", "", "All files (*.*)", options = QFileDialog.DontUseNativeDialog))
			# Use relative path. 
			if file_path:
				file_path = os.path.relpath(file_path)
			
		if file_path != "":
			target_widget.setText(file_path)
	
	def select_dir(self, target_widget):
		dir_path = str(QFileDialog.getExistingDirectory(self, "Select Directory", "", options = QFileDialog.ShowDirsOnly | QFileDialog.DontResolveSymlinks | QFileDialog.DontUseNativeDialog))
		if dir_path != "":
			# Use relative path. 
			target_widget.setText(os.path.relpath(dir_path))
	
	def quit_all_child_applications(self):
		# Quit all child applications
		for child_application in self.child_application_list:
			child_application.kill()
			# child_application.terminate() # This call ends up outputing "Program interrupted" Message and it is not pretty...
	
	"""
#	def show_output_info(self):
#		QMessageBox.information(self, "sx* output","outdir is the name of the output folder specified by the user. If it does not exist, the directory will be created. If it does exist, the program will crash and an error message will come up. Please change the name of directory and restart the program.")
	"""

# ========================================================================================
class SXCmdTab(QWidget):
	def __init__(self, name, parent=None):
		super(SXCmdTab, self).__init__(parent)
		
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.name = name
		self.sxcmdwidget = parent
		
#		# layout parameters
#		self.y1 = 10
#		# self.y2 = self.y1 + 95 #self.y2 = self.y1 + 98
#		
#		self.x1 = 10
#		self.x2 = self.x1 + 500 # self.x2 = self.x1 + 200
#		self.x3 = self.x2 + 135
#		self.x4 = self.x3 + 100
#		self.x5 = 230
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# local constants
		required_cmd_token_restore_tooltip = "please enter the value manually"
		const_cmd_token_restore_tooltip = "retrieve the registed constant value for this parameter"
		default_cmd_token_restore_tooltip = "retrieve this default value"
		
		# Set grid layout
		grid_row_origin = 0; grid_col_origin = 0
		title_row_span = 1; title_col_span = 2
		short_info_row_span = 1; short_info_col_span = 5
		func_btn_row_span = 1; func_btn_col_span = 2
		token_label_row_span = 1; token_label_col_span = 4
		token_widget_row_span = 1; token_widget_col_span = 1
		cmd_frame_row_span = 32; cmd_frame_col_span = 7
		
		title_label_min_width = 150
		title_label_min_height = 80
		short_info_min_width = 360
		short_info_min_height = 80
		func_btn_min_width = 150
		# token_label_min_width = 460 # Best for MAC OSX
		token_label_min_width = 560 # Best for Linux
		token_widget_min_width = 120
		
		grid_layout = QGridLayout(self)
		grid_layout.setMargin(SXLookFeelConst.grid_margin)
		grid_layout.setSpacing(SXLookFeelConst.grid_spacing)
		grid_layout.setColumnMinimumWidth(grid_col_origin + token_label_col_span, token_widget_min_width)
		grid_layout.setColumnMinimumWidth(grid_col_origin + token_label_col_span + token_widget_col_span, token_widget_min_width)
		grid_layout.setColumnMinimumWidth(grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_min_width)
		grid_layout.setColumnMinimumWidth(grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_min_width)
		# # Give the columns of token label a higher priority to stretch relative to the others
		# for col_span in xrange(token_label_col_span):
		# 	grid_layout.setColumnStretch(grid_row_origin + col_span, grid_layout.columnStretch(grid_row_origin+col_span) + 1)
		
		# Define the tab frame within the tab layout
#		tab_frame = QFrame(self)
		tab_frame = QFrame()
		# tab_frame.setFrameStyle(QFrame.StyledPanel)
		grid_layout.addWidget(tab_frame, grid_row_origin, grid_col_origin, cmd_frame_row_span, cmd_frame_col_span)
		
		# Start add command token widgets to the grid layout
		grid_row = grid_row_origin
		
		tab_group = self.name.lower()
		if tab_group == "main":
			# Set a label and its position in this tab
#			temp_label = QLabel("<b>%s</b>" % (self.sxcmdwidget.sxcmd.get_mode_name_for("human")), self)
			temp_label = QLabel("<b>%s</b>" % (self.sxcmdwidget.sxcmd.get_mode_name_for("human")))
#			temp_label.move(self.x1, self.y1)
			temp_label.setMinimumWidth(title_label_min_width)
			temp_label.setMinimumHeight(title_label_min_height)
			grid_layout.addWidget(temp_label, grid_row, grid_col_origin, title_row_span, title_col_span)
			
			# NOTE: 2015/11/17 Toshio Moriya
			# Necessary to separate "<b>%s</b>" from the information for avoiding to invoke the tag interpretations of string
			# e.g. < becomes the escape character
#			temp_label = QLabel("%s" % (self.sxcmdwidget.sxcmd.short_info), self)
			temp_label = QLabel("%s" % (self.sxcmdwidget.sxcmd.short_info))
			temp_label.setWordWrap(True)
			temp_label.setMinimumWidth(short_info_min_width)
			temp_label.setMinimumHeight(short_info_min_height)
#			temp_label.setFixedWidth(600)
#			temp_label.setFixedHeight(80)
#			temp_label.move(self.x1 + 100, self.y1)
			grid_layout.addWidget(temp_label, grid_row, grid_col_origin + title_col_span, short_info_row_span, short_info_col_span)
#			self.y1 += 50
			
			grid_row += short_info_row_span
			
			# Add load paramaters button 
#			self.load_params_btn = QPushButton("Load parameters", self)
			self.load_params_btn = QPushButton("Load parameters")
#			self.load_params_btn.move(self.x1 - 5, self.y1)
			self.load_params_btn.setMinimumWidth(func_btn_min_width)
			self.load_params_btn.setToolTip("load gui parameter settings to retrieve a previously-saved one")
			self.connect(self.load_params_btn, SIGNAL("clicked()"), self.sxcmdwidget.load_params)
			grid_layout.addWidget(self.load_params_btn, grid_row, grid_col_origin, func_btn_row_span, func_btn_col_span)
#			self.y1 += 25
			
		elif tab_group == "advanced":
			# Set a label and its position in this tab
#			temp_label = QLabel("<b>%s</b>" % (self.sxcmdwidget.sxcmd.get_mode_name_for("human")), self)
			temp_label = QLabel("<b>%s</b>" % (self.sxcmdwidget.sxcmd.get_mode_name_for("human")))
#			temp_label.move(self.x1, self.y1)
			temp_label.setMinimumWidth(title_label_min_width)
			temp_label.setMinimumHeight(title_label_min_height)
			grid_layout.addWidget(temp_label, grid_row, grid_col_origin, title_row_span, title_col_span)
			
			temp_label = QLabel("Set advanced parameters", self)
			temp_label.setWordWrap(True)
			temp_label.setMinimumWidth(short_info_min_width)
			temp_label.setMinimumHeight(short_info_min_height)
#			temp_label.setFixedWidth(600)
#			temp_label.move(self.x1 + 100, self.y1)
			grid_layout.addWidget(temp_label, grid_row, grid_col_origin + title_col_span, short_info_row_span, short_info_col_span)
#			self.y1 += 25
		
		# Add space
#		self.y1 = self.y1 + 25 * 1
		grid_row += 2
		
		# Add widget for editing command args and options
		for cmd_token in self.sxcmdwidget.sxcmd.token_list:
			if cmd_token.group == tab_group:
				
				# First, handle very special cases
				if cmd_token.type == "function":
					n_widgets = 2 # function type has two line edit boxes
					cmd_token_widget = [None] * n_widgets
					cmd_token_restore_widget = [None] * n_widgets
					
					# Create widgets for user function name
					widget_index = 0
#					temp_label = QLabel(cmd_token.label[widget_index], self)
					temp_label = QLabel(cmd_token.label[widget_index])
#					temp_label.move(self.x1, self.y1)
					temp_label.setMinimumWidth(token_label_min_width)
					grid_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)
					
					assert(cmd_token.is_required == False)
					cmd_token_restore_widget[widget_index] = QPushButton("[%s]" % cmd_token.restore[widget_index])
					cmd_token_restore_widget[widget_index].setToolTip(default_cmd_token_restore_tooltip)
					grid_layout.addWidget(cmd_token_restore_widget[widget_index], grid_row, grid_col_origin + token_label_col_span, token_widget_row_span, token_widget_col_span)
					
					# cmd_token_widget[widget_index] = QLineEdit(self)
					cmd_token_widget[widget_index] = QLineEdit()
					cmd_token_widget[widget_index].setText(cmd_token.restore[widget_index])
#					cmd_token_widget[widget_index].move(self.x2,self.y1 - 7)
#					cmd_token_widget[widget_index].setMinimumWidth(token_widget_min_width)
					cmd_token_widget[widget_index].setToolTip(cmd_token.help[widget_index])
					grid_layout.addWidget(cmd_token_widget[widget_index], grid_row, grid_col_origin + token_label_col_span + token_widget_col_span, token_widget_row_span, token_widget_col_span)
					
					self.connect(cmd_token_restore_widget[widget_index], SIGNAL("clicked()"), partial(self.handle_restore_widget_event, cmd_token, widget_index))
					
#					self.y1 = self.y1 + 25
					grid_row +=  1
					
					# Create widgets for external file path containing above user function
					widget_index = 1
#					temp_label = QLabel(cmd_token.label[widget_index], self)
					temp_label = QLabel(cmd_token.label[widget_index])
#					temp_label.move(self.x1, self.y1)
#					temp_label.setMinimumWidth(token_label_min_width)
					grid_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)
					
					assert(cmd_token.is_required == False)
					cmd_token_restore_widget[widget_index] = QPushButton("[%s]" % cmd_token.restore[widget_index])
					cmd_token_restore_widget[widget_index].setToolTip(default_cmd_token_restore_tooltip)
					grid_layout.addWidget(cmd_token_restore_widget[widget_index], grid_row, grid_col_origin + token_label_col_span, token_widget_row_span, token_widget_col_span)
					
#					cmd_token_widget[widget_index] = QLineEdit(self)
					cmd_token_widget[widget_index] = QLineEdit()
					cmd_token_widget[widget_index].setText(cmd_token.restore[widget_index]) # Because default user functions is internal
#					cmd_token_widget[widget_index].move(self.x2,self.y1 - 7)
#					cmd_token_widget[widget_index].setMinimumWidth(token_widget_min_width)
					cmd_token_widget[widget_index].setToolTip(cmd_token.help[widget_index])
					grid_layout.addWidget(cmd_token_widget[widget_index], grid_row, grid_col_origin + token_label_col_span + token_widget_col_span, token_widget_row_span, token_widget_col_span)
					
					self.connect(cmd_token_restore_widget[widget_index], SIGNAL("clicked()"), partial(self.handle_restore_widget_event, cmd_token, widget_index))
					
					file_format = "py"
#					temp_btn = QPushButton("Select Script", self)
					temp_btn = QPushButton("Select Script")
#					temp_btn.move(self.x3, self.y1 - 10)
#					temp_btn.setMinimumWidth(token_widget_min_width)
					temp_btn.setToolTip("display open file dailog to select .%s python script file" % file_format)
					grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
					self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget[widget_index], file_format))
					
#					spacer_frame = QFrame()
#					spacer_frame.setMinimumWidth(token_widget_min_width)
#					grid_layout.addWidget(spacer_frame, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
					
					grid_row +=  1
					
#					temp_label = QLabel(cmd_token.help[widget_index], self)
					temp_label = QLabel(cmd_token.help[widget_index])
#					temp_label.move(self.x1, self.y1 + 25)
#					temp_label.setMinimumWidth(token_label_min_width)
					grid_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)
					
#					self.y1 = self.y1 + 25 * 2
					grid_row +=  1
					
				# Then, handle the other cases
				else:
					# Create label widget 
#					temp_label = QLabel(cmd_token.label, self)
					temp_label = QLabel(cmd_token.label)
#					temp_label.move(self.x1, self.y1)
					temp_label.setMinimumWidth(token_label_min_width)
					grid_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)
					
					# Create widget and associate it to this cmd_token
					cmd_token_widget = None
					cmd_token_restore_widget = None
					cmd_token_restore_tooltip = default_cmd_token_restore_tooltip
					if cmd_token.type == "bool":
						btn_name = "NO"
						is_btn_enable = True
						custom_style = "QPushButton {color:gray; }"
						if cmd_token.restore:
							btn_name = "YES"
						if cmd_token.type in parent.sxconst_set.dict.keys():
							custom_style = "QPushButton {color:green; }"
							cmd_token_restore_tooltip = const_cmd_token_restore_tooltip
						elif cmd_token.is_required:
							btn_name = "required"
							custom_style = "QPushButton {color:red; }"
							is_btn_enable = False
							cmd_token_restore_tooltip = required_cmd_token_restore_tooltip
						cmd_token_restore_widget = QPushButton("[%s]" % btn_name)
						cmd_token_restore_widget.setStyleSheet(custom_style)
						cmd_token_restore_widget.setEnabled(is_btn_enable)
						grid_layout.addWidget(cmd_token_restore_widget, grid_row, grid_col_origin + token_label_col_span, token_widget_row_span, token_widget_col_span)
						
						# construct new widget(s) for this command token
#						cmd_token_widget = QCheckBox("", self)
						cmd_token_widget = QCheckBox("")
						if cmd_token.restore == True:
							cmd_token_widget.setCheckState(Qt.Checked)
						else:
							cmd_token_widget.setCheckState(Qt.Unchecked)
						cmd_token_widget.setEnabled(is_btn_enable)
						grid_layout.addWidget(cmd_token_widget, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span, token_widget_row_span, token_widget_col_span)
						
						self.connect(cmd_token_restore_widget, SIGNAL("clicked()"), partial(self.handle_restore_widget_event, cmd_token))
						
					else:
						btn_name = "%s" % cmd_token.restore
						custom_style = "QPushButton {color:gray; }"
						is_btn_enable = True
						if cmd_token.type in parent.sxconst_set.dict.keys():
							custom_style = "QPushButton {color:green; }"
							cmd_token_restore_tooltip = const_cmd_token_restore_tooltip
						elif cmd_token.is_required:
							btn_name = "required"
							custom_style = "QPushButton {color:red; }"
							is_btn_enable = False
							cmd_token_restore_tooltip = required_cmd_token_restore_tooltip
						cmd_token_restore_widget = QPushButton("[%s]" % btn_name)
						cmd_token_restore_widget.setStyleSheet(custom_style)
						cmd_token_restore_widget.setEnabled(is_btn_enable)
						grid_layout.addWidget(cmd_token_restore_widget, grid_row, grid_col_origin + token_label_col_span, token_widget_row_span, token_widget_col_span)
						
#						cmd_token_widget = QLineEdit(self)
						cmd_token_widget = QLineEdit()
						cmd_token_widget.setText(cmd_token.restore)
#						cmd_token_widget.move(self.x2,self.y1 - 7)
#						cmd_token_widget.setMinimumWidth(token_widget_min_width)
						grid_layout.addWidget(cmd_token_widget, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span, token_widget_row_span, token_widget_col_span)
						
						self.connect(cmd_token_restore_widget, SIGNAL("clicked()"), partial(self.handle_restore_widget_event, cmd_token))
						
						if cmd_token.type == "image":
							file_format = "hdf"
#							temp_btn = QPushButton("Select .%s" % file_format, self)
							temp_btn = QPushButton("Select .%s" % file_format)
#							temp_btn.move(self.x3, self.y1 - 12)
#							temp_btn.setMinimumWidth(token_widget_min_width)
							temp_btn.setToolTip("display open file dailog to select .%s format image file" % file_format)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "bdb"
#							temp_btn = QPushButton("Select .%s" % file_format, self)
							temp_btn = QPushButton("Select .%s" % file_format)
#							temp_btn.move(self.x4, self.y1 - 12)
#							temp_btn.setMinimumWidth(token_widget_min_width)
							temp_btn.setToolTip("display open file dailog to select .%s format image file" % file_format)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
						elif cmd_token.type == "any_image":
#							temp_btn = QPushButton("Select Image", self)
							temp_btn = QPushButton("Select Image")
#							temp_btn.move(self.x3, self.y1 - 12)
#							temp_btn.setMinimumWidth(token_widget_min_width)
							temp_btn.setToolTip("display open file dailog to select standard format image file (e.g. .hdf, .mrc)")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget))
							file_format = "bdb"
#							temp_btn = QPushButton("Select .%s" % file_format, self)
							temp_btn = QPushButton("Select .%s" % file_format)
#							temp_btn.move(self.x4 + 40, self.y1 - 12)
#							temp_btn.setMinimumWidth(token_widget_min_width)
							temp_btn.setToolTip("display open file dailog to select .%s format image file" % file_format)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
						elif cmd_token.type == "bdb":
							file_format = "bdb"
#							temp_btn = QPushButton("Select .%s" % file_format, self)
							temp_btn = QPushButton("Select .%s" % file_format)
#							temp_btn.move(self.x3 + 40, self.y1 - 12)
#							temp_btn.setMinimumWidth(token_widget_min_width)
							temp_btn.setToolTip("display open file dailog to select .%s format image file" % file_format)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
#							spacer_frame = QFrame()
#							spacer_frame.setMinimumWidth(token_widget_min_width)
#							grid_layout.addWidget(spacer_frame, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "pdb":
							file_format = "pdb"
#							temp_btn = QPushButton("Select .%s" % file_format, self)
							temp_btn = QPushButton("Select .%s" % file_format)
#							temp_btn.move(self.x3, self.y1 - 12)
#							temp_btn.setMinimumWidth(token_widget_min_width)
							temp_btn.setToolTip("display open file dailog to select .%s format image file" % file_format)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span* 2, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
#							spacer_frame = QFrame()
#							spacer_frame.setMinimumWidth(token_widget_min_width)
#							grid_layout.addWidget(spacer_frame, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "parameters":
#							temp_btn = QPushButton("Select Paramter", self)
							temp_btn = QPushButton("Select Paramter")
#							temp_btn.move(self.x3, self.y1 - 12)
#							temp_btn.setMinimumWidth(token_widget_min_width)
							temp_btn.setToolTip("display open file dailog to select paramter file (e.g. .txt)")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget))
#							spacer_frame = QFrame()
#							spacer_frame.setMinimumWidth(token_widget_min_width)
#							grid_layout.addWidget(spacer_frame, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "any_file":
							temp_btn = QPushButton("Select File")
							temp_btn.setToolTip("display open file dailog to select file (e.g. *.*)")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget))
						elif cmd_token.type == "directory":
#							temp_btn = QPushButton("Select directory", self)
							temp_btn = QPushButton("Select directory")
#							temp_btn.move(self.x3, self.y1 - 12)
#							temp_btn.setMinimumWidth(token_widget_min_width)
							temp_btn.setToolTip("display select directory dailog")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_dir, cmd_token_widget))
#							spacer_frame = QFrame()
#							spacer_frame.setMinimumWidth(token_widget_min_width)
#							grid_layout.addWidget(spacer_frame, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						# elif cmd_token.type == "output":
#							spacer_frame = QFrame()
#							spacer_frame.setMinimumWidth(token_widget_min_width)
#							grid_layout.addWidget(spacer_frame, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
#							spacer_frame = QFrame()
#							spacer_frame.setMinimumWidth(token_widget_min_width)
#							grid_layout.addWidget(spacer_frame, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						# else:
						# 	if cmd_token.type not in ["int", "float", "string", "apix", "wn", "box", "radius", "sym"]: ERROR("Logical Error: Encountered unsupported type (%s). Consult with the developer."  % cmd_token.type, "%s in %s" % (__name__, os.path.basename(__file__)))
#							spacer_frame = QFrame()
#							spacer_frame.setMinimumWidth(token_widget_min_width)
#							grid_layout.addWidget(spacer_frame, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
#							spacer_frame = QFrame()
#							spacer_frame.setMinimumWidth(token_widget_min_width)
#							grid_layout.addWidget(spacer_frame, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
							
					cmd_token_widget.setToolTip(cmd_token.help)
					cmd_token_restore_widget.setToolTip(cmd_token_restore_tooltip)
					
#					self.y1 = self.y1 + 25
					grid_row += 1
				
				# Register this widget
				cmd_token.widget = cmd_token_widget
				cmd_token.restore_widget = cmd_token_restore_widget
		
		if tab_group == "main":
			# Add space
#			self.y1 = self.y1 + 25 * 1
			grid_row += 1
			
			# Add gui components for MPI related paramaters
#			temp_label = QLabel("MPI processors", self)
			temp_label = QLabel("MPI processors")
#			temp_label.move(self.x1, self.y1)
#			temp_label.setMinimumWidth(token_label_min_width)
			grid_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)
			
			# self.mpi_nproc_edit = QLineEdit(self)
			self.mpi_nproc_edit = QLineEdit()
			self.mpi_nproc_edit.setText("1")
#			self.mpi_nproc_edit.move(self.x2, self.y1)
#			self.mpi_nproc_edit.setMinimumWidth(token_widget_min_width)
			self.mpi_nproc_edit.setToolTip("number of processors to use. default is single processor mode")
			grid_layout.addWidget(self.mpi_nproc_edit, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span, token_widget_row_span, token_widget_col_span)
			
			# self.y1 = self.y1 + 25
			grid_row += 1
			
#			temp_label = QLabel("MPI command line template", self)
			temp_label = QLabel("MPI command line template")
#			temp_label.move(self.x1, self.y1)
#			temp_label.setMinimumWidth(token_label_min_width)
			grid_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)
			
#			self.mpi_cmd_line_edit = QLineEdit(self)
			self.mpi_cmd_line_edit = QLineEdit()
			self.mpi_cmd_line_edit.setText("")
#			self.mpi_cmd_line_edit.move(self.x2, self.y1)
#			self.mpi_cmd_line_edit.setMinimumWidth(token_widget_min_width)
			self.mpi_cmd_line_edit.setToolTip("template of MPI command line (e.g. \"mpirun -np XXX_SXMPI_NPROC_XXX --host n0,n1,n2 XXX_SXCMD_LINE_XXX\"). if empty, use \"mpirun -np XXX_SXMPI_NPROC_XXX XXX_SXCMD_LINE_XXX\"")
			grid_layout.addWidget(self.mpi_cmd_line_edit, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span, token_widget_row_span, token_widget_col_span)
			
#			self.y1 = self.y1 + 25
			grid_row += 1
			
			# If MPI is not supported, disable this widget
			self.set_text_entry_widget_enable_state(self.mpi_nproc_edit, self.sxcmdwidget.sxcmd.mpi_support)
			self.set_text_entry_widget_enable_state(self.mpi_cmd_line_edit, self.sxcmdwidget.sxcmd.mpi_support)
			
			# Add gui components for queue submission (qsub)
			is_qsub_enabled = False
#			temp_label = QLabel("submit job to queue", self)
			temp_label = QLabel("submit job to queue")
#			temp_label.move(self.x1, self.y1)
#			temp_label.setMinimumWidth(token_label_min_width)
			grid_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)
			
#			self.qsub_enable_checkbox = QCheckBox("", self)
			self.qsub_enable_checkbox = QCheckBox("")
			if is_qsub_enabled == True:
				self.qsub_enable_checkbox.setCheckState(Qt.Checked)
			else: # assert(is_qsub_enabled == False)
				self.qsub_enable_checkbox.setCheckState(Qt.Unchecked)
#			self.qsub_enable_checkbox.move(self.x2, self.y1)
			self.qsub_enable_checkbox.setToolTip("submit job to queue")
			self.qsub_enable_checkbox.stateChanged.connect(self.set_qsub_enable_state) # To control enable state of the following qsub related widgets
			self.qsub_enable_checkbox.setEnabled(self.sxcmdwidget.sxcmd.is_submittable)
			grid_layout.addWidget(self.qsub_enable_checkbox, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span, token_widget_row_span, token_widget_col_span)
			
			# self.y1 = self.y1 + 25
			grid_row += 1
			
#			temp_label = QLabel("job name", self)
			temp_label = QLabel("job name")
#			temp_label.move(self.x1, self.y1)
#			temp_label.setMinimumWidth(token_label_min_width)
			grid_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)
			
#			self.qsub_job_name_edit = QLineEdit(self)
			self.qsub_job_name_edit = QLineEdit()
			if self.sxcmdwidget.sxcmd.is_submittable == True:
				self.qsub_job_name_edit.setText(self.sxcmdwidget.sxcmd.get_mode_name_for("file_path"))
			else: # assert(self.sxcmdwidget.sxcmd.is_submittable == False)
				assert(self.sxcmdwidget.sxcmd.mpi_support == False)
				self.qsub_job_name_edit.setText("N/A")
#			self.qsub_job_name_edit.move(self.x2, self.y1)
#			self.qsub_job_name_edit.setMinimumWidth(token_widget_min_width)
			self.qsub_job_name_edit.setToolTip("name of this job")
			grid_layout.addWidget(self.qsub_job_name_edit, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span, token_widget_row_span, token_widget_col_span)
			
#			self.y1 = self.y1 + 25
			grid_row += 1
			
#			temp_label = QLabel("submission command", self)
			temp_label = QLabel("submission command")
#			temp_label.move(self.x1, self.y1)
#			temp_label.setMinimumWidth(token_label_min_width)
			grid_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)
			
#			self.qsub_cmd_edit = QLineEdit(self)
			self.qsub_cmd_edit = QLineEdit()
			if self.sxcmdwidget.sxcmd.is_submittable == True:
				self.qsub_cmd_edit.setText("qsub")
			else: # assert(self.sxcmdwidget.sxcmd.is_submittable == False)
				assert(self.sxcmdwidget.sxcmd.mpi_support == False)
				self.qsub_cmd_edit.setText("N/A")
#			self.qsub_cmd_edit.move(self.x2, self.y1)
#			self.qsub_cmd_edit.setMinimumWidth(token_widget_min_width)
			self.qsub_cmd_edit.setToolTip("name of submission command to queue job")
			grid_layout.addWidget(self.qsub_cmd_edit, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span, token_widget_row_span, token_widget_col_span)
			
#			self.y1 = self.y1 + 25
			grid_row += 1
			
#			temp_label = QLabel("submission script template", self)
			temp_label = QLabel("submission script template")
#			temp_label.move(self.x1, self.y1)
#			temp_label.setMinimumWidth(token_label_min_width)
			grid_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)
			
#			self.qsub_script_edit = QLineEdit(self)
			self.qsub_script_edit = QLineEdit()
			if self.sxcmdwidget.sxcmd.is_submittable == True:
				self.qsub_script_edit.setText("msgui_qsub.sh")
			else: # assert(self.sxcmdwidget.sxcmd.is_submittable == False)
				assert(self.sxcmdwidget.sxcmd.mpi_support == False)
				self.qsub_script_edit.setText("N/A")
#			self.qsub_script_edit.move(self.x2, self.y1)
#			self.qsub_script_edit.setMinimumWidth(token_widget_min_width)
			self.qsub_script_edit.setToolTip("file name of submission script template (e.g. $EMAN2DIR/bin/msgui_qsub.sh)")
			grid_layout.addWidget(self.qsub_script_edit, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span, token_widget_row_span, token_widget_col_span)
			
#			self.qsub_script_open_btn = QPushButton("Select Template", self)
			self.qsub_script_open_btn = QPushButton("Select Template")
#			self.qsub_script_open_btn.move(self.x3, self.y1 - 4)
#			self.qsub_script_open_btn.setMinimumWidth(token_widget_min_width)
			self.qsub_script_open_btn.setToolTip("display open file dailog to select job submission script template file")
			self.connect(self.qsub_script_open_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, self.qsub_script_edit))
			grid_layout.addWidget(self.qsub_script_open_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
			
#			self.y1 = self.y1 + 25
			grid_row += 1
			
			# Initialize enable state of qsub related widgets
			self.set_qsub_enable_state()
			
			# Add space
#			self.y1 = self.y1 + 25 * 1
			grid_row += 1
			
			# Add save paramaters button 
#			self.save_params_btn = QPushButton("Save parameters", self)
			self.save_params_btn = QPushButton("Save parameters")
#			self.save_params_btn.move(self.x1-5, self.y1)
			self.save_params_btn.setMinimumWidth(func_btn_min_width)
			self.save_params_btn.setToolTip("save gui parameter settings")
			self.connect(self.save_params_btn, SIGNAL("clicked()"), self.sxcmdwidget.save_params)
			grid_layout.addWidget(self.save_params_btn, grid_row, grid_col_origin, func_btn_row_span, func_btn_col_span)
			
#			self.y1 = self.y1 + 30
			grid_row += 1
			
#			self.cmd_line_btn = QPushButton("Generate command line", self)
			self.cmd_line_btn = QPushButton("Generate command line")
#			self.cmd_line_btn.move(self.x1-5, self.y1)
			self.cmd_line_btn.setMinimumWidth(func_btn_min_width)
			self.cmd_line_btn.setToolTip("generate command line from gui parameter settings and automatically save settings")
			self.connect(self.cmd_line_btn, SIGNAL("clicked()"), self.sxcmdwidget.save_cmd_line)
			grid_layout.addWidget(self.cmd_line_btn, grid_row, grid_col_origin, func_btn_row_span, func_btn_col_span)
			
			# self.y1 = self.y1 + 30
			grid_row += 1
			
			# Add a run button
#			self.execute_btn = QPushButton("Run %s" % self.sxcmdwidget.sxcmd.get_mode_name_for("human"), self)
			self.execute_btn = QPushButton("Run %s" % self.sxcmdwidget.sxcmd.get_mode_name_for("human"))
			# make 3D textured push button look
			custom_style = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px} QPushButton:focus {font: bold; color: #000;border: 2px solid #8D0;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px}"
			self.execute_btn.setStyleSheet(custom_style)
#			self.execute_btn.move(self.x5, self.y1)
			self.execute_btn.setMinimumWidth(func_btn_min_width)
			self.execute_btn.setToolTip("run %s and automatically save gui parameter settings" % self.sxcmdwidget.sxcmd.get_mode_name_for("human"))
			self.connect(self.execute_btn, SIGNAL("clicked()"), self.sxcmdwidget.execute_cmd_line)
			grid_layout.addWidget(self.execute_btn, grid_row, grid_col_origin + func_btn_col_span, func_btn_row_span, func_btn_col_span)
	
	def set_text_entry_widget_enable_state(self, widget, is_enabled):
		# Set enable state and background color of text entry widget according to enable state
		default_palette = QPalette()
		bg_color = default_palette.color(QPalette.Inactive, QPalette.Base) # For text entry widgets
		# bg_color = default_palette.color(QPalette.Inactive, QPalette.Button) # For button
		if is_enabled == False:
			bg_color = default_palette.color(QPalette.Disabled, QPalette.Base) # For text entry widgets
			# bg_color = default_palette.color(QPalette.Disabled, QPalette.Button) # For button
		
		widget.setEnabled(is_enabled)
		palette = widget.palette()
		palette.setColor(widget.backgroundRole(), bg_color)
		widget.setPalette(palette)
	
	def set_qsub_enable_state(self):
		is_enabled = False
		if self.qsub_enable_checkbox.checkState() == Qt.Checked:
			is_enabled = True
		
		# Set enable state and background color of mpi related widgets
		if self.sxcmdwidget.sxcmd.mpi_support:
			self.set_text_entry_widget_enable_state(self.mpi_cmd_line_edit, not is_enabled)
		
		# Set enable state and background color of qsub related widgets
		self.set_text_entry_widget_enable_state(self.qsub_job_name_edit, is_enabled)
		self.set_text_entry_widget_enable_state(self.qsub_cmd_edit, is_enabled)
		self.set_text_entry_widget_enable_state(self.qsub_script_edit, is_enabled)
		self.qsub_script_open_btn.setEnabled(is_enabled)
	
	def handle_restore_widget_event(self, sxcmd_token, widget_index=0):
		if sxcmd_token.type == "function":
			assert(len(sxcmd_token.widget) == 2 and len(sxcmd_token.restore) == 2 and widget_index < 2)
			sxcmd_token.widget[widget_index].setText("%s" % sxcmd_token.restore[widget_index])
		else:
			if sxcmd_token.type == "bool":
				if sxcmd_token.restore == "YES":
					sxcmd_token.widget.setChecked(Qt.Checked)
				else: # sxcmd_token.restore == "NO"
					sxcmd_token.widget.setChecked(Qt.Unchecked)
			else:
				sxcmd_token.widget.setText("%s" % sxcmd_token.restore)

# ========================================================================================
class SXConstSetWindow(QWidget):
	def __init__(self, sxconst_set, sxcmd_list, parent=None):
		super(SXConstSetWindow, self).__init__(parent)
		
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.sxconst_set = sxconst_set
		self.sxcmd_list = sxcmd_list
		
		self.gui_settings_file_path = "%s/gui_settings_project_settings.txt" % (SXLookFeelConst.project_dir)
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# Set the window title and size
		self.setWindowTitle(self.sxconst_set.name)
		self.resize(SXLookFeelConst.sxconst_set_window_width, SXLookFeelConst.sxconst_set_window_height)
		
#		# Set the background color of this widget
#		self.setAutoFillBackground(True)
#		palette = QPalette()
#		palette.setBrush(QPalette.Background, QBrush(SXLookFeelConst.default_bg_color))
#		self.setPalette(palette)
		
		global_row_origin = 0; global_col_origin = 0
		global_row_span = 4; global_col_span = 1
		
		header_row_origin = 0; header_col_origin = 0
		title_row_span = 1; title_col_span = 1
		short_info_row_span = 1; short_info_col_span = 1
		title_min_width = 300
		short_info_min_width = 300
		short_info_min_height = 80
		
		const_set_row_origin = 0; const_set_col_origin = 0
		const_label_row_span = 1; const_label_col_span = 1
		const_register_widget_row_span = 1; const_register_widget_col_span = 1
		const_widget_row_span = 1; const_widget_col_span = 1
		const_label_min_width = 150
		const_register_widget_min_width = const_label_min_width
		const_widget_min_width = const_label_min_width
		
		btn_row_origin = 0; btn_col_origin = 0
		func_btn_row_span = 1; func_btn_col_span = 1
		register_btn_row_span = 1; register_btn_col_span = 2
		func_btn_min_width = 150
		
		box_layout = QVBoxLayout(self)
		box_layout.setContentsMargins(0,0,0,0)
		box_layout.setSpacing(0)
		scroll_area = QScrollArea()
		scroll_area.setWidgetResizable(True)
		scroll_area_widgets = QWidget(scroll_area)
		
		# Set the background color of scroll area
		scroll_area_widgets.setAutoFillBackground(True)
		palette = QPalette()
		palette.setBrush(QPalette.Background, QBrush(SXLookFeelConst.default_bg_color))
		scroll_area_widgets.setPalette(palette)
		
		scroll_area.setWidget(scroll_area_widgets)
		box_layout.addWidget(scroll_area)
		
#		global_layout = QGridLayout()
		global_layout = QGridLayout(scroll_area_widgets)
		global_layout.setMargin(SXLookFeelConst.grid_margin)
		global_layout.setSpacing(SXLookFeelConst.grid_spacing)
		global_layout.setRowStretch(global_row_span - 1, global_layout.rowStretch(global_row_origin) + 1)
		
		header_layout = QGridLayout()
		header_layout.setMargin(SXLookFeelConst.grid_margin)
		header_layout.setSpacing(SXLookFeelConst.grid_spacing)
		
		const_set_layout = QGridLayout()
		const_set_layout.setMargin(SXLookFeelConst.grid_margin)
		const_set_layout.setSpacing(SXLookFeelConst.grid_spacing)
		
		btn_layout = QGridLayout()
		btn_layout.setMargin(SXLookFeelConst.grid_margin)
		btn_layout.setSpacing(SXLookFeelConst.grid_spacing * 2)
		
		# # Define the frame within the widget layout
		# widget_frame = QFrame()
		# # widget_frame.setFrameStyle(QFrame.StyledPanel)
		# global_layout.addWidget(widget_frame, global_row_origin, global_col_origin, global_row_span, global_col_span)
		
		global_grid_row = global_row_origin
		
		# Start add title widgets to the grid layout
		header_grid_row = header_row_origin
		
		# Set a label and its position in this tab
		temp_label = QLabel("<b>%s</b>" % (self.sxconst_set.name))
		temp_label.setMinimumWidth(title_min_width)
		header_layout.addWidget(temp_label, header_grid_row, header_col_origin, title_row_span, title_col_span)
		
		header_grid_row += 1
		
		# NOTE: 2015/11/17 Toshio Moriya
		# Necessary to separate "<b>%s</b>" from the information for avoiding to invoke the tag interpretations of string
		# e.g. < becomes the escape character
		temp_label = QLabel("%s" % (self.sxconst_set.short_info))
		temp_label.setWordWrap(True)
		temp_label.setMinimumWidth(short_info_min_width)
		temp_label.setMinimumHeight(short_info_min_height)
		header_layout.addWidget(temp_label, header_grid_row, header_col_origin, short_info_row_span, short_info_col_span)
		
		# Add const set grid layout to global layout
		global_layout.addLayout(header_layout, global_grid_row, global_col_origin)
		global_grid_row += 1
		
		# Start add project parameter constant widgets to the grid layout
		const_set_grid_row = const_set_row_origin
		
		# Add widget for editing command args and options
		for sxconst in self.sxconst_set.list:
			# Create widget associated to this project constant parameter
			temp_label = QLabel(sxconst.label)
			temp_label.setMinimumWidth(const_label_min_width)
			const_set_layout.addWidget(temp_label, const_set_grid_row, const_set_col_origin, const_label_row_span, const_label_col_span)
			
			sxconst_register_widget = QPushButton("[%s]" % sxconst.register)
			sxconst_register_widget.setMinimumWidth(const_register_widget_min_width)
			custom_style = "QPushButton {color:green; }"
			sxconst_register_widget.setStyleSheet(custom_style)
			const_set_layout.addWidget(sxconst_register_widget, const_set_grid_row, const_set_row_origin + const_label_col_span, const_register_widget_row_span, const_register_widget_col_span)
			sxconst_register_widget.setToolTip("retrieve this registered value to edit box")
			self.connect(sxconst_register_widget, SIGNAL("clicked()"), partial(self.handle_regster_widget_event, sxconst))
			
			sxconst_widget = QLineEdit()
			sxconst_widget.setMinimumWidth(const_widget_min_width)
			sxconst_widget.setText(sxconst.register)
			sxconst_widget.setToolTip(sxconst.help)
			const_set_layout.addWidget(sxconst_widget, const_set_grid_row, const_set_row_origin + const_label_col_span + const_register_widget_col_span, const_widget_row_span, const_widget_col_span)
			
			const_set_grid_row += 1
			
			# Register this widget
			sxconst.register_widget = sxconst_register_widget
			sxconst.widget = sxconst_widget
		
		# Add const set grid layout to global layout
		global_layout.addLayout(const_set_layout, global_grid_row, global_col_origin)
		global_grid_row += 1
		
		# Start add buttons to the grid layout
		btn_grid_row = btn_row_origin
		
		# Add a register button
		self.execute_btn = QPushButton("Register settings")
		# make 3D textured push button look
		custom_style = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px}"
		self.execute_btn.setStyleSheet(custom_style)
		self.execute_btn.setMinimumWidth(func_btn_min_width * register_btn_col_span)
		self.execute_btn.setToolTip("register project constant parameter settings to automatically set values to command arguments and options")
		self.connect(self.execute_btn, SIGNAL("clicked()"), self.register_const_set)
		btn_layout.addWidget(self.execute_btn, btn_grid_row, btn_col_origin, register_btn_row_span, register_btn_col_span)
		
		btn_grid_row += 1
		
		# Add save project constant parameter settings button 
		self.save_consts_btn = QPushButton("Save settings")
		self.save_consts_btn.setMinimumWidth(func_btn_min_width)
		self.save_consts_btn.setToolTip("save project constant parameter settings")
		self.connect(self.save_consts_btn, SIGNAL("clicked()"), self.save_consts)
		btn_layout.addWidget(self.save_consts_btn, btn_grid_row, btn_col_origin, func_btn_row_span, func_btn_col_span)
		
		# Add load project constant parameter settings button 
		self.load_consts_btn = QPushButton("Load settings")
		self.load_consts_btn.setMinimumWidth(func_btn_min_width)
		self.load_consts_btn.setToolTip("load project constant parameter settings to retrieve the previously-saved one")
		self.connect(self.load_consts_btn, SIGNAL("clicked()"), self.load_consts)
		btn_layout.addWidget(self.load_consts_btn, btn_grid_row, btn_col_origin + func_btn_col_span, func_btn_row_span, func_btn_col_span)
		
		btn_grid_row += 1
		
		# Add a close button
		self.close_btn = QPushButton("Close")
		self.close_btn.setMinimumWidth(func_btn_min_width)
		self.close_btn.setToolTip("close this window")
		self.connect(self.close_btn, SIGNAL("clicked()"), self.close)
		btn_layout.addWidget(self.close_btn, btn_grid_row, btn_col_origin + func_btn_col_span, func_btn_row_span, func_btn_col_span)
		
		# Add button grid layout to global layout
		global_layout.addLayout(btn_layout, global_grid_row, global_col_origin)
		
		# Load the previously saved parameter setting of this sx command
		if os.path.exists(self.gui_settings_file_path):
			self.read_consts(self.gui_settings_file_path)
	
	def handle_regster_widget_event(self, sxconst):
		sxconst.widget.setText(sxconst.register)
		
	def register_const_set(self):
		# Loop through all project constant parameters
		for sxconst in self.sxconst_set.list:
			sxconst.register = sxconst.widget.text()
			sxconst.register_widget.setText("[%s]" % sxconst.register)
		
		# Loop through all command tokens
		for sxcmd in self.sxcmd_list:
			for cmd_token in sxcmd.token_list:
				if cmd_token.type in self.sxconst_set.dict.keys():
					sxconst = self.sxconst_set.dict[cmd_token.type]
					cmd_token.restore = sxconst.register
					cmd_token.restore_widget.setText("[%s]" % cmd_token.restore)
					cmd_token.widget.setText(cmd_token.restore)
					# print "MRK_DEBUG: %s, %s, %s, %s, %s" % (sxcmd.name, cmd_token.key_base, cmd_token.type, cmd_token.default, cmd_token.restore)
		
		# Save the current state of GUI settings
		if os.path.exists(SXLookFeelConst.project_dir) == False:
			os.mkdir(SXLookFeelConst.project_dir)
		self.write_consts(self.gui_settings_file_path)
	
	def write_consts(self, file_path_out):
		file_out = open(file_path_out,"w")
		
		# Write script name for consistency check upon loading
		file_out.write("@@@@@ project settings gui settings - ")
		# file_out.write(EMANVERSION + " (CVS" + CVSDATESTAMP[6:-2] +")")
		file_out.write(EMANVERSION + " (GITHUB: " + DATESTAMP +")" )
		file_out.write(" @@@@@ \n")
		
		# Loop through all project constant parameters
		for sxconst in self.sxconst_set.list:
			# The other type has only one line edit box
			val_str = str(sxconst.widget.text())
			file_out.write("<%s> %s (registered %s) == %s \n" % (sxconst.key, sxconst.label, sxconst.register, val_str))
			
		file_out.close()
	
	def read_consts(self, file_path_in):
		file_in = open(file_path_in,"r")
		
		# Check if this parameter file is for this sx script
		line_in = file_in.readline()
		if line_in.find("@@@@@ project settings gui settings") != -1:
			n_function_type_lines = 2
			function_type_line_counter = 0
			# loop through the rest of lines
			for line_in in file_in:
				# Extract label (which should be left of "=="). Also strip the ending spaces
				label_in = line_in.split("==")[0].strip()
				# Extract value (which should be right of "=="). Also strip all spaces
				val_str_in = line_in.split("==")[1].strip() 
				
				# Extract key_base of this command token
				target_operator = "<"
				item_tail = label_in.find(target_operator)
				if item_tail != 0: 
					QMessageBox.warning(self, "Invalid Project Settings File Format", "Project settings entry should start from \"%s\" for entry key in line (%s). The format of this file might be corrupted. Please save the project settings file again." % (target_operator, line_in))
				label_in = label_in[item_tail + len(target_operator):].strip() # Get the rest of line
				target_operator = ">"
				item_tail = label_in.find(target_operator)
				if item_tail == -1: 
					QMessageBox.warning(self, "Invalid Project Settings File Format", "Project settings entry should have \"%s\" closing entry key in line (%s) The format of this file might be corrupted. Please save the project settings file again." % (target_operator, line_in))
				key = label_in[0:item_tail]
				# Get corresponding sxconst
				if key not in self.sxconst_set.dict.keys(): 
					QMessageBox.warning(self, "Invalid Project Settings File Format", "Invalid entry key for project settings \"%s\" is found in line (%s). This project settings file might be imcompatible with the current version. Please save the project settings file again." % (key, line_in))
				sxconst = self.sxconst_set.dict[key]
				sxconst.widget.setText(val_str_in)
						
		else:
			QMessageBox.warning(self, "Fail to load project settings", "The specified file is not project settings file.")
		
		file_in.close()
	
	def save_consts(self):
		file_path_out = str(QFileDialog.getSaveFileName(self, "Save settings", options = QFileDialog.DontUseNativeDialog))
		if file_path_out != "":
			self.write_consts(file_path_out)
	
	def load_consts(self):
		file_path_in = str(QFileDialog.getOpenFileName(self, "Load settings", options = QFileDialog.DontUseNativeDialog))
		if file_path_in != "":
			self.read_consts(file_path_in)

# ========================================================================================
# Layout of the Pop Up window SXInfoWindow; started by the function info of the main window
class SXInfoWindow(QWidget):
	def __init__(self, parent = None):
		super(SXInfoWindow, self).__init__(parent)
		
		#Here we just set the window title and  3 different labels, with their positions in the window
		self.setWindowTitle("SPHIRE-GUI Info Page")
		
		# Set the background color of this widget
		self.setAutoFillBackground(True)
		palette = QPalette()
		palette.setBrush(QPalette.Background, QBrush(SXLookFeelConst.default_bg_color))
		self.setPalette(palette)
		
		label_row_span = 1; label_col_span = 3
		close_row_span = 1; close_col_span = 1
		spacer_min_width = 12
		
		grid_layout = QGridLayout(self)
		grid_layout.setMargin(SXLookFeelConst.grid_margin)
		grid_layout.setSpacing(SXLookFeelConst.grid_spacing)
		
		grid_col = 0
		grid_row = 0; grid_layout.setRowMinimumHeight(grid_row, spacer_min_width)
		grid_row += 1; temp_label=QLabel("<b>SPHIRE GUI Prototype</b>", self)
		# temp_label.move(10,10)
		temp_label.setAlignment(Qt.AlignHCenter)
		grid_layout.addWidget(temp_label, grid_row, grid_col, label_row_span, label_col_span)
		grid_row += 1; temp_label=QLabel("<b>Author: Toshio Moriya</b> ", self)
		# temp_label.move(10,40)
		temp_label.setAlignment(Qt.AlignHCenter)
		grid_layout.addWidget(temp_label, grid_row, grid_col, label_row_span, label_col_span)
		grid_row += 1; grid_layout.setRowMinimumHeight(grid_row, spacer_min_width)
		grid_row += 1; temp_label=QLabel("For more information visit:%s " % SPARX_DOCUMENTATION_WEBSITE, self)
		# temp_label.move(10,70)
		temp_label.setAlignment(Qt.AlignHCenter)
		grid_layout.addWidget(temp_label, grid_row, grid_col, label_row_span, label_col_span)
		grid_row += 1; grid_layout.setRowMinimumHeight(grid_row, spacer_min_width)
		grid_row += 1; grid_layout.setRowMinimumHeight(grid_row, spacer_min_width)
		grid_row += 1; close_btn = QPushButton("Close")
		self.connect(close_btn, SIGNAL("clicked()"), self.close)
		grid_layout.addWidget(close_btn, grid_row, grid_col + 1, close_row_span, close_col_span)

# ========================================================================================
# Base-class for sx command window (inherited by classes SXUtilWindow and SXMainWindow )
# It provides basic functionalities and widgets for sx command related windows 
class SXCmdWindowBase(QWidget):
	def __init__(self, sxcmd_category, sxconst_set, sxcmd_list, parent = None):
		super(SXCmdWindowBase, self).__init__(parent)
		
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.sxcmd_category = sxcmd_category
		self.sxconst_set = sxconst_set
		self.sxcmd_list = sxcmd_list
		self.cur_sxcmd = None
		
		self.expected_cmd_counts = 32
		
		self.menu_widget_row_span = 1
		self.menu_widget_col_span = 2
		
		self.menu_area_row_span = self.menu_widget_row_span * self.expected_cmd_counts
		self.menu_area_col_span = self.menu_widget_col_span
		
		self.sxcmd_widget_row_span = self.menu_area_row_span
		self.sxcmd_widget_col_span = 1
		
		self.sxcmd_area_row_span = self.sxcmd_widget_row_span
		self.sxcmd_area_col_span = self.sxcmd_widget_col_span
		
		self.grid_layout = None # grid layout
		
		self.grid_row_origin = 0
		self.grid_col_origin = 0
		
		self.grid_row = self.grid_row_origin # Keep current row
		self.grid_col = self.grid_col_origin # keep current column
		
		self.cmd_btn_group = None
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		
	def setup_layout(self, scroll_widget_brush):
		# Setup global layout
		global_layout = QVBoxLayout(self)
		global_layout.setContentsMargins(0,0,0,0)
		global_layout.setSpacing(0)
		# Setup scroll area and its widget 
		scroll_area = QScrollArea()
		scroll_area.setWidgetResizable(True)
		scroll_area_widget = QWidget(scroll_area)
		# Setup scroll widget and its background color
		scroll_area_widget.setAutoFillBackground(True)
		palette = QPalette()
		palette.setBrush(QPalette.Background, scroll_widget_brush)
		scroll_area_widget.setPalette(palette)
		# Register the widget to scroll area
		scroll_area.setWidget(scroll_area_widget)
		# Register the scroll area to the global layout
		global_layout.addWidget(scroll_area)
		
		# Setup grid layout in the scroll area
		self.grid_layout = QGridLayout(scroll_area_widget)
		self.grid_layout.setMargin(SXLookFeelConst.grid_margin)
		self.grid_layout.setSpacing(SXLookFeelConst.grid_spacing)
		self.grid_layout.setColumnMinimumWidth(0, SXLookFeelConst.menu_area_min_width // 2)
		self.grid_layout.setColumnMinimumWidth(1, SXLookFeelConst.menu_area_min_width // 2)
		self.grid_layout.setColumnMinimumWidth(2, SXLookFeelConst.sxcmd_area_min_width)
		# Give the column of the command settings area a higher stretch priority so that the other area does not stretch horizontally
		self.grid_layout.setColumnStretch(self.grid_col_origin + self.menu_area_col_span, self.grid_layout.columnStretch(self.grid_col_origin + self.menu_area_col_span) + 1)
		
		# # MRK_DEBUG: Add the menu area frame for checking
		# menu_area_frame = QFrame()
		# menu_area_frame.resize(SXLookFeelConst.menu_area_min_width , SXLookFeelConst.menu_area_min_height)
		# menu_area_frame.setFrameStyle(QFrame.StyledPanel)
		# self.grid_layout.addWidget(menu_area_frame, self.grid_row_origin, self.grid_col_origin, self.menu_area_row_span, self.menu_area_col_span)
		
		# # MRK_DEBUG: Add the sx command area frame for checking
		# cmd_area_frame = QFrame()
		# cmd_area_frame.resize(SXLookFeelConst.sxcmd_area_min_width, SXLookFeelConst.sxcmd_area_min_height)
		# cmd_area_frame.setFrameStyle(QFrame.StyledPanel)
		# self.grid_layout.addWidget(cmd_area_frame, self.grid_row_origin, self.grid_col_origin + self.menu_area_col_span, self.sxcmd_area_row_span, self.sxcmd_area_col_span)
	
	# Add Pipeline SX Commands (sx*.py) associated widgets
	def add_sxcmd_widgets(self):
		self.cmd_btn_group = QButtonGroup()
		# self.cmd_btn_group.setExclusive(True) # NOTE: 2016/02/18 Toshio Moriya: Without QPushButton.setCheckable(True). This does not do anything. Let manually do this
		
		# Add SX Commands (sx*.py) associated widgets
		for sxcmd in self.sxcmd_list:
			if sxcmd.category == self.sxcmd_category:
				# Add buttons for this sx*.py processe
				sxcmd.btn = QPushButton(sxcmd.label)
				# sxcmd.btn.setCheckable(True) # NOTE: 2016/02/18 Toshio Moriya: With this setting, we can not move the focus to the unchecked butttons... PyQt bug?
				sxcmd.btn.setToolTip(sxcmd.short_info)
				self.cmd_btn_group.addButton(sxcmd.btn)
				self.grid_layout.addWidget(sxcmd.btn, self.grid_row, self.grid_col_origin, self.menu_widget_row_span, self.menu_widget_col_span)
				
				# Create SXCmdWidget for this sx*.py processe
				sxcmd.widget = SXCmdWidget(self.sxconst_set, sxcmd)
				sxcmd.widget.hide()
				self.grid_layout.addWidget(sxcmd.widget, self.grid_row_origin, self.grid_col_origin + self.menu_area_col_span, self.sxcmd_widget_row_span, self.sxcmd_widget_col_span)
				
				# connect widget signals
				self.connect(sxcmd.btn, SIGNAL("clicked()"), partial(self.handle_sxcmd_btn_event, sxcmd))
				
				self.grid_row += 1
	
	def handle_sxcmd_btn_event(self, sxcmd):
		modifiers = QApplication.keyboardModifiers()
		if modifiers == Qt.ShiftModifier:
			os.system("python -m webbrowser %s%s" % (SPARX_DOCUMENTATION_WEBSITE, sxcmd.name))
			return
		
		if self.cur_sxcmd == sxcmd: return
		
		if self.cur_sxcmd != None:
			assert(self.cur_sxcmd.widget.isVisible() == True)
			self.cur_sxcmd.widget.hide()
			custom_style = "QPushButton {color:black; }" # custom_style = "QPushButton {color:#000; }"
			self.cur_sxcmd.btn.setStyleSheet(custom_style)
			
		self.cur_sxcmd = sxcmd
		
		if self.cur_sxcmd != None:
			assert(self.cur_sxcmd.widget.isVisible() == False)
			self.cur_sxcmd.widget.show()
			custom_style = "QPushButton {font: bold; color:blue; }" # custom_style = "QPushButton {font: bold; color:#8D0; }"
			self.cur_sxcmd.btn.setStyleSheet(custom_style)
			
	def handle_open_window_btn_event(self, window, position_shift):
		if window.isVisible():
			window.raise_()
			window.activateWindow()
			return
		
		window.move(self.pos() + position_shift)
		window.show()
		window.raise_()
		
# ========================================================================================
# Utility Window (opened by class SXMainWindow)
# This class includes the layout of the utility window
class SXUtilWindow(SXCmdWindowBase):
	def __init__(self, sxconst_set, sxcmd_list, parent = None):
		super(SXUtilWindow, self).__init__("util", sxconst_set, sxcmd_list, parent)
		
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		
		# --------------------------------------------------------------------------------
		# Setup Window Properties
		# --------------------------------------------------------------------------------
		self.setWindowTitle("Utilities")
		self.resize(SXLookFeelConst.sxutil_window_width, SXLookFeelConst.sxutil_window_height)
		
		# --------------------------------------------------------------------------------
		# Setup Window Layout
		# --------------------------------------------------------------------------------
		self.setup_layout(QBrush(SXLookFeelConst.default_bg_color))
		
		# --------------------------------------------------------------------------------
		# Add General Buttons
		# --------------------------------------------------------------------------------
		# Add Close button
		close_btn = QPushButton("Close")
		close_btn.setToolTip("close Utility Window")
		self.grid_layout.addWidget(close_btn, self.grid_row, self.grid_col_origin, self.menu_widget_row_span, self.menu_widget_col_span)
		self.connect(close_btn, SIGNAL("clicked()"), self.close)
		
		self.grid_row += 1
		
		# --------------------------------------------------------------------------------
		# Add Utility SX Commands (sx*.py) associated widgets
		# --------------------------------------------------------------------------------
		# Add title label and set position and font style
		title=QLabel("<span style=\'font-size:18pt; font-weight:600; color:#aa0000;\'><b>UTILITIES </b></span><span style=\'font-size:12pt; font-weight:60; color:#aa0000;\'>(shift-click for wiki)</span>")
		self.grid_layout.addWidget(title, self.grid_row, self.grid_col_origin, self.menu_widget_row_span, self.menu_widget_col_span)
		
		self.grid_row += 1
		
		self.add_sxcmd_widgets()
		
		# --------------------------------------------------------------------------------
		# Register project constant parameter settings upon initialization
		# --------------------------------------------------------------------------------
		self.sxconst_set.window.register_const_set()
		
		# --------------------------------------------------------------------------------
		# Load the previously saved parameter setting of this sx command
		# Override the registration of project constant parameter settings with the previously-saved one
		# --------------------------------------------------------------------------------
		for sxcmd in self.sxcmd_list:
			if os.path.exists(sxcmd.widget.gui_settings_file_path):
				sxcmd.widget.read_params(sxcmd.widget.gui_settings_file_path)
	
#	def closeEvent(self, event):
#		# Quit child applications of all sxcmd widgets
#		for sxcmd in self.sxcmd_list:
#			sxcmd.widget.quit_all_child_applications()
#		
#		event.accept()

# ========================================================================================
# Child Window status, used by only main window
class SXChildStatus:
	def __init__(self, window = None):
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.window = window       # Child window widget
		self.is_minimized = False  # Is this window currently minimized?
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	
# ========================================================================================
# Main Window (started by class SXApplication)
class SXMainWindow(SXCmdWindowBase):
	
	def __init__(self, sxconst_set, sxcmd_list, parent = None):
		super(SXMainWindow, self).__init__("pipe", sxconst_set, sxcmd_list, parent)
		
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.child_status_list = [] 
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		
		# --------------------------------------------------------------------------------
		# Setup Window Properties
		# --------------------------------------------------------------------------------
		self.setWindowTitle("SPHIRE-GUI Main (Alpha Version)")
		self.resize(SXLookFeelConst.sxmain_window_width, SXLookFeelConst.sxmain_window_height)
		self.installEventFilter(self) # Necessary for self.eventFilter()
		# print "MRK_DEBUG: SXMainWindow Window Flags: 0x%08x " % (self.windowFlags())
		
		# --------------------------------------------------------------------------------
		# Setup Window Layout
		# --------------------------------------------------------------------------------
		self.setup_layout(QBrush(QPixmap(get_image_directory() + "sxgui.py_main_window_background_image.png")))
		
		# --------------------------------------------------------------------------------
		# General 
		# --------------------------------------------------------------------------------
		# Create SXInfoWindow
		sxinfo_window = SXInfoWindow()
		sxinfo_window.setWindowFlags((sxinfo_window.windowFlags()| Qt.CustomizeWindowHint) & ~Qt.WindowMinimizeButtonHint) # Disabled minimize icon button in window title bar
		sxinfo_window.setWindowFlags((sxinfo_window.windowFlags()| Qt.CustomizeWindowHint) & ~Qt.WindowMaximizeButtonHint) # Disabled maximize icon button in window title bar
		sxinfo_window.hide()
		self.child_status_list.append(SXChildStatus(sxinfo_window))
		
		# Add Push button to display popup window for info about the application
		info_btn = QPushButton()
		icon = QIcon(get_image_directory() + "sparxicon.png") # Decorates the button with the sphire image
		info_btn.setIcon(icon)
		info_btn.setToolTip("open info page")
		self.grid_layout.addWidget(info_btn, self.grid_row, self.grid_col_origin, self.menu_widget_row_span, self.menu_widget_col_span // 2)
		position_shift = QPoint(0, 0)
		self.connect(info_btn, SIGNAL("clicked()"), partial(self.handle_open_window_btn_event, sxinfo_window, position_shift))
		
		# Add Close button
		quit_btn = QPushButton("Quit")
		quit_btn.setToolTip("quit SPHIRE GUI ")
		self.grid_layout.addWidget(quit_btn, self.grid_row, self.grid_col_origin + self.menu_widget_col_span // 2, self.menu_widget_row_span, self.menu_widget_col_span // 2)
		# self.connect(quit_btn, SIGNAL("clicked()"), qApp, SLOT("quit()")) # NOTE: 2016/03/01 Toshio Moriya; This does not emit close event
		self.connect(quit_btn, SIGNAL("clicked()"), self.close)
		
		self.grid_row += 1
		
		# Add title label and set position and font style
		title=QLabel("<span style=\'font-size:18pt; font-weight:600; color:#aa0000;\'><b>PROJECT SETTINGS </b></span>")
		self.grid_layout.addWidget(title, self.grid_row, self.grid_col_origin, self.menu_widget_row_span, self.menu_widget_col_span)
		
		self.grid_row += 1
		
		# Create SXConstSetWindow
		sxconst_set.window = SXConstSetWindow(self.sxconst_set, self.sxcmd_list)
		sxconst_set.window.setWindowFlags((sxconst_set.window.windowFlags()| Qt.CustomizeWindowHint) & ~Qt.WindowMinimizeButtonHint) # Disabled minimize icon button in window title bar
		sxconst_set.window.hide()
		self.child_status_list.append(SXChildStatus(sxconst_set.window))
		
		# Add project constant parameter set associated button
		self.sxconst_set.btn = QPushButton("Open %s Window" % self.sxconst_set.label)
		self.sxconst_set.btn.setToolTip(self.sxconst_set.short_info)
		self.grid_layout.addWidget(self.sxconst_set.btn, self.grid_row, self.grid_col_origin, self.menu_widget_row_span, self.menu_widget_col_span)
		position_shift = QPoint(-1 * SXLookFeelConst.sxconst_set_window_width, 0)
		self.connect(self.sxconst_set.btn, SIGNAL("clicked()"), partial(self.handle_open_window_btn_event, self.sxconst_set.window, position_shift))
		
		self.grid_row += 1
		
		# --------------------------------------------------------------------------------
		# Add Pipeline SX Commands (sx*.py) associated widgets
		# --------------------------------------------------------------------------------
		# Add title label and set position and font style
		title=QLabel("<span style=\'font-size:18pt; font-weight:600; color:#aa0000;\'><b>PROGRAMS </b></span><span style=\'font-size:12pt; font-weight:60; color:#aa0000;\'>(shift-click for wiki)</span>")
		self.grid_layout.addWidget(title, self.grid_row, self.grid_col_origin, self.menu_widget_row_span, self.menu_widget_col_span)
		self.grid_row += 1
		
		self.add_sxcmd_widgets()
		
		# --------------------------------------------------------------------------------
		# Create Utility SX Commands (sx*.py) related widgets
		# --------------------------------------------------------------------------------
		# Add title label and set position and font style
		title=QLabel("<span style=\'font-size:18pt; font-weight:600; color:#aa0000;\'><b>UTILITIES </b></span>")
		self.grid_layout.addWidget(title, self.grid_row, self.grid_col_origin, self.menu_widget_row_span, self.menu_widget_col_span)
		self.grid_row += 1
		
		# Create Utility Window
		sxutil_window = SXUtilWindow(self.sxconst_set, self.sxcmd_list)
		sxutil_window.setWindowFlags((sxutil_window.windowFlags()| Qt.CustomizeWindowHint) & ~Qt.WindowMinimizeButtonHint) # Disabled minimize icon button in window title bar
		sxutil_window.hide()
		self.child_status_list.append(SXChildStatus(sxutil_window))
		
		# Add all utilities button
		open_util_window_btn = QPushButton("Open Utilities Window")
		open_util_window_btn.setToolTip("open utilities window")
		position_shift = QPoint(SXLookFeelConst.menu_area_min_width + SXLookFeelConst.grid_margin * 2, 0)
		self.connect(open_util_window_btn, SIGNAL("clicked()"), partial(self.handle_open_window_btn_event, sxutil_window, position_shift))
		
		self.grid_layout.addWidget(open_util_window_btn, self.grid_row, self.grid_col_origin, self.menu_widget_row_span, self.menu_widget_col_span)
		self.grid_row += 1
		
		# --------------------------------------------------------------------------------
		# Register project constant parameter settings upon initialization
		# --------------------------------------------------------------------------------
		self.sxconst_set.window.register_const_set()
		
		# --------------------------------------------------------------------------------
		# Load the previously saved parameter setting of this sx command
		# Override the registration of project constant parameter settings with the previously-saved one
		# --------------------------------------------------------------------------------
		for sxcmd in self.sxcmd_list:
			if os.path.exists(sxcmd.widget.gui_settings_file_path):
				sxcmd.widget.read_params(sxcmd.widget.gui_settings_file_path)
	
	def eventFilter(self, source, event):
		if event.type() == QEvent.WindowStateChange:
			# print "MRK_DEBUG: Hello QEvent.WindowStateChange"
			if self.windowState() & Qt.WindowMinimized:
				# print "MRK_DEBUG: sxgui main window has minimized"
				assert(self.isMinimized() == True)
				#
				# NOTE: 2016/03/07 Toshio Moriya
				# Minimize icon button of child window should be disabled
				#
				for child_status in self.child_status_list:
					if child_status.window.isVisible(): 
						# assert(not child_status.window.isMinimized())
						child_status.window.hide()
						child_status.is_minimized = True
						# child_status.window.setWindowState(child_status.window.windowState() | Qt.WindowMinimized)
#				
#				if self.sxinfo_window.isVisible(): 
#					# assert(not self.sxinfo_window.isMinimized())
#					self.sxinfo_window.hide()
#					self.sxinfo_window_minimized = True
#					# self.sxinfo_window.setWindowState(self.sxinfo_window.windowState() | Qt.WindowMinimized)
#				if  self.sxconst_set.window.isVisible(): 
#					# assert(not self.sxconst_set.window.isMinimized())
#					self.sxconst_set.window.hide()
#					self.sxconst_set_window_minimized = True
#					# self.sxconst_set.window.setWindowState(self.sxconst_set.window.windowState() | Qt.WindowMinimized)
#				if self.sxutil_window.isVisible(): 
#					# assert(not self.sxutil_window.isMinimized())
#					self.sxutil_window.hide()
#					self.sxutil_window_minimized = True
#					# self.sxutil_window.setWindowState(self.sxutil_window.windowState() | Qt.WindowMinimized)
#				
			else:
				# print "MRK_DEBUG: sxgui main window has not minimized"
				assert(self.isMinimized() == False)
				#
				# NOTE: 2016/03/07 Toshio Moriya
				# Minimize icon button of child window should be disabled
				#
				for child_status in self.child_status_list:
					# if self.sxinfo_window.isMinimized():
					if child_status.is_minimized == True:
						assert(not child_status.window.isVisible()) 
						# self.child_status.window.setWindowState(child_status.window.windowState() & ~Qt.WindowMinimized)
						child_status.window.show()
						child_status.window.raise_()
						child_status.is_minimized = False
				assert(self.isVisible()) 
				self.raise_()
#				
#				# if self.sxinfo_window.isMinimized():
#				if self.sxinfo_window_minimized == True:
#					assert(not self.sxinfo_window.isVisible()) 
#					# self.sxinfo_window.setWindowState(self.sxinfo_window.windowState() & ~Qt.WindowMinimized)
#					self.sxinfo_window.show()
#					self.sxinfo_window.raise_()
#					self.sxinfo_window_minimized = False
#				# if self.sxconst_set.window.isMinimized(): 
#				if self.sxconst_set_window_minimized == True: 
#					assert(not self.sxconst_set.window.isVisible()) 
#					# self.sxconst_set.window.setWindowState(self.sxconst_set.window.windowState() & ~Qt.WindowMinimized)
#					self.sxconst_set.window.show()
#					self.sxconst_set.window.raise_()
#					self.sxconst_set_window_minimized = False
#				# if self.sxutil_window.isMinimized(): 
#				if self.sxutil_window_minimized == True: 
#					assert(not self.sxutil_window.isVisible()) 
#					# self.sxutil_window.setWindowState(self.sxutil_window.windowState() & ~Qt.WindowMinimized)
#					self.sxutil_window.show()
#					self.sxutil_window.raise_()
#					self.sxutil_window_minimized = False
#				assert(self.isVisible()) 
#				self.raise_()
#				
		elif event.type() == QEvent.WindowActivate:
			# print "MRK_DEBUG: sxgui main window has gained focus (beome active)"
			for child_status in self.child_status_list:
				if child_status.window.isVisible(): 
					child_status.window.raise_()
			assert(self.isVisible()) 
			self.raise_()
			
#			if self.sxinfo_window.isVisible(): 
#				self.sxinfo_window.raise_()
#			if self.sxconst_set.window.isVisible(): 
#				self.sxconst_set.window.raise_()
#			if self.sxutil_window.isVisible(): 
#				self.sxutil_window.raise_()
#			self.raise_()
		# elif event.type()== QEvent.WindowDeactivate:
		# 	print "MRK_DEBUG: sxgui main window has lost focus (beome deactive)"
		# elif event.type()== QEvent.FocusIn:
		# 	print "MRK_DEBUG: sxgui main has gained keyboard focus"
		# elif event.type()== QEvent.FocusOut:
		# 	print "MRK_DEBUG: sxgui main has lost keyboard focus"
		
		return super(SXMainWindow, self).eventFilter(source, event)
	
	def closeEvent(self, event):
		# Quit child applications of all sxcmd widgets
		for sxcmd in self.sxcmd_list:
			sxcmd.widget.quit_all_child_applications()
		
		# close all child windows
		for child_status in self.child_status_list:
			if child_status.window:
				child_status.window.close()
		
#		if self.sxinfo_window:
#			self.sxinfo_window.close()
#		if self.sxutil_window:
#			self.sxutil_window.close()
#		if self.sxconst_set.window:
#			self.sxconst_set.window.close()
		
		event.accept()

# ========================================================================================
class SXApplication(QApplication):
	def __init__(self, *args):
#		QApplication.__init__(self, *args)
		super(SXApplication, self).__init__(*args)
		
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.sxmain_window = None
		
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		"""
		style=QtGui.QStyleFactory.create("Cleanlooks")
		if style==None:
			print "Note: standard Cleanlooks style not available, controls may be distorted. Using ",
			# the first one should work, but we have the loop, just in case
			for s in list(QtGui.QStyleFactory.keys()):
				style=QtGui.QStyleFactory.create(s)
				if style!=None: 
					print s
					break
		if style!=None: self.setStyle(style)
		"""
		
		# Construct list of sxscript objects (extracted from associated wiki documents)
		sxconst_set = construct_sxconst_set()
		
		# Construct list of sxscript objects (extracted from associated wiki documents)
		sxcmd_list = construct_sxcmd_list()
		
		# Define the main window (class SXMainWindow)
		self.sxmain_window = SXMainWindow(sxconst_set, sxcmd_list)
#		# Best setting for MAC OSX
#		# self.sxmain_window.resize(SXLookFeelConst.sxcmd_area_min_width + SXLookFeelConst.menu_area_min_width + SXLookFeelConst.grid_margin * (7 + 1), SXLookFeelConst.sxcmd_area_min_height + SXLookFeelConst.grid_margin * 2)
#		# Best setting for Linux
#		self.sxmain_window.resize(SXLookFeelConst.sxcmd_area_min_width + SXLookFeelConst.menu_area_min_width + SXLookFeelConst.grid_margin * (7 + 7), SXLookFeelConst.sxcmd_area_min_height + SXLookFeelConst.grid_margin * 2)
		
		# main_size = self.sxmain_window.minimumSizeHint();
		# desktop = self.desktop();
		# screen_rect = desktop.screenGeometry();
		# self.sxmain_window.move((screen_rect.width()/2) - (main_size.width()/2), (screen_rect.height()/2) - (main_size.height()/2));
		self.sxmain_window.move(QPoint(SXLookFeelConst.sxconst_set_window_width, 0));
		
		# When all windows are closed, function byebye of class SXApplication will be called
		self.connect(self, SIGNAL("lastWindowClosed()"), self.byebye )
		
		# Show main window
		self.sxmain_window.show()
		self.sxmain_window.raise_()
		
	#function byebye (just quit)  
	def byebye( self ):
		print" bye bye!"
		self.exit(0)

# ========================================================================================
def main(args):
#	from optparse import OptionParser
#	progname = os.path.basename(sys.argv[0])
#	usage = """%prog [options] <image>
#	
#Automatic and manual particle selection. This version is specifically aimed at square boxes
#for single particle analysis."""
#	parser = OptionParser(usage=usage,version=EMANVERSION)
#	
#	parser.add_option("--demo", type="string", default="",   help="Name of the demo whose input parameters should be used to initialize the GUI fields: --demo=mpibdb means the input parameters in demo/mpi_bdb will be used, --demo=mpibdbctf means the input parameters in demo/mpi_bdb_ctf will be used")
#	global options
#	(options, args) = parser.parse_args()
#	global DEMO_mpibdbctf
#	DEMO_mpibdbctf = "mpibdbctf"
#	global DEMO_mpibdb
#	DEMO_mpibdb = "mpibdb"
	
#	global sxapplication
	sxapplication = SXApplication(args)
	sxapplication.setWindowIcon(QIcon(get_image_directory()+"sparxicon.png"))
	
	app_font = sxapplication.font()
	app_font_info = QFontInfo(sxapplication.font())
	new_point_size = app_font_info.pointSize() + 1
	# # MRK_DEBUG: Check the default system font
	# print "MRK_DEBUG: app_font_info.style()      = ", app_font_info.style()
	# print "MRK_DEBUG: app_font_info.styleHint()  = ", app_font_info.styleHint()
	# print "MRK_DEBUG: app_font_info.styleName()  = ", app_font_info.styleName()
	# print "MRK_DEBUG: app_font_info.family()     = ", app_font_info.family()
	# print "MRK_DEBUG: app_font_info.fixedPitch() = ", app_font_info.fixedPitch()
	# print "MRK_DEBUG: app_font_info.pixelSize()  = ", app_font_info.pixelSize()
	# print "MRK_DEBUG: app_font_info.pointSize()  = ", app_font_info.pointSize()
	# print "MRK_DEBUG: app_font_info.pointSizeF() = ", app_font_info.pointSizeF()
	# print "MRK_DEBUG: app_font_info.bold ()      = ", app_font_info.bold()
	# print "MRK_DEBUG: app_font_info.italic()     = ", app_font_info.italic()
	# 
	# NOTE: 2019/02/19 Toshio Moriya
	# The following method of changing font size works with Linux.
	# However, it does not work Mac OSX. The text of widget classes below won't change,
	# still showing the default font size:
	# QPushButton, QLable, Window Title, and QToolTip
	# 
	app_font.setPointSize(new_point_size) # app_font.setPointSize(13) # and setPointSizeF() are device independent, while setPixelSize() is device dependent
	sxapplication.setFont(app_font)
	
	# sxapplication.setStyleSheet("QPushButton {font-size:18pt;}");  # NOTE: 2016/02/19 Toshio Moriya: Doesn't work 
	# sxapplication.setStyleSheet("QLabel {font-size:18pt;}"); # NOTE: 2016/02/19 Toshio Moriya: Doesn't work 
	# sxapplication.setStyleSheet("QToolTip {font-size:14pt; color:white; padding:2px; border-width:2px; border-style:solid; border-radius:20px; background-color: black; border: 1px solid white;}");
	sxapplication.setStyleSheet("QToolTip {font-size:%dpt;}" % (new_point_size));
#	sxapplication.setStyleSheet("QScrollArea {background-color: transparent;}");
#	sxapplication.setStyleSheet("QScrollArea > QWidget > QWidget {background-color: transparent;}");
	
	sxapplication.sxmain_window.show()
	sxapplication.sxmain_window.raise_()
	sxapplication.exec_()

# ========================================================================================
if __name__ == "__main__":
	main(sys.argv)

# ========================================================================================
# END OF SCRIPT
# ========================================================================================

