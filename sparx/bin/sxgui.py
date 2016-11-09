#!/usr/bin/env python
#
# Authors:
# Toshio Moriya, 11/11/2015 (toshio.moriya@mpi-dortmund.mpg.de)
# Markus Stabrin, 09/06/2016 (markus.stabrin@mpi-dortmund.mpg.de)
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
# Helper Functions
# 
# This function is added here because db_convert_path in EMAN2db.py has a bug.
# 
def translate_to_bdb_path(std_path):
	'''
	Translate a standard file path (std_path) to bdb syntax (return value). 
	The path pass must contain at lease EMAN2DB directory and .bdb file name.
	For instance, if the path is particles/EMAN2DB/data.bdb,
	will return bdb:particles#data.
	'''
	
	# Check error conditions
	if not isinstance(std_path,str): 
		raise RuntimeError("Path has to be a string")
	path_tokens = std_path.split("/")
	
	if len(path_tokens) < 2: 
		raise ValueError("Invalid file path. The path pass must contain at least \'EMAN2DB\' directory and \'.bdb\' file name (e.g \'./EMAN2DB/data.bdb\'). ")

	if path_tokens[-2] != "EMAN2DB": 
		raise ValueError("Invalid file path. The path pass must contain \'EMAN2DB\' directory (e.g \'./EMAN2DB/data.bdb\').")
	
	if os.path.splitext(path_tokens[-1])[1] != ".bdb": 
		raise ValueError("Path is invalid. The path pass must contain \'.bdb\' file name (e.g \'./EMAN2DB/data.bdb\').")
	
	# If necessary, compose directory path as a relative path at first
	dir = ""
	if len(path_tokens) > 2:
		for idx in xrange(0, len(path_tokens) - 2):
			if idx != 0:
				dir += "/"
			dir += path_tokens[idx] # accrue the directory
	
	# if the input file path is a absolute path, add '/' at the head of the path
	if std_path[0] == "/" and dir[0] != "/": 
		dir = "/" + dir
	
	# Add '#' before the database name (file basename without extension)
	bdb_path = "bdb:"
	if dir != "":
		bdb_path += dir + "#"
	# Finally, add file basename (without .bdb extension)
	assert(os.path.splitext(path_tokens[-1])[1] == ".bdb")
	bdb_path += os.path.splitext(path_tokens[-1])[0]
	
	return bdb_path

# ========================================================================================
# Inherited by SXcmd_category and SXconst_set
# SXMainWindow use this class to handle events from menu item buttons
class SXmenu_item(object):
	def __init__(self, name = "", label = "", short_info = ""):
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.name = name              # Name of this menu item, used as a key of dictionary
		self.label = label            # User friendly name of this menu item
		self.short_info = short_info  # Short description of this menu item
		self.btn = None               # <Used only in sxgui.py> QPushButton button instance associating with this menu item
		self.widget = None            # <Used only in sxgui.py> SXCmdWidget instance associating with this menu item
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

# ========================================================================================
class SXcmd_token(object):
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
class SXcmd(object):
	def __init__(self, category = "", role = "", is_submittable = True):
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.name = ""                        # Name of this command (i.e. name of sx*.py script but without .py extension), used for generating command line
		self.mode = ""                        # key base name of a command token, defining mode/subset of this command. For fullset command, use empty string
		self.label = ""                       # User friendly name of this command
		self.short_info = ""                  # Short description of this command
		self.mpi_support = False              # Flag to indicate if this command suppors MPI version
		self.mpi_add_flag = False             # DESIGN_NOTE: 2015/11/12 Toshio Moriya. This can be removed when --MPI flag is removed from all sx*.py scripts
		self.category = category              # Category of this command: sxc_movie_micrograph, sxc_ctf, sxc_particle_stack, sxc_2d_clustering, sxc_initial_3d_modeling, sxc_3d_refinement, sxc_3d_clustering, sxc_utilities
		self.role = role                      # Role of this command; sxr_pipe (pipeline), sxr_alt (alternative) sxr_util (utility)
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

	def get_category_dir_path(self, parent_dir_path = ""):
		category_dir_path = self.category.replace("sxc_", "")
		if parent_dir_path != "":
			category_dir_path = os.path.join(parent_dir_path, category_dir_path)

		return category_dir_path

# ========================================================================================
class SXcmd_category(SXmenu_item):
	def __init__(self, name = "", label = "", short_info = ""):
		super(SXcmd_category, self).__init__(name, label, short_info)

		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		# self.name = name              # <Inherit from SXmenu_item> Name of this command category (i.e. sxc_movie_micrograph, sxc_ctf, sxc_particle_stack, sxc_2d_clustering, sxc_initial_3d_modeling, sxc_3d_refinement, sxc_3d_clustering, sxc_utilities), used as a key of dictionary
		# self.label = label            # <Inherit from SXmenu_item> User friendly name of this command category
		# self.short_info = short_info  # <Inherit from SXmenu_item> Short description of this command category
		self.cmd_list = []              # <Used only in sxgui.py> list of commands in this category. Need this to keep the order of commands
#		self.cmd_dict = {}              # <Used only in sxgui.py> dictionary of commands in this category, organised by names of commands. Easy to access a command but looses their order
		# self.btn = None               # <Inherit from SXmenu_item> QPushButton button instance associating with this category
		# self.widget = None            # <Inherit from SXmenu_item> SXCmdWidget instance associating with this category

		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

# ========================================================================================
class SXconst(object):
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
class SXconst_set(SXmenu_item):
	def __init__(self):
		super(SXconst_set, self).__init__()
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		# self.name = ""        # <Inherit from SXmenu_item> Name of this constant parameter set
		# self.label = ""       # <Inherit from SXmenu_item> User friendly name of this set
		# self.short_info = ""  # <Inherit from SXmenu_item> Short description of this set
		self.list = []          # <Used only in sxgui.py> list of constant parameters. Need this to keep the order of constant parameters
		self.dict = {}          # <Used only in sxgui.py> dictionary of constant parameters, organised by keys of constant parameters. Easy to access each constant parameter but looses their order
		# self.btn = None       # <Inherit from SXmenu_item> QPushButton button instance associating with this set
		# self.widget = None    # <Inherit from SXmenu_item> Widget instance associating with this set
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

# ========================================================================================
class SXLookFeelConst(object):
	# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	# static class variables
	# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	default_bg_color = QColor(229, 229, 229, 192) # default_bg_color = QColor(229, 229, 229, 242) # Greyish-White Transparent
	sxinfo_widget_bg_color = QColor(0, 0, 0, 10) # Almost-Completely Transparent
	sxcmd_widget_bg_color = QColor(0, 0, 0, 0) # Completely Transparent
	sxcmd_tab_bg_color = QColor(229, 229, 229, 200) # White Transparent

	# Constants
	project_dir = "sxgui_settings"
	sxmain_window_left = 0
	sxmain_window_top = 0
	sxmain_window_min_width = 1500 # Requirement of specification
	sxmain_window_min_height = 360 # Requirement of specification
	expected_cmd_counts = 32
	grid_margin = 6 # grid_margin = 12
	grid_spacing = 6

	# Constants initialised with invalid values.
	# Valid values should be set by initialise() function
	screen_height = -1
	screen_width = -1
	sxmain_window_width = -1
	sxmain_window_height = -1
	sxmenu_item_btn_width = -1
	grid_distance = -1
	sxmenu_btn_area_min_width = -1
	sxcmd_btn_area_min_width = -1
	sxcmd_widget_area_min_width = -1

	file_dialog_dir = ""

	@staticmethod
	def initialise(sxapp):
		# Set the directory for all file dialogs to script directory
		SXLookFeelConst.file_dialog_dir = os.getcwd()

		monitor_index = 0
		# Search for maximun screen height and set it to SXLookFeelConst singleton class
		max_screen_height = sxapp.desktop().screenGeometry().height()
		for index in range(sxapp.desktop().screenCount()):
			screen_height = sxapp.desktop().screenGeometry(index).height()
			if max_screen_height < screen_height:
				monitor_index = index
				max_screen_height = screen_height
		SXLookFeelConst.screen_height = max_screen_height
		# Search for maximun screen width and set it to SXLookFeelConst singleton class
		SXLookFeelConst.screen_width = sxapp.desktop().screenGeometry(monitor_index).width()

		# Set size of the main window depending on the screen size
		SXLookFeelConst.sxmain_window_height = SXLookFeelConst.screen_height / 2
		if SXLookFeelConst.sxmain_window_height <= SXLookFeelConst.sxmain_window_min_height:
			SXLookFeelConst.sxmain_window_height = SXLookFeelConst.sxmain_window_min_height

		SXLookFeelConst.sxmain_window_width = SXLookFeelConst.sxmain_window_min_width
		if SXLookFeelConst.sxmain_window_width >= SXLookFeelConst.screen_width * 3 / 4:
			SXLookFeelConst.sxmain_window_width = SXLookFeelConst.screen_width * 3 / 4
			if SXLookFeelConst.sxmain_window_width < 960:
				SXLookFeelConst.sxmain_window_width = 960

		# SXLookFeelConst.sxmain_window_height = SXLookFeelConst.screen_height / 2
		# SXLookFeelConst.sxmain_window_width =SXLookFeelConst.sxmain_window_min_width

		SXLookFeelConst.sxmenu_item_btn_width = SXLookFeelConst.sxmain_window_height * 0.125
		SXLookFeelConst.grid_distance = SXLookFeelConst.sxmenu_item_btn_width / 10

		SXLookFeelConst.sxmenu_btn_area_min_width = 2 * SXLookFeelConst.sxmenu_item_btn_width + SXLookFeelConst.grid_distance + 18
		SXLookFeelConst.sxcmd_btn_area_min_width = 240
		SXLookFeelConst.sxcmd_widget_area_min_width = SXLookFeelConst.sxmain_window_width - SXLookFeelConst.sxmenu_btn_area_min_width - SXLookFeelConst.sxcmd_btn_area_min_width

	@staticmethod
	def format_path(path):
		formatted_path = os.path.relpath(path)
		if formatted_path[:len("../")] == "../":
			# if the path is above the project root directory (current directory)
			# use absolute path
			formatted_path = path
		# else:
			# if the path is project subdirectory
			# use relative path

		return formatted_path

# ========================================================================================
class SXLogoButton(QPushButton):
	def __init__(self, logo_file_path, parent = None):
		super(SXLogoButton, self).__init__(parent)

		# print "MRK_DEBUG: logo_file_path = %s" % logo_file_path
		# print "MRK_DEBUG: os.path.exists(logo_file_path) %s" % os.path.exists(logo_file_path)

		# Width of logo image
		logo_width = SXLookFeelConst.sxmenu_item_btn_width * 2 + SXLookFeelConst.grid_distance

		# Style of widget
		self.setFixedSize(logo_width, 0.434 * logo_width)
		self.customButtonStyle = """
			SXLogoButton {{background-color: rgba(0, 0, 0, 0); border: 0px solid black; border-radius: 0px; image: url("{0}");}}
			SXLogoButton:focus {{background-color: rgba(0, 0, 0, 0); border: 0px solid grey; border-radius: 0px; image: url("{0}");}}
			SXLogoButton:pressed {{background-color: rgba(0, 0, 0, 0); border: 0px solid red; border-radius: 0px; image: url("{0}");}}
			""".format(logo_file_path)
		self.customButtonStyleClicked = """
			SXLogoButton {{background-color: rgba(0, 0, 0, 0); border: 0px solid black; border-radius: 0px; image: url("{0}");}}
			SXLogoButton:focus {{background-color: rgba(0, 0, 0, 0); border: 0px solid grey; border-radius: 0px; image: url("{0}");}}
			SXLogoButton:pressed {{background-color: rgba(0, 0, 0, 0); border: 0px solid red; border-radius: 0px; image: url("{0}");}}
			""".format(logo_file_path)

		# Set style and add click event
		self.setStyleSheet(self.customButtonStyle)

		# Add ToolTip
		self.setToolTip('Help')

# ========================================================================================
class SXPictogramButton(QPushButton):
	def __init__(self, pictogram_name, pictogram_file_path, parent = None):
		super(SXPictogramButton, self).__init__(parent)

		# print "MRK_DEBUG: pictogram_file_path = %s" % pictogram_file_path
		# print "MRK_DEBUG: os.path.exists(logo_file_path) %s" % os.path.exists(pictogram_file_path)

		# Width of pictogram image
		pictogram_width = SXLookFeelConst.sxmenu_item_btn_width

		# Style of widget
		self.setFixedSize(pictogram_width, pictogram_width)
		self.customButtonStyle = """
			SXPictogramButton {{background-color: rgba(0, 0, 0, 0); border: 2px solid rgba(0, 0, 0, 0); border-radius: {1}px; image: url("{0}");}}
			SXPictogramButton:focus {{background-color: rgba(0, 0, 0, 0); border: 2px solid grey; border-radius: {1}px; image: url("{0}");}}
			SXPictogramButton:pressed {{background-color: rgba(0, 0, 0, 0); border: 2px solid rgb(153, 153, 153); border-radius: {1}px; image: url("{0}");}}
			""".format(pictogram_file_path, pictogram_width / 6)
		self.customButtonStyleClicked = """
			SXPictogramButton:pressed {{background-color: rgba(0, 0, 0, 0); border: 2px solid rgb(153, 153, 153); border-radius: {1}px; image: url("{0}");}}
			SXPictogramButton {{background-color: rgba(0, 0, 0, 0); border: 2px solid rgb(220, 220, 220); border-radius: {1}px; image: url("{0}");}}
			""".format(pictogram_file_path, pictogram_width / 6)

		# Set style and add click event
		self.setStyleSheet(self.customButtonStyle)

		# Add tooltipp
		self.setToolTip('{0}{1}'.format(pictogram_name[0].upper(), pictogram_name[1:]))

class SXMenuItemBtnAreaWidget(QWidget):
	def __init__(self, sxconst_set, sxcmd_category_list, sxinfo, parent = None):
		super(SXMenuItemBtnAreaWidget, self).__init__(parent)

		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

		# Create widgets for pipeline command category button area and miscellaneous function button area
		sxcmd_category_btn_subarea_widget = self.create_sxmenu_item_btn_subarea_widget()
		misc_func_btn_subarea_widget = self.create_sxmenu_item_btn_subarea_widget()
		for sxcmd_category in sxcmd_category_list:
			if sxcmd_category.name != "sxc_utilities":
				self.add_sxmenu_item_btn_widget(sxcmd_category, sxcmd_category_btn_subarea_widget)
			else: # assert(sxcmd_category.name == "sxc_utilities")
				self.add_sxmenu_item_btn_widget(sxcmd_category, misc_func_btn_subarea_widget)
		self.add_sxmenu_item_btn_widget(sxconst_set, misc_func_btn_subarea_widget)

		global_layout = QVBoxLayout()
		global_layout.setContentsMargins(0, 0, 0, 0)

		sxmenu_item_btn_area_widget = QWidget(self)
		sxmenu_item_btn_area_widget.setObjectName('SXMenuItemBtnAreaWidget')
		sxmenu_item_btn_area_widget.setStyleSheet('QWidget#SXMenuItemBtnAreaWidget {background-color: rgba(0, 0, 0, 153);}')
		sxmenu_item_btn_area_widget.setFixedWidth(SXLookFeelConst.sxmenu_btn_area_min_width)

		sxmenu_item_btn_area_layout = QVBoxLayout()

		# Add widget of pipeline command category button area to layout
		sxmenu_item_btn_area_layout.addWidget(sxcmd_category_btn_subarea_widget)

		# Create and Add separator label
		layout_label = QHBoxLayout()
		line_label = QLabel(sxmenu_item_btn_area_widget)
		line_label.setFixedHeight(1)
		line_label.setFixedWidth(SXLookFeelConst.sxmenu_item_btn_width * 2)
		line_label.setStyleSheet('background-color: rgba(220, 220, 220, 100)')
		layout_label.addWidget(line_label)
		layout_label.setContentsMargins(0, 7, 0, 7)

		sxmenu_item_btn_area_layout.addLayout(layout_label)

		# Add widget of miscellaneous function button area to layout
		sxmenu_item_btn_area_layout.addWidget(misc_func_btn_subarea_widget)

		# Add stretch to make a space and keep sizes of the other widgets to be constant
		sxmenu_item_btn_area_layout.addStretch(1)

		# Add menu item button for application information
		sxmenu_item_btn_pictograph_file_path = '{0}sxgui_logo_sphire.png'.format(get_image_directory())
		sxmenu_item_btn = SXLogoButton(sxmenu_item_btn_pictograph_file_path)
		sxinfo.btn = sxmenu_item_btn

		sxmenu_item_btn_area_layout.addWidget(sxmenu_item_btn)

		# Set menu item button area layout to the widget
		sxmenu_item_btn_area_widget.setLayout(sxmenu_item_btn_area_layout)

		# self related settings
		global_layout.addWidget(sxmenu_item_btn_area_widget)
		self.setLayout(global_layout)

	def create_sxmenu_item_btn_subarea_widget(self):
		sxmenu_item_btn_subarea_widget = QWidget()

		grid_layout = QGridLayout()
		grid_layout.setSpacing(SXLookFeelConst.grid_distance)
		grid_layout.setContentsMargins(0, 0, 0, 0)

		sxmenu_item_btn_subarea_widget.setLayout(grid_layout)

		return sxmenu_item_btn_subarea_widget

	def add_sxmenu_item_btn_widget(self, sxmenu_item, sxmenu_item_btn_subarea_widget):
		assert(isinstance(sxmenu_item, SXmenu_item) == True) # Assuming the sxmenu_item is an instance of class SXmenu_item

		sxmenu_item_btn_pictograph_file_path = "{0}sxgui_pictograph_{1}.png".format(get_image_directory(), sxmenu_item.name.replace("sxc_", ""))
		sxmenu_item.btn = SXPictogramButton(sxmenu_item.name.replace("sxc_", ""), sxmenu_item_btn_pictograph_file_path, self)
		cur_widget_counts = sxmenu_item_btn_subarea_widget.layout().count()
		sxmenu_item_btn_subarea_widget.layout().addWidget(sxmenu_item.btn, cur_widget_counts // 2, cur_widget_counts % 2)

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

		self.gui_settings_file_path = "%s/gui_settings_%s.txt" % (self.sxcmd.get_category_dir_path(SXLookFeelConst.project_dir), self.sxcmd.get_mode_name_for("file_path"))

		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

		# Set grid layout
		grid_layout = QGridLayout(self)
		# grid_layout.setMargin(SXLookFeelConst.grid_margin)
		# grid_layout.setSpacing(SXLookFeelConst.grid_spacing)

		self.setAutoFillBackground(True)
		palette = QPalette()
		palette.setBrush(QPalette.Background, QBrush(SXLookFeelConst.sxcmd_widget_bg_color))
		self.setPalette(palette)

		self.sxcmd_tab_main = SXCmdTab("Main", self)
		self.sxcmd_tab_advance = SXCmdTab("Advanced", self)
		tab_widget = QTabWidget()
		tab_widget.insertTab(0, self.sxcmd_tab_main, self.sxcmd_tab_main.name)
		tab_widget.insertTab(1, self.sxcmd_tab_advance, self.sxcmd_tab_advance.name)
		tab_widget.setAutoFillBackground(True)
		tab_widget.setStyleSheet("""QTabWidget::pane {
			border-top: 2px solid #C2C7CB;
			position: absolute;
			top: -0.5em;
		}

		QTabWidget::tab-bar {
			alignment: center;
		}

		QTabBar::tab {
			background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1,
					stop: 0 #E1E1E1, stop: 0.4 #DDDDDD,
					stop: 0.5 #D8D8D8, stop: 1.0 #D3D3D3);
			border: 2px solid #C4C4C3;
			border-bottom-color: #C2C7CB; /* same as the pane color */
			border-top-left-radius: 4px;
			border-top-right-radius: 4px;
			min-width: 8ex;
			padding: 2px;
		}

		QTabBar::tab:selected, QTabBar::tab:hover {
			background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1,
					stop: 0 #fafafa, stop: 0.4 #f4f4f4,
					stop: 0.5 #e7e7e7, stop: 1.0 #fafafa);
		}

		QTabBar::tab:selected {
			border-color: #9B9B9B;
			border-bottom-color: #C2C7CB; /* same as pane color */
		}""")
		palette = tab_widget.palette()
		palette.setBrush(QPalette.Background, QBrush(SXLookFeelConst.sxcmd_widget_bg_color))
		tab_widget.setPalette(palette)
		grid_layout.addWidget(tab_widget, 0, 0)

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
						QMessageBox.warning(self, "Invalid parameter value", "Exteranl File Path (%s) should include the python script extension (.py)." % (external_file_path))
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
							QMessageBox.warning(self, "Invalid parameter value", "Token (%s) of command (%s) is required. Please set the value for this." % (sxcmd_token.label, self.sxcmd.get_mode_name_for("message_output")))
							return ""
						### if sxcmd_token.widget.text() == sxcmd_token.default and sxcmd_token.is_required == False: # Do not add this token to command line
						# else: # assert(sxcmd_token.is_required == False) # Do not add to this command line
					else: # sxcmd_token.widget.text() != sxcmd_token.default
						### if sxcmd_token.widget.text() != sxcmd_token.default and sxcmd_token.is_required == True:  # Add this token to command line
						### if sxcmd_token.widget.text() != sxcmd_token.default and sxcmd_token.is_required == False: # Add this token to command line

						# For now, using line edit box for the other type
						widget_text = str(sxcmd_token.widget.text())
						if sxcmd_token.type not in ["int", "float", "apix", "wn", "box", "radius", "any_file_list", "any_image_list", "any_directory"]:
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
				#
				# DESIGN_NOTE: 2016/03/17 Toshio Moriya
				# The MPI policy below has changed!!! An example of this exception is sxcter.py.
				# Don't add --MPI flag if np == 1
				#
				# DESIGN_NOTE: 2015/10/27 Toshio Moriya
				# Since we now assume sx*.py exists in only MPI version, always add --MPI flag if necessary
				# This is not elegant but can be removed when --MPI flag is removed from all sx*.py scripts
				#
				if self.sxcmd.mpi_add_flag and np > 1:
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
						QMessageBox.warning(self, "Invalid parameter value", "\"%s\" must be larger than 0. Please check the setting" % (required_label))
						return ""

					valid_np = np
					if valid_np % required_divisor != 0:
						if valid_np < required_divisor:
							valid_np = required_divisor
						else:
							valid_np = valid_np - (valid_np % required_divisor)
						QMessageBox.warning(self, "Invalid parameter value", "The number of \"MPI processes\" (%d) is invalid. It MUST BE multiplicity of \"%s\" (%d). Please check the setting. A close valid number is %d." % (np, required_label, required_divisor,valid_np))
						return ""

			# else: assert(np == 1) # because the "MPI Processes" is disabled for sx*.py process which does not support mpi

			# Generate command line according to the case
			cmd_line = ""
			if self.sxcmd_tab_main.qsub_enable_checkbox.checkState() == Qt.Checked:
				# Case 1: queue submission is enabled (MPI can be supported or unsupported)
				# Create script for queue submission from a give template
				if os.path.exists(self.sxcmd_tab_main.qsub_script_edit.text()) != True:
					QMessageBox.warning(self, "Invalid parameter value", "Invalid file path for qsub script template (%s)." % (self.sxcmd_tab_main.qsub_script_edit.text()))
					return ""

				file_template = open(self.sxcmd_tab_main.qsub_script_edit.text(),"r")
				# Extract command line from qsub script template
				for line in file_template:
					if line.find("XXX_SXCMD_LINE_XXX") != -1:
						if np > 1:
							cmd_line = line
						else:
							cmd_line = "XXX_SXCMD_LINE_XXX"
						cmd_line = cmd_line.replace("XXX_SXCMD_LINE_XXX", sxcmd_line)
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
					if np > 1:
						cmd_line = "mpirun -np XXX_SXMPI_NPROC_XXX XXX_SXCMD_LINE_XXX"
					else:
						cmd_line = "XXX_SXCMD_LINE_XXX"
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
					if os.path.exists(sxcmd_token.widget.text()) or db_check_dict(str(sxcmd_token.widget.text())):
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
					QMessageBox.warning(self, "Invalid parameter value", "Invalid file path for qsub script template (%s). Aborting execution ..." % (template_file_path))
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
			if os.path.exists(self.sxcmd.get_category_dir_path(SXLookFeelConst.project_dir)) == False:
				os.mkdir(self.sxcmd.get_category_dir_path(SXLookFeelConst.project_dir))
			self.write_params(self.gui_settings_file_path)
		# else: SX command line is be empty because an error happens in generate_cmd_line. Let's do nothing

	def print_cmd_line(self):
		# Generate command line
		cmd_line = self.generate_cmd_line()
		if cmd_line:
			message_line = "Generated the following command line:"
			print message_line
			print cmd_line
			QtGui.QMessageBox.information(self, "Information","%s \n\n%s" % (message_line, cmd_line))

			# Save the current state of GUI settings
			if os.path.exists(self.sxcmd.get_category_dir_path(SXLookFeelConst.project_dir)) == False:
				os.mkdir(self.sxcmd.get_category_dir_path(SXLookFeelConst.project_dir))
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
		# Write Qsub parameters
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
						self.sxcmd_tab_main.qsub_enable_checkbox.setChecked(Qt.Checked)
					else: # assert(val_str_in == "NO")
						self.sxcmd_tab_main.qsub_enable_checkbox.setChecked(Qt.Unchecked)
					# self.sxcmd_tab_main.set_qsub_enable_state() # Somehow this place does not paint the text boxes upon application startup
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
						QMessageBox.warning(self, "Invalid Parameter File Format", "Command token entry should start from \"%s\" for key base name in line (%s) of file (%s). The format of this file might be corrupted. Please save the paramater file again." % (target_operator, line_in, file_path_in))
					label_in = label_in[item_tail + len(target_operator):].strip() # Get the rest of line
					target_operator = ">"
					item_tail = label_in.find(target_operator)
					if item_tail == -1:
						QMessageBox.warning(self, "Invalid Parameter File Format", "Command token entry should have \"%s\" closing key base name in line (%s) of file (%s). The format of this file might be corrupted. Please save the paramater file again." % (target_operator, line_in, file_path_in))
					key_base = label_in[0:item_tail]
					# Get corresponding cmd_token
					if key_base not in self.sxcmd.token_dict.keys():
						QMessageBox.warning(self, "Invalid Parameter File Format", "Invalid base name of command token \"%s\" is found in line (%s) of file (%s). This parameter file might be imcompatible with the current version. Please save the paramater file again." % (key_base, line_in, file_path_in))
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
			QMessageBox.warning(self, "Fail to load parameters", "The specified file is not parameter file for %s." % self.sxcmd.get_mode_name_for("human"))

		file_in.close()

	def save_params(self):
		file_path_out = str(QFileDialog.getSaveFileName(self, "Save Parameters", SXLookFeelConst.file_dialog_dir, options = QFileDialog.DontUseNativeDialog))
		if file_path_out != "":
			self.write_params(file_path_out)

	def load_params(self):
		file_path_in = str(QFileDialog.getOpenFileName(self, "Load parameters", SXLookFeelConst.file_dialog_dir, options = QFileDialog.DontUseNativeDialog))
		if file_path_in != "":
			self.read_params(file_path_in)
			self.sxcmd_tab_main.set_qsub_enable_state()

	def select_file(self, target_widget, file_format = ""):
		file_path = ""
		if file_format == "bdb":
			file_path = str(QFileDialog.getOpenFileName(self, "Select BDB File", SXLookFeelConst.file_dialog_dir, "BDB files (*.bdb)", options = QFileDialog.DontUseNativeDialog))
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
				file_path = translate_to_bdb_path(file_path)
		elif file_format == "py":
			file_path = str(QFileDialog.getOpenFileName(self, "Select Python File", SXLookFeelConst.file_dialog_dir, "PY files (*.py)", options = QFileDialog.DontUseNativeDialog))
			# Use full path
		elif file_format == "pdb":
			file_path = str(QFileDialog.getOpenFileName(self, "Select PDB File", SXLookFeelConst.file_dialog_dir, "PDB files (*.pdb *.pdb*)", options = QFileDialog.DontUseNativeDialog))
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "mrc":
			file_path = str(QFileDialog.getOpenFileName(self, "Select MRC File", SXLookFeelConst.file_dialog_dir, "MRC files (*.mrc *.mrcs)", options = QFileDialog.DontUseNativeDialog))
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "exe":
			file_path = str(QFileDialog.getOpenFileName(self, "Select EXE File", SXLookFeelConst.file_dialog_dir, "EXE files (*.exe );; All files (*)", options = QFileDialog.DontUseNativeDialog))
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "any_file_list" or file_format == "any_image_list":
			file_path_list = QFileDialog.getOpenFileNames(self, "Select Files", SXLookFeelConst.file_dialog_dir, "All files (*)", options = QFileDialog.DontUseNativeDialog)
			# Use relative path.
			for a_file_path in file_path_list:
				file_path += SXLookFeelConst.format_path(str(a_file_path)) + " "
		else:
			if file_format:
				file_path = str(QFileDialog.getOpenFileName(self, "Select %s File" % (file_format.upper()), SXLookFeelConst.file_dialog_dir, "%s files (*.%s)"  % (file_format.upper(), file_format), options = QFileDialog.DontUseNativeDialog))
			else:
				file_path = str(QFileDialog.getOpenFileName(self, "Select File", SXLookFeelConst.file_dialog_dir, "All files (*)", options = QFileDialog.DontUseNativeDialog))
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)

		if file_path != "":
			target_widget.setText(file_path)

	def select_dir(self, target_widget):
		dir_path = str(QFileDialog.getExistingDirectory(self, "Select Directory", SXLookFeelConst.file_dialog_dir, options = QFileDialog.ShowDirsOnly | QFileDialog.DontResolveSymlinks | QFileDialog.DontUseNativeDialog))
		if dir_path != "":
			# Use relative path.
			target_widget.setText(SXLookFeelConst.format_path(dir_path))

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

		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# local constants
		required_cmd_token_restore_tooltip = "Please enter the value manually"
		const_cmd_token_restore_tooltip = "Retrieve the registed constant value for this parameter"
		default_cmd_token_restore_tooltip = "Retrieve this default value"

		# Setting for layout
		grid_row_origin = 0; grid_col_origin = 0
		title_row_span = 1; title_col_span = 2
		short_info_row_span = 1; short_info_col_span = 5
		func_btn_row_span = 1; func_btn_col_span = 2
		token_label_row_span = 1; token_label_col_span = 4
		token_widget_row_span = 1; token_widget_col_span = 1
		cmd_frame_row_span = 32; cmd_frame_col_span = 7

		title_label_min_width = 180 # title_label_min_width = 150
		title_label_min_height = 40 #title_label_min_height = 80
		short_info_min_width = 260 # short_info_min_width = 360
		short_info_min_height = 40 # short_info_min_height = 80
		func_btn_min_width = 150
		btn_min_width = 300
		token_label_min_width = 300 # token_label_min_width = 360
		token_widget_min_width = 120
		mpi_label_min_width = 100

		# Setup global layout
		global_layout = QVBoxLayout(self)
		global_layout.setContentsMargins(0,0,0,0)
		global_layout.setSpacing(0)
		# Setup scroll area and its widget
		scroll_area = QScrollArea()
		# scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
		# scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn) # MRK_DEBUG: Useful during designing layout
		scroll_area.setWidgetResizable(True)
		scroll_area_widget = QWidget(scroll_area)
		# Setup scroll widget and its background color
		scroll_area.setStyleSheet("QScrollArea {background-color:transparent;}");
		### scroll_area_widget.setStyleSheet("background-color:transparent;");
		scroll_area_widget.setAutoFillBackground(True)
		scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
		palette = QPalette()
		palette.setBrush(QPalette.Background, QBrush(SXLookFeelConst.sxcmd_tab_bg_color))
		scroll_area_widget.setPalette(palette)
		# Register the widget to scroll area
		scroll_area.setWidget(scroll_area_widget)
		# Register the scroll area to the global layout
		global_layout.addWidget(scroll_area)

		# Setup other layouts
		scroll_layout = QVBoxLayout(scroll_area_widget)
		scroll_layout.setContentsMargins(0,0,0,0)
		title_hbox = QHBoxLayout()
		title_layout = QGridLayout()
		title_layout.setMargin(SXLookFeelConst.grid_margin)
		title_layout.setSpacing(SXLookFeelConst.grid_spacing)
#		title_layout.setColumnMinimumWidth(grid_col_origin + token_label_col_span, token_widget_min_width)
#		title_layout.setColumnMinimumWidth(grid_col_origin + token_label_col_span + token_widget_col_span, token_widget_min_width)
#		title_layout.setColumnMinimumWidth(grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_min_width)
#		title_layout.setColumnMinimumWidth(grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_min_width)
		grid_layout = QGridLayout()
		grid_layout.setMargin(SXLookFeelConst.grid_margin)
		grid_layout.setSpacing(SXLookFeelConst.grid_spacing)
		grid_layout.setColumnMinimumWidth(grid_col_origin + token_label_col_span, token_widget_min_width)
		grid_layout.setColumnMinimumWidth(grid_col_origin + token_label_col_span + token_widget_col_span, token_widget_min_width)
		grid_layout.setColumnMinimumWidth(grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_min_width)
		grid_layout.setColumnMinimumWidth(grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_min_width)
		submit_layout = QGridLayout()
		submit_layout.setMargin(SXLookFeelConst.grid_margin)
		submit_layout.setSpacing(SXLookFeelConst.grid_spacing)
		submit_layout.setColumnMinimumWidth(grid_col_origin + token_label_col_span, token_widget_min_width)
		submit_layout.setColumnMinimumWidth(grid_col_origin + token_label_col_span + token_widget_col_span, token_widget_min_width)
		submit_layout.setColumnMinimumWidth(grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_min_width)
		submit_layout.setColumnMinimumWidth(grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_min_width)
		title_hbox.addLayout(title_layout)
#		title_hbox.addStretch(1)
		title_layout.setColumnStretch(grid_row_origin + token_label_col_span, title_layout.columnStretch(grid_row_origin+token_label_col_span) + 1)
		scroll_layout.addLayout(title_hbox)
		scroll_layout.addLayout(grid_layout)
		scroll_layout.addLayout(submit_layout)
		scroll_layout.addStretch(1)
		# # Give the columns of token label a higher priority to stretch relative to the others
		# for col_span in xrange(token_label_col_span):
		# 	grid_layout.setColumnStretch(grid_row_origin + col_span, grid_layout.columnStretch(grid_row_origin+col_span) + 1)

		# Define the tab frame within the tab layout
		# tab_frame = QFrame()
		# grid_layout.addWidget(tab_frame, grid_row_origin, grid_col_origin, cmd_frame_row_span, cmd_frame_col_span)

		# Start add command token widgets to the grid layout
		grid_row = grid_row_origin

		tab_group = self.name.lower()
		if tab_group == "main":
			# Set a label and its position in this tab
			temp_label = QLabel("<b>%s</b>" % (self.sxcmdwidget.sxcmd.get_mode_name_for("human")))
			temp_label.setMinimumWidth(title_label_min_width)
			temp_label.setMinimumHeight(title_label_min_height)
#			temp_label.setFixedWidth(title_label_min_width)
#			temp_label.setFixedHeight(title_label_min_height)
			title_layout.addWidget(temp_label, grid_row, grid_col_origin, title_row_span, title_col_span)

			#
			# NOTE: 2015/11/17 Toshio Moriya
			# Necessary to separate "<b>%s</b>" from the information for avoiding to invoke the tag interpretations of string
			# e.g. < becomes the escape character
			#
			temp_label = QLabel("%s" % (self.sxcmdwidget.sxcmd.short_info))
			temp_label.setWordWrap(True)
			temp_label.setMinimumWidth(short_info_min_width)
			temp_label.setMinimumHeight(short_info_min_height)
#			temp_label.setFixedHeight(short_info_min_height)
			title_layout.addWidget(temp_label, grid_row, grid_col_origin + title_col_span, short_info_row_span, short_info_col_span)

			grid_row += short_info_row_span

		elif tab_group == "advanced":
			# Set a label and its position in this tab
			temp_label = QLabel("<b>%s</b>" % (self.sxcmdwidget.sxcmd.get_mode_name_for("human")))
			temp_label.setMinimumWidth(title_label_min_width)
			temp_label.setMinimumHeight(title_label_min_height)
#			temp_label.setFixedWidth(title_label_min_width)
#			temp_label.setFixedHeight(title_label_min_height)
			title_layout.addWidget(temp_label, grid_row, grid_col_origin, title_row_span, title_col_span)

			temp_label = QLabel("Set advanced parameters", self)
			temp_label.setWordWrap(True)
			temp_label.setMinimumWidth(short_info_min_width)
			temp_label.setMinimumHeight(short_info_min_height)
#			temp_label.setFixedHeight(short_info_min_height)
			title_layout.addWidget(temp_label, grid_row, grid_col_origin + title_col_span, short_info_row_span, short_info_col_span)

		# Add space
		grid_row += 2

		# Add widget for editing command args and options
		for cmd_token in self.sxcmdwidget.sxcmd.token_list:
			if cmd_token.group == tab_group:

				# First, handle very special cases
				if cmd_token.type == "function":
					n_widgets = 2 # function type has two line edit boxes
					cmd_token_widget = [None] * n_widgets
					cmd_token_restore_widget = [None] * n_widgets

					# Define custom style for restore widgets
					custom_style = "QPushButton {color:gray; }"

					# Create widgets for user function name
					widget_index = 0
					temp_label = QLabel(cmd_token.label[widget_index])
					temp_label.setMinimumWidth(token_label_min_width)
					grid_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)

					assert(cmd_token.is_required == False)
					cmd_token_restore_widget[widget_index] = QPushButton("%s" % cmd_token.restore[widget_index])
					cmd_token_restore_widget[widget_index].setStyleSheet(custom_style)
					cmd_token_restore_widget[widget_index].setToolTip('<FONT>'+default_cmd_token_restore_tooltip+'</FONT>')
					grid_layout.addWidget(cmd_token_restore_widget[widget_index], grid_row, grid_col_origin + token_label_col_span, token_widget_row_span, token_widget_col_span)

					# cmd_token_widget[widget_index] = QLineEdit(self)
					cmd_token_widget[widget_index] = QLineEdit()
					cmd_token_widget[widget_index].setText(cmd_token.restore[widget_index])
					cmd_token_widget[widget_index].setToolTip('<FONT>'+cmd_token.help[widget_index]+'</FONT>')
					grid_layout.addWidget(cmd_token_widget[widget_index], grid_row, grid_col_origin + token_label_col_span + token_widget_col_span, token_widget_row_span, token_widget_col_span)

					self.connect(cmd_token_restore_widget[widget_index], SIGNAL("clicked()"), partial(self.handle_restore_widget_event, cmd_token, widget_index))

					grid_row +=  1

					# Create widgets for external file path containing above user function
					widget_index = 1
					temp_label = QLabel(cmd_token.label[widget_index])
					grid_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)

					assert(cmd_token.is_required == False)
					cmd_token_restore_widget[widget_index] = QPushButton("%s" % cmd_token.restore[widget_index])
					cmd_token_restore_widget[widget_index].setStyleSheet(custom_style)
					cmd_token_restore_widget[widget_index].setToolTip('<FONT>'+default_cmd_token_restore_tooltip+'</FONT>')
					grid_layout.addWidget(cmd_token_restore_widget[widget_index], grid_row, grid_col_origin + token_label_col_span, token_widget_row_span, token_widget_col_span)

					cmd_token_widget[widget_index] = QLineEdit()
					cmd_token_widget[widget_index].setText(cmd_token.restore[widget_index]) # Because default user functions is internal
					cmd_token_widget[widget_index].setToolTip('<FONT>'+cmd_token.help[widget_index]+'</FONT>')
					grid_layout.addWidget(cmd_token_widget[widget_index], grid_row, grid_col_origin + token_label_col_span + token_widget_col_span, token_widget_row_span, token_widget_col_span)

					self.connect(cmd_token_restore_widget[widget_index], SIGNAL("clicked()"), partial(self.handle_restore_widget_event, cmd_token, widget_index))

					file_format = "py"
					temp_btn = QPushButton("Select Script")
					temp_btn.setToolTip('<FONT>'+"Display open file dailog to select .%s python script file</FONT>" % file_format)
					grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
					self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget[widget_index], file_format))

					grid_row +=  1

#					temp_label = QLabel(cmd_token.help[widget_index])
#					grid_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)
#
#					grid_row +=  1

				# Then, handle the other cases
				else:
					# Create label widget
					temp_label = QLabel(cmd_token.label)
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
						cmd_token_restore_widget = QPushButton("%s" % btn_name)
						cmd_token_restore_widget.setStyleSheet(custom_style)
						cmd_token_restore_widget.setEnabled(is_btn_enable)
						grid_layout.addWidget(cmd_token_restore_widget, grid_row, grid_col_origin + token_label_col_span, token_widget_row_span, token_widget_col_span)

						# construct new widget(s) for this command token
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
						cmd_token_restore_widget = QPushButton("%s" % btn_name)
						cmd_token_restore_widget.setStyleSheet(custom_style)
						cmd_token_restore_widget.setEnabled(is_btn_enable)
						grid_layout.addWidget(cmd_token_restore_widget, grid_row, grid_col_origin + token_label_col_span, token_widget_row_span, token_widget_col_span)

						cmd_token_widget = QLineEdit()
						cmd_token_widget.setText(cmd_token.restore)
						grid_layout.addWidget(cmd_token_widget, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span, token_widget_row_span, token_widget_col_span)

						self.connect(cmd_token_restore_widget, SIGNAL("clicked()"), partial(self.handle_restore_widget_event, cmd_token))

						if cmd_token.type == "image":
							file_format = "hdf"
							temp_btn = QPushButton("Select .%s" % file_format)
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select .%s format image file</FONT>" % file_format)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "bdb"
							temp_btn = QPushButton("Select .%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select .%s format image file</FONT>" % file_format)
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
						elif cmd_token.type == "any_image":
							temp_btn = QPushButton("Select Image")
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select standard format image file (e.g. .hdf, .mrc)</FONT>")
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget))
							file_format = "invisible"
							temp_btn = QPushButton("Select .%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select .%s format image file</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "any_micrograph":
							temp_btn = QPushButton("Select Image")
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select standard format image file (e.g. .hdf, .mrc)</FONT>")
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget))
							file_format = "txt"
							temp_btn = QPushButton("Select .%s" % file_format)
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select .%s parameter file</FONT>" % file_format)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
						elif cmd_token.type == "any_file_list":
							temp_btn = QPushButton("Select Files")
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select files (e.g. *.*)</FONT>")
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, cmd_token.type))
							file_format = "bdb"
							temp_btn = QPushButton("Select .%s" % file_format)
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select .%s format image file</FONT>" % file_format)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
						elif cmd_token.type == "any_image_list":
							temp_btn = QPushButton("Select Images")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select standard format image files (e.g. .hdf, .mrc)</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, cmd_token.type))
							file_format = "invisible"
							temp_btn = QPushButton("Select .%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select .%s format image file</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "bdb" or cmd_token.key_base == "makevstack":
							file_format = "bdb"
							temp_btn = QPushButton("Select .%s" % file_format)
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select .%s format image file</FONT>" % file_format)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "invisible"
							temp_btn = QPushButton("Select .%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select .%s format image file</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "pdb":
							file_format = "pdb"
							temp_btn = QPushButton("Select .%s" % file_format)
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select .%s format image file</FONT>" % file_format)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "invisible"
							temp_btn = QPushButton("Select .%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select .%s format image file</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "hdf":
							file_format = "hdf"
							temp_btn = QPushButton("Select .%s" % file_format)
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select .%s format image file</FONT>" % file_format)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "invisible"
							temp_btn = QPushButton("Select .%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select .%s format image file</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "mrc":
							file_format = "mrc"
							temp_btn = QPushButton("Select .%s" % file_format)
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select .%s format image file</FONT>" % file_format)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "invisible"
							temp_btn = QPushButton("Select .%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select .%s format image file</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "parameters":
							temp_btn = QPushButton("Select Parameter")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select parameter file (e.g. .txt)</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget))
							file_format = "invisible"
							temp_btn = QPushButton("Select .%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select .%s format image file</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "txt":
							file_format = "txt"
							temp_btn = QPushButton("Select .%s" % file_format)
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select .%s parameter file</FONT>" % file_format)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "invisible"
							temp_btn = QPushButton("Select .%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select .%s format image file</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "exe":
							file_format = "exe"
							temp_btn = QPushButton("Select .%s" % file_format)
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select .%s parameter file</FONT>" % file_format)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "invisible"
							temp_btn = QPushButton("Select .%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select .%s format image file</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "any_file":
							temp_btn = QPushButton("Select File")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select file (e.g. *.*)</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget))
							file_format = "invisible"
							temp_btn = QPushButton("Select .%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select .%s format image file</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "directory" or cmd_token.type == "any_directory":
							temp_btn = QPushButton("Select directory")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display select directory dailog"+'</FONT>')
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_dir, cmd_token_widget))
							file_format = "invisible"
							temp_btn = QPushButton("Select .%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select .%s format image file</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						# elif cmd_token.type == "output":
						# else:
						# 	if cmd_token.type not in ["int", "float", "string", "apix", "wn", "box", "radius", "sym"]: ERROR("Logical Error: Encountered unsupported type (%s). Consult with the developer."  % cmd_token.type, "%s in %s" % (__name__, os.path.basename(__file__)))

					cmd_token_widget.setToolTip('<FONT>'+cmd_token.help+'</FONT>')
					cmd_token_restore_widget.setToolTip('<FONT>'+cmd_token_restore_tooltip+'</FONT>')

					grid_row += 1

				# Register this widget
				cmd_token.widget = cmd_token_widget
				cmd_token.restore_widget = cmd_token_restore_widget

		if tab_group == "main":
			# Add space
			grid_row += 1

			# Add gui components for MPI related paramaters
			temp_label = QLabel("MPI processors")
			temp_label.setMinimumWidth(token_label_min_width)
			submit_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)

			# self.mpi_nproc_edit = QLineEdit(self)
			self.mpi_nproc_edit = QLineEdit()
			self.mpi_nproc_edit.setText("1")
			self.mpi_nproc_edit.setToolTip('<FONT>'+"Number of processors to use. default is single processor mode"+'</FONT>')
			submit_layout.addWidget(self.mpi_nproc_edit, grid_row, grid_col_origin + token_label_col_span, token_widget_row_span, token_widget_col_span)

			# Add save paramaters button
			self.save_params_btn = QPushButton("Save parameters")
			self.save_params_btn.setMinimumWidth(btn_min_width)
			self.save_params_btn.setToolTip('<FONT>'+"Save gui parameter settings"+'</FONT>')
			self.connect(self.save_params_btn, SIGNAL("clicked()"), self.sxcmdwidget.save_params)
			submit_layout.addWidget(self.save_params_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span*2)

			grid_row += 1

			temp_label = QLabel("MPI command line template")
			temp_label.setMinimumWidth(token_label_min_width)
			submit_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)

			self.mpi_cmd_line_edit = QLineEdit()
			self.mpi_cmd_line_edit.setText("")
			self.mpi_cmd_line_edit.setToolTip('<FONT>'+"Template of MPI command line (e.g. \"mpirun -np XXX_SXMPI_NPROC_XXX --host n0,n1,n2 XXX_SXCMD_LINE_XXX\"). if empty, use \"mpirun -np XXX_SXMPI_NPROC_XXX XXX_SXCMD_LINE_XXX\"</FONT>")
			submit_layout.addWidget(self.mpi_cmd_line_edit, grid_row, grid_col_origin + token_label_col_span, token_widget_row_span, token_widget_col_span)

			# Add load paramaters button
			self.load_params_btn = QPushButton("Load parameters")
			self.load_params_btn.setMinimumWidth(btn_min_width)
			self.load_params_btn.setToolTip('<FONT>'+"Load gui parameter settings to retrieve a previously-saved one"+'</FONT>')
			self.connect(self.load_params_btn, SIGNAL("clicked()"), self.sxcmdwidget.load_params)
			submit_layout.addWidget(self.load_params_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span*2)

			grid_row += 1

			# If MPI is not supported, disable this widget
			self.set_text_entry_widget_enable_state(self.mpi_nproc_edit, self.sxcmdwidget.sxcmd.mpi_support)
			self.set_text_entry_widget_enable_state(self.mpi_cmd_line_edit, self.sxcmdwidget.sxcmd.mpi_support)

			# Add gui components for queue submission (qsub)
			is_qsub_enabled = False
			temp_label = QLabel("Submit job to queue")
			temp_label.setMinimumWidth(token_label_min_width)
			submit_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)

			self.qsub_enable_checkbox = QCheckBox("")
			if is_qsub_enabled == True:
				self.qsub_enable_checkbox.setCheckState(Qt.Checked)
			else: # assert(is_qsub_enabled == False)
				self.qsub_enable_checkbox.setCheckState(Qt.Unchecked)
			self.qsub_enable_checkbox.setToolTip('<FONT>'+"Submit job to queue"+'</FONT>')
			self.qsub_enable_checkbox.stateChanged.connect(self.set_qsub_enable_state) # To control enable state of the following qsub related widgets
			self.qsub_enable_checkbox.setEnabled(self.sxcmdwidget.sxcmd.is_submittable)
			submit_layout.addWidget(self.qsub_enable_checkbox, grid_row, grid_col_origin + token_label_col_span, token_widget_row_span, token_widget_col_span)

			grid_row += 1

			temp_label = QLabel("Job name")
			temp_label.setMinimumWidth(token_label_min_width)
			submit_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)

			self.qsub_job_name_edit = QLineEdit()
			if self.sxcmdwidget.sxcmd.is_submittable == True:
				self.qsub_job_name_edit.setText(self.sxcmdwidget.sxcmd.get_mode_name_for("file_path"))
			else: # assert(self.sxcmdwidget.sxcmd.is_submittable == False)
				assert(self.sxcmdwidget.sxcmd.mpi_support == False)
				self.qsub_job_name_edit.setText("N/A")
			self.qsub_job_name_edit.setToolTip('<FONT>'+"Name of this job"+'</FONT>')
			submit_layout.addWidget(self.qsub_job_name_edit, grid_row, grid_col_origin + token_label_col_span, token_widget_row_span, token_widget_col_span)

			grid_row += 1

			temp_label = QLabel("Submission command")
			temp_label.setMinimumWidth(token_label_min_width)
			submit_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)

			self.qsub_cmd_edit = QLineEdit()
			if self.sxcmdwidget.sxcmd.is_submittable == True:
				self.qsub_cmd_edit.setText("qsub")
			else: # assert(self.sxcmdwidget.sxcmd.is_submittable == False)
				assert(self.sxcmdwidget.sxcmd.mpi_support == False)
				self.qsub_cmd_edit.setText("N/A")
			self.qsub_cmd_edit.setToolTip('<FONT>'+"Name of submission command to queue job"+'</FONT>')
			submit_layout.addWidget(self.qsub_cmd_edit, grid_row, grid_col_origin + token_label_col_span, token_widget_row_span, token_widget_col_span)

			self.cmd_line_btn = QPushButton("Generate command line")
			self.cmd_line_btn.setMinimumWidth(btn_min_width)
			self.cmd_line_btn.setToolTip('<FONT>'+"Generate command line from gui parameter settings and automatically save settings"+'</FONT>')
			self.connect(self.cmd_line_btn, SIGNAL("clicked()"), self.sxcmdwidget.print_cmd_line)
			submit_layout.addWidget(self.cmd_line_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span*2)

			grid_row += 1

			temp_label = QLabel("Submission script template")
			temp_label.setMinimumWidth(token_label_min_width)
			submit_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)

			self.qsub_script_edit = QLineEdit()
			if self.sxcmdwidget.sxcmd.is_submittable == True:
				self.qsub_script_edit.setText("msgui_qsub.sh")
			else: # assert(self.sxcmdwidget.sxcmd.is_submittable == False)
				assert(self.sxcmdwidget.sxcmd.mpi_support == False)
				self.qsub_script_edit.setText("N/A")
			self.qsub_script_edit.setToolTip('<FONT>'+"File name of submission script template (e.g. $EMAN2DIR/bin/msgui_qsub.sh)"+'</FONT>')
			submit_layout.addWidget(self.qsub_script_edit, grid_row, grid_col_origin + token_label_col_span, token_widget_row_span, token_widget_col_span)

			self.qsub_script_open_btn = QPushButton("Select Template")
			self.qsub_script_open_btn.setMinimumWidth(func_btn_min_width)
			self.qsub_script_open_btn.setToolTip('<FONT>'+"Display open file dailog to select job submission script template file"+'</FONT>')
			self.connect(self.qsub_script_open_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, self.qsub_script_edit))
			submit_layout.addWidget(self.qsub_script_open_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span, token_widget_row_span, token_widget_col_span)

			# Add a run button
			self.execute_btn = QPushButton("Run %s" % self.sxcmdwidget.sxcmd.get_mode_name_for("human"))
			# make 3D textured push button look
			custom_style = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0)} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084)} QPushButton:focus {font: bold; color: #000;border: 2px solid #8D0;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0)}"
			self.execute_btn.setStyleSheet(custom_style)
			self.execute_btn.setMinimumWidth(btn_min_width)
			self.execute_btn.setToolTip('<FONT>'+"Run %s and automatically save gui parameter settings</FONT>" % self.sxcmdwidget.sxcmd.get_mode_name_for("human"))
			self.connect(self.execute_btn, SIGNAL("clicked()"), self.sxcmdwidget.execute_cmd_line)
			submit_layout.addWidget(self.execute_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span*2)

			grid_row += 1

			# Initialize enable state of qsub related widgets
			self.set_qsub_enable_state()

	def set_text_entry_widget_enable_state(self, widget, is_enabled):
		# Set enable state and background color of text entry widget according to enable state
		default_palette = QPalette()
		bg_color = default_palette.color(QPalette.Inactive, QPalette.Base)
		if is_enabled == False:
			bg_color = default_palette.color(QPalette.Disabled, QPalette.Base)

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
# Command Category Widget (opened by class SXMainWindow)
class SXCmdCategoryWidget(QWidget):
	def __init__(self, sxconst_set, sxcmd_category, parent = None):
		super(SXCmdCategoryWidget, self).__init__(parent)

		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.sxconst_set = sxconst_set
		self.sxcmd_category = sxcmd_category
		self.cur_sxcmd = None

		# Layout constants
		self.sxcmd_btn_row_span = 1
		self.sxcmd_btn_col_span = 1

		self.sxcmd_btn_area_row_span = self.sxcmd_btn_row_span * SXLookFeelConst.expected_cmd_counts
		self.sxcmd_btn_area_col_span = self.sxcmd_btn_col_span

		self.sxcmd_widget_area_row_span = self.sxcmd_btn_area_row_span
		self.sxcmd_widget_area_col_span = 1

		self.grid_row_origin = 0
		self.grid_col_origin = 0

		# Layout variables
		self.grid_layout = None # grid layout

		self.grid_row = self.grid_row_origin # Keep current row
		self.grid_col = self.grid_col_origin # keep current column

		self.sxcmd_btn_group = None
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

		# --------------------------------------------------------------------------------
		# Setup Window Layout
		# --------------------------------------------------------------------------------
		self.setup_layout(QBrush(SXLookFeelConst.default_bg_color))

		# --------------------------------------------------------------------------------
		# Add SX Commands (sx*.py) associated widgets
		# --------------------------------------------------------------------------------
		self.add_sxcmd_widgets()

#		# --------------------------------------------------------------------------------
#		# Load the previously saved parameter setting of this sx command
#		# Override the registration of project constant parameter settings with the previously-saved one
#		# --------------------------------------------------------------------------------
#		for sxcmd in self.sxcmd_category.cmd_list:
#			if os.path.exists(sxcmd.widget.gui_settings_file_path):
#				sxcmd.widget.read_params(sxcmd.widget.gui_settings_file_path)

		# --------------------------------------------------------------------------------
		# Alway select the 1st entry of the command list upon startup
		# --------------------------------------------------------------------------------
		self.handle_sxcmd_btn_event(self.sxcmd_category.cmd_list[0])

	def setup_layout(self, background_brush):
		# Setup background color of this widget
		self.setAutoFillBackground(True)
		palette = QPalette()
		palette.setBrush(QPalette.Background, background_brush)
		self.setPalette(palette)

		# Setup grid layout in the scroll area
		self.grid_layout = QGridLayout()
		self.grid_layout.setMargin(SXLookFeelConst.grid_margin)
		self.grid_layout.setSpacing(SXLookFeelConst.grid_spacing)
		self.grid_layout.setColumnMinimumWidth(0, SXLookFeelConst.sxcmd_btn_area_min_width)
		# self.grid_layout.setColumnMinimumWidth(1, SXLookFeelConst.sxcmd_widget_area_min_width)
		# Give the column of the command settings area a higher stretch priority so that the other area does not stretch horizontally
		# self.grid_layout.setColumnStretch(self.grid_col_origin + self.sxcmd_btn_area_col_span, self.grid_layout.columnStretch(self.grid_col_origin + self.sxcmd_btn_area_col_span) + 1)

	# Add Pipeline SX Commands (sx*.py) associated widgets
	def add_sxcmd_widgets(self):
		self.sxcmd_btn_group = QButtonGroup()
		# self.sxcmd_btn_group.setExclusive(True) # NOTE: 2016/02/18 Toshio Moriya: Without QPushButton.setCheckable(True). This does not do anything. Let manually do this

		current_role = None
		self.stacked_layout = QStackedLayout()
		grid_box_layout = QVBoxLayout()
		grid_box_layout.addLayout(self.grid_layout)
		grid_box_layout.addStretch(1)
		global_layout = QHBoxLayout()
		global_layout.addLayout(grid_box_layout)
		global_layout.addLayout(self.stacked_layout, stretch=1)
		self.setLayout(global_layout)

		# Add SX Commands (sx*.py) associated widgets
		for sxcmd in self.sxcmd_category.cmd_list:
			if sxcmd.role != current_role:
				# Add title label and set position and font style
				label_text = ""
				if sxcmd.role == "sxr_pipe":
					label_text = "COMMANDS"
				elif sxcmd.role == "sxr_alt":
					label_text = "ALTERNATIVES"
				elif sxcmd.role == "sxr_util":
					label_text = "UTILITIES"
				else:
					label_text = "UNKNOWN"

				if current_role !=  None:
					self.grid_row += 1

				# title=QLabel("<span style=\'font-size:18pt; font-weight:600; color:#aa0000;\'><b>%s </b></span><span style=\'font-size:12pt; font-weight:60; color:#aa0000;\'>(shift-click for wiki)</span>" % label_text)
				title=QLabel("<span style=\'font-size:18pt; font-weight:600; color:#000000;\'><b>%s </b></span><span style=\'font-size:12pt; font-weight:60; color:#000000;\'>(shift-click for wiki)</span>" % label_text)
				self.grid_layout.addWidget(title, self.grid_row, self.grid_col_origin, self.sxcmd_btn_row_span, self.sxcmd_btn_col_span)

				self.grid_row += 1

				current_role = sxcmd.role

			# Add buttons for this sx*.py processe
			sxcmd.btn = QPushButton(sxcmd.label)
			# sxcmd.btn.setCheckable(True) # NOTE: 2016/02/18 Toshio Moriya: With this setting, we can not move the focus to the unchecked butttons... PyQt bug?
			sxcmd.btn.setToolTip('<FONT>'+sxcmd.short_info+'</FONT>')
			self.sxcmd_btn_group.addButton(sxcmd.btn)
			self.grid_layout.addWidget(sxcmd.btn, self.grid_row, self.grid_col_origin, self.sxcmd_btn_row_span, self.sxcmd_btn_col_span)

			# Create SXCmdWidget for this sx*.py processe
			sxcmd.widget = SXCmdWidget(self.sxconst_set, sxcmd)
			self.stacked_layout.addWidget(sxcmd.widget)

			# connect widget signals
			self.connect(sxcmd.btn, SIGNAL("clicked()"), partial(self.handle_sxcmd_btn_event, sxcmd))

			self.grid_row += 1

	def load_previous_session(self):
		for sxcmd in self.sxcmd_category.cmd_list:
			if os.path.exists(sxcmd.widget.gui_settings_file_path):
				sxcmd.widget.read_params(sxcmd.widget.gui_settings_file_path)

	def handle_sxcmd_btn_event(self, sxcmd):
		modifiers = QApplication.keyboardModifiers()
		if modifiers == Qt.ShiftModifier:
			os.system("python -m webbrowser %s%s" % (SPARX_DOCUMENTATION_WEBSITE, sxcmd.name))
			return

		if self.cur_sxcmd == sxcmd: return

		if self.cur_sxcmd != None:
			custom_style = "QPushButton {font: normal; color:black; }" # custom_style = "QPushButton {color:#000; }"
			self.cur_sxcmd.btn.setStyleSheet(custom_style)

		self.cur_sxcmd = sxcmd

		if self.cur_sxcmd != None:
			self.stacked_layout.setCurrentWidget(self.cur_sxcmd.widget)
			custom_style = "QPushButton {font: bold; color:blue; }" # custom_style = "QPushButton {font: bold; color:#8D0; }"
			self.cur_sxcmd.btn.setStyleSheet(custom_style)

	def quit_all_child_applications(self):
		# Quit all child applications
		for sxcmd in self.sxcmd_category.cmd_list:
			sxcmd.widget.quit_all_child_applications()

# ========================================================================================
# Layout of the project constants parameters widget; owned by the main window
class SXConstSetWidget(QWidget):
	def __init__(self, sxconst_set, sxcmd_category_list, parent=None):
		super(SXConstSetWidget, self).__init__(parent)

		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.sxconst_set = sxconst_set
		self.sxcmd_category_list = sxcmd_category_list

		self.gui_settings_file_path = "%s/gui_settings_project_settings.txt" % (SXLookFeelConst.project_dir)

		# Layout constants and variables
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
		func_btn_min_width = 50
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

		# Set the background color of this widget
		self.setAutoFillBackground(True)
		palette = QPalette()
		palette.setBrush(QPalette.Background, QBrush(SXLookFeelConst.default_bg_color))
		self.setPalette(palette)

		global_layout = QGridLayout()
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

		global_grid_row = global_row_origin

		# Start add title widgets to the grid layout
		header_grid_row = header_row_origin

		# Set a label and its position in this tab
		temp_label = QLabel("<b>%s</b>" % (self.sxconst_set.label))
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

			sxconst_register_widget = QPushButton("%s" % sxconst.register)
			sxconst_register_widget.setMinimumWidth(const_register_widget_min_width)
			custom_style = "QPushButton {color:green; }"
			sxconst_register_widget.setStyleSheet(custom_style)
			const_set_layout.addWidget(sxconst_register_widget, const_set_grid_row, const_set_row_origin + const_label_col_span, const_register_widget_row_span, const_register_widget_col_span)
			sxconst_register_widget.setToolTip('<FONT>'+"Retrieve this registered value to edit box"+'</FONT>')
			self.connect(sxconst_register_widget, SIGNAL("clicked()"), partial(self.handle_regster_widget_event, sxconst))

			sxconst_widget = QLineEdit()
			sxconst_widget.setMinimumWidth(const_widget_min_width)
			sxconst_widget.setText(sxconst.register)
			sxconst_widget.setToolTip('<FONT>'+sxconst.help+'</FONT>')
			const_set_layout.addWidget(sxconst_widget, const_set_grid_row, const_set_row_origin + const_label_col_span + const_register_widget_col_span, const_widget_row_span, const_widget_col_span)

			const_set_grid_row += 1

			# Register this widget
			sxconst.register_widget = sxconst_register_widget
			sxconst.widget = sxconst_widget

		# Add const set grid layout to global layout
		global_layout.addLayout(const_set_layout, global_grid_row, global_col_origin)
		# global_grid_row += 1

		# Start add buttons to the grid layout
		btn_grid_row = btn_row_origin

		# Add a register button
		self.execute_btn = QPushButton("Register settings")
		# make 3D textured push button look
		custom_style = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px}"
		self.execute_btn.setStyleSheet(custom_style)
		self.execute_btn.setMinimumWidth(func_btn_min_width * register_btn_col_span)
		self.execute_btn.setToolTip('<FONT>'+"Register project constant parameter settings to automatically set values to command arguments and options"+'</FONT>')
		self.connect(self.execute_btn, SIGNAL("clicked()"), self.register_const_set)
		btn_layout.addWidget(self.execute_btn, btn_grid_row, btn_col_origin, register_btn_row_span, register_btn_col_span)

		btn_grid_row += 1

		# Add save project constant parameter settings button
		self.save_consts_btn = QPushButton("Save settings")
		self.save_consts_btn.setMinimumWidth(func_btn_min_width)
		self.save_consts_btn.setToolTip('<FONT>'+"Save project constant parameter settings"+'</FONT>')
		self.connect(self.save_consts_btn, SIGNAL("clicked()"), self.save_consts)
		btn_layout.addWidget(self.save_consts_btn, btn_grid_row, btn_col_origin, func_btn_row_span, func_btn_col_span)

		# Add load project constant parameter settings button
		self.load_consts_btn = QPushButton("Load settings")
		self.load_consts_btn.setMinimumWidth(func_btn_min_width)
		self.load_consts_btn.setToolTip('<FONT>'+"Load project constant parameter settings to retrieve the previously-saved one"+'</FONT>')
		self.connect(self.load_consts_btn, SIGNAL("clicked()"), self.load_consts)
		btn_layout.addWidget(self.load_consts_btn, btn_grid_row, btn_col_origin + func_btn_col_span, func_btn_row_span, func_btn_col_span)

		btn_grid_row += 1

		# Add button grid layout to global layout
		# global_layout.addLayout(btn_layout, global_grid_row, global_col_origin) # Maybe later :)

		# Load the previously saved parameter setting of this sx command
		if os.path.exists(self.gui_settings_file_path):
			self.read_consts(self.gui_settings_file_path)

		# Layout for a constant size
		constant_height_layout = QVBoxLayout()
		constant_height_layout.addLayout(global_layout)
		constant_height_layout.addLayout(btn_layout)
		constant_height_layout.addStretch(1)
		constant_width_layout = QHBoxLayout(self)
		constant_width_layout.addLayout(constant_height_layout)
		constant_width_layout.addStretch(1)

	def handle_regster_widget_event(self, sxconst):
		sxconst.widget.setText(sxconst.register)

	def register_const_set(self):
		# Loop through all project constant parameters
		for sxconst in self.sxconst_set.list:
			sxconst.register = sxconst.widget.text()
			sxconst.register_widget.setText("%s" % sxconst.register)

		# Loop through all command categories
		for sxcmd_category in self.sxcmd_category_list:
			# Loop through all commands of this category
			for sxcmd in sxcmd_category.cmd_list:
				# Loop through all command tokens of this command
				for cmd_token in sxcmd.token_list:
					if cmd_token.type in self.sxconst_set.dict.keys():
						sxconst = self.sxconst_set.dict[cmd_token.type]
						cmd_token.restore = sxconst.register
						cmd_token.restore_widget.setText("%s" % cmd_token.restore)
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
		file_path_out = str(QFileDialog.getSaveFileName(self, "Save settings", SXLookFeelConst.file_dialog_dir, options = QFileDialog.DontUseNativeDialog))
		if file_path_out != "":
			self.write_consts(file_path_out)

	def load_consts(self):
		file_path_in = str(QFileDialog.getOpenFileName(self, "Load settings", SXLookFeelConst.file_dialog_dir, options = QFileDialog.DontUseNativeDialog))
		if file_path_in != "":
			self.read_consts(file_path_in)

# ========================================================================================
# Layout of the information widget; owned by the main window
class SXInfoWidget(QWidget):
	def __init__(self, parent = None):
		super(SXInfoWidget, self).__init__(parent)

		# Set the background color of this widget
		self.setAutoFillBackground(True)
		palette = QPalette()
		palette.setBrush(QPalette.Background, QBrush(SXLookFeelConst.sxinfo_widget_bg_color))
		self.setPalette(palette)

		label_row_span = 1; label_col_span = 3
		close_row_span = 1; close_col_span = 1
		spacer_min_width = 12

		grid_layout = QGridLayout(self)
		grid_layout.setMargin(SXLookFeelConst.grid_margin)
		grid_layout.setSpacing(SXLookFeelConst.grid_spacing)

		grid_col = 0
		grid_row = 0; grid_layout.setRowMinimumHeight(grid_row, spacer_min_width)
		grid_row += 10; temp_label=QLabel("<span style=\'font-size:18pt; font-weight:600; color:#ffffff;\'><b>SPHIRE GUI Prototype</b></span>")
		temp_label.setAlignment(Qt.AlignHCenter)
		grid_layout.addWidget(temp_label, grid_row, grid_col, label_row_span, label_col_span)
		grid_row += 1; temp_label=QLabel("<span style=\'font-size:18pt; font-weight:600; color:#ffffff;\'><b>Author: Toshio Moriya</b></span>")
		temp_label.setAlignment(Qt.AlignHCenter)
		grid_layout.addWidget(temp_label, grid_row, grid_col, label_row_span, label_col_span)
		grid_row += 1; grid_layout.setRowMinimumHeight(grid_row, spacer_min_width)
		grid_row += 1; temp_label=QLabel("<span style=\'font-size:18pt; font-weight:600; color:#ffffff;\'><b>For more information visit:%s </b></span>" % SPARX_DOCUMENTATION_WEBSITE)
		temp_label.setAlignment(Qt.AlignHCenter)
		grid_layout.addWidget(temp_label, grid_row, grid_col, label_row_span, label_col_span)
		grid_row += 1; grid_layout.setRowMinimumHeight(grid_row, spacer_min_width)
		grid_row += 1; grid_layout.setRowMinimumHeight(grid_row, spacer_min_width)

# ========================================================================================
# Main Window (started by class SXApplication)
class SXMainWindow(QMainWindow): # class SXMainWindow(QWidget):

	def __init__(self, parent = None):
		super(SXMainWindow, self).__init__(parent)

		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		self.sxinfo = None
		self.sxconst_set = None
		self.sxcmd_category_list = None

		self.cur_sxmenu_item = None
		self.sxmenu_item_widget_stacked_layout = None

		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

		# --------------------------------------------------------------------------------
		# Construct menu items
		# --------------------------------------------------------------------------------
		self.construct_sxinfo()              # Construct application information
		self.construct_sxconst_set()         # Construct project constant set for project settings
		self.construct_sxcmd_category_list() # Construct list of categorised sxscript objects (extracted from associated wiki documents)

		# --------------------------------------------------------------------------------
		# Setup Window Layout
		# --------------------------------------------------------------------------------
		background_image_file_path = '{0}sxgui_background.png'.format(get_image_directory())

		# Central widget
		central_widget = QWidget(self)
		central_widget.setObjectName('central')
		central_widget.setStyleSheet(
			'QWidget#central {{background-image: url("{0}")}}'.format(background_image_file_path)
			)
		self.setCentralWidget(central_widget)

		# Layout for central widget
		central_layout = QHBoxLayout(central_widget)
		central_widget.setLayout(central_layout)

		# --------------------------------------------------------------------------------
		# Construct and add a widget for menu item button area (containing all menu item buttons)
		# --------------------------------------------------------------------------------
		sxmenu_item_btn_area_widget = SXMenuItemBtnAreaWidget(self.sxconst_set, self.sxcmd_category_list, self.sxinfo, central_widget)
		central_layout.addWidget(sxmenu_item_btn_area_widget)

		# --------------------------------------------------------------------------------
		# Construct and add widgets for menu item widget area (containing all menu item widgets)
		# --------------------------------------------------------------------------------
		# Stacked layout for sx menu item widgets area
		self.sxmenu_item_widget_stacked_layout = QStackedLayout()
		central_layout.addLayout(self.sxmenu_item_widget_stacked_layout, stretch = 1)

		# Construct and add widgets for sx command categories
		for sxcmd_category in self.sxcmd_category_list:
			# Create SXCmdCategoryWidget for this command category
			sxcmd_category.widget = SXCmdCategoryWidget(self.sxconst_set, sxcmd_category)
			self.sxmenu_item_widget_stacked_layout.addWidget(sxcmd_category.widget)

		# Construct and add a widget for project constants settings
		self.sxconst_set.widget = SXConstSetWidget(self.sxconst_set, self.sxcmd_category_list)
		self.sxmenu_item_widget_stacked_layout.addWidget(self.sxconst_set.widget)

		# Construct and add a widget for GUI application information
		self.sxinfo.widget = SXInfoWidget()
		self.sxmenu_item_widget_stacked_layout.addWidget(self.sxinfo.widget)

		# --------------------------------------------------------------------------------
		# Set up event handler of all menu item buttons
		# --------------------------------------------------------------------------------
		for sxcmd_category in self.sxcmd_category_list:
			sxcmd_category.btn.clicked.connect(partial(self.handle_sxmenu_item_btn_event, sxcmd_category))
		self.sxconst_set.btn.clicked.connect(partial(self.handle_sxmenu_item_btn_event, self.sxconst_set))
		self.sxinfo.btn.clicked.connect(partial(self.handle_sxmenu_item_btn_event, self.sxinfo))

		# --------------------------------------------------------------------------------
		# Register project constant parameter settings upon initialization
		# --------------------------------------------------------------------------------
		self.sxconst_set.widget.register_const_set()

		# --------------------------------------------------------------------------------
		# Load the previously saved parameter setting of all sx commands
		# Override the registration of project constant parameter settings with the previously-saved one
		# --------------------------------------------------------------------------------
		for sxcmd_category in self.sxcmd_category_list:
			sxcmd_category.widget.load_previous_session()

		# --------------------------------------------------------------------------------
		# Start widget
		# --------------------------------------------------------------------------------
		start_widget = QtGui.QWidget()
		logo_container = QtGui.QWidget()
		layout_start_widget = QtGui.QHBoxLayout()
		layout_logo_container = QtGui.QVBoxLayout()
		logo_container.setStyleSheet('border-image: url("{0}sxgui_pictograph_info.png")'.format(get_image_directory()))
		logo_container.setFixedSize(100, 100)
		layout_start_widget.setContentsMargins(0, 0, 0, 20)

		layout_logo_container.addStretch(1)
		layout_logo_container.addWidget(logo_container)
		layout_start_widget.addLayout(layout_logo_container)
		layout_start_widget.addStretch(1)
		start_widget.setLayout(layout_start_widget)
		self.sxmenu_item_widget_stacked_layout.addWidget(start_widget)

		# --------------------------------------------------------------------------------
		# Display application information upon startup
		# --------------------------------------------------------------------------------
		self.sxmenu_item_widget_stacked_layout.setCurrentWidget(start_widget)

		# --------------------------------------------------------------------------------
		# Get focus to main window
		# --------------------------------------------------------------------------------
		self.setFocus()

	def construct_sxinfo(self):
		sxinfo = SXmenu_item(); sxinfo.name = "GUI Information"; sxinfo.label = "GUI Appliation Information"; sxinfo.short_info = "DUMMY STRING"

		# Store GUI application information as a class data member
		self.sxinfo = sxinfo

	def construct_sxconst_set(self):
		sxconst_set = SXconst_set(); sxconst_set.name = "sxc_project_settings"; sxconst_set.label = "Project Settings"; sxconst_set.short_info = "Set constant parameter values for this project. These constants will be used as default values of associated arugments and options in command settings. However, the setting here is not required to run commands."
		sxconst = SXconst(); sxconst.key = "protein"; sxconst.label = "Protein name"; sxconst.help = "a valid string for file names on your OS."; sxconst.register = "MY_PROTEIN"; sxconst.type = "string"; sxconst_set.list.append(sxconst); sxconst_set.dict[sxconst.key] = sxconst
		sxconst = SXconst(); sxconst.key = "apix"; sxconst.label = "Micrograph pixel size [A]"; sxconst.help = ""; sxconst.register = "1.0"; sxconst.type = "float"; sxconst_set.list.append(sxconst); sxconst_set.dict[sxconst.key] = sxconst
		sxconst = SXconst(); sxconst.key = "ctfwin"; sxconst.label = "CTF window size [pixels]"; sxconst.help = "it should be slightly larger than particle box size"; sxconst.register = "512"; sxconst.type = "int"; sxconst_set.list.append(sxconst); sxconst_set.dict[sxconst.key] = sxconst
		sxconst = SXconst(); sxconst.key = "box"; sxconst.label = "Particle box size [pixels]" ; sxconst.help = ""; sxconst.register = "-1"; sxconst.type = "int"; sxconst_set.list.append(sxconst); sxconst_set.dict[sxconst.key] = sxconst
		sxconst = SXconst(); sxconst.key = "radius"; sxconst.label = "Protein particle radius [pixels]"; sxconst.help = ""; sxconst.register = "-1"; sxconst.type = "int"; sxconst_set.list.append(sxconst); sxconst_set.dict[sxconst.key] = sxconst
		sxconst = SXconst(); sxconst.key = "sym"; sxconst.label = "Point-group symmetry"; sxconst.help = "e.g. c1, c4, d5"; sxconst.register = "c1"; sxconst.type = "string"; sxconst_set.list.append(sxconst); sxconst_set.dict[sxconst.key] = sxconst
		sxconst = SXconst(); sxconst.key = "mass"; sxconst.label = "Protein molecular mass [kDa]"; sxconst.help = ""; sxconst.register = "-1.0"; sxconst.type = "float"; sxconst_set.list.append(sxconst); sxconst_set.dict[sxconst.key] = sxconst
		sxconst = SXconst(); sxconst.key = "config"; sxconst.label = "Imaging configurations"; sxconst.help = "a free-style string for your record. please use it to describe the set of imaging configurations used in this project (e.g. types of microscope, detector, enegy filter, abbration corrector, phase plate, and etc."; sxconst.register = "MY_MICROSCOPE"; sxconst.type = "int"; sxconst_set.list.append(sxconst); sxconst_set.dict[sxconst.key] = sxconst

		# Store the project constant parameter set as a class data member
		self.sxconst_set = sxconst_set

	def construct_sxcmd_category_list(self):
		sxcmd_category_list = []
		sxcmd_list = []           # Used only within this function
		sxcmd_category_dict = {}  # Used only within this function

		# Actual configurations of all sx command categories and sx commands are inserted into the following section by wikiparser.py
		# as sxcmd_category_list and sxcmd_list
		# @@@@@ START_INSERTION @@@@@
		sxcmd_category = SXcmd_category(); sxcmd_category.name = "sxc_movie_micrograph"; sxcmd_category.label = "Movie Micrograph"; sxcmd_category.short_info = "movie frame alignemnt, and drift assessment"
		sxcmd_category_list.append(sxcmd_category)
		sxcmd_category = SXcmd_category(); sxcmd_category.name = "sxc_ctf"; sxcmd_category.label = "CTF"; sxcmd_category.short_info = "ctf estinatim, and ctf assessment"
		sxcmd_category_list.append(sxcmd_category)
		sxcmd_category = SXcmd_category(); sxcmd_category.name = "sxc_particle_stack"; sxcmd_category.label = "Particle Stack"; sxcmd_category.short_info = "particle picking, and particle windowing"
		sxcmd_category_list.append(sxcmd_category)
		sxcmd_category = SXcmd_category(); sxcmd_category.name = "sxc_2d_clustering"; sxcmd_category.label = "2D Clustering"; sxcmd_category.short_info = "2d clustering with isac, and post-processing"
		sxcmd_category_list.append(sxcmd_category)
		sxcmd_category = SXcmd_category(); sxcmd_category.name = "sxc_initial_3d_modeling"; sxcmd_category.label = "Initial 3D Modeling"; sxcmd_category.short_info = "initial 3d modeling with viper/rviper"
		sxcmd_category_list.append(sxcmd_category)
		sxcmd_category = SXcmd_category(); sxcmd_category.name = "sxc_3d_refinement"; sxcmd_category.label = "3D Refinement"; sxcmd_category.short_info = "3d refinement and post-processing"
		sxcmd_category_list.append(sxcmd_category)
		sxcmd_category = SXcmd_category(); sxcmd_category.name = "sxc_3d_clustering"; sxcmd_category.label = "3D Clustering"; sxcmd_category.short_info = "3d variability, and 3d clustering protocol I & II"
		sxcmd_category_list.append(sxcmd_category)
		sxcmd_category = SXcmd_category(); sxcmd_category.name = "sxc_localres"; sxcmd_category.label = "Local Resolution"; sxcmd_category.short_info = "local resolution, and local filter"
		sxcmd_category_list.append(sxcmd_category)
		sxcmd_category = SXcmd_category(); sxcmd_category.name = "sxc_utilities"; sxcmd_category.label = "Utilities"; sxcmd_category.short_info = "miscellaneous utlitity commands"
		sxcmd_category_list.append(sxcmd_category)

		sxcmd = SXcmd(); sxcmd.name = "sxunblur"; sxcmd.mode = ""; sxcmd.label = "Micrograph movie alignment"; sxcmd.short_info = "Align micrograph movies with Unblur."; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_movie_micrograph"; sxcmd.role = "sxr_pipe"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "unblur"; token.key_prefix = ""; token.label = "Path to unblur executable"; token.help = ""; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "exe"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "input_mrc_micrograph"; token.key_prefix = ""; token.label = "Input micrographs"; token.help = "You can use the wild cards (e.g. *) to specify a list of micrographs. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "mrc"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "output"; token.key_prefix = ""; token.label = "Output directory"; token.help = "The results will be written here. This directory will be created automatically and it must not exist previously. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "output"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "nr_threads"; token.key_prefix = "--"; token.label = "Number of threads"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "2"; token.restore = "2"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "nr_frames"; token.key_prefix = "--"; token.label = "Number of frames"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "3"; token.restore = "3"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "pixel_size"; token.key_prefix = "--"; token.label = "Pixel size [A]"; token.help = ""; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "apix"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "dose_filter"; token.key_prefix = "--"; token.label = "Apply dose filter"; token.help = "Requires additionally to set the correct exposure per frame, voltage, and pre exposure. "; token.group = "main"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "exposure_per_frame"; token.key_prefix = "--"; token.label = "Per frame exposure [e/A^2]"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "1.0"; token.restore = "1.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "voltage"; token.key_prefix = "--"; token.label = "Microscope voltage [kV]"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "300.0"; token.restore = "300.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "pre_exposure"; token.key_prefix = "--"; token.label = "Pre exposure [e/A^2]"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "0.0"; token.restore = "0.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "save_frames"; token.key_prefix = "--"; token.label = "Save aligned frames"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "expert_mode"; token.key_prefix = "--"; token.label = "Set expert mode"; token.help = "Requires the options initial shift, shift radius, b-factor, fourier_vertical, fourier_horizontal, shift threshold, iterations, restore noise, and verbosity. "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "shift_initial"; token.key_prefix = "--"; token.label = "Minimum shift for initial search [A]"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "2.0"; token.restore = "2.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "shift_radius"; token.key_prefix = "--"; token.label = "Outer radius shift limit [A]"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "200.0"; token.restore = "200.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "b_factor"; token.key_prefix = "--"; token.label = "Apply B-factor to images [A^2]"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "1500.0"; token.restore = "1500.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "fourier_vertical"; token.key_prefix = "--"; token.label = "Half-width of central vertical line of fourier mask"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "1"; token.restore = "1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "fourier_horizontal"; token.key_prefix = "--"; token.label = "Half-width of central horizontal line of fourier mask"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "1"; token.restore = "1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "shift_threshold"; token.key_prefix = "--"; token.label = "Termination shift threshold"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "0.1"; token.restore = "0.1"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "iterations"; token.key_prefix = "--"; token.label = "Maximum iterations"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "10"; token.restore = "10"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "restore_noise"; token.key_prefix = "--"; token.label = "Restore noise power"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "verbose"; token.key_prefix = "--"; token.label = "Verbose"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "filter_sum"; token.key_prefix = "--"; token.label = "Apply filter to micrographs"; token.help = "When selected, an extra set images will be created by applying a band-pass filter to the corrected sums. Useful for picking purposes. "; token.group = "main"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "lowpass"; token.key_prefix = "--"; token.label = "Highest frequency cut-off"; token.help = "Absolute frequency cut-off for the low-pass filter. "; token.group = "main"; token.is_required = False; token.default = "0.033"; token.restore = "0.033"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "highpass"; token.key_prefix = "--"; token.label = "Lowest frequency cut-off"; token.help = "Absolute frequency cut-off for the high-pass filter. "; token.group = "main"; token.is_required = False; token.default = "0.00033"; token.restore = "0.00033"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "remove_sum"; token.key_prefix = "--"; token.label = "Remove corrected sums"; token.help = "Discards the micrographs containing the corrected sums. "; token.group = "main"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxgui_unblur"; sxcmd.mode = ""; sxcmd.label = "Drift Assessment"; sxcmd.short_info = "Assess micrographs based on drift estimation produced by Unblur."; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_movie_micrograph"; sxcmd.role = "sxr_pipe"; sxcmd.is_submittable = False
		token = SXcmd_token(); token.key_base = "inputfile"; token.key_prefix = ""; token.label = "Shift files"; token.help = "A wild card * can be used to process multiple shift files. "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "any_file"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "input_shift_list_file"; token.key_prefix = ""; token.label = "Input shift list file"; token.help = "Extension of input shift list file must be '.txt'. If this is not provided, all files matched with the micrograph name pattern will be processed. "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "txt"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "e2display"; sxcmd.mode = ""; sxcmd.label = "Display Data"; sxcmd.short_info = "Displays images, volumes, or 1D plots."; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_movie_micrograph"; sxcmd.role = "sxr_util"; sxcmd.is_submittable = False
		token = SXcmd_token(); token.key_base = "input_data_list"; token.key_prefix = ""; token.label = "Input files"; token.help = "List of input images, volumes, plots. Wild cards (e.g *) can be used to select a list of files. Not recommended when the list is too large. "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "any_file_list"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "singleimage"; token.key_prefix = "--"; token.label = "Single image view"; token.help = "Displays a stack in a single image view: "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "fullrange"; token.key_prefix = "--"; token.label = "Skip auto-contrast"; token.help = "Disables the auto contrast for the display of particles stacks and 2D images. "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "verbose"; token.key_prefix = "--"; token.label = "Verbose"; token.help = "Accepted values 0-9. "; token.group = "advanced"; token.is_required = False; token.default = "0"; token.restore = "0"; token.type = "int"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxcter"; sxcmd.mode = ""; sxcmd.label = "CTF Estimation"; sxcmd.short_info = "Automated estimation of CTF parameters with error assessment."; sxcmd.mpi_support = True; sxcmd.mpi_add_flag = True; sxcmd.category = "sxc_ctf"; sxcmd.role = "sxr_pipe"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "input_image"; token.key_prefix = ""; token.label = "Input micrographs"; token.help = "For the multi-micrograph mode, specify micrograph list file name or file name pattern with a wild card (i.e. *). A particle stack can also be supplied with --stack_mode, but this mode is not supported by sxgui. Images of bdb format can not be used. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "any_micrograph"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "output_directory"; token.key_prefix = ""; token.label = "Output directory"; token.help = "The CTF parameters (partres file), rotationally averaged power spectra (rotinf), and micrograph thumbnails (thumb files) will be written here. This directory will be created automatically and it must not exist previously. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "output"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "wn"; token.key_prefix = "--"; token.label = "CTF window size [Pixels]"; token.help = "It should be slightly larger than particle box size. This is ignored in the stack mode. "; token.group = "main"; token.is_required = False; token.default = "512"; token.restore = "512"; token.type = "ctfwin"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "apix"; token.key_prefix = "--"; token.label = "Pixel size [A]"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "-1.0"; token.restore = "-1.0"; token.type = "apix"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "Cs"; token.key_prefix = "--"; token.label = "Microscope spherical aberration [mm]"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "2.0"; token.restore = "2.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "voltage"; token.key_prefix = "--"; token.label = "Microscope voltage [kV]"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "300.0"; token.restore = "300.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ac"; token.key_prefix = "--"; token.label = "Amplitude contrast [%]"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "10.0"; token.restore = "10.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "f_start"; token.key_prefix = "--"; token.label = "Starting frequency [1/A]"; token.help = "Lowest frequency to be considered in the CTF estimation. Determined automatically by default. "; token.group = "main"; token.is_required = False; token.default = "-1.0"; token.restore = "-1.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "f_stop"; token.key_prefix = "--"; token.label = "Stop frequency [1/A]"; token.help = "by Highest frequency to be considered in the CTF estimation. Determined automatically by default. "; token.group = "main"; token.is_required = False; token.default = "-1.0"; token.restore = "-1.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "kboot"; token.key_prefix = "--"; token.label = "Number of estimates per micrograph"; token.help = "Used for error assessment. "; token.group = "advanced"; token.is_required = False; token.default = "16"; token.restore = "16"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "overlap_x"; token.key_prefix = "--"; token.label = "X overlap [%]"; token.help = "Overlap between the windows in the x direction. "; token.group = "advanced"; token.is_required = False; token.default = "50"; token.restore = "50"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "overlap_y"; token.key_prefix = "--"; token.label = "Y overlap [%]"; token.help = "Overlap between the windows in the y direction. "; token.group = "advanced"; token.is_required = False; token.default = "50"; token.restore = "50"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "edge_x"; token.key_prefix = "--"; token.label = "Edge x [Pixels]"; token.help = "Defines the edge of the tiling area in the x direction. Normally it does not need to be modified. "; token.group = "advanced"; token.is_required = False; token.default = "0"; token.restore = "0"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "edge_y"; token.key_prefix = "--"; token.label = "Edge y [Pixels]"; token.help = "Defines the edge of the tiling area in the y direction. Normally it does not need to be modified. "; token.group = "advanced"; token.is_required = False; token.default = "0"; token.restore = "0"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "set_ctf_header"; token.key_prefix = "--"; token.label = "Export parameters to header"; token.help = "Exports the estimated CTF parameters to the image header. Not supported in the particle stack mode. "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "debug"; token.key_prefix = "--"; token.label = "Debug info"; token.help = "Prints debug information. "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxgui_cter"; sxcmd.mode = ""; sxcmd.label = "CTF Assessment"; sxcmd.short_info = "GUI tool to assess and sort micrographs based on their estimated CTF by sxcter."; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_ctf"; sxcmd.role = "sxr_pipe"; sxcmd.is_submittable = False
		token = SXcmd_token(); token.key_base = "cter_ctf_file"; token.key_prefix = ""; token.label = "CTF parameter file"; token.help = "This file is produced by sxcter and normally called partres.txt "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "txt"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "e2display"; sxcmd.mode = ""; sxcmd.label = "Display Data"; sxcmd.short_info = "Displays images, volumes, or 1D plots."; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_ctf"; sxcmd.role = "sxr_util"; sxcmd.is_submittable = False
		token = SXcmd_token(); token.key_base = "input_data_list"; token.key_prefix = ""; token.label = "Input files"; token.help = "List of input images, volumes, plots. Wild cards (e.g *) can be used to select a list of files. Not recommended when the list is too large. "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "any_file_list"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "singleimage"; token.key_prefix = "--"; token.label = "Single image view"; token.help = "Displays a stack in a single image view: "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "fullrange"; token.key_prefix = "--"; token.label = "Skip auto-contrast"; token.help = "Disables the auto contrast for the display of particles stacks and 2D images. "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "verbose"; token.key_prefix = "--"; token.label = "Verbose"; token.help = "Accepted values 0-9. "; token.group = "advanced"; token.is_required = False; token.default = "0"; token.restore = "0"; token.type = "int"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "e2boxer_old"; sxcmd.mode = ""; sxcmd.label = "Particle Picking"; sxcmd.short_info = "Manual and automatic particle picking."; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_particle_stack"; sxcmd.role = "sxr_pipe"; sxcmd.is_submittable = False
		token = SXcmd_token(); token.key_base = "input_micrograph_list"; token.key_prefix = ""; token.label = "Input micrographs"; token.help = "Wild cards (e.g. *) can be used to specify a list of micrographs. Not recommended if the number is very large. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "any_image_list"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "boxsize"; token.key_prefix = "--"; token.label = "Box size [Pixels]"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "-1"; token.restore = "-1"; token.type = "box"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "exclude_edges"; token.key_prefix = "--"; token.label = "Exclude edge images"; token.help = "Exclude particles extending outside the micrograph: "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "norm"; token.key_prefix = "--"; token.label = "Particle normalization"; token.help = "Normalization processor applied to the final particle images: "; token.group = "advanced"; token.is_required = False; token.default = "normalize.edgemean"; token.restore = "normalize.edgemean"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "invert"; token.key_prefix = "--"; token.label = "Invert contrast"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "verbose"; token.key_prefix = "--"; token.label = "Verbose"; token.help = "Accepted values 0-9. "; token.group = "advanced"; token.is_required = False; token.default = "0"; token.restore = "0"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "gauss_autoboxer"; token.key_prefix = "--"; token.label = "Sample autoboxed image"; token.help = "Gets the parameters used for autoboxing from the EMANDB. It requires the name of a micrograph used for automatic boxing previously. "; token.group = "advanced"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "string"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxwindow"; sxcmd.mode = ""; sxcmd.label = "Particle Extraction"; sxcmd.short_info = "Windows particles from micrographs using the particles coordinates."; sxcmd.mpi_support = True; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_particle_stack"; sxcmd.role = "sxr_pipe"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "input_micrograph_pattern"; token.key_prefix = ""; token.label = "Input micrograph path pattern"; token.help = "Specify path pattern of input micrographs with a wild card (*). Use the wild card to indicate the place of variable part of the file names (e.g. serial number, time stamp, and etc). The path pattern must be enclosed by single quotes (') or double quotes ('). (Note: sxgui.py automatically adds single quotes (')). The substring at the variable part must be same between the associated pair of input micrograph and coordinates file. bdb files can not be selected as input micrographs. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "any_image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "input_coordinates_pattern"; token.key_prefix = ""; token.label = "Input coordinates path pattern"; token.help = "Specify path pattern of input coordinates files with a wild card (*). Use the wild card to indicate the place of variable part of the file names (e.g. serial number, time stamp, and etc). The path pattern must be enclosed by single quotes (') or double quotes ('). (Note: sxgui.py automatically adds single quotes (')). The substring at the variable part must be same between the associated pair of input micrograph and coordinates file. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "parameters"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "output_directory"; token.key_prefix = ""; token.label = "Output directory"; token.help = "The results will be written here. This directory will be created automatically and it must not exist previously. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "output"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "selection_list"; token.key_prefix = "--"; token.label = "Micrograph selecting list"; token.help = "Specify a name of micrograph selection list text file for Selected Micrographs Mode. The file extension must be '.txt'. Alternatively, the file name of a single micrograph can be specified for Single Micrograph Mode. "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "any_micrograph"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "coordinates_format"; token.key_prefix = "--"; token.label = "Coordinate file format"; token.help = "Allowed values are 'sparx', 'eman1', 'eman2', or 'spider'. The sparx, eman2, and spider formats use the particle center as coordinates. The eman1 format uses the lower left corner of the box as coordinates. "; token.group = "main"; token.is_required = False; token.default = "eman1"; token.restore = "eman1"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "box_size"; token.key_prefix = "--"; token.label = "Particle box size [Pixels]"; token.help = "The x and y dimensions of square area to be windowed. The box size after resampling is assumed when resample_ratio < 1.0. "; token.group = "main"; token.is_required = False; token.default = "256"; token.restore = "256"; token.type = "box"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "skip_invert"; token.key_prefix = "--"; token.label = "Skip invert image contrast"; token.help = "Use this option for negative staining data. By default, the image contrast is inverted for cryo data. "; token.group = "main"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "import_ctf"; token.key_prefix = "--"; token.label = "CTER CTF parameter file"; token.help = "The file produced by sxcter and normally called partres.txt. "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "parameters"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "limit_ctf"; token.key_prefix = "--"; token.label = "Use CTF limit filter"; token.help = "Frequencies whose oscillation can not be properly modeled at the current pixel size are discarded in the images with the appropriate low-pass filter. This option requires --import_ctf. "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "astigmatism_error"; token.key_prefix = "--"; token.label = "Astigmatism error limit [Degrees]"; token.help = "Set astigmatism to zero for all micrographs where the angular error computed by sxcter is larger than the desired value. "; token.group = "advanced"; token.is_required = False; token.default = "360.0"; token.restore = "360.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "resample_ratio"; token.key_prefix = "--"; token.label = "Ratio between new and original pixel size"; token.help = "Use values between 0.0 and 1.0. The new pixel size is automatically recalculated when resample_ratio < 1.0 is used. "; token.group = "advanced"; token.is_required = False; token.default = "1.0"; token.restore = "1.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "check_consistency"; token.key_prefix = "--"; token.label = "Check consistency of inputs"; token.help = "Create the text file containing the list of inconsistent Micrograph ID entries (i.e. inconsist_mic_list_file.txt) "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "e2bdb"; sxcmd.mode = "makevstack"; sxcmd.label = "Particle Stack"; sxcmd.short_info = "Output virtual image stack. Make a 'virtual' BDB image stack with the specified name from one or more other stacks. "; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_particle_stack"; sxcmd.role = "sxr_pipe"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "makevstack"; token.key_prefix = "--"; token.label = "Output virtual image stack"; token.help = "Make a 'virtual' BDB image stack with the specified name from one or more other stacks. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "output"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "input_bdb_stack_pattern"; token.key_prefix = ""; token.label = "Input BDB image stack pattern"; token.help = "Specify file path pattern of stack subsets created in particle extraction using a wild card /'*/' (e.g. /'//sxwindow_output_dir//*/'). The stack subsets are located in the sxwindow output directory."; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "any_directory"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "e2display"; sxcmd.mode = ""; sxcmd.label = "Display Data"; sxcmd.short_info = "Displays images, volumes, or 1D plots."; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_particle_stack"; sxcmd.role = "sxr_util"; sxcmd.is_submittable = False
		token = SXcmd_token(); token.key_base = "input_data_list"; token.key_prefix = ""; token.label = "Input files"; token.help = "List of input images, volumes, plots. Wild cards (e.g *) can be used to select a list of files. Not recommended when the list is too large. "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "any_file_list"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "singleimage"; token.key_prefix = "--"; token.label = "Single image view"; token.help = "Displays a stack in a single image view: "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "fullrange"; token.key_prefix = "--"; token.label = "Skip auto-contrast"; token.help = "Disables the auto contrast for the display of particles stacks and 2D images. "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "verbose"; token.key_prefix = "--"; token.label = "Verbose"; token.help = "Accepted values 0-9. "; token.group = "advanced"; token.is_required = False; token.default = "0"; token.restore = "0"; token.type = "int"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxisac"; sxcmd.mode = ""; sxcmd.label = "ISAC - 2D Clustering"; sxcmd.short_info = "Iterative Stable Alignment and Clustering (ISAC) of a 2D image stack."; sxcmd.mpi_support = True; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_2d_clustering"; sxcmd.role = "sxr_pipe"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "stack_file"; token.key_prefix = ""; token.label = "Input image stack"; token.help = "The images must to be square (''nx''=''ny''). The stack can be either in bdb or hdf format. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "output_directory"; token.key_prefix = ""; token.label = "Output directory"; token.help = "The directory will be automatically created and the results will be written here. If the directory already exists, results will be written there, possibly overwriting previous runs. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "output"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "radius"; token.key_prefix = "--"; token.label = "Particle radius [Pixels]"; token.help = "Radius of the particle (pixels) "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "radius"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "img_per_grp"; token.key_prefix = "--"; token.label = "Images per class"; token.help = "Number of images per class in an ideal situation. In practice, it defines the maximum size of the classes. "; token.group = "main"; token.is_required = False; token.default = "100"; token.restore = "100"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "CTF"; token.key_prefix = "--"; token.label = "Phase-flip"; token.help = "If set, the data will be phase-flipped using CTF information included in the image headers. "; token.group = "main"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "xr"; token.key_prefix = "--"; token.label = "X search range [Pixels]"; token.help = "The translational search range in the x direction. Set by the program by default. "; token.group = "main"; token.is_required = False; token.default = "1"; token.restore = "1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "thld_err"; token.key_prefix = "--"; token.label = "Pixel error threshold [Pixels]"; token.help = "Used for checking stability. It is defined as the root mean square of distances between corresponding pixels from set of found transformations and theirs average transformation, depends linearly on square of radius (parameter ou). units - pixels. "; token.group = "main"; token.is_required = False; token.default = "0.7"; token.restore = "0.7"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "stop_after_candidates"; token.key_prefix = "--"; token.label = "Stop after candidates step"; token.help = "The run stops after the 'candidate_class_averages' section is created. "; token.group = "main"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "restart_section"; token.key_prefix = "--"; token.label = "Restart section"; token.help = "Each iteration contains three sections: 'restart', 'candidate_class_averages', and 'reproducible_class_averages'. To restart, for example, from generation 4 - section 'candidate_class_averages', then set to: '--restart_section=candidate_class_averages,4'. \\ The option requires no white space before or after the comma. By default, the execution restarts from where it stopped. A default restart requires also to provide the name of the directory as argument. Alternatively, the '--use_latest_master_directory' option can be used. "; token.group = "main"; token.is_required = False; token.default = "' '"; token.restore = "' '"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "target_radius"; token.key_prefix = "--"; token.label = "Target particle radius [Pixels]"; token.help = "Particle radius used by isac to process the data. The images will be resized to fit this radius "; token.group = "main"; token.is_required = False; token.default = "29"; token.restore = "29"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "target_nx"; token.key_prefix = "--"; token.label = "Target particle image size [Pixels]"; token.help = "Image size used by isac to process the data. The images will be resized according to target particle radius and then cut/padded to achieve the target image size. When xr > 0, the final image size for isac processing is 'target_nx + xr - 1'  "; token.group = "main"; token.is_required = False; token.default = "76"; token.restore = "76"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ir"; token.key_prefix = "--"; token.label = "Inner ring [Pixels]"; token.help = "Inner of the resampling to polar coordinates. "; token.group = "advanced"; token.is_required = False; token.default = "1"; token.restore = "1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "rs"; token.key_prefix = "--"; token.label = "Ring step [Pixels]"; token.help = "Step of the resampling to polar coordinates. "; token.group = "advanced"; token.is_required = False; token.default = "1"; token.restore = "1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "yr"; token.key_prefix = "--"; token.label = "Y search range [Pixels]"; token.help = "The translational search range in the y direction. Set as xr by default. "; token.group = "advanced"; token.is_required = False; token.default = "-1"; token.restore = "-1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ts"; token.key_prefix = "--"; token.label = "Search step [Pixels]"; token.help = "Translational search step. "; token.group = "advanced"; token.is_required = False; token.default = "1.0"; token.restore = "1.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "maxit"; token.key_prefix = "--"; token.label = "Reference-free alignment iterations"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "30"; token.restore = "30"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "center_method"; token.key_prefix = "--"; token.label = "Centering method"; token.help = "Method to center global 2D average during the initial prealignment of the data (0 : no centering; -1 : average shift method; please see center_2D in utilities.py for methods 1-7). "; token.group = "advanced"; token.is_required = False; token.default = "-1"; token.restore = "-1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "dst"; token.key_prefix = "--"; token.label = "Discrete angle used for within-group alignment"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "90.0"; token.restore = "90.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "FL"; token.key_prefix = "--"; token.label = "Lowest filter frequency [1/Pixel]"; token.help = "Lowest frequency used for the tangent filter. "; token.group = "advanced"; token.is_required = False; token.default = "0.2"; token.restore = "0.2"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "FH"; token.key_prefix = "--"; token.label = "Highest filter frequency [1/Pixel]"; token.help = "Highest frequency used for the tangent filter. "; token.group = "advanced"; token.is_required = False; token.default = "0.3"; token.restore = "0.3"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "FF"; token.key_prefix = "--"; token.label = "Tangent filter fall-off"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "0.2"; token.restore = "0.2"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "init_iter"; token.key_prefix = "--"; token.label = "SAC initialization iterations"; token.help = "Number of ab-initio-within-cluster alignment runs used for stability evaluation during SAC initialization. "; token.group = "advanced"; token.is_required = False; token.default = "3"; token.restore = "3"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "main_iter"; token.key_prefix = "--"; token.label = "SAC main iterations"; token.help = "Number of ab-initio-within-cluster alignment runs used for stability evaluation during the main SAC. "; token.group = "advanced"; token.is_required = False; token.default = "3"; token.restore = "3"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "iter_reali"; token.key_prefix = "--"; token.label = "SAC stability check interval"; token.help = "Defines every how many iterations the SAC stability checking is performed. "; token.group = "advanced"; token.is_required = False; token.default = "1"; token.restore = "1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "match_first"; token.key_prefix = "--"; token.label = "Initial phase 2-way matching iterations"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "1"; token.restore = "1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "max_round"; token.key_prefix = "--"; token.label = "Maximum candidate generation rounds"; token.help = "Maximum number of rounds to generate the candidate class averages in the first phase. "; token.group = "advanced"; token.is_required = False; token.default = "20"; token.restore = "20"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "match_second"; token.key_prefix = "--"; token.label = "Second phae 2- or 3-way matching iterations"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "5"; token.restore = "5"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "stab_ali"; token.key_prefix = "--"; token.label = "Number of alignments for stability check"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "5"; token.restore = "5"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "indep_run"; token.key_prefix = "--"; token.label = "Independent runs for reproducibility tests"; token.help = "By default, perform full ISAC to 4-way matching. Value indep_run=2 will restrict ISAC to 2-way matching and 3 to 3-way matching.  Note the number of used MPI processes requested in mpirun must be a multiplicity of indep_run. "; token.group = "advanced"; token.is_required = False; token.default = "4"; token.restore = "4"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "thld_grp"; token.key_prefix = "--"; token.label = "Minimum size of reproducible class"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "10"; token.restore = "10"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "n_generations"; token.key_prefix = "--"; token.label = "Maximum generations"; token.help = "The program stops when reaching this total number of generations: "; token.group = "advanced"; token.is_required = False; token.default = "10"; token.restore = "10"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "rand_seed"; token.key_prefix = "--"; token.label = "Seed"; token.help = "Useful for testing purposes. By default, isac sets a random seed number. "; token.group = "advanced"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "new"; token.key_prefix = "--"; token.label = "Use new code"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "debug"; token.key_prefix = "--"; token.label = "Debug info"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "use_latest_master_directory"; token.key_prefix = "--"; token.label = "Use latest master directory"; token.help = "When active, the program looks for the latest directory that starts with the word 'master'. "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "skip_prealignment"; token.key_prefix = "--"; token.label = "Skip pre-alignment"; token.help = "Useful if images are already centered. The 2dalignment directory will still be generated but the parameters will be zero. "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxprocess"; sxcmd.mode = "isacselect"; sxcmd.label = "Get ISAC Particles"; sxcmd.short_info = "Create ISAC particle ID list. Retrieve original image numbers listed in ISAC output average stack. see ISAC documentation for details. "; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_2d_clustering"; sxcmd.role = "sxr_pipe"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "isacselect"; token.key_prefix = "--"; token.label = "Create ISAC particle ID list"; token.help = "Retrieve original image numbers listed in ISAC output average stack. see ISAC documentation for details. "; token.group = "main"; token.is_required = True; token.default = True; token.restore = True; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "class_file_name"; token.key_prefix = ""; token.label = "ISAC class file name"; token.help = "File name of the class averages. It is located in the ISAC output directory."; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "output_list"; token.key_prefix = ""; token.label = "Output ISAC particle ID list"; token.help = "Output text file containing retrieved member particle IDs of all ISAC classes."; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "output"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "e2bdb"; sxcmd.mode = "makevstack"; sxcmd.label = "Create Stack Subset"; sxcmd.short_info = "Output virtual image stack. Make a 'virtual' BDB image stack with the specified name from one or more other stacks. "; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_2d_clustering"; sxcmd.role = "sxr_pipe"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "makevstack"; token.key_prefix = "--"; token.label = "Output virtual image stack"; token.help = "Make a 'virtual' BDB image stack with the specified name from one or more other stacks. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "output"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "input_bdb_stack_file"; token.key_prefix = ""; token.label = "Input BDB image stack"; token.help = "File name of the class averages. It is located in the ISAC output directory."; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "bdb"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "list"; token.key_prefix = "--"; token.label = "File containing selection list of images"; token.help = "Input ASCII file containing a list of selected image names to creates a new virtual BDB image stack from an existed virtual stack. "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "parameters"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxisac_post_processing"; sxcmd.mode = ""; sxcmd.label = "ISAC - Postprocessing"; sxcmd.short_info = "Postprocess of the 2D clustering result produced by ISAC."; sxcmd.mpi_support = True; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_2d_clustering"; sxcmd.role = "sxr_pipe"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "stack_file"; token.key_prefix = ""; token.label = "Original image stack"; token.help = "Particles required to create the full-sized class averages. The images must be square (''nx''=''ny'') and the stack must be bdb format "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "bdb"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "isac_directory"; token.key_prefix = ""; token.label = "Isac output directory"; token.help = "Name of the directory where isac was run previously. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "directory"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "class_file_name_no_dir_info"; token.key_prefix = ""; token.label = "Isac class file name"; token.help = "Name of the file (no directory info) that contains the class averages. It is located in the Isac directory. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "radius"; token.key_prefix = "--"; token.label = "Particle radius [Pixels]"; token.help = "There is no default radius. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "radius"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "CTF"; token.key_prefix = "--"; token.label = "Phase-flip"; token.help = "If set, the data will be phase-flipped using CTF information included in the image headers. "; token.group = "main"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "single_stack_output"; token.key_prefix = "--"; token.label = "Single stack output"; token.help = "If set, only one stack for all classes will be generated. Otherwise, which is the default, one file per class will be generated. "; token.group = "main"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "e2display"; sxcmd.mode = ""; sxcmd.label = "Display Data"; sxcmd.short_info = "Displays images, volumes, or 1D plots."; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_2d_clustering"; sxcmd.role = "sxr_util"; sxcmd.is_submittable = False
		token = SXcmd_token(); token.key_base = "input_data_list"; token.key_prefix = ""; token.label = "Input files"; token.help = "List of input images, volumes, plots. Wild cards (e.g *) can be used to select a list of files. Not recommended when the list is too large. "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "any_file_list"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "singleimage"; token.key_prefix = "--"; token.label = "Single image view"; token.help = "Displays a stack in a single image view: "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "fullrange"; token.key_prefix = "--"; token.label = "Skip auto-contrast"; token.help = "Disables the auto contrast for the display of particles stacks and 2D images. "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "verbose"; token.key_prefix = "--"; token.label = "Verbose"; token.help = "Accepted values 0-9. "; token.group = "advanced"; token.is_required = False; token.default = "0"; token.restore = "0"; token.type = "int"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxrviper"; sxcmd.mode = ""; sxcmd.label = "Initial 3D Model - RVIPER"; sxcmd.short_info = "Reproducible ''ab initio'' 3D structure determination. The program is designed to determine a validated initial intermediate resolution structure using a small set (<100?) of class averages produced by ISAC."; sxcmd.mpi_support = True; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_initial_3d_modeling"; sxcmd.role = "sxr_pipe"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "stack"; token.key_prefix = ""; token.label = "Input images stack"; token.help = "A small set (<100) of class averages produced by ISAC. The images must be square and the stack must be hdf format . "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "hdf"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "output_directory"; token.key_prefix = ""; token.label = "Output directory"; token.help = "The directory will be automatically created and the results will be written here. If the directory already exists, results will be written there, possibly overwriting previous runs. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "output"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "radius"; token.key_prefix = "--"; token.label = "Target particle radius [Pixels]"; token.help = "Use the same value as in ISAC. It has to be less than half the box size. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "sym"; token.key_prefix = "--"; token.label = "Point-group symmetry"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "c1"; token.restore = "c1"; token.type = "sym"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "iteration_start"; token.key_prefix = "--"; token.label = "Starting iteration"; token.help = "Iteration from which to start the program. 0 means go to the most recent one. "; token.group = "main"; token.is_required = False; token.default = "0"; token.restore = "0"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ir"; token.key_prefix = "--"; token.label = "Inner rotational search radius [Pixels]"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "1"; token.restore = "1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "rs"; token.key_prefix = "--"; token.label = "Ring step size"; token.help = "Step between rings used for the rotational search. "; token.group = "advanced"; token.is_required = False; token.default = "1"; token.restore = "1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "xr"; token.key_prefix = "--"; token.label = "X search range [Pixels]"; token.help = "The translational search range in the x direction will take place in a +/xr range. "; token.group = "advanced"; token.is_required = False; token.default = "'0'"; token.restore = "'0'"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "yr"; token.key_prefix = "--"; token.label = "Y search range [Pixels]"; token.help = "The translational search range in the y direction. If omitted it will be xr. "; token.group = "advanced"; token.is_required = False; token.default = "'0'"; token.restore = "'0'"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ts"; token.key_prefix = "--"; token.label = "Translational search step [Pixels]"; token.help = "The search will be performed in -xr, -xr+ts, 0, xr-ts, xr, can be fractional. "; token.group = "advanced"; token.is_required = False; token.default = "'1.0'"; token.restore = "'1.0'"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "delta"; token.key_prefix = "--"; token.label = "Projection angular step [Degree]"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "'2.0'"; token.restore = "'2.0'"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "center"; token.key_prefix = "--"; token.label = "Center 3D template"; token.help = "0: no centering; 1: center of gravity. "; token.group = "advanced"; token.is_required = False; token.default = "-1.0"; token.restore = "-1.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "maxit1"; token.key_prefix = "--"; token.label = "Maximum iterations - GA step"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "400"; token.restore = "400"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "maxit2"; token.key_prefix = "--"; token.label = "Maximum iterations - Finish step"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "50"; token.restore = "50"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "L2threshold"; token.key_prefix = "--"; token.label = "GA stop threshold"; token.help = "Defines the maximum relative dispersion of volumes' L2 norms. "; token.group = "advanced"; token.is_required = False; token.default = "0.03"; token.restore = "0.03"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "doga"; token.key_prefix = "--"; token.label = "Threshold to start GA"; token.help = "Do GA when the fraction of orientation that changes less than 1.0 degrees is at least this fraction. "; token.group = "advanced"; token.is_required = False; token.default = "0.1"; token.restore = "0.1"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "n_shc_runs"; token.key_prefix = "--"; token.label = "GA population size"; token.help = "This defines the number of quasi-independent volumes generated. (same as '--nruns' parameter from sxviper.py). "; token.group = "advanced"; token.is_required = False; token.default = "4"; token.restore = "4"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "n_rv_runs"; token.key_prefix = "--"; token.label = "Number of iterations"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "10"; token.restore = "10"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "n_v_runs"; token.key_prefix = "--"; token.label = "Number of viper runs per cycle"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "3"; token.restore = "3"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "outlier_percentile"; token.key_prefix = "--"; token.label = "Percentile for outlier"; token.help = "Threshold above which images are considered outliers and removed if 'percentile' is used as outlier selection method. "; token.group = "advanced"; token.is_required = False; token.default = "95.0"; token.restore = "95.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ref_a"; token.key_prefix = "--"; token.label = "Projection generation method"; token.help = "Method for generating the quasi-uniformly distributed projection directions. S- Saff algorithm, or P - Penczek 1994 algorithm. "; token.group = "advanced"; token.is_required = False; token.default = "S"; token.restore = "S"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "npad"; token.key_prefix = "--"; token.label = "Image padding factor"; token.help = "The images are padded to achieve the original size times this option. "; token.group = "advanced"; token.is_required = False; token.default = "2"; token.restore = "2"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "fl"; token.key_prefix = "--"; token.label = "Low-pass filter frequency [1/Pixel]"; token.help = "Using a hyperbolic tangent low-pass filter. "; token.group = "advanced"; token.is_required = False; token.default = "0.25"; token.restore = "0.25"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "aa"; token.key_prefix = "--"; token.label = "Low-pass filter fall-off"; token.help = "Fall-off of for the hyperbolic tangent low-pass filter. "; token.group = "advanced"; token.is_required = False; token.default = "0.1"; token.restore = "0.1"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "pwreference"; token.key_prefix = "--"; token.label = "Power spectrum reference"; token.help = "Text file containing a 1D reference power spectrum. "; token.group = "advanced"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "parameters"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "mask3D"; token.key_prefix = "--"; token.label = "3D mask"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "sphere"; token.restore = "sphere"; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "moon_elimination"; token.key_prefix = "--"; token.label = "Eliminate disconnected regions"; token.help = "Used to removed disconnected pieces from the model. It requires as argument a comma separated string with the mass in KDa and the pixel size. "; token.group = "advanced"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "criterion_name"; token.key_prefix = "--"; token.label = "Stable projection criterion"; token.help = "Used to decide if the volumes have a set of stable projections. Valid options are - '80th percentile',  or 'fastest increase in the last quartile'. "; token.group = "advanced"; token.is_required = False; token.default = "80th percentile"; token.restore = "80th percentile"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "outlier_index_threshold_method"; token.key_prefix = "--"; token.label = "Outlier selection method"; token.help = "Used to decide which images to keep. Valid options are -  'discontinuity_in_derivative', 'percentile', or 'angle_measure'. "; token.group = "advanced"; token.is_required = False; token.default = "discontinuity_in_derivative"; token.restore = "discontinuity_in_derivative"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "angle_threshold"; token.key_prefix = "--"; token.label = "Angle threshold"; token.help = "Threshold used to remove projections if 'angle_measure' is used to decide the outliers. "; token.group = "advanced"; token.is_required = False; token.default = "30"; token.restore = "30"; token.type = "int"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxprocess"; sxcmd.mode = "changesize"; sxcmd.label = "Resample VIPER Model"; sxcmd.short_info = "Resample 2D or 3D images. resample (decimate or interpolate up) images (2D or 3D) in a stack to change the pixel size. "; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_initial_3d_modeling"; sxcmd.role = "sxr_pipe"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "changesize"; token.key_prefix = "--"; token.label = "Resample 2D or 3D images"; token.help = "resample (decimate or interpolate up) images (2D or 3D) in a stack to change the pixel size. "; token.group = "main"; token.is_required = True; token.default = True; token.restore = True; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "input_stack"; token.key_prefix = ""; token.label = "Input 2D/3D image stack"; token.help = "Input 2D/3D image stack."; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "output_stack"; token.key_prefix = ""; token.label = "Output 2D/3D image stack"; token.help = "Resampled (decimated or interpolated up) 2D/3D image stack."; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "output"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ratio"; token.key_prefix = "--"; token.label = "Ratio of new to old image size"; token.help = "if < 1, the pixel size will increase and image size decrease. if > 1, the other way round. "; token.group = "main"; token.is_required = False; token.default = "1.0"; token.restore = "1.0"; token.type = "float"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "e2proc3d"; sxcmd.mode = "clip"; sxcmd.label = "Pad/Clip VIPER Model"; sxcmd.short_info = "Pad or clip. Make the output volume have this size by padding/clipping. 1, 3 or 6 arguments; '<x>[,<y>,<z>[,<xc>,<yc>,<zc>]]'. "; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_initial_3d_modeling"; sxcmd.role = "sxr_pipe"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "clip"; token.key_prefix = "--"; token.label = "Pad/Clip volume to specified size [Pixels]"; token.help = "Make the output volume have this size by padding/clipping. 1, 3 or 6 arguments; '<x>[,<y>,<z>[,<xc>,<yc>,<zc>]]'. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "input_volume"; token.key_prefix = ""; token.label = "Input volume"; token.help = "Input volume file name."; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "output_file"; token.key_prefix = ""; token.label = "Output clipped/padded volume"; token.help = "Output clipped/padded volume file name."; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "output"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxviper"; sxcmd.mode = ""; sxcmd.label = "Initial 3D Model - VIPER"; sxcmd.short_info = "''ab initio'' 3D structure determination using Validation of Individual Parameter Reproducibility (VIPER). Designed to determine a validated initial model using a small set of class averages produced by ISAC."; sxcmd.mpi_support = True; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_initial_3d_modeling"; sxcmd.role = "sxr_alt"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "stack"; token.key_prefix = ""; token.label = "Input images stack"; token.help = "A small set of Class averages produced by ISAC. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "directory"; token.key_prefix = ""; token.label = "Output directory"; token.help = "The directory will be automatically created and the results will be written here. If the directory already exists, results will be written there, possibly overwriting previous runs. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "output"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ir"; token.key_prefix = "--"; token.label = "Inner rotational search radius [Pixels]"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "1"; token.restore = "1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "radius"; token.key_prefix = "--"; token.label = "Target particle radius [Pixels]"; token.help = "Use the same value as in ISAC. It has to be less than half the box size. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "rs"; token.key_prefix = "--"; token.label = "Ring step size"; token.help = "Step between rings used for the rotational search. "; token.group = "advanced"; token.is_required = False; token.default = "1"; token.restore = "1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "xr"; token.key_prefix = "--"; token.label = "X search range [Pixels]"; token.help = "The translational search range in the x direction will take place in a +/xr range. "; token.group = "advanced"; token.is_required = False; token.default = "'0'"; token.restore = "'0'"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "yr"; token.key_prefix = "--"; token.label = "Y search range [Pixels]"; token.help = "The translational search range in the y direction. If omitted it will be xr. "; token.group = "advanced"; token.is_required = False; token.default = "'0'"; token.restore = "'0'"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "mask3D"; token.key_prefix = "--"; token.label = "3D mask"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "sphere"; token.restore = "sphere"; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "moon_elimination"; token.key_prefix = "--"; token.label = "Eliminate disconnected regions"; token.help = "Used to removed disconnected pieces from the model. It requires as argument a comma separated string with the mass in KDa and the pixel size. "; token.group = "advanced"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ts"; token.key_prefix = "--"; token.label = "Translational search step [Pixels]"; token.help = "The search will be performed in -xr, -xr+ts, 0, xr-ts, xr, can be fractional. "; token.group = "advanced"; token.is_required = False; token.default = "'1.0'"; token.restore = "'1.0'"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "delta"; token.key_prefix = "--"; token.label = "Projection angular step [Degree]"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "'2.0'"; token.restore = "'2.0'"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "center"; token.key_prefix = "--"; token.label = "Center 3D template"; token.help = "0: no centering; 1: center of gravity "; token.group = "advanced"; token.is_required = False; token.default = "-1.0"; token.restore = "-1.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "maxit1"; token.key_prefix = "--"; token.label = "Maximum iterations - GA step"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "400"; token.restore = "400"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "maxit2"; token.key_prefix = "--"; token.label = "Maximum iterations - Finish step"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "50"; token.restore = "50"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "L2threshold"; token.key_prefix = "--"; token.label = "GA stop threshold"; token.help = "Defines the maximum relative dispersion of volumes' L2 norms. "; token.group = "advanced"; token.is_required = False; token.default = "0.03"; token.restore = "0.03"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ref_a"; token.key_prefix = "--"; token.label = "Projection generation method"; token.help = "Method for generating the quasi-uniformly distributed projection directions. S- Saff algorithm, or P - Penczek 1994 algorithm. "; token.group = "advanced"; token.is_required = False; token.default = "S"; token.restore = "S"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "sym"; token.key_prefix = "--"; token.label = "Point-group symmetry"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "c1"; token.restore = "c1"; token.type = "sym"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "nruns"; token.key_prefix = "--"; token.label = "GA population size"; token.help = "This defines the number of quasi-independent volumes generated. "; token.group = "advanced"; token.is_required = False; token.default = "6"; token.restore = "6"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "doga"; token.key_prefix = "--"; token.label = "Threshold to start GA"; token.help = "Do GA when the fraction of orientation that changes less than 1.0 degrees is at least this fraction. "; token.group = "advanced"; token.is_required = False; token.default = "0.1"; token.restore = "0.1"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "fl"; token.key_prefix = "--"; token.label = "Low-pass filter frequency [1/Pixel]"; token.help = "Using a hyperbolic tangent low-pass filter. "; token.group = "advanced"; token.is_required = False; token.default = "0.25"; token.restore = "0.25"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "aa"; token.key_prefix = "--"; token.label = "Low-pass filter fall-off"; token.help = "Fall-off of for the hyperbolic tangent low-pass filter. "; token.group = "advanced"; token.is_required = False; token.default = "0.1"; token.restore = "0.1"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "pwreference"; token.key_prefix = "--"; token.label = "Power spectrum reference"; token.help = "Text file containing a 1D reference power spectrum. "; token.group = "advanced"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "parameters"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "debug"; token.key_prefix = "--"; token.label = "Verbose"; token.help = "Print debug info. "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxpdb2em"; sxcmd.mode = ""; sxcmd.label = "PDB File Conversion"; sxcmd.short_info = "Converts an atomic model into a simulated density map"; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_initial_3d_modeling"; sxcmd.role = "sxr_alt"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "input_pdb"; token.key_prefix = ""; token.label = "Input PDB file"; token.help = "Starting atomic coordinates: "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "pdb"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "output_hdf"; token.key_prefix = ""; token.label = "Output map"; token.help = ""; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "output"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "apix"; token.key_prefix = "--"; token.label = "Pixel size of the output map [A]"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "1.0"; token.restore = "1.0"; token.type = "apix"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "box"; token.key_prefix = "--"; token.label = "Output box size [Voxels]"; token.help = "If not given, the program will find the minimum box size fitting the structure. Be aware that this will most likely result in a rectangular box. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "box"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "het"; token.key_prefix = "--"; token.label = "Include hetero atoms"; token.help = "Otherwise the HETATM entries in the PDB file are ignored. "; token.group = "main"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "center"; token.key_prefix = "--"; token.label = "Center model at the origin"; token.help = ": Specifies whether the atomic model should be moved to the origin before generating density map.  Available options are: c - Use the geometrical center of atoms; a - Use the center of mass (recommended); x,y,z - Vector to be subtracted for the coordinates. Default: no centering, in which case (0,0,0) in the PDB space will map to the center of the EM volume. "; token.group = "main"; token.is_required = False; token.default = "n"; token.restore = "n"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "O"; token.key_prefix = "--"; token.label = "Apply additional rotation"; token.help = "This can be used to modify the orientation of the atomic model. "; token.group = "main"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "tr0"; token.key_prefix = "--"; token.label = "Rotational matrix file"; token.help = ": This file must contain the 3x4 transformation matrix to be applied to the PDB coordinates after centering. The translation vector (last column of the matrix) must be specified in Angstrom. "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "parameters"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "quiet"; token.key_prefix = "--"; token.label = "Silent mode"; token.help = "Does not print any information to the monitor. "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "e2display"; sxcmd.mode = ""; sxcmd.label = "Display Data"; sxcmd.short_info = "Displays images, volumes, or 1D plots."; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_initial_3d_modeling"; sxcmd.role = "sxr_util"; sxcmd.is_submittable = False
		token = SXcmd_token(); token.key_base = "input_data_list"; token.key_prefix = ""; token.label = "Input files"; token.help = "List of input images, volumes, plots. Wild cards (e.g *) can be used to select a list of files. Not recommended when the list is too large. "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "any_file_list"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "singleimage"; token.key_prefix = "--"; token.label = "Single image view"; token.help = "Displays a stack in a single image view: "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "fullrange"; token.key_prefix = "--"; token.label = "Skip auto-contrast"; token.help = "Disables the auto contrast for the display of particles stacks and 2D images. "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "verbose"; token.key_prefix = "--"; token.label = "Verbose"; token.help = "Accepted values 0-9. "; token.group = "advanced"; token.is_required = False; token.default = "0"; token.restore = "0"; token.type = "int"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxprocess"; sxcmd.mode = "adaptive_mask"; sxcmd.label = "Adaptive 3D Mask"; sxcmd.short_info = "Create soft-edged 3D mask from reference volume. "; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_initial_3d_modeling"; sxcmd.role = "sxr_util"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "adaptive_mask"; token.key_prefix = "--"; token.label = "Create soft-edged 3D mask from reference volume"; token.help = ""; token.group = "main"; token.is_required = True; token.default = True; token.restore = True; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "input_volume"; token.key_prefix = ""; token.label = "Input volume"; token.help = "Input reference volume"; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "output_mask3D"; token.key_prefix = ""; token.label = "Output mask"; token.help = "Output 3D mask"; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "output"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "nsigma"; token.key_prefix = "--"; token.label = "Density standard deviations threshold"; token.help = "Defines the threshold used to find the main volume within the data. All voxels with density <= mean + nsigma standard deviations will be included into the main volume. This option will not be used if the option threshold is larger than -9999.0. "; token.group = "main"; token.is_required = False; token.default = "1.0"; token.restore = "1.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "threshold"; token.key_prefix = "--"; token.label = "Binarization threshold"; token.help = "Below this value the data is assumed to not belong to the main volume. With the value lower than the default, the option will be ignored and the mask will be set according to nsigma. "; token.group = "main"; token.is_required = False; token.default = "-9999.0"; token.restore = "-9999.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ndilation"; token.key_prefix = "--"; token.label = "Mask extension cycles"; token.help = "The initial mask will be extended this number of cycles. To keep the size of the main volume, set this to kernel_size/2 "; token.group = "main"; token.is_required = False; token.default = "3"; token.restore = "3"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "kernel_size"; token.key_prefix = "--"; token.label = "Gaussian kernel size"; token.help = "Size of the gaussian kernel used to smooth the binary mask. "; token.group = "main"; token.is_required = False; token.default = "11"; token.restore = "11"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "gauss_standard_dev"; token.key_prefix = "--"; token.label = "Kernel standard deviation "; token.help = "Standard deviation used in the construction of the gaussian smoothing of the mask. "; token.group = "main"; token.is_required = False; token.default = "9"; token.restore = "9"; token.type = "int"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxprocess"; sxcmd.mode = "angular_distribution"; sxcmd.label = "Angular Distribution"; sxcmd.short_info = "Create angular distribution file. The file will contain a 3D representation of the given angular distribution, which can be viewed with UCFS Chimera "; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_initial_3d_modeling"; sxcmd.role = "sxr_util"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "angular_distribution"; token.key_prefix = "--"; token.label = "Create angular distribution file"; token.help = "The file will contain a 3D representation of the given angular distribution, which can be viewed with UCFS Chimera "; token.group = "main"; token.is_required = True; token.default = True; token.restore = True; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "inputfile"; token.key_prefix = ""; token.label = "Alignment Parameter file"; token.help = "Alignment Parameter file created by a previous 3D reconstruction step (e.g. sxmeridien.py)"; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "any_file"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "pixel_size"; token.key_prefix = "--"; token.label = "Pixel size [A]"; token.help = "Pixel size in A. "; token.group = "main"; token.is_required = False; token.default = "1.0"; token.restore = "1.0"; token.type = "apix"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "round_digit"; token.key_prefix = "--"; token.label = "Number precision"; token.help = "Decimal numbers will be rounded to this number of decimal points. "; token.group = "main"; token.is_required = False; token.default = "5"; token.restore = "5"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "box_size"; token.key_prefix = "--"; token.label = "Box size [Pixels]"; token.help = "Box size in pixel used for calculating the center of the particle. "; token.group = "main"; token.is_required = False; token.default = "500"; token.restore = "500"; token.type = "box"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "particle_radius"; token.key_prefix = "--"; token.label = "Particle radius [Pixels]"; token.help = "Used for the representation in Chimera. Defines where the cylinders representing the histogram must start. "; token.group = "main"; token.is_required = False; token.default = "175"; token.restore = "175"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "cylinder_width"; token.key_prefix = "--"; token.label = "Cylinder width"; token.help = "Used for the representation in Chimera. This will define the width of the cylinders representing the histogram."; token.group = "main"; token.is_required = False; token.default = "1"; token.restore = "1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "cylinder_length"; token.key_prefix = "--"; token.label = "Cylinder length"; token.help = "Used for the representation in Chimera. This will define the relative size of the cylinders representing the histogram. "; token.group = "main"; token.is_required = False; token.default = "10000"; token.restore = "10000"; token.type = "int"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxmeridien"; sxcmd.mode = ""; sxcmd.label = "3D Refinement"; sxcmd.short_info = "Performs 3D structure refinement."; sxcmd.mpi_support = True; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_3d_refinement"; sxcmd.role = "sxr_pipe"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "stack"; token.key_prefix = ""; token.label = "Input image stack"; token.help = ""; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "output_directory"; token.key_prefix = ""; token.label = "Output directory"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "current directory"; token.restore = "current directory"; token.type = "output"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "initial_volume"; token.key_prefix = ""; token.label = "Initial 3D reference"; token.help = ""; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "radius"; token.key_prefix = "--"; token.label = "Particle radius [Pixels]"; token.help = "Outer radius [in pixels] of particles < int(nx/2)-1 "; token.group = "main"; token.is_required = False; token.default = "-1"; token.restore = "-1"; token.type = "radius"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "xr"; token.key_prefix = "--"; token.label = "Search range [Pixels]"; token.help = "Range for translation search in both directions. Search is +/-xr. It can be fractional "; token.group = "advanced"; token.is_required = False; token.default = "5.0"; token.restore = "5.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ts"; token.key_prefix = "--"; token.label = "Search step size [Pixels]"; token.help = "Step size of translation search in both directions. Search is within a circle of radius xr on a grid with steps ts. It can be fractional. "; token.group = "advanced"; token.is_required = False; token.default = "2.0"; token.restore = "2.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "inires"; token.key_prefix = "--"; token.label = "Starting resolution [A]"; token.help = "Resolution used to start the refinement. "; token.group = "advanced"; token.is_required = False; token.default = "25.0"; token.restore = "25.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "mask3D"; token.key_prefix = "--"; token.label = "3D mask"; token.help = "Soft mask for the volume. If not given, a hard sphere of radius boxsize/2-1 will be used. "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "hardmask"; token.key_prefix = "--"; token.label = "Apply hard 2D mask"; token.help = "Apply hard mask to 2D data with the particle radius. "; token.group = "main"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "sym"; token.key_prefix = "--"; token.label = "Point-group symmetry"; token.help = "Symmetry of the refined structure. Acceptable values are: cn, dn, where n is multiplicity. "; token.group = "main"; token.is_required = False; token.default = "c1"; token.restore = "c1"; token.type = "sym"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "skip_prealignment"; token.key_prefix = "--"; token.label = "Skip 2D pre-alignment step"; token.help = "Use this option if images are already centered. "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "initialshifts"; token.key_prefix = "--"; token.label = "Read shifts from header"; token.help = "Start with the shift parameters stored in the image headers. Useful to jumpstart the procedure. "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "center_method"; token.key_prefix = "--"; token.label = "Centering method"; token.help = "Method for centering averages during initial 2D prealignment of data (0: no centering; -1: average shift method; For 1-7, see center_2D in utilities.py) "; token.group = "advanced"; token.is_required = False; token.default = "-1"; token.restore = "-1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "target_radius"; token.key_prefix = "--"; token.label = "Target particle radius"; token.help = "For 2D prealignment, images will be shrank/enlarged to this radius. "; token.group = "advanced"; token.is_required = False; token.default = "29"; token.restore = "29"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "delta"; token.key_prefix = "--"; token.label = "Initial angular sampling step"; token.help = "Initial angular sampling step. "; token.group = "advanced"; token.is_required = False; token.default = "15.0"; token.restore = "15.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "shake"; token.key_prefix = "--"; token.label = "Shake"; token.help = "Shake "; token.group = "advanced"; token.is_required = False; token.default = "0.5"; token.restore = "0.5"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "small_memory"; token.key_prefix = "--"; token.label = "Do not keep data in memory"; token.help = "Data will not be kept in memory if small_memory is true "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ref_a"; token.key_prefix = "--"; token.label = "Projection generation method"; token.help = "Method for generating the quasi-uniformly distributed projection directions. S- Saff algorithm, or P - Penczek 1994 algorithm. "; token.group = "advanced"; token.is_required = False; token.default = "S"; token.restore = "S"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ccfpercentage"; token.key_prefix = "--"; token.label = "Correlation peaks to be included [%]"; token.help = "Percentage of correlation peaks to be included. 0.0 corresponds to hard matching. "; token.group = "advanced"; token.is_required = False; token.default = "99.9"; token.restore = "99.9"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "nonorm"; token.key_prefix = "--"; token.label = "Skip image norm correction"; token.help = "Do not apply image norm correction. "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "do_final"; token.key_prefix = "--"; token.label = "Perform final reconstruction"; token.help = "Perform final reconstruction using orientation parameters from the specified iteration. By default, it uses iteration achieved best resolution "; token.group = "main"; token.is_required = False; token.default = "-1"; token.restore = "-1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "memory_per_node"; token.key_prefix = "--"; token.label = "Memory per node [GB]"; token.help = "User provided information about memory per node in GB (NOT per CPU). By default, it uses 2GB * (number of CPUs per node) "; token.group = "main"; token.is_required = False; token.default = "-1.0"; token.restore = "-1.0"; token.type = "float"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxprocess"; sxcmd.mode = "postprocess"; sxcmd.label = "3D Refinement Postprocess"; sxcmd.short_info = "Adjust power spectrum based on B-factor. B-factor is estimated from unfiltered odd/even 3D volumes, one volume, or ISAC 2D averages. "; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_3d_refinement"; sxcmd.role = "sxr_pipe"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "postprocess"; token.key_prefix = "--"; token.label = "Adjust power spectrum based on B-factor"; token.help = "B-factor is estimated from unfiltered odd/even 3D volumes, one volume, or ISAC 2D averages. "; token.group = "main"; token.is_required = True; token.default = True; token.restore = True; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "firstvolume"; token.key_prefix = ""; token.label = "First unfiltered half-volume "; token.help = "Generated by sxmeridien"; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "secondvolume"; token.key_prefix = ""; token.label = "Second unfiltered half-volume "; token.help = "Generated by sxmeridien"; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "mtf"; token.key_prefix = "--"; token.label = "MTF file"; token.help = "File contains the MTF (modulation transfer function) of the detector used. "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "parameters"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "fsc_adj"; token.key_prefix = "--"; token.label = "Apply FSC-based low-pass filter"; token.help = "Applies an FSC-based low-pass filter to the merged volume before the B-factor estimation. This is not applicable to 2D image. "; token.group = "main"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "B_enhance"; token.key_prefix = "--"; token.label = "B-factor enhancement"; token.help = "-1.0: B-factor is not applied; 0.0: program estimates B-factor from options. B_start (usually 10 Angstrom) to the resolution determined by FSC143; 128.0: program use the given value 128.0 [A^2] to enhance map. "; token.group = "main"; token.is_required = False; token.default = "0.0"; token.restore = "0.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "low_pass_filter"; token.key_prefix = "--"; token.label = "Low-pass filter frequency"; token.help = "0.0: low_pass filter to resolution; -1.0: no low_pass filter; 5.8: low_pass filter to 5.8 [A]; 0.2: low_pass filter to 0.2 (absolute frequency). "; token.group = "main"; token.is_required = False; token.default = "0.0"; token.restore = "0.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "aa"; token.key_prefix = "--"; token.label = "Low-pass filter fall-off"; token.help = "Low-pass filter fall-off. "; token.group = "main"; token.is_required = False; token.default = "0.1"; token.restore = "0.1"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "mask"; token.key_prefix = "--"; token.label = "User-provided mask"; token.help = "User-provided mask. "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "output"; token.key_prefix = "--"; token.label = "Output file"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "vol_postrefine_masked.hdf"; token.restore = "vol_postrefine_masked.hdf"; token.type = "output"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "pixel_size"; token.key_prefix = "--"; token.label = "Pixel size [A]"; token.help = "Pixel size in A. "; token.group = "main"; token.is_required = False; token.default = "1.0"; token.restore = "1.0"; token.type = "apix"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "B_start"; token.key_prefix = "--"; token.label = "B-factor lower limit [A]"; token.help = "Lower limit for B-factor estimation. "; token.group = "main"; token.is_required = False; token.default = "10.0"; token.restore = "10.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "B_stop"; token.key_prefix = "--"; token.label = "B-factor higher limit [A]"; token.help = "Higher limit for B-factor estimation. "; token.group = "main"; token.is_required = False; token.default = "0.0"; token.restore = "0.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "do_adaptive_mask"; token.key_prefix = "--"; token.label = "Apply adaptive mask"; token.help = "Program creates mask adaptively with given density threshold. "; token.group = "main"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "mask_threshold"; token.key_prefix = "--"; token.label = "Adaptive mask threshold"; token.help = "Density threshold for creating adaptive surface mask. "; token.group = "main"; token.is_required = False; token.default = "0.02"; token.restore = "0.02"; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "consine_edge"; token.key_prefix = "--"; token.label = "Cosine edge width [Pixels]"; token.help = "Width of cosine transition area for soft-edge masking. "; token.group = "main"; token.is_required = False; token.default = "6.0"; token.restore = "6.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "dilation"; token.key_prefix = "--"; token.label = "surface dilation size [Pixels]"; token.help = "Size of surface dilation or erosion. "; token.group = "main"; token.is_required = False; token.default = "3.0"; token.restore = "3.0"; token.type = "float"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "e2display"; sxcmd.mode = ""; sxcmd.label = "Display Data"; sxcmd.short_info = "Displays images, volumes, or 1D plots."; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_3d_refinement"; sxcmd.role = "sxr_util"; sxcmd.is_submittable = False
		token = SXcmd_token(); token.key_base = "input_data_list"; token.key_prefix = ""; token.label = "Input files"; token.help = "List of input images, volumes, plots. Wild cards (e.g *) can be used to select a list of files. Not recommended when the list is too large. "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "any_file_list"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "singleimage"; token.key_prefix = "--"; token.label = "Single image view"; token.help = "Displays a stack in a single image view: "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "fullrange"; token.key_prefix = "--"; token.label = "Skip auto-contrast"; token.help = "Disables the auto contrast for the display of particles stacks and 2D images. "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "verbose"; token.key_prefix = "--"; token.label = "Verbose"; token.help = "Accepted values 0-9. "; token.group = "advanced"; token.is_required = False; token.default = "0"; token.restore = "0"; token.type = "int"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxprocess"; sxcmd.mode = "adaptive_mask"; sxcmd.label = "Adaptive 3D Mask"; sxcmd.short_info = "Create soft-edged 3D mask from reference volume. "; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_3d_refinement"; sxcmd.role = "sxr_util"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "adaptive_mask"; token.key_prefix = "--"; token.label = "Create soft-edged 3D mask from reference volume"; token.help = ""; token.group = "main"; token.is_required = True; token.default = True; token.restore = True; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "input_volume"; token.key_prefix = ""; token.label = "Input volume"; token.help = "Input reference volume"; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "output_mask3D"; token.key_prefix = ""; token.label = "Output mask"; token.help = "Output 3D mask"; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "output"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "nsigma"; token.key_prefix = "--"; token.label = "Density standard deviations threshold"; token.help = "Defines the threshold used to find the main volume within the data. All voxels with density <= mean + nsigma standard deviations will be included into the main volume. This option will not be used if the option threshold is larger than -9999.0. "; token.group = "main"; token.is_required = False; token.default = "1.0"; token.restore = "1.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "threshold"; token.key_prefix = "--"; token.label = "Binarization threshold"; token.help = "Below this value the data is assumed to not belong to the main volume. With the value lower than the default, the option will be ignored and the mask will be set according to nsigma. "; token.group = "main"; token.is_required = False; token.default = "-9999.0"; token.restore = "-9999.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ndilation"; token.key_prefix = "--"; token.label = "Mask extension cycles"; token.help = "The initial mask will be extended this number of cycles. To keep the size of the main volume, set this to kernel_size/2 "; token.group = "main"; token.is_required = False; token.default = "3"; token.restore = "3"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "kernel_size"; token.key_prefix = "--"; token.label = "Gaussian kernel size"; token.help = "Size of the gaussian kernel used to smooth the binary mask. "; token.group = "main"; token.is_required = False; token.default = "11"; token.restore = "11"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "gauss_standard_dev"; token.key_prefix = "--"; token.label = "Kernel standard deviation "; token.help = "Standard deviation used in the construction of the gaussian smoothing of the mask. "; token.group = "main"; token.is_required = False; token.default = "9"; token.restore = "9"; token.type = "int"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxprocess"; sxcmd.mode = "angular_distribution"; sxcmd.label = "Angular Distribution"; sxcmd.short_info = "Create angular distribution file. The file will contain a 3D representation of the given angular distribution, which can be viewed with UCFS Chimera "; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_3d_refinement"; sxcmd.role = "sxr_util"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "angular_distribution"; token.key_prefix = "--"; token.label = "Create angular distribution file"; token.help = "The file will contain a 3D representation of the given angular distribution, which can be viewed with UCFS Chimera "; token.group = "main"; token.is_required = True; token.default = True; token.restore = True; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "inputfile"; token.key_prefix = ""; token.label = "Alignment Parameter file"; token.help = "Alignment Parameter file created by a previous 3D reconstruction step (e.g. sxmeridien.py)"; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "any_file"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "pixel_size"; token.key_prefix = "--"; token.label = "Pixel size [A]"; token.help = "Pixel size in A. "; token.group = "main"; token.is_required = False; token.default = "1.0"; token.restore = "1.0"; token.type = "apix"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "round_digit"; token.key_prefix = "--"; token.label = "Number precision"; token.help = "Decimal numbers will be rounded to this number of decimal points. "; token.group = "main"; token.is_required = False; token.default = "5"; token.restore = "5"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "box_size"; token.key_prefix = "--"; token.label = "Box size [Pixels]"; token.help = "Box size in pixel used for calculating the center of the particle. "; token.group = "main"; token.is_required = False; token.default = "500"; token.restore = "500"; token.type = "box"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "particle_radius"; token.key_prefix = "--"; token.label = "Particle radius [Pixels]"; token.help = "Used for the representation in Chimera. Defines where the cylinders representing the histogram must start. "; token.group = "main"; token.is_required = False; token.default = "175"; token.restore = "175"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "cylinder_width"; token.key_prefix = "--"; token.label = "Cylinder width"; token.help = "Used for the representation in Chimera. This will define the width of the cylinders representing the histogram."; token.group = "main"; token.is_required = False; token.default = "1"; token.restore = "1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "cylinder_length"; token.key_prefix = "--"; token.label = "Cylinder length"; token.help = "Used for the representation in Chimera. This will define the relative size of the cylinders representing the histogram. "; token.group = "main"; token.is_required = False; token.default = "10000"; token.restore = "10000"; token.type = "int"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sx3dvariability"; sxcmd.mode = "symmetrize"; sxcmd.label = "3D Variability Preprocess"; sxcmd.short_info = "Symmetrise input stack. Please skip this preparation step if the structure is asymmetrical (i.e. c1), since it is required only when the structure has internal symmetry. "; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_3d_clustering"; sxcmd.role = "sxr_pipe"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "symmetrize"; token.key_prefix = "--"; token.label = "Symmetrise input stack"; token.help = "Please skip this preparation step if the structure is asymmetrical (i.e. c1), since it is required only when the structure has internal symmetry. "; token.group = "main"; token.is_required = True; token.default = True; token.restore = True; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "prj_stack"; token.key_prefix = ""; token.label = "Input image stack"; token.help = "The images must containt the 3D orientation parameters in the header and optionally CTF information. The output image stack is bdb:sdata. Please use it as an input image stack of sx3dvariability."; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "sym"; token.key_prefix = "--"; token.label = "Point-group symmetry"; token.help = "Only required if the structure has symmetry higher than c1. Notice this step can be run with only one CPU and there is no MPI version for it. "; token.group = "main"; token.is_required = False; token.default = "c1"; token.restore = "c1"; token.type = "sym"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sx3dvariability"; sxcmd.mode = ""; sxcmd.label = "3D Variablity Estimation"; sxcmd.short_info = "Calculate 3D variability using a set of aligned projection images as an input."; sxcmd.mpi_support = True; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_3d_clustering"; sxcmd.role = "sxr_pipe"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "prj_stack"; token.key_prefix = ""; token.label = "Input image stack"; token.help = "The images must containt the 3D orientation parameters in the header and optionally CTF information. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ave2D"; token.key_prefix = "--"; token.label = "Output 2D averages"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "var2D"; token.key_prefix = "--"; token.label = "Output 2D variances"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ave3D"; token.key_prefix = "--"; token.label = "Output 3D average"; token.help = "3D reconstruction computed from projections averaged within respective angular neighbourhood. It should be used to assess the resolvability and possible artifacts of the variability map. "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "var3D"; token.key_prefix = "--"; token.label = "Output 3D variability"; token.help = "It creates a volume containing, for each voxel, a measure of the variability in the dataset. Careful, time consuming! "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "img_per_grp"; token.key_prefix = "--"; token.label = "Number of projections"; token.help = "Images from the angular neighbourhood that will be used to estimate 2D variance for each projection data. The larger the number the less noisy the estimate, but the lower the resolution. Usage of large number also results in rotational artifacts in variances that will be visible in 3D variability volume. "; token.group = "main"; token.is_required = False; token.default = "10"; token.restore = "10"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "no_norm"; token.key_prefix = "--"; token.label = "Skip normalization"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "radius"; token.key_prefix = "--"; token.label = "Radius for 3D variability [Pixels]"; token.help = "Normally, use the particle radius. "; token.group = "advanced"; token.is_required = False; token.default = "-1"; token.restore = "-1"; token.type = "radius"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "npad"; token.key_prefix = "--"; token.label = "Image padding factor"; token.help = "The images are padded to achieve the original size times this option. "; token.group = "advanced"; token.is_required = False; token.default = "2"; token.restore = "2"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "sym"; token.key_prefix = "--"; token.label = "Point-group symmetry"; token.help = "Only required if the structure has symmetry higher than c1. Notice this step can be run with only one CPU and there is no MPI version for it. "; token.group = "main"; token.is_required = False; token.default = "c1"; token.restore = "c1"; token.type = "sym"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "fl"; token.key_prefix = "--"; token.label = "Low-pass filter frequency [1/Pixel]"; token.help = "Filter to be applied to the images prior to variability calculation. By default, no filtering. "; token.group = "main"; token.is_required = False; token.default = "0.0"; token.restore = "0.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "aa"; token.key_prefix = "--"; token.label = "Low-pass filter fall-off"; token.help = "Filter to be applied to the images prior to variability calculation. By default, no filtering. "; token.group = "main"; token.is_required = False; token.default = "0.0"; token.restore = "0.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "CTF"; token.key_prefix = "--"; token.label = "Use CTF"; token.help = "If set to true, the a CTF correction will be applied using the parameters found in the image headers. "; token.group = "main"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "VERBOSE"; token.key_prefix = "--"; token.label = "Verbose"; token.help = "Useful for debugging. "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "VAR"; token.key_prefix = "--"; token.label = "Stack on input consists of 2D variances"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "SND"; token.key_prefix = "--"; token.label = "Squared normalized differences"; token.help = "If set to true, it computes the squared normalised differences. "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxsort3d"; sxcmd.mode = ""; sxcmd.label = "3D Clustering Protocol I"; sxcmd.short_info = "Sort 3D heterogeneity based on the reproducible members of K-means and Equal K-means classification. It runs after 3D refinement where the alignment parameters are determined."; sxcmd.mpi_support = True; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_3d_clustering"; sxcmd.role = "sxr_pipe"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "stack"; token.key_prefix = ""; token.label = "Input images stack"; token.help = ""; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "outdir"; token.key_prefix = ""; token.label = "Output directory"; token.help = "There is a log.txt that describes the sequences of computations in the program. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "output"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "mask"; token.key_prefix = ""; token.label = "3D mask"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "focus"; token.key_prefix = "--"; token.label = "Focus 3D mask"; token.help = "Mask used for focused clustering "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ir"; token.key_prefix = "--"; token.label = "Inner radius for rotational correlation [Pixels]"; token.help = "Must be bigger than 1. "; token.group = "advanced"; token.is_required = False; token.default = "1"; token.restore = "1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "radius"; token.key_prefix = "--"; token.label = "Outer radius for rotational correlation [Pixels]"; token.help = "Must be smaller than half the box size. Please set to the radius of the particle. "; token.group = "main"; token.is_required = False; token.default = "-1"; token.restore = "-1"; token.type = "radius"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "maxit"; token.key_prefix = "--"; token.label = "Maximum iterations"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "25"; token.restore = "25"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "rs"; token.key_prefix = "--"; token.label = "Step between rings in rotational correlation"; token.help = "Must be bigger than 0. "; token.group = "advanced"; token.is_required = False; token.default = "1"; token.restore = "1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "xr"; token.key_prefix = "--"; token.label = "X search range [Pixels]"; token.help = "The translational search range in the x direction will take place in a +/xr range. "; token.group = "advanced"; token.is_required = False; token.default = "1"; token.restore = "1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "yr"; token.key_prefix = "--"; token.label = "Y search range [Pixels]"; token.help = "The translational search range in the y direction. If omitted it will be set as xr. "; token.group = "advanced"; token.is_required = False; token.default = "-1"; token.restore = "-1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ts"; token.key_prefix = "--"; token.label = "Translational search step [Pixels]"; token.help = "The search will be performed in -xr, -xr+ts, 0, xr-ts, xr, can be fractional. "; token.group = "advanced"; token.is_required = False; token.default = "0.25"; token.restore = "0.25"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "delta"; token.key_prefix = "--"; token.label = "Angular step for projections"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "'2'"; token.restore = "'2'"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "an"; token.key_prefix = "--"; token.label = "Local angular search width [Degrees]"; token.help = "This defines the neighbourhood where the local angular search will be performed. "; token.group = "advanced"; token.is_required = False; token.default = "'-1'"; token.restore = "'-1'"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "center"; token.key_prefix = "--"; token.label = "Centering method"; token.help = "0 - if you do not want the volume to be centered, 1 - center the volume using the center of gravity. "; token.group = "advanced"; token.is_required = False; token.default = "0"; token.restore = "0"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "nassign"; token.key_prefix = "--"; token.label = "Number of reassignment iterations"; token.help = "Performed for each angular step. "; token.group = "advanced"; token.is_required = False; token.default = "1"; token.restore = "1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "nrefine"; token.key_prefix = "--"; token.label = "Number of alignment iterations"; token.help = "Performed for each angular step. "; token.group = "advanced"; token.is_required = False; token.default = "0"; token.restore = "0"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "CTF"; token.key_prefix = "--"; token.label = "Use CTF"; token.help = "Do a full CTF correction during the alignment. "; token.group = "main"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "stoprnct"; token.key_prefix = "--"; token.label = "Assignment convergence threshold [%]"; token.help = "Used to asses convergence of the run. It is the minimum percentage of assignment change required to stop the run.  "; token.group = "advanced"; token.is_required = False; token.default = "3.0"; token.restore = "3.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "sym"; token.key_prefix = "--"; token.label = "Point-group symmetry"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "c1"; token.restore = "c1"; token.type = "sym"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "function"; token.key_prefix = "--"; token.label = "Reference preparation function"; token.help = "Function used to prepare the reference volume. "; token.group = "advanced"; token.is_required = False; token.default = "do_volume_mrk02"; token.restore = "do_volume_mrk02"; token.type = "function"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "independent"; token.key_prefix = "--"; token.label = "Number of independent runs"; token.help = "Number of independent equal-Kmeans"; token.group = "advanced"; token.is_required = False; token.default = "3"; token.restore = "3"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "number_of_images_per_group"; token.key_prefix = "--"; token.label = "Images per group"; token.help = "Critical number of images per group, defined by user. "; token.group = "main"; token.is_required = False; token.default = "1000"; token.restore = "1000"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "low_pass_filter"; token.key_prefix = "--"; token.label = "Low-pass filter frequency [1/Pixel]"; token.help = "Low-pass filter used for the 3D sorting on the original image size. "; token.group = "advanced"; token.is_required = False; token.default = "-1.0"; token.restore = "-1.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "nxinit"; token.key_prefix = "--"; token.label = "Initial image size for sorting"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "64"; token.restore = "64"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "unaccounted"; token.key_prefix = "--"; token.label = "Reconstruct unaccounted images"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "seed"; token.key_prefix = "--"; token.label = "Seed"; token.help = "Seed used for the initial random assignment for EQ Kmeans. The program generates a random integer by default. "; token.group = "advanced"; token.is_required = False; token.default = "-1"; token.restore = "-1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "smallest_group"; token.key_prefix = "--"; token.label = "Smallest group size"; token.help = "Minimum members for identified group. "; token.group = "main"; token.is_required = False; token.default = "500"; token.restore = "500"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "sausage"; token.key_prefix = "--"; token.label = "Use sausage filter"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "chunk0"; token.key_prefix = "--"; token.label = "Chunk file name for 1st subset"; token.help = "Name of chunk file containing particle IDs of 1st subset (chunk0) for computing margin of error. "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "parameters"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "chunk1"; token.key_prefix = "--"; token.label = "Chunk file name for 2nd subset"; token.help = "Name of chunk file containing particle IDs of 2nd subset (chunk0) for computing margin of error. "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "parameters"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "PWadjustment"; token.key_prefix = "--"; token.label = "Power spectrum reference"; token.help = "Text file containing a 1D reference power spectrum. "; token.group = "advanced"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "parameters"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "protein_shape"; token.key_prefix = "--"; token.label = "Protein Shape"; token.help = "It defines protein preferred orientation angles. 'g' is for globular proteins and 'f' is for filament proteins. "; token.group = "advanced"; token.is_required = False; token.default = "'g'"; token.restore = "'g'"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "upscale"; token.key_prefix = "--"; token.label = "Power spectrum adjustment strength"; token.help = "This parameters adjusts how strongly the power spectrum of the volume should be modified to match the reference. A value of 1 brings the volume's power spectrum completely to the reference, while a value of 0 means no modification. "; token.group = "advanced"; token.is_required = False; token.default = "0.5"; token.restore = "0.5"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "wn"; token.key_prefix = "--"; token.label = "Target image size [Pixels]"; token.help = "If different than 0, then the images will be rescaled to fit this size. "; token.group = "advanced"; token.is_required = False; token.default = "0"; token.restore = "0"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "interpolation"; token.key_prefix = "--"; token.label = "3D interpolation method"; token.help = "Method interpolation in 3D. Options are tr1 or 4nn. "; token.group = "advanced"; token.is_required = False; token.default = "'4nn'"; token.restore = "'4nn'"; token.type = "string"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxrsort3d"; sxcmd.mode = ""; sxcmd.label = "3D Clustering Protocol II"; sxcmd.short_info = "Sort out 3D heterogeneity of 2D data whose 3D reconstruction parameters (xform.projection) have been determined already using 3D sorting protocol I (P1)."; sxcmd.mpi_support = True; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_3d_clustering"; sxcmd.role = "sxr_pipe"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "stack"; token.key_prefix = ""; token.label = "Input images stack"; token.help = ""; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "outdir"; token.key_prefix = ""; token.label = "Output directory"; token.help = "There is a log.txt that describes the sequences of computations in the program. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "output"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "mask"; token.key_prefix = ""; token.label = "3D mask"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "previous_run1"; token.key_prefix = "--"; token.label = "Directory of first sxsort3d.py run"; token.help = ""; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "directory"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "previous_run2"; token.key_prefix = "--"; token.label = "Directory of second sxsort3d.py run"; token.help = ""; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "directory"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "focus"; token.key_prefix = "--"; token.label = "Focus 3D mask"; token.help = "Mask used for focused clustering "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ir"; token.key_prefix = "--"; token.label = "Inner radius for rotational correlation [Pixels]"; token.help = "Must be bigger than 1. "; token.group = "advanced"; token.is_required = False; token.default = "1"; token.restore = "1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "radius"; token.key_prefix = "--"; token.label = "Outer radius for rotational correlation [Pixels]"; token.help = "Must be smaller than half the box size. Please set to the radius of the particle. "; token.group = "main"; token.is_required = False; token.default = "-1"; token.restore = "-1"; token.type = "radius"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "maxit"; token.key_prefix = "--"; token.label = "Maximum iterations"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "25"; token.restore = "25"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "rs"; token.key_prefix = "--"; token.label = "Step between rings in rotational correlation"; token.help = "Must be bigger than 0. "; token.group = "advanced"; token.is_required = False; token.default = "1"; token.restore = "1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "xr"; token.key_prefix = "--"; token.label = "X search range [Pixels]"; token.help = "The translational search range in the x direction will take place in a +/xr range. "; token.group = "advanced"; token.is_required = False; token.default = "1"; token.restore = "1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "yr"; token.key_prefix = "--"; token.label = "Y search range [Pixels]"; token.help = "The translational search range in the y direction. If omitted it will be set as xr. "; token.group = "advanced"; token.is_required = False; token.default = "-1"; token.restore = "-1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ts"; token.key_prefix = "--"; token.label = "Translational search step [Pixels]"; token.help = "The search will be performed in -xr, -xr+ts, 0, xr-ts, xr, can be fractional. "; token.group = "advanced"; token.is_required = False; token.default = "0.25"; token.restore = "0.25"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "delta"; token.key_prefix = "--"; token.label = "Angular step for projections"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "'2'"; token.restore = "'2'"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "an"; token.key_prefix = "--"; token.label = "Local angular search width [Degrees]"; token.help = "This defines the neighbourhood where the local angular search will be performed. "; token.group = "advanced"; token.is_required = False; token.default = "'-1'"; token.restore = "'-1'"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "center"; token.key_prefix = "--"; token.label = "Centering method"; token.help = "0 - if you do not want the volume to be centered, 1 - center the volume using the center of gravity. "; token.group = "advanced"; token.is_required = False; token.default = "0"; token.restore = "0"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "nassign"; token.key_prefix = "--"; token.label = "Number of reassignment iterations"; token.help = "Performed for each angular step. "; token.group = "advanced"; token.is_required = False; token.default = "1"; token.restore = "1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "nrefine"; token.key_prefix = "--"; token.label = "Number of alignment iterations"; token.help = "Performed for each angular step. "; token.group = "advanced"; token.is_required = False; token.default = "0"; token.restore = "0"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "CTF"; token.key_prefix = "--"; token.label = "Use CTF"; token.help = "Do a full CTF correction during the alignment. "; token.group = "main"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "stoprnct"; token.key_prefix = "--"; token.label = "Assignment convergence threshold [%]"; token.help = "Used to asses convergence of the run. It is the minimum percentage of assignment change required to stop the run.  "; token.group = "advanced"; token.is_required = False; token.default = "3.0"; token.restore = "3.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "sym"; token.key_prefix = "--"; token.label = "Point-group symmetry"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "c1"; token.restore = "c1"; token.type = "sym"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "independent"; token.key_prefix = "--"; token.label = "Number of independent runs"; token.help = "Number of independent equal-Kmeans"; token.group = "advanced"; token.is_required = False; token.default = "3"; token.restore = "3"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "number_of_images_per_group"; token.key_prefix = "--"; token.label = "Images per group"; token.help = "Critical number of images per group, defined by user. "; token.group = "main"; token.is_required = False; token.default = "1000"; token.restore = "1000"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "low_pass_filter"; token.key_prefix = "--"; token.label = "Low-pass filter frequency [1/Pixel]"; token.help = "Low-pass filter used for the 3D sorting on the original image size. "; token.group = "advanced"; token.is_required = False; token.default = "-1.0"; token.restore = "-1.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "nxinit"; token.key_prefix = "--"; token.label = "Initial image size for sorting"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "64"; token.restore = "64"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "unaccounted"; token.key_prefix = "--"; token.label = "Reconstruct unaccounted images"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "function"; token.key_prefix = "--"; token.label = "Reference preparation function"; token.help = "Function used to prepare the reference volume. "; token.group = "advanced"; token.is_required = False; token.default = "do_volume_mrk02"; token.restore = "do_volume_mrk02"; token.type = "function"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "seed"; token.key_prefix = "--"; token.label = "Seed"; token.help = "Seed used for the initial random assignment for EQ Kmeans. The program generates a random integer by default. "; token.group = "advanced"; token.is_required = False; token.default = "-1"; token.restore = "-1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "smallest_group"; token.key_prefix = "--"; token.label = "Smallest group size"; token.help = "Minimum members for identified group. "; token.group = "main"; token.is_required = False; token.default = "500"; token.restore = "500"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "group_size_for_unaccounted"; token.key_prefix = "--"; token.label = "Unaccounted particles group size"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "sausage"; token.key_prefix = "--"; token.label = "Use sausage filter"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "chunkdir"; token.key_prefix = "--"; token.label = "Chunkdir"; token.help = "Used for computing margin of error. "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "parameters"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "PWadjustment"; token.key_prefix = "--"; token.label = "Power spectrum reference"; token.help = "Text file containing a 1D reference power spectrum. "; token.group = "advanced"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "parameters"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "protein_shape"; token.key_prefix = "--"; token.label = "Protein Shape"; token.help = "It defines protein preferred orientation angles. 'g' is for globular proteins and 'f' is for filament proteins. "; token.group = "advanced"; token.is_required = False; token.default = "'g'"; token.restore = "'g'"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "upscale"; token.key_prefix = "--"; token.label = "Power spectrum adjustment strength"; token.help = "This parameters adjusts how strongly the power spectrum of the volume should be modified to match the reference. A value of 1 brings the volume's power spectrum completely to the reference, while a value of 0 means no modification.  "; token.group = "advanced"; token.is_required = False; token.default = "0.5"; token.restore = "0.5"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "wn"; token.key_prefix = "--"; token.label = "Target image size [Pixels]"; token.help = "If different than 0, then the images will be rescaled to fit this size. "; token.group = "advanced"; token.is_required = False; token.default = "0"; token.restore = "0"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "interpolation"; token.key_prefix = "--"; token.label = "3D interpolation method"; token.help = "Method interpolation in 3D. Options are tr1 or 4nn. "; token.group = "advanced"; token.is_required = False; token.default = "'4nn'"; token.restore = "'4nn'"; token.type = "string"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "e2display"; sxcmd.mode = ""; sxcmd.label = "Display Data"; sxcmd.short_info = "Displays images, volumes, or 1D plots."; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_3d_clustering"; sxcmd.role = "sxr_util"; sxcmd.is_submittable = False
		token = SXcmd_token(); token.key_base = "input_data_list"; token.key_prefix = ""; token.label = "Input files"; token.help = "List of input images, volumes, plots. Wild cards (e.g *) can be used to select a list of files. Not recommended when the list is too large. "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "any_file_list"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "singleimage"; token.key_prefix = "--"; token.label = "Single image view"; token.help = "Displays a stack in a single image view: "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "fullrange"; token.key_prefix = "--"; token.label = "Skip auto-contrast"; token.help = "Disables the auto contrast for the display of particles stacks and 2D images. "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "verbose"; token.key_prefix = "--"; token.label = "Verbose"; token.help = "Accepted values 0-9. "; token.group = "advanced"; token.is_required = False; token.default = "0"; token.restore = "0"; token.type = "int"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxprocess"; sxcmd.mode = "adaptive_mask"; sxcmd.label = "Adaptive 3D Mask"; sxcmd.short_info = "Create soft-edged 3D mask from reference volume. "; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_3d_clustering"; sxcmd.role = "sxr_util"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "adaptive_mask"; token.key_prefix = "--"; token.label = "Create soft-edged 3D mask from reference volume"; token.help = ""; token.group = "main"; token.is_required = True; token.default = True; token.restore = True; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "input_volume"; token.key_prefix = ""; token.label = "Input volume"; token.help = "Input reference volume"; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "output_mask3D"; token.key_prefix = ""; token.label = "Output mask"; token.help = "Output 3D mask"; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "output"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "nsigma"; token.key_prefix = "--"; token.label = "Density standard deviations threshold"; token.help = "Defines the threshold used to find the main volume within the data. All voxels with density <= mean + nsigma standard deviations will be included into the main volume. This option will not be used if the option threshold is larger than -9999.0. "; token.group = "main"; token.is_required = False; token.default = "1.0"; token.restore = "1.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "threshold"; token.key_prefix = "--"; token.label = "Binarization threshold"; token.help = "Below this value the data is assumed to not belong to the main volume. With the value lower than the default, the option will be ignored and the mask will be set according to nsigma. "; token.group = "main"; token.is_required = False; token.default = "-9999.0"; token.restore = "-9999.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ndilation"; token.key_prefix = "--"; token.label = "Mask extension cycles"; token.help = "The initial mask will be extended this number of cycles. To keep the size of the main volume, set this to kernel_size/2 "; token.group = "main"; token.is_required = False; token.default = "3"; token.restore = "3"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "kernel_size"; token.key_prefix = "--"; token.label = "Gaussian kernel size"; token.help = "Size of the gaussian kernel used to smooth the binary mask. "; token.group = "main"; token.is_required = False; token.default = "11"; token.restore = "11"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "gauss_standard_dev"; token.key_prefix = "--"; token.label = "Kernel standard deviation "; token.help = "Standard deviation used in the construction of the gaussian smoothing of the mask. "; token.group = "main"; token.is_required = False; token.default = "9"; token.restore = "9"; token.type = "int"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxprocess"; sxcmd.mode = "binary_mask"; sxcmd.label = "Binary 3D Mask"; sxcmd.short_info = "Create binary 3D mask from reference volume. "; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_3d_clustering"; sxcmd.role = "sxr_util"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "binary_mask"; token.key_prefix = "--"; token.label = "Create binary 3D mask from reference volume"; token.help = ""; token.group = "main"; token.is_required = True; token.default = True; token.restore = True; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "input_volume"; token.key_prefix = ""; token.label = "Input volume"; token.help = "Input reference volume"; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "output_mask3D"; token.key_prefix = ""; token.label = "Output mask"; token.help = "Output 3D mask"; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "output"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "bin_threshold"; token.key_prefix = "--"; token.label = "Binarization threshold"; token.help = "Below this value the data is assumed to not belong to the main volume. "; token.group = "main"; token.is_required = False; token.default = "0.0"; token.restore = "0.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ne"; token.key_prefix = "--"; token.label = "Erosion cycles"; token.help = "After initial binarization the volume is eroded to remove fragmented pieces of the volume. "; token.group = "main"; token.is_required = False; token.default = "0"; token.restore = "0"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "nd"; token.key_prefix = "--"; token.label = "Dilation cycles"; token.help = "After erosing the binary volume is dilated back to smooth the surface and match the original size. "; token.group = "main"; token.is_required = False; token.default = "0"; token.restore = "0"; token.type = "int"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxprocess"; sxcmd.mode = "angular_distribution"; sxcmd.label = "Angular Distribution"; sxcmd.short_info = "Create angular distribution file. The file will contain a 3D representation of the given angular distribution, which can be viewed with UCFS Chimera "; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_3d_clustering"; sxcmd.role = "sxr_util"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "angular_distribution"; token.key_prefix = "--"; token.label = "Create angular distribution file"; token.help = "The file will contain a 3D representation of the given angular distribution, which can be viewed with UCFS Chimera "; token.group = "main"; token.is_required = True; token.default = True; token.restore = True; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "inputfile"; token.key_prefix = ""; token.label = "Alignment Parameter file"; token.help = "Alignment Parameter file created by a previous 3D reconstruction step (e.g. sxmeridien.py)"; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "any_file"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "pixel_size"; token.key_prefix = "--"; token.label = "Pixel size [A]"; token.help = "Pixel size in A. "; token.group = "main"; token.is_required = False; token.default = "1.0"; token.restore = "1.0"; token.type = "apix"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "round_digit"; token.key_prefix = "--"; token.label = "Number precision"; token.help = "Decimal numbers will be rounded to this number of decimal points. "; token.group = "main"; token.is_required = False; token.default = "5"; token.restore = "5"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "box_size"; token.key_prefix = "--"; token.label = "Box size [Pixels]"; token.help = "Box size in pixel used for calculating the center of the particle. "; token.group = "main"; token.is_required = False; token.default = "500"; token.restore = "500"; token.type = "box"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "particle_radius"; token.key_prefix = "--"; token.label = "Particle radius [Pixels]"; token.help = "Used for the representation in Chimera. Defines where the cylinders representing the histogram must start. "; token.group = "main"; token.is_required = False; token.default = "175"; token.restore = "175"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "cylinder_width"; token.key_prefix = "--"; token.label = "Cylinder width"; token.help = "Used for the representation in Chimera. This will define the width of the cylinders representing the histogram."; token.group = "main"; token.is_required = False; token.default = "1"; token.restore = "1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "cylinder_length"; token.key_prefix = "--"; token.label = "Cylinder length"; token.help = "Used for the representation in Chimera. This will define the relative size of the cylinders representing the histogram. "; token.group = "main"; token.is_required = False; token.default = "10000"; token.restore = "10000"; token.type = "int"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxlocres"; sxcmd.mode = ""; sxcmd.label = "Local Resolution"; sxcmd.short_info = "Compute local resolution of a map."; sxcmd.mpi_support = True; sxcmd.mpi_add_flag = True; sxcmd.category = "sxc_localres"; sxcmd.role = "sxr_pipe"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "firstvolume"; token.key_prefix = ""; token.label = "First half-volume"; token.help = ""; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "secondvolume"; token.key_prefix = ""; token.label = "Second half-volume"; token.help = ""; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "maskfile"; token.key_prefix = ""; token.label = "3D mask"; token.help = "Defines the region where the local filtering should be applied. "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "outputfile"; token.key_prefix = ""; token.label = "Output volume"; token.help = "Each voxel contains the resolution for this area in absolute frequency units. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "output"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "wn"; token.key_prefix = "--"; token.label = "FSC window size"; token.help = "Defines the size of window where the local real-space FSC is computed. "; token.group = "main"; token.is_required = False; token.default = "7"; token.restore = "7"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "step"; token.key_prefix = "--"; token.label = "Fourier shell step size [Pixels]"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "1.0"; token.restore = "1.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "cutoff"; token.key_prefix = "--"; token.label = "Resolution cut-off [1/Pixel]"; token.help = "Cut-off for the FSC. By default it goes to the Nyquist frequencies. "; token.group = "main"; token.is_required = False; token.default = "0.5"; token.restore = "0.5"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "res_overall"; token.key_prefix = "--"; token.label = "Overall resolution [1/Pixel]"; token.help = "A value between 0 to 0.5 for calibration of the average local resolution. See Description section in the wiki page for details. "; token.group = "main"; token.is_required = False; token.default = "-1.0"; token.restore = "-1.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "radius"; token.key_prefix = "--"; token.label = "Mask radius [Pixels]"; token.help = "In case no mask is provided, a hard sphere of this radius will be used. "; token.group = "main"; token.is_required = False; token.default = "-1"; token.restore = "-1"; token.type = "radius"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "fsc"; token.key_prefix = "--"; token.label = "FSC output file"; token.help = "Contains the overall FSC curve computed by rotational averaging of local resolution values. It is truncated to res_overall. "; token.group = "main"; token.is_required = False; token.default = "no curve"; token.restore = "no curve"; token.type = "string"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxfilterlocal"; sxcmd.mode = ""; sxcmd.label = "3D Local Filter"; sxcmd.short_info = "Locally filter maps according to the local resolution determined by sxlocres.py."; sxcmd.mpi_support = True; sxcmd.mpi_add_flag = True; sxcmd.category = "sxc_localres"; sxcmd.role = "sxr_pipe"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "input_volume"; token.key_prefix = ""; token.label = "Input volume"; token.help = ""; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "locres_volume"; token.key_prefix = ""; token.label = "Local resolution file"; token.help = "Volume file containing the local resolution estimate produced by sxlocres. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "maskfile"; token.key_prefix = ""; token.label = "3D mask"; token.help = "Defines the region where the local filtering should be applied. "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "outputfile"; token.key_prefix = ""; token.label = "Output volume"; token.help = ""; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "output"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "radius"; token.key_prefix = "--"; token.label = "Mask radius [Pixels]"; token.help = "In case no mask is provided, a hard sphere of this radius will be used. "; token.group = "main"; token.is_required = False; token.default = "-1"; token.restore = "-1"; token.type = "radius"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "falloff"; token.key_prefix = "--"; token.label = "Low-pass filter fall-off"; token.help = "The program uses a tangent low-pass filter. "; token.group = "main"; token.is_required = False; token.default = "0.1"; token.restore = "0.1"; token.type = "float"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "e2display"; sxcmd.mode = ""; sxcmd.label = "Display Data"; sxcmd.short_info = "Displays images, volumes, or 1D plots."; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_localres"; sxcmd.role = "sxr_util"; sxcmd.is_submittable = False
		token = SXcmd_token(); token.key_base = "input_data_list"; token.key_prefix = ""; token.label = "Input files"; token.help = "List of input images, volumes, plots. Wild cards (e.g *) can be used to select a list of files. Not recommended when the list is too large. "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "any_file_list"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "singleimage"; token.key_prefix = "--"; token.label = "Single image view"; token.help = "Displays a stack in a single image view: "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "fullrange"; token.key_prefix = "--"; token.label = "Skip auto-contrast"; token.help = "Disables the auto contrast for the display of particles stacks and 2D images. "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "verbose"; token.key_prefix = "--"; token.label = "Verbose"; token.help = "Accepted values 0-9. "; token.group = "advanced"; token.is_required = False; token.default = "0"; token.restore = "0"; token.type = "int"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxprocess"; sxcmd.mode = "adaptive_mask"; sxcmd.label = "Adaptive 3D Mask"; sxcmd.short_info = "Create soft-edged 3D mask from reference volume. "; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_localres"; sxcmd.role = "sxr_util"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "adaptive_mask"; token.key_prefix = "--"; token.label = "Create soft-edged 3D mask from reference volume"; token.help = ""; token.group = "main"; token.is_required = True; token.default = True; token.restore = True; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "input_volume"; token.key_prefix = ""; token.label = "Input volume"; token.help = "Input reference volume"; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "output_mask3D"; token.key_prefix = ""; token.label = "Output mask"; token.help = "Output 3D mask"; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "output"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "nsigma"; token.key_prefix = "--"; token.label = "Density standard deviations threshold"; token.help = "Defines the threshold used to find the main volume within the data. All voxels with density <= mean + nsigma standard deviations will be included into the main volume. This option will not be used if the option threshold is larger than -9999.0. "; token.group = "main"; token.is_required = False; token.default = "1.0"; token.restore = "1.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "threshold"; token.key_prefix = "--"; token.label = "Binarization threshold"; token.help = "Below this value the data is assumed to not belong to the main volume. With the value lower than the default, the option will be ignored and the mask will be set according to nsigma. "; token.group = "main"; token.is_required = False; token.default = "-9999.0"; token.restore = "-9999.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ndilation"; token.key_prefix = "--"; token.label = "Mask extension cycles"; token.help = "The initial mask will be extended this number of cycles. To keep the size of the main volume, set this to kernel_size/2 "; token.group = "main"; token.is_required = False; token.default = "3"; token.restore = "3"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "kernel_size"; token.key_prefix = "--"; token.label = "Gaussian kernel size"; token.help = "Size of the gaussian kernel used to smooth the binary mask. "; token.group = "main"; token.is_required = False; token.default = "11"; token.restore = "11"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "gauss_standard_dev"; token.key_prefix = "--"; token.label = "Kernel standard deviation "; token.help = "Standard deviation used in the construction of the gaussian smoothing of the mask. "; token.group = "main"; token.is_required = False; token.default = "9"; token.restore = "9"; token.type = "int"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxprocess"; sxcmd.mode = "angular_distribution"; sxcmd.label = "Angular Distribution"; sxcmd.short_info = "Create angular distribution file. The file will contain a 3D representation of the given angular distribution, which can be viewed with UCFS Chimera "; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_localres"; sxcmd.role = "sxr_util"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "angular_distribution"; token.key_prefix = "--"; token.label = "Create angular distribution file"; token.help = "The file will contain a 3D representation of the given angular distribution, which can be viewed with UCFS Chimera "; token.group = "main"; token.is_required = True; token.default = True; token.restore = True; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "inputfile"; token.key_prefix = ""; token.label = "Alignment Parameter file"; token.help = "Alignment Parameter file created by a previous 3D reconstruction step (e.g. sxmeridien.py)"; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "any_file"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "pixel_size"; token.key_prefix = "--"; token.label = "Pixel size [A]"; token.help = "Pixel size in A. "; token.group = "main"; token.is_required = False; token.default = "1.0"; token.restore = "1.0"; token.type = "apix"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "round_digit"; token.key_prefix = "--"; token.label = "Number precision"; token.help = "Decimal numbers will be rounded to this number of decimal points. "; token.group = "main"; token.is_required = False; token.default = "5"; token.restore = "5"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "box_size"; token.key_prefix = "--"; token.label = "Box size [Pixels]"; token.help = "Box size in pixel used for calculating the center of the particle. "; token.group = "main"; token.is_required = False; token.default = "500"; token.restore = "500"; token.type = "box"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "particle_radius"; token.key_prefix = "--"; token.label = "Particle radius [Pixels]"; token.help = "Used for the representation in Chimera. Defines where the cylinders representing the histogram must start. "; token.group = "main"; token.is_required = False; token.default = "175"; token.restore = "175"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "cylinder_width"; token.key_prefix = "--"; token.label = "Cylinder width"; token.help = "Used for the representation in Chimera. This will define the width of the cylinders representing the histogram."; token.group = "main"; token.is_required = False; token.default = "1"; token.restore = "1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "cylinder_length"; token.key_prefix = "--"; token.label = "Cylinder length"; token.help = "Used for the representation in Chimera. This will define the relative size of the cylinders representing the histogram. "; token.group = "main"; token.is_required = False; token.default = "10000"; token.restore = "10000"; token.type = "int"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "e2display"; sxcmd.mode = ""; sxcmd.label = "Display Data"; sxcmd.short_info = "Displays images, volumes, or 1D plots."; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_utilities"; sxcmd.role = "sxr_util"; sxcmd.is_submittable = False
		token = SXcmd_token(); token.key_base = "input_data_list"; token.key_prefix = ""; token.label = "Input files"; token.help = "List of input images, volumes, plots. Wild cards (e.g *) can be used to select a list of files. Not recommended when the list is too large. "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "any_file_list"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "singleimage"; token.key_prefix = "--"; token.label = "Single image view"; token.help = "Displays a stack in a single image view: "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "fullrange"; token.key_prefix = "--"; token.label = "Skip auto-contrast"; token.help = "Disables the auto contrast for the display of particles stacks and 2D images. "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "verbose"; token.key_prefix = "--"; token.label = "Verbose"; token.help = "Accepted values 0-9. "; token.group = "advanced"; token.is_required = False; token.default = "0"; token.restore = "0"; token.type = "int"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxpdb2em"; sxcmd.mode = ""; sxcmd.label = "PDB File Conversion"; sxcmd.short_info = "Converts an atomic model into a simulated density map"; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_utilities"; sxcmd.role = "sxr_util"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "input_pdb"; token.key_prefix = ""; token.label = "Input PDB file"; token.help = "Starting atomic coordinates: "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "pdb"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "output_hdf"; token.key_prefix = ""; token.label = "Output map"; token.help = ""; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "output"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "apix"; token.key_prefix = "--"; token.label = "Pixel size of the output map [A]"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "1.0"; token.restore = "1.0"; token.type = "apix"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "box"; token.key_prefix = "--"; token.label = "Output box size [Voxels]"; token.help = "If not given, the program will find the minimum box size fitting the structure. Be aware that this will most likely result in a rectangular box. "; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "box"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "het"; token.key_prefix = "--"; token.label = "Include hetero atoms"; token.help = "Otherwise the HETATM entries in the PDB file are ignored. "; token.group = "main"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "center"; token.key_prefix = "--"; token.label = "Center model at the origin"; token.help = ": Specifies whether the atomic model should be moved to the origin before generating density map.  Available options are: c - Use the geometrical center of atoms; a - Use the center of mass (recommended); x,y,z - Vector to be subtracted for the coordinates. Default: no centering, in which case (0,0,0) in the PDB space will map to the center of the EM volume. "; token.group = "main"; token.is_required = False; token.default = "n"; token.restore = "n"; token.type = "string"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "O"; token.key_prefix = "--"; token.label = "Apply additional rotation"; token.help = "This can be used to modify the orientation of the atomic model. "; token.group = "main"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "tr0"; token.key_prefix = "--"; token.label = "Rotational matrix file"; token.help = ": This file must contain the 3x4 transformation matrix to be applied to the PDB coordinates after centering. The translation vector (last column of the matrix) must be specified in Angstrom. "; token.group = "main"; token.is_required = False; token.default = "none"; token.restore = "none"; token.type = "parameters"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "quiet"; token.key_prefix = "--"; token.label = "Silent mode"; token.help = "Does not print any information to the monitor. "; token.group = "advanced"; token.is_required = False; token.default = False; token.restore = False; token.type = "bool"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxprocess"; sxcmd.mode = "adaptive_mask"; sxcmd.label = "Adaptive 3D Mask"; sxcmd.short_info = "Create soft-edged 3D mask from reference volume. "; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_utilities"; sxcmd.role = "sxr_util"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "adaptive_mask"; token.key_prefix = "--"; token.label = "Create soft-edged 3D mask from reference volume"; token.help = ""; token.group = "main"; token.is_required = True; token.default = True; token.restore = True; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "input_volume"; token.key_prefix = ""; token.label = "Input volume"; token.help = "Input reference volume"; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "output_mask3D"; token.key_prefix = ""; token.label = "Output mask"; token.help = "Output 3D mask"; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "output"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "nsigma"; token.key_prefix = "--"; token.label = "Density standard deviations threshold"; token.help = "Defines the threshold used to find the main volume within the data. All voxels with density <= mean + nsigma standard deviations will be included into the main volume. This option will not be used if the option threshold is larger than -9999.0. "; token.group = "main"; token.is_required = False; token.default = "1.0"; token.restore = "1.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "threshold"; token.key_prefix = "--"; token.label = "Binarization threshold"; token.help = "Below this value the data is assumed to not belong to the main volume. With the value lower than the default, the option will be ignored and the mask will be set according to nsigma. "; token.group = "main"; token.is_required = False; token.default = "-9999.0"; token.restore = "-9999.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ndilation"; token.key_prefix = "--"; token.label = "Mask extension cycles"; token.help = "The initial mask will be extended this number of cycles. To keep the size of the main volume, set this to kernel_size/2 "; token.group = "main"; token.is_required = False; token.default = "3"; token.restore = "3"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "kernel_size"; token.key_prefix = "--"; token.label = "Gaussian kernel size"; token.help = "Size of the gaussian kernel used to smooth the binary mask. "; token.group = "main"; token.is_required = False; token.default = "11"; token.restore = "11"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "gauss_standard_dev"; token.key_prefix = "--"; token.label = "Kernel standard deviation "; token.help = "Standard deviation used in the construction of the gaussian smoothing of the mask. "; token.group = "main"; token.is_required = False; token.default = "9"; token.restore = "9"; token.type = "int"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxprocess"; sxcmd.mode = "binary_mask"; sxcmd.label = "Binary 3D Mask"; sxcmd.short_info = "Create binary 3D mask from reference volume. "; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_utilities"; sxcmd.role = "sxr_util"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "binary_mask"; token.key_prefix = "--"; token.label = "Create binary 3D mask from reference volume"; token.help = ""; token.group = "main"; token.is_required = True; token.default = True; token.restore = True; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "input_volume"; token.key_prefix = ""; token.label = "Input volume"; token.help = "Input reference volume"; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "image"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "output_mask3D"; token.key_prefix = ""; token.label = "Output mask"; token.help = "Output 3D mask"; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "output"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "bin_threshold"; token.key_prefix = "--"; token.label = "Binarization threshold"; token.help = "Below this value the data is assumed to not belong to the main volume. "; token.group = "main"; token.is_required = False; token.default = "0.0"; token.restore = "0.0"; token.type = "float"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "ne"; token.key_prefix = "--"; token.label = "Erosion cycles"; token.help = "After initial binarization the volume is eroded to remove fragmented pieces of the volume. "; token.group = "main"; token.is_required = False; token.default = "0"; token.restore = "0"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "nd"; token.key_prefix = "--"; token.label = "Dilation cycles"; token.help = "After erosing the binary volume is dilated back to smooth the surface and match the original size. "; token.group = "main"; token.is_required = False; token.default = "0"; token.restore = "0"; token.type = "int"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		sxcmd = SXcmd(); sxcmd.name = "sxprocess"; sxcmd.mode = "angular_distribution"; sxcmd.label = "Angular Distribution"; sxcmd.short_info = "Create angular distribution file. The file will contain a 3D representation of the given angular distribution, which can be viewed with UCFS Chimera "; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False; sxcmd.category = "sxc_utilities"; sxcmd.role = "sxr_util"; sxcmd.is_submittable = True
		token = SXcmd_token(); token.key_base = "angular_distribution"; token.key_prefix = "--"; token.label = "Create angular distribution file"; token.help = "The file will contain a 3D representation of the given angular distribution, which can be viewed with UCFS Chimera "; token.group = "main"; token.is_required = True; token.default = True; token.restore = True; token.type = "bool"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "inputfile"; token.key_prefix = ""; token.label = "Alignment Parameter file"; token.help = "Alignment Parameter file created by a previous 3D reconstruction step (e.g. sxmeridien.py)"; token.group = "main"; token.is_required = True; token.default = ""; token.restore = ""; token.type = "any_file"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "pixel_size"; token.key_prefix = "--"; token.label = "Pixel size [A]"; token.help = "Pixel size in A. "; token.group = "main"; token.is_required = False; token.default = "1.0"; token.restore = "1.0"; token.type = "apix"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "round_digit"; token.key_prefix = "--"; token.label = "Number precision"; token.help = "Decimal numbers will be rounded to this number of decimal points. "; token.group = "main"; token.is_required = False; token.default = "5"; token.restore = "5"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "box_size"; token.key_prefix = "--"; token.label = "Box size [Pixels]"; token.help = "Box size in pixel used for calculating the center of the particle. "; token.group = "main"; token.is_required = False; token.default = "500"; token.restore = "500"; token.type = "box"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "particle_radius"; token.key_prefix = "--"; token.label = "Particle radius [Pixels]"; token.help = "Used for the representation in Chimera. Defines where the cylinders representing the histogram must start. "; token.group = "main"; token.is_required = False; token.default = "175"; token.restore = "175"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "cylinder_width"; token.key_prefix = "--"; token.label = "Cylinder width"; token.help = "Used for the representation in Chimera. This will define the width of the cylinders representing the histogram."; token.group = "main"; token.is_required = False; token.default = "1"; token.restore = "1"; token.type = "int"; sxcmd.token_list.append(token)
		token = SXcmd_token(); token.key_base = "cylinder_length"; token.key_prefix = "--"; token.label = "Cylinder length"; token.help = "Used for the representation in Chimera. This will define the relative size of the cylinders representing the histogram. "; token.group = "main"; token.is_required = False; token.default = "10000"; token.restore = "10000"; token.type = "int"; sxcmd.token_list.append(token)

		sxcmd_list.append(sxcmd)

		# @@@@@ END_INSERTION @@@@@

		# Create dictionaries from the constructed lists
		for sxcmd_category in sxcmd_category_list:
			sxcmd_category_dict[sxcmd_category.name] = sxcmd_category

		# Create command token dictionary for each SXcmd instance
		# Then, register SXcmd instance to an associated SXcmd_category
		for sxcmd in sxcmd_list:
			for sxcmd_token in sxcmd.token_list:
				# Handle very special cases
				if sxcmd_token.type == "function":
					n_widgets = 2 # function type has two line edit boxes
					sxcmd_token.label = [sxcmd_token.label, "Python script for user function"]
					sxcmd_token.help = [sxcmd_token.help, "Please leave it blank if file is not external to sphire"]
					sxcmd_token.default = [sxcmd_token.default, "none"]
					sxcmd_token.restore = sxcmd_token.default
				# else: Do nothing for the other types

				# Register this to command token dictionary
				sxcmd.token_dict[sxcmd_token.key_base] = sxcmd_token

			# Register this to command to command category dictionary
			assert sxcmd_category_dict.has_key(sxcmd.category), "sxcmd.category %s" % (sxcmd.category)
			sxcmd_category_dict[sxcmd.category].cmd_list.append(sxcmd)

		# Store the constructed lists and dictionary as a class data member
		self.sxcmd_category_list = sxcmd_category_list

	def update_qsub_enable_states(self):
		# Construct and add widgets for sx command categories
		for sxcmd_category in self.sxcmd_category_list:
			# Create SXCmdCategoryWidget for this command category
			for sxcmd in sxcmd_category.cmd_list:
				sxcmd.widget.sxcmd_tab_main.set_qsub_enable_state()

	def handle_sxmenu_item_btn_event(self, sxmenu_item):
		assert(isinstance(sxmenu_item, SXmenu_item) == True) # Assuming the sxmenu_item is an instance of class SXmenu_item

		if self.cur_sxmenu_item == sxmenu_item: return

		if self.cur_sxmenu_item != None:
			self.cur_sxmenu_item.btn.setStyleSheet(self.cur_sxmenu_item.btn.customButtonStyle)

		self.cur_sxmenu_item = sxmenu_item

		if self.cur_sxmenu_item != None:
			self.cur_sxmenu_item.btn.setStyleSheet(self.cur_sxmenu_item.btn.customButtonStyleClicked)
			self.sxmenu_item_widget_stacked_layout.setCurrentWidget(self.cur_sxmenu_item.widget)

	def closeEvent(self, event):
		event.ignore() # event.accept()

		# Quit child applications of all sxcmd widgets
		for sxcmd_category in self.sxcmd_category_list:
			sxcmd_category.widget.quit_all_child_applications()

		print("bye bye")
		QtCore.QCoreApplication.instance().quit()

	# def changeEvent(self, event):
	# 	print(self.frameGeometry())

# ========================================================================================
def main(args):
	sxapp = QApplication(args)
	# The valid keys can be retrieved using the keys() function.
	# Typically they include "windows", "motif", "cde", "plastique" and "cleanlooks".
	# Depending on the platform, "windowsxp", "windowsvista" and "macintosh" may be available. Note that keys are case insensitive.
	# sxapp.setStyle("macintosh")
	sxapp.setStyle("cleanlooks")
	# sxapp.setStyle("plastique")

	# print "MRK_DEBUG:"
	# print "MRK_DEBUG: sxapp.style().metaObject().className() == %s" % (str(sxapp.style().metaObject().className()))
	# for key in QStyleFactory.keys():
	# 	print "MRK_DEBUG: str(key) == %s" % str(key)
	# 	print "MRK_DEBUG: QStyleFactory.create(key) = %s" % (str(QStyleFactory.create(key).metaObject().className()))
	# 	if sxapp.style().metaObject().className() == QStyleFactory.create(key).metaObject().className():
	# 		print "MRK_DEBUG: !!!USING THE STYLE: %s!!!" % str(key)
	# print "MRK_DEBUG:"

	sxapp.setWindowIcon(QIcon(get_image_directory()+"sxgui_icon_sphire.png"))

	sxapp_font = sxapp.font()
	sxapp_font_info = QFontInfo(sxapp.font())
#	new_point_size = sxapp_font_info.pointSize() + 1
	new_point_size = sxapp_font_info.pointSize()
	# # MRK_DEBUG: Check the default system font
	# print "MRK_DEBUG: sxapp_font_info.style()      = ", sxapp_font_info.style()
	# print "MRK_DEBUG: sxapp_font_info.styleHint()  = ", sxapp_font_info.styleHint()
	# print "MRK_DEBUG: sxapp_font_info.styleName()  = ", sxapp_font_info.styleName()
	# print "MRK_DEBUG: sxapp_font_info.family()     = ", sxapp_font_info.family()
	# print "MRK_DEBUG: sxapp_font_info.fixedPitch() = ", sxapp_font_info.fixedPitch()
	# print "MRK_DEBUG: sxapp_font_info.pixelSize()  = ", sxapp_font_info.pixelSize()
	# print "MRK_DEBUG: sxapp_font_info.pointSize()  = ", sxapp_font_info.pointSize()
	# print "MRK_DEBUG: sxapp_font_info.pointSizeF() = ", sxapp_font_info.pointSizeF()
	# print "MRK_DEBUG: sxapp_font_info.bold ()      = ", sxapp_font_info.bold()
	# print "MRK_DEBUG: sxapp_font_info.italic()     = ", sxapp_font_info.italic()
	#
	# NOTE: 2019/02/19 Toshio Moriya
	# The following method of changing font size works with Linux.
	# However, it does not work Mac OSX. The text of widget classes below won't change,
	# still showing the default font size:
	# QPushButton, QLable, Window Title, and QToolTip
	#
#	sxapp_font.setPointSize(new_point_size) # and setPointSizeF() are device independent, while setPixelSize() is device dependent
#	sxapp.setFont(sxapp_font)

	# sxapp.setStyleSheet("QPushButton {font-size:18pt;}");  # NOTE: 2016/02/19 Toshio Moriya: Doesn't work
	# sxapp.setStyleSheet("QLabel {font-size:18pt;}"); # NOTE: 2016/02/19 Toshio Moriya: Doesn't work
	# sxapp.setStyleSheet("QToolTip {font-size:14pt; color:white; padding:2px; border-width:2px; border-style:solid; border-radius:20px; background-color: black; border: 1px solid white;}");
	sxapp.setStyleSheet("QToolTip {font-size:%dpt;}" % (new_point_size));

	# Initialise a singleton class for look & feel constants
	SXLookFeelConst.initialise(sxapp)

	# Define the main window (class SXMainWindow)
	sxmain_window = SXMainWindow()
	sxmain_window.setWindowTitle("SPHIRE-GUI Main (Alpha Version)")
	sxmain_window.setMinimumWidth(SXLookFeelConst.sxmain_window_width)
	sxmain_window.setMinimumHeight(SXLookFeelConst.sxmain_window_height)
	sxmain_window.resize(SXLookFeelConst.sxmain_window_width, SXLookFeelConst.sxmain_window_height)
	sxmain_window.move(QPoint(SXLookFeelConst.sxmain_window_left, SXLookFeelConst.sxmain_window_top));

	# Show main window
	sxmain_window.show()
	sxmain_window.raise_()

	# Update qsub enable state of all sx command category widgets after window is displayed and raised
	sxmain_window.update_qsub_enable_states()

	# Start event handling loop
	sxapp.exec_()

# ========================================================================================
if __name__ == "__main__":
	main(sys.argv)

# ========================================================================================
# END OF SCRIPT
# ========================================================================================
