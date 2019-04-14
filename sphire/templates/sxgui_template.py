#!/usr/bin/env python
from __future__ import print_function
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

import os
import subprocess
import sys
import collections
from builtins import range
from builtins import object
try:
    from PyQt4.Qt import *
    from PyQt4 import QtGui
    from PyQt4 import QtCore
except ImportError:
    from PyQt5.Qt import *
    from PyQt5 import QtWidgets as QtGui
    from PyQt5 import QtCore
from EMAN2 import *
from EMAN2_cppwrap import *
from sp_global_def import *
from sp_sparx import *
from optparse import OptionParser
from functools import partial  # Use to connect event-source widget and event handler
from subprocess import *
import re

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
		for idx in range(0, len(path_tokens) - 2):
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
# Inherited by SXcmd_category, SXconst_set, and SXoperand_set
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
		self.key_base = ""            # key base name of command token (argument or option) in command line
		self.key_prefix = ""          # key prefix of of command token. None for argument, "--" or "-" for option
		self.label = ""               # User friendly name of argument or option
		self.help = ""                # Help info
		self.group = ""               # Tab group: main or advanced
		self.dependency_group = [["", "", ""]]    # Depencency group: Disables or enables widgets
		self.is_required = False      # Required argument or options. No default values are available
		self.is_locked = False        # The restore value will be used as the locked value.
		self.is_reversed = False      # Reversed default value of bool. The flag will be added if the value is same as default 
		self.default = ""             # Default value
		self.restore = [[""], [""]]             # Restore value
		self.type = ""                # Type of value
		self.filament_tab = ""                # Type of value
		# NOTE: Toshio Moriya 2018/01/19
		# self.is_in_io should be removed after cleaning up MoinMoin related codes.
		self.is_in_io = False         # <Used only in wikiparser.py> To check consistency between "usage in command line" and list in "== Input ==" and "== Output ==" sections
		self.restore_widget = None    # <Used only in sxgui.py> Restore widget instance associating with this command token
		self.widget = None            # <Used only in sxgui.py> Main widget instance associating with this command token
		self.subwidget_left = None    # <Used only in sxgui.py> Subwidget instance at the left associating with the helper utility of this command token (e.g. conversion calculator)
		self.subwidget_right = None   # <Used only in sxgui.py> SubWidget instance at the right associating with the helper utility of this command token (e.g. conversion calculator)
		self.calculator_dialog = None # <Used only in sxgui.py> Calculator dialog instance associating with the helper utility of this command token (e.g. conversion calculator)
		self.other_dialog_list = []   # <Used only in sxgui.py> List of the other calculator dialog instances associating with this command (e.g. conversion calculator)
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

	def initialize_edit(self, key_base):
		self.key_base = key_base
		self.key_prefix = None
		self.label = None
		self.help = None
		self.group = None
		self.is_required = None
		self.is_locked = None
		self.is_reversed = None
		self.default = None
		self.restore = None
		self.type = None

# ========================================================================================
class SXcmd(object):
	def __init__(self, category = "", role = "", is_submittable = True):
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.name = ""                        # Name of this command (i.e. name of sx*.py script but without .py extension), used for generating command line
		self.subname = ""                     # Subname of this command (i.e. sub-command name), used for generating command line. For fullset command, use empty string
		self.mode = ""                        # Key base name of a command option token, defining mode/subset of this command (i.e. option mode name). For fullset command, use empty string
		self.subset_config = ""               # Unique name to differentiate subset configuration of this command. For example, to name a command argument mode, which dependes on the number of input arguments. If not necessary, use empty string 
		self.label = ""                       # User friendly name of this command
		self.short_info = ""                  # Short description of this command
		self.mpi_support = False              # Flag to indicate if this command suppors MPI version
		self.mpi_add_flag = False             # DESIGN_NOTE: 2015/11/12 Toshio Moriya. This can be removed when --MPI flag is removed from all sx*.py scripts
		self.category = category              # Category of this command: sxc_movie, sxc_cter, sxc_window, sxc_isac, sxc_viper, sxc_meridien, sxc_sort3d, sxc_localres, sxc_utilities
		self.role = role                      # Role of this command; sxr_pipe (pipeline), sxr_alt (alternative) sxr_util (utility)
		self.is_submittable = is_submittable  # External GUI Application (e.g. sxgui_cter.py) should not be submitted to job queue
		self.token_list = []                  # List of command tokens. Need this to keep the order of command tokens
		self.token_dict = {}                  # Dictionary of command tokens, organised by key base name of command token. Easy to access a command token but looses their order
		self.dependency_dict = collections.OrderedDict()             # Dictionary of command tokens, containing the dependencies
		self.btn = None                       # <Used only in sxgui.py> QPushButton button instance associating with this command
		self.widget = None                    # <Used only in sxgui.py> SXCmdWidget instance associating with this command
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

	def get_mode_name_for(self, target_name):
		mode_name = self.name
		
		if self.subname != "":
			if target_name in ["file_path"]:
				mode_name = "%s_%s" % (mode_name, self.subname)
			elif target_name in ["human"]:
				mode_name = "%s %s" % (mode_name, self.subname)
		
		if self.mode != "":
			if target_name in ["file_path"]:
				mode_name = "%s_%s" % (mode_name, self.mode)
			elif target_name in ["human"]:
				mode_name = "%s %s%s" % (mode_name, self.token_dict[self.mode].key_prefix, self.mode)
		
		if self.subset_config != "":
			if target_name in ["file_path"]:
				mode_name = "%s_%s" % (mode_name, self.subset_config.replace(" ", "_"))
			elif target_name in ["human"]:
				mode_name = "%s (%s)" % (mode_name, self.subset_config)
		
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
		# self.name = name              # <Inherit from SXmenu_item> Name of this command category (i.e. sxc_movie, sxc_cter, sxc_window, sxc_isac, sxc_viper, sxc_meridien, sxc_sort3d, sxc_localres, sxc_utilities), used as a key of dictionary
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
		super(SXmenu_item, self).__init__()
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
class SXoperand(object):
	def __init__(self):
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.key = ""                  # <Used only in sxgui.py> key of calculator operand
		self.label = ""                # <Used only in sxgui.py> User friendly name of calculator operand
		self.help = ""                 # <Used only in sxgui.py> Help info
		self.register = ""             # <Used only in sxgui.py> Default value
		self.type = ""                 # <Used only in sxgui.py> Type of value
		self.is_input = True           # <Used only in sxgui.py> Flag to indicate if this operand is input or output
		self.validated_register = None # <Used only in sxgui.py> Contain validated numerical value of register only if the value is valid. If not, it is None. This is used as a flag to indicate if the registered value is valid or not
		self.validated = None          # <Used only in sxgui.py> Contain Validated numerical value of command widget only if the value is valid. If not, it is None. This is used as a flag to indicate if the registered value is valid or not
		self.register_widget = None    # <Used only in sxgui.py> Restore widget instance associating with this calculator operand
		self.widget = None             # <Used only in sxgui.py> Widget instance associating with this calculator operand
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

# ========================================================================================
class SXoperand_set(SXmenu_item):
	def __init__(self):
		super(SXmenu_item, self).__init__()
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		# self.name = ""        # <Inherit from SXmenu_item> Name of this calculator operand set
		# self.label = ""       # <Inherit from SXmenu_item> User friendly name of this set
		# self.short_info = ""  # <Inherit from SXmenu_item> Short description of this set
		self.list = []          # <Used only in sxgui.py> list of calculator operands. Need this to keep the order of calculator operands
		self.dict = {}          # <Used only in sxgui.py> dictionary of calculator operands, organised by keys of calculator operands. Easy to access each calculator operand but looses their order
		# self.btn = None       # <Inherit from SXmenu_item> QPushButton button instance associating with this set
		# self.widget = None    # <Inherit from SXmenu_item> Widget instance associating with this set
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

# ========================================================================================
class SXLookFeelConst(object):
	# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	# static class variables
	# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	default_bg_color = QColor(229, 229, 229, 192) # default_bg_color = QColor(229, 229, 229, 242) # Greyish-White Transparent
	default_bg_color_string = 'rgba(229, 229, 229, 192)' # default_bg_color = QColor(229, 229, 229, 242) # Greyish-White Transparent
	default_bg_color_calculator = QColor(214, 214, 214) # Greyish-White
	sxinfo_widget_bg_color = QColor(0, 0, 0, 10) # Almost-Completely Transparent
	sxcmd_widget_bg_color = QColor(0, 0, 0, 0) # Completely Transparent
	sxcmd_tab_bg_color = QColor(229, 229, 229, 200) # White Transparent
	sxcmd_tab_bg_color_string = 'rgba(229, 229, 229, 200)' # White Transparent

	# Constants
	project_dir_raw = "sxgui_settings"
	project_dir = project_dir_raw
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
	def initialise(sxapp, version):
		# Set the directory for all file dialogs to script directory
		SXLookFeelConst.project_dir += '_{0}'.format(version)
		SXLookFeelConst.file_dialog_dir = os.getcwd()
		
		# Check whether current settings directory exists
		if os.path.exists(SXLookFeelConst.project_dir):
			current_settings_exist = True
		else:
			current_settings_exist = False
			
		do_longer_warning = False
		older_settings_exist = False
		for file_name in os.listdir(SXLookFeelConst.file_dialog_dir):
			if SXLookFeelConst.project_dir_raw in file_name and file_name != SXLookFeelConst.project_dir:
				print('WARNING: Directory %s/ detected that belongs to another version of SPHIRE.' % file_name)
				older_settings_exist = True
				
				if not current_settings_exist:
					do_longer_warning = True
				
		if do_longer_warning:
			print('To load old settings, please load the gui settings manually.')
			print('Backwards compatibility cannot be guaranteed.')
		
		if not current_settings_exist and not older_settings_exist:
			print("\nSettings directory for current SPHIRE version %s doesn't exist\nWould you like to create a new project directory and continue?\nYou need to run the sphire command in the foreground (without &) to answer this question\n[y/n] "% SXLookFeelConst.project_dir)
			answer = raw_input()
			while answer.lower() not in ('y', 'n'):
				print('Answer needs to be y or n')
				answer = raw_input()
			if answer.lower() == 'n':
				print("\nbye bye")
				exit()
		
		monitor_index = 0
		# Search for maximum screen height and set it to SXLookFeelConst singleton class
		max_screen_height = sxapp.desktop().screenGeometry().height()
		for index in range(sxapp.desktop().screenCount()):
			screen_height = sxapp.desktop().screenGeometry(index).height()
			if max_screen_height < screen_height:
				monitor_index = index
				max_screen_height = screen_height
		SXLookFeelConst.screen_height = max_screen_height
		# Search for maximum screen width and set it to SXLookFeelConst singleton class
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
		formatted_path = os.path.relpath(str(path))
		if formatted_path[:len("../")] == "../":
			# if the path is above the project root directory (current directory)
			# use absolute path
			formatted_path = path
		# else:
			# if the path is project subdirectory
			# use relative path

		return formatted_path

	@staticmethod
	def generate_sxcmd_wiki_url(sxcmd, wiki_type = "SPHIRE"):
		if wiki_type == "SPHIRE":
			# First, handle exceptional cases
			if sxcmd.name in ["e2display", "sxpdb2em", "sxrelion2sphire", "sxprocess", "e2proc3d", "sxheader", "e2bdb"]:
				sxcmd_category_name = "utilities"
			elif sxcmd.name in ["sxpipe"] and sxcmd.subname in ["resample_micrographs", "organize_micrographs", "restacking", "moon_eliminator"]:
				sxcmd_category_name = "utilities"
			# elif sxcmd.name in ["sxpipe"] and sxcmd.subname in ["reboxing"]:
			# 	sxcmd_category_name = "meridien"
			elif sxcmd.name in ["sxmeridien", "sxmeridien_20171120"]:
				sxcmd_category_name = "meridien"
			elif sxcmd.name in ["sxunblur", "sxsummovie"]:
				sxcmd_category_name = "movie"
			else:
				sxcmd_category_name = sxcmd.category.replace("sxc_", "")
			# URL Format: "http://sphire.mpg.de/wiki/doku.php?id=pipeline:CMD_CATEGORY:CMD_BASE
			sxcmd_wiki_url = "http://sphire.mpg.de/wiki/doku.php?id=pipeline:%s:%s" % (sxcmd_category_name, sxcmd.name)
			if sxcmd.subname != "":
				sxcmd_wiki_url = "%s_%s" % (sxcmd_wiki_url, sxcmd.subname)
		else:
			assert (wiki_type == "SPARX")
			sxcmd_wiki_url = "%s%s" % (SPARX_DOCUMENTATION_WEBSITE, sxcmd.name)
			if sxcmd.subname != "":
				sxcmd_wiki_url = "%s_%s" % (sxcmd_wiki_url, sxcmd.subname)
		
		return sxcmd_wiki_url

	@staticmethod
	def generate_sxmenu_item_wiki_url(sxmenu_item):
		# First, handle exceptional cases
		sxmenu_item_wiki_url = None
		if sxmenu_item.name in ["sxc_gui_info"]:
			sxmenu_item_wiki_url = "http://sphire.mpg.de/wiki/doku.php?id=start"
		else:
			# URL Format: "http://sphire.mpg.de/wiki/doku.php?id=pipeline:CMD_CATEGORY:start"
			sxmenu_item_wiki_url = "http://sphire.mpg.de/wiki/doku.php?id=pipeline:%s:start" % (sxmenu_item.name.replace("sxc_", ""))
		assert (sxmenu_item_wiki_url is not None)
		
		return sxmenu_item_wiki_url

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
		self.setToolTip('HELP')

# ========================================================================================
class SXPictogramButton(QPushButton):
	def __init__(self, pictogram_name, pictogram_file_path, parent = None):
		super(SXPictogramButton, self).__init__(parent)

		# print "MRK_DEBUG: pictogram_file_path = %s" % pictogram_file_path
		# print "MRK_DEBUG: os.path.exists(logo_file_path) %s" % os.path.exists(pictogram_file_path)

		# Width of pictogram image
		self.parent = parent
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
		self.setToolTip(pictogram_name.upper())

class SXMenuItemBtnAreaWidget(QWidget):
	def __init__(self, sxconst_set, sxcmd_category_list, sxinfo, parent = None):
		super(SXMenuItemBtnAreaWidget, self).__init__(parent)

		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

		# Create widgets for pipeline command category button area and miscellaneous function button area
		sxcmd_category_btn_subarea_widget = self.create_sxmenu_item_btn_subarea_widget()
		misc_func_btn_subarea_widget = self.create_sxmenu_item_btn_subarea_widget()
		self.add_sxmenu_item_btn_widget(sxconst_set, sxcmd_category_btn_subarea_widget)
		self.sxconst_set = sxconst_set
		self.sxcmd_category_list = sxcmd_category_list
		for sxcmd_category in sxcmd_category_list:
			if sxcmd_category.name != "sxc_utilities" and sxcmd_category.name != "sxc_movie":
				self.add_sxmenu_item_btn_widget(sxcmd_category, sxcmd_category_btn_subarea_widget)
			else: # assert(sxcmd_category.name == "sxc_utilities")
				self.add_sxmenu_item_btn_widget(sxcmd_category, misc_func_btn_subarea_widget)

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
	### process_started = pyqtSignal()

	def __init__(self, sxconst_set, sxcmd, helical, parent = None):
		super(SXCmdWidget, self).__init__(parent)

		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.sxconst_set = sxconst_set
		self.sxcmd = sxcmd

		self.sxcmd_tab_main = None
		self.sxcmd_tab_advance = None

		self.child_application_list = []

		self.gui_settings_file_path = "%s/gui_settings_%s_%d.txt" % (self.sxcmd.get_category_dir_path(SXLookFeelConst.project_dir), self.sxcmd.get_mode_name_for("file_path"), int(helical))

		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

		# Set grid layout
		grid_layout = QGridLayout(self)
		# grid_layout.setContentsMargins(SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin)
		# grid_layout.setSpacing(SXLookFeelConst.grid_spacing)

		self.setAutoFillBackground(True)
		palette = QPalette()
		palette.setBrush(QPalette.Background, QBrush(SXLookFeelConst.sxcmd_widget_bg_color))
		self.setPalette(palette)

		self.sxcmd_tab_main = SXCmdTab("Main", helical, self)
		self.sxcmd_tab_advance = SXCmdTab("Advanced", helical, self)
		for key in self.sxcmd.dependency_dict:
			if key not in self.sxcmd.token_dict:
				continue
			self.sxcmd.token_dict[key].widget.my_changed_state.connect(self.change_state)
			if isinstance(self.sxcmd.token_dict[key].widget, SXCheckBox):
				idx = self.sxcmd.token_dict[key].widget.checkState()
				if idx == 0:
					self.sxcmd.token_dict[key].widget.setCheckState(2)
				elif idx == 1:
					self.sxcmd.token_dict[key].widget.setCheckState(2)
				elif idx == 2:
					self.sxcmd.token_dict[key].widget.setCheckState(0)
				self.sxcmd.token_dict[key].widget.setCheckState(idx)
			elif isinstance(self.sxcmd.token_dict[key].widget, SXLineEdit):
				text = self.sxcmd.token_dict[key].widget.text()
				self.sxcmd.token_dict[key].widget.setText('0000000000000000099999')
				self.sxcmd.token_dict[key].widget.setText(text)


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

	@QtCore.pyqtSlot(str, str)
	@QtCore.pyqtSlot(str, int)
	def change_state(self, dependency, state, old_dependency=None, sender=None, prev_sender=None):
		"""
		Change the state of widgets based on the choice of the corresponding combo box

		name - Name of the group to change status (Emitted by the combo box)

		Returns:
		None
		"""

		#print(dependency, state, old_dependency)
		if sender is None:
			sender = self.sender()
		else:
			sender = sender

		dependency = str(dependency)
		if old_dependency is None:
			old_dependency = [dependency]
		else:
			old_dependency.append(dependency)
		parent = sender.parent()
		while True:
			if hasattr(parent, 'sxcmd'):
				break
			else:
				parent = parent.parent()

		final_state = True
		for name, exp_state, inverse in parent.sxcmd.token_dict[dependency].dependency_group:
			if not name:
				continue

			try:
				widget = parent.sxcmd.token_dict[name].widget
			except KeyError:
				continue

			try:
				curr_state = widget.checkState()
			except AttributeError:
				curr_state = widget.text()

			if isinstance(widget, SXCheckBox):
				type_check = bool
			elif isinstance(widget, SXLineEdit):
				type_check = str
			else:
				assert False

			if str(type_check(curr_state)) != exp_state and inverse == 'True':
				is_enabled = True
			elif str(type_check(curr_state)) == exp_state and inverse == 'True':
				is_enabled = False
			elif str(type_check(curr_state)) == exp_state:
				is_enabled = True
			else:
				is_enabled = False

			if len(parent.sxcmd.token_dict[dependency].dependency_group) == 2 and is_enabled:
				if prev_sender is not None:
					is_enabled = prev_sender.isEnabled()

			if not is_enabled:
				final_state = False
		sender.setEnabled(final_state)

		try:
			chain_deps = parent.sxcmd.dependency_dict[dependency]
		except KeyError:
			pass
		else:
			for name, exp_state, inverse in chain_deps:
				if not name:
					continue

				try:
					widget = parent.sxcmd.token_dict[name].widget
				except KeyError:
					continue

				if not name in old_dependency:
					self.change_state(name, exp_state, old_dependency=old_dependency, sender=widget, prev_sender=sender)

		#try:
		#	for name, exp_state, inverse in parent.sxcmd.dependency_dict[str(dependency)]:
		#		if not parent.sxcmd.token_dict[str(dependency)].widget.isEnabled():
		#			is_enabled = False
		#		elif str(type_check(state)) != exp_state and inverse == 'True':
		#			is_enabled = True
		#		elif str(type_check(state)) == exp_state and inverse == 'True':
		#			is_enabled = False
		#		elif str(type_check(state)) == exp_state:
		#			is_enabled = True
		#		else:
		#			is_enabled = False
		#		print('STATE', name, exp_state, state, is_enabled, type_check, sender)
		#		SXCmdTab.set_text_entry_widget_enable_state(parent.sxcmd.token_dict[name].widget, is_enabled)
		#	child_dict = {}
		#	for name, exp_state, inverse in parent.sxcmd.dependency_dict[str(dependency)]:
		#		for entry in parent.sxcmd.token_dict[name].dependency_group:
		#			if not entry[0]:
		#				continue
		#			if name not in old_dependency and str(dependency) != entry[0]:
		#				try:
		#					new_state = parent.sxcmd.token_dict[name].widget.checkState()
		#				except AttributeError:
		#					new_state = parent.sxcmd.token_dict[name].widget.text()
		#				child_dict[name] = new_state
		#	for key, value in child_dict.items():
		#		print(key, value)
		#		self.change_state(key, value, old_dependency, parent.sxcmd.token_dict[key].widget)
		#except KeyError:
		#	return None
		#print(state)
		#try:
		#	for entry in self.cmd[name]:
		#		widget = entry[0]
		#		state = entry[1]
		#		sub_name = entry[2]
		#		if not self.content[name].isEnabled():
		#			widget.setEnabled(False)
		#		elif self.content[name].edit.currentText() == state:
		#			widget.setEnabled(True)
		#		else:
		#			widget.setEnabled(False)
		#		self.change_state(name=sub_name)
		#except KeyError:
		#	return None


	def map_widgets_to_sxcmd_line(self):
		# Add program name to command line
		sxcmd_line = "%s.py" % (self.sxcmd.name)
		
		if self.sxcmd.subname != "":
			sxcmd_line += " %s" % (self.sxcmd.subname)

		# Loop through all command tokens
		for sxcmd_token in self.sxcmd.token_list:
			# First, handle very special cases
			if not isinstance(sxcmd_token.widget, list):
				widgets = [sxcmd_token.widget]
			else:
				widgets = sxcmd_token.widget

			do_continue = False
			if not sxcmd_token.is_locked:
				for widget in widgets:
					if not widget.isEnabled():
						do_continue = True
			if do_continue:
				continue
			elif sxcmd_token.type == "user_func":
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
				if sxcmd_token.type == "bool_ignore":
					### Possbile cases:
					### if not sxcmd_token.is_reversed  and  (sxcmd_token.widget.checkState() == Qt.Checked) != sxcmd_token.default  and  sxcmd_token.is_required == True:  # Add this token to command line
					### if not sxcmd_token.is_reversed  and  (sxcmd_token.widget.checkState() == Qt.Checked) == sxcmd_token.default  and  sxcmd_token.is_required == True:  # Add this token to command line
					### if not sxcmd_token.is_reversed  and  (sxcmd_token.widget.checkState() == Qt.Checked) != sxcmd_token.default  and  sxcmd_token.is_required == False: # Add this token to command line
					### if not sxcmd_token.is_reversed  and  (sxcmd_token.widget.checkState() == Qt.Checked) == sxcmd_token.default  and  sxcmd_token.is_required == False: # Do not add this token to command line
					### 
					### if     sxcmd_token.is_reversed  and  (sxcmd_token.widget.checkState() == Qt.Checked) != sxcmd_token.default  and  sxcmd_token.is_required == True:  # Add this token to command line
					### if     sxcmd_token.is_reversed  and  (sxcmd_token.widget.checkState() == Qt.Checked) == sxcmd_token.default  and  sxcmd_token.is_required == True:  # Add this token to command line
					### if     sxcmd_token.is_reversed  and  (sxcmd_token.widget.checkState() == Qt.Checked) != sxcmd_token.default  and  sxcmd_token.is_required == False: # Do not add this token to command line
					### if     sxcmd_token.is_reversed  and  (sxcmd_token.widget.checkState() == Qt.Checked) == sxcmd_token.default  and  sxcmd_token.is_required == False: # Add this token to command line
					### 
					#else: # Do not add this token to command line
					continue
					
				elif sxcmd_token.type == "bool":
					### Possbile cases:
					### if not sxcmd_token.is_reversed  and  (sxcmd_token.widget.checkState() == Qt.Checked) != sxcmd_token.default  and  sxcmd_token.is_required == True:  # Add this token to command line
					### if not sxcmd_token.is_reversed  and  (sxcmd_token.widget.checkState() == Qt.Checked) == sxcmd_token.default  and  sxcmd_token.is_required == True:  # Add this token to command line
					### if not sxcmd_token.is_reversed  and  (sxcmd_token.widget.checkState() == Qt.Checked) != sxcmd_token.default  and  sxcmd_token.is_required == False: # Add this token to command line
					### if not sxcmd_token.is_reversed  and  (sxcmd_token.widget.checkState() == Qt.Checked) == sxcmd_token.default  and  sxcmd_token.is_required == False: # Do not add this token to command line
					### 
					### if     sxcmd_token.is_reversed  and  (sxcmd_token.widget.checkState() == Qt.Checked) != sxcmd_token.default  and  sxcmd_token.is_required == True:  # Add this token to command line
					### if     sxcmd_token.is_reversed  and  (sxcmd_token.widget.checkState() == Qt.Checked) == sxcmd_token.default  and  sxcmd_token.is_required == True:  # Add this token to command line
					### if     sxcmd_token.is_reversed  and  (sxcmd_token.widget.checkState() == Qt.Checked) != sxcmd_token.default  and  sxcmd_token.is_required == False: # Do not add this token to command line
					### if     sxcmd_token.is_reversed  and  (sxcmd_token.widget.checkState() == Qt.Checked) == sxcmd_token.default  and  sxcmd_token.is_required == False: # Add this token to command line
					### 
					is_flag_required = False
					if sxcmd_token.is_required:
						is_flag_required = True
					else:
						assert (not sxcmd_token.is_required)
						if not sxcmd_token.is_reversed and (sxcmd_token.widget.checkState() == Qt.Checked) != sxcmd_token.default:
							is_flag_required = True
						elif sxcmd_token.is_reversed and (sxcmd_token.widget.checkState() == Qt.Checked) == sxcmd_token.default:
							is_flag_required = True
						else:
							assert (not is_flag_required)
					
					if is_flag_required:
						sxcmd_line += " %s%s" % (sxcmd_token.key_prefix, sxcmd_token.key_base)
					#else: # Do not add this token to command line
					
				else:
					if sxcmd_token.widget.text() == sxcmd_token.default:
						### if sxcmd_token.widget.text() == sxcmd_token.default and sxcmd_token.is_required == True:  # Error case
						if sxcmd_token.is_required == True:
							QMessageBox.warning(self, "Invalid parameter value", "Token (%s) of command (%s) is required. Please set the value for this." % (sxcmd_token.label, self.sxcmd.get_mode_name_for("human")))
							return ""
						### if sxcmd_token.widget.text() == sxcmd_token.default and sxcmd_token.is_required == False: # Do not add this token to command line
						# else: # assert(sxcmd_token.is_required == False) # Do not add to this command line
					else: # sxcmd_token.widget.text() != sxcmd_token.default
						### if sxcmd_token.widget.text() != sxcmd_token.default and sxcmd_token.is_required == True:  # Add this token to command line
						### if sxcmd_token.widget.text() != sxcmd_token.default and sxcmd_token.is_required == False: # Add this token to command line

						# For now, using line edit box for the other type
						widget_text = str(sxcmd_token.widget.text())
###						if sxcmd_token.type not in ["int", "float", "apix", "ctfwin", "box", "radius", "mass", "displayable_list", "mic_one_list", "any_file_list", "any_image_list", "dir_list"]:
						if sxcmd_token.type not in ["int", "float", "abs_freq", "apix", "ctfwin", "box", "radius", "mass", "displayable_list", "mic_one_list", "dir_list", "filament_width"]:
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

	@QtCore.pyqtSlot()
	def execute_cmd_line(self, execute=True, output_dir='.', number=None):
		# Disable the run command button
		execute_btn = self.sender()
		execute_btn.setEnabled(False)
		QtCore.QTimer.singleShot(5000, lambda: execute_btn.setEnabled(True))

		# Generate command line
		cmd_line = self.generate_cmd_line()

		if cmd_line:
			# Command line is not empty
			# First, check existence of outputs
			for sxcmd_token in self.sxcmd.token_list:
				if sxcmd_token.type == "output" or sxcmd_token.type == "output_continue" or sxcmd_token.type == "output_bdb2d_stack":
					if os.path.exists(sxcmd_token.widget.text()) or db_check_dict(str(sxcmd_token.widget.text())):
						# DESIGN_NOTE: 2015/11/24 Toshio Moriya
						# This special case needs to be handled with more general method...
						if sxcmd_token.type == "output_continue":
							reply = QMessageBox.question(self, "Output Directory/File", "Output Directory/File (%s) already exists. Do you really want to run the program with continue mode?" % (sxcmd_token.widget.text()), QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
							if reply == QMessageBox.No:
								return
							# else: # Do nothing
						else:
							assert(sxcmd_token.type == "output" or sxcmd_token.type == "output_bdb2d_stack")
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
				if number is not None:
					file_name_qsub_script = os.path.join(output_dir, "{0:04d}_".format(number) + "qsub_" + str(self.sxcmd_tab_main.qsub_job_name_edit.text()) + ".sh")
				else:
					file_name_qsub_script = os.path.join(output_dir, "qsub_" + str(self.sxcmd_tab_main.qsub_job_name_edit.text()) + ".sh")
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
				print("Wrote the following command line in the queue submission script: ")
				print(cmd_line_in_script)
				if execute:
					print("Submitted a job by the following command: ")
					print(cmd_line)
				else:
					print("Saved command to: ")
					print(file_name_qsub_script)
			elif self.sxcmd_tab_main.qsub_enable_checkbox.checkState() == Qt.Unchecked and not execute:
				QMessageBox.warning(self, "Qsub template needs to be specified for pipeline option", "Qsub template needs to be specified for pipeline option")
				return
			else:
				# Case 2: queue submission is disabled (MPI can be supported or unsupported)
				if self.sxcmd_tab_main.qsub_enable_checkbox.checkState() == Qt.Checked: ERROR("Logical Error: Encountered unexpected condition for sxcmd_tab_main.qsub_enable_checkbox.checkState. Consult with the developer.", "%s in %s" % (__name__, os.path.basename(__file__)))
				print("Executed the following command: ")
				print(cmd_line)

			# Execute the generated command line
			if execute:
				process = subprocess.Popen(cmd_line, shell=True)
				### self.process_started.emit(process.pid)
				if self.sxcmd.is_submittable == False:
					assert(self.sxcmd.mpi_support == False)
					# Register to This is a GUI application
					self.child_application_list.append(process)

			# Save the current state of GUI settings
			if os.path.exists(self.sxcmd.get_category_dir_path(SXLookFeelConst.project_dir)) == False:
				os.mkdir(self.sxcmd.get_category_dir_path(SXLookFeelConst.project_dir))
			self.write_params(self.gui_settings_file_path)
		# else: SX command line is be empty because an error happens in generate_cmd_line. Let's do nothing

	def add_to_pipeline(self):
		name = QtGui.QFileDialog.getExistingDirectory(self, 'Pipeline output directory', SXLookFeelConst.file_dialog_dir, options = QFileDialog.DontUseNativeDialog)
		if isinstance(name, tuple):
			dir_name = str(name[0])
		else:
			dir_name = str(name)
		if dir_name:
			number_match = re.compile('(\d{4})_')
			number = -1
			for file_name in sorted(os.listdir(dir_name)):
				try:
					number = number_match.match(file_name).group(1)
				except AttributeError:
					QMessageBox.warning(self, "Invalid file in pipeline directory", "There is a file in the pipeline directory not starting with four digits. Abort...")
					return

			self.execute_cmd_line(execute=False, output_dir=dir_name, number=int(number)+1)

	def print_cmd_line(self):
		# Generate command line
		cmd_line = self.generate_cmd_line()
		if cmd_line:
			message_line = "Generated the following command line:"
			print(message_line)
			print(cmd_line)
			QtGui.QMessageBox.information(self, "Information","%s \n\n%s" % (message_line, cmd_line))

			# Save the current state of GUI settings
			if os.path.exists(self.sxcmd.get_category_dir_path(SXLookFeelConst.project_dir)) == False:
				os.mkdir(self.sxcmd.get_category_dir_path(SXLookFeelConst.project_dir))
			self.write_params(self.gui_settings_file_path)
		# else: Do nothing

	def show_dialog_calculated_res(self):
		# Generate command line
		calculate = FSCPlot(self)
		fsc_plot.plot(*args)

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
					if cmd_token.type == "user_func":
						# This type has two line edit boxes as a list of widget
						n_widgets = 2
						for widget_index in range(n_widgets):
							val_str = str(cmd_token.widget[widget_index].text())
							file_out.write("<%s> %s (default %s) == %s \n" % (cmd_token.key_base, cmd_token.label[widget_index], cmd_token.default[widget_index], val_str))
					# Then, handle the other cases
					else:
						val_str = ""
						if cmd_token.type in ("bool", "bool_ignore"):
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
			cmd_token_apix = None
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
						QMessageBox.warning(self, "Invalid Parameter File Format", "Command token entry should start from \"%s\" for key base name in line (%s) of file (%s). The format of this file might be corrupted. Please save the parameter file again." % (target_operator, line_in, file_path_in))
					label_in = label_in[item_tail + len(target_operator):].strip() # Get the rest of line
					target_operator = ">"
					item_tail = label_in.find(target_operator)
					if item_tail == -1:
						QMessageBox.warning(self, "Invalid Parameter File Format", "Command token entry should have \"%s\" closing key base name in line (%s) of file (%s). The format of this file might be corrupted. Please save the parameter file again." % (target_operator, line_in, file_path_in))
					key_base = label_in[0:item_tail]
					# Get corresponding cmd_token
					if key_base not in list(self.sxcmd.token_dict.keys()):
						QMessageBox.warning(self, "Invalid Parameter File Format", "Invalid base name of command token \"%s\" is found in line (%s) of file (%s). This parameter file might be incompatible with the current version. Please save the parameter file again." % (key_base, line_in, file_path_in))
					cmd_token = self.sxcmd.token_dict[key_base]
					if not cmd_token.is_locked: 
						# First, handle very special cases
						if cmd_token.type == "user_func":
							cmd_token.widget[function_type_line_counter].setText(val_str_in)
							function_type_line_counter += 1
							function_type_line_counter %= n_function_type_lines # function have two line edit boxes
						else:
							if cmd_token.type in ("bool", "bool_ignore"):
								# construct new widget(s) for this command token
								if val_str_in == "YES":
									cmd_token.widget.setChecked(Qt.Checked)
								else: # val_str_in == "NO"
									cmd_token.widget.setChecked(Qt.Unchecked)
							# Then, handle the other cases
							else:
								# For now, use line edit box for the other type
								cmd_token.widget.setText(val_str_in)
								if cmd_token.type == "apix":
									cmd_token_apix = cmd_token 
			if cmd_token_apix is not None:
				assert (cmd_token_apix.type == "apix")
				# if len(cmd_token_apix.other_dialog_list) > 0:
				# 	print("MRK_DEBUG: ")
				# 	print("MRK_DEBUG: ----- SXCmdWidget.read_params() ----- ")
				# 	print("MRK_DEBUG: cmd_token_apix.widget.text() := \"{}\"".format(cmd_token_apix.widget.text()))
				# 	print("MRK_DEBUG: len(cmd_token_apix.other_dialog_list) := \"{}\"".format(len(cmd_token_apix.other_dialog_list)))
				for sxcmd_token_apix_other_dialog in cmd_token_apix.other_dialog_list:
				# 	print("MRK_DEBUG: BEFORE sxcmd_token_apix_other_dialog.sxconst_register_widget_apix.text() := \"{}\"".format(sxcmd_token_apix_other_dialog.sxconst_register_widget_apix.text()))
				# 	print("MRK_DEBUG: BEFORE sxcmd_token_apix_other_dialog.sxcmd_token_widget_apix.text() := \"{}\"".format(sxcmd_token_apix_other_dialog.sxcmd_token_widget_apix.text()))
				# 	print("MRK_DEBUG: BEFORE sxcmd_token_apix_other_dialog.sxcmd_token_widget_abs_freq.text() := \"{}\"".format(sxcmd_token_apix_other_dialog.sxcmd_token_widget_abs_freq.text()))
					sxcmd_token_apix_other_dialog.reflect_external_local_update_apix_and_abs_freq()
				# 	print("MRK_DEBUG: AFTER sxcmd_token_apix_other_dialog.sxconst_register_widget_apix.text() := \"{}\"".format(sxcmd_token_apix_other_dialog.sxconst_register_widget_apix.text()))
				# 	print("MRK_DEBUG: AFTER sxcmd_token_apix_other_dialog.sxcmd_token_widget_apix.text() := \"{}\"".format(sxcmd_token_apix_other_dialog.sxcmd_token_widget_apix.text()))
				# 	print("MRK_DEBUG: AFTER sxcmd_token_apix_other_dialog.sxcmd_token_widget_abs_freq.text() := \"{}\"".format(sxcmd_token_apix_other_dialog.sxcmd_token_widget_abs_freq.text()))
		else:
			QMessageBox.warning(self, "Fail to load parameters", "The specified file parameter file for %s could not be read." % self.sxcmd.get_mode_name_for("human"))

		file_in.close()

	def save_params(self):
		name = QFileDialog.getSaveFileName(self, "Save Parameters", SXLookFeelConst.file_dialog_dir, options = QFileDialog.DontUseNativeDialog)
		if isinstance(name, tuple):
			file_path = str(name[0])
		else:
			file_path = str(name)
		if file_path != "":
			self.write_params(file_path)

	def load_params(self):
		name = QFileDialog.getOpenFileName(self, "Load parameters", SXLookFeelConst.file_dialog_dir, options = QFileDialog.DontUseNativeDialog)
		if isinstance(name, tuple):
			file_path = str(name[0])
		else:
			file_path = str(name)
		if file_path != "":
			self.read_params(file_path)
			self.sxcmd_tab_main.set_qsub_enable_state()

	def select_file(self, target_widget, file_format = ""):
		file_path = ""
		# NOTE: Toshio Moriya 2018/01/25
		# All supported image/volume formats according to http://blake.bcm.edu/emanwiki/EMAN2ImageFormats
		# ;; HDF (*.hdf);; MRC (*.mrc);; MRCS (*.mrcs);; Spider (*.spi);; Imagic (*.img *hed);; TIFF (*.tif *.tiff);; PNG (*.png);; JPEG (*.jpg *.jpeg);; BDB (*.bdb);; Gatan (*.dm2 *.dm3 *.dm4);; FEI (*.ser);; EM (*.em);; ICOS (*.icos);; Amira (*.am);; DF3 (*.d3);; FITS (*.fts);; LST (*.lst);; LSTFAST (*.lsx *.lst);; OMAP (*.omap);; PGM (*.pgm);; PIF (*.pif);; SAL (*.hdr *.img);; SITUS (*.situs);; V4L (*.v4l);; VTK (*.vtk);; XPLOR (*.xplor)
		# 
		if file_format == "displayable_list":
			# Read not supported: ;; JPEG (*.jpg *.jpeg)
			name = QFileDialog.getOpenFileNames(self, "Select any displayable files", SXLookFeelConst.file_dialog_dir, "Typical displayable files (*.hdf *.bdb *.mrc *.mrcs *.spi *.img *.tif *.tiff *.png *.txt);; HDF (*.hdf);; BDB (*.bdb);; MRC (*.mrc);; MRCS (*.mrcs);; Spider (*.spi);; Imagic (*.img *hed);; TIFF (*.tif *.tiff);; PNG (*.png);; Text (*.txt);; Gatan (*.dm2 *.dm3 *.dm4);; FEI (*.ser);; EM (*.em);; ICOS (*.icos);; Amira (*.am);; DF3 (*.d3);; FITS (*.fts);; LST (*.lst);; LSTFAST (*.lsx *.lst);; OMAP (*.omap);; PGM (*.pgm);; PIF (*.pif);; SAL (*.hdr *.img);; SITUS (*.situs);; V4L (*.v4l);; VTK (*.vtk);; XPLOR (*.xplor);; All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path_list = name[0]
			else:
				file_path_list = name
			for a_file_path in file_path_list:
				# Use relative path.
				a_file_path = str(SXLookFeelConst.format_path(a_file_path))
				try: # Check if the path is bdb
					a_file_path = translate_to_bdb_path(a_file_path) # Convert the standard path to bdb key if possible.
				except ValueError:  # If the path is not bdb, we will receive this exception
					pass # This is not bdb path. Then, use standard path
				file_path += a_file_path + " "
		elif file_format == "data2d3d_both":
			# Read not supported: ;; JPEG (*.jpg *.jpeg)
			name = QFileDialog.getOpenFileName(self, "Select any image/volume file", SXLookFeelConst.file_dialog_dir, "Typical image & volume files (*.hdf *.bdb *.mrc *.mrcs *.spi *.img *.tif *.tiff *.png);; HDF (*.hdf);; BDB (*.bdb);; MRC (*.mrc);; MRCS (*.mrcs);; Spider (*.spi);; Imagic (*.img *hed);; TIFF (*.tif *.tiff);; PNG (*.png);; Gatan (*.dm2 *.dm3 *.dm4);; FEI (*.ser);; EM (*.em);; ICOS (*.icos);; Amira (*.am);; DF3 (*.d3);; FITS (*.fts);; LST (*.lst);; LSTFAST (*.lsx *.lst);; OMAP (*.omap);; PGM (*.pgm);; PIF (*.pif);; SAL (*.hdr *.img);; SITUS (*.situs);; V4L (*.v4l);; VTK (*.vtk);; XPLOR (*.xplor);; All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
				try: # Check if the path is bdb
					file_path = translate_to_bdb_path(file_path) # Convert the standard path to bdb key if possible.
				except ValueError:  # If the path is not bdb, we will receive this exception
					pass # This is not bdb path. Then, use standard path
		elif file_format == "mrc2d_mic_both":
			name = QFileDialog.getOpenFileName(self, "Select MRC micrograph/movie file", SXLookFeelConst.file_dialog_dir, "MRC micrograph & movie files (*.mrc *.mrcs);; MRC (*.mrc);; MRCS (*.mrcs)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "mic_both":
			# Read not supported: ;; JPEG (*.jpg *.jpeg)
			name = QFileDialog.getOpenFileName(self, "Select any micrograph/movie file", SXLookFeelConst.file_dialog_dir, "Typical micrograph & movie files (*.mrc *.mrcs *.tif *.tiff *.hdf *.bdb *.spi *.img);; MRC (*.mrc);; MRCS (*.mrcs);; TIFF (*.tif *.tiff);; HDF (*.hdf);; BDB (*.bdb);; Spider (*.spi);; Imagic (*.img);; PNG (*.png);; Gatan (*.dm2 *.dm3 *.dm4);; FEI (*.ser);; EM (*.em);; ICOS (*.icos);; Amira (*.am);; DF3 (*.d3);; FITS (*.fts);; LST (*.lst);; LSTFAST (*.lsx *.lst);; OMAP (*.omap);; PGM (*.pgm);; PIF (*.pif);; SAL (*.img );; SITUS (*.situs);; V4L (*.v4l);; VTK (*.vtk);; XPLOR (*.xplor);; All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
				try: # Check if the path is bdb
					file_path = translate_to_bdb_path(file_path) # Convert the standard path to bdb key if possible.
				except ValueError:  # If the path is not bdb, we will receive this exception
					pass # This is not bdb path. Then, use standard path
		elif file_format == "mrc2d_mic_one":
			name = QFileDialog.getOpenFileName(self, "Select MRC micrograph file", SXLookFeelConst.file_dialog_dir, "MRC micrograph files (*.mrc);; MRCS (*.mrcs)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "mic_one":
			# NOTE: Toshio Moriya 2018/01/25
			# Currently, the distinction between MRC and MRCS is not always used, and 
			# MRC can be also micrograph stack dependes on external programs (e.g. unblur, summovie)...
			# 
			# Read not supported: ;; JPEG (*.jpg *.jpeg)
			# Only stack: ;; MRCS (*.mrcs)
			name = QFileDialog.getOpenFileName(self, "Select any micrograph file", SXLookFeelConst.file_dialog_dir, "Typical micrograph files (*.mrc *.mrcs *.tif *.tiff *.hdf *.bdb *.spi *.img);; MRC (*.mrc);; MRCS (*.mrcs);; TIFF (*.tif *.tiff);; HDF (*.hdf);; BDB (*.bdb);; Spider (*.spi);; Imagic (*.img);; PNG (*.png);; Gatan (*.dm2 *.dm3 *.dm4);; FEI (*.ser);; EM (*.em);; ICOS (*.icos);; Amira (*.am);; DF3 (*.d3);; FITS (*.fts);; LST (*.lst);; LSTFAST (*.lsx *.lst);; OMAP (*.omap);; PGM (*.pgm);; PIF (*.pif);; SAL (*.img );; SITUS (*.situs);; V4L (*.v4l);; VTK (*.vtk);; XPLOR (*.xplor);; All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
				try: # Check if the path is bdb
					file_path = translate_to_bdb_path(file_path) # Convert the standard path to bdb key if possible.
				except ValueError:  # If the path is not bdb, we will receive this exception
					pass # This is not bdb path. Then, use standard path
		elif file_format == "mrc2d_mic_one_list":
			# NOTE: Toshio Moriya 2018/01/25
			# Currently, the distinction between MRC and MRCS is not always used, and 
			# MRC can be also micrograph stack dependes on external programs (e.g. unblur, summovie)...
			# 
			name = QFileDialog.getOpenFileNames(self, "Select MRC micrograph files", SXLookFeelConst.file_dialog_dir, "MRC files (*.mrc);; MRCS (*.mrcs)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path_list = name[0]
			else:
				file_path_list = name
			# Use relative path.
			for a_file_path in file_path_list:
				file_path += str(SXLookFeelConst.format_path(a_file_path)) + " "
		elif file_format == "mic_one_list":
			# NOTE: Toshio Moriya 2018/01/25
			# Currently, the distinction between MRC and MRCS is not always used, and 
			# MRC can be also micrograph stack dependes on external programs (e.g. unblur, summovie)...
			# 
			# Read not supported: ;; JPEG (*.jpg *.jpeg)
			# Only stack: ;; MRCS (*.mrcs)
			name = QFileDialog.getOpenFileNames(self, "Select any micrograph files", SXLookFeelConst.file_dialog_dir, "Typical micrograph files (*.mrc *.mrcs *.tif *.tiff *.hdf *.bdb *.spi *.img);; MRC (*.mrc);; MRCS (*.mrcs);; TIFF (*.tif *.tiff);; HDF (*.hdf);; BDB (*.bdb);; Spider (*.spi);; Imagic (*.img);; PNG (*.png);; Gatan (*.dm2 *.dm3 *.dm4);; FEI (*.ser);; EM (*.em);; ICOS (*.icos);; Amira (*.am);; DF3 (*.d3);; FITS (*.fts);; LST (*.lst);; LSTFAST (*.lsx *.lst);; OMAP (*.omap);; PGM (*.pgm);; PIF (*.pif);; SAL (*.img );; SITUS (*.situs);; V4L (*.v4l);; VTK (*.vtk);; XPLOR (*.xplor);; All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path_list = name[0]
			else:
				file_path_list = name
			# Use relative path.
			for a_file_path in file_path_list:
				a_file_path = str(SXLookFeelConst.format_path(a_file_path))
				try: # Check if the path is bdb
					a_file_path = translate_to_bdb_path(a_file_path) # Convert the standard path to bdb key if possible.
				except ValueError:  # If the path is not bdb, we will receive this exception
					pass # This is not bdb path. Then, use standard path
				file_path += a_file_path + " "
		elif file_format == "mrc2d_mic_stack":
			name = QFileDialog.getOpenFileName(self, "Select MRC movie file", SXLookFeelConst.file_dialog_dir, "MRC movie files (*.mrcs);; MRC (*.mrc)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "mic_stack":
			# NOTE: Toshio Moriya 2018/01/25
			# Currently, the distinction between MRC and MRCS is not always used, and 
			# MRC can be also micrograph stack dependes on external programs (e.g. unblur, summovie)...
			# 
			# Read not supported: ;; JPEG (*.jpg *.jpeg)
			# 2D image stack not supported: ;; Gatan (*.dm2 *.dm3);; EM (*.em);; ICOS (*.icos);; Amira (*.am);; DF3 (*.d3);; FITS (*.fts);; OMAP (*.omap);; PGM (*.pgm);; PNG (*.png);; SAL (*.hdr *.img);; SITUS (*.situs);; TIFF (*.tif *.tiff);; V4L (*.v4l);; VTK (*.vtk);; XPLOR (*.xplor)
			# Maybe only single 2D image: ;; MRC (*.mrc)
			name = QFileDialog.getOpenFileName(self, "Select any movie file", SXLookFeelConst.file_dialog_dir, "Typical movie files (*.mrcs *.mrc *.bdb *.hdf *.spi *.img );; MRCS (*.mrcs);; MRC (*.mrc);; BDB (*.bdb);; HDF (*.hdf);; Spider (*.spi);; Imagic (*.img *hed);; Gatan (*.dm4);; FEI (*.ser);; LST (*.lst);; LSTFAST (*.lsx *.lst);; PIF (*.pif);; All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
				try: # Check if the path is bdb
					file_path = translate_to_bdb_path(file_path) # Convert the standard path to bdb key if possible.
				except ValueError:  # If the path is not bdb, we will receive this exception
					pass # This is not bdb path. Then, use standard path
		elif file_format == "hdf2d_one":
			name = QFileDialog.getOpenFileName(self, "Select HDF image file", SXLookFeelConst.file_dialog_dir, "HDF image files (*.hdf)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "data2d_one":
			# Read not supported: ;; JPEG (*.jpg *.jpeg)
			# Maybe only 2D image stack: ;; MRCS (*.mrcs)
			name = QFileDialog.getOpenFileName(self, "Select any image file", SXLookFeelConst.file_dialog_dir, "Typical image files (*.hdf *.bdb *.mrc *.spi *.img *.tif *.tiff *.png);; HDF (*.hdf);; BDB (*.bdb);; MRC (*.mrc);; Spider (*.spi);; Imagic (*.img);; TIFF (*.tif *.tiff);; PNG (*.png);; Gatan (*.dm2 *.dm3 *.dm4);; FEI (*.ser);; EM (*.em);; ICOS (*.icos);; Amira (*.am);; DF3 (*.d3);; FITS (*.fts);; LST (*.lst);; LSTFAST (*.lsx *.lst);; OMAP (*.omap);; PGM (*.pgm);; PIF (*.pif);; SAL (*.img );; SITUS (*.situs);; V4L (*.v4l);; VTK (*.vtk);; XPLOR (*.xplor);; All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
				try: # Check if the path is bdb
					file_path = translate_to_bdb_path(file_path) # Convert the standard path to bdb key if possible.
				except ValueError:  # If the path is not bdb, we will receive this exception
					pass # This is not bdb path. Then, use standard path
		elif file_format == "bdb2d_stack":
			name = QFileDialog.getOpenFileName(self, "Select BDB image stack file", SXLookFeelConst.file_dialog_dir, "BDB files (*.bdb)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
				file_path = translate_to_bdb_path(file_path)
		elif file_format == "data2d_stack":
			# NOTE: Toshio Moriya 2018/01/25
			# Currently, this case is not used. Instead, using bdb2d_stack
			# 
			# Read not supported: ;; JPEG (*.jpg *.jpeg)
			# 2D image stack not supported: ;; Gatan (*.dm2 *.dm3);; EM (*.em);; ICOS (*.icos);; Amira (*.am);; DF3 (*.d3);; FITS (*.fts);; OMAP (*.omap);; PGM (*.pgm);; PNG (*.png);; SAL (*.hdr *.img);; SITUS (*.situs);; TIFF (*.tif *.tiff);; V4L (*.v4l);; VTK (*.vtk);; XPLOR (*.xplor)
			# Maybe only single 2D image: ;; MRC (*.mrc)
			name = QFileDialog.getOpenFileName(self, "Select any image stack file", SXLookFeelConst.file_dialog_dir, "Typical image stack files (*.bdb *.hdf *.mrcs *.spi *.img );; BDB (*.bdb);; HDF (*.hdf);; MRCS (*.mrcs);; Spider (*.spi);; Imagic (*.img *hed);; Gatan (*.dm4);; FEI (*.ser);; LST (*.lst);; LSTFAST (*.lsx *.lst);; PIF (*.pif);; All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
				try: # Check if the path is bdb
					file_path = translate_to_bdb_path(file_path) # Convert the standard path to bdb key if possible.
				except ValueError:  # If the path is not bdb, we will receive this exception
					pass # This is not bdb path. Then, use standard path
		elif file_format == "hdf3d_one":
			name = QFileDialog.getOpenFileName(self, "Select HDF volume file", SXLookFeelConst.file_dialog_dir, "HDF volume files (*.hdf)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "data3d_one":
			# Read not supported: ;; JPEG (*.jpg *.jpeg)
			# 3D volume not supported: ;; Gatan (*.dm2 *.dm3);; FEI (*.ser);; SAL (*.hdr *.img);; PGM (*.pgm);; PNG (*.png);; TIFF (*.tif *.tiff);; V4L (*.v4l)
			# Maybe only 3D volume stack: ;; MRCS (*.mrcs)
			name = QFileDialog.getOpenFileName(self, "Select any volume file", SXLookFeelConst.file_dialog_dir, "Typical volume files (*.hdf *.bdb *.mrc *.spi *.img);; HDF (*.hdf);; BDB (*.bdb);; MRC (*.mrc);; Spider (*.spi);; Imagic (*.img);; Gatan (*.dm4);; EM (*.em);; ICOS (*.icos);; Amira (*.am);; DF3 (*.d3);; FITS (*.fts);; LST (*.lst);; LSTFAST (*.lsx *.lst);; OMAP (*.omap);; PIF (*.pif);; SITUS (*.situs);; VTK (*.vtk);; XPLOR (*.xplor);; All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
				try: # Check if the path is bdb
					file_path = translate_to_bdb_path(file_path) # Convert the standard path to bdb key if possible.
				except ValueError:  # If the path is not bdb, we will receive this exception
					pass # This is not bdb path. Then, use standard path
		elif file_format == "hdf3d_stack":
			# NOTE: Toshio Moriya 2018/01/25
			# Currently, this case is not used.
			# 
			name = QFileDialog.getOpenFileName(self, "Select HDF volume stack file", SXLookFeelConst.file_dialog_dir, "HDF volume stack files (*.hdf)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "data3d_stack":
			# NOTE: Toshio Moriya 2018/01/25
			# Currently, this case is not used.
			# 
			# Read not supported: ;; JPEG (*.jpg *.jpeg)
			# 3D volume stack not supported: ;; Gatan (*.dm2)
			# Maybe 3D volume stack not supported: ;; Gatan (*.dm3 *.dm4);; FEI (*.ser);; EM (*.em);; ICOS (*.icos);; Spider (*.spi);; Amira (*.am);; DF3 (*.d3);; FITS (*.fts);; LST (*.lst);; LSTFAST (*.lsx *.lst);; OMAP (*.omap);; PGM (*.pgm);; PIF (*.pif);; PNG (*.png);; SAL (*.hdr *.img);; SITUS (*.situs);; TIFF (*.tif *.tiff);; V4L (*.v4l);; VTK (*.vtk);; XPLOR (*.xplor)
			# Maybe only sigle 3D volume: ;; MRC (*.mrc)
			name = QFileDialog.getOpenFileName(self, "Select any volume stack file", SXLookFeelConst.file_dialog_dir, "Typical volume stack files (*.hdf *.bdb *.mrcs *.img);; HDF (*.hdf);; BDB (*.bdb);; MRCS (*.mrcs);; Imagic (*.img *hed);; All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
				try: # Check if the path is bdb
					file_path = translate_to_bdb_path(file_path) # Convert the standard path to bdb key if possible.
				except ValueError:  # If the path is not bdb, we will receive this exception
					pass # This is not bdb path. Then, use standard path
		elif file_format == "select_mic_both":
			name = QFileDialog.getOpenFileName(self, "Select micrograph/movie selection file", SXLookFeelConst.file_dialog_dir, "Micrograph/Movie selection files (*.txt);; All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "select_mic_one":
			name = QFileDialog.getOpenFileName(self, "Select micrograph selection file", SXLookFeelConst.file_dialog_dir, "Micrograph selection files (*.txt);; All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "select_mic_stack":
			name = QFileDialog.getOpenFileName(self, "Select micrograph movie selection file", SXLookFeelConst.file_dialog_dir, "Micrograph movie selection files (*.txt);; All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "select_data2d_stack":
			name = QFileDialog.getOpenFileName(self, "Select image selection file", SXLookFeelConst.file_dialog_dir, "Image selection files (*.txt);; All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "select_drift_params":
			name = QFileDialog.getOpenFileName(self, "Select drift shift params selection file", SXLookFeelConst.file_dialog_dir, "Drift shift params selection files (*.txt);; All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "params_any_txt":
			name = QFileDialog.getOpenFileName(self, "Select parameters file", SXLookFeelConst.file_dialog_dir, "Parameters text files (*.txt);; All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "params_proj_txt":
			name = QFileDialog.getOpenFileName(self, "Select projection parameters file", SXLookFeelConst.file_dialog_dir, "Projection parameters files (*.txt);; All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "params_coords_box":
			name = QFileDialog.getOpenFileName(self, "Select EMAN BOX coordinates file", SXLookFeelConst.file_dialog_dir, "EMAN BOX coordinates files (*.box)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "params_coords_any":
			name = QFileDialog.getOpenFileName(self, "Select any coordinates file", SXLookFeelConst.file_dialog_dir, "Typical coordinates files (*.box *.json *.dat *.txt);; EMAN BOX (*.box);; EMAN2 JSON (*.json);; SPIDER DAT (*.dat);; SPHIRE TXT (*.txt);; All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "params_cter_txt":
			name = QFileDialog.getOpenFileName(self, "Select CTER partres parameters file", SXLookFeelConst.file_dialog_dir, "CTER partres parameters files (*.txt)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "params_rebox_rbx":
			name = QFileDialog.getOpenFileName(self, "Select SPHIRE rebox file", SXLookFeelConst.file_dialog_dir, "SPHIRE rebox files (*.rbx)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "params_drift_txt":
			name = QFileDialog.getOpenFileName(self, "Select drift shift parameters file", SXLookFeelConst.file_dialog_dir, "Drift shift parameters files (*.txt);; All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "rot_matrix":
			name = QFileDialog.getOpenFileName(self, "Select rotational matrix file", SXLookFeelConst.file_dialog_dir, "Rotational matrix files (*.txt);; All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "params_relion_star":
			name = QFileDialog.getOpenFileName(self, "Select RELION STAR file", SXLookFeelConst.file_dialog_dir, "RELION STAR files (*.star);; All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "params_any_json":
			name = QFileDialog.getOpenFileName(self, "Select JSON file", SXLookFeelConst.file_dialog_dir, "JSON files (*.json);; All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "params_any_h5":
			name = QFileDialog.getOpenFileName(self, "Select h5 file", SXLookFeelConst.file_dialog_dir, "h5 files (*.h5);; All files (*)", options = QFileDialog.DontUseNativeDialog)
			# Use relative path.
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "spectrum1d":
			name = QFileDialog.getOpenFileName(self, "Select 1D power spectrum file", SXLookFeelConst.file_dialog_dir, "1D power spectrum files (*.txt);; All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "mtf":
			name = QFileDialog.getOpenFileName(self, "Select MTF data file", SXLookFeelConst.file_dialog_dir, "MTF data files (*.txt);; All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "pdb":
			name = QFileDialog.getOpenFileName(self, "Select PDB data file", SXLookFeelConst.file_dialog_dir, "PDB data files (*.pdb *.pdb*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "exe":
			name = QFileDialog.getOpenFileName(self, "Select executable file", SXLookFeelConst.file_dialog_dir, "All files (*);; Python (*.py) ;; Executable files (*.exe)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)
		elif file_format == "py":
			name = QFileDialog.getOpenFileName(self, "Select Python script file", SXLookFeelConst.file_dialog_dir, "PY files (*.py)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use full path
###		elif file_format == "bdb":
###			file_path = QFileDialog.getOpenFileName(self, "Select BDB image file", SXLookFeelConst.file_dialog_dir, "BDB files (*.bdb)", options = QFileDialog.DontUseNativeDialog)
###			# Use relative path.
###			if file_path:
###				file_path = SXLookFeelConst.format_path(file_path)
###				file_path = translate_to_bdb_path(file_path)
###		elif file_format == "any_file_list" or file_format == "any_image_list":
###			file_path_list = QFileDialog.getOpenFileNames(self, "Select files", SXLookFeelConst.file_dialog_dir, "All files (*)", options = QFileDialog.DontUseNativeDialog)
###			# Use relative path.
###			for a_file_path in file_path_list:
###				file_path += SXLookFeelConst.format_path(str(a_file_path)) + " "
		elif file_format == 'submission_template':
			if 'SPHIRE_SUBMISSION_SCRIPT_TEMPLATE_FOLDER' in os.environ:
				template_folder = os.environ['SPHIRE_SUBMISSION_SCRIPT_TEMPLATE_FOLDER']
			else:
				template_folder = SXLookFeelConst.file_dialog_dir
			name = QFileDialog.getOpenFileName(self, "Select any file", template_folder, "All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
		else:
			if file_format:
				name = QFileDialog.getOpenFileName(self, "Select %s file" % (file_format.upper()), SXLookFeelConst.file_dialog_dir, "%s files (*.%s)"  % (file_format.upper(), file_format), options = QFileDialog.DontUseNativeDialog)
			else:
				name = QFileDialog.getOpenFileName(self, "Select any file", SXLookFeelConst.file_dialog_dir, "All files (*)", options = QFileDialog.DontUseNativeDialog)
			if isinstance(name, tuple):
				file_path = str(name[0])
			else:
				file_path = str(name)
			# Use relative path.
			if file_path:
				file_path = SXLookFeelConst.format_path(file_path)

		if file_path != "":
			target_widget.setText(file_path)

	def select_dir(self, target_widget):
		dir_path = QFileDialog.getExistingDirectory(self, "Select directory", SXLookFeelConst.file_dialog_dir, options = QFileDialog.ShowDirsOnly | QFileDialog.DontResolveSymlinks | QFileDialog.DontUseNativeDialog)
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

class SXLineEdit(QtGui.QLineEdit):
	my_changed_state = QtCore.pyqtSignal(str, str)

	def __init__(self, dependency, parent=None):
		super(SXLineEdit, self).__init__(parent)
		self.dependency = dependency
		self.textChanged.connect(lambda x: self.my_changed_state.emit(self.dependency, x))


class SXCheckBox(QtGui.QCheckBox):
	my_changed_state = QtCore.pyqtSignal(str, int)

	def __init__(self, dependency, parent=None):
		super(SXCheckBox, self).__init__(parent)
		self.dependency = dependency
		self.stateChanged.connect(lambda x: self.my_changed_state.emit(self.dependency, x))


# ========================================================================================
class SXCmdTab(QWidget):
	def __init__(self, name, helical, parent=None):
		super(SXCmdTab, self).__init__(parent)

		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.name = name
		self.sxcmdwidget = parent

		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# local constants
		required_cmd_token_restore_tooltip = "Please enter the value manually"
		locked_cmd_token_restore_tooltip = "This value is locked"
		const_cmd_token_restore_tooltip = "Retrieve the registed constant value for this parameter"
		default_cmd_token_restore_tooltip = "Retrieve this default value"

		# Setting for layout
		grid_row_origin = 0; grid_col_origin = 0
		title_row_span = 1; title_col_span = 2
		short_info_row_span = 1; short_info_col_span = 5
		func_btn_row_span = 1; func_btn_col_span = 2
		token_label_row_span = 1; token_label_col_span = 4
		token_widget_row_span = 1; token_widget_col_span = 1.5
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
		title_layout.setContentsMargins(SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin)
		title_layout.setSpacing(SXLookFeelConst.grid_spacing)
#		title_layout.setColumnMinimumWidth(grid_col_origin + token_label_col_span, token_widget_min_width)
#		title_layout.setColumnMinimumWidth(grid_col_origin + token_label_col_span + token_widget_col_span, token_widget_min_width)
#		title_layout.setColumnMinimumWidth(grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_min_width)
#		title_layout.setColumnMinimumWidth(grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_min_width)
		grid_layout = QGridLayout()
		grid_layout.setContentsMargins(SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin)
		grid_layout.setSpacing(SXLookFeelConst.grid_spacing)
		grid_layout.setColumnMinimumWidth(grid_col_origin + token_label_col_span, token_widget_min_width)
		grid_layout.setColumnMinimumWidth(grid_col_origin + token_label_col_span + token_widget_col_span, token_widget_min_width)
		grid_layout.setColumnMinimumWidth(grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_min_width)
		grid_layout.setColumnMinimumWidth(grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_min_width)
		submit_layout = QGridLayout()
		submit_layout.setContentsMargins(SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin)
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
		cmd_token_apix = None
		cmd_token_abs_freq_list = []
		for cmd_token in self.sxcmdwidget.sxcmd.token_list:
			# Keep some command tokens for further setting after the loop
			if cmd_token.type == "apix":
				cmd_token_apix = cmd_token
			elif cmd_token.type == "abs_freq":
				cmd_token_abs_freq_list.append(cmd_token)
			# else: # Do nothing
			curr_group = cmd_token.group
			if cmd_token.filament_tab and helical:
				curr_group = cmd_token.filament_tab

			if curr_group == tab_group:
				# First, handle very special cases
				if cmd_token.type == "user_func":
					n_widgets = 2 # function type has two line edit boxes
					cmd_token_restore_widget = [None] * n_widgets
					cmd_token_widget = [None] * n_widgets
					cmd_token_subwidget_left = None
					cmd_token_subwidget_right = None
					cmd_token_calculator_dialog = None

					# Define custom style for restore widgets
					custom_style = "QPushButton {color:gray; }"

					# Create widgets for user function name
					widget_index = 0
					temp_label = QLabel(cmd_token.label[widget_index])
					temp_label.setMinimumWidth(token_label_min_width)
					grid_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)

					assert(cmd_token.is_required == False)
					cmd_token_restore_widget[widget_index] = QPushButton("%s" % cmd_token.restore[widget_index][helical][0])
					cmd_token_restore_widget[widget_index].setStyleSheet(custom_style)
					cmd_token_restore_widget[widget_index].setToolTip('<FONT>'+default_cmd_token_restore_tooltip+'</FONT>')
					grid_layout.addWidget(cmd_token_restore_widget[widget_index], grid_row, grid_col_origin + token_label_col_span, token_widget_row_span, token_widget_col_span)

					# cmd_token_widget[widget_index] = SXLineEdit(self)
					cmd_token_widget[widget_index] = SXLineEdit(cmd_token.key_base)
					cmd_token_widget[widget_index].setText(cmd_token.restore[widget_index][helical][0])
					cmd_token_widget[widget_index].setToolTip('<FONT>'+cmd_token.help[widget_index]+'</FONT>')
					grid_layout.addWidget(cmd_token_widget[widget_index], grid_row, grid_col_origin + token_label_col_span + token_widget_col_span, token_widget_row_span, token_widget_col_span)

					cmd_token_restore_widget[widget_index].clicked.connect(partial(self.handle_restore_widget_event, cmd_token, widget_index))

					grid_row +=  1

					# Create widgets for external file path containing above user function
					widget_index = 1
					temp_label = QLabel(cmd_token.label[widget_index])
					grid_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)

					assert(cmd_token.is_required == False)
					cmd_token_restore_widget[widget_index] = QPushButton("%s" % cmd_token.restore[widget_index][helical][0])
					cmd_token_restore_widget[widget_index].setStyleSheet(custom_style)
					cmd_token_restore_widget[widget_index].setToolTip('<FONT>'+default_cmd_token_restore_tooltip+'</FONT>')
					grid_layout.addWidget(cmd_token_restore_widget[widget_index], grid_row, grid_col_origin + token_label_col_span, token_widget_row_span, token_widget_col_span)

					cmd_token_widget[widget_index] = SXLineEdit(cmd_token.key_base)
					cmd_token_widget[widget_index].setText(cmd_token.restore[widget_index][helical][0]) # Because default user functions is internal
					cmd_token_widget[widget_index].setToolTip('<FONT>'+cmd_token.help[widget_index]+'</FONT>')
					grid_layout.addWidget(cmd_token_widget[widget_index], grid_row, grid_col_origin + token_label_col_span + token_widget_col_span, token_widget_row_span, token_widget_col_span)

					cmd_token_restore_widget[widget_index].clicked.connect(partial(self.handle_restore_widget_event, cmd_token, widget_index))

					file_format = "py"
					temp_btn = QPushButton("Select Python script")
					temp_btn.setToolTip('<FONT>'+"Display open file dialog to select .%s python script file</FONT>" % file_format)
					grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
					temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget[widget_index], file_format))

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
					cmd_token_restore_widget = None
					cmd_token_restore_tooltip = default_cmd_token_restore_tooltip
					cmd_token_widget = None
					cmd_token_subwidget_left = None
					cmd_token_subwidget_right = None
					cmd_token_calculator_dialog = None

					restores = cmd_token.restore[helical]
					if cmd_token.type in ("bool", "bool_ignore"):
						if len(restores) == 1:
							btn_name = "NO"
							is_btn_enable = True
							custom_style = "QPushButton {color:gray; }"
							if restores[0]:
								btn_name = "YES"
							if cmd_token.type in list(parent.sxconst_set.dict.keys()):
								custom_style = "QPushButton {color:green; }"
								cmd_token_restore_tooltip = const_cmd_token_restore_tooltip
							elif cmd_token.is_required:
								if cmd_token.is_locked:
									btn_name = "locked"
									custom_style = "QPushButton {color:blue; }"
									is_btn_enable = False
									cmd_token_restore_tooltip = locked_cmd_token_restore_tooltip
								else:
									btn_name = "required"
									custom_style = "QPushButton {color:red; }"
									is_btn_enable = False
									cmd_token_restore_tooltip = required_cmd_token_restore_tooltip
								
							cmd_token_restore_widget = QPushButton("%s" % btn_name)
							cmd_token_restore_widget.setStyleSheet(custom_style)
							cmd_token_restore_widget.setEnabled(is_btn_enable)
							grid_layout.addWidget(
								cmd_token_restore_widget,
								grid_row,
								grid_col_origin + token_label_col_span,
								token_widget_row_span,
								token_widget_col_span
								)
							cmd_token_restore_widget.clicked.connect(partial(self.handle_restore_widget_event, cmd_token))
						else:
							custom_style = "QComboBox {color:gray; }"
							is_btn_enable = True
							if cmd_token.type in list(parent.sxconst_set.dict.keys()):
								custom_style = "QComboBox {color:green; }"
								cmd_token_restore_tooltip = const_cmd_token_restore_tooltip
							elif cmd_token.is_required:
								if cmd_token.is_locked:
									btn_name = "locked"
									custom_style = "QComboBox {color:blue; }"
									is_btn_enable = False
									cmd_token_restore_tooltip = locked_cmd_token_restore_tooltip
								else:
									btn_name = "required"
									custom_style = "QComboBox {color:red; }"
									is_btn_enable = False
									cmd_token_restore_tooltip = required_cmd_token_restore_tooltip
							cmd_token_restore_widget = QComboBox()
							cmd_token_restore_widget.setEditable(True)
							cmd_token_restore_widget.lineEdit().setReadOnly(True)
							cmd_token_restore_widget.lineEdit().setAlignment(QtCore.Qt.AlignCenter)

							cmd_token_restore_widget.addItems(restores)
							cmd_token_restore_widget.setStyleSheet(custom_style)
							cmd_token_restore_widget.setEnabled(is_btn_enable)
							grid_layout.addWidget(
								cmd_token_restore_widget,
								grid_row,
								grid_col_origin + token_label_col_span,
								token_widget_row_span,
								token_widget_col_span
								)
							cmd_token_restore_widget.activated.connect(partial(self.handle_restore_widget_event, cmd_token))

						# construct new widget(s) for this command token
						cmd_token_widget = SXCheckBox(cmd_token.key_base)
						if restores[0] == True:
							cmd_token_widget.setCheckState(Qt.Checked)
						else:
							cmd_token_widget.setCheckState(Qt.Unchecked)
						cmd_token_widget.setEnabled(not cmd_token.is_locked)
						grid_layout.addWidget(cmd_token_widget, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span, token_widget_row_span, token_widget_col_span)


					else:

						if len(restores) == 1:
							btn_name = "%s" % restores[0]
							custom_style = "QPushButton {color:gray; }"
							is_btn_enable = True
							if cmd_token.type in list(parent.sxconst_set.dict.keys()):
								custom_style = "QPushButton {color:green; }"
								cmd_token_restore_tooltip = const_cmd_token_restore_tooltip
							elif cmd_token.is_required:
								if cmd_token.is_locked:
									btn_name = "locked"
									custom_style = "QPushButton {color:blue; }"
									is_btn_enable = False
									cmd_token_restore_tooltip = locked_cmd_token_restore_tooltip
								else:
									btn_name = "required"
									custom_style = "QPushButton {color:red; }"
									is_btn_enable = False
									cmd_token_restore_tooltip = required_cmd_token_restore_tooltip
							cmd_token_restore_widget = QPushButton("%s" % btn_name)
							cmd_token_restore_widget.setStyleSheet(custom_style)
							cmd_token_restore_widget.setEnabled(is_btn_enable)
							grid_layout.addWidget(
								cmd_token_restore_widget,
								grid_row,
								grid_col_origin + token_label_col_span,
								token_widget_row_span,
								token_widget_col_span
								)
							cmd_token_restore_widget.clicked.connect(partial(self.handle_restore_widget_event, cmd_token))
						else:
							custom_style = "QComboBox {color:gray; }"
							is_btn_enable = True
							if cmd_token.type in list(parent.sxconst_set.dict.keys()):
								custom_style = "QComboBox {color:green; }"
								cmd_token_restore_tooltip = const_cmd_token_restore_tooltip
							elif cmd_token.is_required:
								if cmd_token.is_locked:
									btn_name = "locked"
									custom_style = "QComboBox {color:blue; }"
									is_btn_enable = False
									cmd_token_restore_tooltip = locked_cmd_token_restore_tooltip
								else:
									btn_name = "required"
									custom_style = "QComboBox {color:red; }"
									is_btn_enable = False
									cmd_token_restore_tooltip = required_cmd_token_restore_tooltip
							cmd_token_restore_widget = QComboBox()
							cmd_token_restore_widget.setEditable(True)
							cmd_token_restore_widget.lineEdit().setReadOnly(True)
							cmd_token_restore_widget.lineEdit().setAlignment(QtCore.Qt.AlignCenter)

							cmd_token_restore_widget.addItems(restores)
							cmd_token_restore_widget.setStyleSheet(custom_style)
							cmd_token_restore_widget.setEnabled(is_btn_enable)
							grid_layout.addWidget(
								cmd_token_restore_widget,
								grid_row,
								grid_col_origin + token_label_col_span,
								token_widget_row_span,
								token_widget_col_span
								)
							cmd_token_restore_widget.activated.connect(partial(self.handle_restore_widget_event, cmd_token))

						cmd_token_widget = SXLineEdit(cmd_token.key_base)
						cmd_token_widget.setText(str(restores[0]))
						cmd_token_widget.setEnabled(not cmd_token.is_locked)
						grid_layout.addWidget(cmd_token_widget, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span, token_widget_row_span, token_widget_col_span)


						if cmd_token.type == "displayable_list":
							file_format = cmd_token.type
							temp_btn = QPushButton("Select any displayables")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select displayable data files of any supported formats</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "INVISIBLE"
							temp_btn = QPushButton("%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "data2d3d_both":
							file_format = cmd_token.type
							temp_btn = QPushButton("Select any image/volume")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select an image/volume file of any supported format</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "INVISIBLE"
							temp_btn = QPushButton("%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "mic_both":
							file_format = "mrc2d_mic_both"
							temp_btn = QPushButton("Select MRC mic/movie")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a MRC format micrograph or movie file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = cmd_token.type
							temp_btn = QPushButton("Select any mic/movie")
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a micrograph or movie file of any supported format</FONT>")
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
						elif cmd_token.type == "mic_one":
							file_format = "mrc2d_mic_one"
							temp_btn = QPushButton("Select MRC micrograph")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a MRC format micrograph file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = cmd_token.type
							temp_btn = QPushButton("Select any micrograph")
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a micrograph file of any supported format</FONT>")
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
						elif cmd_token.type == "mic_one_list":
							file_format = "mrc2d_mic_one_list"
							temp_btn = QPushButton("Select MRC micrographs")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select MRC format micrograph files</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = cmd_token.type
							temp_btn = QPushButton("Select any micrographs")
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select micrograph files of any supported formats</FONT>")
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
						elif cmd_token.type == "mic_stack":
							file_format = "mrc2d_mic_stack"
							temp_btn = QPushButton("Select MRC movie")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a MRC format movie file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = cmd_token.type
							temp_btn = QPushButton("Select any movie")
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a movie file of any supported format</FONT>")
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
						elif cmd_token.type == "data2d_one":
							file_format = "hdf2d_one"
							temp_btn = QPushButton("Select HDF image")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a HDF format image file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = cmd_token.type
							temp_btn = QPushButton("Select any image")
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a image file of any supported format</FONT>")
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
						elif cmd_token.type == "bdb2d_stack" or cmd_token.type == "output_bdb2d_stack":
							file_format = "bdb2d_stack"
							temp_btn = QPushButton("Select BDB image stack")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a BDB format image stack file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "INVISIBLE"
							temp_btn = QPushButton("%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "data2d_stack":
							file_format = "bdb2d_stack"
							temp_btn = QPushButton("Select BDB image stack")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a BDB format image stack file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = cmd_token.type
							temp_btn = QPushButton("Select any image stack")
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a image stack file of any supported format</FONT>")
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
						elif cmd_token.type == "data3d_one":
							file_format = "hdf3d_one"
							temp_btn = QPushButton("Select HDF volume")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a HDF format volume file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = cmd_token.type
							temp_btn = QPushButton("Select any volume")
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a volume file of any supported format</FONT>")
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
						elif cmd_token.type == "data3d_stack":
							file_format = "hdf3d_stack"
							temp_btn = QPushButton("Select HDF volume stack")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a HDF format volume stack file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = cmd_token.type
							temp_btn = QPushButton("Select any volume stack")
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a volume stack file of any supported format</FONT>")
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
						elif cmd_token.type == "select_mic_both":
							file_format = cmd_token.type
							temp_btn = QPushButton("Select mic/movie list")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a micrograph/movie selection file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "INVISIBLE"
							temp_btn = QPushButton("%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "select_mic_one":
							file_format = cmd_token.type
							temp_btn = QPushButton("Select micrograph list")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a micrograph selection file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "INVISIBLE"
							temp_btn = QPushButton("%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "select_mic_one_ext":
							file_format = "select_mic_one"
							temp_btn = QPushButton("Select micrograph list")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a micrograph selection file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "mic_one"
							temp_btn = QPushButton("Select any micrograph")
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a micrograph file of any supported format</FONT>")
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
						elif cmd_token.type == "select_mic_stack":
							file_format = cmd_token.type
							temp_btn = QPushButton("Select movie list")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a micrograph movie selection file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "INVISIBLE"
							temp_btn = QPushButton("%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "select_data2d_stack":
							file_format = cmd_token.type
							temp_btn = QPushButton("Select image list")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a image selection file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "INVISIBLE"
							temp_btn = QPushButton("%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "select_drift_params":
							file_format = cmd_token.type
							temp_btn = QPushButton("Select drift params list")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a drift shift parameters selection file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "INVISIBLE"
							temp_btn = QPushButton("%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "params_any_txt":
							file_format = cmd_token.type
							temp_btn = QPushButton("Select parameters text")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a parameters text file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "INVISIBLE"
							temp_btn = QPushButton("%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "params_proj_txt":
							file_format = cmd_token.type
							temp_btn = QPushButton("Select projection params")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a projection parameters file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "INVISIBLE"
							temp_btn = QPushButton("%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "params_coords_any":
							file_format = "params_coords_box"
							temp_btn = QPushButton("Select BOX coordinates")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a EMAN BOX coordinates file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = cmd_token.type
							temp_btn = QPushButton("Select any coordinates")
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a coordinates file of any supported format</FONT>")
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
						elif cmd_token.type == "params_cter_txt":
							file_format = cmd_token.type
							temp_btn = QPushButton("Select CTER partres")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a CTER partres parameters file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "INVISIBLE"
							temp_btn = QPushButton("%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "params_rebox_rbx":
							file_format = cmd_token.type
							temp_btn = QPushButton("Select SPHIRE rebox")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a SPHIRE rebox file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "INVISIBLE"
							temp_btn = QPushButton("%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "params_drift_txt":
							file_format = cmd_token.type
							temp_btn = QPushButton("Select drift shift params")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a drift shift parameters file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "INVISIBLE"
							temp_btn = QPushButton("%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "rot_matrix":
							file_format = cmd_token.type
							temp_btn = QPushButton("Select matrix file")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a rotational matrix file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "INVISIBLE"
							temp_btn = QPushButton("%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "params_relion_star":
							file_format = cmd_token.type
							temp_btn = QPushButton("Select RELION STAR file")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a RELION STAR file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "INVISIBLE"
							temp_btn = QPushButton("%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "params_any_json":
							file_format = cmd_token.type
							temp_btn = QPushButton("Select JSON file")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select a JSON file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "INVISIBLE"
							temp_btn = QPushButton("%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "params_any_h5":
							file_format = cmd_token.type
							temp_btn = QPushButton("Select h5 file")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dailog to select a h5 file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "INVISIBLE"
							temp_btn = QPushButton("%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "spectrum1d":
							file_format = cmd_token.type
							temp_btn = QPushButton("Select power spectrum")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a 1D power spectrum file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "INVISIBLE"
							temp_btn = QPushButton("%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "mtf":
							file_format = cmd_token.type
							temp_btn = QPushButton("Select MTF data")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a MTF data file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "INVISIBLE"
							temp_btn = QPushButton("%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "pdb":
							file_format = cmd_token.type
							temp_btn = QPushButton("Select PDB data")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a PDB data file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "INVISIBLE"
							temp_btn = QPushButton("%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "exe":
							file_format = cmd_token.type
							temp_btn = QPushButton("Select executable")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select an executable file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
							file_format = "INVISIBLE"
							temp_btn = QPushButton("%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "py":
							file_format = cmd_token.type
							temp_btn = QPushButton("Select python file")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip(
								'<FONT>' + "Display open file dialog to select an executable file</FONT>")
							grid_layout.addWidget(temp_btn, grid_row,
												  grid_col_origin + token_label_col_span + token_widget_col_span * 2,
												  token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(
								partial(self.sxcmdwidget.select_file, cmd_token_widget,
										file_format))
							file_format = "INVISIBLE"
							temp_btn = QPushButton("%s" % file_format)
							temp_btn.setToolTip('<FONT>' + "This is %s button</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet(
								'background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row,
												  grid_col_origin + token_label_col_span + token_widget_col_span * 3,
												  token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "dir" or cmd_token.type == "dir_list" or cmd_token.type == "output_continue":
							temp_btn = QPushButton("Select directory")
							temp_btn.setMinimumWidth(func_btn_min_width)
							temp_btn.setToolTip('<FONT>'+"Display select directory dialog"+'</FONT>')
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							temp_btn.clicked.connect(partial(self.sxcmdwidget.select_dir, cmd_token_widget))
							file_format = "INVISIBLE"
							temp_btn = QPushButton("%s" % file_format)
							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
							temp_btn.setEnabled(False)
							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
							temp_btn.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						elif cmd_token.type == "abs_freq":
							# Create button for resolution display. Here, I use button to keep the look & feel.
							cmd_token_subwidget_left = QPushButton("Label for Resolution [A]")
							cmd_token_subwidget_left.setEnabled(False)
							cmd_token_subwidget_left.setMinimumWidth(func_btn_min_width)
							cmd_token_subwidget_left.setToolTip('<FONT>'+"Resolution [A] corresponding to absolute frequency [1/Pixel]"+'</FONT>')
							grid_layout.addWidget(cmd_token_subwidget_left, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
							# Create button to show the associated calculator dialog.
							cmd_token_subwidget_right = QPushButton("Use resolution [A]")
							cmd_token_subwidget_right.setToolTip('<FONT>'+"Display calculator dialog to use the resolution [A] instead of absolute frequency [1/Pixel]. It calculates absolute frequency [1/Pixel] (abs_freq) From resolution [A] (ares) using a give pixel size [A/Pixel] (apix), where abs_freq = apix/ares. </FONT>")
							cmd_token_subwidget_right.setMinimumWidth(func_btn_min_width)
							grid_layout.addWidget(cmd_token_subwidget_right, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
							# Associated this subwidget to open the calculator dialog
							const_register_widget_apix = self.sxcmdwidget.sxconst_set.dict["apix"].register_widget
							cmd_token_calculator_dialog = SXDialogCalculator(const_register_widget_apix, cmd_token_widget, cmd_token_subwidget_left, self)
							cmd_token_calculator_dialog.reflect_external_local_update_abs_freq()
							# Connect the main widget "editing finished" event to the calculator dialog
							cmd_token_widget.editingFinished.connect(cmd_token_calculator_dialog.reflect_external_local_update_abs_freq)
							# Connect the right subwidget "clicked" event to the calculator dialog
							cmd_token_subwidget_right.clicked.connect(cmd_token_calculator_dialog.show)
###						elif cmd_token.type == "any_micrograph":
###							temp_btn = QPushButton("Select Image")
###							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select standard format image file (e.g. .hdf, .mrc)</FONT>")
###							temp_btn.setMinimumWidth(func_btn_min_width)
###							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
###							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget))
###							file_format = "txt"
###							temp_btn = QPushButton("Select text file")
###							temp_btn.setMinimumWidth(func_btn_min_width)
###							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a parameters text file</FONT>" )
###							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
###							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
###						elif cmd_token.type == "any_image":
###							temp_btn = QPushButton("Select Image")
###							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select standard format image file (e.g. .hdf, .mrc)</FONT>")
###							temp_btn.setMinimumWidth(func_btn_min_width)
###							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
###							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget))
###							file_format = "INVISIBLE"
###							temp_btn = QPushButton("%s" % file_format)
###							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
###							temp_btn.setEnabled(False)
###							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
###							temp_btn.setMinimumWidth(func_btn_min_width)
###							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
###						elif cmd_token.type == "any_image_list":
###							temp_btn = QPushButton("Select Images")
###							temp_btn.setMinimumWidth(func_btn_min_width)
###							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select standard format image files (e.g. .hdf, .mrc)</FONT>")
###							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
###							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, cmd_token.type))
###							file_format = "INVISIBLE"
###							temp_btn = QPushButton("%s" % file_format)
###							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
###							temp_btn.setEnabled(False)
###							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
###							temp_btn.setMinimumWidth(func_btn_min_width)
###							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
###						elif cmd_token.type == "any_file":
###							temp_btn = QPushButton("Select File")
###							temp_btn.setMinimumWidth(func_btn_min_width)
###							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select file (e.g. *.*)</FONT>")
###							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
###							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget))
###							file_format = "INVISIBLE"
###							temp_btn = QPushButton("%s" % file_format)
###							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
###							temp_btn.setEnabled(False)
###							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
###							temp_btn.setMinimumWidth(func_btn_min_width)
###							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
###						elif cmd_token.type == "any_file_list":
###							temp_btn = QPushButton("Select Files")
###							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select files (e.g. *.*)</FONT>")
###							temp_btn.setMinimumWidth(func_btn_min_width)
###							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
###							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, cmd_token.type))
###							file_format = "bdb"
###							temp_btn = QPushButton("Select .%s" % file_format)
###							temp_btn.setMinimumWidth(func_btn_min_width)
###							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select .%s format image file</FONT>" % file_format)
###							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
###							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
###						elif cmd_token.type == "image":
###							file_format = "hdf"
###							temp_btn = QPushButton("Select .%s" % file_format)
###							temp_btn.setMinimumWidth(func_btn_min_width)
###							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select .%s format image file</FONT>" % file_format)
###							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
###							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
###							file_format = "bdb"
###							temp_btn = QPushButton("Select .%s" % file_format)
###							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select .%s format image file</FONT>" % file_format)
###							temp_btn.setMinimumWidth(func_btn_min_width)
###							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
###							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
###						elif cmd_token.type == "bdb":
###							file_format = "bdb"
###							temp_btn = QPushButton("Select .%s" % file_format)
###							temp_btn.setMinimumWidth(func_btn_min_width)
###							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select .%s format image file</FONT>" % file_format)
###							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
###							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
###							file_format = "INVISIBLE"
###							temp_btn = QPushButton("%s" % file_format)
###							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
###							temp_btn.setEnabled(False)
###							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
###							temp_btn.setMinimumWidth(func_btn_min_width)
###							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
###						elif cmd_token.type == "hdf":
###							file_format = cmd_token.type
###							temp_btn = QPushButton("Select .%s" % file_format)
###							temp_btn.setMinimumWidth(func_btn_min_width)
###							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select .%s format image file</FONT>" % file_format)
###							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
###							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
###							file_format = "INVISIBLE"
###							temp_btn = QPushButton("%s" % file_format)
###							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
###							temp_btn.setEnabled(False)
###							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
###							temp_btn.setMinimumWidth(func_btn_min_width)
###							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
###						elif cmd_token.type == "mrc":
###							file_format = cmd_token.type
###							temp_btn = QPushButton("Select .%s" % file_format)
###							temp_btn.setMinimumWidth(func_btn_min_width)
###							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select .%s format image file</FONT>" % file_format)
###							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
###							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
###							file_format = "INVISIBLE"
###							temp_btn = QPushButton("%s" % file_format)
###							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
###							temp_btn.setEnabled(False)
###							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
###							temp_btn.setMinimumWidth(func_btn_min_width)
###							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
###						elif cmd_token.type == "txt":
###							file_format = cmd_token.type
###							temp_btn = QPushButton("Select text file")
###							temp_btn.setMinimumWidth(func_btn_min_width)
###							temp_btn.setToolTip('<FONT>'+"Display open file dialog to select a parameters text file</FONT>")
###							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span)
###							self.connect(temp_btn, SIGNAL("clicked()"), partial(self.sxcmdwidget.select_file, cmd_token_widget, file_format))
###							file_format = "INVISIBLE"
###							temp_btn = QPushButton("%s" % file_format)
###							temp_btn.setToolTip('<FONT>'+"This is %s button</FONT>" % file_format)
###							temp_btn.setEnabled(False)
###							temp_btn.setStyleSheet('background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid')
###							temp_btn.setMinimumWidth(func_btn_min_width)
###							grid_layout.addWidget(temp_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 3, token_widget_row_span, token_widget_col_span)
						else:
							if cmd_token.type not in ["int", "float", "string", "output", "apix", "ctfwin", "box", "radius", "sym", "mass", "filament_width"]: ERROR("Logical Error: Encountered unsupported type (%s). Consult with the developer."  % cmd_token.type, "%s in %s" % (__name__, os.path.basename(__file__)))
						
						# if cmd_token.type in ["output", "output_continue", "output_bdb2d_stack"]:
						# 	# Need to add output info button in future

					cmd_token_widget.setToolTip('<FONT>'+cmd_token.help+'</FONT>')
					cmd_token_restore_widget.setToolTip('<FONT>'+cmd_token_restore_tooltip+'</FONT>')

					grid_row += 1
					
				# Register this widget
				cmd_token.restore_widget = cmd_token_restore_widget
				cmd_token.widget = cmd_token_widget
				cmd_token.subwidget_left = cmd_token_subwidget_left
				cmd_token.subwidget_right = cmd_token_subwidget_right
				cmd_token.calculator_dialog = cmd_token_calculator_dialog
		
		is_necessary_to_connet_now = True
		if len(cmd_token_abs_freq_list) == 0:
			is_necessary_to_connet_now = False
		if is_necessary_to_connet_now and len(cmd_token_abs_freq_list) > 0:
			for cmd_token_abs_freq in cmd_token_abs_freq_list:
				if cmd_token_abs_freq.calculator_dialog is None:
					is_necessary_to_connet_now = False
					break
		if is_necessary_to_connet_now and cmd_token_apix is None:
			is_necessary_to_connet_now = False
		if is_necessary_to_connet_now and cmd_token_apix.widget is None:
			is_necessary_to_connet_now = False
		# else: # Do nothing
		
		if is_necessary_to_connet_now:
			# Associated the widget of apix command token to the calculator dialog each other
			assert (cmd_token_apix.type == "apix")
			assert (cmd_token_apix.widget is not None)
			# Loop through all absolute frequency tokens of this command
			# print("MRK_DEBUG: ")
			# print("MRK_DEBUG: ----- SXConstSetWidget Constructor ----- ")
			# print("MRK_DEBUG: self.sxcmdwidget.sxcmd.name          := {}".format(self.sxcmdwidget.sxcmd.name))
			# print("MRK_DEBUG: self.sxcmdwidget.sxcmd.subname       := {}".format(self.sxcmdwidget.sxcmd.subname))
			# print("MRK_DEBUG: self.sxcmdwidget.sxcmd.mode          := {}".format(self.sxcmdwidget.sxcmd.mode))
			# print("MRK_DEBUG: self.sxcmdwidget.sxcmd.subset_config := {}".format(self.sxcmdwidget.sxcmd.subset_config))
			for cmd_token_abs_freq in cmd_token_abs_freq_list:
				assert (cmd_token_abs_freq.type == "abs_freq")
				# Register the calculator dialogs of all the other absolute frequency tokens
				# by looping through all the other absolute frequency tokens
				for other_cmd_token_abs_freq in cmd_token_abs_freq_list:
					assert (other_cmd_token_abs_freq.type == "abs_freq")
					# Exclude itself
					if cmd_token_abs_freq.key_base != other_cmd_token_abs_freq.key_base:
						cmd_token_abs_freq.other_dialog_list.append(other_cmd_token_abs_freq.calculator_dialog)
						cmd_token_abs_freq.calculator_dialog.sxcmd_token_other_dialog_list_abs_freq.append(other_cmd_token_abs_freq.calculator_dialog)
				# Register pixel size token widget to this absolut frequency calculator dialog
				# if cmd_token_abs_freq is None:
				# 	print("MRK_DEBUG: cmd_token_abs_freq is None")
				# if cmd_token_abs_freq.calculator_dialog is None:
				# 	print("MRK_DEBUG: cmd_token_abs_freq.calculator_dialog is None")
				# 	print("MRK_DEBUG: cmd_token_abs_freq.key_base := {}".format(cmd_token_abs_freq.key_base))
				# 	print("MRK_DEBUG: cmd_token_abs_freq.type     := {}".format(cmd_token_abs_freq.type))
				cmd_token_abs_freq.calculator_dialog.sxcmd_token_widget_apix = cmd_token_apix.widget
				# Register this absolut frequency calculator dialog to pixel size token
				cmd_token_apix.other_dialog_list.append(cmd_token_abs_freq.calculator_dialog)
			
			# print("MRK_DEBUG: len(cmd_token_apix.other_dialog_list) := \"{}\"".format(len(cmd_token_apix.other_dialog_list)))
			
			# Initialise pixel size of all calculate dialogs
			self.handle_apix_token_editing_finished_event(cmd_token_apix)
			# Connect the apix command token widget "editing finished" event to the calculator dialog
			cmd_token_apix.widget.editingFinished.connect(partial(self.handle_apix_token_editing_finished_event, cmd_token_apix))
		
		if tab_group == "main":
			# Add space
			grid_row += 1

			# Add gui components for MPI related parameters
			temp_label = QLabel("MPI processors")
			temp_label.setMinimumWidth(token_label_min_width)
			submit_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)

			# self.mpi_nproc_edit = QLineEdit(self)
			self.mpi_nproc_edit = QLineEdit()
			self.mpi_nproc_edit.setText("1")
			self.mpi_nproc_edit.setToolTip('<FONT>'+"Number of processors to use. default is single processor mode"+'</FONT>')
			submit_layout.addWidget(self.mpi_nproc_edit, grid_row, grid_col_origin + token_label_col_span, token_widget_row_span, token_widget_col_span)

			# Add save parameters button
			self.save_params_btn = QPushButton("Save parameters")
			self.save_params_btn.setMinimumWidth(btn_min_width)
			self.save_params_btn.setToolTip('<FONT>'+"Save gui parameter settings"+'</FONT>')
			self.save_params_btn.clicked.connect(self.sxcmdwidget.save_params)
			submit_layout.addWidget(self.save_params_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span*2)

			grid_row += 1

			temp_label = QLabel("MPI command line template")
			temp_label.setMinimumWidth(token_label_min_width)
			submit_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)

			self.mpi_cmd_line_edit = QLineEdit()
			if 'SPHIRE_MPI_COMMAND_LINE_TEMPLATE' in os.environ:
				mpi_command_line_template = os.environ['SPHIRE_MPI_COMMAND_LINE_TEMPLATE']
			else:
				mpi_command_line_template = ""
			self.mpi_cmd_line_edit.setText(mpi_command_line_template)
			self.mpi_cmd_line_edit.setToolTip('<FONT>'+"Template of MPI command line (e.g. \"mpirun -np XXX_SXMPI_NPROC_XXX --host n0,n1,n2 XXX_SXCMD_LINE_XXX\"). if empty, use \"mpirun -np XXX_SXMPI_NPROC_XXX XXX_SXCMD_LINE_XXX\"</FONT>")
			submit_layout.addWidget(self.mpi_cmd_line_edit, grid_row, grid_col_origin + token_label_col_span, token_widget_row_span, token_widget_col_span)

			# Add load parameters button
			self.load_params_btn = QPushButton("Load parameters")
			self.load_params_btn.setMinimumWidth(btn_min_width)
			self.load_params_btn.setToolTip('<FONT>'+"Load gui parameter settings to retrieve a previously-saved one"+'</FONT>')
			self.load_params_btn.clicked.connect(self.sxcmdwidget.load_params)
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
				if 'SPHIRE_SUBMISSION_COMMAND' in os.environ:
					submission_command = os.environ['SPHIRE_SUBMISSION_COMMAND']
				else:
					submission_command = "qsub"
				self.qsub_cmd_edit.setText(submission_command)
			else: # assert(self.sxcmdwidget.sxcmd.is_submittable == False)
				assert(self.sxcmdwidget.sxcmd.mpi_support == False)
				self.qsub_cmd_edit.setText("N/A")
			self.qsub_cmd_edit.setToolTip('<FONT>'+"Name of submission command to queue job"+'</FONT>')
			submit_layout.addWidget(self.qsub_cmd_edit, grid_row, grid_col_origin + token_label_col_span, token_widget_row_span, token_widget_col_span)

			self.cmd_line_btn = QPushButton("Generate command line")
			self.cmd_line_btn.setMinimumWidth(btn_min_width)
			self.cmd_line_btn.setToolTip('<FONT>'+"Generate command line from gui parameter settings and automatically save settings"+'</FONT>')
			self.cmd_line_btn.clicked.connect(self.sxcmdwidget.print_cmd_line)
			submit_layout.addWidget(self.cmd_line_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span*2)

			grid_row += 1

			self.pipe_line_btn = QPushButton("Add to pipeline folder")
			self.pipe_line_btn.setMinimumWidth(btn_min_width)
			self.pipe_line_btn.setToolTip('<FONT>'+"Generate executable files that and add them to the queue folder."+'</FONT>')
			self.pipe_line_btn.clicked.connect(self.sxcmdwidget.add_to_pipeline)
			submit_layout.addWidget(self.pipe_line_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span*2)

			grid_row += 1

			temp_label = QLabel("Submission script template")
			temp_label.setMinimumWidth(token_label_min_width)
			submit_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)

			self.qsub_script_edit = QLineEdit()
			if self.sxcmdwidget.sxcmd.is_submittable == True:
				if 'SPHIRE_SUBMISSION_SCRIPT_TEMPLATE' in os.environ:
					submission_script_template = os.environ['SPHIRE_SUBMISSION_SCRIPT_TEMPLATE']
				else:
					submission_script_template = "msgui_qsub.sh"
				self.qsub_script_edit.setText(submission_script_template)
			else: # assert(self.sxcmdwidget.sxcmd.is_submittable == False)
				assert(self.sxcmdwidget.sxcmd.mpi_support == False)
				self.qsub_script_edit.setText("N/A")
			self.qsub_script_edit.setToolTip('<FONT>'+"File name of submission script template (e.g. $PROJECT_DIR/msgui_qsub.sh)"+'</FONT>')
			submit_layout.addWidget(self.qsub_script_edit, grid_row, grid_col_origin + token_label_col_span, token_widget_row_span, token_widget_col_span)

			self.qsub_script_open_btn = QPushButton("Select template")
			self.qsub_script_open_btn.setMinimumWidth(func_btn_min_width)
			self.qsub_script_open_btn.setToolTip('<FONT>'+"Display open file dialog to select job submission script template file"+'</FONT>')
			self.qsub_script_open_btn.clicked.connect(partial(self.sxcmdwidget.select_file, self.qsub_script_edit, 'submission_template'))
			submit_layout.addWidget(self.qsub_script_open_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span, token_widget_row_span, token_widget_col_span)

			# Add a run button
			# self.execute_btn = QPushButton("Run %s" % self.sxcmdwidget.sxcmd.get_mode_name_for("human"))
			self.execute_btn = QPushButton("Run command")
			# make 3D textured push button look
			custom_style = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0)} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084)} QPushButton:focus {font: bold; color: #000;border: 2px solid #8D0;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0)} QPushButton:disabled {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #ff0000)}"
			self.execute_btn.setStyleSheet(custom_style)
			self.execute_btn.setMinimumWidth(btn_min_width)
			self.execute_btn.setToolTip('<FONT>'+"Run %s and automatically save gui parameter settings</FONT>" % self.sxcmdwidget.sxcmd.get_mode_name_for("human"))
			self.execute_btn.clicked.connect(self.sxcmdwidget.execute_cmd_line)
			submit_layout.addWidget(self.execute_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span*2)

			grid_row += 1

###			# Add a debug button
###			sxoperand_set = self.construct_sxoperand_set()
###			self.res_calculator = SXDialogCalculator(sxoperand_set, self)
###			self.cmd_line_btn = QPushButton("Debug")
###			self.cmd_line_btn.setMinimumWidth(btn_min_width)
###			self.cmd_line_btn.setToolTip('<FONT>'+"Button for debug"+'</FONT>')
###			self.connect(self.cmd_line_btn, SIGNAL("clicked()"), self.res_calculator.show)
###			submit_layout.addWidget(self.cmd_line_btn, grid_row, grid_col_origin + token_label_col_span + token_widget_col_span * 2, token_widget_row_span, token_widget_col_span*2)
###
###			grid_row += 1

			# Initialize enable state of qsub related widgets
			self.set_qsub_enable_state()

	@staticmethod
	def set_text_entry_widget_enable_state(widget, is_enabled):
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

	def handle_restore_widget_event(self, sxcmd_token, widget_index=0, restore_idx=0):
		assert(not sxcmd_token.is_locked)
		if sxcmd_token.type == "user_func":
			assert(len(sxcmd_token.widget) == 2 and len(sxcmd_token.restore) == 2 and widget_index < 2)
			sxcmd_token.widget[widget_index].setText("%s" % sxcmd_token.restore[widget_index][restore_idx][0])
		else:
			if sxcmd_token.type in ("bool", "bool_ignore"):
				try:
					text = str(self.sender().text())
				except AttributeError:
					if restore_idx == 0:
						text = str(self.sender().currentText())
					else:
						text = str(self.sender().itemText(restore_idx))
					text = text.strip()
					self.sender().blockSignals(True)
					self.sender().setCurrentIndex(0)
					self.sender().blockSignals(False)
				if text == 'YES':
					sxcmd_token.widget.setChecked(Qt.Checked)
				else: # sxcmd_token.restore == False
					sxcmd_token.widget.setChecked(Qt.Unchecked)
			else:
				try:
					text = str(self.sender().text())
				except AttributeError:
					if restore_idx == 0:
						text = str(self.sender().currentText())
					else:
						text = str(self.sender().itemText(restore_idx))
					text = text.strip()
					self.sender().blockSignals(True)
					self.sender().setCurrentIndex(0)
					self.sender().blockSignals(False)
				sxcmd_token.widget.setText(text)
				if sxcmd_token.type == "abs_freq":
					sxcmd_token.calculator_dialog.reflect_external_local_update_abs_freq()
					for sxcmd_token_other_dialog in sxcmd_token.other_dialog_list:
						sxcmd_token_other_dialog.reflect_external_local_update_abs_freq()
				elif sxcmd_token.type == "apix":
					for sxcmd_token_other_dialog in sxcmd_token.other_dialog_list:
						sxcmd_token_other_dialog.reflect_external_local_update_apix()


	def handle_apix_token_editing_finished_event(self, sxcmd_token_apix):
		assert (sxcmd_token_apix.type == "apix")
		# print("MRK_DEBUG: ")
		# print("MRK_DEBUG: ----- SXCmdTab.handle_apix_token_editing_finished_event() ----- ")
		# print("MRK_DEBUG: len(sxcmd_token_apix.other_dialog_list) := \"{}\"".format(len(sxcmd_token_apix.other_dialog_list)))
		for cmd_token_calculator_dialog_abs_freq in sxcmd_token_apix.other_dialog_list:
			# print("MRK_DEBUG: cmd_token_calculator_dialog_abs_freq := \"{}\"".format(cmd_token_calculator_dialog_abs_freq))
			if cmd_token_calculator_dialog_abs_freq is not None:
				cmd_token_calculator_dialog_abs_freq.reflect_external_local_update_apix()

###	def handle_abs_freq_editing_finished_event(self, sxcmd_token_widget_abs_freq, sxcmd_token_subwidget_left_ares):
###		apix_str = self.sxcmdwidget.sxconst_set.dict["apix"].register
###		abs_freq_str = sxcmd_token_widget_abs_freq.text()
###		ares_str, is_valid_ares = SXDialogCalculator.convert_abs_freq_to_ares(apix_str, abs_freq_str)
###		print("MRK_DEBUG: ----- handle_abs_freq_editing_finished_event() ----- ")
###		print("MRK_DEBUG: Input Abs. Freq. [1/Pixel] string ; abs_freq_str := {}".format(abs_freq_str))
###		print("MRK_DEBUG: Input Pixel Size [A/Pixel] string ; apix_str := {}".format(apix_str))
###		print("MRK_DEBUG: Output Resolution [A] string      ; ares_str:= {}".format(ares_str))
###		if is_valid_ares:
###			sxcmd_token_subwidget_left_ares.setText("{}[A]@{}[A/Pix]".format(ares_str, apix_str))
###		else:
###			sxcmd_token_subwidget_left_ares.setText("Mode {}".format(ares_str))


# ========================================================================================
# Command Category Widget (opened by class SXMainWindow)
class SXCmdCategoryWidget(QWidget):
	def __init__(self, sxconst_set, sxcmd_category, helical=False, parent = None):
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
		self.add_sxcmd_widgets(helical)

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
		self.grid_layout.setContentsMargins(SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin)
		self.grid_layout.setSpacing(SXLookFeelConst.grid_spacing)
		self.grid_layout.setColumnMinimumWidth(0, SXLookFeelConst.sxcmd_btn_area_min_width)
		# self.grid_layout.setColumnMinimumWidth(1, SXLookFeelConst.sxcmd_widget_area_min_width)
		# Give the column of the command settings area a higher stretch priority so that the other area does not stretch horizontally
		# self.grid_layout.setColumnStretch(self.grid_col_origin + self.sxcmd_btn_area_col_span, self.grid_layout.columnStretch(self.grid_col_origin + self.sxcmd_btn_area_col_span) + 1)

	# Add Pipeline SX Commands (sx*.py) associated widgets
	def add_sxcmd_widgets(self, helical):
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
			sxcmd.widget = SXCmdWidget(self.sxconst_set, sxcmd, helical, self)
			self.stacked_layout.addWidget(sxcmd.widget)

			# connect widget signals
			sxcmd.btn.clicked.connect(partial(self.handle_sxcmd_btn_event, sxcmd))

			self.grid_row += 1

	def load_previous_session(self):
		for sxcmd in self.sxcmd_category.cmd_list:
			if os.path.exists(sxcmd.widget.gui_settings_file_path):
				sxcmd.widget.read_params(sxcmd.widget.gui_settings_file_path)

	def handle_sxcmd_btn_event(self, sxcmd):
		modifiers = QApplication.keyboardModifiers()
		if modifiers == Qt.ShiftModifier:
			# os.system("python -m webbrowser %s%s" % (SPARX_DOCUMENTATION_WEBSITE, sxcmd.name))
			# sxcmd_wiki_url = SXLookFeelConst.generate_sxcmd_wiki_url(sxcmd, wiki_type = "SPARX")
			sxcmd_wiki_url = SXLookFeelConst.generate_sxcmd_wiki_url(sxcmd)
			print("Opening Wiki Page ...")
			print(sxcmd_wiki_url)
			os.system("python -m webbrowser %s" % (sxcmd_wiki_url))
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
	def __init__(self, sxconst_set, sxcmd_category_list, helical=False, parent=None):
		super(SXConstSetWidget, self).__init__(parent)

		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.parent = parent
		self.sxconst_set = sxconst_set
		self.sxcmd_category_list = sxcmd_category_list

		self.gui_settings_file_path = "%s/gui_settings_project_%d.txt" % (SXLookFeelConst.project_dir, int(helical))

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
		global_layout.setContentsMargins(SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin)
		global_layout.setSpacing(SXLookFeelConst.grid_spacing)
		global_layout.setRowStretch(global_row_span - 1, global_layout.rowStretch(global_row_origin) + 1)

		header_layout = QGridLayout()
		header_layout.setContentsMargins(SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin)
		header_layout.setSpacing(SXLookFeelConst.grid_spacing)

		const_set_layout = QGridLayout()
		const_set_layout.setContentsMargins(SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin)
		const_set_layout.setSpacing(SXLookFeelConst.grid_spacing)

		btn_layout = QGridLayout()
		btn_layout.setContentsMargins(SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin)
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
			sxconst_register_widget.clicked.connect(partial(self.handle_regster_widget_event, sxconst))

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
		self.execute_btn.clicked.connect(self.register_const_set)
		btn_layout.addWidget(self.execute_btn, btn_grid_row, btn_col_origin, register_btn_row_span, register_btn_col_span)

		btn_grid_row += 1

		# Add save project constant parameter settings button
		self.save_consts_btn = QPushButton("Save settings")
		self.save_consts_btn.setMinimumWidth(func_btn_min_width)
		self.save_consts_btn.setToolTip('<FONT>'+"Save project constant parameter settings"+'</FONT>')
		self.save_consts_btn.clicked.connect(self.save_consts)
		btn_layout.addWidget(self.save_consts_btn, btn_grid_row, btn_col_origin, func_btn_row_span, func_btn_col_span)

		# Add load project constant parameter settings button
		self.load_consts_btn = QPushButton("Load settings")
		self.load_consts_btn.setMinimumWidth(func_btn_min_width)
		self.load_consts_btn.setToolTip('<FONT>'+"Load project constant parameter settings to retrieve the previously-saved one"+'</FONT>')
		self.load_consts_btn.clicked.connect(self.load_consts)
		btn_layout.addWidget(self.load_consts_btn, btn_grid_row, btn_col_origin + func_btn_col_span, func_btn_row_span, func_btn_col_span)

		btn_grid_row += 1
		btn_layout.addWidget(QLabel(), btn_grid_row, btn_col_origin, register_btn_row_span, register_btn_col_span)

		btn_grid_row += 1

		# Add a register button
		custom_style = "QPushButton {color: #68a2c3}"
		if helical:
			self.helical_spa_btn = QPushButton("Switch to SPA GUI")
			self.helical_spa_btn.setToolTip('Switch to the SPA GUI')
			self.helical_spa_btn.setStyleSheet(custom_style)
		else:
			self.helical_spa_btn = QPushButton("Switch to helical GUI")
			self.helical_spa_btn.setToolTip('Switch to the helical GUI')
			self.helical_spa_btn.setStyleSheet(custom_style)
		# make 3D textured push button look
		self.helical_spa_btn.setMinimumWidth(func_btn_min_width * register_btn_col_span)
		self.helical_spa_btn.clicked.connect(self.switch_gui)
		btn_layout.addWidget(self.helical_spa_btn, btn_grid_row, btn_col_origin, register_btn_row_span, register_btn_col_span)

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

	def switch_gui(self):
		self.parent.helical = not self.parent.helical
		self.parent.sxmenu_item_widget_stacked_layout_global.setCurrentIndex(self.parent.helical)
		self.parent.sxmenu_item_widget_stacked_layout[self.parent.helical].setCurrentWidget(self.parent.sxconst_set[self.parent.helical].widget)

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
					if not cmd_token.is_locked and cmd_token.type in list(self.sxconst_set.dict.keys()):
						sxconst = self.sxconst_set.dict[cmd_token.type]
						cmd_token.restore = sxconst.register
						try:
							cmd_token.restore_widget.setText("%s" % cmd_token.restore)
						except:
							cmd_token.restore_widget.setCurrentIndex(0)
						cmd_token.widget.setText(cmd_token.restore)
						# print "MRK_DEBUG: %s, %s, %s, %s, %s, %s" % (sxcmd.name, sxcmd.subname, cmd_token.key_base, cmd_token.type, cmd_token.default, cmd_token.restore)
					elif cmd_token.type == "abs_freq":
						assert("apix" in list(self.sxconst_set.dict.keys()))
###						print("MRK_DEBUG: ----- register_const_set() ----- ")
###						print("MRK_DEBUG: cmd_token.type = {}".format(cmd_token.type))
###						sxconst_apix = self.sxconst_set.dict["apix"]
###						print("MRK_DEBUG: sxconst_apix.type = {}".format(sxconst_apix.type))
###						print("MRK_DEBUG: sxconst_apix.register = {}".format(sxconst_apix.register))
###						assert (cmd_token.subwidget_left is not None)
###						assert (cmd_token.subwidget_right is not None)
						assert (cmd_token.calculator_dialog is not None)
###						sxoperand_apix = cmd_token.calculator_dialog.sxoperand_set.dict["apix"]
###						print("MRK_DEBUG: BEFORE sxoperand_apix.register = {}".format(sxoperand_apix.register))
###						print("MRK_DEBUG: BEFORE sxoperand_apix.register_widget.text() = {}".format(sxoperand_apix.register_widget.text()))
###						print("MRK_DEBUG: BEFORE sxoperand_apix.widget.text() = {}".format(sxoperand_apix.widget.text()))
###						sxoperand_apix.register = sxconst_apix.register
###						sxoperand_apix.register_widget.setText("%s" % sxoperand_apix.register)
###						sxoperand_apix.widget.setText(sxoperand_apix.register)
###						print("MRK_DEBUG: AFTER sxoperand_apix.register = {}".format(sxoperand_apix.register))
###						print("MRK_DEBUG: AFTER sxoperand_apix.register_widget.text() = {}".format(sxoperand_apix.register_widget.text()))
###						print("MRK_DEBUG: AFTER sxoperand_apix.widget.text() = {}".format(sxoperand_apix.widget.text()))
						cmd_token.calculator_dialog.reflect_external_global_update_apix()

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
				if key not in list(self.sxconst_set.dict.keys()):
					QMessageBox.warning(self, "Invalid Project Settings File Format", "Invalid entry key for project settings \"%s\" is found in line (%s). This project settings file might be incompatible with the current version. Please save the project settings file again." % (key, line_in))
				sxconst = self.sxconst_set.dict[key]
				sxconst.widget.setText(val_str_in)
		else:
			QMessageBox.warning(self, "Fail to load project settings", "The specified file is not project settings file.")

		file_in.close()

	def save_consts(self):
		name = QFileDialog.getSaveFileName(self, "Save settings", SXLookFeelConst.file_dialog_dir, options = QFileDialog.DontUseNativeDialog)
		if isinstance(name, tuple):
			file_path = str(name[0])
		else:
			file_path = str(name)
		if file_path != "":
			self.write_consts(file_path)

	def load_consts(self):
		name = QFileDialog.getOpenFileName(self, "Load settings", SXLookFeelConst.file_dialog_dir, options = QFileDialog.DontUseNativeDialog)
		if isinstance(name, tuple):
			file_path = str(name[0])
		else:
			file_path = str(name)
		if file_path != "":
			self.read_consts(file_path)

# ========================================================================================
# Layout of the information widget; owned by the main window
class SXInfoWidget(QWidget):
	def __init__(self, parent = None):
		super(SXInfoWidget, self).__init__(parent)

		self.setStyleSheet("background-color: {0}".format(SXLookFeelConst.default_bg_color_string))
		widget = QWidget(self)

		# Get the picture name
		pic_name = '{0}sxgui_info.png'.format(get_image_directory())
		# Import the picture as pixmap to get the right dimensions
		self.pixmap = QPixmap(pic_name)
		width = self.pixmap.width()
		height = self.pixmap.height()

		# Scrol widget
		scroll_widget = QWidget()
		scroll_widget.setStyleSheet('background-color: transparent')

		label1 = QLabel()
		label1.setFixedHeight(40)
		label2 = QLabel()
		label2.setFixedHeight(40)

		# Create a QLabel and show the picture
		self.label = QLabel()
		self.label.setFixedSize(width, height)
		self.label.setStyleSheet('border-image: url("{0}"); background-color: transparent'.format(pic_name))

		# Layout for the scroll widet vert
		label3 = QLabel()
		label3.setFixedWidth(40)
		label4 = QLabel()
		label4.setFixedWidth(40)

		# Layout for the scroll widget hor
		layout_vert = QHBoxLayout()
		layout_vert.addWidget(label3)
		layout_vert.addWidget(self.label)
		layout_vert.addWidget(label4)

		# Layout for the scroll widget hor
		layout = QVBoxLayout(scroll_widget)
		layout.addWidget(label1)
		layout.addLayout(layout_vert)
		layout.addWidget(label2)

		# Add a scroll area for vertical scrolling
		scroll_area = QScrollArea(widget)
		scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
		scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
		scroll_area.setWidget(scroll_widget)
		scroll_area.setStyleSheet("background-color: {0}".format(SXLookFeelConst.sxcmd_tab_bg_color_string))

		layout = QHBoxLayout(widget)
		layout.addWidget(scroll_area, stretch=1)

		layout = QHBoxLayout(self)
		layout.addWidget(widget)
		layout.setContentsMargins(0, 0, 0, 0)

"""
# ========================================================================================
class SXDialogCalculator(QDialog):
	def __init__(self, parent = None):
		super(QDialog, self).__init__(parent)
		
		self.setWindowModality(Qt.ApplicationModal)
		
		# self.setWindowTitle()
		self.setWindowTitle("Absolute Frequency Calculator")

		temp_label = QLabel("Calculate absolute frequency [1/Pixel] from resolution [A]", self)
		temp_label.move(50,50)
		
		# Create label widget
		temp_label = QLabel("Resolution [A]", self)
		# temp_label.setMinimumWidth(token_label_min_width)
		# grid_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)
		temp_label.move(50,100)
		self.edit_res = QLineEdit(self)
		self.edit_res.setText('Enter Resolution Here')
		self.edit_res.move(200,100)

		temp_label = QLabel("Pixel Size [A/Pixel]", self)
		# temp_label.setMinimumWidth(token_label_min_width)
		# grid_layout.addWidget(temp_label, grid_row, grid_col_origin, token_label_row_span, token_label_col_span)
		temp_label.move(50,200)
		self.edit_apix = QLineEdit(self)
		self.edit_apix.setText('Enter Pixel Size Here')
		self.edit_apix.move(200,200)
		
		self.btn_apply = QPushButton("Apply", self)
		self.btn_apply.move(50,300)
		self.btn_cancel = QPushButton("Cancel", self)
		self.btn_cancel.move(200,300)
		# self.connect(cmd_token_restore_widget[widget_index], SIGNAL("clicked()"), partial(self.handle_restore_widget_event, cmd_token, widget_index))
		
		### self.show()
"""

class SXDialogCalculator(QDialog):
	# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	# static class variables
	# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	precision_apix = 7
	precision_abs_freq = 7
	precision_ares = 2

	def __init__(self, sxconst_register_widget_apix, sxcmd_token_widget_abs_freq, sxcmd_token_subwidget_left_ares, parent = None):
		super(QDialog, self).__init__(parent)
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.sxoperand_set = self.construct_sxoperand_set()
		self.sxconst_register_widget_apix = sxconst_register_widget_apix # SXconst.register_widget of apix type
		self.sxcmd_token_widget_abs_freq = sxcmd_token_widget_abs_freq
		self.sxcmd_token_subwidget_left_ares = sxcmd_token_subwidget_left_ares
		# This have to be set upon the construction order...
		self.sxcmd_token_widget_apix = None                              # SXcmd_token.widget of apix type if this command has one
		self.sxcmd_token_other_dialog_list_abs_freq = []                 # SXcmd_token.other_dialog_list of the associated abs_freq type token if this command has more than one abs_freq type tokens
		
		# This should be a modal dialog
		self.setWindowModality(Qt.ApplicationModal)
		self.setWindowTitle("%s" % (self.sxoperand_set.label))

		# Layout constants and variables
		global_row_origin = 0; global_col_origin = 0
		global_row_span = 4; global_col_span = 1

		header_row_origin = 0; header_col_origin = 0
		title_row_span = 1; title_col_span = 1
		short_info_row_span = 1; short_info_col_span = 1
		title_min_width = 300
		short_info_min_width = 300
		short_info_min_height = 80

		operand_set_row_origin = 0; operand_set_col_origin = 0
		operand_label_row_span = 1; operand_label_col_span = 1
		operand_register_widget_row_span = 1; operand_register_widget_col_span = 1
		operand_widget_row_span = 1; operand_widget_col_span = 1
		operand_label_min_width = 150
		operand_register_widget_min_width = operand_label_min_width
		operand_widget_min_width = operand_label_min_width

		btn_row_origin = 0; btn_col_origin = 0
		func_btn_row_span = 1; func_btn_col_span = 1
		register_btn_row_span = 1; register_btn_col_span = 2
		func_btn_min_width = 50
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# background_image_file_path = '{0}sxgui_background.png'.format(get_image_directory())
		# self.setStyleSheet('background-image: url("{0}")'.format(background_image_file_path))

		# Set the background color of this widget
		self.setAutoFillBackground(True)
		palette = QPalette()
		palette.setBrush(QPalette.Background, QBrush(SXLookFeelConst.default_bg_color_calculator))
		self.setPalette(palette)

		global_layout = QGridLayout()
		global_layout.setContentsMargins(SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin)
		global_layout.setSpacing(SXLookFeelConst.grid_spacing)
		global_layout.setRowStretch(global_row_span - 1, global_layout.rowStretch(global_row_origin) + 1)

#		header_layout = QGridLayout()
#		header_layout.setContentsMargins(SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin)
#		header_layout.setSpacing(SXLookFeelConst.grid_spacing)

		operand_layout = QGridLayout()
		operand_layout.setContentsMargins(SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin)
		operand_layout.setSpacing(SXLookFeelConst.grid_spacing)

		btn_layout = QGridLayout()
		btn_layout.setContentsMargins(SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin, SXLookFeelConst.grid_margin)
		btn_layout.setSpacing(SXLookFeelConst.grid_spacing * 2)

		global_grid_row = global_row_origin

#		# Start add title widgets to the grid layout
#		header_grid_row = header_row_origin
#
#		# Set a label and its position in this dialog
#		temp_label = QLabel("<b>%s</b>" % (self.sxoperand_set.label))
#		temp_label.setMinimumWidth(title_min_width)
#		header_layout.addWidget(temp_label, header_grid_row, header_col_origin, title_row_span, title_col_span)
#
#		header_grid_row += 1

#		# NOTE: 2015/11/17 Toshio Moriya
#		# Necessary to separate "<b>%s</b>" from the information for avoiding to invoke the tag interpretations of string
#		# e.g. < becomes the escape character
#		temp_label = QLabel("%s" % (self.sxoperand_set.short_info))
#		temp_label.setWordWrap(True)
#		temp_label.setMinimumWidth(short_info_min_width)
#		temp_label.setMinimumHeight(short_info_min_height)
#		header_layout.addWidget(temp_label, header_grid_row, header_col_origin, short_info_row_span, short_info_col_span)
#
#		# Add const set grid layout to global layout
#		global_layout.addLayout(header_layout, global_grid_row, global_col_origin)
#		global_grid_row += 1

		# Start add project parameter constant widgets to the grid layout
		operand_grid_row = operand_set_row_origin

		# Add widget for editing command args and options
		for sxoperand in self.sxoperand_set.list:
			# Create widget associated to this project constant parameter
			temp_label = QLabel(sxoperand.label)
			temp_label.setMinimumWidth(operand_label_min_width)
			operand_layout.addWidget(temp_label, operand_grid_row, operand_set_col_origin, operand_label_row_span, operand_label_col_span)

			sxoperand_register_widget = QPushButton("%s" % sxoperand.register)
			sxoperand_register_widget.setMinimumWidth(operand_register_widget_min_width)
			custom_style = ""
			if sxoperand.is_input:
				if sxoperand.type in ["apix"]:
					custom_style = "QPushButton {color:green; }" # custom_style = "QPushButton {font: normal; color:black;}"
					sxoperand_register_widget.setStyleSheet(custom_style)
				sxoperand_register_widget.setToolTip('<FONT>'+"Retrieve this default value to edit box"+'</FONT>')
				sxoperand_register_widget.clicked.connect(partial(self.handle_regster_widget_event, sxoperand))
			else:
				assert (not sxoperand.is_input)
				custom_style = "QPushButton {background: rgba(0, 0, 0, 0); color: rgba(0, 0, 0, 0); border: 0px rgba(0, 0, 0, 0) solid;}"
				sxoperand_register_widget.setStyleSheet(custom_style)
				sxoperand_register_widget.setEnabled(False)
			operand_layout.addWidget(sxoperand_register_widget, operand_grid_row, operand_set_row_origin + operand_label_col_span, operand_register_widget_row_span, operand_register_widget_col_span)

			sxoperand_widget = QLineEdit()
			sxoperand_widget.setMinimumWidth(operand_widget_min_width)
			sxoperand_widget.setText(sxoperand.register)
			sxoperand_widget.setToolTip('<FONT>'+sxoperand.help+'</FONT>')
			operand_layout.addWidget(sxoperand_widget, operand_grid_row, operand_set_row_origin + operand_label_col_span + operand_register_widget_col_span, operand_widget_row_span, operand_widget_col_span)
			if not sxoperand.is_input:
				sxoperand_widget.setEnabled(False)
				
			operand_grid_row += 1

			# Register this widget
			sxoperand.register_widget = sxoperand_register_widget
			sxoperand.widget = sxoperand_widget

		# Add const set grid layout to global layout
		global_layout.addLayout(operand_layout, global_grid_row, global_col_origin)
		# global_grid_row += 1

		# Start add buttons to the grid layout
		btn_grid_row = btn_row_origin

		# Add save project constant parameter settings button
		self.convert_units_btn = QPushButton("Convert units")
		self.convert_units_btn.setMinimumWidth(func_btn_min_width)
		self.convert_units_btn.setToolTip('<FONT>'+"Calculate unit conversion"+'</FONT>')
		self.convert_units_btn.setDefault(True);
		self.convert_units_btn.clicked.connect(self.handle_convert_units)
		btn_layout.addWidget(self.convert_units_btn, btn_grid_row, btn_col_origin, func_btn_row_span, func_btn_col_span)

		# Add load project constant parameter settings button
		self.cancel_btn = QPushButton("Cancel")
		self.cancel_btn.setMinimumWidth(func_btn_min_width)
		self.cancel_btn.setToolTip('<FONT>'+"Cancel the unit convertion. Close the calculator dialog without applying the converted value to the associated edit box"+'</FONT>')
		self.cancel_btn.clicked.connect(self.close)
		btn_layout.addWidget(self.cancel_btn, btn_grid_row, btn_col_origin + func_btn_col_span, func_btn_row_span, func_btn_col_span)

		btn_grid_row += 1

###		# For now, initialise with Nyquist resolution of 1.0 [A] pixel size
###		self.apix_str = self.sxconst_register_widget_apix.text()
###		self.abs_freq_str = self.sxcmd_token_widget_abs_freq.text()
###		self.convert_abs_freq_to_ares()
###		sxoperand_apix = self.sxoperand_set.dict["apix"]
###		sxoperand_apix.register = apix_str
###		sxoperand_apix.widget.setText(sxoperand_apix.register)
###		sxoperand_ares = self.sxoperand_set.dict["ares"]
###		sxoperand_ares.register = ares_str
###		sxoperand_ares.widget.setText(sxoperand_ares.register)
		
		# Add a apply button
		self.apply_btn = QPushButton("Apply")
		# make 3D textured push button look
		custom_style = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0)} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084)} QPushButton:focus {font: bold; color: #000;border: 2px solid #8D0;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0)} QPushButton:disabled {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #ff0000)}"
		self.apply_btn.setStyleSheet(custom_style)
		self.apply_btn.setMinimumWidth(func_btn_min_width * register_btn_col_span)
		self.apply_btn.setToolTip('<FONT>'+"Apply converted value to the corresponding command arguments and options"+'</FONT>')
		self.apply_btn.clicked.connect(self.handle_apply_unit_conversion)
		btn_layout.addWidget(self.apply_btn, btn_grid_row, btn_col_origin, register_btn_row_span, register_btn_col_span)

		btn_grid_row += 1

		# Add button grid layout to global layout
		# global_layout.addLayout(btn_layout, global_grid_row, global_col_origin) # Maybe later :)

#		# Load the previously saved parameter setting of this sx command
#		if os.path.exists(self.gui_settings_file_path):
#			self.read_consts(self.gui_settings_file_path)

		# Layout for a constant size
		constant_height_layout = QVBoxLayout()
		constant_height_layout.addLayout(global_layout)
		constant_height_layout.addLayout(btn_layout)
		constant_height_layout.addStretch(1)
		constant_width_layout = QHBoxLayout(self)
		constant_width_layout.addLayout(constant_height_layout)
		constant_width_layout.addStretch(1)
		
		# Synchronise the initial status of operands in this dialog
		self.handle_convert_units(is_enable_message = False)
		
	def construct_sxoperand_set(self):
		sxoperand_set = SXoperand_set(); sxoperand_set.name = "sxc_calculator"; sxoperand_set.label = "Absolute Frequency Calculator"; sxoperand_set.short_info = "From resolution [A] (ares), calculate absolute frequency [1/Pixel] (abs_freq) using a give pixel size [A/Pixel] (apix), where abs_freq = apix/ares."
		sxoperand = SXoperand(); sxoperand.key = "ares"; sxoperand.label = "Resolution [A]"; sxoperand.help = "Resolution [A] to calculate the corresponding absolute frequency [1/Pixel] using a given pixel size [A/Pixel]."; sxoperand.register = " -9999.9"; sxoperand.type = "float"; sxoperand_set.list.append(sxoperand); sxoperand_set.dict[sxoperand.key] = sxoperand
		sxoperand = SXoperand(); sxoperand.key = "apix"; sxoperand.label = "Pixel size [A/Pixel]"; sxoperand.help = "Pixel size of dataset [A/Pixel]."; sxoperand.register = " -9999.9"; sxoperand.type = "apix"; sxoperand_set.list.append(sxoperand); sxoperand_set.dict[sxoperand.key] = sxoperand
		sxoperand = SXoperand(); sxoperand.key = "abs_freq"; sxoperand.label = "Absolute Frequency [1/Pixel]"; sxoperand.help = "Absolute frquency [1/Pixel] corresponding to the given resolution [A/Pixel]."; sxoperand.register = " -9999.9"; sxoperand.type = "abs_freq"; sxoperand.is_input = False; sxoperand_set.list.append(sxoperand); sxoperand_set.dict[sxoperand.key] = sxoperand

		# Store the project constant parameter set as a class data member
		return sxoperand_set
	
	def handle_regster_widget_event(self, sxoperand):
		sxoperand.widget.setText(sxoperand.register)
		self.handle_convert_units()

###	def handle_show_event(self):
###		# Use current value of resolution [A] widget as a register value of resolution [A] operand upon showing this dialog
###		sxoperand_ares = self.sxoperand_set.dict["ares"]
###		sxoperand_ares.register = self.sxcmd_token_subwidget_left_ares.text()
###		sxoperand_ares.register_widget.setText("%s" % sxoperand_ares.register)
###		sxoperand_ares.widget.setText(sxoperand_ares.register)
###		
###		# Show this dialog
###		self.show()

###	def handle_abs_freq_editing_finished_event(self):
###		apix_str = self.sxconst_register_widget_apix.text()
###		abs_freq_str = sxcmd_token_widget_abs_freq.text()
###		ares_str, is_valid_ares = selft.convert_abs_freq_to_ares(apix_str, abs_freq_str)
###		sxoperand_
###		print("MRK_DEBUG: ----- handle_abs_freq_editing_finished_event() ----- ")
###		print("MRK_DEBUG: Input Abs. Freq. [1/Pixel] string ; abs_freq_str := {}".format(abs_freq_str))
###		print("MRK_DEBUG: Input Pixel Size [A/Pixel] string ; apix_str := {}".format(apix_str))
###		print("MRK_DEBUG: Output Resolution [A] string      ; ares_str:= {}".format(ares_str))
###		if is_valid_ares:
###			sxcmd_token_subwidget_left_ares.setText("{}[A]@{}[A/Pix]".format(ares_str, apix_str))
###		else:
###			sxcmd_token_subwidget_left_ares.setText("Mode {}".format(ares_str))
	
###	def convert_abs_freq_to_ares():
###		# Initialise resolution [A] string to absolute frequency [1/Pixel] string
###		self.ares_str = self.abs_freq_str
###		self.is_valid_ares = False
###		
###		# Initialise resolution [A] and pixel size [A/Pixel] with invalid values
###		self.abs_freq_str abs_freq = 0.0
###		self.abs_freq_strapix = 0.0
###		# Check data type of arguments
###		try:
###			self.apix = float(self.apix_str)
###		except ValueError:
###			# The text string can not be converted to float. This must be special value of this option
###			assert (self.ares_str == self.abs_freq_str)
###			assert (not self.is_valid_ares)
###			return
###		try:
###			self.abs_freq = float(self.abs_freq_str)
###		except ValueError:
###			# The text string can not be converted to float. This must be special value of this option
###			assert (self.ares_str == self.abs_freq_str)
###			assert (not self.is_valid_ares)
###			return
###		# Check invalid values of arguments
###		if self.apix <= 0.0:
###			# The valid range of pixel size [A] is apix > 0.0. The valid value must have not been set yet.
###			assert (self.ares_str == self.abs_freq_str)
###			assert (not self.is_valid_ares)
###			return ares_str, self.is_valid_ares
###		assert (apix > 0.0)
###		if abs_freq <= 0.0 or abs_freq > 0.5:
###			# The valid range of absolut frequency is 0.0 < abs_freq <= 0.5. This must be special value of this option
###			assert (ares_str == abs_freq_str)
###			assert (not is_valid_ares)
###			return ares_str, is_valid_ares
###		assert (abs_freq > 0.0 or abs_freq <= 0.5)
###		
###		ares = round(apix/abs_freq, SXDialogCalculator.ares_precision)
###		ares_str = "{}".format(ares)
###		is_valid_ares = True
###		
###		print("MRK_DEBUG: ----- convert_abs_freq_to_ares() ----- ")
###		print("MRK_DEBUG: Input Abs. Freq. [1/Pixel]   ; abs_freq := {}".format(abs_freq))
###		print("MRK_DEBUG: Input Pixel Size [A/Pixel]   ; apix := {}".format(apix))
###		print("MRK_DEBUG: Output Resolution [A]        ; ares:= {}".format(ares))
###		print("MRK_DEBUG: Output Resolution [A] string ; ares_str:= {}".format(ares_str))
###		
###		return ares_str, is_valid_ares

	def handle_convert_units(self, is_enable_message = True):
		sxoperand_apix = self.sxoperand_set.dict["apix"]
		apix_str = sxoperand_apix.widget.text()
		sxoperand_ares = self.sxoperand_set.dict["ares"]
		ares_str = sxoperand_ares.widget.text()
		sxoperand_abs_freq = self.sxoperand_set.dict["abs_freq"]
		
		apix = self.convert_str_to_float_apix(apix_str)
		
		if apix is None:
			sxoperand_apix.validated = None
			sxoperand_ares.validated = None
			sxoperand_abs_freq.validated = None
			sxoperand_abs_freq.widget.setText("Invalid Pixel Size {}".format(apix_str))
			self.apply_btn.setEnabled(False)
			if is_enable_message:
				QMessageBox.warning(self, "Invalid Pixel Size [A/Pixel]", "Invalid Value {} for Pixel Size [A/Pixel] is provided. It must be a non-zero positive float value...".format(apix_str))
		else:
			assert (apix is not None)
			sxoperand_apix.validated = apix
			nyquist_res = 2.0 * apix
			nyquist_res_str = "{}".format(nyquist_res)
			ares = self.convert_str_to_float_ares(ares_str, nyquist_res)
			if ares is None:
				sxoperand_ares.validated = None
				sxoperand_abs_freq.validated = None
				sxoperand_abs_freq.widget.setText("Invalid Resolution {}".format(ares_str))
				self.apply_btn.setEnabled(False)
				if is_enable_message:
					QMessageBox.warning(self, "Invalid Resolution [A]", "Invalid Value {} for Resolution [A] is provided. It must be a float value, and larger than or equal to Nyquist resolution {} [A] (2.0 * Pixel Size [A/Pixel])...".format(ares_str, nyquist_res_str))
			else:
				assert (ares >= nyquist_res)
				sxoperand_ares.validated = ares
				abs_freq = round(apix/ares, self.precision_abs_freq)
				# The valid range of absolute frequency [1/Pixel] is 0.0 < abs_freq <= 0.5
				# If both pixel size and resolution values are valid, the absolute frequency must be always valid.
				assert (abs_freq > 0.0 or abs_freq <= 0.5)
				sxoperand_abs_freq.validated = abs_freq
				abs_freq_str = "{}".format(abs_freq)
				# Update widget associating to absolute frequency [1/Pixel] in this dialog
				sxoperand_abs_freq.widget.setText(abs_freq_str)
				self.apply_btn.setEnabled(True)
		
		
#		if sxoperand_apix.validated is None:
#			sxoperand_abs_freq.validated = None
#			sxoperand_abs_freq.widget.setText("Invalid Pixel Size {}".format(apix_str))
#			self.apply_btn.setEnabled(False)
#			if is_enable_message:
#				QMessageBox.warning(self, "Invalid Pixel Size [A/Pixel]", "Invalid Value {} for Pixel Size [A/Pixel] is provided. It must be a non-zero positive float value...".format(apix_str))
#		else:
#			assert (sxoperand_apix.validated is not None)
#			assert (apix > 0.0)
#			nyquist_res = 2.0 * apix
#			nyquist_res_str = "{}".format(nyquist_res)
#			if sxoperand_ares.validated is None:
#				sxoperand_abs_freq.validated = None
#				sxoperand_abs_freq.widget.setText("Invalid Resolution {}".format(ares_str))
#				self.apply_btn.setEnabled(False)
#				if is_enable_message:
#					QMessageBox.warning(self, "Invalid Resolution [A]", "Invalid Value {} for Resolution [A] is provided. It must be a float value, and larger than or equal to Nyquist resolution {} [A] (2.0 * Pixel Size [A/Pixel])...".format(ares_str, nyquist_res_str))
#			else:
#				assert (sxoperand_apix.validated is not None)
#				assert (nyquist_res * 2.0 )
#				abs_freq = round(apix/ares, self.precision_abs_freq)
#				# The valid range of absolute frequency [1/Pixel] is 0.0 < abs_freq <= 0.5
#				# If both pixel size and resolution values are valid, the absolute frequency must be always valid.
#				assert (abs_freq > 0.0 or abs_freq <= 0.5)
#				sxoperand_abs_freq.validated = abs_freq
#				abs_freq_str = "{}".format(abs_freq)
#				# Update widget associating to absolute frequency [1/Pixel] in this dialog
#				sxoperand_abs_freq.widget.setText(abs_freq_str)
#				self.apply_btn.setEnabled(True)
		
		# print("MRK_DEBUG: ")
		# print("MRK_DEBUG: ----- SXDialogCalculator.handle_convert_units() ----- ")
		# print("MRK_DEBUG: Edit Widget Resolution [A]              ; sxoperand_ares.widget.text()              := \"{}\"".format(sxoperand_ares.widget.text()))
		# print("MRK_DEBUG: Edit Validated Resolution [A]           ; sxoperand_ares.validated                  := \"{}\"".format(sxoperand_ares.validated))
		# print("MRK_DEBUG: Register Widget Resolution [A]          ; sxoperand_ares.register_widget.text()     := \"{}\"".format(sxoperand_ares.register_widget.text()))
		# print("MRK_DEBUG: Register Validated Resolution [A]       ; sxoperand_ares.validated_register         := \"{}\"".format(sxoperand_ares.validated_register))
		# print("MRK_DEBUG: Register Resolution [A]                 ; sxoperand_ares.register                   := \"{}\"".format(sxoperand_ares.register))
		# 
		# print("MRK_DEBUG: Edit Widget Pixel Size [A/Pixel]        ; sxoperand_apix.widget.text()              := \"{}\"".format(sxoperand_apix.widget.text()))
		# print("MRK_DEBUG: Edit Validated Abs. Freq. [1/Pixel]     ; sxoperand_apix.validated                  := \"{}\"".format(sxoperand_apix.validated))
		# print("MRK_DEBUG: Register Widget Pixel Size [A/Pixel]    ; sxoperand_apix.register_widget.text()     := \"{}\"".format(sxoperand_apix.register_widget.text()))
		# print("MRK_DEBUG: Register Validated Pixel Size [A/Pixel] ; sxoperand_apix.validated_register         := \"{}\"".format(sxoperand_apix.validated_register))
		# print("MRK_DEBUG: Register Pixel Size [A/Pixel]           ; sxoperand_apix.register                   := \"{}\"".format(sxoperand_apix.register))
		# 
		# print("MRK_DEBUG: Edit Widget bs. Freq. [1/Pixel]         ; sxoperand_abs_freq.widget.text()          := \"{}\"".format(sxoperand_abs_freq.widget.text()))
		# print("MRK_DEBUG: Edit Validated Abs. Freq. [1/Pixel]     ; sxoperand_abs_freq.validated              := \"{}\"".format(sxoperand_abs_freq.validated))
		# print("MRK_DEBUG: Register Widget Abs. Freq. [1/Pixel]    ; sxoperand_abs_freq.register_widget.text() := \"{}\"".format(sxoperand_abs_freq.register_widget.text()))
		# print("MRK_DEBUG: Register Validated Abs. Freq. [1/Pixel] ; sxoperand_abs_freq.validated_register     := \"{}\"".format(sxoperand_abs_freq.validated_register))
		# print("MRK_DEBUG: Register Abs. Freq. [1/Pixel]           ; sxoperand_abs_freq.register               := \"{}\"".format(sxoperand_abs_freq.register))

	def handle_apply_unit_conversion(self):
		# Do the unit conversion first
		self.handle_convert_units()
		sxoperand_abs_freq = self.sxoperand_set.dict["abs_freq"]
		sxoperand_ares = self.sxoperand_set.dict["ares"]
		sxoperand_apix = self.sxoperand_set.dict["apix"]
		
		abs_freq = sxoperand_abs_freq.validated
		ares = sxoperand_ares.validated
		apix = sxoperand_apix.validated
		# This button should have been disabled if all operands are valid
		if abs_freq is None:
			ERROR("Logical Error: Encountered unexpected condition in SXDialogCalculator.handle_apply_unit_conversion(). Consult with the developer.", "%s in %s" % (__name__, os.path.basename(__file__)))
			self.sxcmd_token_widget_abs_freq.setText("Logical Error")
			self.sxcmd_token_subwidget_left_ares.setText("Logical Error")
		else:
			assert (abs_freq is not None and ares is not None and apix is not None)
			# Update register, validated_register, and register_widget of absolute frequency [1/Pixel]
			sxoperand_abs_freq.validated_register = abs_freq
			sxoperand_abs_freq.register = "{}".format(sxoperand_abs_freq.validated_register)
			sxoperand_abs_freq.register_widget.setText(sxoperand_abs_freq.register)

			# Update register, validated_register, and register_widget of resolution [A]
			sxoperand_ares.validated_register = ares
			sxoperand_ares.register = "{}".format(sxoperand_ares.validated_register)
			sxoperand_ares.register_widget.setText(sxoperand_ares.register)

			# DO NOT update register, validated_register, and register_widget of pixel size [A/Pixel]
			# These should be updated through only project constants settings
			
			# Update command token subwidget associating to resolution [A] accordingly
			self.sxcmd_token_widget_abs_freq.setText(sxoperand_abs_freq.register)
			self.sxcmd_token_subwidget_left_ares.setText("{}[A]@{}[A/Pix]".format(ares, apix))
			
			# Update the associated pixel size token of this command if there is any
			if self.sxcmd_token_widget_apix is not None:
				self.sxcmd_token_widget_apix.setText("{}".format(apix))
			for sxcmd_token_other_dialog_abs_freq in self.sxcmd_token_other_dialog_list_abs_freq:
				assert (sxcmd_token_other_dialog_abs_freq is not self)
				sxcmd_token_other_dialog_abs_freq.reflect_external_local_update_apix()
		
		# print("MRK_DEBUG: ")
		# print("MRK_DEBUG: ----- SXDialogCalculator.handle_apply_unit_conversion() ----- ")
		# print("MRK_DEBUG: Edit Widget Resolution [A]              ; sxoperand_ares.widget.text()              := \"{}\"".format(sxoperand_ares.widget.text()))
		# print("MRK_DEBUG: Edit Validated Resolution [A]           ; sxoperand_ares.validated                  := \"{}\"".format(sxoperand_ares.validated))
		# print("MRK_DEBUG: Register Widget Resolution [A]          ; sxoperand_ares.register_widget.text()     := \"{}\"".format(sxoperand_ares.register_widget.text()))
		# print("MRK_DEBUG: Register Validated Resolution [A]       ; sxoperand_ares.validated_register         := \"{}\"".format(sxoperand_ares.validated_register))
		# print("MRK_DEBUG: Register Resolution [A]                 ; sxoperand_ares.register                   := \"{}\"".format(sxoperand_ares.register))
		# 
		# print("MRK_DEBUG: Edit Widget Pixel Size [A/Pixel]        ; sxoperand_apix.widget.text()              := \"{}\"".format(sxoperand_apix.widget.text()))
		# print("MRK_DEBUG: Edit Validated Abs. Freq. [1/Pixel]     ; sxoperand_apix.validated                  := \"{}\"".format(sxoperand_apix.validated))
		# print("MRK_DEBUG: Register Widget Pixel Size [A/Pixel]    ; sxoperand_apix.register_widget.text()     := \"{}\"".format(sxoperand_apix.register_widget.text()))
		# print("MRK_DEBUG: Register Validated Pixel Size [A/Pixel] ; sxoperand_apix.validated_register         := \"{}\"".format(sxoperand_apix.validated_register))
		# print("MRK_DEBUG: Register Pixel Size [A/Pixel]           ; sxoperand_apix.register                   := \"{}\"".format(sxoperand_apix.register))
		# 
		# print("MRK_DEBUG: Edit Widget bs. Freq. [1/Pixel]         ; sxoperand_abs_freq.widget.text()          := \"{}\"".format(sxoperand_abs_freq.widget.text()))
		# print("MRK_DEBUG: Edit Validated Abs. Freq. [1/Pixel]     ; sxoperand_abs_freq.validated              := \"{}\"".format(sxoperand_abs_freq.validated))
		# print("MRK_DEBUG: Register Widget Abs. Freq. [1/Pixel]    ; sxoperand_abs_freq.register_widget.text() := \"{}\"".format(sxoperand_abs_freq.register_widget.text()))
		# print("MRK_DEBUG: Register Validated Abs. Freq. [1/Pixel] ; sxoperand_abs_freq.validated_register     := \"{}\"".format(sxoperand_abs_freq.validated_register))
		# print("MRK_DEBUG: Register Abs. Freq. [1/Pixel]           ; sxoperand_abs_freq.register               := \"{}\"".format(sxoperand_abs_freq.register))
		
		self.close()

	def convert_str_to_float_apix(self, apix_str):
		# Initialise with easy-to-recognise invalid value, which won't be used.
		apix = -9999.9 
		# Check data type of string
		try:
			apix = float(apix_str)
		except ValueError:
			# The text string can not be converted to float.
			return None
		# round by predefined precision so that calculated resolution [A] is likely to be a simple float value (e.g. 3.0[A])
		apix = round(apix, self.precision_apix)
		# Check value validity. The valid range of pixel size [A] is apix > 0.0.
		if apix <= 0.0:
			return None
		assert (apix > 0.0)
		return apix

	def register_operand_apix_str(self, apix_str):
		# Validated the operand value
		sxoperand_apix = self.sxoperand_set.dict["apix"]
		sxoperand_apix.register = apix_str
		apix = self.convert_str_to_float_apix(sxoperand_apix.register)
		if apix is None:
			sxoperand_apix.validated_register = None
		else:
			sxoperand_apix.validated_register = apix
			assert (sxoperand_apix.validated_register is not None)
		sxoperand_apix.validated = sxoperand_apix.validated_register
		# Set the string regardless of the validity
		# Update widgets associating to pixel size [A/Pixel] in this dialog
		sxoperand_apix.register_widget.setText(sxoperand_apix.register)
		sxoperand_apix.widget.setText(sxoperand_apix.register)

	def edit_operand_apix_str(self, apix_str):
		# Validated the operand value
		sxoperand_apix = self.sxoperand_set.dict["apix"]
		apix = self.convert_str_to_float_apix(apix_str)
		if apix is None:
			sxoperand_apix.validated = None
		else:
			sxoperand_apix.validated = apix
			assert (sxoperand_apix.validated is not None)
		# Set the string regardless of the validity
		# Update widgets associating to pixel size [A/Pixel] in this dialog
		sxoperand_apix.widget.setText(apix_str)

	def convert_str_to_float_abs_freq(self, abs_freq_str):
		# Initialise with easy-to-recognise invalid value, which won't be used.
		abs_freq = -9999.9 
		# Check data type of arguments
		try:
			abs_freq = float(abs_freq_str)
		except ValueError:
			# The text string can not be converted to float.
			return None
		# round by predefined precision so that calculated resolution [A] is likely to be a simple float value (e.g. 3.0[A])
		abs_freq = round(abs_freq, self.precision_abs_freq)
		# Check value validity. The valid range of absolut frequency is 0.0 < abs_freq <= 0.5. 
		if abs_freq <= 0.0 or abs_freq > 0.5:
			return None
		assert (abs_freq > 0.0 or abs_freq <= 0.5)
		return abs_freq

	def register_operand_abs_freq_str(self, abs_freq_str):
		# Validated the operand value
		sxoperand_abs_freq = self.sxoperand_set.dict["abs_freq"]
		sxoperand_abs_freq.register = abs_freq_str
		abs_freq = self.convert_str_to_float_abs_freq(sxoperand_abs_freq.register)
		if abs_freq is None:
			sxoperand_abs_freq.validated_register = None
		else:
			sxoperand_abs_freq.validated_register = abs_freq
			assert (sxoperand_abs_freq.validated_register is not None)
		sxoperand_abs_freq.validated = sxoperand_abs_freq.validated_register
		# Set the string regardless of the validity
		# Update widgets associating to absolute frequency [1/Pixel] in this dialog
		sxoperand_abs_freq.register_widget.setText(sxoperand_abs_freq.register)
		sxoperand_abs_freq.widget.setText(sxoperand_abs_freq.register)

	def convert_str_to_float_ares(self, ares_str, nyquist_res = 0.0):
		# Initialise with easy-to-recognise invalid value, which won't be used.
		ares = -9999.9 
		# Check data type of arguments
		try:
			ares = float(ares_str)
		except ValueError:
			# The text string can not be converted to float. This must be special value of this option
			return None
		# Check value validity. The valid range of resolution is 0.0 < ares
		# Here, we ignore Nyquist with pixel size.
		ares = round(ares, self.precision_ares)
		if ares < nyquist_res:
			return None
		assert (ares >= nyquist_res)
		return ares

	# def rigister_operand_ares_str(self, ares_str, nyquist_res = 0.0):
	# 	sxoperand_ares = self.sxoperand_set.dict["ares"]
	# 	sxoperand_ares.register = ares_str
	# 	ares = self.convert_str_to_float_ares(sxoperand_ares.register, nyquist_res)
	# 	if ares is None:
	# 		sxoperand_ares.validated_register = None
	# 		return
	# 	sxoperand_ares.validated_register = ares
	# 	assert (sxoperand_ares.validated_register is not None)
	# 	return
	
	def synchronize_external_update_to_ares(self):
		# print("MRK_DEBUG: ----- synchronize_external_update_to_ares ----- ")
		sxoperand_apix = self.sxoperand_set.dict["apix"]
		sxoperand_abs_freq = self.sxoperand_set.dict["abs_freq"]
		sxoperand_ares = self.sxoperand_set.dict["ares"]
		if sxoperand_abs_freq.validated is None:
			# Use registered string of absolute frequency. This can be special value of this option (e.g. indicate lpf modes)
			sxoperand_abs_freq.widget.setText("Mode {}".format(sxoperand_abs_freq.register))
			# sxoperand_ares.register = sxoperand_abs_freq.widget.text()
			sxoperand_ares.register = "No Default"
			sxoperand_ares.validated_register = None
		else:
			assert (sxoperand_abs_freq.validated is not None)
			if sxoperand_apix.validated is None:
				# Invalid value to indicate the invalid pixel size but invalid absolute frequency value...
				sxoperand_ares.validated_register = None
				sxoperand_ares.register = "Invalid Pixel Size {}".format(sxoperand_apix.widget.text())
				sxoperand_abs_freq.widget.setText(sxoperand_ares.register)
			else:
				assert (sxoperand_apix.validated is not None)
				apix = sxoperand_apix.validated
				abs_freq = sxoperand_abs_freq.validated
				ares = round(apix/abs_freq, self.precision_ares)
				# The valid range of resolution [A] is ares <= apix * 2.0 (Nyquist resolution)
				# If both pixel size and absolute frequency values are valid, this must be always true.
				assert (ares >= apix * 2.0 ) 
				sxoperand_ares.validated_register = ares
				sxoperand_ares.register = "{}".format(sxoperand_ares.validated_register)
		sxoperand_ares.validated = sxoperand_ares.validated_register
		
		# Update widgets associating to resolution [A] in this dialog
		# and command token subwidget associating to resolution [A]
		if sxoperand_ares.validated_register is not None:
			self.apply_btn.setEnabled(True)
			sxoperand_ares.register_widget.setEnabled(True)
			sxoperand_ares.register_widget.setText(sxoperand_ares.register)
			sxoperand_ares.widget.setText(sxoperand_ares.register)
			# Update widget of this command token for absolute frequency [1/Pixel] associating to this dialog
			self.sxcmd_token_widget_abs_freq.setText("{}".format(sxoperand_abs_freq.register))
			# Update left subwidget of this command token for absolute frequency [1/Pixel] associating to this dialog
			self.sxcmd_token_subwidget_left_ares.setText("{}[A]@{}[A/Pix]".format(sxoperand_ares.validated, sxoperand_apix.validated))
		else:
			assert (sxoperand_ares.validated_register is None)
			self.apply_btn.setEnabled(False)
			sxoperand_ares.register_widget.setEnabled(False)
			# sxoperand_ares.register_widget.setText("No Default")
			sxoperand_ares.register_widget.setText(sxoperand_ares.register)
			sxoperand_ares.widget.setText("")
			if sxoperand_abs_freq.validated is None:
				self.sxcmd_token_subwidget_left_ares.setText("Mode {}".format(sxoperand_abs_freq.register))
			else:
				assert (sxoperand_abs_freq.validated is not None)
				self.sxcmd_token_subwidget_left_ares.setText("{}".format(sxoperand_ares.register))
		
		# print("MRK_DEBUG: ")
		# print("MRK_DEBUG: ----- SXDialogCalculator.synchronize_external_update_to_ares() ----- ")
		# print("MRK_DEBUG: Edit Widget Resolution [A]              ; sxoperand_ares.widget.text()              := \"{}\"".format(sxoperand_ares.widget.text()))
		# print("MRK_DEBUG: Edit Validated Resolution [A]           ; sxoperand_ares.validated                  := \"{}\"".format(sxoperand_ares.validated))
		# print("MRK_DEBUG: Register Widget Resolution [A]          ; sxoperand_ares.register_widget.text()     := \"{}\"".format(sxoperand_ares.register_widget.text()))
		# print("MRK_DEBUG: Register Validated Resolution [A]       ; sxoperand_ares.validated_register         := \"{}\"".format(sxoperand_ares.validated_register))
		# print("MRK_DEBUG: Register Resolution [A]                 ; sxoperand_ares.register                   := \"{}\"".format(sxoperand_ares.register))
		# 
		# print("MRK_DEBUG: Edit Widget Pixel Size [A/Pixel]        ; sxoperand_apix.widget.text()              := \"{}\"".format(sxoperand_apix.widget.text()))
		# print("MRK_DEBUG: Edit Validated Abs. Freq. [1/Pixel]     ; sxoperand_apix.validated                  := \"{}\"".format(sxoperand_apix.validated))
		# print("MRK_DEBUG: Register Widget Pixel Size [A/Pixel]    ; sxoperand_apix.register_widget.text()     := \"{}\"".format(sxoperand_apix.register_widget.text()))
		# print("MRK_DEBUG: Register Validated Pixel Size [A/Pixel] ; sxoperand_apix.validated_register         := \"{}\"".format(sxoperand_apix.validated_register))
		# print("MRK_DEBUG: Register Pixel Size [A/Pixel]           ; sxoperand_apix.register                   := \"{}\"".format(sxoperand_apix.register))
		# 
		# print("MRK_DEBUG: Edit Widget bs. Freq. [1/Pixel]         ; sxoperand_abs_freq.widget.text()          := \"{}\"".format(sxoperand_abs_freq.widget.text()))
		# print("MRK_DEBUG: Edit Validated Abs. Freq. [1/Pixel]     ; sxoperand_abs_freq.validated              := \"{}\"".format(sxoperand_abs_freq.validated))
		# print("MRK_DEBUG: Register Widget Abs. Freq. [1/Pixel]    ; sxoperand_abs_freq.register_widget.text() := \"{}\"".format(sxoperand_abs_freq.register_widget.text()))
		# print("MRK_DEBUG: Register Validated Abs. Freq. [1/Pixel] ; sxoperand_abs_freq.validated_register     := \"{}\"".format(sxoperand_abs_freq.validated_register))
		# print("MRK_DEBUG: Register Abs. Freq. [1/Pixel]           ; sxoperand_abs_freq.register               := \"{}\"".format(sxoperand_abs_freq.register))
		
	def reflect_external_global_update_apix(self):
		# print("MRK_DEBUG: ")
		# print("MRK_DEBUG: ----- SXDialogCalculator.reflect_external_global_update_apix() ----- ")
		# print("MRK_DEBUG: self.sxconst_register_widget_apix.text() := \"{}\"".format(self.sxconst_register_widget_apix.text()))
		apix_str = self.sxconst_register_widget_apix.text()
		self.register_operand_apix_str(apix_str)
		self.synchronize_external_update_to_ares()

	def reflect_external_local_update_apix(self):
		# print("MRK_DEBUG: ")
		# print("MRK_DEBUG: ----- SXDialogCalculator.reflect_external_local_update_apix() ----- ")
		# print("MRK_DEBUG: self.sxcmd_token_widget_apix.text() := \"{}\"".format(self.sxcmd_token_widget_apix.text()))
		apix_str = self.sxcmd_token_widget_apix.text()
		self.edit_operand_apix_str(apix_str)
		self.synchronize_external_update_to_ares()

	def reflect_external_local_update_abs_freq(self):
		# print("MRK_DEBUG: ")
		# print("MRK_DEBUG: ----- SXDialogCalculator.reflect_external_local_update_abs_freq() ----- ")
		# print("MRK_DEBUG: self.sxcmd_token_widget_abs_freq.text() := \"{}\"".format(self.sxcmd_token_widget_abs_freq.text()))
		abs_freq_str = self.sxcmd_token_widget_abs_freq.text()
		self.register_operand_abs_freq_str(abs_freq_str)
		self.synchronize_external_update_to_ares()

	def reflect_external_local_update_apix_and_abs_freq(self):
		# print("MRK_DEBUG: ")
		# print("MRK_DEBUG: ----- SXDialogCalculator.reflect_external_local_update_apix_and_abs_freq() ----- ")
		# print("MRK_DEBUG: self.sxcmd_token_widget_apix.text() := \"{}\"".format(self.sxcmd_token_widget_apix.text()))
		# print("MRK_DEBUG: self.sxcmd_token_widget_abs_freq.text() := \"{}\"".format(self.sxcmd_token_widget_abs_freq.text()))
		apix_str = self.sxcmd_token_widget_apix.text()
		self.edit_operand_apix_str(apix_str)
		abs_freq_str = self.sxcmd_token_widget_abs_freq.text()
		self.register_operand_abs_freq_str(abs_freq_str)
		self.synchronize_external_update_to_ares()

###	@staticmethod
###	def convert_ares_to_abs_freq(apix_str, ares_str):
###		# Initialise resolution [A] string to absolute frequency [1/Pixel] string
###		abs_freq_str = ares_str
###		is_valid_abs_freq = False
###		
###		# Initialise resolution [A] and pixel size [A/Pixel] with invalid values
###		ares = 0.0
###		apix = 0.0
###		# Check data type of arguments
###		try:
###			apix = float(apix_str)
###		except ValueError:
###			# The text string can not be converted to float. This must be special value of this option
###			assert (abs_freq_str == ares_str)
###			assert (not is_valid_abs_freq)
###			return abs_freq_str, is_valid_abs_freq
###		try:
###			ares = float(ares_str)
###		except ValueError:
###			# The text string can not be converted to float. This must be special value of this option
###			assert (abs_freq_str == ares_str)
###			assert (not is_valid_abs_freq)
###			return abs_freq_str, is_valid_abs_freq
###		if apix <= 0.0:
###			# The valid range of pixel size [A] is apix > 0.0. The valid value must have not been set yet.
###			assert (abs_freq_str == ares_str)
###			assert (not is_valid_abs_freq)
###			return abs_freq_str, is_valid_abs_freq
###		assert (apix > 0.0)
###		nyquist_res = apix * 2.0
###		# Check invalid values of arguments
###		if ares < nyquist_res:
###			# The valid range of absolut frequency is 0.0 < abs_freq <= 0.5. This must be special value of this option
###			assert (abs_freq_str == ares_str)
###			assert (not is_valid_abs_freq)
###			return abs_freq_str, is_valid_abs_freq
###		assert (ares >= nyquist_res)
###		
###		abs_freq = round(apix/ares, SXDialogCalculator.abs_freq_precision)
###		abs_freq_str = "{}".format(abs_freq)
###		is_valid_abs_freq = True
###		
###		print("MRK_DEBUG: ----- convert_ares_to_abs_freq() ----- ")
###		print("MRK_DEBUG: Input Abs. Freq. [1/Pixel]   ; abs_freq := {}".format(abs_freq))
###		print("MRK_DEBUG: Input Pixel Size [A/Pixel]   ; apix := {}".format(apix))
###		print("MRK_DEBUG: Output Resolution [A]        ; ares:= {}".format(ares))
###		print("MRK_DEBUG: Output Resolution [A] string ; ares_str:= {}".format(ares_str))
###		
###		return abs_freq_str, is_valid_abs_freq

###	def convert_units(self):
###		# Initialise resolution [A] and pixel size [A/Pixel] with invalid values
###		apix = 0.0
###		ares = 0.0
###		# Get resolution [A] and pixel size [A/Pixel] values from the edit boxes
###		apix_str = self.sxoperand_set.dict["apix"].widget.text()
###		ares_str = self.sxoperand_set.dict["ares"].widget.text()
###		try:
###			apix = float(apix_str)
###		except ValueError:
###			# Check invalid values make sure the text string can be converted to float
###			QMessageBox.warning(self, "Invalid Value for Pixel Size [A/Pixel]", "Invalid Value {} for Pixel Size [A/Pixel] is provided. It must be a non-zero positive float value...".format(apix_str))
###			return
###		try:
###			ares = float(ares_str)
###		except ValueError:
###			# Check invalid values make sure the text string can be converted to float
###			QMessageBox.warning(self, "Invalid Value for Resolution [A]", "Invalid Value {} for Resolution [A] is provided. It must be a non-zero positive float value...".format(ares_str))
###			return
###		# Check invalid values
###		if apix <= 0.0:
###			QMessageBox.warning(self, "Invalid Value for Pixel Size [A/Pixel]", "Invalid Value {} for Pixel Size [A/Pixel] is provided. It must be a non-zero positive float value...".format(apix))
###			return
###		assert (apix > 0.0)
###		nyquist_res = apix * 2.0
###		if ares < nyquist_res:
###			QMessageBox.warning(self, "Invalid Value for Resolution [A]", "Invalid Value {} for Resolution [A] is provided. It must be larger than or equal to the Nyquest resolution {} [A] of the provided pixel size {} [A/Pixel]...".format(ares, nyquist_res, apix))
###			return
###		assert (ares >= nyquist_res)
###		
###		self.abs_freq = round(apix/ares, self.abs_freq_precision)
###		self.apply_btn.setText("Apply Abs. Freq. {} [1/Pixel]".format(self.abs_freq))
###		
###		print("MRK_DEBUG: ----- convert_units() ----- ")
###		print("MRK_DEBUG: Input Resolution [A]       ; ares := {}".format(ares))
###		print("MRK_DEBUG: Input Pixel Size [A/Pixel] ; apix := {}".format(apix))
###		print("MRK_DEBUG: Output Abs. Freq. [1/Pixel]; self.abs_freq:= {}".format(self.abs_freq))
###		print("MRK_DEBUG: Output Resolution [A]      ; self.abs_freq:= {}".format(round(apix/self.abs_freq, self.ares_precision)))

# ========================================================================================
# Main Window (started by class SXApplication)
class SXMainWindow(QMainWindow): # class SXMainWindow(QWidget):

	def __init__(self, helical, parent = None):
		super(SXMainWindow, self).__init__(parent)

		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		self.helical = helical
		self.sxinfo = []
		self.sxconst_set = []
		self.sxcmd_category_list = []

		self.cur_sxmenu_item = None
		self.sxmenu_item_widget_stacked_layout = []

		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

		# --------------------------------------------------------------------------------
		# Construct menu items
		# --------------------------------------------------------------------------------
		self.sxinfo.append(self.construct_sxinfo())              # Construct application information
		self.sxinfo.append(self.construct_sxinfo())              # Construct application information
		self.sxconst_set.append(self.construct_sxconst_set(0))         # Construct project constant set for project settings
		self.sxconst_set.append(self.construct_sxconst_set(1))         # Construct project constant set for project settings
		self.sxcmd_category_list.append(self.construct_sxcmd_category_list()) # Construct list of categorised sxscript objects (extracted from associated wiki documents)
		self.sxcmd_category_list.append(self.construct_sxcmd_category_list()) # Construct list of categorised sxscript objects (extracted from associated wiki documents)

		# --------------------------------------------------------------------------------
		# Setup Window Layout
		# --------------------------------------------------------------------------------
		background_image_file_path = '{0}sxgui_background.png'.format(get_image_directory())

		central_widget_global = QWidget(self)
		central_widget_global.setObjectName('central')
		central_widget_global.setStyleSheet(
			'QWidget#central {{background-image: url("{0}")}}'.format(background_image_file_path)
			)
		self.setCentralWidget(central_widget_global)
		self.sxmenu_item_widget_stacked_layout_global = QStackedLayout()
		central_widget_global.setLayout(self.sxmenu_item_widget_stacked_layout_global)

		for idx in range(2):
			# Central widget
			central_widget = QWidget(self)
			central_widget.setObjectName('central')
			central_widget.setStyleSheet(
				'QWidget#central {{background-image: url("{0}")}}'.format(background_image_file_path)
				)

			# Layout for central widget
			central_layout = QHBoxLayout(central_widget)
			central_widget.setLayout(central_layout)

			# --------------------------------------------------------------------------------
			# Construct and add a widget for menu item button area (containing all menu item buttons)
			# --------------------------------------------------------------------------------
			sxmenu_item_btn_area_widget = SXMenuItemBtnAreaWidget(self.sxconst_set[idx], self.sxcmd_category_list[idx], self.sxinfo[idx], central_widget)
			central_layout.addWidget(sxmenu_item_btn_area_widget)

			# --------------------------------------------------------------------------------
			# Construct and add widgets for menu item widget area (containing all menu item widgets)
			# --------------------------------------------------------------------------------
			# Stacked layout for sx menu item widgets area
			self.sxmenu_item_widget_stacked_layout.append(QStackedLayout())
			central_layout.addLayout(self.sxmenu_item_widget_stacked_layout[idx], stretch = 1)

			# Construct and add a widget for project constants settings
			self.sxconst_set[idx].widget = SXConstSetWidget(self.sxconst_set[idx], self.sxcmd_category_list[idx], helical=idx, parent=self)
			self.sxmenu_item_widget_stacked_layout[idx].addWidget(self.sxconst_set[idx].widget)

			# Construct and add widgets for sx command categories
			for sxcmd_category in self.sxcmd_category_list[idx]:
				# Create SXCmdCategoryWidget for this command category
				sxcmd_category.widget = SXCmdCategoryWidget(self.sxconst_set[idx], sxcmd_category, idx)
				self.sxmenu_item_widget_stacked_layout[idx].addWidget(sxcmd_category.widget)

			# Construct and add a widget for GUI application information
			self.sxinfo[idx].widget = SXInfoWidget()
			self.sxmenu_item_widget_stacked_layout[idx].addWidget(self.sxinfo[idx].widget)

			# --------------------------------------------------------------------------------
			# Set up event handler of all menu item buttons
			# --------------------------------------------------------------------------------
			for i, sxcmd_category in enumerate(self.sxcmd_category_list[idx]):
				sxcmd_category.btn.clicked.connect(partial(self.handle_sxmenu_item_btn_event, sxcmd_category, i))
			self.sxconst_set[idx].btn.clicked.connect(partial(self.handle_sxmenu_item_btn_event, self.sxconst_set[idx]))
			self.sxinfo[idx].btn.clicked.connect(partial(self.handle_sxmenu_item_btn_event, self.sxinfo[idx]))

			# --------------------------------------------------------------------------------
			# Register project constant parameter settings upon initialization
			# --------------------------------------------------------------------------------
			self.sxconst_set[idx].widget.register_const_set()

			# --------------------------------------------------------------------------------
			# Load the previously saved parameter setting of all sx commands
			# Override the registration of project constant parameter settings with the previously-saved one
			# --------------------------------------------------------------------------------
			for sxcmd_category in self.sxcmd_category_list[idx]:
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
			self.sxmenu_item_widget_stacked_layout[idx].addWidget(start_widget)

			# --------------------------------------------------------------------------------
			# Display application information upon startup
			# --------------------------------------------------------------------------------
			self.sxmenu_item_widget_stacked_layout[idx].setCurrentWidget(start_widget)
			self.sxmenu_item_widget_stacked_layout_global.addWidget(central_widget)
		self.sxmenu_item_widget_stacked_layout_global.setCurrentIndex(self.helical)

		# --------------------------------------------------------------------------------
		# Get focus to main window
		# --------------------------------------------------------------------------------
		self.setFocus()

	def construct_sxinfo(self):
###		sxinfo = SXmenu_item(); sxinfo.name = "GUI Information"; sxinfo.label = "GUI Appliation Information"; sxinfo.short_info = "DUMMY STRING"
		sxinfo = SXmenu_item(); sxinfo.name = "sxc_gui_info"; sxinfo.label = "GUI Appliation Information"; sxinfo.short_info = "DUMMY STRING"

		# Store GUI application information as a class data member
		return sxinfo

	def construct_sxconst_set(self, helical):
		sxconst_set = SXconst_set(); sxconst_set.name = "sxc_project"; sxconst_set.label = "Project Settings"; sxconst_set.short_info = "Set constant parameter values for this project. These constants will be used as default values of associated arguments and options in command settings. However, the project settings here are not required to run commands."
		sxconst = SXconst(); sxconst.key = "protein"; sxconst.label = "Protein name"; sxconst.help = "a valid string for file names on your OS."; sxconst.register = "MY_PROTEIN"; sxconst.type = "string"; sxconst_set.list.append(sxconst); sxconst_set.dict[sxconst.key] = sxconst
		sxconst = SXconst(); sxconst.key = "apix"; sxconst.label = "Micrograph pixel size [A]"; sxconst.help = ""; sxconst.register = "1.0"; sxconst.type = "float"; sxconst_set.list.append(sxconst); sxconst_set.dict[sxconst.key] = sxconst
		sxconst = SXconst(); sxconst.key = "ctfwin"; sxconst.label = "CTF window size [pixels]"; sxconst.help = "it should be slightly larger than particle box size"; sxconst.register = "512"; sxconst.type = "int"; sxconst_set.list.append(sxconst); sxconst_set.dict[sxconst.key] = sxconst
		sxconst = SXconst(); sxconst.key = "box"; sxconst.label = "Particle box size [pixels]" ; sxconst.help = ""; sxconst.register = "-1"; sxconst.type = "int"; sxconst_set.list.append(sxconst); sxconst_set.dict[sxconst.key] = sxconst
		sxconst = SXconst(); sxconst.key = "radius"; sxconst.label = "Protein particle radius [pixels]"; sxconst.help = ""; sxconst.register = "-1"; sxconst.type = "int"; sxconst_set.list.append(sxconst); sxconst_set.dict[sxconst.key] = sxconst
		sxconst = SXconst(); sxconst.key = "sym"; sxconst.label = "Point-group symmetry"; sxconst.help = "e.g. c1, c4, d5"; sxconst.register = "c1"; sxconst.type = "string"; sxconst_set.list.append(sxconst); sxconst_set.dict[sxconst.key] = sxconst
		sxconst = SXconst(); sxconst.key = "mass"; sxconst.label = "Protein molecular mass [kDa]"; sxconst.help = ""; sxconst.register = "-1.0"; sxconst.type = "float"; sxconst_set.list.append(sxconst); sxconst_set.dict[sxconst.key] = sxconst
		if helical:
			sxconst = SXconst(); sxconst.key = "filament_width"; sxconst.label = "Filament width [pixels]"; sxconst.help = "Width of the filament in case of filamentous processing."; sxconst.register = "-1.0"; sxconst.type = "float"; sxconst_set.list.append(sxconst); sxconst_set.dict[sxconst.key] = sxconst
		sxconst = SXconst(); sxconst.key = "config"; sxconst.label = "Imaging configurations"; sxconst.help = "a free-style string for your record. please use it to describe the set of imaging configurations used in this project (e.g. types of microscope, detector, enegy filter, abbration corrector, phase plate, and etc."; sxconst.register = "MY_MICROSCOPE"; sxconst.type = "int"; sxconst_set.list.append(sxconst); sxconst_set.dict[sxconst.key] = sxconst

		# Store the project constant parameter set as a class data member
		return sxconst_set

	def construct_sxcmd_category_list(self):
		sxcmd_category_list = []
		sxcmd_list = []           # Used only within this function
		sxcmd_category_dict = {}  # Used only within this function

		# Actual configurations of all sx command categories and sx commands are inserted into the following section by wikiparser.py
		# as sxcmd_category_list and sxcmd_list
		# @@@@@ START_INSERTION @@@@@
		# @@@@@ END_INSERTION @@@@@

		# Create dictionaries from the constructed lists
		for sxcmd_category in sxcmd_category_list:
			sxcmd_category_dict[sxcmd_category.name] = sxcmd_category

		# Create command token dictionary for each SXcmd instance
		# Then, register SXcmd instance to an associated SXcmd_category
		for sxcmd in sxcmd_list:
			for sxcmd_token in sxcmd.token_list:
				# Handle very special cases
				if sxcmd_token.type == "user_func":
					n_widgets = 2 # function type has two line edit boxes
					sxcmd_token.label = [sxcmd_token.label, "Python script for user function"]
					sxcmd_token.help = [sxcmd_token.help, "Please leave it none/blank if file is not external to SPHIRE"]
					sxcmd_token.default = sxcmd_token.default.split('|||')
					if not sxcmd_token.is_locked:
						sxcmd_token.restore = [[sxcmd_token.default, sxcmd_token.default], [["none"], ["none"]]]
					else:
						sxcmd_token.restore = [["", ""], ["", ""]]
					sxcmd_token.default = [sxcmd_token.default[0], "none"]
				# else: Do nothing for the other types

				# Register this to command token dictionary
				sxcmd.token_dict[sxcmd_token.key_base] = sxcmd_token

			# Register this to command to command category dictionary
			assert sxcmd.category in sxcmd_category_dict, "sxcmd.category %s" % (sxcmd.category)
			sxcmd_category_dict[sxcmd.category].cmd_list.append(sxcmd)

		# Store the constructed lists and dictionary as a class data member
		return sxcmd_category_list

	def update_qsub_enable_states(self):
		# Construct and add widgets for sx command categories
		for idx in range(2):
			for sxcmd_category in self.sxcmd_category_list[idx]:
				# Create SXCmdCategoryWidget for this command category
				for sxcmd in sxcmd_category.cmd_list:
					sxcmd.widget.sxcmd_tab_main.set_qsub_enable_state()

	def handle_sxmenu_item_btn_event(self, sxmenu_item, i=0):
		assert(isinstance(sxmenu_item, SXmenu_item) == True) # Assuming the sxmenu_item is an instance of class SXmenu_item

		modifiers = QApplication.keyboardModifiers()
		if modifiers == Qt.ShiftModifier:
			sxmenu_item_wiki_url = SXLookFeelConst.generate_sxmenu_item_wiki_url(sxmenu_item)
			print("Opening Wiki Page ...")
			print(sxmenu_item_wiki_url)
			os.system("python -m webbrowser %s" % (sxmenu_item_wiki_url))
			return

		if self.cur_sxmenu_item == sxmenu_item: return

		if self.cur_sxmenu_item != None:
			self.cur_sxmenu_item.btn.setStyleSheet(self.cur_sxmenu_item.btn.customButtonStyle)

		if isinstance(sxmenu_item, SXconst_set):
			self.cur_sxmenu_item = self.sender().parent.sxconst_set
		elif isinstance(sxmenu_item, SXcmd_category):
			self.cur_sxmenu_item = self.sender().parent.sxcmd_category_list[i]
		else:
			self.cur_sxmenu_item = sxmenu_item

		if self.cur_sxmenu_item != None:
			self.cur_sxmenu_item.btn.setStyleSheet(self.cur_sxmenu_item.btn.customButtonStyleClicked)
			self.sxmenu_item_widget_stacked_layout[self.helical].setCurrentWidget(self.cur_sxmenu_item.widget)

	def closeEvent(self, event):
		event.ignore() # event.accept()

		# Quit child applications of all sxcmd widgets
		for idx in range(2):
			for sxcmd_category in self.sxcmd_category_list[idx]:
				sxcmd_category.widget.quit_all_child_applications()

		print("bye bye")
		QtCore.QCoreApplication.instance().quit()

	# def changeEvent(self, event):
	# 	print(self.frameGeometry())

# ========================================================================================

def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + """ 
	The main SPHIRE GUI application. It is designed as the command generator for the SPHIRE single particle analysis pipeline.
	"""
	parser = OptionParser(usage, version=SPARXVERSION)
	parser.add_option('--helical', action='store_true', default=False, help='Start the GUI in helical mode. This can be changed after the start. (default False)')
	# No options!!! Does not need to call parser.add_option()
	
	(options, args) = parser.parse_args(sys.argv[1:])
	
	if len(args) > 1:
		print("see usage " + usage)
		sys.exit()
	
	sxapp = QApplication(sys.argv)
	# The valid keys can be retrieved using the keys() function.
	# Typically they include "windows", "motif", "cde", "plastique" and "cleanlooks".
	# Depending on the platform, "windowsxp", "windowsvista" and "macintosh" may be available. Note that keys are case insensitive.
	# sxapp.setStyle("macintosh")
	if sys.platform.startswith('darwin'):
		sxapp.setStyle("macintosh")
	else:
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
	version_string = '1.2_rc9'
	SXLookFeelConst.initialise(sxapp, version_string)

	# Define the main window (class SXMainWindow)
	sxmain_window = SXMainWindow(options.helical)
	sxmain_window.setWindowTitle("SPHIRE-GUI Main Version {0}".format(version_string))
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
	main()

# ========================================================================================
# END OF SCRIPT
# ========================================================================================
