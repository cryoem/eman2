#!/usr/bin/env python
#
# Author: Toshio Moriya, 11/11/2015 (toshio.moriya@mpi-dortmund.mpg.de)
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
from global_def import *
from PyQt4.Qt import *
from PyQt4 import QtGui
from PyQt4 import QtCore
from subprocess import *
from EMAN2 import *
from sparx import *
from EMAN2_cppwrap import *

# ========================================================================================
class SXcmd_token:
	def __init__(self):
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.key_base = ""		    # key base name of command token (argument or option) in command line
		self.key_prefix = ""		# key prefix of of command token. None for argument, '--' or '-' for option
		self.label = ""				# User friendly name of argument or option
		self.help = ""		    	# Help info
		self.group = ""				# Tab group: main or advanced
		self.is_required = False	# Required argument or options. No default value are available 
		self.default = ""			# Default value
		self.type = ""		    	# Type of value
		self.is_in_io = False		# To check consistency between 'usage in command line' and list in '== Input ==' and '== Output ==' sections
		self.widget_list = []		# The list of associated widget instances to this command token
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

# ========================================================================================
class SXcmd:
	def __init__(self, wiki_file_path):
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.wiki_file_path = wiki_file_path  # File path to wiki documents of this command
		self.name = ""                        # Name of this command (i.e. name of sx*.py script but without .py extension)
		self.short_info = ""                  # Short description of this command
		self.mpi_support = False              # Flag to indicate if this command suppors MPI version
		self.mpi_add_flag = False             # NOTE: 2015/11/12 Toshio Moriya. This can be removed when --MPI flag is removed from all sx*.py scripts 
		self.token_list = []                  # list of command tokens. Need this to keep the order of command tokens
		self.token_dict = {}                  # dictionary of command tokens, organised by key base name of command token. Easy to access a command token but looses their order
		self.line = ""
		# self.button_widget = None		      # The list of associated widget instances to this command token
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		
		self.construct_token_list_from_wiki()

	# Set self.mpi_support and self.sxcmd_token_list by extracting info from option parser
	def construct_token_list_from_wiki(self):
		# Private helper class used only in this function
		class SXkeyword_map:
			def __init__(self, priority, token_type):
				if priority >= 100: ERROR("Priority should be lower than 100", "SXkeyword_map::__init__() in wikiparser")
				# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
				# class variables
				self.priority = priority      # Priority of this keyword. Highest priority is 0. The type of higher priority will be used to avoid the conflict among keywords
				self.token_type = token_type  # Token value type 
				# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
				
		# Define dictionary of keywords:
		# The dictionary maps command token to special data types
		# If a command token extracted from 'usage in command line' contains the keyword defined here
		# the associated special data type will be assigned to this command token.
		# 
		# - output      : Line edit box and output info button
		#                 GUI also checks the existence of output directory/file before execution of the sx*.py
		#                 GUI abort the execution if the directory/file exists already
		# - directory   : Line edit box and open directory button
		# - image       : Line edit box and open file buttons for .hdf and .bdb 
		# - parameters  : Line edit box and open file button for all file types 
		# - pdb         : Line edit box and open file button for .pdb 
		# - function    : Two line edit boxes (function name & file path of the container script)
		#                 and open file button for .py
		# 
		keyword_dict = {}
		# Use priority 0 to overrule the exceptional cases (This is a reason why priority is introduced...)
		keyword_dict["use_latest_master_directory"] = SXkeyword_map(0, "")           # --use_latest_master_directory (contains keyworkd 'directory' but this should be bool type)
		# Use priority 1 for output
		keyword_dict["output"]                      = SXkeyword_map(1, "output")     # output.hdf, output_directory, outputfile
		keyword_dict["outdir"]                      = SXkeyword_map(1, "output")     # outdir1, outdir2, outdir, --outdir=output_directory
		keyword_dict["locresvolume"]                = SXkeyword_map(1, "output")     # locresvolume (this contained keyword "volume" also... This is another reason why priority is introduced...)
		keyword_dict["directory"]                   = SXkeyword_map(1, "output")     # directory
		# Use priority 2 for the others
		keyword_dict["indir"]                       = SXkeyword_map(2, "directory")  # --indir=input_directory
		keyword_dict["coords_dir"]                  = SXkeyword_map(2, "directory")  # --coords_dir=coords_directory
		keyword_dict["stack"]                       = SXkeyword_map(2, "image")      # stack, stack_file, prj_stack
		keyword_dict["volume"]                      = SXkeyword_map(2, "image")      # initial_volume, firstvolume, secondvolume, inputvolume
		keyword_dict["mask"]                        = SXkeyword_map(2, "image")      # --mask3D=mask3D, maskfile, mask
		keyword_dict["focus"]                       = SXkeyword_map(2, "image")      # --focus=3Dmask
		keyword_dict["importctf"]                   = SXkeyword_map(2, "paramters")  # --importctf=ctf_file
		keyword_dict["pwreference"]                 = SXkeyword_map(2, "paramters")  # --pwreference=pwreference_file
		keyword_dict["pdb"]                         = SXkeyword_map(2, "pdb")        # input.pdb
		keyword_dict["function"]                    = SXkeyword_map(2, "function")   # --function=user_function
				
		# Define list of target sections for GUI and set current
		section_lists = []		
		section_lists.append("= Name ="); section_name = len(section_lists) - 1; 
		section_lists.append("= Usage ="); section_usage = len(section_lists) - 1; 
		section_lists.append("=== Typical usage ==="); section_typical = len(section_lists) - 1; 
		section_lists.append("== Input =="); section_input = len(section_lists) - 1; 
		section_lists.append("== Output =="); section_output = len(section_lists) - 1; 
		current_section = section_name
		
		# Define list of subsections of input section and set current		
		group_main = "main"
		group_advanced = "advanced"
		group = group_main
				
		# Define States and set current
		state_searching  = 0
		state_processing = 1
		state_done = 1
		current_state = state_searching
				
		# NOTE: 2015/11/11 Toshio Moriya
		# This should be exception. Need to decide if this should be skipped or exit system.
		if os.path.exists(self.wiki_file_path) == False: ERROR("Rutime Error: Wiki document is not found.", "main() in wikiparser")
		
		file_wiki = open(self.wiki_file_path,'r')
		
		# Loop through all lines in the wiki document file
		for line_wiki in file_wiki:
			# make sure spaces & new line are not included at head and tail of this line
			line_wiki = line_wiki.strip()  

			if not line_wiki:
				# This is empty line. Always ignore it regardless of state
				continue
				
			if current_state == state_searching:
				if line_wiki.find(section_lists[current_section]) != -1:
					# Found the current target section
					current_state = state_processing
				# else: just ignore this line 
			else:
				if current_state != state_processing: ERROR("Logical Error: This condition should not happen! State setting must be incorrect.", "main() in wikiparser")
				if line_wiki[0] == "=": # Assuming the section always starts with "="
					# Reached the next section (might be not target)
					current_section += 1 # Update current target section
					group = group_main   # reset group (subsection) for '== Input ==' and '== Output ==' sections
					if current_section == len(section_lists):
						# All target sections are handled
						current_state = state_done
						break
						
					if line_wiki.find(section_lists[current_section]) != -1:
						# Found the current target section
						if current_section >= len(section_lists): ERROR("Logical Error: This condition should not happen! Section setting must be incorrect.", "main() in wikiparser")
						current_state = state_processing
					else:
						# This section was not the current target. Go back to searching state
						current_state = state_searching											
				else:
					# We are in a target section
					if current_section == section_name:
						# Extract the name of sxscript
						target_operator = "-"
						item_tail = line_wiki.find(target_operator)
						if item_tail == -1: ERROR("Wiki Format Error: '= Name =' section should contain only one valid line, and the line should starts from 'sx* - ': %s" % line_wiki, "main() in wikiparser")
						self.name = line_wiki[0:item_tail].strip()
						# Extract the short info about this sxscript (can be empty)
						self.short_info = line_wiki[item_tail + len(target_operator):].strip()							
					elif current_section == section_usage:
						# Extract 'usage in command line' to identify each command token is either an argument (no-prefix) or option ('--' prefix)
						# This information is also used to check consistency between 'usage in command line' and list in '== Input ==' and '== Output ==' sections
						if line_wiki[0:len("sx")] == "sx":
							usage_token_list = line_wiki.split()
							if usage_token_list[0] != self.name + ".py": ERROR("Wiki Format Error: First token should be script name with .py (sx*.py)", "main() in wikiparser")
							# Register arguments and options
							for usage_token in usage_token_list[1:]:
								# Check if --MPI is used in this script
								# NOTE: 2015/11/12 Toshio Moriya
								# The following can be removed when --MPI flag is removed from all sx*.py scripts 
								if usage_token == "--MPI":
									# ERROR("Warning: The 'usage in command line' contains --MPI flag. The flag will be removed in near future, so ignoring this line...'.", "main() in wikiparser", action = 0)
									self.mpi_support = True
									self.mpi_add_flag = True
									continue
								# Allocate memory for new command token
								token = SXcmd_token()
								# Extract key of command token. 
								key = usage_token.split("=")[0] # Remove characters after '=' if token contains it (i.e. some options)
								token.key_base = key.strip("-") # Get key base name by removing prefix ('--' or '-' for option)
								token.key_prefix = key[0:len(key) - len(token.key_base)]
								# Try to set the special type base on the keyword dictionary
								best_keyword_map = SXkeyword_map(99, "")
								for keyword in keyword_dict.keys():
									if token.key_base.find(keyword) != -1:
										# command token contains keyword
										keyword_map = keyword_dict[keyword]
										if best_keyword_map.priority > keyword_map.priority:
											# Update keyword_map to one with a higher priority
											best_keyword_map = keyword_map
								token.type = best_keyword_map.token_type # If command token does not contains any keywords, its type stays with ""
								# Register this command token to the list (ordered) and dictionary (unordered)		
								self.token_list.append(token)
								self.token_dict[token.key_base] = token
						# else: Ignore this line (must be comments).
					elif current_section == section_typical:
						target_operator = "mpirun"
						if self.mpi_support == False and line_wiki.find(target_operator) > 1:
							self.mpi_support = True
						# else: Ignore this line
					elif current_section == section_input or current_section == section_output:
						if line_wiki[0] == "*" and line_wiki.find("optional"):
							# Reached the option subsection (argument subsection is done)
							group = group_advanced
						else:
							line_buffer = line_wiki
							# Extract key base name of command token
							target_operator = "::"
							item_tail = line_buffer.find(target_operator)
							if item_tail == -1: 
								# ERROR("Warning: This line (%s) is missing key base name (maybe comment line?). Ignoring this line...'."  % (line_wiki), "main() in wikiparser", action = 0)
								continue
							key_base = line_buffer[0:item_tail]
							if key_base == "MPI":
								# ERROR("Warning: This line (%s) contains MPI flag. The flag will be removed in near future, so ignoring this line...'."  % (line_wiki), "main() in wikiparser", action = 0)
								if self.mpi_support == False or self.mpi_add_flag == False: ERROR("Logical Error: Since MPI flag is found and the command should support MPI.", "main() in wikiparser")
								continue
							line_buffer = line_buffer[item_tail + len(target_operator):].strip() # Get the rest of line											
							# check consistency between 'usage in command line' and this
							if key_base not in self.token_dict.keys(): ERROR("Wiki Format Error: Key base (%s) is missing from 'usage in command line' in '= Usage ='." % key_base, "main() in wikiparser")
							# Get the reference to the command token object associated with this key base name
							token = self.token_dict[key_base]
							if token.key_base != key_base: ERROR("Logical Error: Registered command token with wrong key base name into the dictionary.", "main() in wikiparser")
							token.is_in_io = True # Set flag to tell this command token is find in input or output section
							token.group = group  # Set group of command token according to the current subsection
							# Extract label of command token
							target_operator = ":"
							item_tail = line_buffer.find(target_operator)
							if item_tail == -1: ERROR("Wiki Format Error: This line (%s) is missing label. Please check the format in Wiki document." % line_wiki, "main() in wikiparser")
							token.label = line_buffer[0:item_tail]
							line_buffer = line_buffer[item_tail + len(target_operator):].strip() # Get the rest of line	
							# Extract help of command token before default value
							target_operator = "(default"
							item_tail = line_buffer.find(target_operator)
							if item_tail == -1: ERROR("Wiki Format Error: This line (%s) is missing default setting. Please check the format in Wiki document." % line_wiki, "main() in wikiparser")
							token.help = line_buffer[0:item_tail]
							line_buffer = line_buffer[item_tail + len(target_operator):].strip() # Get the rest of line	
							# Extract default value of command token
							target_operator = ")"
							item_tail = line_buffer.find(target_operator)
							if item_tail == -1: ERROR("Wiki Format Error: This line (%s) is missing ')' for default setting. Please check the format in Wiki document." % line_wiki, "main() in wikiparser")
							default_value = line_buffer[0:item_tail].strip() # make sure spaces & new line are not included at head and tail
							if default_value.find("required") != -1:
								# This is a required command token and should have value type instead of default value
								token.is_required = True
								token.default = ""
								if not token.type:
									# Type is still empty, meaning no special type is assigned
									# Extract the data type (the rest of line)
									token.type = default_value.replace("required", "").strip()
							else: 
								# This is not required command token and should have default value
								token.is_required = False
								token.default = default_value
								
								if not token.type:
									# Type is still empty, meaning no special type is assigned
									# Find out the data type from default value
									try: 
										int(token.default)
										token.type = "int"
									except:
										try:  	 
											float(token.default)
											token.type = "float"
										except:  
											if token.default == "True":
												token.default = True # convert the default value to boolean
												token.type = "bool"
											elif token.default == "False":
												token.default = False # convert the default value to boolean
												token.type = "bool"
											else:
												token.type = "string"
								# else: keep the special type
							# Ignore the rest of line ...
					else:
						ERROR("Logical Error: This section is invalid. Did you assigne an invalid section?", "main() in wikiparser")
					
		if current_state != state_done: ERROR("Wiki Format Error: parser could not extract all information necessary. Please check if the Wiki format has all required sections.", "main() in wikiparser")

		# Make sure there are no extra arguments or options in 'usage in command line' of '= Usage ='
		for token in self.token_list:
			if token.is_in_io == False: ERROR("Wiki Format Error: An extra argument or option (%s) is found in 'usage in command line' of '= Usage ='." % token.key_base, "main() in wikiparser")
				
		file_wiki.close()
		
		print "Succeed to parsing Wiki document (%s)" % self.wiki_file_path
		"""
		# For DEBUG
		if self.name == "sxwindow": 
			print "><><>< DEBUG OUTPUT ><><><"
			print ""
			print "------"
			print "GLOBAL"
			print "------"
			print "name            : %s" % self.name 
			print "short_info      : %s" % self.short_info 
			print "mpi_support     : %s" % self.mpi_support 
			print "mpi_add_flag    : %s" % self.mpi_add_flag 
			print "len(token_list) : %d" % len(self.token_list)
			print "len(token_dict) : %d" % len(self.token_dict)
			print ""
			print "--------------"
			print "cmd_token_list"
			print "--------------"
			for token in self.token_list:
				print "%s%s (group=%s, required=%s, default=%s, type=%s) <%s>" % (token.key_prefix, token.key_base, token.group, token.is_required, token.default, token.type, token.label),token.help
			print ""
		"""

	def generate_cmd_line(self):
	
		self.map_gui_to_cmd_line()
		
		file_name_out = QtGui.QFileDialog.getSaveFileName(self, "Generate Command Line", options = QtGui.QFileDialog.DontUseNativeDialog)
		if file_name_out != '':
			file_out = open(file_name_out,'w')
			file_out.write(self.sxcmd_line + '\n')
			file_out.close()
			print 'Saved the following command to %s:' % file_name_out
			print self.sxcmd_line

	
# ========================================================================================
# Provides all necessary functionarity
# tabs only contains gui and knows how to layout them
class SXPopup(QWidget):
	def __init__(self, sxcmd):
		QWidget.__init__(self)
		
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.sxcmd = sxcmd
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		
		self.setWindowTitle(self.sxcmd.name)
		self.tab_main = SXTab_main(self)
		self.tab_advance = SXTab_advance(self)
		self.tab_main.w1 = self.tab_advance
		self.TabWidget = QtGui.QTabWidget()
		self.TabWidget.insertTab(0,self.tab_main,'Main')
		self.TabWidget.insertTab(1,self.tab_advance,'Advanced')
		self.TabWidget.resize(730,1080) # self.TabWidget.resize(730,860)
		self.TabWidget.show()
		
	def map_gui_to_cmd_line(self):
		# Add program name to command line
		sxcmd_line = "%s.py" % self.name
		
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
	
	"""
	def save_params(self):		
		file_name_out = QtGui.QFileDialog.getSaveFileName(self, "Save Parameters", options = QtGui.QFileDialog.DontUseNativeDialog)
		if file_name_out != '':
			file_out = open(file_name_out,'w')
		
			# Write script name for consistency check upon loading
			file_out.write('%s \n' % (self.name))	

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
	"""

	"""
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
				QMessageBox.warning(self, 'Fail to load paramters', 'The specified file is not paramter file for %s.' % self.name)
	
			file_in.close()
	"""

	"""
#	def choose_file(self):
#		#opens a file browser, showing files only in .hdf format
#		file_name = QtGui.QFileDialog.getOpenFileName(self, "Open Data File", "", "HDF files (*.hdf)")
#		#after the user selected a file, we obtain this filename as a Qstring
#		a=QtCore.QString(file_name)
#		print a
#		#we convert this Qstring to a string and send it to line edit classed stackname edit of the Poptwodali window
#		self.stacknameedit.setText(str(a))
	"""
		
	"""
#		#Function choose_file started when  the  open_file of the  Poptwodali window is clicked (same as above but for bdb files(maybe we can combine these two into one function)
#	def choose_file1(self):
#		file_name1 = QtGui.QFileDialog.getOpenFileName(self, "Open Data File", "EMAN2DB/", "BDB FILES (*.bdb)" )
#		a=QtCore.QString(file_name1)
#		b=os.path.basename(str(a))
#		c=os.path.splitext(b)[0]
#		d="bdb:"+c
#		print d
#		self.stacknameedit.setText(d)
	"""
	
	"""
#	def show_output_info(self):
#		QMessageBox.information(self, "sx* output",'outdir is the name of the output folder specified by the user. If it does not exist, the directory will be created. If it does exist, the program will crash and an error message will come up. Please change the name of directory and restart the program.')
	"""
		
# ========================================================================================
class SXTab_main(QWidget):

	def __init__(self, parent = None):
		QWidget.__init__(self, parent)
		
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.sxcmd = parent.sxcmd
		
		# layout parameters
		self.y1 = 10
		self.y2 = self.y1 + 95 #self.y2 = self.y1 + 98
		self.y4 = self.y2 + 450
		self.y5 = self.y4 + 95
		self.yspc = 4
		
		self.x1 = 10
		self.x2 = self.x1 + 500 # self.x2 = self.x1 + 200
		self.x3 = self.x2 + 145
		self.x4 = self.x3 + 100
		self.x5 = 230
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		
		# # Set the window title
		# self.setWindowTitle(self.sxcmd.name)
		# Set a label and its position in this tab
		temp_label = QtGui.QLabel('<b>%s</b>' % (self.sxcmd.name), self)
		temp_label.move(self.x1,self.y1)
		# NOTE: 2015/11/17 Toshio Moriya
		# Necessary to separate '<b>%s</b>' from the information for avoiding to invoke the tag interpretations of string
		# e.g. < becomes the escape character
		temp_label = QtGui.QLabel('%s' % (self.sxcmd.short_info), self)
		temp_label.setWordWrap(True)
		temp_label.setFixedWidth(600)
		temp_label.move(self.x1 + 100, self.y1)
		self.y1 += 50
		
		# Add load paramater button 
		self.load_params_btn = QPushButton("Load parameters", self)
		self.load_params_btn.move(self.x1-5,self.y1)
		self.load_params_btn.setToolTip('Load gui parameter settings to retrieve a previously-saved one')
#		self.connect(self.load_params_btn, SIGNAL("clicked()"), parent.load_params)
		
		# Add widget for editing command args and options
		for cmd_token in self.sxcmd.token_list:
			if cmd_token.group == 'main':			
				# Create label widget 
				label_widget = QtGui.QLabel(cmd_token.label, self)
				label_widget.move(self.x1,self.y2)
				# Create widget and associate it to this cmd_token 
				cmd_token_widget = None
				if cmd_token.type == "bool":
					# construct new widget(s) for this command token
					cmd_token_widget = QtGui.QCheckBox("", self)
					cmd_token_widget.setCheckState(cmd_token.default)
				else:
					# For now, use line edit box for the other type
					cmd_token_widget = QtGui.QLineEdit(self)
					cmd_token_widget.setText(cmd_token.default)
					# if cmd_token.type == "output":
					# elif cmd_token.type == "directory":
					# elif cmd_token.type == "image":
					# elif cmd_token.type == "parameters":
					# elif cmd_token.type == "pdb":
					# elif cmd_token.type == "function":
					# else:
					#	if cmd_token.type not in ["int", "float", "string"]: ERROR("Logical Error: Encountered unsupported type (%s). Consult with the developer."  % 	line_wiki, "SXTab_main::__init__() in sxgui.py")
				cmd_token_widget.move(self.x2,self.y2 - 7)
				cmd_token_widget.setToolTip(cmd_token.help)		
				cmd_token.widget_list.append(cmd_token_widget)
				
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

		temp_label = QtGui.QLabel('MPI command line template', self)
		temp_label.move(self.x1,self.y2)
		self.mpi_cmd_line_edit = QtGui.QLineEdit(self)
		self.mpi_cmd_line_edit.setText('')
		self.mpi_cmd_line_edit.move(self.x2,self.y2)
		self.mpi_cmd_line_edit.setToolTip('The template of MPI command line (e.g. "mpirun -np XXX_SXMPI_NPROC_XXX --host n0,n1,n2 XXX_SXCMD_LINE_XXX"). If empty, use "mpirun -np XXX_SXMPI_NPROC_XXX XXX_SXCMD_LINE_XXX"')

		self.y2 =self.y2+25
		
		# If MPI is not supported, disable this widget
		self.set_widget_enable_state(self.mpi_nproc_edit, self.sxcmd.mpi_support)
		self.set_widget_enable_state(self.mpi_cmd_line_edit, self.sxcmd.mpi_support)

		# Add gui components for queue submission (qsub)
		is_qsub_enabled = False
		temp_label = QtGui.QLabel('submit job to queue', self)
		temp_label.move(self.x1,self.y2)
		self.qsub_enable_checkbox = QtGui.QCheckBox("", self)
		self.qsub_enable_checkbox.setCheckState(is_qsub_enabled)
		self.qsub_enable_checkbox.stateChanged.connect(self.set_qsub_enable_state) # To control enable state of the following qsub related widgets
		self.qsub_enable_checkbox.move(self.x2,self.y2)
		self.qsub_enable_checkbox.setToolTip('submit job to queue')
		
		self.y2 =self.y2+25
		
		temp_label = QtGui.QLabel('job name', self)
		temp_label.move(self.x1,self.y2)
		self.qsub_job_name_edit = QtGui.QLineEdit(self)
		self.qsub_job_name_edit.setText(self.sxcmd.name)
		self.qsub_job_name_edit.move(self.x2,self.y2)
		self.qsub_job_name_edit.setToolTip('name of this job')

		self.y2 =self.y2+25

		temp_label = QtGui.QLabel('submission command', self)
		temp_label.move(self.x1,self.y2)
		self.qsub_cmd_edit = QtGui.QLineEdit(self)
		self.qsub_cmd_edit.setText('qsub')
		self.qsub_cmd_edit.move(self.x2,self.y2)
		self.qsub_cmd_edit.setToolTip('name of submission command to queue job')

		self.y2 =self.y2+25

		temp_label = QtGui.QLabel('submission script template', self)
		temp_label.move(self.x1,self.y2)
		self.qsub_script_edit = QtGui.QLineEdit(self)
		self.qsub_script_edit.setText('msgui_qsub.sh')
		self.qsub_script_edit.move(self.x2,self.y2)
		self.qsub_script_edit.setToolTip('file name of submission script template (e.g. $EMAN2DIR/bin/msgui_qsub.sh')

		self.y2 =self.y2+25
		
		# Initialize enable state of qsub related widgets
		self.set_qsub_enable_state()
		
		# Add space
		self.y2 = self.y2+25*1

		# Add save paramater button 
		self.save_params_btn = QPushButton("Save parameters", self)
		# self.save_params_btn.move(self.x1-5,  self.y4)
		self.save_params_btn.move(self.x1-5,  self.y2)
		self.save_params_btn.setToolTip('Save gui parameter settings')
#		self.connect(self.save_params_btn, SIGNAL("clicked()"), parent.save_params)
		
		# self.y4 = self.y4+30
		self.y2 = self.y2+30

		self.cmd_line_btn = QPushButton("Generate command line", self)
		# self.cmd_line_btn.move(self.x1-5,  self.y4)
		self.cmd_line_btn.move(self.x1-5,  self.y2)
		self.cmd_line_btn.setToolTip('Generate command line from gui parameter settings')
#		self.connect(self.cmd_line_btn, SIGNAL("clicked()"), parent.generate_cmd_line)
		
		self.y2 = self.y2+30

		# Add a run button
		self.execute_btn = QtGui.QPushButton('Run %s' % self.sxcmd.name, self)
		# make 3D textured push button look
		s = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px}"
		self.execute_btn.setStyleSheet(s)
		# self.execute_btn.move(self.x5,  self.y5)
		self.execute_btn.move(self.x5,  self.y2)
#		self.connect(self.execute_btn, SIGNAL("clicked()"), parent.execute_cmd)

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
		if self.sxcmd.mpi_support:
			self.set_widget_enable_state(self.mpi_cmd_line_edit, not is_enabled)
		
		# Set enable state and background color of qsub related widgets
		self.set_widget_enable_state(self.qsub_job_name_edit, is_enabled)
		self.set_widget_enable_state(self.qsub_cmd_edit, is_enabled)
		self.set_widget_enable_state(self.qsub_script_edit, is_enabled)
		
# ========================================================================================
class SXTab_advance(QWidget):
	def __init__(self, parent = None):
		QWidget.__init__(self)
				
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.sxcmd = parent.sxcmd
		
		# layout parameters
		self.y1=10
		self.yspc = 4
		
		self.x1 = 20
		self.x2 = self.x1 + 500 # self.x2 = self.x1+280
		self.x3 = self.x2 + 145
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		
		# Set the window title
		#self.setWindowTitle('%s advanced parameter selection' % self.sxcmd.name)
		# Set a label and its position in this tab
		temp_label = QtGui.QLabel('<b>%s</b>' % (self.sxcmd.name), self)
		temp_label.move(self.x1,self.y1)
		temp_label = QtGui.QLabel('Set advanced parameters', self)
		temp_label.setWordWrap(True)
		temp_label.setFixedWidth(600)
		temp_label.move(self.x1 + 100, self.y1)
		self.y1 = self.y1+25
		
		# Add gui components for editing command args and options
		for cmd_token in self.sxcmd.token_list:
			if cmd_token.group == 'advanced':
				# Create label widget
				label_widget = QtGui.QLabel(cmd_token.label, self)
				label_widget.move(self.x1,self.y1)		
				# Create widget and associate it to this cmd_token 
				cmd_token_widget = None
				if cmd_token.type == "bool":
					# construct new widget(s) for this command token
					cmd_token_widget = QtGui.QCheckBox("", self)
					cmd_token_widget.setCheckState(cmd_token.default)
				else:
					# For now, use line edit box for the other type
					cmd_token_widget = QtGui.QLineEdit(self)
					cmd_token_widget.setText(cmd_token.default)
					# if cmd_token.type == "output":
					# elif cmd_token.type == "directory":
					# elif cmd_token.type == "image":
					# elif cmd_token.type == "parameters":
					# elif cmd_token.type == "pdb":
					# elif cmd_token.type == "function":
					# else:
					#	if cmd_token.type not in ["int", "float", "string"]: ERROR("Logical Error: Encountered unsupported type (%s). Consult with the developer."  % 	line_wiki, "SXTab_main::__init__() in sxgui.py")
				cmd_token_widget.move(self.x2,self.y1)
				cmd_token_widget.setToolTip(cmd_token.help)		
				cmd_token.widget_list.append(cmd_token_widget)
					
				self.y1 = self.y1+25
# sx_end

# ========================================================================================
# Layout of the Pop Up window SXPopup_info; started by the function info of the main window
class SXPopup_info(QWidget):
	def __init__(self):
		QWidget.__init__(self)
		#Here we just set the window title and  3 different labels, with their positions in the window
		self.setWindowTitle('Sparx GUI Info Page')
		title1=QtGui.QLabel('<b>Sparx GUI Prototype</b>', self)
		title1.move(10,10)
		title2=QtGui.QLabel('<b>Authors: Toshio Moriya</b> ', self)
		title2.move(10,40)
		title3=QtGui.QLabel('For more information visit:\n%s ' % SPARX_DOCUMENTATION_WEBSITE, self)
		title3.move(10,70)

# ========================================================================================
# Main Window (started by class App)
# This class includes the layout of the main window			
class MainWindow(QtGui.QWidget):
	def __init__(self):
		QtGui.QWidget.__init__(self)
		
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
		# Get all necessary informations from wiki documents of sx*.py scripts
		# --------------------------------------------------------------------------------
		# NOTE: 2015/11/17 Toshio Moriya
		# (1) The path should be corrected, needs to follow the installation
		# (2) It is desirable to read the wiki documents file path from a setting file for more flexibility 
		work_dir = "%s/src/eman2/sparx/doc" % os.environ['EMAN2DIR']
		wiki_file_path_list = []
		
		wiki_file_path_list.append("%s/cter.txt" % work_dir)
		wiki_file_path_list.append("%s/window.txt" % work_dir)
		wiki_file_path_list.append("%s/isac.txt" % work_dir)
		wiki_file_path_list.append("%s/viper.txt" % work_dir)
		wiki_file_path_list.append("%s/meridien.txt" % work_dir)
		wiki_file_path_list.append("%s/3dvariability.txt" % work_dir)
		wiki_file_path_list.append("%s/locres.txt" % work_dir)
		wiki_file_path_list.append("%s/filterlocal.txt" % work_dir)
		wiki_file_path_list.append("%s/sort3d.txt" % work_dir)
		
		self.y2 = 95
		
		for wiki_file_path in wiki_file_path_list:
			# Construct sxscript object associated with this wiki document
			sxcmd = SXcmd(wiki_file_path)
			
			# Add buttons for this sx*.py processe
			temp_btn = QPushButton(sxcmd.name, self)
			temp_btn.move(10, self.y2)
			temp_btn.setToolTip(sxcmd.short_info)
			from functools import partial
			self.connect(temp_btn, SIGNAL("clicked()"), partial(self.handle_sxcmd_btn_event, sxcmd))

			self.y2 += 30
			
		# Set the width and height of the main window
		self.resize(300,400)

	# Click actions: The following functions are associated with the click event of push buttons (btn##) on the main window. 
	def handle_sxcmd_btn_event(self, sxcmd):
		modifiers = QtGui.QApplication.keyboardModifiers()
		if modifiers == QtCore.Qt.ShiftModifier:
			os.system("python -m webbrowser %s%s" % (SPARX_DOCUMENTATION_WEBSITE, sxcmd.name))
			return
			
		self.w = SXPopup(sxcmd)
		
	#This is the function info, which is being started when the Pushbutton btn_info of the main window is being clicked
	def info(self):
		# print "Opening a new popup window..."
		# Opens the window SXPopup_info, and defines its width and height
		# The layout of the SXPopup_info window is defined in class SXPopup_info(QWidget Window)
		self.w = SXPopup_info()
		self.w.resize(250,200)
		self.w.show()

# ========================================================================================
#  This is the main class of the program
#  Here we provide the necessary imports. The basic GUI widgets are located in QtGui module.
class App(QApplication):
	def __init__(self, *args):
		QApplication.__init__(self, *args)
		# Define the main window (class MainWindow)
		self.main = MainWindow()
		# self.main.resize(400,450)
		self.main.resize(1000, 755)
		# Define that when all windows are closed, function byebye of class App will be started
		self.connect(self, SIGNAL("lastWindowClosed()"), self.byebye )
		# Show main window
		self.main.show()
		
	#function byebye (just quit)  
	def byebye( self ):
		print' bye bye!'
		self.exit(0)

# ========================================================================================
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

# ========================================================================================
if __name__ == "__main__":
	main(sys.argv)

# ========================================================================================
# END OF SCRIPT
# ========================================================================================

