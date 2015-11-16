#!/usr/bin/env python
# 
#
# Author: Toshio Moriya 12/11/2015 (toshio.moriya@mpi-dortmund.mpg.de)
# ========================================================================================

# from EMAN2 import *
# from sparx import *
import os
import sys
from   global_def import ERROR

class SXcmd_token:
	def __init__(self):
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

class SXkeyword_map:
	def __init__(self, priority, cmd_token_type):
		if priority >= 100: ERROR("Priority should be lower than 100", "SXkeyword_map::__init__() in wikiparser")
		# class variables
		self.priority = priority              # Priority of this keyword. Highest priority is 0. The type of higher priority will be used to avoid the conflict among keywords
		self.cmd_token_type = cmd_token_type  # Command token value type 
		
def main():
		# The following should be data members of the popupwindow
		sxscript = ""
		short_info = ""
		mpi_support = False
		mpi_add_flag = False # NOTE: 2015/11/12 Toshio Moriya. This can be removed when --MPI flag is removed from all sx*.py scripts 
		cmd_token_list = [] # The list of command tokens. Need this to keep the order of command tokens
		cmd_token_dict = {} # The dictionary of command tokens, organised by key base name of command token. Easy to access a command token but looses their order
		
		# Open wiki document file 
		# This should be passed as an arguments
		work_dir = "%s/src/eman2/sparx/doc" % os.environ['EMAN2DIR']
		# file_path_wiki = "%s/viper.txt" % work_dir
		# file_path_wiki = "%s/window.txt" % work_dir
		# file_path_wiki = "%s/cter.txt" % work_dir
		# file_path_wiki = "%s/isac.txt" % work_dir
		# file_path_wiki = "%s/meridien.txt" % work_dir
		# file_path_wiki = "%s/3dvariability.txt" % work_dir
		# file_path_wiki = "%s/locres.txt" % work_dir
		# file_path_wiki = "%s/filterlocal.txt" % work_dir
		file_path_wiki = "%s/sort3d.txt" % work_dir
		
		# NOTE: 2015/11/11 Toshio Moriya
		# This should be exception. Need to decide if this should be skipped or exit system.
		if os.path.exists(file_path_wiki) == False: ERROR("Rutime Error: Wiki document is not found.", "main() in wikiparser")
			
		file_wiki = open(file_path_wiki,'r')
		
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
		# Use priority 2 for others
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
						sxscript = line_wiki[0:item_tail].strip()
						# Extract the short info about this sxscript (can be empty)
						short_info = line_wiki[item_tail + len(target_operator):].strip()							
					elif current_section == section_usage:
						# Extract 'usage in command line' to identify each command token is either an argument (no-prefix) or option ('--' prefix)
						# This information is also used to check consistency between 'usage in command line' and list in '== Input ==' and '== Output ==' sections
						if line_wiki[0:len("sx")] == "sx":
							usage_token_list = line_wiki.split()
							if usage_token_list[0] != sxscript + ".py": ERROR("Wiki Format Error: First token should be script name with .py (sx*.py)", "main() in wikiparser")
							# Register arguments and options
							for usage_token in usage_token_list[1:]:
								# Allocate memory for new command token
								cmd_token = SXcmd_token()
								# Extract key of command token. 
								key = usage_token.split("=")[0] # Remove characters after '=' if token contains it (i.e. some options)
								cmd_token.key_base = key.strip("-") # Get key base name by removing prefix ('--' or '-' for option)
								cmd_token.key_prefix = key[0:len(key) - len(cmd_token.key_base)]
								# Try to set the special type base on the keyword dictionary
								best_keyword_map = SXkeyword_map(99, "")
								for keyword in keyword_dict.keys():
									if cmd_token.key_base.find(keyword) != -1:
										# command token contains keyword
										keyword_map = keyword_dict[keyword]
										if best_keyword_map.priority > keyword_map.priority:
											# Update keyword_map to one with a higher priority
											best_keyword_map = keyword_map
								cmd_token.type = best_keyword_map.cmd_token_type # If command token does not contains any keywords, its type stays with ""
								# Register this command token to the list (ordered) and dictionary (unordered)		
								cmd_token_list.append(cmd_token)
								cmd_token_dict[cmd_token.key_base] = cmd_token
							# Check if --MPI is used in this script
							# NOTE: 2015/11/12 Toshio Moriya
							# The following can be removed when --MPI flag is removed from all sx*.py scripts 
							if "--MPI" in cmd_token_dict.keys():
								mpi_support = True
								mpi_add_flag = True
							# else: Do nothing.
						# else: Ignore this line (must be comments).
					elif current_section == section_typical:
						target_operator = "mpirun"
						if mpi_support == False and line_wiki.find(target_operator) != 1:
							mpi_support = True
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
								ERROR("Warning: This line (%s) is missing key base name (maybe comment line?). Ignoring this line...'."  % 	line_wiki, "main() in wikiparser", action = 0)
								continue
							key_base = line_buffer[0:item_tail]
							line_buffer = line_buffer[item_tail + len(target_operator):].strip() # Get the rest of line											
							# check consistency between 'usage in command line' and this
							if key_base not in cmd_token_dict.keys(): ERROR("Wiki Format Error: Key base (%s) is missing from 'usage in command line' in '= Usage ='." % key_base, "main() in wikiparser")
							# Get the reference to the command token object associated with this key base name
							cmd_token = cmd_token_dict[key_base]
							if cmd_token.key_base != key_base: ERROR("Logical Error: Registered command token with wrong key base name into the dictionary.", "main() in wikiparser")
							cmd_token.is_in_io = True # Set flag to tell this command token is find in input or output section
							cmd_token.group = group  # Set group of command token according to the current subsection
							# Extract label of command token
							target_operator = ":"
							item_tail = line_buffer.find(target_operator)
							if item_tail == -1: ERROR("Wiki Format Error: This line (%s) is missing label. Please check the format in Wiki document." % line_wiki, "main() in wikiparser")
							cmd_token.label = line_buffer[0:item_tail]
							line_buffer = line_buffer[item_tail + len(target_operator):].strip() # Get the rest of line	
							# Extract help of command token before default value
							target_operator = "(default"
							item_tail = line_buffer.find(target_operator)
							if item_tail == -1: ERROR("Wiki Format Error: This line (%s) is missing default setting. Please check the format in Wiki document." % line_wiki, "main() in wikiparser")
							cmd_token.help = line_buffer[0:item_tail]
							line_buffer = line_buffer[item_tail + len(target_operator):].strip() # Get the rest of line	
							# Extract default value of command token
							target_operator = ")"
							item_tail = line_buffer.find(target_operator)
							if item_tail == -1: ERROR("Wiki Format Error: This line (%s) is missing ')' for default setting. Please check the format in Wiki document." % line_wiki, "main() in wikiparser")
							default_value = line_buffer[0:item_tail].strip() # make sure spaces & new line are not included at head and tail
							if default_value.find("required") != -1:
								# This is a required command token and should have value type instead of default value
								cmd_token.is_required = True
								cmd_token.default = ""
								if not cmd_token.type:
									# Type is still empty, meaning no special type is assigned
									# Extract the data type (the rest of line)
									cmd_token.type = default_value.replace("required", "").strip()
							else: 
								# This is not required command token and should have default value
								cmd_token.is_required = False
								cmd_token.default = default_value
								
								if not cmd_token.type:
									# Type is still empty, meaning no special type is assigned
									# Find out the data type from default value
									try: 
										int(cmd_token.default)
										cmd_token.type = "int"
									except:
										try:  	 
											float(cmd_token.default)
											cmd_token.type = "float"
										except:  
											if cmd_token.default == "True" or cmd_token.default == "False":
												cmd_token.type = "bool"
											else:
												cmd_token.type = "string"
								# else: keep the special type
							# Ignore the rest of line ...
					else:
						ERROR("Logical Error: This section is invalid. Did you assigne an invalid section?", "main() in wikiparser")
					
		if current_state != state_done: ERROR("Wiki Format Error: parser could not extract all information necessary. Please check if the Wiki format has all required sections.", "main() in wikiparser")

		# Make sure there are no extra arguments or options in 'usage in command line' of '= Usage ='
		for cmd_token in cmd_token_list:
			if cmd_token.is_in_io == False: ERROR("Wiki Format Error: An extra argument or option (%s) is found in 'usage in command line' of '= Usage ='." % cmd_token.key_base, "main() in wikiparser")
					
		file_wiki.close()
		
		print "Succeed to parsing Wiki document (%s)" % file_path_wiki

		
		# For DEBUG
		print "><><>< DEBUG OUTPUT ><><><"
		print "------"
		print "GLOBAL"
		print "------"
		print "sxscript           : %s" % sxscript 
		print "short_info         : %s" % short_info 
		print "mpi_support        : %s" % mpi_support 
		print "mpi_add_flag       : %s" % mpi_add_flag 
		print "len(cmd_token_list): %d" % len(cmd_token_list)
		print "len(cmd_token_dict): %d" % len(cmd_token_dict)
		print "--------------"
		print "cmd_token_list"
		print "--------------"
		for cmd_token in cmd_token_list:
			print "%s%s (group=%s, required=%s, default=%s, type=%s) <%s>" % (cmd_token.key_prefix, cmd_token.key_base, cmd_token.group, cmd_token.is_required, cmd_token.default, cmd_token.type, cmd_token.label),cmd_token.help
		
if __name__ == '__main__':
	main()
