#!/usr/bin/env python
# 
#
# Author: Toshio Moriya 12/11/2015 (toshio.moriya@mpi-dortmund.mpg.de)
# ========================================================================================

# from EMAN2 import *
# from sparx import *
import os
import sys
import copy
from global_def import ERROR

from sxgui_template import SXcmd_token, SXcmd

# ========================================================================================
class SXsubcmd_config:
	def __init__(self, label = "", token_edit_list = [], mpi_support = None):
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.label = label                              # User friendly name of command subset
		self.token_edit_list = token_edit_list          # To edit some attributes of tokens (label, help, group, is_required, default). If the original value should be kept, set to None. First entry should be mode token.
		self.mpi_support = mpi_support                  # Flag to indicate if this command suppors MPI version
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

# ========================================================================================
class SXcmd_config:
	def __init__(self, wiki, category, is_submittable = True, exclude_list = [], subconfig = None):
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.wiki = wiki                      # Wiki document file path
		self.category = category              # Category of this command: pipe (pipeline), util (utility)
		self.is_submittable = is_submittable  # External GUI Application (e.g. sxgui_cter.py) should not be submitted to job queue
		self.exclude_list = exclude_list      # token key base list to be excluded
		self.subconfig = subconfig            # Subset configuration of this command (e.g. sxprocess and sxlocres). None includes all command tokens, and do not make subcmd
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

# ========================================================================================
def construct_token_list_from_wiki(sxcmd_config):
	# Private helper class used only in this function
	class SXkeyword_map:
		def __init__(self, priority, token_type):
			if priority >= 100: ERROR("Priority should be lower than 100", "%s in %s" % (__name__, os.path.basename(__file__)))
			# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
			# class variables
			self.priority = priority      # Priority of this keyword. Highest priority is 0. The type of higher priority will be used to avoid the conflict among keywords
			self.token_type = token_type  # Token value type 
			# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	
	print "Start parsing Wiki document (%s as %s command) " % (sxcmd_config.wiki, sxcmd_config.category)
	
	# Allocate memory for new SXcmd instance
	sxcmd = SXcmd(sxcmd_config.category, sxcmd_config.is_submittable)
	
	# Define dictionary of keywords:
	# The dictionary maps command token to special data types
	# If a command token extracted from 'usage in command line' contains the keyword defined here
	# the associated special data type will be assigned to this command token.
	# 
	# - output      : Line edit box for formatted string type, and output info button. 
	#                 GUI also checks the existence of output directory/file before execution of the sx*.py
	#                 GUI abort the execution if the directory/file exists already
	# - image       : Line edit box for formatted string type, and open file buttons for .hdf and .bdb 
	# - any_image   : Line edit box for formatted string type, and open file buttons for all file types (also mrc, tiff, and etc) and .bdb 
	# - parameters  : Line edit box for formatted string type, and open file button for all file types 
	# - any_file    : Line edit box for formatted string type, and open file button for all file types 
	# - bdb         : Line edit box for formatted string type, and open file button for .bdb 
	# - pdb         : Line edit box for formatted string type, and open file button for .pdb 
	# - function    : Two line edit boxes for formatted string type (function name & file path of the container script), 
	#                 and open file button for .py
	# - directory   : Line edit box for formatted string type, and open directory button (NOTE: Toshio Moriya 2016/03/03: Not used at this point)
	# 
	# - apix        : Project constant - float type
	# - wn          : Project constant - int type
	# - box         : Project constant - int type
	# - radius      : Project constant - int type
	# - sym         : Project constant - formatted string type
	#
	
	keyword_dict = {}
	
	# Use priority 0 to overrule the exceptional cases (This is a reason why priority is introduced...)
	keyword_dict["--use_latest_master_directory"] = SXkeyword_map(0, "")           # --use_latest_master_directory (contains keyworkd 'directory' but this should be bool type)
	keyword_dict["stack_file"]                    = SXkeyword_map(0, "bdb")        # stack_file (contains keyworkd 'stack' but this should be bdb type)
	keyword_dict["--stack_mode"]                  = SXkeyword_map(0, "")           # stack_mode (contains keyworkd 'stack' but this should be bool type)
	keyword_dict["--adaptive_mask"]               = SXkeyword_map(0, "")           # --adaptive_mask (contains keyworkd 'mask' but this should be bool type)
	keyword_dict["--symmetrize"]                  = SXkeyword_map(0, "")           # --symmetrize (contains keyworkd '--sym' but this should be bool type)
	# Use priority 1 for output
	keyword_dict["output"]                        = SXkeyword_map(1, "output")     # output_hdf, output_directory, outputfile, outputfile, --output=OUTPUT
	keyword_dict["outdir"]                        = SXkeyword_map(1, "output")     # outdir
	keyword_dict["locres_volume"]                 = SXkeyword_map(1, "output")     # locres_volume (this contained keyword "volume" also... This is another reason why priority is introduced...)
	keyword_dict["directory"]                     = SXkeyword_map(1, "output")     # directory
	keyword_dict["rotpw"]                         = SXkeyword_map(1, "output")     # rotpw
	keyword_dict["output_mask3D"]                 = SXkeyword_map(1, "output")     # output_mask3D
	# Use priority 2 for the others
	keyword_dict["stack"]                         = SXkeyword_map(2, "image")      # stack, prj_stack
	keyword_dict["volume"]                        = SXkeyword_map(2, "image")      # initial_volume, firstvolume, secondvolume, input_volume
	keyword_dict["mask"]                          = SXkeyword_map(2, "image")      # --mask3D=mask3D, maskfile, mask, --mask=MASK
	keyword_dict["--focus"]                       = SXkeyword_map(2, "image")      # --focus=3Dmask
	keyword_dict["--input"]                       = SXkeyword_map(2, "image")      # --input=INPUT
	keyword_dict["input_micrograph"]              = SXkeyword_map(2, "any_image")  # input_micrograph_pattern
	keyword_dict["input_image"]                   = SXkeyword_map(2, "any_image")  # input_image
	keyword_dict["--tr0"]                         = SXkeyword_map(2, "parameters") # --tr0=matrix_file
	keyword_dict["input_coordinates"]             = SXkeyword_map(2, "parameters") # input_coordinates_pattern
	keyword_dict["--import_ctf"]                  = SXkeyword_map(2, "parameters") # --import_ctf=ctf_file
	keyword_dict["--importctf"]                   = SXkeyword_map(2, "parameters") # --importctf=IMPORTCTF
	keyword_dict["cter_ctf_file"]                 = SXkeyword_map(2, "parameters") # cter_ctf_file 
	keyword_dict["--pwreference"]                 = SXkeyword_map(2, "parameters") # --pwreference=pwreference 
	keyword_dict["inputfile"]                     = SXkeyword_map(2, "any_file")   # inputfile
	keyword_dict["input_pdb"]                     = SXkeyword_map(2, "pdb")        # input_pdb
	keyword_dict["--function"]                    = SXkeyword_map(2, "function")   # --function=user_function
	
	keyword_dict["--apix"]                        = SXkeyword_map(2, "apix")       # --apix=pixel_size, --apix 
	keyword_dict["--pixel_size"]                  = SXkeyword_map(2, "apix")       # --pixel_size=PIXEL_SIZE
	keyword_dict["--wn"]                          = SXkeyword_map(2, "ctfwin")     # --wn
	keyword_dict["--box"]                         = SXkeyword_map(2, "box")        # --box=box_size, --box_size=box_size
	keyword_dict["--radius"]                      = SXkeyword_map(2, "radius")     # --radius=particle_radius, --radius=outer_radius, --radius=outer_radius, --radius=particle_radius, --radius=outer_radius, --radius=outer_radius
	keyword_dict["--sym"]                         = SXkeyword_map(2, "sym")        # --sym=c1, --sym=c1, --sym=c1, --sym=symmetry, --sym=c1, --sym=c4
	
	# NOTE: 2016/02/23 Toshio Moriya
	# Below might be useful to include
	# reference power spectrum? --pwreference of viper, --pwreference of rviper, --PWadjustment of sort3d, --PWadjustment of rsort3d
	#
	# Below must be exceptional cases
	# --wn of locres, sort3d, & rsort3d; same as ctfwin?
	# --radiusvar of 3dvariability; same as radius?
	# --radius of locres & filterlocal; same as radius? 
	#
	
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
	state_done = 2
	current_state = state_searching
	
	# NOTE: 2015/11/11 Toshio Moriya
	# This should be exception. Need to decide if this should be skipped or exit system.
	if os.path.exists(sxcmd_config.wiki) == False: ERROR("Rutime Error: Wiki document is not found.", "%s in %s" % (__name__, os.path.basename(__file__)))
	
	file_wiki = open(sxcmd_config.wiki,'r')
	
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
			if current_state != state_processing: ERROR("Logical Error: This condition should not happen! State setting must be incorrect.", "%s in %s" % (__name__, os.path.basename(__file__)))
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
					if current_section >= len(section_lists): ERROR("Logical Error: This condition should not happen! Section setting must be incorrect.", "%s in %s" % (__name__, os.path.basename(__file__)))
					current_state = state_processing
				else:
					# This section was not the current target. Go back to searching state
					current_state = state_searching
			else:
				# We are in a target section
				if current_section == section_name:
					line_buffer = line_wiki
					# Extract the name of sxscript
					target_operator = "-"
					item_tail = line_buffer.find(target_operator)
					if item_tail == -1: ERROR("Wiki Format Error: '= Name =' section should contain only one valid line, and the line should starts from 'sx* - ': %s" % line_wiki, "%s in %s" % (__name__, os.path.basename(__file__)))
					sxcmd.name = line_buffer[0:item_tail].strip()
					line_buffer = line_buffer[item_tail + len(target_operator):].strip() # Get the rest of line
					# Extract the label of this sxscript
					target_operator = ":"
					item_tail = line_buffer.find(target_operator)
					if item_tail == -1: ERROR("Wiki Format Error: '= Name =' section should contain a label ended with ':' after 'sx* - ': %s" % line_wiki, "%s in %s" % (__name__, os.path.basename(__file__)))
					sxcmd.label = line_buffer[0:item_tail].strip()
					# Extract the short info about this sxscript (can be empty)
					sxcmd.short_info = line_buffer[item_tail + len(target_operator):].strip() # Get the rest of line
				elif current_section == section_usage:
					# Extract 'usage in command line' to identify each command token is either an argument (no-prefix) or option ('--' prefix)
					# This information is also used to check consistency between 'usage in command line' and list in '== Input ==' and '== Output ==' sections
					if line_wiki[0:len("sx")] == "sx":
						usage_token_list = line_wiki.split()
						if usage_token_list[0] != sxcmd.name + ".py": ERROR("Wiki Format Error: First token should be script name with .py (sx*.py)", "%s in %s" % (__name__, os.path.basename(__file__)))
						# Register arguments and options
						for usage_token in usage_token_list[1:]:
							# Check if --MPI is used in this script
							# NOTE: 2015/11/12 Toshio Moriya
							# The following can be removed when --MPI flag is removed from all sx*.py scripts 
							if usage_token == "--MPI":
								# ERROR("Warning: The 'usage in command line' contains --MPI flag. The flag will be removed in near future, so ignoring this line...'.","%s in %s" % (__name__, os.path.basename(__file__)), action = 0)
								sxcmd.mpi_support = True
								sxcmd.mpi_add_flag = True
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
								if key.find(keyword) != -1:
									# command token contains keyword
									keyword_map = keyword_dict[keyword]
									if best_keyword_map.priority > keyword_map.priority:
										# Update keyword_map to one with a higher priority
										best_keyword_map = keyword_map
							token.type = best_keyword_map.token_type # If command token does not contains any keywords, its type stays with ""
							# Register this command token to the list (ordered) and dictionary (unordered)		
							sxcmd.token_list.append(token)
							sxcmd.token_dict[token.key_base] = token
					# else: Ignore this line (must be comments).
				elif current_section == section_typical:
					target_operator = "mpirun"
					if sxcmd.mpi_support == False and line_wiki.find(target_operator) > 1 and line_wiki.find(sxcmd.name + ".py") > 1:
						sxcmd.mpi_support = True
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
							# ERROR("Warning: This line (%s) is missing key base name (maybe comment line?). Ignoring this line...'."  % (line_wiki),"%s in %s" % (__name__, os.path.basename(__file__)), action = 0)
							continue
						key_base = line_buffer[0:item_tail]
						if key_base == "MPI":
							# ERROR("Warning: This line (%s) contains MPI flag. The flag will be removed in near future, so ignoring this line...'."  % (line_wiki), "%s in %s" % (__name__, os.path.basename(__file__)), action = 0)
							if sxcmd.mpi_support == False or sxcmd.mpi_add_flag == False: ERROR("Logical Error: Since MPI flag is found, the command should support MPI.", "%s in %s" % (__name__, os.path.basename(__file__)))
							continue
						line_buffer = line_buffer[item_tail + len(target_operator):].strip() # Get the rest of line
						# check consistency between 'usage in command line' and this
						if key_base not in sxcmd.token_dict.keys(): ERROR("Wiki Format Error: Key base (%s) is missing from 'usage in command line' in '= Usage ='." % key_base, "%s in %s" % (__name__, os.path.basename(__file__)))
						# Get the reference to the command token object associated with this key base name
						token = sxcmd.token_dict[key_base]
						if token.key_base != key_base: ERROR("Logical Error: Registered command token with wrong key base name into the dictionary.", "%s in %s" % (__name__, os.path.basename(__file__)))
						token.is_in_io = True # Set flag to tell this command token is find in input or output section
						token.group = group  # Set group of command token according to the current subsection
						# Extract label of command token
						target_operator = ":"
						item_tail = line_buffer.find(target_operator)
						if item_tail == -1: ERROR("Wiki Format Error: This line (%s) is missing label. Please check the format in Wiki document." % line_wiki, "%s in %s" % (__name__, os.path.basename(__file__)))
						token.label = line_buffer[0:item_tail]
						line_buffer = line_buffer[item_tail + len(target_operator):].strip() # Get the rest of line	
						# Extract help of command token before default value
						target_operator = "(default"
						item_tail = line_buffer.find(target_operator)
						if item_tail == -1: ERROR("Wiki Format Error: This line (%s) is missing default setting. Please check the format in Wiki document." % line_wiki, "%s in %s" % (__name__, os.path.basename(__file__)))
						token.help = line_buffer[0:item_tail]
						line_buffer = line_buffer[item_tail + len(target_operator):].strip() # Get the rest of line	
						# Extract default value of command token
						target_operator = ")"
						item_tail = line_buffer.find(target_operator)
						if item_tail == -1: ERROR("Wiki Format Error: This line (%s) is missing ')' for default setting. Please check the format in Wiki document." % line_wiki,"%s in %s" % (__name__, os.path.basename(__file__)))
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
						# Initialise restore value with default value
						token.restore = token.default
						# Ignore the rest of line ...
				else:
					ERROR("Logical Error: This section is invalid. Did you assigne an invalid section?", "%s in %s" % (__name__, os.path.basename(__file__)))
				
	if current_state != state_done: ERROR("Wiki Format Error: parser could not extract all information necessary. Please check if the Wiki format has all required sections.", "%s in %s" % (__name__, os.path.basename(__file__)))

	# Make sure there are no extra arguments or options in 'usage in command line' of '= Usage ='
	for token in sxcmd.token_list:
		if token.is_in_io == False: ERROR("Wiki Format Error: An extra argument or option (%s) is found in 'usage in command line' of '= Usage ='." % token.key_base, "%s in %s" % (__name__, os.path.basename(__file__)))
	
	file_wiki.close()
	
	# DESIGN_NOTE: 2016/02/05 Toshio Moriya
	# Handle exceptional cases due to the limitation of software design 
	# In future, we should remove these exception handling by reviewing the design
	if sxcmd.name == "sxfilterlocal":
		assert(sxcmd.token_dict["locres_volume"].key_base == "locres_volume")
		assert(sxcmd.token_dict["locres_volume"].type == "output")
		sxcmd.token_dict["locres_volume"].type = "image"
	elif sxcmd.name in ["sxlocres",  "sxsort3d", "sxrsort3d"]:
		assert(sxcmd.token_dict["wn"].key_base == "wn")
		assert(sxcmd.token_dict["wn"].type == "ctfwin")
		sxcmd.token_dict["wn"].type = "int"
	
	print "Succeed to parse Wiki document (%s as %s command)" % (sxcmd_config.wiki, sxcmd_config.category)
	
	"""
	# For DEBUG
	if sxcmd.name == "sxwindow": 
		print "><><>< DEBUG OUTPUT ><><><"
		print ""
		print "------"
		print "GLOBAL"
		print "------"
		print "name            : %s" % sxcmd.name 
		print "label           : %s" % sxcmd.label 
		print "short_info      : %s" % sxcmd.short_info 
		print "mpi_support     : %s" % sxcmd.mpi_support 
		print "mpi_add_flag    : %s" % sxcmd.mpi_add_flag 
		print "type            : %s" % sxcmd.type 
		print "len(token_list) : %d" % len(sxcmd.token_list)
		print "len(token_dict) : %d" % len(sxcmd.token_dict)
		print ""
		print "--------------"
		print "cmd_token_list"
		print "--------------"
		for token in sxcmd.token_list:
			print "%s%s (group=%s, required=%s, default=%s, type=%s, restore=%s) <%s>" % (token.key_prefix, token.key_base, token.group, token.is_required, token.default, token.type, token.restore, token.label, token.help)
		print ""
	"""
	
	return sxcmd

def apply_exclude_list(sxcmd_config_exclude_list, sxcmd):
	assert(len(sxcmd_config_exclude_list) > 0)
	assert(len(sxcmd.token_list) == len(sxcmd.token_dict))
	
	for token_key_base in sxcmd_config_exclude_list:
		sxcmd.token_list.remove(sxcmd.token_dict[token_key_base])
		del sxcmd.token_dict[token_key_base]

def apply_sxsubcmd_config(sxsubcmd_config, sxcmd):
	assert(sxsubcmd_config != None)
	assert(len(sxsubcmd_config.token_edit_list) > 0)
	
	# Copy command token dictionary, then clear the command token list and dictionary
	fullset_token_dict = copy.deepcopy(sxcmd.token_dict)
	sxcmd.token_list = []
	sxcmd.token_dict = {}
	
	# Using the first entry in token edit list as command mode token of this subset, 
	# get mode token from sxcmd (having a fullset of tokens)
	mode_token_edit = sxsubcmd_config.token_edit_list[0]
	if mode_token_edit.key_base not in fullset_token_dict.keys(): ERROR("Logical Error: This condition should not happen! Subset command configuration must be incorrect.", "%s in %s" % (__name__, os.path.basename(__file__)))
	mode_token = fullset_token_dict[mode_token_edit.key_base]
	
	# Create mode name of this subset, append key base of mode token to mode_name of this command
	sxcmd.mode = mode_token.key_base
	# print "MRK_DEBUG: sxcmd.mode = %s" % (sxcmd.mode)
	
	# Set command label of this subset
	sxcmd.label = sxsubcmd_config.label
	# print "MRK_DEBUG: sxcmd.label = %s" % (sxcmd.label)
	
	# Set command mpi support of this subset if necessary
	if sxsubcmd_config.mpi_support != None:
		sxcmd.mpi_support = sxsubcmd_config.mpi_support
		if sxcmd.mpi_support == False:
			sxcmd.mpi_add_flag = False
	# print "MRK_DEBUG: sxcmd.mpi_support = %s" % (sxcmd.mpi_support)
	
	# Use label of mode token as a short info of subset command
	sxcmd.short_info = "%s. %s" % (mode_token.label, mode_token.help)
	# print "MRK_DEBUG: sxcmd.short_info = %s" % (sxcmd.short_info)
	
	# Reconstruct token list
	for token_edit in sxsubcmd_config.token_edit_list:
		# print "MRK_DEBUG: token_edit.key_base = %s" % (token_edit.key_base)
		token = None
		if token_edit.key_base not in fullset_token_dict.keys():
			# token key base is not found in fullset. This must be an argument to be added
			if token_edit.key_prefix != "": ERROR("Logical Error: This condition should not happen! Subset command configuration must be incorrect.", "%s in %s" % (__name__, os.path.basename(__file__)))
			token = token_edit
		else:
			# token key base is found in fullset. This must be an option.
			token = fullset_token_dict[token_edit.key_base]
			assert(token_edit.key_prefix == None)
			if token_edit.label != None:
				token.label = token_edit.label
			if token_edit.help != None:
				token.help = token_edit.help
			if token_edit.group != None:
				token.group = token_edit.group
			if token_edit.is_required != None:
				token.is_required = token_edit.is_required
			if token_edit.default != None:
				token.default = token_edit.default
			if token_edit.type != None:
				token.type = token_edit.type
		assert(token != None)
		token.restore = token.default
		sxcmd.token_list.append(token)
		sxcmd.token_dict[token_edit.key_base] = (token)
		
	assert(len(sxcmd.token_list) == len(sxsubcmd_config.token_edit_list))
	assert(len(sxcmd.token_dict) == len(sxsubcmd_config.token_edit_list))

def insert_sxcmd_to_file(sxcmd, output_file, sxcmd_variable_name):
	output_file.write("\t")
	output_file.write("%s = SXcmd()" % sxcmd_variable_name)
	output_file.write("; %s.name = \"%s\"" % (sxcmd_variable_name, sxcmd.name)) 
	output_file.write("; %s.mode = \"%s\"" % (sxcmd_variable_name, sxcmd.mode)) 
	output_file.write("; %s.label = \"%s\"" % (sxcmd_variable_name, sxcmd.label)) 
	output_file.write("; %s.short_info = \"%s\"" % (sxcmd_variable_name, sxcmd.short_info.replace("\"", "'")))
	output_file.write("; %s.mpi_support = %s" % (sxcmd_variable_name, sxcmd.mpi_support))
	output_file.write("; %s.mpi_add_flag = %s" % (sxcmd_variable_name, sxcmd.mpi_add_flag))
	output_file.write("; %s.category = \"%s\"" % (sxcmd_variable_name, sxcmd.category))
	output_file.write("; %s.is_submittable = %s" % (sxcmd_variable_name, sxcmd.is_submittable))
	output_file.write("\n")
	
	for token in sxcmd.token_list:
		output_file.write("\t")
		output_file.write("token = SXcmd_token()")
		output_file.write("; token.key_base = \"%s\"" % token.key_base)
		output_file.write("; token.key_prefix = \"%s\"" % token.key_prefix)
		output_file.write("; token.label = \"%s\"" % token.label)
		output_file.write("; token.help = \"%s\"" % token.help.replace("\"", "'")) 
		output_file.write("; token.group = \"%s\"" % token.group)
		output_file.write("; token.is_required = %s" % token.is_required)
		if token.type == "bool":
			output_file.write("; token.default = %s" % token.default)
			output_file.write("; token.restore = %s" % token.restore)
		else:
			if token.is_required:
				output_file.write("; token.default = \"\"")
				output_file.write("; token.restore = \"\"")
			else:
				output_file.write("; token.default = \"%s\"" % token.default)
				output_file.write("; token.restore = \"%s\"" % token.restore)
		output_file.write("; token.type = \"%s\"" % token.type)
		# output_file.write("; token.is_in_io = %s" % token.is_in_io)
		
		output_file.write("; %s.token_list.append(token)" % sxcmd_variable_name)
		output_file.write("\n")
	
	return 
		
# ========================================================================================
def main():
	# --------------------------------------------------------------------------------
	# Get all necessary informations from wiki documents of sx*.py scripts
	# and create gui generation paramter
	# --------------------------------------------------------------------------------
	sxcmd_config_list = []
	
	# --------------------------------------------------------------------------------
	# Define pipeline command settings
	# --------------------------------------------------------------------------------
	sxcmd_config_list.append(SXcmd_config("../doc/cter.txt", "pipe", exclude_list=["stack_mode"]))
	
	sxcmd_config_list.append(SXcmd_config("../doc/gui_cter.txt", "pipe", is_submittable = False))
	
	sxcmd_config_list.append(SXcmd_config("../doc/window.txt", "pipe"))
	
	sxcmd_config_list.append(SXcmd_config("../doc/isac.txt", "pipe"))
	# sxcmd_config_list.append(SXcmd_config("../doc/isac_snr4.txt", "pipe"))
	
	sxcmd_config_list.append(SXcmd_config("../doc/isac_post_processing.txt", "pipe"))
	
	sxcmd_config_list.append(SXcmd_config("../doc/viper.txt", "pipe"))

	# NOTE: Toshio Moriya 2016/03/11
	# Temporarily disabled sxrviper for the 03/07/2016 release
	# sxcmd_config_list.append(SXcmd_config("../doc/rviper.txt", "pipe"))
	
	sxcmd_config_list.append(SXcmd_config("../doc/meridien.txt", "pipe"))
	
	
	token_edit_list = []
	token_edit = SXcmd_token(); token_edit.initialize_edit("postprocess"); token_edit.is_required = True; token_edit.default = True; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("firstvolume"); token_edit.key_prefix = ""; token_edit.label = "first unfiltered half-volume "; token_edit.help = ""; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "image"; token_edit_list.append(token_edit) 
	token_edit = SXcmd_token(); token_edit.initialize_edit("secondvolume"); token_edit.key_prefix = ""; token_edit.label = "second unfiltered half-volume "; token_edit.help = ""; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "image"; token_edit_list.append(token_edit) 
	token_edit = SXcmd_token(); token_edit.initialize_edit("fsc_weighted"); token_edit.help = "main"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("low_pass_filter"); token_edit.help = "main"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("ff"); token_edit.help = "main"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("aa"); token_edit.help = "main"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("mask"); token_edit.help = "main"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("output"); token_edit.help = "main"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("pixel_size"); token_edit.help = "main"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("B_start"); token_edit.help = "main"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("FSC_cutoff"); token_edit.help = "main"; token_edit_list.append(token_edit)
	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("3D Refinement Postprocess", token_edit_list, sxsubcmd_mpi_support)
	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "pipe", subconfig = sxcmd_subconfig))

	token_edit_list = []
	token_edit = SXcmd_token(); token_edit.initialize_edit("symmetrize"); token_edit.is_required = True; token_edit.default = True; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("input_volume"); token_edit.key_prefix = ""; token_edit.label = "input volume"; token_edit.help = ""; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "image"; token_edit_list.append(token_edit) 
	token_edit = SXcmd_token(); token_edit.initialize_edit("sym"); token_edit.help = "main"; token_edit_list.append(token_edit)
	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("3D Variability Preprocess", token_edit_list, sxsubcmd_mpi_support)
	sxcmd_config_list.append(SXcmd_config("../doc/3dvariability.txt", "pipe", subconfig = sxcmd_subconfig))
	
	sxcmd_config_list.append(SXcmd_config("../doc/3dvariability.txt", "pipe", exclude_list=["symmetrize"]))
	
	sxcmd_config_list.append(SXcmd_config("../doc/locres.txt", "pipe"))

	sxcmd_config_list.append(SXcmd_config("../doc/sort3d.txt", "pipe"))

	sxcmd_config_list.append(SXcmd_config("../doc/rsort3d.txt", "pipe"))

	# --------------------------------------------------------------------------------
	# Define utility command settings
	# --------------------------------------------------------------------------------
	sxcmd_config_list.append(SXcmd_config("../doc/pdb2em.txt", "util"))

	token_edit_list = []
	token_edit = SXcmd_token(); token_edit.initialize_edit("adaptive_mask"); token_edit.is_required = True; token_edit.default = True; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("input_volume"); token_edit.key_prefix = ""; token_edit.label = "input volume"; token_edit.help = ""; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "image"; token_edit_list.append(token_edit) 
	token_edit = SXcmd_token(); token_edit.initialize_edit("output_mask3D"); token_edit.key_prefix = ""; token_edit.label = "output 3D mask"; token_edit.help = ""; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "output"; token_edit_list.append(token_edit) 
	token_edit = SXcmd_token(); token_edit.initialize_edit("nsigma"); token_edit.help = "main"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("ndilation"); token_edit.help = "main"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("kernel_size"); token_edit.help = "main"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("gauss_standard_dev"); token_edit.help = "main"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("threshold"); token_edit.help = "main"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("ne"); token_edit.help = "main"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("nd"); token_edit.help = "main"; token_edit_list.append(token_edit)
	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("Adaptive 3D Mask", token_edit_list, sxsubcmd_mpi_support)
	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "util", subconfig = sxcmd_subconfig))
	
	sxcmd_config_list.append(SXcmd_config("../doc/filterlocal.txt", "util"))

	# sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "util"))
	
#	token_edit_list = []
#	token_edit = SXcmd_token(); token_edit.initialize_edit("stack_mode"); token_edit.group = "main"; token_edit.is_required = True; token_edit.default = True; token_edit_list.append(token_edit)
#	token_edit = SXcmd_token(); token_edit.initialize_edit("stack"); token_edit.key_prefix = ""; token_edit.label = "2D images in a stack file (bdb or hdf)"; token_edit.help = ""; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "image"; token_edit_list.append(token_edit) 
#	token_edit = SXcmd_token(); token_edit.initialize_edit("output_directory"); token_edit_list.append(token_edit) 
#	token_edit = SXcmd_token(); token_edit.initialize_edit("apix"); token_edit_list.append(token_edit)
#	token_edit = SXcmd_token(); token_edit.initialize_edit("Cs"); token_edit_list.append(token_edit)
#	token_edit = SXcmd_token(); token_edit.initialize_edit("voltage"); token_edit_list.append(token_edit)
#	token_edit = SXcmd_token(); token_edit.initialize_edit("ac"); token_edit_list.append(token_edit)
#	token_edit = SXcmd_token(); token_edit.initialize_edit("f_start"); token_edit_list.append(token_edit)
#	token_edit = SXcmd_token(); token_edit.initialize_edit("f_stop"); token_edit_list.append(token_edit)
#	# token_edit = SXcmd_token(); token_edit.initialize_edit("kboot"); token_edit_list.append(token_edit)
#	# token_edit = SXcmd_token(); token_edit.initialize_edit("overlap_x"); token_edit_list.append(token_edit)
#	# token_edit = SXcmd_token(); token_edit.initialize_edit("overlap_y"); token_edit_list.append(token_edit)
#	# token_edit = SXcmd_token(); token_edit.initialize_edit("edge_x"); token_edit_list.append(token_edit)
#	# token_edit = SXcmd_token(); token_edit.initialize_edit("edge_y"); token_edit_list.append(token_edit)
#	token_edit = SXcmd_token(); token_edit.initialize_edit("debug"); token_edit_list.append(token_edit)
#	sxsubcmd_mpi_support = False
#	sxcmd_subconfig = SXsubcmd_config("CTF Estimation (Stack Mode)", token_edit_list, sxsubcmd_mpi_support)
#	sxcmd_config_list.append(SXcmd_config("../doc/cter.txt", "util", subconfig = sxcmd_subconfig))

	# --------------------------------------------------------------------------------
	# Generate sxgui.py
	# --------------------------------------------------------------------------------
	sxgui_template_file_path = "sxgui_template.py"
	
	output_file_path = "../bin/sxgui.py" # output_file_path = "sxgui_trial.py"
	# remove the previous output
	if os.path.exists(output_file_path):
		os.remove(output_file_path)
	
	sxgui_template_file = open(sxgui_template_file_path,'r')
	output_file = open(output_file_path,'w')
	
	# Define States and set current
	state_template  = 0
	state_insertion = 1
	current_state = state_template
	
	for line in sxgui_template_file:
		output_file.write(line)
		if current_state == state_template:
			if line.find("# @@@@@ START_INSERTION @@@@@") != -1:
				current_state = state_insertion
				sxcmd_variable_name = "sxcmd"
				for sxcmd_config in sxcmd_config_list:
					# Construct sxscript object associated with this wiki document
					sxcmd = construct_token_list_from_wiki(sxcmd_config)
					if sxcmd_config.subconfig != None:
						apply_sxsubcmd_config(sxcmd_config.subconfig, sxcmd)
					if len(sxcmd_config.exclude_list) > 0:
						apply_exclude_list(sxcmd_config.exclude_list, sxcmd)
					insert_sxcmd_to_file(sxcmd, output_file, sxcmd_variable_name)
					output_file.write("\n")
					output_file.write("\tsxcmd_list.append(%s)\n" % sxcmd_variable_name)
					output_file.write("\n")
			# else: do nothing
		else:
			if current_state != state_insertion: ERROR("Logical Error: This condition should not happen! State setting must be incorrect.", "%s in %s" % (__name__, os.path.basename(__file__)))
			if line.find("# @@@@@ END_INSERTION @@@@@") != -1:
				current_state = state_template
			# else: do nothing
			
	if current_state == state_insertion: ERROR("Script Template Format Error: START_INSERTION and END_INSERTION must be paired.", "%s in %s" % (__name__, os.path.basename(__file__)))
	
	output_file.close()
	
	os.system("chmod +x %s" % output_file_path)
	
# ========================================================================================
if __name__ == '__main__':
	main()

# ========================================================================================
# END OF SCRIPT
# ========================================================================================
