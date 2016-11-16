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

from sxgui_template import SXcmd_token, SXcmd, SXcmd_category

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
	def __init__(self, wiki, format, category, role, is_submittable = True, exclude_list = [], subconfig = None):
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.wiki = wiki                      # Wiki document file path
		self.format = format                  # Wiki document format: MoinMoinWiki or DokuWiki
		self.category = category              # Category of this command: sxc_movie_micrograph, sxc_ctf, sxc_particle_stack, sxc_2d_clustering, sxc_initial_3d_modeling, sxc_3d_refinement, sxc_3d_clustering, sxc_utilities
		self.role = role                      # Role of this command; sxr_pipe (pipeline), sxr_alt (alternative) sxr_util (utility)
		self.is_submittable = is_submittable  # External GUI Application (e.g. sxgui_cter.py) should not be submitted to job queue. If it is true, this command will be registered to child_application_list of SXCmdWidget
		self.exclude_list = exclude_list      # token key base list to be excluded
		self.subconfig = subconfig            # Subset configuration of this command (e.g. sxprocess and sxlocres). None includes all command tokens, and do not make subcmd
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><


# ========================================================================================
# Helper class used only in construct_token_list_from_*() functions
class SXkeyword_map:
	def __init__(self, priority, token_type):
		if priority >= 100: ERROR("Priority should be lower than 100", "%s in %s" % (__name__, os.path.basename(__file__)))
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.priority = priority      # Priority of this keyword. Highest priority is 0. The type of higher priority will be used to avoid the conflict among keywords
		self.token_type = token_type  # Token value type
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		
# ----------------------------------------------------------------------------------------
def construct_keyword_dict():
	# Define dictionary of keywords:
	# The dictionary maps command token to special data types
	# If a command token extracted from 'usage in command line' contains the keyword defined here
	# the associated special data type will be assigned to this command token.
	#
	# - output          : Line edit box for formatted string type, and output info button.
	#                     GUI also checks the existence of output directory/file before execution of the sx*.py
	#                     GUI abort the execution if the directory/file exists already
	# - image           : Line edit box for formatted string type, and open file buttons for .hdf and .bdb
	# - any_image       : Line edit box for formatted string type, and open file buttons for all file types (also mrc, tiff, and etc) and .bdb
	# - any_micrograph  : Line edit box for formatted string type, and open file buttons for all file types (also mrc, tiff, and etc) and .txt
	# - parameters      : Line edit box for formatted string type, and open file button for all file types
	# - any_file        : Line edit box for formatted string type, and open file button for all file types
	# - bdb             : Line edit box for formatted string type, and open file button for .bdb
	# - hdf             : Line edit box for formatted string type, and open file button for .hdf
	# - pdb             : Line edit box for formatted string type, and open file button for .pdb
	# - mrc             : Line edit box for formatted string type, and open file button for .mrc
	# - txt             : Line edit box for formatted string type, and open file button for .txt
	# - exe             : Line edit box for formatted string type, and open file button for .exe and no extension
	# - any_file_list   : Line edit box for formatted string type, and open file button for all file types and .bdb
	#                     The string with space is interpreted as a list of any image file names upon command generation. (i.e. does not enclose the string with single quotes)
	# - any_image_list  : Line edit box for formatted string type, and open file button for all file types (also mrc, tiff, and etc) and .bdb.
	#                     The string with space is interpreted as a list of any image file names upon command generation. (i.e. does not enclose the string with single quotes)
	# - function        : Two line edit boxes for formatted string type (function name & file path of the container script),
	#                     and open file button for .py
	# - directory       : Line edit box for formatted string type, and open directory button
	# - any_directory   : Line edit box for formatted string type, and open directory button
	#                     The string with space is interpreted as a list of any directory names upon command generation. (i.e. does not enclose the string with single quotes)
	#
	# - apix            : Project constant - float typeJ78
	# - ctfwin          : Project constant - int type
	# - box             : Project constant - int type
	# - radius          : Project constant - int type
	# - sym             : Project constant - formatted string type
	#

	keyword_dict = {}

	# Use priority 0 to overrule the exceptional cases (This is a reason why priority is introduced...)
	keyword_dict["--use_latest_master_directory"] = SXkeyword_map(0, "")               # --use_latest_master_directory (contains keyworkd 'directory' but this should be bool type)
	keyword_dict["--stack_mode"]                  = SXkeyword_map(0, "")               # stack_mode (contains keyworkd 'stack' but this should be bool type)
	keyword_dict["--adaptive_mask"]               = SXkeyword_map(0, "")               # --adaptive_mask (contains keyworkd 'mask' but this should be bool type)
	keyword_dict["--binary_mask"]                 = SXkeyword_map(0, "")               # --binary_mask (contains keyworkd 'mask' but this should be bool type)
	keyword_dict["--symmetrize"]                  = SXkeyword_map(0, "")               # --symmetrize (contains keyworkd '--sym' but this should be bool type)
#	keyword_dict["input_micrograph_list_file"]    = SXkeyword_map(0, "txt")            # input_micrograph_list_file (contains keyworkd 'input_micrograph_list' but this should be parameters type)
	keyword_dict["isac_directory"]                = SXkeyword_map(0, "directory")      # isac_directory (contains keyworkd 'directory' but this should be directory type)
	keyword_dict["--single_stack_output"]         = SXkeyword_map(0, "bool")           # --single_stack_output (contains keyworkd 'output' but this should be bool type)
	keyword_dict["--hardmask"]                    = SXkeyword_map(0, "bool")           # --hardmask (contains keyworkd 'mask' but this should be bool type)
	keyword_dict["--do_adaptive_mask"]            = SXkeyword_map(0, "bool")           # --do_adaptive_mask (contains keyworkd 'mask' but this should be bool type)
	# Use priority 1 for output
	keyword_dict["output"]                        = SXkeyword_map(1, "output")         # output_hdf, output_directory, outputfile, outputfile, --output=OUTPUT, output_stack, output_file
	keyword_dict["outdir"]                        = SXkeyword_map(1, "output")         # outdir
	keyword_dict["locres_volume"]                 = SXkeyword_map(1, "output")         # locres_volume (this contained keyword "volume" also... This is another reason why priority is introduced...)
	keyword_dict["directory"]                     = SXkeyword_map(1, "output")         # directory
	keyword_dict["rotpw"]                         = SXkeyword_map(1, "output")         # rotpw
	keyword_dict["output_mask3D"]                 = SXkeyword_map(1, "output")         # output_mask3D
	keyword_dict["--makevstack"]                  = SXkeyword_map(1, "output")         # --makevstack
	keyword_dict["input_micrograph_list"]         = SXkeyword_map(1, "any_image_list") # input_micrograph_list (contains keyworkd 'input_micrograph' but this should be image_list type)
	# Use priority 2 for the others
	keyword_dict["stack"]                         = SXkeyword_map(2, "image")          # stack, prj_stack, input_stack
	keyword_dict["volume"]                        = SXkeyword_map(2, "image")          # initial_volume, firstvolume, secondvolume, input_volume
	keyword_dict["mask"]                          = SXkeyword_map(2, "image")          # --mask3D=mask3D, maskfile, mask, --mask=MASK
	keyword_dict["--focus"]                       = SXkeyword_map(2, "image")          # --focus=3Dmask
	keyword_dict["--input"]                       = SXkeyword_map(2, "image")          # --input=INPUT
	keyword_dict["class_file_name_no_dir_info"]   = SXkeyword_map(2, "image")          # class_file_name_no_dir_info
	keyword_dict["input_micrograph"]              = SXkeyword_map(2, "any_image")      # input_micrograph_pattern
	keyword_dict["input_image"]                   = SXkeyword_map(2, "any_micrograph") # input_image
	keyword_dict["selection_list"]                = SXkeyword_map(2, "any_micrograph") # selection_list
	keyword_dict["--tr0"]                         = SXkeyword_map(2, "parameters")     # --tr0=matrix_file
	keyword_dict["input_coordinates"]             = SXkeyword_map(2, "parameters")     # input_coordinates_pattern
	keyword_dict["input_ctf_params_source"]       = SXkeyword_map(2, "parameters")     # input_ctf_params_source
	keyword_dict["--importctf"]                   = SXkeyword_map(2, "parameters")     # --importctf=IMPORTCTF
	keyword_dict["--pwreference"]                 = SXkeyword_map(2, "parameters")     # --pwreference=pwreference
	keyword_dict["--PWadjustment"]                = SXkeyword_map(2, "parameters")     # --PWadjustment=PWadjustment
	keyword_dict["--mtf"]                         = SXkeyword_map(2, "parameters")     # --mtf=MTF_FILE_NAME
	keyword_dict["--chunk"]                       = SXkeyword_map(2, "parameters")     # --chunk0=CHUNK0_FILE_NAME, --chunk1=CHUNK1_FILE_NAME
	keyword_dict["--list"]                        = SXkeyword_map(2, "parameters")     # --list
	keyword_dict["inputfile"]                     = SXkeyword_map(2, "any_file")       # inputfile
	keyword_dict["unblur"]                        = SXkeyword_map(2, "exe")            # unblur
	keyword_dict["input_pdb"]                     = SXkeyword_map(2, "pdb")            # input_pdb
	keyword_dict["input_mrc_micrograph"]          = SXkeyword_map(2, "mrc")            # input_mrc_micrograph
	keyword_dict["input_bdb_stack_file"]          = SXkeyword_map(2, "bdb")            # input_bdb_stack_file
	keyword_dict["input_shift_list_file"]         = SXkeyword_map(2, "txt")            # input_shift_list_file
	keyword_dict["cter_ctf_file"]                 = SXkeyword_map(2, "txt")            # cter_ctf_file
	keyword_dict["input_data_list"]               = SXkeyword_map(2, "any_file_list")  # input_data_list
	keyword_dict["--function"]                    = SXkeyword_map(2, "function")       # --function=user_function
	keyword_dict["--previous_run"]                = SXkeyword_map(2, "directory")      # --previous_run1=run1_directory, --previous_run2=run2_directory
	keyword_dict["input_bdb_stack_pattern"]       = SXkeyword_map(2, "any_directory")  # input_bdb_stack_pattern

	keyword_dict["--apix"]                        = SXkeyword_map(2, "apix")           # --apix=pixel_size, --apix
	keyword_dict["--pixel_size"]                  = SXkeyword_map(2, "apix")           # --pixel_size=PIXEL_SIZE
	keyword_dict["--wn"]                          = SXkeyword_map(2, "ctfwin")         # --wn
	keyword_dict["--box"]                         = SXkeyword_map(2, "box")            # --box=box_size, --box_size=box_size
	keyword_dict["--radius"]                      = SXkeyword_map(2, "radius")         # --radius=particle_radius, --radius=outer_radius, --radius=outer_radius, --radius=particle_radius, --radius=outer_radius, --radius=outer_radius
	keyword_dict["--sym"]                         = SXkeyword_map(2, "sym")            # --sym=c1, --sym=c1, --sym=c1, --sym=symmetry, --sym=c1, --sym=c4

	# NOTE: 2016/02/23 Toshio Moriya
	# Below might be useful to include
	# reference power spectrum? --pwreference of viper, --pwreference of rviper, --PWadjustment of sort3d, --PWadjustment of rsort3d
	#
	# Below must be exceptional cases
	# --wn of locres, sort3d, & rsort3d; same as ctfwin?
	# --radius of locres & filterlocal; same as radius?
	#
	
	return keyword_dict

# ----------------------------------------------------------------------------------------
def handle_exceptional_keywords(sxcmd):
	# DESIGN_NOTE: 2016/02/05 Toshio Moriya
	# Handle exceptional cases due to the limitation of software design
	# In future, we should remove these exception handling by reviewing the design
	if sxcmd.name == "sxisac_post_processing":
		assert(sxcmd.token_dict["stack_file"].key_base == "stack_file")
		assert(sxcmd.token_dict["stack_file"].type == "image")
		sxcmd.token_dict["stack_file"].type = "bdb"
	elif sxcmd.name == "sxfilterlocal":
		assert(sxcmd.token_dict["locres_volume"].key_base == "locres_volume")
		assert(sxcmd.token_dict["locres_volume"].type == "output")
		sxcmd.token_dict["locres_volume"].type = "image"
	elif sxcmd.name in ["sxlocres",  "sxsort3d", "sxrsort3d"]:
		assert(sxcmd.token_dict["wn"].key_base == "wn")
		assert(sxcmd.token_dict["wn"].type == "ctfwin")
		sxcmd.token_dict["wn"].type = "int"
	elif sxcmd.name in ["sxrviper"]:
		assert(sxcmd.token_dict["stack"].key_base == "stack")
		assert(sxcmd.token_dict["stack"].type == "image")
		sxcmd.token_dict["stack"].type = "hdf"
		# Typically, this is target particle radius used by ISAC.
		assert(sxcmd.token_dict["radius"].key_base == "radius")
		assert(sxcmd.token_dict["radius"].type == "radius")
		sxcmd.token_dict["radius"].type = "int"
	elif sxcmd.name in ["sxviper"]:
		# Typically, this is target particle radius used by ISAC.
		assert(sxcmd.token_dict["radius"].key_base == "radius")
		assert(sxcmd.token_dict["radius"].type == "radius")
		sxcmd.token_dict["radius"].type = "int"

# ----------------------------------------------------------------------------------------
def remove_MoinMoinWiki_makeup(target_text):
	# makeup for link
	# [[URL|DISPLAY_TEXT]]
	makeup_begin = "[["
	makeup_end = "]]"
	makeup_separator = "|"

	item_head = target_text.find(makeup_begin)
	while item_head != -1:
		# Found a start of wiki makeup
		item_tail = target_text.find(makeup_end)
		if item_tail == -1:
			ERROR("Wiki Format Warning: The string \"%s\" contains \"%s\" but not \"%s\". Removing \"%s\", but please check the format in Wiki document." % (target_text, makeup_begin, makeup_end, makeup_begin), "%s in %s" % (__name__, os.path.basename(__file__)))
			target_text = target_text.replace(makeup_begin, "", 1)
		else: # assert (item_tail > -1)
			makeup_token = target_text[item_head:item_tail+len(makeup_end)]
			display_item = makeup_token.replace(makeup_begin, "")
			display_item = display_item.replace(makeup_end, "")
			if display_item.find(makeup_separator) != -1:
				item_tokens = display_item.split(makeup_separator)
				assert (len(item_tokens) == 2)
				display_item = item_tokens[1] # 2nd one should be display text
			print "### Found a wiki makeup token \"%s\". Changed to \"%s\"" % (makeup_token, display_item)
			target_text = target_text.replace(makeup_token, display_item, 1)

		# Try to find the next
		item_head = target_text.find(makeup_begin)

	return target_text


# ----------------------------------------------------------------------------------------
def construct_token_list_from_MoinMoinWiki(sxcmd_config):

	print "Start parsing MoinMoinWiki document (%s as %s %s command) " % (sxcmd_config.wiki, sxcmd_config.category, sxcmd_config.role)

	if sxcmd_config.format != "MoinMoinWiki": ERROR("Logical Error: Incorrect Wiki format %s! Check the sxcmd_config setting in this script." % (sxcmd_config.format), "%s in %s" % (__name__, os.path.basename(__file__)))
					
	# Allocate memory for new SXcmd instance
	sxcmd = SXcmd(sxcmd_config.category, sxcmd_config.role, sxcmd_config.is_submittable)
	
	# Define dictionary of keywords:
	# The dictionary maps command token to special data types
	keyword_dict = construct_keyword_dict()

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
					if item_tail == -1: ERROR("Wiki Format Error: '= Name =' section should contain only one valid line, and the line should starts from 'sx* - ' or 'e2* - ': %s" % line_wiki, "%s in %s" % (__name__, os.path.basename(__file__)))
					sxcmd.name = line_buffer[0:item_tail].strip()
					line_buffer = line_buffer[item_tail + len(target_operator):].strip() # Get the rest of line
					# Extract the label of this sxscript
					target_operator = ":"
					item_tail = line_buffer.find(target_operator)
					if item_tail == -1: ERROR("Wiki Format Error: '= Name =' section should contain a label ended with ':' after 'sx* - ' or 'e2* - ': %s" % line_wiki, "%s in %s" % (__name__, os.path.basename(__file__)))
					sxcmd.label = line_buffer[0:item_tail].strip()
					# Extract the short info about this sxscript (can be empty)
					sxcmd.short_info = remove_MoinMoinWiki_makeup(line_buffer[item_tail + len(target_operator):].strip()) # Get the rest of line
				elif current_section == section_usage:
					# Extract 'usage in command line' to identify each command token is either an argument (no-prefix) or option ('--' prefix)
					# This information is also used to check consistency between 'usage in command line' and list in '== Input ==' and '== Output ==' sections
					if line_wiki[0:len("sx")] == "sx" or line_wiki[0:len("e2")] == "e2":
						usage_token_list = line_wiki.split()
						if usage_token_list[0] != sxcmd.name + ".py": ERROR("Wiki Format Error: First token should be script name with .py (sx*.py or e2*.py)", "%s in %s" % (__name__, os.path.basename(__file__)))
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
						token.help = remove_MoinMoinWiki_makeup(line_buffer[0:item_tail])
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
	
	handle_exceptional_keywords(sxcmd)

	print "Succeed to parse MoinMoinWiki document (%s as %s %s command)" % (sxcmd_config.wiki, sxcmd_config.category, sxcmd_config.role)

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


# ----------------------------------------------------------------------------------------
def remove_DokuWiki_makeup(target_text):
	# makeup for link
	# [[URL|DISPLAY_TEXT]]
	makeup_begin = "[["
	makeup_end = "]]"
	makeup_separator = "|"

	item_head = target_text.find(makeup_begin)
	while item_head != -1:
		# Found a start of wiki makeup
		item_tail = target_text.find(makeup_end)
		if item_tail == -1:
			ERROR("Wiki Format Warning: The string \"%s\" contains \"%s\" but not \"%s\". Removing \"%s\", but please check the format in Wiki document." % (target_text, makeup_begin, makeup_end, makeup_begin), "%s in %s" % (__name__, os.path.basename(__file__)))
			target_text = target_text.replace(makeup_begin, "", 1)
		else: # assert (item_tail > -1)
			makeup_token = target_text[item_head:item_tail+len(makeup_end)]
			display_item = makeup_token.replace(makeup_begin, "")
			display_item = display_item.replace(makeup_end, "")
			if display_item.find(makeup_separator) != -1:
				item_tokens = display_item.split(makeup_separator)
				assert (len(item_tokens) == 2)
				display_item = item_tokens[1] # 2nd one should be display text
			print "### Found a wiki makeup token \"%s\". Changed to \"%s\"" % (makeup_token, display_item)
			target_text = target_text.replace(makeup_token, display_item, 1)

		# Try to find the next
		item_head = target_text.find(makeup_begin)

	return target_text


# ----------------------------------------------------------------------------------------
def construct_token_list_from_DokuWiki(sxcmd_config):

	print "Start parsing DokuWiki document (%s as %s %s command) " % (sxcmd_config.wiki, sxcmd_config.category, sxcmd_config.role)

	if sxcmd_config.format != "DokuWiki": ERROR("Logical Error: Incorrect Wiki format %s! Check the sxcmd_config setting in this script." % (sxcmd_config.format), "%s in %s" % (__name__, os.path.basename(__file__)))
	
	# Allocate memory for new SXcmd instance
	sxcmd = SXcmd(sxcmd_config.category, sxcmd_config.role, sxcmd_config.is_submittable)

	# Define dictionary of keywords:
	# The dictionary maps command token to special data types
	keyword_dict = construct_keyword_dict()

	# Define list of target sections for GUI and set current
	section_lists = []
	section_lists.append("====== Usage ======"); section_usage = len(section_lists) - 1;
	section_lists.append("====== Typical usage ======"); section_typical = len(section_lists) - 1;
	section_lists.append("====== Input ======"); section_input = len(section_lists) - 1;
	section_lists.append("====== Output ======"); section_output = len(section_lists) - 1;
	current_section = section_usage

	# Define list of subsections of input section and set current
	group_lists = []
	group_lists.append("=== Main Options ==="); group_main = len(group_lists) - 1;
	group_lists.append("=== Advanced Options ==="); group_advanced = len(group_lists) - 1;
	current_group_name = group_lists[group_main].replace("===", "").replace("Options", "").strip().lower() 
	
	# Define States and set current
	state_searching_header  = 0
	state_processing_header  = 1
	state_searching  = 2
	state_processing = 3
	state_done = 4
	current_state = state_searching_header # Assuming the first line starts from ====== COMMAND_NAME =======

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
		# print "MRK_DEBUG: line_wiki = %s" % (line_wiki)

		if current_state == state_searching_header:
			# Extract the name of sxscript
			target_operator = "======"
			if line_wiki[:len(target_operator)] != target_operator or line_wiki[-len(target_operator):] != target_operator:
				ERROR("Wiki Format Error: The Wiki document must start from header section title defined by ====== COMMAND_NAME =======.", "%s in %s" % (__name__, os.path.basename(__file__)))
			sxcmd.name = line_wiki.replace(target_operator, "").strip()
			current_state = state_processing_header
		elif current_state == state_processing_header:
			# Extract the label of this sxscript
			line_buffer = line_wiki
			target_operator = ":"
			item_tail = line_buffer.find(target_operator)
			if item_tail == -1: ERROR("Wiki Format Error: Header section body should contain a label ended with ':': %s" % line_wiki, "%s in %s" % (__name__, os.path.basename(__file__)))
			sxcmd.label = line_buffer[0:item_tail].strip()
			# Extract the short info about this sxscript (can be empty)
			sxcmd.short_info = remove_DokuWiki_makeup(line_buffer[item_tail + len(target_operator):].strip()) # Get the rest of line
			current_state = state_searching
		elif current_state == state_searching:
			if line_wiki.find(section_lists[current_section]) != -1:
				# Found the current target section
				current_state = state_processing
			# else: just ignore this line
		else:
			if current_state != state_processing: ERROR("Logical Error: This condition should not happen! State setting must be incorrect.", "%s in %s" % (__name__, os.path.basename(__file__)))
			target_operator = "======"
			if line_wiki[:len(target_operator)] == target_operator: # Assuming the section always starts with "======"
				# Reached the next section (might be not target)
				current_section += 1 # Update current target section
				current_group_name = group_lists[group_main].replace("===", "").replace("Options", "").strip().lower()    # reset group (subsection) for '====== Input ======' and '====== Output ======' sections
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
				if current_section == section_usage:
					# Extract 'usage in command line' to identify each command token is either an argument (no-prefix) or option ('--' prefix)
					# This information is also used to check consistency between 'usage in command line' and list in '== Input ==' and '== Output ==' sections
					if line_wiki[0:len("sx")] == "sx" or line_wiki[0:len("e2")] == "e2":
						usage_token_list = line_wiki.split()
						if usage_token_list[0] != sxcmd.name + ".py": ERROR("Wiki Format Error: First token should be script name with .py (sx*.py or e2*.py)", "%s in %s" % (__name__, os.path.basename(__file__)))
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
					if line_wiki.find(group_lists[group_advanced]) > -1:
						# Reached the option subsection (argument subsection is done)
						current_group_name = group_lists[group_advanced].replace("===", "").replace("Options", "").strip().lower() 
					elif line_wiki[0] == ";":
						# Option entry must start with ";"
						line_buffer = line_wiki[1:]  # remove ";"
						# Extract key base name of command token
						target_operator = ":"
						item_tail = line_buffer.find(target_operator)
						if item_tail == -1:
							# ERROR("Warning: This line (%s) is missing key base name (maybe comment line?). Ignoring this line...'."  % (line_wiki),"%s in %s" % (__name__, os.path.basename(__file__)), action = 0)
							continue
						key_base = line_buffer[0:item_tail].strip()
						if key_base == "MPI":
							# ERROR("Warning: This line (%s) contains MPI flag. The flag will be removed in near future, so ignoring this line...'."  % (line_wiki), "%s in %s" % (__name__, os.path.basename(__file__)), action = 0)
							if sxcmd.mpi_support == False or sxcmd.mpi_add_flag == False: ERROR("Logical Error: Since MPI flag is found, the command should support MPI.", "%s in %s" % (__name__, os.path.basename(__file__)))
							continue
						line_buffer = line_buffer[item_tail + len(target_operator):].strip() # Get the rest of line
						# check consistency between 'usage in command line' and this
						if key_base not in sxcmd.token_dict.keys(): ERROR("Wiki Format Error: Key base (%s) is missing from 'usage in command line' in '====== Usage ======'." % key_base, "%s in %s" % (__name__, os.path.basename(__file__)))
						# Get the reference to the command token object associated with this key base name
						token = sxcmd.token_dict[key_base]
						if token.key_base != key_base: ERROR("Logical Error: Registered command token with wrong key base name into the dictionary.", "%s in %s" % (__name__, os.path.basename(__file__)))
						token.is_in_io = True # Set flag to tell this command token is find in input or output section
						token.group = current_group_name # Set group of command token according to the current subsection
						# Extract label of command token
						target_operator = ":"
						item_tail = line_buffer.find(target_operator)
						if item_tail == -1: ERROR("Wiki Format Error: This line (%s) is missing label. Please check the format in Wiki document." % line_wiki, "%s in %s" % (__name__, os.path.basename(__file__)))
						token.label = line_buffer[0:item_tail].strip()
						line_buffer = line_buffer[item_tail + len(target_operator):].strip() # Get the rest of line
						# Extract help of command token before default value
						target_operator = "(default"
						item_tail = line_buffer.find(target_operator)
						if item_tail == -1: ERROR("Wiki Format Error: This line (%s) is missing default setting. Please check the format in Wiki document." % line_wiki, "%s in %s" % (__name__, os.path.basename(__file__)))
						token.help = remove_DokuWiki_makeup(line_buffer[0:item_tail])
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
					# else: 
						# This is not option entry. Ignore this line
				else:
					ERROR("Logical Error: This section is invalid. Did you assigne an invalid section?", "%s in %s" % (__name__, os.path.basename(__file__)))

	if current_state != state_done: ERROR("Wiki Format Error: parser could not extract all information necessary. Please check if the Wiki format has all required sections.", "%s in %s" % (__name__, os.path.basename(__file__)))

	# Make sure there are no extra arguments or options in 'usage in command line' of '====== Usage ======'
	for token in sxcmd.token_list:
		if token.is_in_io == False: ERROR("Wiki Format Error: An extra argument or option (%s) is found in 'usage in command line' of '====== Usage ======'." % token.key_base, "%s in %s" % (__name__, os.path.basename(__file__)))

	file_wiki.close()
	
	handle_exceptional_keywords(sxcmd)

	print "Succeed to parse MoinMoinWiki document (%s as %s %s command)" % (sxcmd_config.wiki, sxcmd_config.category, sxcmd_config.role)

	"""
	# For DEBUG
	if sxcmd.name == "sxmeridien":
		print "><><>< DEBUG OUTPUT ><><><"
		print ""
		print "------"
		print "GLOBAL"
		print "------"
		print "name            : %s" % sxcmd.name
		print "mode            : %s" % sxcmd.mode
		print "label           : %s" % sxcmd.label
		print "short_info      : %s" % sxcmd.short_info
		print "mpi_support     : %s" % sxcmd.mpi_support
		print "mpi_add_flag    : %s" % sxcmd.mpi_add_flag
		print "category        : %s" % sxcmd.category
		print "role            : %s" % sxcmd.role
		print "is_submittable  : %s" % sxcmd.is_submittable
		print "len(token_list) : %d" % len(sxcmd.token_list)
		print "len(token_dict) : %d" % len(sxcmd.token_dict)
		print ""
		print "--------------"
		print "cmd_token_list"
		print "--------------"
		for token in sxcmd.token_list:
			print "%s%s (group=%s, required=%s, default=%s, type=%s, restore=%s label=%s help=%s" % (token.key_prefix, token.key_base, token.group, token.is_required, token.default, token.type, token.restore, token.label, token.help)
		print ""
	"""
	
	return sxcmd


# ========================================================================================
def apply_exclude_list(sxcmd_config_exclude_list, sxcmd):
	assert(len(sxcmd_config_exclude_list) > 0)
	assert(len(sxcmd.token_list) == len(sxcmd.token_dict))

	for token_key_base in sxcmd_config_exclude_list:
		sxcmd.token_list.remove(sxcmd.token_dict[token_key_base])
		del sxcmd.token_dict[token_key_base]

# ----------------------------------------------------------------------------------------
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
	if mode_token_edit.key_base not in fullset_token_dict.keys(): ERROR("Logical Error: This condition should not happen! Subset command configuration must be incorrect. Key (%s) should not exists." % (mode_token_edit.key_base), "%s in %s" % (__name__, os.path.basename(__file__)))
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
			if token_edit.key_prefix != "": ERROR("Logical Error: This condition should not happen! Subset command configuration must be incorrect. Key (%s) should be argument." % (token_edit.key_base) , "%s in %s" % (__name__, os.path.basename(__file__)))
			token = token_edit
		else:
			# token key base is found in fullset. This must be an option.
			token = fullset_token_dict[token_edit.key_base]
			if token_edit.key_prefix != None:
				token.key_prefix = token_edit.key_prefix
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

# ----------------------------------------------------------------------------------------
def insert_sxcmd_category_list_to_file(sxcmd_category_list, output_file):
	sxcmd_category_variable_name = "sxcmd_category"
	for sxcmd_category in sxcmd_category_list:
		output_file.write("\t\t")
		output_file.write("%s = SXcmd_category()" % sxcmd_category_variable_name)
		output_file.write("; %s.name = \"%s\"" % (sxcmd_category_variable_name, sxcmd_category.name))
		output_file.write("; %s.label = \"%s\"" % (sxcmd_category_variable_name, sxcmd_category.label))
		output_file.write("; %s.short_info = \"%s\"" % (sxcmd_category_variable_name, sxcmd_category.short_info.replace("\"", "'")))
		output_file.write("\n")
		output_file.write("\t\t")
		output_file.write("%s_list.append(%s)" % (sxcmd_category_variable_name, sxcmd_category_variable_name))
		output_file.write("\n")

	output_file.write("\n")
	return

# ----------------------------------------------------------------------------------------
def insert_sxcmd_to_file(sxcmd, output_file, sxcmd_variable_name):
	output_file.write("\t\t")
	output_file.write("%s = SXcmd()" % sxcmd_variable_name)
	output_file.write("; %s.name = \"%s\"" % (sxcmd_variable_name, sxcmd.name))
	output_file.write("; %s.mode = \"%s\"" % (sxcmd_variable_name, sxcmd.mode))
	output_file.write("; %s.label = \"%s\"" % (sxcmd_variable_name, sxcmd.label))
	output_file.write("; %s.short_info = \"%s\"" % (sxcmd_variable_name, sxcmd.short_info.replace("\"", "'")))
	output_file.write("; %s.mpi_support = %s" % (sxcmd_variable_name, sxcmd.mpi_support))
	output_file.write("; %s.mpi_add_flag = %s" % (sxcmd_variable_name, sxcmd.mpi_add_flag))
	output_file.write("; %s.category = \"%s\"" % (sxcmd_variable_name, sxcmd.category))
	output_file.write("; %s.role = \"%s\"" % (sxcmd_variable_name, sxcmd.role))
	output_file.write("; %s.is_submittable = %s" % (sxcmd_variable_name, sxcmd.is_submittable))
	output_file.write("\n")

	for token in sxcmd.token_list:
		output_file.write("\t\t")
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

	output_file.write("\n")
	output_file.write("\t\t%s_list.append(%s)\n" % (sxcmd_variable_name, sxcmd_variable_name))
	output_file.write("\n")

	return

# ========================================================================================
def create_sxcmd_subconfig_window_makevstack():
	token_edit_list = []
	token_edit = SXcmd_token(); token_edit.initialize_edit("makevstack"); token_edit.is_required = True; token_edit.default = "none"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("input_bdb_stack_pattern"); token_edit.key_prefix = ""; token_edit.label = "Input BDB image stack pattern"; token_edit.help = "Specify file path pattern of stack subsets created in particle extraction using a wild card /'*/' (e.g. /'//sxwindow_output_dir//*/'). The stack subsets are located in the sxwindow output directory."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "any_directory"; token_edit_list.append(token_edit)
	
	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("Particle Stack", token_edit_list, sxsubcmd_mpi_support)

	return sxcmd_subconfig

def create_sxcmd_subconfig_isacselect():
	token_edit_list = []
	token_edit = SXcmd_token(); token_edit.initialize_edit("isacselect"); token_edit.is_required = True; token_edit.default = True; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("class_file_name"); token_edit.key_prefix = ""; token_edit.label = "ISAC class file name"; token_edit.help = "File name of the class averages. It is located in the ISAC output directory."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "image"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("output_list"); token_edit.key_prefix = ""; token_edit.label = "Output ISAC particle ID list"; token_edit.help = "Output text file containing retrieved member particle IDs of all ISAC classes."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "output"; token_edit_list.append(token_edit)
	
	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("Get ISAC Particles", token_edit_list, sxsubcmd_mpi_support)

	return sxcmd_subconfig

def create_sxcmd_subconfig_isac_makevstack():
	token_edit_list = []
	token_edit = SXcmd_token(); token_edit.initialize_edit("makevstack"); token_edit.is_required = True; token_edit.default = "none"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("input_bdb_stack_file"); token_edit.key_prefix = ""; token_edit.label = "Input BDB image stack"; token_edit.help = "File name of the class averages. It is located in the ISAC output directory."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "bdb"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("list"); token_edit_list.append(token_edit)
	
	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("Create Stack Subset", token_edit_list, sxsubcmd_mpi_support)

	return sxcmd_subconfig

def create_sxcmd_subconfig_changesize():
	token_edit_list = []
	token_edit = SXcmd_token(); token_edit.initialize_edit("changesize"); token_edit.is_required = True; token_edit.default = True; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("input_stack"); token_edit.key_prefix = ""; token_edit.label = "Input 2D/3D image stack"; token_edit.help = "Input 2D/3D image stack."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "image"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("output_stack"); token_edit.key_prefix = ""; token_edit.label = "Output 2D/3D image stack"; token_edit.help = "Resampled (decimated or interpolated up) 2D/3D image stack."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "output"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("ratio"); token_edit_list.append(token_edit)
	
	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("Resample VIPER Model", token_edit_list, sxsubcmd_mpi_support)

	return sxcmd_subconfig

def create_sxcmd_subconfig_clip():
	token_edit_list = []
	token_edit = SXcmd_token(); token_edit.initialize_edit("clip"); token_edit.label = "Pad/Clip volume to specified size [Pixels]"; token_edit.is_required = True; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("input_volume"); token_edit.key_prefix = ""; token_edit.label = "Input volume"; token_edit.help = "Input volume file name."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = "none"; token_edit.type = "image"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("output_file"); token_edit.key_prefix = ""; token_edit.label = "Output clipped/padded volume"; token_edit.help = "Output clipped/padded volume file name."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = "none"; token_edit.type = "output"; token_edit_list.append(token_edit)
	
	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("Pad/Clip VIPER Model", token_edit_list, sxsubcmd_mpi_support)

	return sxcmd_subconfig

# def create_sxcmd_subconfig_scale_clip():
# 	token_edit_list = []
# 	token_edit = SXcmd_token(); token_edit.initialize_edit("scale"); token_edit.label = "Resample ratio"; token_edit.help = "Rescale the volume by the specified ratio before padding/clipping."; token_edit.is_required = True; token_edit_list.append(token_edit)
# 	token_edit = SXcmd_token(); token_edit.initialize_edit("clip"); token_edit.label = "Pad/Clip volume to specified size [Pixels]"; token_edit.is_required = True; token_edit_list.append(token_edit)
# 	token_edit = SXcmd_token(); token_edit.initialize_edit("input_volume"); token_edit.key_prefix = ""; token_edit.label = "Input volume"; token_edit.help = "Input volume file name."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = "none"; token_edit.type = "image"; token_edit_list.append(token_edit)
# 	token_edit = SXcmd_token(); token_edit.initialize_edit("output_file"); token_edit.key_prefix = ""; token_edit.label = "Output resampled volume"; token_edit.help = "Output resampled volume file name."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = "none"; token_edit.type = "output"; token_edit_list.append(token_edit)
# 	
# 	sxsubcmd_mpi_support = False
# 	sxcmd_subconfig = SXsubcmd_config("Resample", token_edit_list, sxsubcmd_mpi_support)
# 
# 	return sxcmd_subconfig

def create_sxcmd_subconfig_adaptive_mask3d():
	token_edit_list = []
	token_edit = SXcmd_token(); token_edit.initialize_edit("adaptive_mask"); token_edit.is_required = True; token_edit.default = True; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("input_volume"); token_edit.key_prefix = ""; token_edit.label = "Input volume"; token_edit.help = "Input reference volume"; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "image"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("output_mask3D"); token_edit.key_prefix = ""; token_edit.label = "Output mask"; token_edit.help = "Output 3D mask"; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "output"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("nsigma"); token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("threshold"); token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("ndilation"); token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("kernel_size"); token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("gauss_standard_dev"); token_edit_list.append(token_edit)
	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("Adaptive 3D Mask", token_edit_list, sxsubcmd_mpi_support)

	return sxcmd_subconfig

def create_sxcmd_subconfig_binary_mask3d():
	token_edit_list = []
	token_edit = SXcmd_token(); token_edit.initialize_edit("binary_mask"); token_edit.is_required = True; token_edit.default = True; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("input_volume"); token_edit.key_prefix = ""; token_edit.label = "Input volume"; token_edit.help = "Input reference volume"; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "image"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("output_mask3D"); token_edit.key_prefix = ""; token_edit.label = "Output mask"; token_edit.help = "Output 3D mask"; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "output"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("bin_threshold"); token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("ne"); token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("nd"); token_edit_list.append(token_edit)
	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("Binary 3D Mask", token_edit_list, sxsubcmd_mpi_support)

	return sxcmd_subconfig

def create_sxcmd_subconfig_refine3d_postprocess():
	token_edit_list = []
	token_edit = SXcmd_token(); token_edit.initialize_edit("postprocess"); token_edit.is_required = True; token_edit.default = True; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("firstvolume"); token_edit.key_prefix = ""; token_edit.label = "First unfiltered half-volume "; token_edit.help = "Generated by sxmeridien"; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "image"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("secondvolume"); token_edit.key_prefix = ""; token_edit.label = "Second unfiltered half-volume "; token_edit.help = "Generated by sxmeridien"; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "image"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("mtf"); token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("fsc_adj"); token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("B_enhance"); token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("low_pass_filter"); token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("aa"); token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("mask"); token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("output"); token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("pixel_size"); token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("B_start"); token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("B_stop"); token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("do_adaptive_mask"); token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("mask_threshold"); token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("consine_edge"); token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("dilation"); token_edit_list.append(token_edit)
	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("3D Refinement Postprocess", token_edit_list, sxsubcmd_mpi_support)

	return sxcmd_subconfig

def create_sxcmd_subconfig_variability_preprocess():
	token_edit_list = []
	token_edit = SXcmd_token(); token_edit.initialize_edit("symmetrize"); token_edit.is_required = True; token_edit.default = True; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("prj_stack"); token_edit.key_prefix = ""; token_edit.label = "Input image stack"; token_edit.help = "The images must containt the 3D orientation parameters in the header and optionally CTF information. The output image stack is bdb:sdata. Please use it as an input image stack of sx3dvariability."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "image"; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("sym"); token_edit_list.append(token_edit)
	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("3D Variability Preprocess", token_edit_list, sxsubcmd_mpi_support)

	return sxcmd_subconfig

def create_sxcmd_subconfig_refine3d_angular_distribution():
	token_edit_list = []
	token_edit = SXcmd_token(); token_edit.initialize_edit("angular_distribution"); token_edit.is_required = True; token_edit.default = True; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("inputfile"); token_edit.label = "Alignment Parameter file"; token_edit.help = "Alignment Parameter file created by a previous 3D reconstruction step (e.g. sxmeridien.py)"; token_edit.is_required = True; token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("pixel_size"); token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("round_digit"); token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("box_size"); token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("particle_radius"); token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("cylinder_width"); token_edit_list.append(token_edit)
	token_edit = SXcmd_token(); token_edit.initialize_edit("cylinder_length"); token_edit_list.append(token_edit)
	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("Angular Distribution", token_edit_list, sxsubcmd_mpi_support)

	return sxcmd_subconfig

def create_exclude_list_boxer():
	exclude_list = []

	exclude_list.append("write_dbbox")
	exclude_list.append("write_ptcls")
	exclude_list.append("force")
	exclude_list.append("format")
	exclude_list.append("suffix")
	exclude_list.append("dbls")
	exclude_list.append("autoboxer")
	exclude_list.append("ppid")
	exclude_list.append("gui")
	exclude_list.append("do_ctf")
	exclude_list.append("cter")
	exclude_list.append("indir")
	exclude_list.append("nameroot")
	exclude_list.append("micsuffix")
	exclude_list.append("wn")
	exclude_list.append("Cs")
	exclude_list.append("voltage")
	exclude_list.append("ac")
	exclude_list.append("kboot")
	exclude_list.append("debug")
	exclude_list.append("apix")

	return exclude_list

def create_exclude_list_display():
	exclude_list = []

	exclude_list.append("classmx")
	exclude_list.append("classes")
	exclude_list.append("pdb")
	exclude_list.append("plot")
	exclude_list.append("plot3")
	exclude_list.append("newwidget")
	exclude_list.append("ppid")

	return exclude_list

# ========================================================================================
def main():
	# --------------------------------------------------------------------------------
	# Define command categories used in GUI
	# --------------------------------------------------------------------------------
	sxcmd_category_list = []
	sxcmd_category_list.append(SXcmd_category("sxc_movie", "Movie Micrograph", "movie frame alignemnt, and drift assessment"))
	sxcmd_category_list.append(SXcmd_category("sxc_cter", "CTF", "ctf estinatim, and ctf assessment"))
	sxcmd_category_list.append(SXcmd_category("sxc_window", "Particle Stack", "particle picking, and particle windowing"))
	sxcmd_category_list.append(SXcmd_category("sxc_isac", "2D Clustering", "2d clustering with isac, and post-processing"))
	sxcmd_category_list.append(SXcmd_category("sxc_viper", "Initial 3D Modeling", "initial 3d modeling with viper/rviper"))
	sxcmd_category_list.append(SXcmd_category("sxc_meridien", "3D Refinement", "3d refinement and post-processing"))
	sxcmd_category_list.append(SXcmd_category("sxc_sort3d", "3D Clustering", "3d variability, and 3d clustering protocol I & II"))
	sxcmd_category_list.append(SXcmd_category("sxc_localres", "Local Resolution", "local resolution, and local filter"))
	sxcmd_category_list.append(SXcmd_category("sxc_utilities", "Utilities", "miscellaneous utlitity commands"))

	# --------------------------------------------------------------------------------
	# Get all necessary informations from wiki documents of sx*.py scripts
	# and create gui generation parameter
	# --------------------------------------------------------------------------------
	sxcmd_config_list = []
	
	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_movie"

	sxcmd_role = "sxr_pipe"
	sxcmd_config_list.append(SXcmd_config("../doc/unblur.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/gui_unblur.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, is_submittable = False))

	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/e2display.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))

	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_cter"

	sxcmd_role = "sxr_pipe"
	sxcmd_config_list.append(SXcmd_config("../doc/cter.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list=["stack_mode"]))
	sxcmd_config_list.append(SXcmd_config("../doc/gui_cter.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, is_submittable = False))

	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/e2display.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))

	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_window"

	sxcmd_role = "sxr_pipe"
	sxcmd_config_list.append(SXcmd_config("../doc/e2boxer_old.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_boxer(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/window.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/e2bdb.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_window_makevstack()))

	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/e2display.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))

	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_isac"

	sxcmd_role = "sxr_pipe"
	sxcmd_config_list.append(SXcmd_config("../doc/isac.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_isacselect()))
	sxcmd_config_list.append(SXcmd_config("../doc/e2bdb.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_isac_makevstack()))
	sxcmd_config_list.append(SXcmd_config("../doc/isac_post_processing.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))

	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/e2display.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))

	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_viper"

	sxcmd_role = "sxr_pipe"
	sxcmd_config_list.append(SXcmd_config("../doc/rviper.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_changesize()))
	sxcmd_config_list.append(SXcmd_config("../doc/e2proc3d.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_clip()))
# 	sxcmd_config_list.append(SXcmd_config("../doc/e2proc3d.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_scale_clip()))

	sxcmd_role = "sxr_alt"
	sxcmd_config_list.append(SXcmd_config("../doc/viper.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/pdb2em.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))

	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/e2display.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_adaptive_mask3d()))
	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_refine3d_angular_distribution()))

	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_meridien"

	sxcmd_role = "sxr_pipe"
	# sxcmd_config_list.append(SXcmd_config("../doc/meridien.doku.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	# sxcmd_config_list.append(SXcmd_config("../doc/meridien.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
	# sxcmd_config_list.append(SXcmd_config("../doc/meridien-08-08-2016.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/meridien-09-09-2016.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_refine3d_postprocess()))

	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/e2display.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_adaptive_mask3d()))
	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_refine3d_angular_distribution()))

	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_sort3d"

	sxcmd_role = "sxr_pipe"
	sxcmd_config_list.append(SXcmd_config("../doc/3dvariability.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_variability_preprocess()))
	sxcmd_config_list.append(SXcmd_config("../doc/3dvariability.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list=["symmetrize"]))
	sxcmd_config_list.append(SXcmd_config("../doc/sort3d.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/rsort3d.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))

	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/e2display.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_adaptive_mask3d()))
	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_binary_mask3d()))
	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_refine3d_angular_distribution()))

	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_localres"

	sxcmd_role = "sxr_pipe"
	sxcmd_config_list.append(SXcmd_config("../doc/locres.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/filterlocal.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))

	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/e2display.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_adaptive_mask3d()))
	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_refine3d_angular_distribution()))

	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_utilities"

	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/e2display.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/pdb2em.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_adaptive_mask3d()))
	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_binary_mask3d()))
	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_refine3d_angular_distribution()))
	# sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))

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
#	sxcmd_config_list.append(SXcmd_config("../doc/cter.txt", "MoinMoinWiki", "sxr_util", subconfig = sxcmd_subconfig))

	# --------------------------------------------------------------------------------
	# Check consistency between sxcmd_category_list and sxcmd_config_list
	# --------------------------------------------------------------------------------
	sxcmd_category_names = []
	for sxcmd_category in sxcmd_category_list:
		sxcmd_category_names.append(sxcmd_category.name)

	for sxcmd_config in sxcmd_config_list:
		if not sxcmd_config.category in sxcmd_category_names:
			ERROR("Logical Error: sxcmd_config for %s is using invalid category %s." % (sxcmd_config.wiki, sxcmd_config.category), "%s in %s" % (__name__, os.path.basename(__file__)))

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
				# Insert Command Category
				insert_sxcmd_category_list_to_file(sxcmd_category_list, output_file)

				# Insert Command List
				sxcmd_variable_name = "sxcmd"
				for sxcmd_config in sxcmd_config_list:
					# Construct sxscript object associated with this wiki document
					if sxcmd_config.format == "MoinMoinWiki":
						sxcmd = construct_token_list_from_MoinMoinWiki(sxcmd_config)
					elif sxcmd_config.format == "DokuWiki":
						sxcmd = construct_token_list_from_DokuWiki(sxcmd_config)
					else:
						ERROR("Logical Error: Invalid Wiki format %s! Check the sxcmd_config setting in this script." % (sxcmd_config.format), "%s in %s" % (__name__, os.path.basename(__file__)))
					
					if sxcmd_config.subconfig != None:
						apply_sxsubcmd_config(sxcmd_config.subconfig, sxcmd)
					if len(sxcmd_config.exclude_list) > 0:
						apply_exclude_list(sxcmd_config.exclude_list, sxcmd)
					insert_sxcmd_to_file(sxcmd, output_file, sxcmd_variable_name)
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
