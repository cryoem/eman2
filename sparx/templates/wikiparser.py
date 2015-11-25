#!/usr/bin/env python
# 
#
# Author: Toshio Moriya 12/11/2015 (toshio.moriya@mpi-dortmund.mpg.de)
# ========================================================================================

# from EMAN2 import *
# from sparx import *
import os
import sys
from global_def import ERROR

from sxgui_template import SXcmd_token, SXcmd

# ========================================================================================
def construct_token_list_from_wiki(wiki_file_path):
	# Private helper class used only in this function
	class SXkeyword_map:
		def __init__(self, priority, token_type):
			if priority >= 100: ERROR("Priority should be lower than 100", "%s in %s" % (__name__, os.path.basename(__file__)))
			# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
			# class variables
			self.priority = priority      # Priority of this keyword. Highest priority is 0. The type of higher priority will be used to avoid the conflict among keywords
			self.token_type = token_type  # Token value type 
			# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	
	print "Start parsing Wiki document (%s)" % wiki_file_path

	# Allocate memory for new SXcmd instance
	sxcmd = SXcmd()
	
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
	keyword_dict["importctf"]                   = SXkeyword_map(2, "parameters") # --importctf=ctf_file
	keyword_dict["pwreference"]                 = SXkeyword_map(2, "parameters") # --pwreference=pwreference_file
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
	state_done = 2
	current_state = state_searching
			
	# NOTE: 2015/11/11 Toshio Moriya
	# This should be exception. Need to decide if this should be skipped or exit system.
	if os.path.exists(wiki_file_path) == False: ERROR("Rutime Error: Wiki document is not found.", "%s in %s" % (__name__, os.path.basename(__file__)))
	
	file_wiki = open(wiki_file_path,'r')
	
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
					# Extract the name of sxscript
					target_operator = "-"
					item_tail = line_wiki.find(target_operator)
					if item_tail == -1: ERROR("Wiki Format Error: '= Name =' section should contain only one valid line, and the line should starts from 'sx* - ': %s" % line_wiki, "%s in %s" % (__name__, os.path.basename(__file__)))
					sxcmd.name = line_wiki[0:item_tail].strip()
					# Extract the short info about this sxscript (can be empty)
					sxcmd.short_info = line_wiki[item_tail + len(target_operator):].strip()							
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
								if token.key_base.find(keyword) != -1:
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
					if sxcmd.mpi_support == False and line_wiki.find(target_operator) > 1:
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
						# Ignore the rest of line ...
				else:
					ERROR("Logical Error: This section is invalid. Did you assigne an invalid section?", "%s in %s" % (__name__, os.path.basename(__file__)))
				
	if current_state != state_done: ERROR("Wiki Format Error: parser could not extract all information necessary. Please check if the Wiki format has all required sections.", "%s in %s" % (__name__, os.path.basename(__file__)))

	# Make sure there are no extra arguments or options in 'usage in command line' of '= Usage ='
	for token in sxcmd.token_list:
		if token.is_in_io == False: ERROR("Wiki Format Error: An extra argument or option (%s) is found in 'usage in command line' of '= Usage ='." % token.key_base, "%s in %s" % (__name__, os.path.basename(__file__)))
			
	file_wiki.close()
	
	print "Succeed to parse Wiki document (%s)" % wiki_file_path
	
	"""
	# For DEBUG
	if sxcmd.name == "sxwindow": 
		print "><><>< DEBUG OUTPUT ><><><"
		print ""
		print "------"
		print "GLOBAL"
		print "------"
		print "name            : %s" % sxcmd.name 
		print "short_info      : %s" % sxcmd.short_info 
		print "mpi_support     : %s" % sxcmd.mpi_support 
		print "mpi_add_flag    : %s" % sxcmd.mpi_add_flag 
		print "len(token_list) : %d" % len(sxcmd.token_list)
		print "len(token_dict) : %d" % len(sxcmd.token_dict)
		print ""
		print "--------------"
		print "cmd_token_list"
		print "--------------"
		for token in sxcmd.token_list:
			print "%s%s (group=%s, required=%s, default=%s, type=%s) <%s>" % (token.key_prefix, token.key_base, token.group, token.is_required, token.default, token.type, token.label),token.help
		print ""
	"""
	
	return sxcmd

def insert_sxcmd_to_file(sxcmd, output_file, sxcmd_variable_name):
	output_file.write("\t")
	output_file.write("%s = SXcmd()" % sxcmd_variable_name)
	output_file.write("; %s.name = \"%s\"" % (sxcmd_variable_name, sxcmd.name)) 
	output_file.write("; %s.short_info = \"%s\"" % (sxcmd_variable_name, sxcmd.short_info.replace("\"", "'")))
	output_file.write("; %s.mpi_support = %s" % (sxcmd_variable_name, sxcmd.mpi_support))
	output_file.write("; %s.mpi_add_flag = %s" % (sxcmd_variable_name, sxcmd.mpi_add_flag))
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
		if token.is_required:
			output_file.write("; token.default = \"\"")
		elif token.type == "bool":
			output_file.write("; token.default = %s" % token.default)
		else:
			output_file.write("; token.default = \"%s\"" % token.default)
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
	wiki_file_path_list = []
	wiki_file_path_list.append("../doc/cter.txt")
	wiki_file_path_list.append("../doc/window.txt")
	wiki_file_path_list.append("../doc/isac.txt")
	wiki_file_path_list.append("../doc/viper.txt")
	wiki_file_path_list.append("../doc/rviper.txt")
	wiki_file_path_list.append("../doc/meridien.txt")
	wiki_file_path_list.append("../doc/3dvariability.txt")
	wiki_file_path_list.append("../doc/locres.txt")
	wiki_file_path_list.append("../doc/filterlocal.txt")
	wiki_file_path_list.append("../doc/sort3d.txt")
	
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
				for wiki_file_path in wiki_file_path_list:
					# Construct sxscript object associated with this wiki document
					sxcmd = construct_token_list_from_wiki(wiki_file_path)
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
