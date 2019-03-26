#!/usr/bin/env python
from __future__ import print_function
#
#
# Author: Toshio Moriya 12/11/2015 (toshio.moriya@mpi-dortmund.mpg.de)
# ========================================================================================

import copy
import sp_global_def
import os
import sxgui_template
from builtins import object

# ========================================================================================
class SXsubcmd_config(object):
	def __init__(self, label = "", short_info = None, token_edit_list = [], mpi_support = None, is_modeless = False, subset_config=""):
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.label = label                              # User friendly name of command subset
		self.short_info = short_info                    # Short description of command subset
		self.token_edit_list = token_edit_list          # To edit some attributes of tokens (label, help, group, is_required, default). If the original value should be kept, set to None. First entry should be mode token.
		self.is_modeless = is_modeless                  # For special command subset, set this to True to suppress setting of SXcmd.mode.
		self.subset_config = subset_config              # Unique name to differentiate subset configuration of this command. For example, to name a command argument mode, which dependes on the number of input arguments. If not necessary, use empty string 
		self.mpi_support = mpi_support                  # Flag to indicate if this command suppors MPI version
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

# ========================================================================================
class SXcmd_config(object):
	def __init__(self, wiki, format, category, role, is_submittable = True, exclude_list = [], subconfig = None):
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.wiki = wiki                      # Wiki document file path
		self.format = format                  # Wiki document format: MoinMoinWiki or DokuWiki
		self.category = category              # Category of this command: sxc_movie, sxc_cter, sxc_window, sxc_isac, sxc_viper, sxc_meridien, sxc_sort3d, sxc_localres, sxc_utilities
		self.role = role                      # Role of this command; sxr_pipe (pipeline), sxr_alt (alternative) sxr_util (utility)
		self.is_submittable = is_submittable  # External GUI Application (e.g. sxgui_cter.py) should not be submitted to job queue. If it is true, this command will be registered to child_application_list of SXCmdWidget
		self.exclude_list = exclude_list      # token key base list to be excluded
		self.subconfig = subconfig            # Subset configuration of this command (e.g. sxprocess and sxlocres). None includes all command tokens, and do not make subcmd
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><


# ========================================================================================
# Helper class used only in construct_token_list_from_*() functions
class SXkeyword_map(object):
	def __init__(self, priority, token_type):
		if priority >= 100: sp_global_def.ERROR("Priority should be lower than 100", "%s in %s" % (__name__, os.path.basename(__file__)))
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
	# - output              : Line edit box for formatted string type, and output info button.
	#                         GUI also checks the existence of output directory/file before execution of the sx*.py
	#                         GUI abort the execution if the directory/file exists already
	# - output_continue     : Line edit box for formatted string type, open directory button, and output info button.
	#                         GUI also checks the existence of output directory/file before execution of the sx*.py
	#                         GUI displays a YES/NO dialog to run continuous mode if the directory/file exists already
	# - output_bdb2d_stack  : Line edit box for formatted string type, and output info button. In addition, open file buttons for .bdb (bdb2d_stack)
	#                         GUI also checks the existence of output directory/file before execution of the sx*.py
	#                         GUI abort the execution if the directory/file exists already
	# 
	# - displayable_list    : Line edit box for formatted string type, and open file button  for all supported data formats: 2D image, 2D stack, 3D volume, and 3D stack as well as 1D graph (displayable_list with bdb path parser)
	#                         The string with space is interpreted as a list of any data file names upon command generation. (i.e. does not enclose the string with single quotes)
	# - data2d3d_both       : Line edit box for formatted string type, and open file button  for all supported data formats: 2D image, 2D stack, 3D volume, and 3D stack (data2d3d_both with bdb path parser)
	# - mic_both            : Line edit box for formatted string type, and open file buttons for .mrc/.mrcs (mrc2d_mic_both) and all supported 2D image and 2D stack formats (mic_both with bdb path parser)
	# - mic_one             : Line edit box for formatted string type, and open file buttons for .mrc/.mrcs (mrc2d_mic_one) and all supported 2D image formats (mic_one with bdb path parser)
	# - mic_one_list        : Line edit box for formatted string type, and open file buttons for .mrc/.mrcs (mrc2d_mic_one_list) and all supported 2D image formats (mic_one_list with bdb path parser)
	#                         The string with space is interpreted as a list of any 2D image file names upon command generation. (i.e. does not enclose the string with single quotes)
	# - mic_stack           : Line edit box for formatted string type, and open file buttons for .mrcs/.mrc (mrc2d_mic_stack) and all supported 2D stack formats (mic_stack with bdb path parser)
	# - data2d_one          : Line edit box for formatted string type, and open file buttons for .hdf (hdf2d_one) and all supported 2D image formats (data2d_one with bdb path parser)
	# - bdb2d_stack         : Line edit box for formatted string type, and open file buttons for .bdb (bdb2d_stack)
	# - data2d_stack        : Line edit box for formatted string type, and open file buttons for .bdb (bdb2d_stack) and all supported 2D stack formats (data2d_stack with bdb path parser) NOTE: Toshio Moriya 2018/01/25 Currently, this case is not used. Instead, using bdb2d_stack
	# - data3d_one          : Line edit box for formatted string type, and open file buttons for .hdf (hdf3d_one) and all supported 3D volume formats (data3d_one with bdb path parser)
	# - data3d_stack        : Line edit box for formatted string type, and open file buttons for .hdf (hdf3d_stack) and all supported 3D stack formats (with bdb path parser) NOTE: Toshio Moriya 2018/01/25 Currently, this case is not used. Instead, using bdb2d_stack
	#
	# - select_mic_both     : Line edit box for formatted string type, and open file button  for micrograph/movie selection file (*.txt) with all files (*) (select_mic_both)
	# - select_mic_one      : Line edit box for formatted string type, and open file button  for micrograph selection file (*.txt) with all files (*) (select_mic_one)
	# - select_mic_one_ext  : Line edit box for formatted string type, and open file buttons for micrograph selection file (*.txt) with all files (*) (select_mic_one) and all supported 2D image formats (mic_one with bdb path parser)
	# - select_mic_stack    : Line edit box for formatted string type, and open file button  for movie selection file (*.txt) with all files (*) (select_mic_one)
	# - select_data2d_stack : Line edit box for formatted string type, and open file button  for 2D image selection file (*.txt) with all files (*) (select_data2d_stack)
	# - select_drift_params : Line edit box for formatted string type, and open file button  for drift shit parameters selection file (*.txt) with all files (*) (select_drift_params)
	#
	# - params_any_txt      : Line edit box for formatted string type, and open file button  for parameters text file (*.txt) with all files (*) (params_any_txt)
	# - params_proj_txt     : Line edit box for formatted string type, and open file button  for projection parameters text file (*.txt) with all files (*) (params_proj_txt)
	# - params_coords_any   : Line edit box for formatted string type, and open file buttons for EMAN1 box coordinates file (*.box) (params_coords_box) and all coordinates formats (params_coords_any)
	# - params_cter_txt     : Line edit box for formatted string type, and open file button  for CTER partres parameters text file (*.txt) (params_cter_txt)
	# - params_rebox_rbx    : Line edit box for formatted string type, and open file buttons for SPHIRE rebox file (*.rbx) (params_rebox_rbx)
	# - params_drift_txt    : Line edit box for formatted string type, and open file button  for drift shit parameters text file (*.txt) with all files (*) (params_drift_txt)
	# - params_relion_star  : Line edit box for formatted string type, and open file button  for RELION STAR file (*.star) with all files (*) (params_relion_star)
	# 
	# - spectrum1d          : Line edit box for formatted string type, and open file button  for 1D power spectrum file (*.txt) with all files (*) (spectrum1d) 
	# - mtf                 : Line edit box for formatted string type, and open file button  for MTF data file (*.txt) with all files (*) (mtf) 
	# - rot_matrix          : Line edit box for formatted string type, and open file button  for rotational matrix file (*.txt) with all files (*) (rot_matrix)
	# - pdb                 : Line edit box for formatted string type, and open file button  for PDB data file (*.pdb *.pdb*) (pdb)
	# - exe                 : Line edit box for formatted string type, and open file button  for executable file (*.exe) with all files (*) (exe)
	# 
	# - dir                 : Line edit box for formatted string type, and open directory button
	# - dir_list            : Line edit box for formatted string type, and open directory button
	#                         The string with space is interpreted as a list of any directory names upon command generation. (i.e. does not enclose the string with single quotes)
	# 
	# - user_func           : Two line edit boxes for formatted string type (function name & file path of the container script),
	#                         and open file button for .py
	# - abs_freq            : Line edit box for float type, label for float type of corresponding resolution [A], and convert button to show the calculator dialog 
	#
	# - apix                : Project constant - float type
	# - ctfwin              : Project constant - int type
	# - box                 : Project constant - int type
	# - radius              : Project constant - int type
	# - sym                 : Project constant - formatted string type
	# - mass                : Project constant - int type
	#
###	# - any_micrograph  : Line edit box for formatted string type, and open file buttons for all file types (also mrc, tiff, and etc) and .txt      <= Should be removed in near future (2018/01/25 Toshio Moriya)
###	# - any_image       : Line edit box for formatted string type, and open file buttons for all file types (also mrc, tiff, and etc) and .bdb      <= Should be removed in near future (2018/01/25 Toshio Moriya)
###	# - any_image_list  : Line edit box for formatted string type, and open file button  for all file types (also mrc, tiff, and etc) and .bdb.    <= Should be removed in near future (2018/01/25 Toshio Moriya)
###	#                     The string with space is interpreted as a list of any image file names upon command generation. (i.e. does not enclose the string with single quotes)
###	# - any_file        : Line edit box for formatted string type, and open file button  for all file types                                         <= Should be removed in near future (2018/01/25 Toshio Moriya)
###	# - any_file_list   : Line edit box for formatted string type, and open file button  for all file types and .bdb                                <= Should be removed in near future (2018/01/25 Toshio Moriya)
###	#                     The string with space is interpreted as a list of any image file names upon command generation. (i.e. does not enclose the string with single quotes)
###	# - image           : Line edit box for formatted string type, and open file buttons for .hdf and .bdb                                          <= Should be removed in near future (2018/01/25 Toshio Moriya)
###	# - bdb             : Line edit box for formatted string type, and open file button  for .bdb                                                   <= Should be removed in near future (2018/01/25 Toshio Moriya)
###	# - hdf             : Line edit box for formatted string type, and open file button  for .hdf                                                   <= Should be removed in near future (2018/01/25 Toshio Moriya)
###	# - mrc             : Line edit box for formatted string type, and open file button  for .mrc                                                   <= Should be removed in near future (2018/01/25 Toshio Moriya)
###	# - txt             : Line edit box for formatted string type, and open file button  for .txt
	# 

	keyword_dict = {}

	# Use priority 0 to overrule the exceptional cases (This is a reason why priority is introduced...)
	keyword_dict["--chunkdir"]                    = SXkeyword_map(0, "dir")                 # --chunkdir=chunkdir (contains keyworkd 'chunk' but this should be directory type)
	keyword_dict["destination_directory"]         = SXkeyword_map(0, "dir")                 # destination_directory (contains keyworkd 'directory' but this should be directory type)
	keyword_dict["main000"]         = SXkeyword_map(0, "dir")                 # destination_directory (contains keyworkd 'directory' but this should be directory type)
###	keyword_dict["--use_latest_master_directory"] = SXkeyword_map(0, "bool")                # --use_latest_master_directory (contains keyworkd 'directory' but this should be bool type)
	keyword_dict["--stack_mode"]                  = SXkeyword_map(0, "bool")                # stack_mode (contains keyworkd 'stack' but this should be bool type)
	keyword_dict["--generate_mask"]               = SXkeyword_map(0, "bool")                # --generate_mask (contains keyworkd 'mask' but this should be bool type)
	keyword_dict["--binary_mask"]                 = SXkeyword_map(0, "bool")                # --binary_mask (contains keyworkd 'mask' but this should be bool type)
	keyword_dict["--symmetrize"]                  = SXkeyword_map(0, "bool")                # --symmetrize (contains keyworkd '--sym' but this should be bool type)
	keyword_dict["--do_adaptive_mask"]            = SXkeyword_map(0, "bool")                # --do_adaptive_mask (contains keyworkd 'mask' but this should be bool type)
	keyword_dict["--adaptive_mask"]               = SXkeyword_map(0, "bool")                # --adaptive_mask (contains keyworkd 'mask' but this should be bool type)
###	keyword_dict["--skip_create_substack"]        = SXkeyword_map(0, "bool")                # --skip_create_substack (contains keyworkd 'stack' but this should be bool type)
	keyword_dict["--no_virtual_stack"]            = SXkeyword_map(0, "bool")                # --no_virtual_stack (contains keyworkd 'stack' but this should be bool type)
	keyword_dict["--create_stack"]                = SXkeyword_map(0, "bool")                # --create_stack (contains keyworkd 'stack' but this should be bool type)
	keyword_dict["--do_not_create_stack"]         = SXkeyword_map(0, "bool")                # --do_not_create_stack (contains keyworkd 'stack' but this should be bool type)
	keyword_dict["--save_vstack"]                 = SXkeyword_map(0, "bool")                # --save_vstack (contains keyworkd 'stack' but this should be bool type)
	keyword_dict["--use_mol_mass"]                = SXkeyword_map(0, "bool_ignore")                # if one wants to use the --mol_mass=KILODALTON
	keyword_dict["--s_use_mol_mass"]                = SXkeyword_map(0, "bool_ignore")                # if one wants to use the --mol_mass=KILODALTON
	keyword_dict["--use_second_mask"]             = SXkeyword_map(0, "bool_ignore")                # if one wants to use the --mol_mass=KILODALTON
	keyword_dict["--mask_threshold"]              = SXkeyword_map(0, "float")               # --mask_threshold=MASK_THRESHOLD (contains keyworkd 'mask' but this should be bool type)
	# Use priority 1 for output
	keyword_dict["output"]                        = SXkeyword_map(1, "output")              # output_hdf, output_directory, outputfile, outputfile, --output=OUTPUT, output_stack, output_file
	keyword_dict["outdir"]                        = SXkeyword_map(1, "output")              # outdir
	keyword_dict["--ave3D"]                       = SXkeyword_map(1, "output")              # --ave3D
	keyword_dict["--var3D"]                       = SXkeyword_map(1, "output")              # --var3D
	keyword_dict["--ave2D"]                       = SXkeyword_map(1, "output")              # --ave2D
	keyword_dict["--var2D"]                       = SXkeyword_map(1, "output")              # --var2D
###	keyword_dict["--masterdir"]                   = SXkeyword_map(1, "output")              # --masterdir=master_dir
	keyword_dict["locres_volume"]                 = SXkeyword_map(1, "output")              # locres_volume (this contained keyword "volume" also... This is another reason why priority is introduced...)
	keyword_dict["directory"]                     = SXkeyword_map(1, "output")              # directory
	keyword_dict["rotpw"]                         = SXkeyword_map(1, "output")              # rotpw
	keyword_dict["--makevstack"]                  = SXkeyword_map(1, "output_bdb2d_stack")  # --makevstack
	keyword_dict["input_micrograph_list"]         = SXkeyword_map(1, "mic_one_list")        # input_micrograph_list (contains keyword 'input_micrograph' but this should be image_list type)
	keyword_dict["--ctref_orgstack"]              = SXkeyword_map(1, "bdb2d_stack")         # --ctref_orgstack=stack_for_continuation
	keyword_dict["input_bdb_stack_path"]          = SXkeyword_map(1, "bdb2d_stack")         # input_bdb_stack_path (contains keyword 'stack' but this should be bdb type)
	keyword_dict["--substack_basename"]           = SXkeyword_map(1, "string")              # --substack_basename=SUBSTACK_BASENAME (contains keyword 'volume' but this should be string type)
	keyword_dict["--importctf"]                   = SXkeyword_map(1, "params_cter_txt")     # --importctf=IMPORTCTF (contains keyword '--import' but this should be params_cter_txt type) # NOTE: Toshio Moriya 2018/01/27 This is used by sxprocess.py but this function might be obsoleted at this point
	# Use priority 2 for the others
	keyword_dict["input_data_list"]               = SXkeyword_map(2, "displayable_list")    # input_data_list
	keyword_dict["--input"]                       = SXkeyword_map(2, "data2d3d_both")       # --input=INPUT
	keyword_dict["input_image_path"]              = SXkeyword_map(2, "mic_one")             # input_image_path
	keyword_dict["source_micrograph"]             = SXkeyword_map(2, "mic_both")            # source_micrograph_pattern
	keyword_dict["input_micrograph"]              = SXkeyword_map(2, "mic_stack")           # input_micrograph_pattern 
	keyword_dict["isac_class_avgs_path"]          = SXkeyword_map(2, "data2d_one")          # isac_class_avgs_path
	keyword_dict["stack"]                         = SXkeyword_map(2, "bdb2d_stack")         # stack, prj_stack, input_stack, --instack=input_stack_file
	keyword_dict["volume"]                        = SXkeyword_map(2, "data3d_one")          # initial_volume, firstvolume, secondvolume, input_volume
	keyword_dict["mask"]                          = SXkeyword_map(2, "data3d_one")          # --mask3D=mask3D, maskfile, mask, --mask=MASK, --mask3D=mask3d_file
	keyword_dict["--focus"]                       = SXkeyword_map(2, "data3d_one")          # --focus=3D_focus_mask, --focus=focus3d_file
###	keyword_dict["isac_averages"]                 = SXkeyword_map(2, "data2d_one")          # isac_averages
	keyword_dict["--ctref_initvol"]               = SXkeyword_map(2, "data3d_one")          # --ctref_initvol=restarting_initial_volume
	keyword_dict["selection_list"]                = SXkeyword_map(2, "select_mic_stack")    # selection_list
	keyword_dict["--chunk"]                       = SXkeyword_map(2, "select_data2d_stack") # --chunk0=CHUNK0_FILE_NAME, --chunk1=CHUNK1_FILE_NAME
	keyword_dict["--ctref_subset"]                = SXkeyword_map(2, "select_data2d_stack") # --ctref_subset=selection_file_path
	keyword_dict["--list"]                        = SXkeyword_map(2, "select_data2d_stack") # --list
###	keyword_dict["--subset"]                      = SXkeyword_map(2, "select_data2d_stack") # --subset=subset_file_path
	keyword_dict["input_shift_list_file"]         = SXkeyword_map(2, "select_drift_params") # input_shift_list_file
###	keyword_dict["--resample_ratio_source"]       = SXkeyword_map(2, "params_any_txt")      # --resample_ratio_source
	keyword_dict["params_file"]                      = SXkeyword_map(2, "params_any_txt")      # --import=INPUT_PARAMS_PATH
	keyword_dict["--import"]                      = SXkeyword_map(2, "params_any_txt")      # --import=INPUT_PARAMS_PATH
	keyword_dict["--export"]                      = SXkeyword_map(2, "params_any_txt")      # --export==OUTPUT_PARAMS_FILE
	keyword_dict["--params_2d_file"]			  = SXkeyword_map(2, "params_any_txt")
	keyword_dict["--params_3d_file"]			  = SXkeyword_map(2, "params_any_txt")
	keyword_dict["--params_3d_index_file"]		  = SXkeyword_map(2, "params_any_txt")
	keyword_dict["--params_3d_chunk_file_0"] 	  = SXkeyword_map(2, "params_any_txt")
	keyword_dict["--params_3d_chunk_file_1"] 	  = SXkeyword_map(2, "params_any_txt")
	keyword_dict["input_coordinates_pattern"]     = SXkeyword_map(2, "params_coords_any")   # input_coordinates_pattern
	keyword_dict["input_rebox_pattern"]           = SXkeyword_map(2, "params_rebox_rbx")    # input_rebox_pattern
	keyword_dict["input_ctf_params_source"]       = SXkeyword_map(2, "params_cter_txt")     # input_ctf_params_source
	keyword_dict["cter_ctf_file"]                 = SXkeyword_map(2, "params_cter_txt")     # cter_ctf_file
	keyword_dict["swap_ctf_params"]               = SXkeyword_map(2, "params_cter_txt")     # swap_ctf_params
	keyword_dict["inputfile"]                     = SXkeyword_map(2, "params_drift_txt")    # inputfile
	keyword_dict["input_shift_pattern"]           = SXkeyword_map(2, "params_drift_txt")    # input_shift_pattern
	keyword_dict["input_star_file"]               = SXkeyword_map(2, "params_relion_star")  # input_star_file
	keyword_dict["--modelpw"]                     = SXkeyword_map(2, "spectrum1d")          # --modelpw=pw2_model_txt
	keyword_dict["--pwreference"]                 = SXkeyword_map(2, "spectrum1d")          # --pwreference=pwreference
	keyword_dict["--PWadjustment"]                = SXkeyword_map(2, "spectrum1d")          # --PWadjustment=PWadjustment, --PWadjustment=ref_pwspectrum1d_file 
	keyword_dict["--mtf"]                         = SXkeyword_map(2, "mtf")                 # --mtf=MTF_FILE_NAME
	keyword_dict["--tr0"]                         = SXkeyword_map(2, "rot_matrix")          # --tr0=matrix_file
	keyword_dict["input_pdb"]                     = SXkeyword_map(2, "pdb")                 # input_pdb
	keyword_dict["unblur_path"]                   = SXkeyword_map(2, "exe")                 # unblur_path
	keyword_dict["summovie_path"]                 = SXkeyword_map(2, "exe")                 # summovie_path
	keyword_dict["--isac_dir"]                    = SXkeyword_map(2, "dir")                 # --isac_dir
	keyword_dict["input_run_dir"]                 = SXkeyword_map(2, "dir")                 # input_run_dir
###	keyword_dict["--oldrefdir"]                   = SXkeyword_map(2, "dir")                 # --oldrefdir=refine_dir_path
	keyword_dict["--ctref_oldrefdir"]             = SXkeyword_map(2, "dir")                 # --ctref_oldrefdir=refine_dir_path
	keyword_dict["--refinement_dir"]              = SXkeyword_map(2, "dir")                 # --refinement_dir=refinemen_out_dir
	keyword_dict["--previous_run"]                = SXkeyword_map(2, "dir")                 # --previous_run1=sort3d_run1_directory, --previous_run2=sort3d_run2_directory (need for sort3d.txt)
	keyword_dict["--relion_project_dir"]          = SXkeyword_map(2, "dir")                 # --relion_project_dir=DIR_PATH
	keyword_dict["--function"]                    = SXkeyword_map(2, "user_func")           # --function=user_function
	keyword_dict["--function_ai"]                    = SXkeyword_map(2, "user_func")           # --function=user_function
	# keyword_dict["--FL"]                        = SXkeyword_map(2, "abs_freq")            # (--FL=FL); This is ISAC2 & advanced that We do NOT supported at this point because of pixel size reduction by target radius
	# keyword_dict["--FH"]                        = SXkeyword_map(2, "abs_freq")            # (--FH=FH); This is ISAC2 & advanced that We do NOT supported at this point because of pixel size reduction by target radius
	# keyword_dict["--FF"]                        = SXkeyword_map(2, "abs_freq")            # (--FF=FF); This is ISAC2 & falloff that we do NOT supported at this point (see below NOTE: 2018/02/06 Toshio Moriya)
	keyword_dict["--fl"]                          = SXkeyword_map(2, "abs_freq")            # --fl=cutoff_frequency, (--fl=fl), (--fl=fl), (--fl=lpf_cutoff_freq), --fl=LPF_CUTOFF_FREQ
	# keyword_dict["--aa"]                        = SXkeyword_map(2, "abs_freq")            # (--aa=aa), (--aa=aa), (--aa=lpf_falloff), (--aa=LPF_FALLOFF_WIDTH); This is falloff that we do NOT supported at this point (see below NOTE: 2018/02/06 Toshio Moriya)
	keyword_dict["--fh"]                          = SXkeyword_map(2, "abs_freq")            # --fh=fh
	keyword_dict["--res_overall"]                 = SXkeyword_map(2, "abs_freq")            # --res_overall=overall_resolution
	# keyword_dict["--falloff"]                   = SXkeyword_map(2, "abs_freq")            # --falloff=falloff; This is falloff that we do NOT supported at this point (see below NOTE: 2018/02/06 Toshio Moriya)

	keyword_dict["--apix"]                        = SXkeyword_map(2, "apix")                # --apix=pixel_size, --apix, --apix=PIXEL_SIZE
	keyword_dict["--pixel_size"]                  = SXkeyword_map(2, "apix")                # --pixel_size=PIXEL_SIZE
	keyword_dict["--wn"]                          = SXkeyword_map(2, "ctfwin")              # --wn
	keyword_dict["--box"]                         = SXkeyword_map(2, "box")                 # --box=box_size, --box_size=box_size, --boxsize=BOX_SIZE
	keyword_dict["--rb_box_size"]                 = SXkeyword_map(2, "box")                 # --rb_box_size=BOX_SIZE
	keyword_dict["--radius"]                      = SXkeyword_map(2, "radius")              # --radius=particle_radius, --radius=outer_radius, --radius=outer_radius, --radius=particle_radius, --radius=outer_radius, --radius=outer_radius
	keyword_dict["--sym"]                         = SXkeyword_map(2, "sym")                 # --sym=c1, --sym=symmetry, --sym=c4
	keyword_dict["--mol_mass"]                    = SXkeyword_map(2, "mass")                # --mol_mass=KILODALTON
	#keyword_dict["--symmetry"]                    = SXkeyword_map(2, "sym")                 # --symmetry
	#keyword_dict["--max_occupy"]                  = SXkeyword_map(2, "int")                 # --max_occupy

	#Added keywords for crYOLO
	keyword_dict["--input_size"] = SXkeyword_map(0, "int")
	keyword_dict["training_dir"] = SXkeyword_map(2, "dir")
	keyword_dict["annot_dir"] = SXkeyword_map(2, "dir")
	keyword_dict["valid_image_dir"] = SXkeyword_map(2, "dir")
	keyword_dict["valid_annot_dir"] = SXkeyword_map(2, "dir")
	keyword_dict["cryolo_train_path"] = SXkeyword_map(2, "exe")
	keyword_dict["cryolo_predict_path"] = SXkeyword_map(2, "exe")
	keyword_dict["target_dir"] = SXkeyword_map(2, "dir")
	keyword_dict["config_path"] = SXkeyword_map(2, "params_any_json")  # --import=INPUT_PARAMS_PATH
	keyword_dict["model_path"] = SXkeyword_map(2, "params_any_h5")
	keyword_dict["pretrained_weights_name"] = SXkeyword_map(2, "params_any_h5")
	keyword_dict["box_dir"] = SXkeyword_map(2, "dir")

	# Added keywords for ctf refine
	keyword_dict["inputstack"] = SXkeyword_map(2, "data2d_stack")
	keyword_dict["path_to_input_stack"] = SXkeyword_map(2, "data2d_stack")
	keyword_dict["outstack"] = SXkeyword_map(2, "dir")
	keyword_dict["output_directory"] = SXkeyword_map(2, "dir")
	keyword_dict["refinement_dir"] = SXkeyword_map(2, "dir")
	keyword_dict["volume_path"] = SXkeyword_map(2, "data3d_one")


	keyword_dict["--range"] = SXkeyword_map(2, "float")
	keyword_dict["--delta"] = SXkeyword_map(2, "float")
	keyword_dict["--resolution"] = SXkeyword_map(2, "float")
	keyword_dict["--number_part"] = SXkeyword_map(2, "int")
	# NOTE: 2018/02/06 Toshio Moriya
	# Low-pass filter fall-off width does not make sense to convert to resolution [A] directly. 
	# It might make more sense to compute Angstrom range from the given cutoff, falloff width, and pixel size
	# For example, fall-off range of 0.1 [1/pix] fall-off width with 2.5[A] cutoff is
	# 1.0[A/Pix]/2.5[A] + 0.1 [1/pix]/2 = 0.45 [1/pix] => 2.2222[A]
	# 1.0[A/Pix]/2.5[A] - 0.1 [1/pix]/2 = 0.35 [1/pix] => 2.8571[A]
	# 
	# Reversely, we can also convert from Angstrom range (e.g. Width from 2.0 to 4.0[A]@1.0[A/Pix]) to absolute frequenyc [A]
	# 1.0[A/Pix]/2.0[A] - 1.0[A/Pix]/4.0[A] = 0.5[1/pix] - 0.25[1/pix] = 0.25[1/pix]
	# In this case, the cut off will be around 
	# 1.0[A/Pix]/4.0[A] + 0.25[1/pix]/2 = 0.25[1/pix] + 0.25[1/pix]/2 = 0.375 [1/Pix] =>  2.6667[A]
	# 
	# NOTE: 2018/02/06 Toshio Moriya
	# At this point, the options related to abs_freq are:
	# sxisac2.py (with target_radius): [MAIN] [ADVINCED] (--FL), (--FH), (--FF fall off)
	# sxcompute_isac_avg.py (with apix): [MAIN] --fl [ADVINCED] --fh
	# sxrviper.py (with target_radius & apix for moon?): [MAIN] [ADVINCED] (--fl), (--aa fall off)
	# sxviper.py (with target_radius & apix for moon?): [MAIN] [ADVINCED] (--fl), (--aa fall off)
	# sxprocess.py --combinemap (with apix): [MAIN] (--fl in [A]) [ADVINCED] (--aa fall off)
	# sx3dvariability.py (with decimate): [MAIN] --fl (--aa fall off) [ADVINCED]
	# sxlocres.py (with apix): [MAIN] --res_overall [ADVINCED]
	# sxfilterlocal.py: [MAIN] (--falloff fall off) [ADVINCED]
	# 
	# NOTE: 2016/02/23 Toshio Moriya
	# Below might be useful to include for project settings
	# reference power spectrum? --pwreference of viper, --pwreference of rviper, --PWadjustment of sort3d, --PWadjustment of rsort3d
	#
	# Below must be exceptional cases
	# --wn of locres, sort3d, & rsort3d; same as ctfwin?
	# --radius of locres & filterlocal; same as radius?
	#
	
	return keyword_dict

# ----------------------------------------------------------------------------------------
def handle_exceptional_cases(sxcmd):
	# DESIGN_NOTE: 2016/02/05 Toshio Moriya
	# Handle exceptional cases due to the limitation of software design
	# In future, we should remove these exception handling by reviewing the design
	if sxcmd.name == "sp_window":
		assert(sxcmd.token_dict["input_micrograph_pattern"].key_base == "input_micrograph_pattern")
		assert(sxcmd.token_dict["input_micrograph_pattern"].type == "mic_stack")
		sxcmd.token_dict["input_micrograph_pattern"].type = "mic_one"
		assert(sxcmd.token_dict["selection_list"].key_base == "selection_list")
		assert(sxcmd.token_dict["selection_list"].type == "select_mic_stack")
		sxcmd.token_dict["selection_list"].type = "select_mic_one"
	elif sxcmd.name == "sp_rewindow":
		assert(sxcmd.token_dict["input_micrograph_pattern"].key_base == "input_micrograph_pattern")
		assert(sxcmd.token_dict["input_micrograph_pattern"].type == "mic_stack")
		sxcmd.token_dict["input_micrograph_pattern"].type = "mic_one"
		assert(sxcmd.token_dict["selection_list"].key_base == "selection_list")
		assert(sxcmd.token_dict["selection_list"].type == "select_mic_stack")
		sxcmd.token_dict["selection_list"].type = "select_mic_one"
	elif sxcmd.name == "sp_cter":
		assert(sxcmd.token_dict["selection_list"].key_base == "selection_list")
		assert(sxcmd.token_dict["selection_list"].type == "select_mic_stack")
		sxcmd.token_dict["selection_list"].type = "select_mic_one_ext"
	elif sxcmd.name == "sp_isac2":
		assert(sxcmd.token_dict["output_directory"].key_base == "output_directory")
		assert(sxcmd.token_dict["output_directory"].type == "output")
		sxcmd.token_dict["output_directory"].type = "output_continue"
	elif sxcmd.name == "sp_isac":
		assert(sxcmd.token_dict["output_directory"].key_base == "output_directory")
		assert(sxcmd.token_dict["output_directory"].type == "output")
		sxcmd.token_dict["output_directory"].type = "output_continue"
	elif sxcmd.name == "sp_mask":
		sxcmd.token_dict["output_directory"].type = "dir"
	elif sxcmd.name == "sp_rviper":
		assert(sxcmd.token_dict["stack"].key_base == "stack")
		assert(sxcmd.token_dict["stack"].type == "bdb2d_stack")
		sxcmd.token_dict["stack"].type = "data2d_one"
		# Typically, this is target particle radius used by ISAC2.
		assert(sxcmd.token_dict["radius"].key_base == "radius")
		assert(sxcmd.token_dict["radius"].type == "radius")
		sxcmd.token_dict["radius"].type = "int"
###		assert(sxcmd.token_dict["output_directory"].key_base == "output_directory")
###		assert(sxcmd.token_dict["output_directory"].type == "output")
##		sxcmd.token_dict["output_directory"].type = "output_continue"
		assert(sxcmd.token_dict["directory"].key_base == "directory")
		assert(sxcmd.token_dict["directory"].type == "output")
		sxcmd.token_dict["directory"].type = "output_continue"
		assert(sxcmd.token_dict["fl"].key_base == "fl")
		assert(sxcmd.token_dict["fl"].type == "abs_freq")
		sxcmd.token_dict["fl"].type = "float"
	elif sxcmd.name == "sp_separate_class":
		sxcmd.token_dict["input_class_avgs"].type = "data2d_one"
	elif sxcmd.name == "sp_proj_compare":
		sxcmd.token_dict["stack"].type = "data2d_one"
		sxcmd.token_dict["classangles"].type = "params_proj_txt"
		sxcmd.token_dict["classselect"].type = "select_data2d_stack"
		sxcmd.token_dict["partangles"].type = "params_proj_txt"
		sxcmd.token_dict["partselect"].type = "select_data2d_stack"
	elif sxcmd.name == "sp_viper":
		assert(sxcmd.token_dict["stack"].key_base == "stack")
		assert(sxcmd.token_dict["stack"].type == "bdb2d_stack")
		sxcmd.token_dict["stack"].type = "data2d_one"
		# Typically, this is target particle radius used by ISAC2.
		assert(sxcmd.token_dict["radius"].key_base == "radius")
		assert(sxcmd.token_dict["radius"].type == "radius")
		sxcmd.token_dict["radius"].type = "int"
		assert(sxcmd.token_dict["directory"].key_base == "directory")
		assert(sxcmd.token_dict["directory"].type == "output")
		sxcmd.token_dict["directory"].type = "output_continue"
		assert(sxcmd.token_dict["fl"].key_base == "fl")
		assert(sxcmd.token_dict["fl"].type == "abs_freq")
		sxcmd.token_dict["fl"].type = "float"
	elif sxcmd.name == "sp_meridien":
		assert(sxcmd.token_dict["output_directory"].key_base == "output_directory")
		assert(sxcmd.token_dict["output_directory"].type == "output")
		sxcmd.token_dict["output_directory"].type = "output_continue"
	elif sxcmd.name == "sp_meridien_20171120":
		assert(sxcmd.token_dict["output_directory"].key_base == "output_directory")
		assert(sxcmd.token_dict["output_directory"].type == "output")
		sxcmd.token_dict["output_directory"].type = "output_continue"
	elif sxcmd.name == "sp_sort3d_depth":
		assert(sxcmd.token_dict["output_dir"].key_base == "output_dir")
		assert(sxcmd.token_dict["output_dir"].type == "output")
		sxcmd.token_dict["output_dir"].type = "output_continue"
	elif sxcmd.name == "sp_sort3d":
		assert(sxcmd.token_dict["wn"].key_base == "wn")
		assert(sxcmd.token_dict["wn"].type == "ctfwin")
		sxcmd.token_dict["wn"].type = "int"
		assert(sxcmd.token_dict["outdir"].key_base == "outdir")
		assert(sxcmd.token_dict["outdir"].type == "output")
		sxcmd.token_dict["outdir"].type = "output_continue"
	elif sxcmd.name == "sp_locres":
		assert(sxcmd.token_dict["wn"].key_base == "wn")
		assert(sxcmd.token_dict["wn"].type == "ctfwin")
		sxcmd.token_dict["wn"].type = "int"
		assert(sxcmd.token_dict["directory"].key_base == "directory")
		assert(sxcmd.token_dict["directory"].type == "output")
		sxcmd.token_dict["directory"].type = "output_continue"
	elif sxcmd.name == "sp_filterlocal":
		assert(sxcmd.token_dict["locres_volume"].key_base == "locres_volume")
		assert(sxcmd.token_dict["locres_volume"].type == "output")
		sxcmd.token_dict["locres_volume"].type = "data3d_one"
	elif sxcmd.name == "sp_pipe":
		if sxcmd.subname == "organize_micrographs":
			assert(sxcmd.token_dict["selection_list"].key_base == "selection_list")
			assert(sxcmd.token_dict["selection_list"].type == "select_mic_stack")
			sxcmd.token_dict["selection_list"].type = "select_mic_both"
		elif sxcmd.subname == "moon_eliminator":
			assert(sxcmd.token_dict["resample_ratio"].key_base == "resample_ratio")
			assert(sxcmd.token_dict["resample_ratio"].type == "string")
			sxcmd.token_dict["resample_ratio"].type = "dir"
			assert(sxcmd.token_dict["fl"].key_base == "fl")
			assert(sxcmd.token_dict["fl"].type == "abs_freq")
			sxcmd.token_dict["fl"].type = "float"

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
			sp_global_def.ERROR("Wiki Format Warning: The string \"%s\" contains \"%s\" but not \"%s\". Removing \"%s\", but please check the format in Wiki document." % (target_text, makeup_begin, makeup_end, makeup_begin), "%s in %s" % (__name__, os.path.basename(__file__)))
			target_text = target_text.replace(makeup_begin, "", 1)
		else: # assert (item_tail > -1)
			makeup_token = target_text[item_head:item_tail+len(makeup_end)]
			display_item = makeup_token.replace(makeup_begin, "")
			display_item = display_item.replace(makeup_end, "")
			if display_item.find(makeup_separator) != -1:
				item_tokens = display_item.split(makeup_separator)
				assert (len(item_tokens) == 2)
				display_item = item_tokens[1] # 2nd one should be display text
			print("### Found a wiki makeup token \"%s\". Changed to \"%s\"" % (makeup_token, display_item))
			target_text = target_text.replace(makeup_token, display_item, 1)

		# Try to find the next
		item_head = target_text.find(makeup_begin)
	
	target_text = target_text.replace("\'\'", "")
	target_text = target_text.replace("\'\'\'", "")
	target_text = target_text.replace("<", "&lt;")  # rich text tag
	target_text = target_text.replace(">", "&gt;")  # rich text tag
	
	return target_text

# ----------------------------------------------------------------------------------------
def construct_token_list_from_MoinMoinWiki(sxcmd_config):

	print("Start parsing MoinMoinWiki document (%s as %s %s command) " % (sxcmd_config.wiki, sxcmd_config.category, sxcmd_config.role))

	if sxcmd_config.format != "MoinMoinWiki": sp_global_def.ERROR("Logical Error: Incorrect Wiki format %s! Check the sxcmd_config setting in this script." % (sxcmd_config.format), "%s in %s" % (__name__, os.path.basename(__file__)))

	# Allocate memory for new SXcmd instance
	sxcmd = sxgui_template.SXcmd(sxcmd_config.category, sxcmd_config.role, sxcmd_config.is_submittable)

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
	if os.path.exists(sxcmd_config.wiki) == False: sp_global_def.ERROR("Rutime Error: Wiki document is not found.", "%s in %s" % (__name__, os.path.basename(__file__)))

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
			if current_state != state_processing: sp_global_def.ERROR("Logical Error: This condition should not happen! State setting must be incorrect.", "%s in %s" % (__name__, os.path.basename(__file__)))
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
					if current_section >= len(section_lists): sp_global_def.ERROR("Logical Error: This condition should not happen! Section setting must be incorrect.", "%s in %s" % (__name__, os.path.basename(__file__)))
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
					if item_tail == -1: sp_global_def.ERROR("Wiki Format Error: '= Name =' section should contain only one valid line, and the line should starts from 'sx* - ' or 'e2* - ': %s" % line_wiki, "%s in %s" % (__name__, os.path.basename(__file__)))
					sxcmd_string = line_buffer[0:item_tail].strip()
					sxcmd_string_token_list = sxcmd_string.split()
					n_sxcmd_string_token_list = len(sxcmd_string_token_list)
					sxcmd.name = sxcmd_string_token_list[0]
					if n_sxcmd_string_token_list == 2:
						sxcmd.subname = sxcmd_string_token_list[1]
					elif n_sxcmd_string_token_list > 2:
						sp_global_def.ERROR("Wiki Format Error: Command string should be only script file name or script file name plus subcommand", "%s in %s" % (__name__, os.path.basename(__file__)))
					line_buffer = line_buffer[item_tail + len(target_operator):].strip() # Get the rest of line
					# Extract the label of this sxscript
					target_operator = ":"
					item_tail = line_buffer.find(target_operator)
					if item_tail == -1: sp_global_def.ERROR("Wiki Format Error: '= Name =' section should contain a label ended with ':' after 'sx* - ' or 'e2* - ': %s" % line_wiki, "%s in %s" % (__name__, os.path.basename(__file__)))
					sxcmd.label = line_buffer[0:item_tail].strip()
					# Extract the short info about this sxscript (can be empty)
					sxcmd.short_info = remove_MoinMoinWiki_makeup(line_buffer[item_tail + len(target_operator):].strip()) # Get the rest of line
				elif current_section == section_usage:
					# Extract 'usage in command line' to identify each command token is either an argument (no-prefix) or option ('--' prefix)
					# This information is also used to check consistency between 'usage in command line' and list in '== Input ==' and '== Output ==' sections
					if line_wiki[0:len("sx_")] == "sx_" or line_wiki[0:len("e2")] == "e2":
						usage_token_list = line_wiki.split()
						head_token_idx = 1
						if usage_token_list[0] != sxcmd.name + ".py": sp_global_def.ERROR("Wiki Format Error: First token should be script name with .py (sx*.py or e2*.py)", "%s in %s" % (__name__, os.path.basename(__file__)))
						if sxcmd.subname != "":
							head_token_idx = 2
							if usage_token_list[1] != sxcmd.subname: sp_global_def.ERROR("Wiki Format Error: Second token of this command should be subname", "%s in %s" % (__name__, os.path.basename(__file__)))
						# Register arguments and options
						for usage_token in usage_token_list[head_token_idx:]:
							# Check if --MPI is used in this script
							# NOTE: 2015/11/12 Toshio Moriya
							# The following can be removed when --MPI flag is removed from all sx*.py scripts
							if usage_token == "--MPI":
								# ERROR("Warning: The 'usage in command line' contains --MPI flag. The flag will be removed in near future, so ignoring this line...'.","%s in %s" % (__name__, os.path.basename(__file__)), action = 0)
								sxcmd.mpi_support = True
								sxcmd.mpi_add_flag = True
								continue
							# Allocate memory for new command token
							token = sxgui_template.SXcmd_token()
							# Extract key of command token.
							key = usage_token.split("=")[0] # Remove characters after '=' if token contains it (i.e. some options)
							token.key_base = key.strip("-") # Get key base name by removing prefix ('--' or '-' for option)
							token.key_prefix = key[0:len(key) - len(token.key_base)]
							# Try to set the special type base on the keyword dictionary
							best_keyword_map = SXkeyword_map(99, "")
							for keyword in list(keyword_dict.keys()):
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
					if sxcmd.mpi_support == False and line_wiki.find(target_operator) > -1 and line_wiki.find(sxcmd.name + ".py") > -1:
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
							if sxcmd.mpi_support == False or sxcmd.mpi_add_flag == False: sp_global_def.ERROR("Logical Error: Since MPI flag is found, the command should support MPI.", "%s in %s" % (__name__, os.path.basename(__file__)))
							continue
						line_buffer = line_buffer[item_tail + len(target_operator):].strip() # Get the rest of line
						# check consistency between 'usage in command line' and this
						if key_base not in list(sxcmd.token_dict.keys()): sp_global_def.ERROR("Wiki Format Error: Key base (%s) is missing from 'usage in command line' in '= Usage ='." % key_base, "%s in %s" % (__name__, os.path.basename(__file__)))
						# Get the reference to the command token object associated with this key base name
						token = sxcmd.token_dict[key_base]
						if token.key_base != key_base: sp_global_def.ERROR("Logical Error: Registered command token with wrong key base name into the dictionary.", "%s in %s" % (__name__, os.path.basename(__file__)))
						token.is_in_io = True # Set flag to tell this command token is find in input or output section
						token.group = group  # Set group of command token according to the current subsection
						# Extract label of command token
						target_operator = ":"
						item_tail = line_buffer.find(target_operator)
						if item_tail == -1: sp_global_def.ERROR("Wiki Format Error: This line (%s) is missing label. Please check the format in Wiki document." % line_wiki, "%s in %s" % (__name__, os.path.basename(__file__)))
						token.label = line_buffer[0:item_tail]
						line_buffer = line_buffer[item_tail + len(target_operator):].strip() # Get the rest of line
						# Extract help of command token before default value
						target_operator = "(default"
						item_tail = line_buffer.find(target_operator)
						if item_tail == -1: sp_global_def.ERROR("Wiki Format Error: This line (%s) is missing default setting. Please check the format in Wiki document." % line_wiki, "%s in %s" % (__name__, os.path.basename(__file__)))
						token.help = remove_MoinMoinWiki_makeup(line_buffer[0:item_tail])
						line_buffer = line_buffer[item_tail + len(target_operator):].strip() # Get the rest of line
						# Extract default value of command token
						target_operator = ")"
						item_tail = line_buffer.find(target_operator)
						if item_tail == -1: sp_global_def.ERROR("Wiki Format Error: This line (%s) is missing ')' for default setting. Please check the format in Wiki document." % line_wiki,"%s in %s" % (__name__, os.path.basename(__file__)))
						default_value = line_buffer[0:item_tail].strip() # make sure spaces & new line are not included at head and tail
						if default_value.find("required") != -1:
							# This is a required command token and should have value type instead of default value
							token.is_required = True
							token.default = ""
							if not token.type:
								# Type is still empty, meaning no special type is assigned
								# Extract the data type (the rest of line)
								token.type = default_value.replace("required", "").strip()
							assert (token.is_required)
							assert (not token.is_locked)
							assert (not token.is_reversed)
						elif default_value.find("question reversed in GUI") != -1:
							# token.is_required = False
							token.default = default_value.replace("question reversed in GUI", "").strip()
							token.type = "bool"
							if token.default == "True":
								token.default = False # convert the default value to opposite boolean
							elif token.default == "False":
								token.default = True # convert the default value to opposite boolean
							else:
								assert (False)
							assert (not token.is_required)
							assert (not token.is_locked)
							assert (not token.is_reversed)
						elif default_value.find("value reversed in GUI") != -1:
							# token.is_required = False
							token.is_reversed = True
							token.default = default_value.replace("value reversed in GUI", "").strip()
							token.type = "bool"
							if token.default == "True":
								token.default = False # convert the default value to opposite boolean
							elif token.default == "False":
								token.default = True # convert the default value to opposite boolean
							else:
								assert (False)
							assert (not token.is_required)
							assert (not token.is_locked)
							assert (token.is_reversed)
						else:
							# This is not required command token and should have default value
							# token.is_required = False
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
							assert (not token.is_required)
							assert (not token.is_locked)
							assert (not token.is_reversed)
						### Initialise restore value with default value
						### if not token.is_locked:
						### 	token.restore = token.default
						# 
						# NOTE: Toshio Moriya 2018/02/07
						# Currently, lock status can be set only through the subcommand configuration
						# Therefore, the lock should be always False yet at this point
						assert (not token.is_locked)
						token.restore = token.default
						# Ignore the rest of line ...
				else:
					sp_global_def.ERROR("Logical Error: This section is invalid. Did you assign an invalid section?", "%s in %s" % (__name__, os.path.basename(__file__)))

	if current_state != state_done: sp_global_def.ERROR("Wiki Format Error: parser could not extract all information necessary (current state %d). Please check if the Wiki format has all required sections." % (current_state), "%s in %s" % (__name__, os.path.basename(__file__)))

	# Make sure there are no extra arguments or options in 'usage in command line' of '= Usage ='
	for token in sxcmd.token_list:
		if token.is_in_io == False: sp_global_def.ERROR("Wiki Format Error: An extra argument or option (%s) is found in 'usage in command line' of '= Usage ='." % token.key_base, "%s in %s" % (__name__, os.path.basename(__file__)))

	file_wiki.close()

	handle_exceptional_cases(sxcmd)

	print("Succeeded to parse MoinMoinWiki document (%s as %s %s command)" % (sxcmd_config.wiki, sxcmd_config.category, sxcmd_config.role))

	"""
	# For DEBUG
	if sxcmd.name == "sxwindow":
		print "><><>< DEBUG OUTPUT ><><><"
		print ""
		print "------"
		print "GLOBAL"
		print "------"
		print "name            : %s" % sxcmd.name
		print "subname         : %s" % sxcmd.subname
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
			sp_global_def.ERROR("Wiki Format Warning: The string \"%s\" contains \"%s\" but not \"%s\". Removing \"%s\", but please check the format in Wiki document." % (target_text, makeup_begin, makeup_end, makeup_begin), "%s in %s" % (__name__, os.path.basename(__file__)))
			target_text = target_text.replace(makeup_begin, "", 1)
		else: # assert (item_tail > -1)
			makeup_token = target_text[item_head:item_tail+len(makeup_end)]
			display_item = makeup_token.replace(makeup_begin, "")
			display_item = display_item.replace(makeup_end, "")
			if display_item.find(makeup_separator) != -1:
				item_tokens = display_item.split(makeup_separator)
				assert (len(item_tokens) == 2)
				display_item = item_tokens[1] # 2nd one should be display text
			print("### Found a wiki makeup token \"%s\". Changed to \"%s\"" % (makeup_token, display_item))
			target_text = target_text.replace(makeup_token, display_item, 1)

		# Try to find the next
		item_head = target_text.find(makeup_begin)

	target_text = target_text.replace("\'\'", "")   # Might not be necessary for DokuWiki (only for MoinMoinWiki?)
	target_text = target_text.replace("\'\'\'", "") # Might not be necessary for DokuWiki (only for MoinMoinWiki?)
	target_text = target_text.replace("**", "")     # Bold
	target_text = target_text.replace("//", "")     # Italic
	target_text = target_text.replace("%%", "")     # Literal
	target_text = target_text.replace("<", "&lt;")  # rich text tag
	target_text = target_text.replace(">", "&gt;")  # rich text tag

	return target_text

# ----------------------------------------------------------------------------------------
def construct_token_list_from_DokuWiki(sxcmd_config):

	print("Start parsing DokuWiki document (%s as %s %s command) " % (sxcmd_config.wiki, sxcmd_config.category, sxcmd_config.role))

	if sxcmd_config.format != "DokuWiki": sp_global_def.ERROR("Logical Error: Incorrect Wiki format %s! Check the sxcmd_config setting in this script." % (sxcmd_config.format), "%s in %s" % (__name__, os.path.basename(__file__)))

	# Allocate memory for new SXcmd instance
	sxcmd = sxgui_template.SXcmd(sxcmd_config.category, sxcmd_config.role, sxcmd_config.is_submittable)
	usage_token_list = []

	# Define dictionary of keywords:
	# The dictionary maps command token to special data types
	keyword_dict = construct_keyword_dict()

	# Define list of target sections for GUI and set current
	section_lists = []
	section_lists.append("===== Usage ====="); section_usage = len(section_lists) - 1;
	section_lists.append("===== Typical usage ====="); section_typical = len(section_lists) - 1;
	section_lists.append("===== Input ====="); section_input = len(section_lists) - 1;
	current_section = section_usage

	# Define list of subsections of input section and set current
	group_lists = []
	group_lists.append("=== Main Parameters ==="); group_main = len(group_lists) - 1;
	group_lists.append("=== Advanced Parameters ==="); group_advanced = len(group_lists) - 1;
	current_group_name = group_lists[group_main].replace("===", "").replace("Parameters", "").strip().lower()

	# Define States and set current
	state_searching_header  = 0
	state_processing_header  = 1
	state_searching  = 2
	state_processing = 3
	state_done = 4
	current_state = state_searching_header # Assuming the first line starts from ====== COMMAND_NAME =======

	# NOTE: 2015/11/11 Toshio Moriya
	# This should be exception. Need to decide if this should be skipped or exit system.
	if os.path.exists(sxcmd_config.wiki) == False: sp_global_def.ERROR("Runtime Error: Wiki document inot found.", "%s in %s" % (__name__, os.path.basename(__file__)))

	file_wiki = open(sxcmd_config.wiki,'r')

	# Loop through all lines in the wiki document file
	for line_wiki in file_wiki:
		# make sure spaces & new line are not included at head and tail of this line
		line_wiki = line_wiki.strip()

		if not line_wiki:
			# This is empty line. Always ignore it regardless of state
			continue
		if line_wiki == "\\":
			# This is a line contains only a DokuWiki linebreak. Always ignore it regardless of state
			continue
		if line_wiki.find("~~NOTOC~~") != -1:
			# This is DokuWiki control macro in the fist line of SPHIRE Wiki Document
			continue

		if current_state == state_searching_header:
			# Extract the name of sxscript
			target_operator = "====="
			if line_wiki[:len(target_operator)] != target_operator or line_wiki[-len(target_operator):] != target_operator:
#				print("MRK_DEBUG: line_wiki                         := %s" % (line_wiki))
#				print("MRK_DEBUG: target_operator                   := %s" % (target_operator))
#				print("MRK_DEBUG: line_wiki[:len(target_operator)]  := %s" % (line_wiki[:len(target_operator)]))
#				print("MRK_DEBUG: line_wiki[-len(target_operator):] := %s" % (line_wiki[-len(target_operator):]))
				sp_global_def.ERROR("Wiki Format Error: The Wiki document must start from header section title defined by ===== COMMAND_NAME =====.", "%s in %s" % (__name__, os.path.basename(__file__)))
			sxcmd_string = line_wiki.replace(target_operator, "").strip()
			sxcmd_string_token_list = sxcmd_string.split()
			n_sxcmd_string_token_list = len(sxcmd_string_token_list)
			assert(len(sxcmd_string_token_list) > 0)
			sxcmd.name = sxcmd_string_token_list[0]
			if n_sxcmd_string_token_list == 2:
				sxcmd.subname = sxcmd_string_token_list[1]
			elif n_sxcmd_string_token_list > 2:
				sp_global_def.ERROR("Wiki Format Error: Command string should be only a script file name or a script file name plus a subcommand", "%s in %s" % (__name__, os.path.basename(__file__)))
			current_state = state_processing_header
		elif current_state == state_processing_header:
			# Extract the label of this sxscript
			line_buffer = line_wiki
			target_operator = ":"
			item_tail = line_buffer.find(target_operator)
			if item_tail == -1: sp_global_def.ERROR("Wiki Format Error: Header section body should contain a label ending with ':': %s" % line_wiki, "%s in %s" % (__name__, os.path.basename(__file__)))
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
			if current_state != state_processing: sp_global_def.ERROR("Logical Error: This condition should not happen! State setting must be incorrect.", "%s in %s" % (__name__, os.path.basename(__file__)))
			target_operator = "====="
			if line_wiki[:len(target_operator)] == target_operator: # Assuming the section always starts with "======"
#				print("MRK_DEBUG: line_wiki                          := %s" % (line_wiki))
#				print("MRK_DEBUG: Old current_section                := %s" % (current_section))
#				print("MRK_DEBUG: Old section_lists[current_section] := %s" % (section_lists[current_section]))
				# Reached the next section (might be not target)
				current_section += 1 # Update current target section
				if current_section == len(section_lists):
					# All target sections are handled
					current_state = state_done
					break
				current_group_name = group_lists[group_main].replace("===", "").replace("Parameters", "").strip().lower()    # reset group (subsection) for '====== Input ======' and '====== Output ======' sections
#				print("MRK_DEBUG: New current_section                := %s" % (current_section))
#				print("MRK_DEBUG: New section_lists[current_section] := %s" % (section_lists[current_section]))
				if line_wiki.find(section_lists[current_section]) != -1:
					# Found the current target section
					if current_section >= len(section_lists): sp_global_def.ERROR("Logical Error: This condition should not happen! Section setting must be incorrect.", "%s in %s" % (__name__, os.path.basename(__file__)))
					current_state = state_processing
				else:
					# This section was not the current target. Go back to searching state
					current_state = state_searching
			else:
				# We are in a target section
				if current_section == section_usage:
					# Extract 'usage in command line' to identify each command token is either an argument (no-prefix) or option ('--' prefix)
					# This information is also used to check consistency between 'usage in command line' and list in '== Input ==' and '== Output ==' sections
					head_token_idx = 1
					if line_wiki[0:len("sp")] == "sp" or line_wiki[0:len("e2")] == "e2":
						usage_token_list_line = line_wiki.split()
						if usage_token_list_line[0] != sxcmd.name + ".py": sp_global_def.ERROR("Wiki Format Error: First token should be script name with .py (sx*.py or e2*.py)", "%s in %s" % (__name__, os.path.basename(__file__)))
						if sxcmd.subname != "":
							head_token_idx = 2
							if usage_token_list_line[1] != sxcmd.subname: sp_global_def.ERROR("Wiki Format Error: Second token of this command should be subname", "%s in %s" % (__name__, os.path.basename(__file__)))
						# Register arguments and options
						for usage_token in usage_token_list_line[head_token_idx:]:
							# Allocate memory for new command token
							token = sxgui_template.SXcmd_token()
							# Extract key of command token.
							key = usage_token.split("=")[0] # Remove characters after '=' if token contains it (i.e. some options)
							if key == "--MPI":
								print("### Found a --MPI flag \"%s\". Ignoring this key..." % (key))
								continue
							token.key_base = key.strip("-") # Get key base name by removing prefix ('--' or '-' for option)
							token.key_prefix = key[0:len(key) - len(token.key_base)]
							# Register this command token to the usage token list (ordered)
							usage_token_list.append(token)
					# else: Ignore this line (must be comments).
				elif current_section == section_typical:
					target_operator = "mpirun"
					if sxcmd.mpi_support == False and line_wiki.find(target_operator) > -1 and line_wiki.find(sxcmd.name + ".py") > -1:
						sxcmd.mpi_support = True
					# else: Ignore this line
				elif current_section == section_input:
					if line_wiki.find(group_lists[group_advanced]) > -1:
						# Reached the option subsection (argument subsection is done)
						current_group_name = group_lists[group_advanced].replace("===", "").replace("Parameters", "").strip().lower()
					elif line_wiki[0] == ";":
						# Option entry must start with "; "
						line_buffer = line_wiki[1:]  # remove "; "
						# Extract key base name of command token
						target_operator = ":"
						item_tail = line_buffer.find(target_operator)
						if item_tail == -1:
							# ERROR("Warning: This line (%s) is missing key base name (maybe comment line?). Ignoring this line...'."  % (line_wiki),"%s in %s" % (__name__, os.path.basename(__file__)), action = 0)
							continue
						key = line_buffer[0:item_tail].strip()
						key = key.replace("%%", "")
						# Check if --MPI is used in this script
						# NOTE: 2015/11/12 Toshio Moriya
						# The following can be removed when --MPI flag is removed from all sx*.py scripts
						if key == "--MPI":
							# ERROR("Warning: This line (%s) contains MPI flag. The flag will be removed in near future, so ignoring this line...'."  % (line_wiki), "%s in %s" % (__name__, os.path.basename(__file__)), action = 0)
							if sxcmd.mpi_support == False: sp_global_def.ERROR("Logical Error: Since MPI flag is found, the command should support MPI.", "%s in %s" % (__name__, os.path.basename(__file__)))
							sxcmd.mpi_add_flag = True
							continue
						# Allocate memory for new command token
						token = sxgui_template.SXcmd_token()
						# Extract key of command token.
						token.key_base = key.strip("-") # Get key base name by removing prefix ('--' or '-' for option)
						token.key_prefix = key[0:len(key) - len(token.key_base)]
						# Try to set the special type base on the keyword dictionary
						best_keyword_map = SXkeyword_map(99, "")
						for keyword in list(keyword_dict.keys()):
							if key.find(keyword) != -1:
								# command token contains keyword
								keyword_map = keyword_dict[keyword]
								if best_keyword_map.priority > keyword_map.priority:
									# Update keyword_map to one with a higher priority
									best_keyword_map = keyword_map
						token.type = best_keyword_map.token_type # If command token does not contains any keywords, its type stays with ""
						line_buffer = line_buffer[item_tail + len(target_operator):].strip() # Get the rest of line
###						token.is_in_io = True # Set flag to tell this command token is found in input section
						token.group = current_group_name # Set group of command token according to the current subsection
						# Extract label of command token
						target_operator = ":"
						item_tail = line_buffer.find(target_operator)
						if item_tail == -1: sp_global_def.ERROR("Wiki Format Error: This line (%s) is missing label. Please check the format in Wiki document." % line_wiki, "%s in %s" % (__name__, os.path.basename(__file__)))
						token.label = line_buffer[0:item_tail].strip()
						line_buffer = line_buffer[item_tail + len(target_operator):].strip() # Get the rest of line
						# Extract help of command token before default value
						target_operator = "(default"
						item_tail = line_buffer.find(target_operator)
						if item_tail == -1: sp_global_def.ERROR("Wiki Format Error: This line (%s) is missing default setting. Please check the format in Wiki document." % line_wiki, "%s in %s" % (__name__, os.path.basename(__file__)))
						token.help = remove_DokuWiki_makeup(line_buffer[0:item_tail])
						line_buffer = line_buffer[item_tail + len(target_operator):].strip() # Get the rest of line
						# Extract default value of command token
						target_operator = ")"
						item_tail = line_buffer.find(target_operator)
						if item_tail == -1: sp_global_def.ERROR("Wiki Format Error: This line (%s) is missing ')' for default setting. Please check the format in Wiki document." % line_wiki,"%s in %s" % (__name__, os.path.basename(__file__)))
						default_value = line_buffer[0:item_tail].strip() # make sure spaces & new line are not included at head and tail
						if default_value.find("required") != -1:
							# This is a required command token and should have value type instead of default value
							token.is_required = True
							token.default = ""
							if not token.type:
								# Type is still empty, meaning no special type is assigned
								# Extract the data type (the rest of line)
								token.type = default_value.replace("required", "").strip()
							assert (token.is_required)
							assert (not token.is_locked)
							assert (not token.is_reversed)
						elif default_value.find("question reversed in GUI") != -1:
							# token.is_required = False
							token.default = default_value.replace("question reversed in GUI", "").strip()
							token.type = "bool"
							if token.default == "True":
								token.default = False # convert the default value to opposite boolean
							elif token.default == "False":
								token.default = True # convert the default value to opposite boolean
							else:
								assert (False)
							assert (not token.is_required)
							assert (not token.is_locked)
							assert (not token.is_reversed)
						elif default_value.find("value reversed in GUI") != -1:
							# token.is_required = False
							token.is_reversed = True
							token.default = default_value.replace("value reversed in GUI", "").strip()
							token.type = "bool"
							if token.default == "True":
								token.default = False # convert the default value to opposite boolean
							elif token.default == "False":
								token.default = True # convert the default value to opposite boolean
							else:
								assert (False)
							assert (not token.is_required)
							assert (not token.is_locked)
							assert (token.is_reversed)
						else:
							# This is not required command token and should have default value
							# token.is_required = False
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
							assert (not token.is_required)
							assert (not token.is_locked)
							assert (not token.is_reversed)
						if not token.is_locked:
							token.restore = token.default
						target_operator = ":"
						item_tail = line_buffer.find(target_operator)
						if item_tail != -1:
							line_buffer = line_buffer[item_tail + len(target_operator):].strip() # Get the rest of line
							target_operator = ":"
							item_tail = line_buffer.find(target_operator)
							if item_tail == -1:
								entry = line_buffer.replace("%%", "").strip()
							else:
								entry = line_buffer[:item_tail].replace("%%", "").strip()

							for dependency in entry.split():
								try:
									group_key, state = dependency.split('==')
								except:
									try:
										group_key, state = dependency.split('!=')
									except Exception as e:
										print(dependency)
										raise
									inverse = True
								else:
									inverse = False
								token.dependency_group.append([group_key.strip('-'), state, inverse])
								try:
									sxcmd.dependency_dict[group_key].append([token.key_base, state, inverse])
								except KeyError:
									sxcmd.dependency_dict[group_key] = [[token.key_base, state, inverse]]
						# Initialise restore value with default value
						# Ignore the rest of line ...
						# Register this command token to the list (ordered) and dictionary (unordered)
						sxcmd.token_list.append(token)
						sxcmd.token_dict[token.key_base] = token
					# else:
						# This is not option entry. Ignore this line
				else:
					sp_global_def.ERROR("Logical Error: This section is invalid. Did you assigne an invalid section?", "%s in %s" % (__name__, os.path.basename(__file__)))

	if current_state != state_done: sp_global_def.ERROR("Wiki Format Error: parser could not extract all information necessary (current state %d). Please check if the Wiki format has all required sections." % (current_state), "%s in %s" % (__name__, os.path.basename(__file__)))

	for sxcmd_token_key_base in list(sxcmd.token_dict.keys()):
		# Make sure there are no extra arguments or options in command token dictionary compared with command token list.
		is_found = False
		for sxcmd_token in sxcmd.token_list:
			if sxcmd_token.key_base == sxcmd_token_key_base:
				is_found = True
				break
		if not is_found: sp_global_def.ERROR("Logical Error: Registered key base name in the command token dictionary is not found in the command token list.", "%s in %s" % (__name__, os.path.basename(__file__)))
		# Make sure there are no extra arguments or options in command token list compared with usage token list extracted from "usage in command line" of "====== Usage ======".
		is_found = False
		for usage_token in usage_token_list:
			if usage_token.key_base == sxcmd_token_key_base:
				is_found = True
				break
			elif sxcmd.token_dict[sxcmd_token_key_base].type == 'bool_ignore':
				is_found = True
				break
		if not is_found: sp_global_def.ERROR("Wiki Format Error: An extra argument or option (%s) is found in the command token dictionary extracted from '===== Input =====' compared with 'usage in command line' of '====== Usage ======'." % sxcmd_token_key_base, "%s in %s" % (__name__, os.path.basename(__file__)))

	for sxcmd_token in sxcmd.token_list:
		# Make sure there are no extra arguments or options in command token list compared with command token dictionary.
		if sxcmd_token.key_base not in list(sxcmd.token_dict.keys()): sp_global_def.ERROR("Logical Error: Registered key base name in the command token list is not registered as key base name in the command token dictionary.", "%s in %s" % (__name__, os.path.basename(__file__)))
		# Make sure there are no extra arguments or options in command token list compared with usage token list extracted from "usage in command line" of "====== Usage ======".
		is_found = False
		for usage_token in usage_token_list:
			if usage_token.key_base == sxcmd_token.key_base:
				is_found = True
				break
			elif sxcmd_token.type == 'bool_ignore':
				is_found = True
				break
		if not is_found: sp_global_def.ERROR("Wiki Format Error: An extra argument or option (%s) is found in the command token list extracted from '===== Input =====' compared with 'usage in command line' of '====== Usage ======'." % sxcmd_token.key_base, "%s in %s" % (__name__, os.path.basename(__file__)))

	for usage_token in usage_token_list:
		# Make sure there are no extra arguments or options in usage token list extracted from "usage in command line" of "====== Usage ======" compared with command token dictionary. 
		if usage_token.key_base not in list(sxcmd.token_dict.keys()): sp_global_def.ERROR("Wiki Format Error: An extra argument or option (%s) is found in 'usage in command line' of '====== Usage ======' compared with the command token dictionary extracted from '===== Input =====' ." % usage_token.key_base, "%s in %s" % (__name__, os.path.basename(__file__)))
		# Make sure there are no extra arguments or options in usage token list extracted from "usage in command line" of "====== Usage ======" compared with command token list. 
		is_found = False
		for sxcmd_token in sxcmd.token_list:
			if sxcmd_token.key_base == sxcmd_token_key_base:
				is_found = True
				break
			elif sxcmd_token.type == 'bool_ignore':
				is_found = True
				break
		if not is_found: sp_global_def.ERROR("Wiki Format Error: An additional argument or option (%s) is found in 'usage in command line' of '====== Usage ======' compared with the command token list extracted from '===== Input =====' ." % usage_token.key_base, "%s in %s" % (__name__, os.path.basename(__file__)))

	file_wiki.close()

	handle_exceptional_cases(sxcmd)

	print("Succeed to parse DokuWiki document (%s as %s %s command)" % (sxcmd_config.wiki, sxcmd_config.category, sxcmd_config.role))

	"""
	# For DEBUG
	if sxcmd.name in ["sxrelion2sphire", "sxpipe"]:
		print("><><>< DEBUG OUTPUT ><><><")
		print("")
		print("------")
		print("GLOBAL")
		print("------")
		print("name            : \'{}\'".format(sxcmd.name))
		print("subname         : \'{}\'".format(sxcmd.subname))
		print("mode            : \'{}\'".format(sxcmd.mode))
		print("subset_config   : \'{}\'".format(sxcmd.subset_config))
		print("label           : \'{}\'".format(sxcmd.label))
		print("short_info      : \'{}\'".format(sxcmd.short_info))
		print("mpi_support     : \'{}\'".format(sxcmd.mpi_support))
		print("mpi_add_flag    : \'{}\'".format(sxcmd.mpi_add_flag))
		print("category        : \'{}\'".format(sxcmd.category))
		print("role            : \'{}\'".format(sxcmd.role))
		print("is_submittable  : \'{}\'".format(sxcmd.is_submittable))
		print("len(token_list) : {}".format(len(sxcmd.token_list)))
		print("len(token_dict) : {}".format(len(sxcmd.token_dict)))
		print("")
		print("--------------")
		print("cmd_token_list")
		print("--------------")
		for token in sxcmd.token_list:
			print("\'{}{}\' (group=\'{}\', is_required=\'{}\', is_locked=\'{}\', is_reversed=\'{}\', default=\'{}\', restore=\'{}\', type=\'{}\', label=\'{}\', help=\'{}\'".format(token.key_prefix, token.key_base, token.group, token.is_required, token.is_locked, token.is_reversed, token.default, token.restore, token.type, token.label, token.help))
		print("")
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
	if mode_token_edit.key_base not in list(fullset_token_dict.keys()): sp_global_def.ERROR("Logical Error: This condition should not happen! Subset command configuration must be incorrect. Key (%s) should exists." % (mode_token_edit.key_base), "%s in %s" % (__name__, os.path.basename(__file__)))
	mode_token = fullset_token_dict[mode_token_edit.key_base]

	# Create mode name of this subset, append key base of mode token to mode_name of this command
	if not sxsubcmd_config.is_modeless:
		sxcmd.mode = mode_token.key_base
	else:
		assert sxsubcmd_config.is_modeless, "MRK_DEBUG: Must be always True"
		assert sxcmd.mode == "", "MRK_DEBUG: This subcommand has no mode! sxcmd.mode must be empty string!"
	# print "MRK_DEBUG: sxcmd.mode = %s" % (sxcmd.mode)
	
	# Set unique name of this subset configuration
	sxcmd.subset_config = sxsubcmd_config.subset_config
	
	# Set command label of this subset
	sxcmd.label = sxsubcmd_config.label
	# print "MRK_DEBUG: sxcmd.label = %s" % (sxcmd.label)

	# Use label of mode token as a short info of subset command
	if sxsubcmd_config.short_info is not None:
		sxcmd.short_info = sxsubcmd_config.short_info
	else:
		sxcmd.short_info = "%s" % (mode_token.help)
		# sxcmd.short_info = "%s. %s" % (mode_token.label, mode_token.help)
	# print "MRK_DEBUG: sxcmd.short_info = %s" % (sxcmd.short_info)

	# Set command mpi support of this subset if necessary
	if sxsubcmd_config.mpi_support != None:
		sxcmd.mpi_support = sxsubcmd_config.mpi_support
		if sxcmd.mpi_support == False:
			sxcmd.mpi_add_flag = False
	# print "MRK_DEBUG: sxcmd.mpi_support = %s" % (sxcmd.mpi_support)

	# Reconstruct token list
	for token_edit in sxsubcmd_config.token_edit_list:
		# print "MRK_DEBUG: token_edit.key_base = %s" % (token_edit.key_base)
		if token_edit.key_base is None: sp_global_def.ERROR("Logical Error: This condition should not happen! Subset command configuration must be incorrect. Invalid None Key (%s)." % (token_edit.key_base) , "%s in %s" % (__name__, os.path.basename(__file__)))
		if token_edit.key_base == "": sp_global_def.ERROR("Logical Error: This condition should not happen! Subset command configuration must be incorrect. Invalid empty string Key (%s)." % (token_edit.key_base) , "%s in %s" % (__name__, os.path.basename(__file__)))
		
		token = None
		if token_edit.key_base not in list(fullset_token_dict.keys()):
			# token key base is not found in fullset. This must be an argument to be added
			if token_edit.key_prefix is None: sp_global_def.ERROR("Logical Error: This condition should not happen! Subset command configuration must be incorrect. Prefix (%s) for Key (%s) should NOT be None." % (token_edit.key_prefix, token_edit.key_base) , "%s in %s" % (__name__, os.path.basename(__file__)))
			if token_edit.key_prefix != "": sp_global_def.ERROR("Logical Error: This condition should not happen! Subset command configuration must be incorrect. Key (%s) should be argument (Prefix (%s) should be empty string)." % (token_edit.key_base, token_edit.key_prefix) , "%s in %s" % (__name__, os.path.basename(__file__)))
			# token = token_edit
			token = sxgui_template.SXcmd_token()
			token.key_base = token_edit.key_base
		else:
			# token key base is found in fullset. This must be an option.
			token = fullset_token_dict[token_edit.key_base]
		assert(token is not None)
		
		if token_edit.key_prefix is not None:
			token.key_prefix = token_edit.key_prefix
		if token_edit.label is not None:
			token.label = token_edit.label
		if token_edit.help is not None:
			token.help = token_edit.help
		if token_edit.group is not None:
			token.group = token_edit.group
		if token_edit.is_required is not None:
			token.is_required = token_edit.is_required
		if token_edit.is_locked is not None:
			token.is_locked = token_edit.is_locked
		if token_edit.is_reversed is not None:
			token.is_reversed = token_edit.is_reversed
		if token_edit.default is not None:
			token.default = token_edit.default
		if token_edit.restore is not None:
			token.restore = token_edit.restore
		if token_edit.type is not None:
			token.type = token_edit.type
		
		if not token.is_locked:
			token.restore = token.default
		
		# Make sure all fields of token are not None
		if token.key_base is None: sp_global_def.ERROR("Logical Error: This condition should not happen! Subset command configuration must be incorrect. token.key_base should NOT be None.", "%s in %s" % (__name__, os.path.basename(__file__)))
		if token.key_prefix is None: sp_global_def.ERROR("Logical Error: This condition should not happen! Subset command configuration must be incorrect. token.key_prefix should NOT be None.", "%s in %s" % (__name__, os.path.basename(__file__)))
		if token.label is None: sp_global_def.ERROR("Logical Error: This condition should not happen! Subset command configuration must be incorrect. token.label should NOT be None.", "%s in %s" % (__name__, os.path.basename(__file__)))
		if token.help is None: sp_global_def.ERROR("Logical Error: This condition should not happen! Subset command configuration must be incorrect. token.help should NOT be None.", "%s in %s" % (__name__, os.path.basename(__file__)))
		if token.group is None: sp_global_def.ERROR("Logical Error: This condition should not happen! Subset command configuration must be incorrect. token.group should NOT be None.", "%s in %s" % (__name__, os.path.basename(__file__)))
		if token.is_required is None: sp_global_def.ERROR("Logical Error: This condition should not happen! Subset command configuration must be incorrect. token.is_required should NOT be None.", "%s in %s" % (__name__, os.path.basename(__file__)))
		if token.is_locked is None: sp_global_def.ERROR("Logical Error: This condition should not happen! Subset command configuration must be incorrect. token.is_locked should NOT be None.", "%s in %s" % (__name__, os.path.basename(__file__)))
		if token.is_reversed is None: sp_global_def.ERROR("Logical Error: This condition should not happen! Subset command configuration must be incorrect. token.is_reversed should NOT be None.", "%s in %s" % (__name__, os.path.basename(__file__)))
		if token.default is None: sp_global_def.ERROR("Logical Error: This condition should not happen! Subset command configuration must be incorrect. token.default should NOT be None.", "%s in %s" % (__name__, os.path.basename(__file__)))
		if token.restore is None: sp_global_def.ERROR("Logical Error: This condition should not happen! Subset command configuration must be incorrect. token.restore should NOT be None.", "%s in %s" % (__name__, os.path.basename(__file__)))
		if token.type is None: sp_global_def.ERROR("Logical Error: This condition should not happen! Subset command configuration must be incorrect. token.type should NOT be None.", "%s in %s" % (__name__, os.path.basename(__file__)))

		sxcmd.token_list.append(token)
		assert token_edit.key_base not in sxcmd.token_dict, token_edit.key_base
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
	output_file.write("; %s.subname = \"%s\"" % (sxcmd_variable_name, sxcmd.subname))
	output_file.write("; %s.mode = \"%s\"" % (sxcmd_variable_name, sxcmd.mode))
	output_file.write("; %s.subset_config = \"%s\"" % (sxcmd_variable_name, sxcmd.subset_config))
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
		output_file.write("; token.is_locked = %s" % token.is_locked)
		output_file.write("; token.is_reversed = %s" % token.is_reversed)
		output_file.write("; token.dependency_group = %s" % str([[str(entry) for entry in group_dep] for group_dep in token.dependency_group]))
		if token.type in ("bool", "bool_ignore"):
			output_file.write("; token.default = %s" % token.default)
			output_file.write("; token.restore = %s" % token.restore)
		else:
			if token.is_required:
				if token.is_locked:
					output_file.write("; token.default = \"\"")
					output_file.write("; token.restore = \"%s\"" % token.restore)
				else:
					output_file.write("; token.default = \"\"")
					output_file.write("; token.restore = \"\"")
			else:
				output_file.write("; token.default = \"%s\"" % token.default)
				output_file.write("; token.restore = \"%s\"" % token.restore)
		output_file.write("; token.type = \"%s\"" % token.type)
		# output_file.write("; token.is_in_io = %s" % token.is_in_io)

		output_file.write("; %s.token_list.append(token)" % sxcmd_variable_name)
		output_file.write("; %s.token_dict[token.key_base] = token" % sxcmd_variable_name)
		for entry in token.dependency_group:
			if entry[0]:
				output_file.write("; %s.dependency_dict.setdefault('%s', []).append([token.key_base, '%s', '%s'])" % (sxcmd_variable_name, entry[0], entry[1], entry[2]))
		output_file.write("\n")

	output_file.write("\n")
	output_file.write("\t\t%s_list.append(%s)\n" % (sxcmd_variable_name, sxcmd_variable_name))
	output_file.write("\n")

	return

# ========================================================================================
def create_sxcmd_subconfig_window_makevstack():
	token_edit_list = []
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("makevstack"); token_edit.is_required = True; token_edit.is_locked = False; token_edit.default = "none"; token_edit.restore = "none"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("input_bdb_stack_pattern"); token_edit.key_prefix = ""; token_edit.label = "Input BDB image stack pattern"; token_edit.help = "Specify file path pattern of stack subsets created in particle extraction using a wild card /'*/' (e.g. /'//sxwindow_output_dir//*/'). The stack subsets are located in the sxwindow output directory."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "dir_list"; token_edit_list.append(token_edit)

	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("Particle Stack", None, token_edit_list, sxsubcmd_mpi_support, subset_config="fullset")

	return sxcmd_subconfig

### def create_sxcmd_subconfig_isacselect():
### 	token_edit_list = []
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("isacselect"); token_edit.is_required = True; token_edit.is_locked = True; token_edit.default = True; token_edit.restore = True; token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("class_file_name"); token_edit.key_prefix = ""; token_edit.label = "ISAC2 class file name"; token_edit.help = "File name of the class averages. It is located in the ISAC2 output directory."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "data2d_one"; token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("output_list"); token_edit.key_prefix = ""; token_edit.label = "Output ISAC2 particle ID list"; token_edit.help = "Output text file containing retrieved member particle IDs of all ISAC2 classes."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "output"; token_edit_list.append(token_edit)
### 
### 	sxsubcmd_mpi_support = False
### 	sxcmd_subconfig = SXsubcmd_config("Get ISAC2 Particles", None, token_edit_list, sxsubcmd_mpi_support)
### 
### 	return sxcmd_subconfig

### def create_sxcmd_subconfig_isac_makevstack():
### 	token_edit_list = []
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("makevstack"); token_edit.is_required = True; token_edit.is_locked = False; token_edit.default = "none"; token_edit.restore = "none"; token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("input_bdb_stack_file"); token_edit.key_prefix = ""; token_edit.label = "Input BDB image stack"; token_edit.help = "Specify path to input BDB stack file which is used for the input of the associated ISAC2 run."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "bdb2d_stack"; token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("list"); token_edit_list.append(token_edit)
### 
### 	sxsubcmd_mpi_support = False
### 	sxcmd_subconfig = SXsubcmd_config("Create Stack Subset", None, token_edit_list, sxsubcmd_mpi_support)
### 
### 	return sxcmd_subconfig

def add_sxcmd_subconfig_isac_beautifier_shared(token_edit_list):
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("stack"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("isac_dir"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("output_dir"); token_edit_list.append(token_edit)
	
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("pixel_size"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("radius"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("noctf"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("skip_local_alignment"); token_edit_list.append(token_edit)

	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("fl"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("xr"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("ts"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("fh"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("maxit"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("navg"); token_edit_list.append(token_edit)

# def create_sxcmd_subconfig_isac_beautifier_to_bfactor():
# 	token_edit_list = []
# 	token_edit = SXcmd_token(); token_edit.initialize_edit("B_enhance"); token_edit.is_required = True; token_edit.is_locked = True; token_edit.default = True; token_edit.restore = True; token_edit_list.append(token_edit)
# 	
# 	token_edit = SXcmd_token(); token_edit.initialize_edit("B_start"); token_edit.group = "main"; token_edit_list.append(token_edit)
# 	token_edit = SXcmd_token(); token_edit.initialize_edit("Bfactor"); token_edit.group = "main"; token_edit_list.append(token_edit)
# 
# 	add_sxcmd_subconfig_isac_beautifier_shared(token_edit_list)
# 
# 	sxsubcmd_mpi_support = True
# 	sxcmd_subconfig = SXsubcmd_config("Beautifier - Adjust to B-factor", "Beautify the ISAC2 2D clustering result with the original pixel size. In addition, adjust the power spectrum of resultant average images using B-factor to enhance averages.", token_edit_list, sxsubcmd_mpi_support)
# 
# 	return sxcmd_subconfig


# def create_sxcmd_subconfig_isac_beautifier_to_rot_avg():
# 	token_edit_list = []
# 	token_edit = SXcmd_token(); token_edit.initialize_edit("adjust_to_given_pw2"); token_edit.is_required = True; token_edit.is_locked = True; token_edit.default = True; token_edit.restore = True; token_edit_list.append(token_edit)
# 	
# 	token_edit = SXcmd_token(); token_edit.initialize_edit("modelpw"); token_edit.group = "main"; token_edit_list.append(token_edit)
# 
# 	add_sxcmd_subconfig_isac_beautifier_shared(token_edit_list)
# 
# 	sxsubcmd_mpi_support = True
# 	sxcmd_subconfig = SXsubcmd_config("Beautifier - Adjust to Rot. Avgs.", "Beautify the ISAC2 2D clustering result with the original pixel size. In addition, adjust the power spectrum of resultant average images to the user-provided 1-D reference power spectrum.", token_edit_list, sxsubcmd_mpi_support)
# 
# 	return sxcmd_subconfig

# def create_sxcmd_subconfig_isac_beautifier_to_model():
# 	token_edit_list = []
# 	token_edit = SXcmd_token(); token_edit.initialize_edit("adjust_to_analytic_model"); token_edit.is_required = True; token_edit.is_locked = True; token_edit.default = True; token_edit.restore = True; token_edit_list.append(token_edit)
# 
# 	add_sxcmd_subconfig_isac_beautifier_shared(token_edit_list)
# 
# 	sxsubcmd_mpi_support = True
# 	sxcmd_subconfig = SXsubcmd_config("Beautifier - Adjust to Model", "Beautify the ISAC2 2D clustering result with the original pixel size. In addition, adjust the power spectrum of resultant average images to an analytical model.", token_edit_list, sxsubcmd_mpi_support)
# 
# 	return sxcmd_subconfig

# def create_sxcmd_subconfig_isac_beautifier_no_adjust():
# 	token_edit_list = []
# 	token_edit = SXcmd_token(); token_edit.initialize_edit("no_adjustment"); token_edit.is_required = True; token_edit.is_locked = True; token_edit.default = True; token_edit.restore = True; token_edit_list.append(token_edit)
# 	
# 	add_sxcmd_subconfig_isac_beautifier_shared(token_edit_list)
# 
# 	sxsubcmd_mpi_support = True
# 	sxcmd_subconfig = SXsubcmd_config("Beautifier - No Adjust", "Beautify the ISAC2 2D clustering result using the original pixel size, without adjusting the power spectrum of resultant average images. Use this option to skip all power spectrum adjustment methods and simply compute class averages with the original pixel size.", token_edit_list, sxsubcmd_mpi_support)
# 
# 	return sxcmd_subconfig

def create_sxcmd_subconfig_viper_changesize():
	token_edit_list = []
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("changesize"); token_edit.label = "Change size of VIPER model"; token_edit.is_required = True; token_edit.is_locked = True; token_edit.default = True; token_edit.restore = True; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("input_stack"); token_edit.key_prefix = ""; token_edit.label = "Input Viper Model"; token_edit.help = "Input Viper Model."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "data3d_one"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("output_stack"); token_edit.key_prefix = ""; token_edit.label = "Output Resized Viper Model"; token_edit.help = "Output resized (decimated or interpolated up) Viper Model."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "output"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("ratio"); token_edit_list.append(token_edit)

	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("Change Size of VIPER Model", None, token_edit_list, sxsubcmd_mpi_support)

	return sxcmd_subconfig

def create_sxcmd_subconfig_viper_window():
	token_edit_list = []
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("clip"); token_edit.label = "Window to specified size [Pixels]"; token_edit.is_required = True; token_edit.is_locked = False; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("input_volume"); token_edit.key_prefix = ""; token_edit.label = "Input volume"; token_edit.help = "Input volume file name."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = "none"; token_edit.type = "data3d_one"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("output_file"); token_edit.key_prefix = ""; token_edit.label = "Output windowed volume"; token_edit.help = "Output windowed (clipped/padded) volume file name."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = "none"; token_edit.type = "output"; token_edit_list.append(token_edit)

	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("Window VIPER Model", None, token_edit_list, sxsubcmd_mpi_support)

	return sxcmd_subconfig

# def create_sxcmd_subconfig_scale_clip():
# 	token_edit_list = []
# 	token_edit = SXcmd_token(); token_edit.initialize_edit("scale"); token_edit.label = "Resample ratio"; token_edit.help = "Rescale the volume by the specified ratio before padding/clipping."; token_edit.is_required = True; token_edit_list.append(token_edit)
# 	token_edit = SXcmd_token(); token_edit.initialize_edit("clip"); token_edit.label = "Pad/Clip volume to specified size [Pixels]"; token_edit.is_required = True; token_edit_list.append(token_edit)
# 	token_edit = SXcmd_token(); token_edit.initialize_edit("input_volume"); token_edit.key_prefix = ""; token_edit.label = "Input volume"; token_edit.help = "Input volume file name."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = "none"; token_edit.type = "data3d_one"; token_edit_list.append(token_edit)
# 	token_edit = SXcmd_token(); token_edit.initialize_edit("output_file"); token_edit.key_prefix = ""; token_edit.label = "Output resampled volume"; token_edit.help = "Output resampled volume file name."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = "none"; token_edit.type = "output"; token_edit_list.append(token_edit)
#
# 	sxsubcmd_mpi_support = False
# 	sxcmd_subconfig = SXsubcmd_config("Resample", None, token_edit_list, sxsubcmd_mpi_support)
#
# 	return sxcmd_subconfig

def create_sxcmd_subconfig_adaptive_mask3d():
	token_edit_list = []
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("adaptive_mask"); token_edit.is_required = True; token_edit.is_locked = True; token_edit.default = True; token_edit.restore = True; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("input_volume"); token_edit.key_prefix = ""; token_edit.label = "Input volume"; token_edit.help = "Input reference volume"; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "data3d_one"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("output_mask3D"); token_edit.key_prefix = ""; token_edit.label = "Output mask"; token_edit.help = "Output 3D mask"; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "output"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("use_mol_mass"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("mol_mass"); token_edit.group = "advanced";  token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("threshold"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("nsigma"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("ndilation"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("edge_width"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("edge_type"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("Adaptive 3D Mask", None, token_edit_list, sxsubcmd_mpi_support)

	return sxcmd_subconfig

def create_sxcmd_subconfig_balance_angles():
	token_edit_list = []
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("balance_angular_distribution"); token_edit.is_required = True; token_edit.is_locked = True; token_edit.default = True; token_edit.restore = True; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("params_any_txt"); token_edit.key_prefix = ""; token_edit.label = "Projection parameters"; token_edit.help = "Typically from MERIDIEN with a filename in the form of Refine3D/final_params_0##.txt"; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "params_any_txt"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("selection_list"); token_edit.key_prefix = ""; token_edit.label = "Output selection list"; token_edit.help = "Text file with a list of retained particle images"; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "output"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("max_occupy"); token_edit.is_required = True; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("angstep"); token_edit.default = 15; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("symmetry"); token_edit.default = "c1"; token_edit_list.append(token_edit)
	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("Balance Angular Distribution", None, token_edit_list, sxsubcmd_mpi_support)

	return sxcmd_subconfig

def create_sxcmd_subconfig_binary_mask3d():
	token_edit_list = []
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("binary_mask"); token_edit.is_required = True; token_edit.is_locked = True; token_edit.default = True; token_edit.restore = True; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("input_volume"); token_edit.key_prefix = ""; token_edit.label = "Input volume"; token_edit.help = "Input reference volume"; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "data3d_one"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("output_mask3D"); token_edit.key_prefix = ""; token_edit.label = "Output mask"; token_edit.help = "Output 3D mask"; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "output"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("use_mol_mass"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("mol_mass"); token_edit.group = "advanced";  token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("threshold"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("nsigma"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("nerosion"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("ndilation"); token_edit_list.append(token_edit)
	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("Binary 3D Mask", None, token_edit_list, sxsubcmd_mpi_support)

	return sxcmd_subconfig

def create_sxcmd_subconfig_postrefiner_halfset_vol():
	token_edit_list = []
	help_string = "Post-refine a pair of unfiltered odd & even halfset maps by combining them, then enhancing the high frequencies (Halfset Maps Mode). B-factor can be automatically estimated from the unfiltered halfset maps. This mode requires two arguments; use unfiltered halfset maps produced by MERIDIEN."
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("combinemaps"); token_edit.label = "Post-refine halfset volumes"; token_edit.help = help_string; token_edit.is_required = True; token_edit.is_locked = True; token_edit.default = True; token_edit.restore = True; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("firstvolume"); token_edit.key_prefix = ""; token_edit.label = "First unfiltered halfset map"; token_edit.help = "As generated by sxmeridien."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "data3d_one"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("secondvolume"); token_edit.key_prefix = ""; token_edit.label = "Second unfiltered halfset map"; token_edit.help = "As generated by sxmeridien."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "data3d_one"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("output_dir"); token_edit_list.append(token_edit)
	
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("pixel_size"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("mask"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("do_adaptive_mask"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("use_mol_mass"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("mol_mass");  token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("threshold"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("nsigma"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("edge_width"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("edge_type"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("do_approx"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("ndilation"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("mtf"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("fsc_adj"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("B_enhance"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("B_start"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("B_stop"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
#	token_edit = SXcmd_token(); token_edit.initialize_edit("randomphasesafter"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("fl"); token_edit.type = "float"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("aa"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("output"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	
	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("PostRefiner", help_string, token_edit_list, sxsubcmd_mpi_support, subset_config="halfset maps")

	return sxcmd_subconfig

def create_sxcmd_subconfig_postrefiner_single_vols():
	token_edit_list = []
	help_string = "Post-refine a single map by enhancing the high frequencies (Single Map Mode). Only ad-hoc low-pass filter cutoff and B-factor can be used. This mode requires one argument; path pattern with wild card \'*\' to specify a list of  volumes or a path to a volume (without wild card \'*\')."
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("combinemaps"); token_edit.label = "Post-refine single map"; token_edit.help = help_string; token_edit.is_required = True; token_edit.is_locked = True; token_edit.default = True; token_edit.restore = True; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("input_volume_pattern"); token_edit.key_prefix = ""; token_edit.label = "Input volume pattern"; token_edit.help = "Specify path pattern of input single volumes with a wild card \'*\' or path to single volume file (without wild card \'*\'). Use the wild card to indicate the place of variable part of the file names. The path pattern must be enclosed by single quotes (\') or double quotes (\") (Note: sxgui.py automatically adds single quotes (\'))."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "data3d_one"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("output_dir"); token_edit_list.append(token_edit)
	
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("pixel_size"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("mask"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("do_adaptive_mask"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("use_mol_mass"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("mol_mass"); token_edit.group = "advanced";  token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("threshold"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("nsigma"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("ndilation"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("edge_width"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("edge_type"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("do_approx"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("mtf"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("B_enhance"); token_edit.is_required = True; token_edit.help = "Non-zero positive value: program use provided B-factor [A^2] to enhance the map; -1.0: B-factor is not applied."; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("fl"); token_edit.is_required = True; token_edit.help = "A value larger than 0.5: low-pass filter to the resolution in Angstrom; -1.0: no low-pass filter."; token_edit.type = "float"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("aa"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("output"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	
	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("PostRefiner (Single Map)", help_string, token_edit_list, sxsubcmd_mpi_support, subset_config="single map")
	
	return sxcmd_subconfig

### NOTE: Toshio Moriya 2018/02/12
### Because SORT3D cluster volumes is heavily low-pass filtered at this point, 
### we did not see any improvement by PostRefiner.
### Therefore, we moved this from main pipeline to utilities for now.
### 
### def create_sxcmd_subconfig_postrefiner_cluster_vol():
### 	token_edit_list = []
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("combinemaps"); token_edit.label = "Enhance cluster volumes"; token_edit.help = "Enhance the power spectrum of cluster volumes, produced by SORT3D_DEPTH, at high frequencies (Cluster Volumes Mode). Only ad-hoc low-pass filter cutoff and B-factor can be used. This mode requires one argument; path pattern with wild card '*' can be used to specify a list of volumes. It is mainly used with SORT3D_DEPTH outputs."; token_edit.is_required = True; token_edit.is_locked = True; token_edit.default = True; token_edit.restore = True; token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("input_cluster_volume_pattern"); token_edit.key_prefix = ""; token_edit.label = "Input cluster volume pattern"; token_edit.help = "Specify path pattern of input cluster volumes, created by SORT3D_DEPTH, with a wild card (*). Use the wild card to indicate the place of variable part of the file names (typically, cluster ID; e.g. \'outdir_sort3d_depth/vol_cluster*.hdf\'). The path pattern must be enclosed by single quotes (\') or double quotes (\") (Note: sxgui.py automatically adds single quotes (\'))."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "data3d_one"; token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("output_dir"); token_edit_list.append(token_edit)
### 	
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("pixel_size"); token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("mask"); token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("do_adaptive_mask"); token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("mask_threshold"); token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("edge_width"); token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("dilation"); token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("mtf"); token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("B_enhance"); token_edit.is_required = True; token_edit.help = "Non-zero positive value: program use the given value [A^2] to enhance map; -1.0: B-factor is not applied."; token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("fl"); token_edit.is_required = True; token_edit.help = "A value larger than 0.5: low-pass filter to the value in Angstrom; -1.0: no low-pass filter."; token_edit.type = "float"; token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("aa"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("output"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
### 	
### 	sxsubcmd_mpi_support = False
### 	sxcmd_subconfig = SXsubcmd_config("PostRefiner of Cluster Volumes", "Enhance the power spectrum of the cluster volumes, produced by SORT3D_DEPTH, at high frequencies (Cluster Volumes Mode). Only ad-hoc low-pass filter cutoff and B-factor can be used. This mode is mainly used with SORT3D_DEPTH outputs, but also can enhance the power spectrum of ANY volumes.", token_edit_list, sxsubcmd_mpi_support, subset_config="cluster volumes")
### 	
### 	return sxcmd_subconfig
### 
### def create_sxcmd_subconfig_postrefiner_single_vol():
### 	token_edit_list = []
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("combinemaps"); token_edit.label = "Enhance single volume"; token_edit.help = "Enhance the power spectrum of any single volume at high frequencies (Single Volume Mode). Only ad-hoc low-pass filter cutoff and B-factor can be used. This mode requires one argument; path to the single volume."; token_edit.is_required = True; token_edit.is_locked = True; token_edit.default = True; token_edit.restore = True; token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("input_volume"); token_edit.key_prefix = ""; token_edit.label = "Input volume"; token_edit.help = "Specify path to single volume file."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "data3d_one"; token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("output_dir"); token_edit_list.append(token_edit)
### 	
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("pixel_size"); token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("mask"); token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("do_adaptive_mask"); token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("mask_threshold"); token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("edge_width"); token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("dilation"); token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("mtf"); token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("B_enhance"); token_edit.is_required = True; token_edit.help = "Non-zero positive value: program use the given value [A^2] to enhance map; -1.0: B-factor is not applied."; token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("fl"); token_edit.is_required = True; token_edit.help = "A value larger than 0.5: low-pass filter to the value in Angstrom; -1.0: no low-pass filter."; token_edit.type = "float"; token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("aa"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("output"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
### 	
### 	sxsubcmd_mpi_support = False
### 	sxcmd_subconfig = SXsubcmd_config("PostRefiner of Single Volume", "Enhance the power spectrum of any single volume at high frequencies (Single Volume Mode). Only ad-hoc low-pass filter cutoff and B-factor can be used.", token_edit_list, sxsubcmd_mpi_support, subset_config="single volume")
### 	
### 	return sxcmd_subconfig


def create_sxcmd_subconfig_variability_preprocess():
	token_edit_list = []
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("symmetrize"); token_edit.is_required = True; token_edit.is_locked = True; token_edit.default = True; token_edit.restore = True; token_edit_list.append(token_edit)

	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("prj_stack"); token_edit.key_prefix = ""; token_edit.label = "Input image stack"; token_edit.help = "The images must contain the 3D orientation parameters in headers and optionally also CTF information. The output image stack is bdb:sdata. Please use it as an input image stack of sx3dvariability."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "bdb2d_stack"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("output_dir"); token_edit_list.append(token_edit)

	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("sym"); token_edit_list.append(token_edit)
	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("3D Variability Preprocess", None, token_edit_list, sxsubcmd_mpi_support)

	return sxcmd_subconfig

def create_sxcmd_subconfig_meridien_20171120_local():
	token_edit_list = []
	# token_edit = SXcmd_token(); token_edit.initialize_edit("continue_from_subset"); token_edit.is_required = True; token_edit.is_locked = True; token_edit.default = True; token_edit.restore = True; token_edit_list.append(token_edit)
	# token_edit = SXcmd_token(); token_edit.initialize_edit("ctrefromsort3d"); token_edit.is_required = True; token_edit.is_locked = True; token_edit.default = True; token_edit.restore = True; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("ctref"); token_edit.is_required = True; token_edit.is_locked = True; token_edit.default = True; token_edit.restore = True; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("output_directory"); token_edit_list.append(token_edit)
	# token_edit = SXcmd_token(); token_edit.initialize_edit("subset"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("ctref_subset"); token_edit_list.append(token_edit)
	# token_edit = SXcmd_token(); token_edit.initialize_edit("oldrefdir"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("ctref_oldrefdir"); token_edit_list.append(token_edit)
	# token_edit = SXcmd_token(); token_edit.initialize_edit("continue_from_iter"); token_edit_list.append(token_edit)
	# token_edit = SXcmd_token(); token_edit.initialize_edit("ctrefromiter"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("ctref_iter"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("ctref_initvol"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("ctref_orgstack"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("ctref_smearing"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("ctref_an"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("memory_per_node"); token_edit_list.append(token_edit)
	
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("radius"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("mask3D"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("symmetry"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("inires"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("delta"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("xr"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("ts"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
###	token_edit = SXcmd_token(); token_edit.initialize_edit("center_method"); token_edit.group = "advanced"; token_edit_list.append(token_edit) # 207/03/10 Toshio Moriya: For now disable 2D alignment related options
###	token_edit = SXcmd_token(); token_edit.initialize_edit("target_radius"); token_edit.group = "advanced"; token_edit_list.append(token_edit) # 207/03/10 Toshio Moriya: For now disable 2D alignment related options
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("shake"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("small_memory"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("ref_a"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("ccfpercentage"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("nonorm"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
###	token_edit = SXcmd_token(); token_edit.initialize_edit("function"); token_edit.group = "advanced"; token_edit_list.append(token_edit) # 207/03/10 Toshio Moriya: Not included for continue run at this point
	
	sxsubcmd_mpi_support = True
	sxcmd_subconfig = SXsubcmd_config("Subset Refinement (OLD)", None, token_edit_list, sxsubcmd_mpi_support)

	return sxcmd_subconfig

# NOTE: Toshio Moriya 2017/12/18
# Default values of some options are different between standard refinement (fresh run & simple restart) and local refinement (stack & iteration mode)
# 
# Standard: parser.add_option("--inires",  type="float",  default= 25.0,  help="Resolution of the initial_volume volume (default 25A)")
# Local   : parser.add_option("--inires",  type="float",  default= -1.0,  help="Resolution of the initial_volume volume (default 25A)")
# 
# Standard: parser.add_option("--delta",   type="float",  default= 7.50,  help="Initial angular sampling step (default 7.5)")
# Local   : parser.add_option("--delta",   type="float",  default= 3.75,  help="Initial angular sampling step (default 7.5)")
# 
def add_sxcmd_subconfig_meridien_shared(token_edit_list):
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("radius"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("mask3D"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("symmetry"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("xr"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("ts"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("memory_per_node"); token_edit_list.append(token_edit)

	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("theta_min"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("theta_max"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("even_angle_method"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("an"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("shake"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("small_memory"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("ccfpercentage"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("nonorm"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("group_by"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("limit_improvement"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("a_criterion"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("function"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("function_ai"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("group_id"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("filament_width"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("helical_rise"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("plot_ang_dist"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("main000"); token_edit_list.append(token_edit)

def add_sxcmd_subconfig_meridien_standard_shared(token_edit_list):
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("inires"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("delta"); token_edit_list.append(token_edit)
	add_sxcmd_subconfig_meridien_shared(token_edit_list)


def add_sxcmd_subconfig_meridien_local_shared_refine(token_edit_list):
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("inires"); token_edit.help = "Resolution of the initial_volume structure. For local refinement, the program automatically calculates the initial resolution using provided orientation parameters."; token_edit.default = -1.0; token_edit.restore = -1.0; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("delta"); token_edit.help = "Initial angular sampling step. For local refinement, the value has to be less than or equal to 3.75."; token_edit.default = 3.75; token_edit.restore = 3.75; token_edit_list.append(token_edit)
	add_sxcmd_subconfig_meridien_shared(token_edit_list)

def create_sxcmd_subconfig_meridien_standard_fresh():
	token_edit_list = []
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("stack"); token_edit.is_required = True; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("output_directory"); token_edit.type = "output"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("initial_volume"); token_edit.is_required = True; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("initialshifts"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("skip_prealignment"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("center_method"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("target_radius"); token_edit_list.append(token_edit)

	add_sxcmd_subconfig_meridien_standard_shared(token_edit_list)

	sxsubcmd_mpi_support = True
	sxcmd_subconfig = SXsubcmd_config("3D Refinement", "Performs 3D structure refinement using a quasi-Maximum Likelihood approach.", token_edit_list, sxsubcmd_mpi_support, is_modeless = True, subset_config="new")

	return sxcmd_subconfig

def create_sxcmd_subconfig_meridien_standard_continuation():
	token_edit_list = []
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("output_directory"); token_edit.label = "Meridien Output Directory"; token_edit.help = "This directory must exist. The results will be written here."; token_edit.is_required = True; token_edit_list.append(token_edit)

###	add_sxcmd_subconfig_meridien_standard_shared(token_edit_list)
	# All options should be in advanced tab
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("radius"); token_edit.group = "advanced"; token_edit.type = "int"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("mask3D"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("symmetry"); token_edit.group = "advanced"; token_edit.type = "string"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("inires"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("delta"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("xr"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("ts"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	# NOTE: Toshio Moriya 200148/02/13 These initialshifts and skip_prealignment options are not used in this mode. They are used only for initial 2D alignment
	# token_edit = SXcmd_token(); token_edit.initialize_edit("initialshifts"); token_edit.group = "advanced"; token_edit_list.append(token_edit)
	# token_edit = SXcmd_token(); token_edit.initialize_edit("skip_prealignment"); token_edit.group ="advanced"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("memory_per_node"); token_edit.group = "advanced"; token_edit_list.append(token_edit)

	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("an"); token_edit_list.append(token_edit)
	# NOTE: Toshio Moriya 200148/02/13 These center_method and target_radius options is not used in this mode. It is used only for initial 2D alignment
	# token_edit = SXcmd_token(); token_edit.initialize_edit("center_method"); token_edit_list.append(token_edit)
	# token_edit = SXcmd_token(); token_edit.initialize_edit("target_radius"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("shake"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("small_memory"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("ccfpercentage"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("nonorm"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("limit_improvement"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("a_criterion"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("function"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("function_ai"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("plot_ang_dist"); token_edit_list.append(token_edit)

	sxsubcmd_mpi_support = True
	sxcmd_subconfig = SXsubcmd_config("3D Refinement Restart", "Restart 3D refinement after the last fully finished iteration of meridien run or local refinement run. One can change some parameters, but MPI settings have to be the same.", token_edit_list, sxsubcmd_mpi_support, is_modeless = True, subset_config="restart")

	return sxcmd_subconfig

def create_sxcmd_subconfig_meridien_header_import_xform_projection():
	token_edit_list = []
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("import"); token_edit.is_required = True; token_edit.is_locked = False; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("stack"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("params"); token_edit.is_required = True; token_edit.is_locked = True; token_edit.default = "xform.projection"; token_edit.restore = "xform.projection"; token_edit_list.append(token_edit)

	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("Import Projection Parameters", "Import projection orientation parameters to the header of the input stack. (Five columns: phi, theta, pshi, sx, sy). These parameters are required by \'Local Refinement from Stack\' mode.", token_edit_list, sxsubcmd_mpi_support)

	return sxcmd_subconfig

def create_sxcmd_subconfig_meridien_local_stack():
	token_edit_list = []
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("local_refinement"); token_edit.is_required = True; token_edit.is_locked = True; token_edit.default = True; token_edit.restore = True; token_edit_list.append(token_edit)

	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("stack"); token_edit.help = "The stack must 3D orientation parameters (xform.projection) in image headers. They can be imporeted using sxheader.py."; token_edit.is_required = True; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("output_directory"); token_edit_list.append(token_edit)

	add_sxcmd_subconfig_meridien_local_shared_refine(token_edit_list)

	sxsubcmd_mpi_support = True
	sxcmd_subconfig = SXsubcmd_config("Local Refinement from Stack", "Perform local refinement in which the restricted search begins from the user-provided orientation parameters stored in image headers. Note delta has to be less than or equal to 3.75[A].", token_edit_list, sxsubcmd_mpi_support, subset_config="stack")
	
	return sxcmd_subconfig

def create_sxcmd_subconfig_meridien_local_iteration():
	token_edit_list = []
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("local_refinement"); token_edit.is_required = True; token_edit.is_locked = True; token_edit.default = True; token_edit.restore = True; token_edit_list.append(token_edit)

	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("output_directory"); token_edit.label = "Meridien Directory"; token_edit.help = "This directory must exist. The results will be written there."; token_edit.is_required = True; token_edit_list.append(token_edit)

	add_sxcmd_subconfig_meridien_local_shared_refine(token_edit_list)

	sxsubcmd_mpi_support = True
	sxcmd_subconfig = SXsubcmd_config("Restart Local Refinement", "Restart local refinement after the last fully finished iteration of meridien run or local refinement run. One can change some parameters, but MPI settings have to be the same as in the original run.", token_edit_list, sxsubcmd_mpi_support, subset_config="iteration")

	return sxcmd_subconfig

def create_sxcmd_subconfig_meridien_final():
	token_edit_list = []
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("do_final"); token_edit.is_required = True; token_edit.is_locked = False; token_edit.default = -1; token_edit.restore = -1; token_edit_list.append(token_edit)
	
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("output_directory"); token_edit.label = "Meridien Output Directory"; token_edit.help = "This directory must exist. In this mode information is read from files in this directory and the results will be written there."; token_edit.is_required = True; token_edit_list.append(token_edit)
	
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("memory_per_node"); token_edit_list.append(token_edit)

	sxsubcmd_mpi_support = True
	sxcmd_subconfig = SXsubcmd_config("Final 3D Reconstruction Only", "Compute a final 3D reconstruction using either select or best resolution iteration of meridien.", token_edit_list, sxsubcmd_mpi_support)

	return sxcmd_subconfig

def create_sxcmd_subconfig_sort3d_header_import_xform_projection():
	token_edit_list = []
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("import"); token_edit.is_required = True; token_edit.is_locked = False; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("stack"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("params"); token_edit.is_required = True; token_edit.is_locked = True; token_edit.default = "xform.projection"; token_edit.restore = "xform.projection"; token_edit_list.append(token_edit)

	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("Import Projection Parameters", "Import projection orientation parameters from a file (for example created by sxmeridien) to header of the input stack; they are required by 3DVARIABILITY.py and SORT3D_DEPTH Stack Mode.", token_edit_list, sxsubcmd_mpi_support)

	return sxcmd_subconfig

def add_sxcmd_subconfig_sort3d_depth_shared_sorting(token_edit_list):
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("mask3D"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("focus"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("radius"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("sym"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("img_per_grp"); token_edit_list.append(token_edit)
	#token_edit = SXcmd_token(); token_edit.initialize_edit("img_per_grp_split_rate"); token_edit_list.append(token_edit)
	#token_edit = SXcmd_token(); token_edit.initialize_edit("minimum_grp_size"); token_edit_list.append(token_edit)
	#token_edit = SXcmd_token(); token_edit.initialize_edit("do_swap_au"); token_edit_list.append(token_edit)
	#token_edit = SXcmd_token(); token_edit.initialize_edit("swap_ratio"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("memory_per_node"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("overhead"); token_edit_list.append(token_edit)

	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("depth_order"); token_edit_list.append(token_edit)
	#token_edit = SXcmd_token(); token_edit.initialize_edit("stop_mgskmeans_percentage"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("nsmear"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("orientation_groups"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("not_include_unaccounted"); token_edit_list.append(token_edit)
###	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("shake"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("notapplybckgnoise"); token_edit_list.append(token_edit)
	#token_edit = SXcmd_token(); token_edit.initialize_edit("random_group_elimination_threshold"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("num_core_set"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("nstep"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("use_umat"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("not_freeze_groups"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("compute_on_the_fly"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("check_smearing"); token_edit_list.append(token_edit)


### # NOTE: 2018/01/08 Toshio Moriya
### # post-refiner embedded sort3d_depth is removed recently.
### def add_sxcmd_subconfig_sort3d_depth_shared_postrefiner(token_edit_list):
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("mtf"); token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("B_enhance"); token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("fl"); token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("aa"); token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("B_start"); token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("B_stop"); token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("nofsc_adj"); token_edit_list.append(token_edit)

def create_sxcmd_subconfig_sort3d_depth_iteration():
	token_edit_list = []
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("refinement_dir"); token_edit.is_required = True; token_edit.is_locked = False; token_edit.default = "none"; token_edit.restore = "none"; token_edit_list.append(token_edit)

	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("output_dir"); token_edit_list.append(token_edit)

	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("niter_for_sorting"); token_edit_list.append(token_edit)
	add_sxcmd_subconfig_sort3d_depth_shared_sorting(token_edit_list)
###	# NOTE: 2018/01/08 Toshio Moriya
###	# post-refiner embedded sort3d_depth is removed recently.
###	add_sxcmd_subconfig_sort3d_depth_shared_postrefiner(token_edit_list)

	sxsubcmd_mpi_support = True
	sxcmd_subconfig = SXsubcmd_config("3D Clustering from Iteration - SORT3D_DEPTH", "Initialize from a given iteration of meridien run using the associated parameters, i.e., full set of orientation parameters per image, including orientation probabilities, normalizations and so on.", token_edit_list, sxsubcmd_mpi_support)

	return sxcmd_subconfig

def create_sxcmd_subconfig_sort3d_depth_stack():
	token_edit_list = []
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("instack"); token_edit.is_required = True; token_edit.is_locked = False; token_edit.default = "none"; token_edit.restore = "none"; token_edit_list.append(token_edit)

	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("output_dir"); token_edit_list.append(token_edit)

	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("nxinit"); token_edit_list.append(token_edit)
	add_sxcmd_subconfig_sort3d_depth_shared_sorting(token_edit_list)
###	# NOTE: 2018/01/08 Toshio Moriya
###	# post-refiner embedded sort3d_depth is removed recently.
###	add_sxcmd_subconfig_sort3d_depth_shared_postrefiner(token_edit_list)

	sxsubcmd_mpi_support = True
	sxcmd_subconfig = SXsubcmd_config("3D Clustering from Stack - SORT3D_DEPTH", "Run from user-provided orientation parameters stored in stack header.  This mode uses only highest probability orientation parameters per image, ths it should be avoided.  This program is used only for testing.", token_edit_list, sxsubcmd_mpi_support)

	return sxcmd_subconfig

def create_sxcmd_subconfig_sort3d_makevstack():
	token_edit_list = []
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("makevstack"); token_edit.label = "Output stack subset"; token_edit.is_required = True; token_edit.is_locked = False; token_edit.default = "none"; token_edit.restore = "none"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("input_bdb_stack_file"); token_edit.key_prefix = ""; token_edit.label = "Original image stack"; token_edit.help = "Specify path to input BDB stack file used for the associated SORT3D_DEPTH run."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "bdb2d_stack"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("list"); token_edit.label = "Image selection file"; token_edit.help = "Input selection text file containing a list of selected image IDs (or indexes of the data subset) to create a new virtual BDB image stack from an existed stack or virtual stack. Typically, Cluster#.txt created by sxsort3d_depth (e.g. Cluster1.txt)."; token_edit_list.append(token_edit)

	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("SORT3D_DEPTH Stack Subset", None, token_edit_list, sxsubcmd_mpi_support, subset_config="subset")

	return sxcmd_subconfig

### # NOTE: 2018/01/08 Toshio Moriya
### # post-refiner embedded sort3d_depth is removed recently.
### def create_sxcmd_subconfig_sort3d_depth_postrefiner():
### 	token_edit_list = []
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("post_sorting_sharpen"); token_edit.is_required = True; token_edit.is_locked = True; token_edit.default = True; token_edit.restore = True; token_edit_list.append(token_edit)
### 
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("output_dir"); token_edit_list.append(token_edit)
### 
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("niter_for_sorting"); token_edit_list.append(token_edit)
### 	token_edit = SXcmd_token(); token_edit.initialize_edit("memory_per_node"); token_edit_list.append(token_edit)
### 	add_sxcmd_subconfig_sort3d_depth_shared_postrefiner(token_edit_list)
### 
### 	sxsubcmd_mpi_support = True
### 	sxcmd_subconfig = SXsubcmd_config("PostRefiner of Cluster Volume", "Reconstruct unfiltered maps from a list of selection files and merge the two halves maps, so that users can adjust B-factors and low-pass filter to have proper visualization.", token_edit_list, sxsubcmd_mpi_support)
### 
### 	return sxcmd_subconfig

def create_sxcmd_subconfig_utility_makevstack():
	token_edit_list = []
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("makevstack"); token_edit.is_required = True; token_edit.is_locked = False; token_edit.default = "none"; token_edit.restore = "none"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("input_bdb_stack_file"); token_edit.key_prefix = ""; token_edit.label = "Input BDB image stack"; token_edit.help = "Specify path to input BDB stack file. "; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "bdb2d_stack"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("list"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("exlist"); token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("step"); token_edit_list.append(token_edit)

	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("Create Virtual Stack", None, token_edit_list, sxsubcmd_mpi_support)

	return sxcmd_subconfig

def create_sxcmd_subconfig_utility_changesize():
	token_edit_list = []
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("changesize"); token_edit.is_required = True; token_edit.is_locked = True; token_edit.default = True; token_edit.restore = True; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("input_stack"); token_edit.key_prefix = ""; token_edit.label = "Input image or volume"; token_edit.help = "Input image or volume."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "data2d3d_both"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("output_stack"); token_edit.key_prefix = ""; token_edit.label = "Output resized image or volume"; token_edit.help = "Resized (decimated or interpolated up) image or volume."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = ""; token_edit.type = "output"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("ratio"); token_edit_list.append(token_edit)

	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("Change Size of Image or Volume", None, token_edit_list, sxsubcmd_mpi_support)

	return sxcmd_subconfig

def create_sxcmd_subconfig_utility_window():
	token_edit_list = []
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("clip"); token_edit.label = "Window to specified size [Pixels]"; token_edit.is_required = True; token_edit.is_locked = False; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("input_volume"); token_edit.key_prefix = ""; token_edit.label = "Input volume"; token_edit.help = "Path to input volume file."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = "none"; token_edit.type = "data3d_one"; token_edit_list.append(token_edit)
	token_edit = sxgui_template.SXcmd_token(); token_edit.initialize_edit("output_file"); token_edit.key_prefix = ""; token_edit.label = "Output windowed volume"; token_edit.help = "Path to output windowed (clipped/padded) volume file."; token_edit.group = "main"; token_edit.is_required = True; token_edit.default = "none"; token_edit.type = "output"; token_edit_list.append(token_edit)

	sxsubcmd_mpi_support = False
	sxcmd_subconfig = SXsubcmd_config("Window Volume", None, token_edit_list, sxsubcmd_mpi_support)

	return sxcmd_subconfig

# ========================================================================================
def create_exclude_list_isac2():
	exclude_list = []

	exclude_list.append("yr")

	return exclude_list

# def create_exclude_list_meridien_20171120():
# 	exclude_list = []
# 
# 	# exclude_list.append("continue_from_subset")
# 	# exclude_list.append("ctrefromsort3d")
# 	exclude_list.append("ctref")
# 	# exclude_list.append("subset")
# 	exclude_list.append("ctref_subset")
# 	# exclude_list.append("oldrefdir")
# 	exclude_list.append("ctref_oldrefdir")
# 	# exclude_list.append("continue_from_iter")
# 	# exclude_list.append("ctrefromiter")
# 	exclude_list.append("ctref_iter")
# 	exclude_list.append("ctref_initvol")
# 	exclude_list.append("ctref_orgstack")
# 	exclude_list.append("ctref_smearing")
#	exclude_list.append("ctref_an")
# 
# 	return exclude_list

def create_exclude_list_sort3d():
	exclude_list = []
	exclude_list.append("seed")
	exclude_list.append("sausage")
	
	return exclude_list

### def create_exclude_list_sort3d_new():
###	exclude_list = []
###
###	exclude_list.append("instack")
###
###	return exclude_list

def create_exclude_list_boxer_old():
	exclude_list = []

	exclude_list.append("write_dbbox")
	exclude_list.append("write_ptcls")
	exclude_list.append("exclude_edges")
	exclude_list.append("force")
	exclude_list.append("format")
	exclude_list.append("norm")
	exclude_list.append("suffix")
	exclude_list.append("dbls")
	exclude_list.append("autoboxer")
	exclude_list.append("ppid")
	exclude_list.append("gui")
	# exclude_list.append("verbose")
	exclude_list.append("gauss_autoboxer")
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

def create_exclude_list_boxer():
	exclude_list = []

	exclude_list.append("ppid")
	exclude_list.append("verbose")

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
def build_config_list_MoinMoinWiki():
	# --------------------------------------------------------------------------------
	# Get all necessary informations from wiki documents of sx*.py scripts in MoinMoinWiki format
	# and create gui generation parameter
	# --------------------------------------------------------------------------------
	sxcmd_config_list = []

	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_cter"

	sxcmd_role = "sxr_pipe"
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/cter.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list=["stack_mode"]))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/gui_cter.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, is_submittable = False))

	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/e2display.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/pipe_organize_micrographs.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_window"

	sxcmd_role = "sxr_pipe"
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/e2boxer_old.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_boxer_old(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/window.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/e2bdb.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_window_makevstack()))

	sxcmd_role = "sxr_alt"
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/e2boxer.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_boxer(), is_submittable = False))

	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/e2display.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/pipe_organize_micrographs.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))

	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_isac"

	sxcmd_role = "sxr_pipe"
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/isac2.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_isac2()))
### 	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_isacselect()))
### 	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/e2bdb.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_isac_makevstack()))
### 	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/isac_post_processing.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
###	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/compute_isac_avg.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
###	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/compute_isac_avg.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_isac_beautifier_to_model()))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/compute_isac_avg.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/pipe_isac_substack.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))

	sxcmd_role = "sxr_alt"
###	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/compute_isac_avg.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_isac_beautifier_no_adjust()))
###	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/compute_isac_avg.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_isac_beautifier_to_rot_avg()))
##	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/compute_isac_avg.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_isac_beautifier_to_bfactor()))

	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/e2display.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))

	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_viper"

	sxcmd_role = "sxr_pipe"
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/rviper.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
#	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_viper_changesize()))
#	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/e2proc3d.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_viper_window()))
#	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/e2proc3d.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_scale_clip()))

	sxcmd_role = "sxr_alt"
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/viper.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/pdb2em.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))

	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/e2display.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_adaptive_mask3d()))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/pipe_angular_distribution.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))

	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_meridien"

	sxcmd_role = "sxr_pipe"
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/meridien.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_meridien_standard_fresh()))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/gui_meridien.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_postrefiner_halfset_vol()))

	sxcmd_role = "sxr_alt"
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/meridien.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_meridien_local_iteration()))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/meridien.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_meridien_local_stack()))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/meridien.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_meridien_standard_continuation()))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/meridien.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_meridien_final()))
	
	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/e2display.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_adaptive_mask3d()))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/pipe_angular_distribution.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))

	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_sort3d"

	sxcmd_role = "sxr_pipe"
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/header.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_sort3d_header_import_xform_projection()))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/3dvariability.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_variability_preprocess()))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/3dvariability.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list=["symmetrize"]))
###	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/sort3d.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_sort3d()))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/sort3d_depth.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_sort3d_depth_iteration()))
###	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_postrefiner_cluster_vol()))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/e2bdb.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_sort3d_makevstack()))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/meridien.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_meridien_local_stack()))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_postrefiner_halfset_vol()))

	sxcmd_role = "sxr_alt"
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/sort3d_depth.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_sort3d_depth_stack()))
### # NOTE: 2018/01/08 Toshio Moriya
### # post-refiner embedded sort3d_depth is removed recently.
###	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/sort3d_depth.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_sort3d_depth_postrefiner()))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/meridien.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_meridien_final()))

	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/e2display.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_adaptive_mask3d()))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_binary_mask3d()))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/pipe_angular_distribution.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_postrefiner_single_vols()))

	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_localres"

	sxcmd_role = "sxr_pipe"
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/locres.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/filterlocal.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))

	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/e2display.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_adaptive_mask3d()))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/pipe_angular_distribution.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))

	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_movie"

	sxcmd_role = "sxr_pipe"
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/unblur.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/gui_unblur.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, is_submittable = False))

	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/e2display.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/summovie.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/pipe_organize_micrographs.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))

	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_utilities"

	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/e2display.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/pdb2em.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/relion2sphire.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/sphire2relion.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))

	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_adaptive_mask3d()))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_binary_mask3d()))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_utility_changesize()))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/e2proc3d.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_utility_window()))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/pipe_angular_distribution.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
###	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_postrefiner_single_vol()))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_postrefiner_single_vols()))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/unblur.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/summovie.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/header.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/e2bdb.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_utility_makevstack()))
	sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/pipe_organize_micrographs.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
	# sxcmd_config_list.append(SXcmd_config("../doc/MoinMoinWiki/process.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))

	return sxcmd_config_list

# ========================================================================================
def build_config_list_DokuWiki(is_dev_mode = False):
	# --------------------------------------------------------------------------------
	# Get all necessary informations from wiki documents of sx*.py scripts in DokuWiki format
	# and create gui generation parameter
	# --------------------------------------------------------------------------------
	sxcmd_config_list = []

	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_cter"

	sxcmd_role = "sxr_pipe"
	sxcmd_config_list.append(SXcmd_config("../doc/cter.txt", "DokuWiki", sxcmd_category, sxcmd_role, exclude_list=["stack_mode"]))
	sxcmd_config_list.append(SXcmd_config("../doc/gui_cter.txt", "DokuWiki", sxcmd_category, sxcmd_role, is_submittable = False))

	sxcmd_role = "sxr_alt"
	sxcmd_config_list.append(SXcmd_config("../doc/pipe_resample_micrographs.txt", "DokuWiki", sxcmd_category, sxcmd_role))

	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/e2display.txt", "DokuWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/pipe_organize_micrographs.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/batch.txt", "DokuWiki", sxcmd_category, sxcmd_role, is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/ctf_refine.txt", "DokuWiki", sxcmd_category, sxcmd_role, is_submittable=True))
	sxcmd_config_list.append(
		SXcmd_config("../doc/ctf_refine_stack.txt", "DokuWiki", sxcmd_category, sxcmd_role,
					 is_submittable=True))
	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_window"

	sxcmd_role = "sxr_pipe"
	sxcmd_config_list.append(SXcmd_config("../doc/cryolo_predict.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/window.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/e2bdb.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_window_makevstack()))

	sxcmd_role = "sxr_alt"
	sxcmd_config_list.append(SXcmd_config("../doc/e2boxer_old.txt", "DokuWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_boxer_old(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/e2boxer.txt", "DokuWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_boxer(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/pipe_restacking.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/rewindow.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/cryolo_train.txt", "DokuWiki", sxcmd_category, sxcmd_role))

	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/e2display.txt", "DokuWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/pipe_organize_micrographs.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/cryolo_boxmanager.txt", "DokuWiki", sxcmd_category, sxcmd_role, is_submittable=False))
	sxcmd_config_list.append(SXcmd_config("../doc/batch.txt", "DokuWiki", sxcmd_category, sxcmd_role, is_submittable = False))

	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_isac"

	sxcmd_role = "sxr_pipe"
	sxcmd_config_list.append(SXcmd_config("../doc/isac2.txt", "DokuWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_isac2()))
### 	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_isacselect()))
### 	sxcmd_config_list.append(SXcmd_config("../doc/e2bdb.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_isac_makevstack()))
### 	sxcmd_config_list.append(SXcmd_config("../doc/isac_post_processing.txt", "MoinMoinWiki", sxcmd_category, sxcmd_role))
###	sxcmd_config_list.append(SXcmd_config("../doc/compute_isac_avg.txt", "DokuWiki", sxcmd_category, sxcmd_role))
###	sxcmd_config_list.append(SXcmd_config("../doc/compute_isac_avg.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_isac_beautifier_to_model()))
	sxcmd_config_list.append(SXcmd_config("../doc/compute_isac_avg.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/pipe_isac_substack.txt", "DokuWiki", sxcmd_category, sxcmd_role))

###	sxcmd_role = "sxr_alt"
###	sxcmd_config_list.append(SXcmd_config("../doc/compute_isac_avg.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_isac_beautifier_no_adjust()))
###	sxcmd_config_list.append(SXcmd_config("../doc/compute_isac_avg.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_isac_beautifier_to_rot_avg()))
###	sxcmd_config_list.append(SXcmd_config("../doc/compute_isac_avg.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_isac_beautifier_to_bfactor()))

	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/e2bdb.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_utility_makevstack()))
	sxcmd_config_list.append(SXcmd_config("../doc/separate_class.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/e2display.txt", "DokuWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/batch.txt", "DokuWiki", sxcmd_category, sxcmd_role, is_submittable = False))

	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_viper"

	sxcmd_role = "sxr_pipe"
	sxcmd_config_list.append(SXcmd_config("../doc/rviper.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/proj_compare.txt", "DokuWiki", sxcmd_category, sxcmd_role))
#	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_viper_changesize()))
#	sxcmd_config_list.append(SXcmd_config("../doc/e2proc3d.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_viper_window()))
#	sxcmd_config_list.append(SXcmd_config("../doc/e2proc3d.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_scale_clip()))
	sxcmd_config_list.append(SXcmd_config("../doc/pipe_moon_eliminator.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/mask.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	#sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_adaptive_mask3d()))

	sxcmd_role = "sxr_alt"
	sxcmd_config_list.append(SXcmd_config("../doc/viper.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/pdb2em.txt", "DokuWiki", sxcmd_category, sxcmd_role))

	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/e2display.txt", "DokuWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/pipe_angular_distribution.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/batch.txt", "DokuWiki", sxcmd_category, sxcmd_role, is_submittable = False))

	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_meridien"

	sxcmd_role = "sxr_pipe"
	sxcmd_config_list.append(SXcmd_config("../doc/meridien.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_meridien_standard_fresh()))
	sxcmd_config_list.append(SXcmd_config("../doc/gui_meridien.txt", "DokuWiki", sxcmd_category, sxcmd_role, is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_postrefiner_halfset_vol()))

	sxcmd_role = "sxr_alt"
	sxcmd_config_list.append(SXcmd_config("../doc/header.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_meridien_header_import_xform_projection()))
	sxcmd_config_list.append(SXcmd_config("../doc/meridien.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_meridien_local_stack()))
	sxcmd_config_list.append(SXcmd_config("../doc/meridien.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_meridien_local_iteration()))
	sxcmd_config_list.append(SXcmd_config("../doc/meridien.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_meridien_standard_continuation()))
	sxcmd_config_list.append(SXcmd_config("../doc/meridien.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_meridien_final()))
	
	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/e2display.txt", "DokuWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/pipe_moon_eliminator.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/mask.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	#sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_adaptive_mask3d()))
	sxcmd_config_list.append(SXcmd_config("../doc/pipe_angular_distribution.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_balance_angles()))
	sxcmd_config_list.append(SXcmd_config("../doc/batch.txt", "DokuWiki", sxcmd_category, sxcmd_role, is_submittable = False))
	sxcmd_config_list.append(
		SXcmd_config("../doc/ctf_refine.txt", "DokuWiki", sxcmd_category, sxcmd_role,
					 is_submittable=True))
	sxcmd_config_list.append(
		SXcmd_config("../doc/ctf_refine_stack.txt", "DokuWiki", sxcmd_category, sxcmd_role,
					 is_submittable=True))


	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_sort3d"

	sxcmd_role = "sxr_pipe"
	sxcmd_config_list.append(SXcmd_config("../doc/header.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_sort3d_header_import_xform_projection()))
	sxcmd_config_list.append(SXcmd_config("../doc/3dvariability.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_variability_preprocess()))
	sxcmd_config_list.append(SXcmd_config("../doc/3dvariability.txt", "DokuWiki", sxcmd_category, sxcmd_role, exclude_list=["symmetrize"]))
###	sxcmd_config_list.append(SXcmd_config("../doc/sort3d.txt", "DokuWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_sort3d()))
	sxcmd_config_list.append(SXcmd_config("../doc/sort3d_depth.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_sort3d_depth_iteration()))
###	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_postrefiner_cluster_vol()))
	sxcmd_config_list.append(SXcmd_config("../doc/e2bdb.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_sort3d_makevstack()))
	sxcmd_config_list.append(SXcmd_config("../doc/meridien.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_meridien_local_stack()))
	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_postrefiner_halfset_vol()))

	sxcmd_role = "sxr_alt"
	if is_dev_mode:
		sxcmd_config_list.append(SXcmd_config("../doc/sort3d.txt", "DokuWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_sort3d()))
	sxcmd_config_list.append(SXcmd_config("../doc/sort3d_depth.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_sort3d_depth_stack()))
### # NOTE: 2018/01/08 Toshio Moriya
### # post-refiner embedded sort3d_depth is removed recently.
###	sxcmd_config_list.append(SXcmd_config("../doc/sort3d_depth.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_sort3d_depth_postrefiner()))
	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_postrefiner_single_vols()))
	sxcmd_config_list.append(SXcmd_config("../doc/meridien.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_meridien_final()))

	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/e2display.txt", "DokuWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/pipe_moon_eliminator.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/mask.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	#sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_adaptive_mask3d()))
	#sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_binary_mask3d()))
	sxcmd_config_list.append(SXcmd_config("../doc/pipe_angular_distribution.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/batch.txt", "DokuWiki", sxcmd_category, sxcmd_role, is_submittable = False))

	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_localres"

	sxcmd_role = "sxr_pipe"
	sxcmd_config_list.append(SXcmd_config("../doc/locres.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/filterlocal.txt", "DokuWiki", sxcmd_category, sxcmd_role))

	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/e2display.txt", "DokuWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/pipe_moon_eliminator.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/mask.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	#sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_adaptive_mask3d()))
	sxcmd_config_list.append(SXcmd_config("../doc/pipe_angular_distribution.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/batch.txt", "DokuWiki", sxcmd_category, sxcmd_role, is_submittable = False))

	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_movie"

	sxcmd_role = "sxr_pipe"
	sxcmd_config_list.append(SXcmd_config("../doc/unblur.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/gui_unblur.txt", "DokuWiki", sxcmd_category, sxcmd_role, is_submittable = False))

	sxcmd_role = "sxr_alt"
	sxcmd_config_list.append(SXcmd_config("../doc/unblur_old.txt", "DokuWiki", sxcmd_category, sxcmd_role))

	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/e2display.txt", "DokuWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/summovie.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/pipe_organize_micrographs.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/batch.txt", "DokuWiki", sxcmd_category, sxcmd_role, is_submittable = False))

	# --------------------------------------------------------------------------------
	sxcmd_category = "sxc_utilities"

	sxcmd_role = "sxr_util"
	sxcmd_config_list.append(SXcmd_config("../doc/batch.txt", "DokuWiki", sxcmd_category, sxcmd_role, is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/e2display.txt", "DokuWiki", sxcmd_category, sxcmd_role, exclude_list = create_exclude_list_display(), is_submittable = False))
	sxcmd_config_list.append(SXcmd_config("../doc/pdb2em.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/relion2sphire.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(		SXcmd_config("../doc/sphire2relion.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/separate_class.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/proj_compare.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/pipe_moon_eliminator.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/mask.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	#sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_adaptive_mask3d()))
	#sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_binary_mask3d()))
	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_utility_changesize()))
	sxcmd_config_list.append(SXcmd_config("../doc/e2proc3d.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_utility_window()))
	sxcmd_config_list.append(SXcmd_config("../doc/pipe_angular_distribution.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/cryolo_boxmanager.txt", "DokuWiki", sxcmd_category, sxcmd_role, is_submittable=False))
	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_balance_angles()))
###	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_postrefiner_single_vol()))
	sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig = create_sxcmd_subconfig_postrefiner_single_vols()))
	sxcmd_config_list.append(SXcmd_config("../doc/summovie.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/header.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/e2bdb.txt", "DokuWiki", sxcmd_category, sxcmd_role, subconfig=create_sxcmd_subconfig_utility_makevstack()))
	sxcmd_config_list.append(SXcmd_config("../doc/pipe_organize_micrographs.txt", "DokuWiki", sxcmd_category, sxcmd_role))
	sxcmd_config_list.append(SXcmd_config("../doc/ctf_refine.txt", "DokuWiki", sxcmd_category, sxcmd_role, is_submittable=True))
	sxcmd_config_list.append(SXcmd_config("../doc/ctf_refine_stack.txt", "DokuWiki", sxcmd_category, sxcmd_role, is_submittable=True))


	# sxcmd_config_list.append(SXcmd_config("../doc/process.txt", "DokuWiki", sxcmd_category, sxcmd_role))

	return sxcmd_config_list

# ========================================================================================
def main(is_dev_mode = False, is_MoinMoinWiki_mode = False):

	if is_MoinMoinWiki_mode and not os.path.exists("../doc/MoinMoinWiki"):
		print("")
		print("This mode is not available for public users." )
		print("")
		return
	
	# --------------------------------------------------------------------------------
	# Define command categories used in GUI
	# --------------------------------------------------------------------------------
	sxcmd_category_list = []
	sxcmd_category_list.append(sxgui_template.SXcmd_category("sxc_cter", "CTF", "CTF Estimation and CTF Assessment"))
	sxcmd_category_list.append(sxgui_template.SXcmd_category("sxc_window", "Particle Stack", "Particle Coordinates and Particle Extraction"))
	sxcmd_category_list.append(sxgui_template.SXcmd_category("sxc_isac", "2D Clustering", "ISAC2 2D Clustering and Beautifier"))
	sxcmd_category_list.append(sxgui_template.SXcmd_category("sxc_viper", "Initial 3D Modeling", "Initial 3D modeling with VIPER/RVIPER"))
	sxcmd_category_list.append(sxgui_template.SXcmd_category("sxc_meridien", "3D Refinement", "MERIDIEN 3d Refinement and PostRefiner"))
	sxcmd_category_list.append(sxgui_template.SXcmd_category("sxc_sort3d", "3D Clustering", "3D Variability and SROT3D_DEPTH 3D Clustering"))
	sxcmd_category_list.append(sxgui_template.SXcmd_category("sxc_localres", "Local Resolution", "Local Resolution, and Local Filtering"))
	sxcmd_category_list.append(sxgui_template.SXcmd_category("sxc_movie", "Movie Micrograph", "Micrograph Movie Alignemnt and Drift Assessment"))
	sxcmd_category_list.append(sxgui_template.SXcmd_category("sxc_utilities", "Utilities", "Miscellaneous Utlitity Commands"))

	# --------------------------------------------------------------------------------
	# Build configuration list
	# --------------------------------------------------------------------------------
	if is_MoinMoinWiki_mode:
		sxcmd_config_list = build_config_list_MoinMoinWiki()
	else:
		assert (not is_MoinMoinWiki_mode)
		sxcmd_config_list = build_config_list_DokuWiki(is_dev_mode)
		
	# --------------------------------------------------------------------------------
	# Check consistency between sxcmd_category_list and sxcmd_config_list
	# --------------------------------------------------------------------------------
	sxcmd_category_names = []
	for sxcmd_category in sxcmd_category_list:
		sxcmd_category_names.append(sxcmd_category.name)

	for sxcmd_config in sxcmd_config_list:
		if not sxcmd_config.category in sxcmd_category_names:
			sp_global_def.ERROR("Logical Error: sxcmd_config for %s is using invalid category %s." % (sxcmd_config.wiki, sxcmd_config.category), "%s in %s" % (__name__, os.path.basename(__file__)))

	# --------------------------------------------------------------------------------
	# Generate sxgui.py
	# --------------------------------------------------------------------------------
	sxgui_template_file_path = "sxgui_template.py"
	
#	output_file_path = "../bin/sxgui.py"
	output_file_path = "./sxgui_auto.py"
	if is_dev_mode:
		output_file_path = "./sxgui_dev.py"
	elif is_MoinMoinWiki_mode:
		output_file_path = "./sxgui_MoinMoinWiki.py"
	
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
						sp_global_def.ERROR("Logical Error: Invalid Wiki format %s! Check the sxcmd_config setting in this script." % (sxcmd_config.format), "%s in %s" % (__name__, os.path.basename(__file__)))

					if sxcmd_config.subconfig != None:
						apply_sxsubcmd_config(sxcmd_config.subconfig, sxcmd)
					if len(sxcmd_config.exclude_list) > 0:
						apply_exclude_list(sxcmd_config.exclude_list, sxcmd)
					insert_sxcmd_to_file(sxcmd, output_file, sxcmd_variable_name)
			# else: do nothing
		else:
			if current_state != state_insertion: sp_global_def.ERROR("Logical Error: This condition should not happen! State setting must be incorrect.", "%s in %s" % (__name__, os.path.basename(__file__)))
			if line.find("# @@@@@ END_INSERTION @@@@@") != -1:
				current_state = state_template
			# else: do nothing

	if current_state == state_insertion: sp_global_def.ERROR("Script Template Format Error: START_INSERTION and END_INSERTION must be paired.", "%s in %s" % (__name__, os.path.basename(__file__)))

	output_file.close()

	os.system("chmod +x %s" % output_file_path)

# ========================================================================================
if __name__ == '__main__':
	print("")
	print("==================== Creating release version ==================== " )
	print("")
	main()
	print("")
	print("==================== Creating development version ==================== " )
	print("")
	main(is_dev_mode = True)
	print("")
# 	print("==================== Creating MoindMoinWiki version ==================== " )
# 	print("")
# 	main(is_MoinMoinWiki_mode = True)

# ========================================================================================
# END OF SCRIPT
# ========================================================================================
