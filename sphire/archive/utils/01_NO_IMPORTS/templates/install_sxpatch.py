#!/usr/bin/env python
from __future__ import print_function
#
# Author: Toshio Moriya 09/06/2017 (toshio.moriya@mpi-dortmund.mpg.de)
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete SPHIRE and EMAN2 software packages have some GPL dependencies,
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
# ========================================================================================
# Imports
# ========================================================================================
# Python Standard Libraries
from __future__ import print_function
pass#IMPORTIMPORTIMPORT import sys
pass#IMPORTIMPORTIMPORT import os
pass#IMPORTIMPORTIMPORT import argparse

# ========================================================================================
# Helper Functions
# ========================================================================================
# ----------------------------------------------------------------------------------------
# Generate command line
# ----------------------------------------------------------------------------------------
def get_cmd_line():
	cmd_line = ""
	for arg in sys.argv:
		cmd_line += arg + "  "
	cmd_line = "Shell line command is \'%s\'" % cmd_line.strip()
	return cmd_line

# ----------------------------------------------------------------------------------------
# Print progress message with time stamp
# ----------------------------------------------------------------------------------------
def print_progress(message):
	pass#IMPORTIMPORTIMPORT from time import strftime, localtime
	time_stamp = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	print(time_stamp, message)

# ----------------------------------------------------------------------------------------
# Generate backup root directory name with current local time
# ----------------------------------------------------------------------------------------
def get_backup_root_dir_name():
	pass#IMPORTIMPORTIMPORT from time import strftime, localtime
	return strftime("backup_%Y%m%d_%H%M%S", localtime())

# ----------------------------------------------------------------------------------------
# Do the actual installation
# ----------------------------------------------------------------------------------------
def install_patch_src_files(patch_src_root_dir_path, patch_src_subdir_rpath, backup_root_dir_path, install_dir_path, install_subdir_rpath, python_script_directive=None):
	pass#IMPORTIMPORTIMPORT import glob
	pass#IMPORTIMPORTIMPORT import shutil
	pass#IMPORTIMPORTIMPORT import stat
	
	print_progress("Checking \'%s\' subdirectory..." % patch_src_subdir_rpath)
	patch_src_subdir_rpath = os.path.join(patch_src_root_dir_path, patch_src_subdir_rpath)
	
	if (os.path.exists(patch_src_subdir_rpath)):
		print_progress("Making the list of patch files in \'%s\'..." % patch_src_subdir_rpath)
		# print_progress("MRK_DEBUG: patch_src_subdir_rpath := %s" % patch_src_subdir_rpath)
		patch_src_file_path_pattern = os.path.join(patch_src_subdir_rpath, "*.py")
		patch_src_file_path_list = glob.glob(patch_src_file_path_pattern)
		print_progress("Found %d Python script files in \'%s\'." % (len(patch_src_file_path_list), patch_src_subdir_rpath))
	
		if len(patch_src_file_path_list) > 0:
			backup_subdir_path = os.path.join(backup_root_dir_path, patch_src_subdir_rpath)
			# print_progress("MRK_DEBUG: backup_subdir_path := %s" % backup_subdir_path)
			if not os.path.exists(backup_subdir_path):
				print_progress("Making backup subdirectory \'%s\'..." % backup_subdir_path)
				os.makedirs(backup_subdir_path)
		else:
			print_progress("Nothing to do..." )
			return
		assert(os.path.exists(backup_subdir_path))
	
		install_subdir_path = os.path.join(install_dir_path, install_subdir_rpath)
		
		for patch_src_file_path in patch_src_file_path_list:
			patch_file_name = os.path.basename(patch_src_file_path)
			install_file_path = os.path.join(install_subdir_path, patch_file_name)
			backup_file_path = os.path.join(backup_subdir_path, patch_file_name)
			# print_progress("MRK_DEBUG: patch_file_name := %s " % patch_file_name)
			# print_progress("MRK_DEBUG: patch_src_file_path := %s " % patch_src_file_path)
			# print_progress("MRK_DEBUG: install_file_path := %s " % install_file_path)
			# print_progress("MRK_DEBUG: backup_file_path := %s " % backup_file_path)
			print_progress("Replacing \'%s\' with \'%s\'..." % (install_file_path, patch_src_file_path))
		
			# Backup the originally-installed file
			if os.path.exists(install_file_path):
				shutil.copy2(install_file_path, backup_file_path)
				
				# Check the consistency of the python script directive if necessary
				install_file = open(install_file_path, "r")
				install_1st_line = install_file.readline().strip()
				# print_progress("MRK_DEBUG: install_1st_line := %s " % install_1st_line)
				if install_1st_line != python_script_directive:
					print_progress("WARNING!!! Found unexpected Python script directive \'%s\' in \'%s\'. However, this should not cause any problems..." % (install_1st_line, install_file_path))
				install_file.close()
			else:
				print_progress("%s does not exist in \'%s\'. The program will not make the backup..." % (patch_file_name, install_subdir_path))
			# assert(not os.path.exists(install_file_path))
		
			# Copy the patch file to install subdirectory
			if python_script_directive is not None:
				# Overwrite python script directive if necessary
				patch_src_file = open(patch_src_file_path, "r")
				patch_src_1st_line = patch_src_file.readline().strip()
				patch_src_remainder_lines = patch_src_file.read()
				patch_src_file.close()
				if patch_src_1st_line != "#!/usr/bin/env python":
					print_progress("WARNING!!! Found unexpected Python script directive \'%s\' in \'%s\'. However, this should not cause any problems..." % (patch_src_1st_line, patch_src_file_path))
			
				install_file = open(install_file_path, "w")
				install_file.write(python_script_directive + "\n")
				install_file.write(patch_src_remainder_lines)
				install_file.close()
			else: 
				# Do simply copy in this case
				shutil.copy2(patch_src_file_path, install_file_path)
				# assert(os.path.exists(install_file_path))
				
			install_file_stat = os.stat(install_file_path)
			os.chmod(install_file_path, install_file_stat.st_mode | 0o111) # Make all executable, Use octal system representation for bits 0o111 = 0b001001001 = Linux 111 = 0d73
	else:
		print_progress("Nothing to do..." )
		return

# ========================================================================================
# Command functions
# ========================================================================================
# ----------------------------------------------------------------------------------------
# TEST COMMAND
# cd /home/moriya/mrk_develop/sxdemo/sxdemo07_20160908/mpi_bdb_ctf
# rm -r mrkout_sxpipe_isac_substack; sxpipe.py isac_substack bdb:beta20161216_pa03b_sxwindow01#data beta20161216_pa04a_sxisac01/class_averages.hdf mrkout_sxpipe_isac_substack
# 
# ----------------------------------------------------------------------------------------
def install_sxpatch(args):
	pass#IMPORTIMPORTIMPORT from distutils.spawn import find_executable
	
	install_dir_path = args.install_dir_path
	if install_dir_path is None:
		sphire_gui_abs_path = find_executable("sphire")
		# print("MRK_DEBUG: sphire_gui_abs_path = ", sphire_gui_abs_path)
		if sphire_gui_abs_path is None:
			print_progress("Program can not find the SPHIRE installation directory. Please manually specify the directory using --install_dir_path options.")
			return
		# assert(os.path.exists(sphire_gui_abs_path))
		install_dir_path = os.path.dirname(os.path.dirname(sphire_gui_abs_path))
	else:
		if not os.path.exists(install_dir_path):
			print_progress("Provided SPHIRE installation directory does not exist. Please check the directory path provided with --install_dir_path options.")
			return
	assert(os.path.exists(install_dir_path))
	print_progress("SPHIRE installation directory is \'%s\'" % (install_dir_path))
	print(" ")
	
	# Define constants
	patch_src_root_dir_path = "src"
	backup_root_dir_path = get_backup_root_dir_name()
	#!/work/software/sphire/beta_20170901/EMAN2/bin/python
	python_script_directive = '#!' +  os.path.join(install_dir_path, "bin/python")
	# print_progress("MRK_DEBUG: python_script_directive := %s " % python_script_directive)
	install_patch_src_files(patch_src_root_dir_path, "eman2/sparx/bin", backup_root_dir_path, install_dir_path, "bin", python_script_directive)
	print(" ")
	install_patch_src_files(patch_src_root_dir_path, "eman2/sparx/lib", backup_root_dir_path, install_dir_path, "lib")
	print(" ")

# ========================================================================================
# Main function
# ========================================================================================
def main():
	# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	# Set up argument parser (supports subcommand)
	# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	parser = argparse.ArgumentParser(description='Python script to install patch for a SPHIRE release.')
	parser.add_argument('--version', action='version', version="Version 0.0.0.0")

	# create the parser for the 'isac_substack' command
	parser.add_argument('--install_dir_path',  type=str,  default=None,  help='Specify path to SPHIRE installation directory. By default, the program try to detect it automatically by extracting the path from \'which sphire\'. (default none)')
	parser.set_defaults(func=install_sxpatch)
	
	# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	# Run specified subcommand
	# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	args = parser.parse_args() # Get namespace object from parser
	# args_dict = vars(parser.parse_args()) # convert it to dictionary object
	# print (args_dict)
	
	print_progress(get_cmd_line())
	print(" ")
	
	# Call the associated function of the specified subcommand
	args.func(args)

	print_progress("DONE!!!")
	print(" ")

# ----------------------------------------------------------------------------------------
if __name__ == '__main__':
	main()

# ========================================================================================
# END OF SCRIPT
# ========================================================================================
