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
		self.key_base = ""          # key base name of command token (argument or option) in command line
		self.key_prefix = ""        # key prefix of of command token. None for argument, '--' or '-' for option
		self.label = ""             # User friendly name of argument or option
		self.help = ""              # Help info
		self.group = ""             # Tab group: main or advanced
		self.is_required = False    # Required argument or options. No default value are available 
		self.default = ""           # Default value
		self.type = ""              # Type of value
		self.is_in_io = False       # <Used only here> To check consistency between 'usage in command line' and list in '== Input ==' and '== Output ==' sections
		self.widget = None          # <Used only in sxgui.py> Widget instances Associating with this command token
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

# ========================================================================================
class SXcmd:
	def __init__(self):
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.name = ""               # Name of this command (i.e. name of sx*.py script but without .py extension)
		self.short_info = ""         # Short description of this command
		self.mpi_support = False     # Flag to indicate if this command suppors MPI version
		self.mpi_add_flag = False    # NOTE: 2015/11/12 Toshio Moriya. This can be removed when --MPI flag is removed from all sx*.py scripts 
		self.token_list = []         # list of command tokens. Need this to keep the order of command tokens
		self.token_dict = {}         # dictionary of command tokens, organised by key base name of command token. Easy to access a command token but looses their order
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		
# ========================================================================================
def construct_sxcmd_list():
	sxcmd_list = []
	
	# Actual sx commands are inserted into the following section by wikiparser.py.
	# @@@@@ START_INSERTION @@@@@
	sxcmd = SXcmd(); sxcmd.name = "sxcter"; sxcmd.short_info = "Automated estimation of CTF parameters with error assessment."; sxcmd.mpi_support = True; sxcmd.mpi_add_flag = True
	token = SXcmd_token(); token.key_base = "stack"; token.key_prefix = ""; token.label = "2D images in a stack file"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "none"; token.type = "image"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "outdir1"; token.key_prefix = ""; token.label = "output directory for rotinf**** files"; token.help = ""; token.group = "main"; token.is_required = True; token.default = ""; token.type = "output"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "outdir2"; token.key_prefix = ""; token.label = "output directory for partres file"; token.help = ""; token.group = "main"; token.is_required = True; token.default = ""; token.type = "output"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "indir"; token.key_prefix = "--"; token.label = "directory containing micrographs to be processed"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "current directory"; token.type = "directory"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "nameroot"; token.key_prefix = "--"; token.label = "prefix of micrographs to be processed"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "none"; token.type = "string"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "micsuffix"; token.key_prefix = "--"; token.label = "a string denoting micrograph type"; token.help = "For example 'mrc', 'hdf', 'ser' ... "; token.group = "main"; token.is_required = False; token.default = "none"; token.type = "string"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "wn"; token.key_prefix = "--"; token.label = "size of window to use"; token.help = "should be slightly larger than particle box size "; token.group = "main"; token.is_required = False; token.default = "512"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "apix"; token.key_prefix = "--"; token.label = "pixel size in angstroms"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "1.0"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "Cs"; token.key_prefix = "--"; token.label = "microscope Cs (spherical aberration)"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "2.0"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "voltage"; token.key_prefix = "--"; token.label = "microscope voltage in KV"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "300.0"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "ac"; token.key_prefix = "--"; token.label = "amplitude contrast in percentage"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "10.0"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "kboot"; token.key_prefix = "--"; token.label = "number of defocus estimates for micrograph"; token.help = "used for error assessment "; token.group = "advanced"; token.is_required = False; token.default = "16"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "debug"; token.key_prefix = "--"; token.label = "debug"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = False; token.type = "bool"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "overlap_x"; token.key_prefix = "--"; token.label = "overlap x in percentage"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "50"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "overlap_y"; token.key_prefix = "--"; token.label = "overlap y in percentage"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "50"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "edge_x"; token.key_prefix = "--"; token.label = "edge x in pixels"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "0"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "edge_y"; token.key_prefix = "--"; token.label = "edge y in pixels"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "0"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "f_start"; token.key_prefix = "--"; token.label = "starting frequency in 1/A"; token.help = "by default determined automatically "; token.group = "advanced"; token.is_required = False; token.default = "-1.0"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "f_stop"; token.key_prefix = "--"; token.label = "stop frequency in 1/A"; token.help = "by default determined automatically "; token.group = "advanced"; token.is_required = False; token.default = "-1.0"; token.type = "float"; sxcmd.token_list.append(token)

	sxcmd_list.append(sxcmd)

	sxcmd = SXcmd(); sxcmd.name = "sxwindow"; sxcmd.short_info = "Window out particles with known coordinates from a micrograph."; sxcmd.mpi_support = False; sxcmd.mpi_add_flag = False
	token = SXcmd_token(); token.key_base = "micrographs_list"; token.key_prefix = ""; token.label = "Input micrographs"; token.help = "if none is specified all files in option --indir are used "; token.group = "main"; token.is_required = False; token.default = "none"; token.type = "string"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "coords_dir"; token.key_prefix = "--"; token.label = "directory containing input coordinates files"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "current directory"; token.type = "directory"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "coords_suffix"; token.key_prefix = "--"; token.label = "suffix of input coordinates files"; token.help = "For example '_ptcls. "; token.group = "main"; token.is_required = False; token.default = "none"; token.type = "string"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "coords_extension"; token.key_prefix = "--"; token.label = "file extension of input coordinates files"; token.help = "e.g 'box' for eman1, 'json' for eman2, ... "; token.group = "main"; token.is_required = False; token.default = "box"; token.type = "string"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "coords_format"; token.key_prefix = "--"; token.label = "format of input coordinates file"; token.help = "'sparx', 'eman1', 'eman2', or 'spider'. The coordinates of sparx, eman2, and spider format is particle center. The coordinates of eman1 format is particle box conner associated with the original box size. "; token.group = "main"; token.is_required = False; token.default = "eman1"; token.type = "string"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "indir"; token.key_prefix = "--"; token.label = "directory containing input micrographs"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "current directory"; token.type = "directory"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "importctf"; token.key_prefix = "--"; token.label = "file name of sxcter output"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "none"; token.type = "paramters"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "limitctf"; token.key_prefix = "--"; token.label = "Filter micrographs based on the CTF limit"; token.help = "It requires --importctf. "; token.group = "advanced"; token.is_required = False; token.default = False; token.type = "bool"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "resample_ratio"; token.key_prefix = "--"; token.label = "Ratio of new to old image size (or old to new pixel size) for resampling"; token.help = "Valid range is 0.0 < resample_ratio <= 1.0. "; token.group = "advanced"; token.is_required = False; token.default = "1.0"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "box_size"; token.key_prefix = "--"; token.label = "xy dimension of square area to be windowed (in pixels)"; token.help = "Pixel size after resampling is assumed when resample_ratio < 1.0 "; token.group = "main"; token.is_required = False; token.default = "256"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "outdir"; token.key_prefix = "--"; token.label = "output directory"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "current directory"; token.type = "output"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "outsuffix"; token.key_prefix = "--"; token.label = "suffix for output stack"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "_ptcls"; token.type = "string"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "micsuffix"; token.key_prefix = "--"; token.label = "file extension of input micrographs"; token.help = "A string denoting micrograph type. "; token.group = "main"; token.is_required = False; token.default = "hdf"; token.type = "string"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "nameroot"; token.key_prefix = "--"; token.label = "root name (Prefix) of input micrographs"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "none"; token.type = "string"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "invert"; token.key_prefix = "--"; token.label = "invert image contrast"; token.help = "(recommended for cryo data) "; token.group = "main"; token.is_required = False; token.default = False; token.type = "bool"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "defocuserror"; token.key_prefix = "--"; token.label = "Defocus errror limit"; token.help = "exclude micrographs whose relative defocus error as estimated by sxcter is larger than defocuserror percent. The error is computed as (std dev defocus)/defocus*100%. "; token.group = "advanced"; token.is_required = False; token.default = "1000000.0"; token.type = "image"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "astigmatismerror"; token.key_prefix = "--"; token.label = "Astigmatism error limit"; token.help = "Set to zero astigmatism for micrographs whose astigmatism angular error as estimated by sxcter is larger than astigmatismerror degrees. "; token.group = "advanced"; token.is_required = False; token.default = "360.0"; token.type = "float"; sxcmd.token_list.append(token)

	sxcmd_list.append(sxcmd)

	sxcmd = SXcmd(); sxcmd.name = "sxisac"; sxcmd.short_info = "Perform Iterative Stable Alignment and Clustering (ISAC) on a 2-D image stack.  NEW VERSION 07/21/2015."; sxcmd.mpi_support = True; sxcmd.mpi_add_flag = False
	token = SXcmd_token(); token.key_base = "stack_file"; token.key_prefix = ""; token.label = "2-D images in a stack file (format must be bdb)"; token.help = "images have to be square (''nx''=''ny'') "; token.group = "main"; token.is_required = True; token.default = ""; token.type = "image"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "output_directory"; token.key_prefix = ""; token.label = "output directory name"; token.help = "into which the results will be written (if it does not exist, it will be created, if it does exist, the results will be written possibly overwriting previous results) "; token.group = "main"; token.is_required = True; token.default = ""; token.type = "output"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "radius"; token.key_prefix = "--"; token.label = "particle radius"; token.help = "there is no default, a sensible number has to be provided, units - pixels "; token.group = "main"; token.is_required = True; token.default = ""; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "img_per_grp"; token.key_prefix = "--"; token.label = "number of images per class"; token.help = "in the ideal case (essentially maximum size of class) "; token.group = "main"; token.is_required = False; token.default = "100"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "CTF"; token.key_prefix = "--"; token.label = "perform phase flipping of the data"; token.help = "using CTF information included in image headers "; token.group = "main"; token.is_required = False; token.default = False; token.type = "bool"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "restart_section"; token.key_prefix = "--"; token.label = "restart section"; token.help = "each generation (iteration) contains three sections: 'restart', 'candidate_class_averages', and 'reproducible_class_averages'. To restart from a particular step, for example, generation 4 and section 'candidate_class_averages' the following option is needed: '--restart_section=candidate_class_averages,4'. The option requires no white space before or after the comma. The default behavior is to restart execution from where it stopped intentionally or unintentionally. For default restart, it is assumed that the name of the directory is provided as argument. Alternatively, the '--use_latest_master_directory' option can be used. "; token.group = "main"; token.is_required = False; token.default = "none"; token.type = "string"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "ir"; token.key_prefix = "--"; token.label = "inner ring"; token.help = "of the resampling to polar coordinates. units - pixels "; token.group = "advanced"; token.is_required = False; token.default = "1"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "rs"; token.key_prefix = "--"; token.label = "ring step"; token.help = "of the resampling to polar coordinates. units - pixels "; token.group = "advanced"; token.is_required = False; token.default = "1"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "xr"; token.key_prefix = "--"; token.label = "x range"; token.help = "of translational search. By default, set by the program. "; token.group = "advanced"; token.is_required = False; token.default = "-1"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "yr"; token.key_prefix = "--"; token.label = "y range"; token.help = "of translational search. By default, same as xr. "; token.group = "advanced"; token.is_required = False; token.default = "-1"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "ts"; token.key_prefix = "--"; token.label = "search step"; token.help = "of translational search: units - pixels "; token.group = "advanced"; token.is_required = False; token.default = "1.0"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "maxit"; token.key_prefix = "--"; token.label = "number of iterations for reference-free alignment"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "30"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "center_method"; token.key_prefix = "--"; token.label = "method for centering"; token.help = "of global 2D average during initial prealignment of data (0 : no centering; -1 : average shift method; please see center_2D in utilities.py for methods 1-7) "; token.group = "advanced"; token.is_required = False; token.default = "7"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "dst"; token.key_prefix = "--"; token.label = "discrete angle used in within group alignment"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "90.0"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "FL"; token.key_prefix = "--"; token.label = "lowest stopband"; token.help = "frequency used in the tangent filter "; token.group = "advanced"; token.is_required = False; token.default = "0.2"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "FH"; token.key_prefix = "--"; token.label = "highest stopband"; token.help = "frequency used in the tangent filter "; token.group = "advanced"; token.is_required = False; token.default = "0.3"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "FF"; token.key_prefix = "--"; token.label = "fall-off of the tangent filter"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "0.2"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "init_iter"; token.key_prefix = "--"; token.label = "SAC initialization iterations"; token.help = "number of runs of ab-initio within-cluster alignment for stability evaluation in SAC initialization "; token.group = "advanced"; token.is_required = False; token.default = "3"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "main_iter"; token.key_prefix = "--"; token.label = "SAC main iterations"; token.help = "number of runs of ab-initio within-cluster alignment for stability evaluation in SAC "; token.group = "advanced"; token.is_required = False; token.default = "3"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "iter_reali"; token.key_prefix = "--"; token.label = "SAC stability check interval"; token.help = "every iter_reali iterations of SAC stability checking is performed "; token.group = "advanced"; token.is_required = False; token.default = "1"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "match_first"; token.key_prefix = "--"; token.label = "number of iterations to run 2-way matching in the first phase"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "1"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "max_round"; token.key_prefix = "--"; token.label = "maximum rounds"; token.help = "of generating candidate class averages in the first phase "; token.group = "advanced"; token.is_required = False; token.default = "20"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "match_second"; token.key_prefix = "--"; token.label = "number of iterations to run 2-way (or 3-way) matching in the second phase"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "5"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "stab_ali"; token.key_prefix = "--"; token.label = "number of alignments when checking stability"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "5"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "thld_err"; token.key_prefix = "--"; token.label = "threshold of pixel error when checking stability"; token.help = "equals root mean square of distances between corresponding pixels from set of found transformations and theirs average transformation, depends linearly on square of radius (parameter ou). units - pixels. "; token.group = "advanced"; token.is_required = False; token.default = "0.7"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "indep_run"; token.key_prefix = "--"; token.label = "level of m-way matching for reproducibility tests"; token.help = "By default, perform full ISAC to 4-way matching. Value indep_run=2 will restrict ISAC to 2-way matching and 3 to 3-way matching.  Note the number of used MPI processes requested in mpirun must be a multiplicity of indep_run. "; token.group = "advanced"; token.is_required = False; token.default = "4"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "thld_grp"; token.key_prefix = "--"; token.label = "threshold of the size of reproducible class"; token.help = "essentially minimum size of class "; token.group = "advanced"; token.is_required = False; token.default = "10"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "n_generations"; token.key_prefix = "--"; token.label = "maximum number of generations"; token.help = "program stops when reaching this total number of generations: "; token.group = "advanced"; token.is_required = False; token.default = "100"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "rand_seed"; token.key_prefix = "--"; token.label = "random seed set before calculations"; token.help = "useful for testing purposes "; token.group = "advanced"; token.is_required = False; token.default = "none - total randomness"; token.type = "string"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "new"; token.key_prefix = "--"; token.label = "use new code"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = False; token.type = "bool"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "debug"; token.key_prefix = "--"; token.label = "debug info printout"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = False; token.type = "bool"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "use_latest_master_directory"; token.key_prefix = "--"; token.label = "use latest master directory"; token.help = "when active, the program looks for the latest directory that starts with the word 'master', so the user does not need to provide a directory name. "; token.group = "advanced"; token.is_required = False; token.default = False; token.type = "bool"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "stop_after_candidates"; token.key_prefix = "--"; token.label = "stop after candidates"; token.help = "stops after the 'candidate_class_averages' section. "; token.group = "advanced"; token.is_required = False; token.default = False; token.type = "bool"; sxcmd.token_list.append(token)

	sxcmd_list.append(sxcmd)

	sxcmd = SXcmd(); sxcmd.name = "sxviper"; sxcmd.short_info = "Validated ''ab initio'' 3D structure determination, aka Validation of Individual Parameter Reproducibility. The program is designed to determine a validated initial intermediate resolution structure using a small set (<100?) of class averages produced by ISAC [[sxisac]]."; sxcmd.mpi_support = True; sxcmd.mpi_add_flag = False
	token = SXcmd_token(); token.key_base = "stack"; token.key_prefix = ""; token.label = "2D images in a stack file"; token.help = ""; token.group = "main"; token.is_required = True; token.default = ""; token.type = "image"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "directory"; token.key_prefix = ""; token.label = "output directory name"; token.help = "into which the results will be written (if it does not exist, it will be created, if it does exist, the results will be written possibly overwriting previous results) "; token.group = "main"; token.is_required = True; token.default = ""; token.type = "output"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "ir"; token.key_prefix = "--"; token.label = "inner radius for rotational search"; token.help = "> 0 "; token.group = "advanced"; token.is_required = False; token.default = "1"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "radius"; token.key_prefix = "--"; token.label = "radius of the particle"; token.help = "has to be less than < int(nx/2)-1 "; token.group = "main"; token.is_required = True; token.default = ""; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "rs"; token.key_prefix = "--"; token.label = "step between rings in rotational search"; token.help = ">0 "; token.group = "advanced"; token.is_required = False; token.default = "1"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "xr"; token.key_prefix = "--"; token.label = "range for translation search in x direction"; token.help = "search is +/xr in pixels "; token.group = "advanced"; token.is_required = False; token.default = "0"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "yr"; token.key_prefix = "--"; token.label = "range for translation search in y direction"; token.help = "if omitted will be set to xr, search is +/yr in pixels "; token.group = "advanced"; token.is_required = False; token.default = "0"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "mask3D"; token.key_prefix = "--"; token.label = "3D mask file"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "sphere"; token.type = "image"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "moon_elimination"; token.key_prefix = "--"; token.label = "elimination of disconnected pieces"; token.help = "two arguments: mass in KDa and resolution in px/A separated by comma, no space "; token.group = "advanced"; token.is_required = False; token.default = "none"; token.type = "string"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "ts"; token.key_prefix = "--"; token.label = "step size of the translation search in x-y directions"; token.help = "search is -xr, -xr+ts, 0, xr-ts, xr, can be fractional "; token.group = "advanced"; token.is_required = False; token.default = "1.0"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "delta"; token.key_prefix = "--"; token.label = "angular step of reference projections"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "2.0"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "center"; token.key_prefix = "--"; token.label = "centering of 3D template"; token.help = "average shift method; 0: no centering; 1: center of gravity "; token.group = "advanced"; token.is_required = False; token.default = "-1"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "maxit1"; token.key_prefix = "--"; token.label = "maximum number of iterations performed for the GA part"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "400"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "maxit2"; token.key_prefix = "--"; token.label = "maximum number of iterations performed for the finishing up part"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "50"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "L2threshold"; token.key_prefix = "--"; token.label = "stopping criterion of GA"; token.help = "given as a maximum relative dispersion of volumes' L2 norms: "; token.group = "advanced"; token.is_required = False; token.default = "0.03"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "ref_a"; token.key_prefix = "--"; token.label = "method for generating the quasi-uniformly distributed projection directions"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "S"; token.type = "string"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "sym"; token.key_prefix = "--"; token.label = "point-group symmetry of the structure"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "c1"; token.type = "string"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "nruns"; token.key_prefix = "--"; token.label = "GA population"; token.help = "aka number of quasi-independent volumes "; token.group = "advanced"; token.is_required = False; token.default = "6"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "doga"; token.key_prefix = "--"; token.label = "do GA when fraction of orientation changes less than 1.0 degrees is at least doga"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "0.1"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "fl"; token.key_prefix = "--"; token.label = "cut-off frequency applied to the template volume"; token.help = "using a hyperbolic tangent low-pass filter "; token.group = "advanced"; token.is_required = False; token.default = "0.25"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "aa"; token.key_prefix = "--"; token.label = "fall-off of hyperbolic tangent low-pass filter"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "0.1"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "pwreference"; token.key_prefix = "--"; token.label = "text file with a reference power spectrum"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "none"; token.type = "paramters"; sxcmd.token_list.append(token)

	sxcmd_list.append(sxcmd)

	sxcmd = SXcmd(); sxcmd.name = "sxmeridien"; sxcmd.short_info = "Performs 3D structure refinement."; sxcmd.mpi_support = True; sxcmd.mpi_add_flag = False
	token = SXcmd_token(); token.key_base = "stack"; token.key_prefix = ""; token.label = "name of input stack"; token.help = ""; token.group = "main"; token.is_required = True; token.default = ""; token.type = "image"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "output_directory"; token.key_prefix = ""; token.label = "output folder"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "current directory"; token.type = "output"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "initial_volume"; token.key_prefix = ""; token.label = "initial 3D structure"; token.help = ""; token.group = "main"; token.is_required = True; token.default = ""; token.type = "image"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "radius"; token.key_prefix = "--"; token.label = "particle radius"; token.help = "radius of the structure in pixels. if not sure, set to boxsize/2-2 "; token.group = "main"; token.is_required = False; token.default = "-1"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "outlier_percentile"; token.key_prefix = "--"; token.label = "percentile above which outliers"; token.help = "are removed every iteration. "; token.group = "main"; token.is_required = False; token.default = "95.0"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "ref_a"; token.key_prefix = "--"; token.label = "method for generating the quasi-uniformly distributed projection directions"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "S"; token.type = "string"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "sym"; token.key_prefix = "--"; token.label = "point-group symmetry of the structure"; token.help = "cn, dn, where n is multiplicity (for example c5 or d3). "; token.group = "main"; token.is_required = False; token.default = "c1"; token.type = "string"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "startangles"; token.key_prefix = "--"; token.label = "Use orientation parameters in the input file header"; token.help = "to jumpstart the procedure "; token.group = "main"; token.is_required = False; token.default = False; token.type = "bool"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "restrict_shifts"; token.key_prefix = "--"; token.label = "Restrict initial searches for translation"; token.help = "unit - original size pixel. By default, no restriction. "; token.group = "main"; token.is_required = False; token.default = "-1"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "local_filter"; token.key_prefix = "--"; token.label = "Use local filtration"; token.help = "By default, uses generic tangent filter. "; token.group = "main"; token.is_required = False; token.default = False; token.type = "bool"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "smear"; token.key_prefix = "--"; token.label = "Use rotational smear"; token.help = ""; token.group = "main"; token.is_required = False; token.default = False; token.type = "bool"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "sausage"; token.key_prefix = "--"; token.label = "Sausage-making filter"; token.help = ""; token.group = "main"; token.is_required = False; token.default = False; token.type = "bool"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "inires"; token.key_prefix = "--"; token.label = "initial resolution"; token.help = "of the initial_volume: unit - angstroms."; token.group = "main"; token.is_required = False; token.default = "25.0"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "mask3D"; token.key_prefix = "--"; token.label = "3D mask"; token.help = "that defines outline of the structure, preferable with soft edges if not given, set to spherical mask with radius boxsize/2-1. "; token.group = "main"; token.is_required = False; token.default = "none"; token.type = "image"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "CTF"; token.key_prefix = "--"; token.label = "Use CTF"; token.help = ""; token.group = "main"; token.is_required = False; token.default = False; token.type = "bool"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "function"; token.key_prefix = "--"; token.label = "name of the reference preparation function"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = "do_volume_mrk02"; token.type = "function"; sxcmd.token_list.append(token)

	sxcmd_list.append(sxcmd)

	sxcmd = SXcmd(); sxcmd.name = "sx3dvariability"; sxcmd.short_info = "Calculate 3D variability field using a set of aligned 2D projection images as an input. The structures with symmetry require preparing data before calculating variability. The data preparation step would symmetrize the data and output a bdb:sdata for variability calculation."; sxcmd.mpi_support = True; sxcmd.mpi_add_flag = False
	token = SXcmd_token(); token.key_base = "prj_stack"; token.key_prefix = ""; token.label = "stack of 2D images"; token.help = "with 3D orientation parameters in header and (optionally) CTF information "; token.group = "main"; token.is_required = True; token.default = ""; token.type = "image"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "ave2D"; token.key_prefix = "--"; token.label = "write to the disk a stack of 2D averages"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "No"; token.type = "string"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "var2D"; token.key_prefix = "--"; token.label = "write to the disk a stack of 2D variances"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "No"; token.type = "string"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "ave3D"; token.key_prefix = "--"; token.label = "write to the disk reconstructed 3D average"; token.help = "3D reconstruction computed from projections averaged within respective angular neighborhood. It should be used to assess the resolvability and possible artifacts of the variability map. "; token.group = "main"; token.is_required = False; token.default = "No"; token.type = "string"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "var3D"; token.key_prefix = "--"; token.label = "compute 3D variability"; token.help = "time consuming! "; token.group = "main"; token.is_required = False; token.default = "No"; token.type = "string"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "img_per_grp"; token.key_prefix = "--"; token.label = "number of projections"; token.help = "from the angular neighborhood that will be used to estimate 2D variance for each projection data. The larger the number the less noisy the estimate, but the lower the resolution. Usage of large number also results in rotational artifacts in variances that will be visible in 3D variability volume. "; token.group = "main"; token.is_required = False; token.default = "10"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "no_norm"; token.key_prefix = "--"; token.label = "do not use normalization"; token.help = ""; token.group = "main"; token.is_required = False; token.default = False; token.type = "bool"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "radiusvar"; token.key_prefix = "--"; token.label = "radius for 3D var"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "-1"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "npad"; token.key_prefix = "--"; token.label = "number of time to pad the original images"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "2"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "sym"; token.key_prefix = "--"; token.label = "point-group symmetry of the structure"; token.help = "specified in case the input structure has symmetry higher than c1. It is specified together with option --sym in the first step for preparing data. Notice this step can be run with only one CPU and there is no MPI version for it. "; token.group = "main"; token.is_required = False; token.default = "c1"; token.type = "string"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "fl"; token.key_prefix = "--"; token.label = "stop-band frequency of the low pass filter"; token.help = "to be applied to 2D data prior to variability calculation By default, no filtration. "; token.group = "main"; token.is_required = False; token.default = "0.0"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "aa"; token.key_prefix = "--"; token.label = "fall-off frequency of the low pass filter"; token.help = "to be applied to 2D data prior to variability calculation By default, no filtration. "; token.group = "main"; token.is_required = False; token.default = "0.0"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "CTF"; token.key_prefix = "--"; token.label = "use CFT correction"; token.help = ""; token.group = "main"; token.is_required = False; token.default = False; token.type = "bool"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "VERBOSE"; token.key_prefix = "--"; token.label = "Long output for debugging"; token.help = ""; token.group = "advanced"; token.is_required = False; token.default = False; token.type = "bool"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "VAR"; token.key_prefix = "--"; token.label = "stack on input consists of 2D variances"; token.help = ""; token.group = "main"; token.is_required = False; token.default = False; token.type = "bool"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "SND"; token.key_prefix = "--"; token.label = "compute squared normalized differences"; token.help = ""; token.group = "main"; token.is_required = False; token.default = False; token.type = "bool"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "symmetrize"; token.key_prefix = "--"; token.label = "Prepare input stack for handling symmetry"; token.help = ""; token.group = "main"; token.is_required = False; token.default = False; token.type = "bool"; sxcmd.token_list.append(token)

	sxcmd_list.append(sxcmd)

	sxcmd = SXcmd(); sxcmd.name = "sxlocres"; sxcmd.short_info = "Compute local resolution in real space within are outlined by the maskfile and within regions wn x wn x wn."; sxcmd.mpi_support = True; sxcmd.mpi_add_flag = True
	token = SXcmd_token(); token.key_base = "firstvolume"; token.key_prefix = ""; token.label = "first half-volume"; token.help = ""; token.group = "main"; token.is_required = True; token.default = ""; token.type = "image"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "secondvolume"; token.key_prefix = ""; token.label = "second half-volume"; token.help = ""; token.group = "main"; token.is_required = True; token.default = ""; token.type = "image"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "maskfile"; token.key_prefix = ""; token.label = "mask volume"; token.help = "outlining the region within which local resolution values will be computed (optional). "; token.group = "main"; token.is_required = False; token.default = "none"; token.type = "image"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "outputfile"; token.key_prefix = ""; token.label = "output local resolution volume"; token.help = "contains, for each voxel, an [[absolute_frequency_units|absolute frequency]] value for which local resolution at this location drops below the specified cut-off FSC value (only regions specified by the mask film or within a sphere have non-zero values). "; token.group = "main"; token.is_required = True; token.default = ""; token.type = "output"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "wn"; token.key_prefix = "--"; token.label = "size of window within which local real-space FSC is computed"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "7"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "step"; token.key_prefix = "--"; token.label = "shell step in Fourier size in pixels"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "1.0"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "cutoff"; token.key_prefix = "--"; token.label = "resolution cut-off for FSC"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "0.5"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "radius"; token.key_prefix = "--"; token.label = "radius for the mask in pixels"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "-1"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "fsc"; token.key_prefix = "--"; token.label = "name output file"; token.help = "that will contain the overall FSC curve computed by rotational averaging of local resolution values (migh be truncated) "; token.group = "main"; token.is_required = False; token.default = "no curve"; token.type = "string"; sxcmd.token_list.append(token)

	sxcmd_list.append(sxcmd)

	sxcmd = SXcmd(); sxcmd.name = "sxfilterlocal"; sxcmd.short_info = "Locally filter input volume based on values within the associated local resolution volume ([[sxlocres.py]]) within area outlined by the maskfile."; sxcmd.mpi_support = True; sxcmd.mpi_add_flag = True
	token = SXcmd_token(); token.key_base = "inputvolume"; token.key_prefix = ""; token.label = "input volume"; token.help = ""; token.group = "main"; token.is_required = True; token.default = ""; token.type = "image"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "locresvolume"; token.key_prefix = ""; token.label = "local resolution volume"; token.help = "as produced by [[sxlocres.py]]. "; token.group = "main"; token.is_required = True; token.default = ""; token.type = "output"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "maskfile"; token.key_prefix = ""; token.label = "mask volume"; token.help = "outlining the region within which local filtration will be applied (optional). "; token.group = "main"; token.is_required = False; token.default = "none"; token.type = "image"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "outputfile"; token.key_prefix = ""; token.label = "locally filtered volume"; token.help = ""; token.group = "main"; token.is_required = True; token.default = ""; token.type = "output"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "radius"; token.key_prefix = "--"; token.label = "radius for the mask in pixels"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "-1"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "falloff"; token.key_prefix = "--"; token.label = "fall-off of low-pass filter"; token.help = "program uses [[filt_tanl|tangent low-pass filter]]. unit - [[absolute_frequency_units|absolute frequency units]]. "; token.group = "main"; token.is_required = False; token.default = "0.1"; token.type = "float"; sxcmd.token_list.append(token)

	sxcmd_list.append(sxcmd)

	sxcmd = SXcmd(); sxcmd.name = "sxsort3d"; sxcmd.short_info = "Sort out 3-D heterogeneity of 2-D data whose 3-D reconstruction parameters (xform.projection) have been determined already using 3-D sorting protocol I (P1)."; sxcmd.mpi_support = True; sxcmd.mpi_add_flag = False
	token = SXcmd_token(); token.key_base = "stack"; token.key_prefix = ""; token.label = "2D images in a stack file"; token.help = ""; token.group = "main"; token.is_required = True; token.default = ""; token.type = "image"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "outdir"; token.key_prefix = ""; token.label = "master output directory"; token.help = "will contain multiple subdirectories. There is a log.txt that describes the sequences of computations in the program. "; token.group = "main"; token.is_required = True; token.default = ""; token.type = "output"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "mask"; token.key_prefix = ""; token.label = "3D mask"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "none"; token.type = "image"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "focus"; token.key_prefix = "--"; token.label = "3D mask for focused clustering"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "none"; token.type = "image"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "ir"; token.key_prefix = "--"; token.label = "inner radius for rotational correlation"; token.help = "> 0. "; token.group = "main"; token.is_required = False; token.default = "1"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "radius"; token.key_prefix = "--"; token.label = "outer radius for rotational correlation"; token.help = "< nx - 1. Please set to the radius of the particle. "; token.group = "main"; token.is_required = False; token.default = "-1"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "maxit"; token.key_prefix = "--"; token.label = "maximum number of iteration"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "25"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "rs"; token.key_prefix = "--"; token.label = "step between rings in rotational correlation"; token.help = "> 0. "; token.group = "main"; token.is_required = False; token.default = "1"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "xr"; token.key_prefix = "--"; token.label = "range for translation search in x direction"; token.help = "search is +/-xr. "; token.group = "main"; token.is_required = False; token.default = "1"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "yr"; token.key_prefix = "--"; token.label = "range for translation search in y direction"; token.help = "search is +/-yr By default, same as xr. "; token.group = "main"; token.is_required = False; token.default = "-1"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "ts"; token.key_prefix = "--"; token.label = "step size of the translation search"; token.help = "in both directions direction. search is -xr, -xr+ts, 0, xr-ts, xr. "; token.group = "main"; token.is_required = False; token.default = "0.25"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "delta"; token.key_prefix = "--"; token.label = "angular step of reference projections"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "2"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "an"; token.key_prefix = "--"; token.label = "angular neighborhood for local searches"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "-1"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "center"; token.key_prefix = "--"; token.label = "centering method"; token.help = "0 - if you do not want the volume to be centered, 1 - center the volume using cog. "; token.group = "main"; token.is_required = False; token.default = "0"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "nassign"; token.key_prefix = "--"; token.label = "number of reassignment iterations"; token.help = "performed for each angular step. "; token.group = "main"; token.is_required = False; token.default = "1"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "nrefine"; token.key_prefix = "--"; token.label = "number of alignment iterations"; token.help = "performed for each angular step. "; token.group = "main"; token.is_required = False; token.default = "0"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "CTF"; token.key_prefix = "--"; token.label = "Consider CTF correction"; token.help = "during the alignment. "; token.group = "main"; token.is_required = False; token.default = False; token.type = "bool"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "stoprnct"; token.key_prefix = "--"; token.label = "Minimum percentage of assignment change to stop the program"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "3.0"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "sym"; token.key_prefix = "--"; token.label = "point-group symmetry of the structure"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "c1"; token.type = "string"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "function"; token.key_prefix = "--"; token.label = "name of the reference preparation function"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "ref_sort3d"; token.type = "function"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "independent"; token.key_prefix = "--"; token.label = "number of independent run"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "3"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "number_of_images_per_group"; token.key_prefix = "--"; token.label = "number of images per group"; token.help = "critical number defined by user. "; token.group = "main"; token.is_required = False; token.default = "1000"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "low_pass_filter"; token.key_prefix = "--"; token.label = "absolute frequency of low-pass filter"; token.help = "for 3d sorting on the original image size. "; token.group = "main"; token.is_required = False; token.default = "-1.0"; token.type = "float"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "nxinit"; token.key_prefix = "--"; token.label = "initial image size for sorting"; token.help = ""; token.group = "main"; token.is_required = False; token.default = "64"; token.type = "int"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "unaccounted"; token.key_prefix = "--"; token.label = "reconstruct the unaccounted images"; token.help = ""; token.group = "main"; token.is_required = False; token.default = False; token.type = "bool"; sxcmd.token_list.append(token)
	token = SXcmd_token(); token.key_base = "seed"; token.key_prefix = "--"; token.label = "random seed"; token.help = "for create initial random assignment for EQ Kmeans "; token.group = "main"; token.is_required = False; token.default = "-1"; token.type = "int"; sxcmd.token_list.append(token)

	sxcmd_list.append(sxcmd)

	# @@@@@ END_INSERTION @@@@@
	
	# Create command token dictionary for each SXcmd instance
	for sxcmd in sxcmd_list:
		for token in sxcmd.token_list:
				sxcmd.token_dict[token.key_base] = token
	
	return sxcmd_list

# ========================================================================================
# Provides all necessary functionarity
# tabs only contains gui and knows how to layout them
class SXPopup(QWidget):
	def __init__(self, sxcmd):
		QWidget.__init__(self)
		
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.sxcmd = sxcmd
		
		self.projct_dir = "sxgui_settings"
		self.gui_settings_file_path = "%s/gui_settings_%s.txt" % (self.projct_dir, self.sxcmd.name)
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
		
		# Load the previously saved parameter setting of this sx command
		if os.path.exists(self.gui_settings_file_path):
			self.read_params(self.gui_settings_file_path)
		
	def map_widgets_to_sxcmd_line(self):
		# Add program name to command line
		sxcmd_line = "%s.py" % self.sxcmd.name
		
		# Loop through all command tokens
		for token in self.sxcmd.token_list:
			if token.type == 'bool':
				if token.is_required == True and self.key_prefix == "--": ERROR("Logical Error: Encountered unexpected condition for bool type token (%s) of command (%s). Consult with the developer." % (token.key_base, self.name), "%s in %s" % (__name__, os.path.basename(__file__)))
				if (token.widget.checkState() == Qt.Checked) != token.default:
					sxcmd_line += " %s%s" % (token.key_prefix, token.key_base)
			else:
				if token.is_required == True and token.widget.text() == token.default:
					ERROR("Warning: Token (%s) of command (%s) is required. Please set the value for this token." % (token.key_base, self.sxcmd.name), "%s in %s" % (__name__, os.path.basename(__file__)), action = 0)
					return ""
				
				if token.widget.text() != token.default:
					# For now, using line edit box for the other type
					if token.key_prefix == "":
						sxcmd_line += " %s" % (token.widget.text())
					elif token.key_prefix == "--":
						sxcmd_line += " %s%s=%s" % (token.key_prefix, token.key_base, token.widget.text())
					else:
						ERROR("Logical Error: Encountered unexpected prefix for token (%s) of command (%s). Consult with the developer." % (token.key_base, self.name), "%s in %s" % (__name__, os.path.basename(__file__)))
				
				# if token.type == "output":
				# elif token.type == "directory":
				# elif token.type == "image":
				# elif token.type == "parameters":
				# elif token.type == "pdb":
				# elif token.type == "function":
				# else:
				#	if token.type not in ["int", "float", "string"]: ERROR("Logical Error: Encountered unsupported type (%s). Consult with the developer." % line_wiki, "%s in %s" % (__name__, os.path.basename(__file__)))
		
		return sxcmd_line
	
	def generate_cmd_line(self):
		sxcmd_line = self.map_widgets_to_sxcmd_line()
		
		if sxcmd_line:
			# SX command line is not empty
			# If mpi is not supported set number of MPI processer (np) to 1
			np = 1
			if self.sxcmd.mpi_support:
				# mpi is supported
				np = int(str(self.tab_main.mpi_nproc_edit.text()))
				# NOTE: 2015/10/27 Toshio Moriya
				# Since we now assume sx*.py exists in only MPI version, always add --MPI flag if necessary
				# This is not elegant but can be removed when --MPI flag is removed from all sx*.py scripts 
				if self.sxcmd.mpi_add_flag:
					sxcmd_line += ' --MPI'
		
			# Generate command line according to the case
			cmd_line = ""
			# Case 1: queue submission is enabled (MPI must be supported)
			if self.tab_main.qsub_enable_checkbox.checkState() == Qt.Checked:
				if self.sxcmd.mpi_support == False: ERROR("Logical Error: Encountered unexpected condition for self.sxcmd.mpi_support. Consult with the developer.", "%s in %s" % (__name__, os.path.basename(__file__)))
				# Create script for queue submission from a give template
				if os.path.exists(self.tab_main.qsub_script_edit.text()) != True: ERROR("Run Time Error: Invalid file path for qsub script template.", "%s in %s" % (__name__, os.path.basename(__file__)))
				file_template = open(self.tab_main.qsub_script_edit.text(),'r')
				# Extract command line from qsub script template 
				for line in file_template:
					if line.find('XXX_SXCMD_LINE_XXX') != -1:
						cmd_line = line.replace('XXX_SXCMD_LINE_XXX', sxcmd_line)
						if cmd_line.find('XXX_SXMPI_NPROC_XXX') != -1:
							cmd_line = cmd_line.replace('XXX_SXMPI_NPROC_XXX', str(np))
						if cmd_line.find('XXX_SXMPI_JOB_NAME_XXX') != -1:
							cmd_line = cmd_line.replace('XXX_SXMPI_JOB_NAME_XXX', str(self.tab_main.qsub_job_name_edit.text()))
				file_template.close()
			# Case 2: queue submission is disabled, but MPI is supported
			elif self.sxcmd.mpi_support:
				if self.tab_main.qsub_enable_checkbox.checkState() == Qt.Checked: ERROR("Logical Error: Encountered unexpected condition for tab_main.qsub_enable_checkbox.checkState. Consult with the developer.", "%s in %s" % (__name__, os.path.basename(__file__)))
				# Add MPI execution to command line
				cmd_line = str(self.tab_main.mpi_cmd_line_edit.text())
				# If empty string is entered, use a default template
				if cmd_line == "":
					cmd_line = 'mpirun -np XXX_SXMPI_NPROC_XXX XXX_SXCMD_LINE_XXX'
				if cmd_line.find('XXX_SXMPI_NPROC_XXX') != -1:
					cmd_line = cmd_line.replace('XXX_SXMPI_NPROC_XXX', str(np))
				if cmd_line.find('XXX_SXCMD_LINE_XXX') != -1:
					cmd_line = cmd_line.replace('XXX_SXCMD_LINE_XXX', sxcmd_line)
			# Case 3: MPI is not supported
			else: 
				if self.sxcmd.mpi_support == True: ERROR("Logical Error: Encountered unexpected condition for self.sxcmd.mpi_support. Consult with the developer.", "%s in %s" % (__name__, os.path.basename(__file__)))
				# Use sx command as it is
				cmd_line = sxcmd_line
		else:
			# SX command line is be empty because an error happens in map_widgets_to_sxcmd_line
			cmd_line = ""
		
		return cmd_line
	
	def execute_cmd_line(self):
		# out_dir = str(self.cmd_token_dict['Output Directory'].gui.text())
		# if os.path.exists(out_dir):
		# 	print "Output directory " + out_dir + " already exists!"
		# 	return
		
		cmd_line = self.generate_cmd_line()
		
		if cmd_line:
			# Command line is not empty
			# If mpi is not supported set number of MPI processer (np) to 1
			np = 1
			if self.sxcmd.mpi_support:
				np = int(str(self.tab_main.mpi_nproc_edit.text()))
		
			# Case 1: queue submission is enabled (MPI must be supported)
			if self.tab_main.qsub_enable_checkbox.checkState() == Qt.Checked:
				if self.sxcmd.mpi_support == False: ERROR("Logical Error: Encountered unexpected condition for self.sxcmd.mpi_support. Consult with the developer.", "%s in %s" % (__name__, os.path.basename(__file__)))
				# Create script for queue submission from a give template
				template_file_path = self.tab_main.qsub_script_edit.text()
				if os.path.exists(template_file_path) == False: 
					ERROR("WARNING: Invalid file path for qsub script template (%s)." % (template_file_path), "%s in %s" % (__name__, os.path.basename(__file__)), action = 0)
					return
				file_template = open(self.tab_main.qsub_script_edit.text(),'r')
				file_name_qsub_script = 'qsub_' + str(self.tab_main.qsub_job_name_edit.text()) + '.sh'
				file_qsub_script = open(file_name_qsub_script,'w')
				for line_io in file_template:
					if line_io.find('XXX_SXCMD_LINE_XXX') != -1:
						line_io = cmd_line
					else:
						if line_io.find('XXX_SXMPI_NPROC_XXX') != -1:
							line_io = line_io.replace('XXX_SXMPI_NPROC_XXX', str(np))
						if line_io.find('XXX_SXMPI_JOB_NAME_XXX') != -1:
							line_io = line_io.replace('XXX_SXMPI_JOB_NAME_XXX', str(self.tab_main.qsub_job_name_edit.text()))
					file_qsub_script.write(line_io)
				file_template.close()
				file_qsub_script.close()
				# Generate command line for queue submission
				cmd_line_in_script = cmd_line
				cmd_line = str(self.tab_main.qsub_cmd_edit.text()) + ' ' + file_name_qsub_script
				print 'Wrote the following command line in the queue submission script: '
				print cmd_line_in_script
				print 'Submitted a job by the following command: '
				print cmd_line
			# Case 2: queue submission is disabled, but MPI is supported
			else:
				if self.tab_main.qsub_enable_checkbox.checkState() == Qt.Checked or self.sxcmd.mpi_support == False: ERROR("Logical Error: Encountered unexpected condition for tab_main.qsub_enable_checkbox.checkState. Consult with the developer.", "%s in %s" % (__name__, os.path.basename(__file__)))
				print 'Executed the following command: '
				print cmd_line
		
			# Execute the generated command line
			process = subprocess.Popen(cmd_line, shell=True)
			self.emit(QtCore.SIGNAL("process_started"), process.pid)
			
			# Save the current state of GUI settings
			if os.path.exists(self.projct_dir) == False:
				os.mkdir(self.projct_dir)
			self.write_params(self.gui_settings_file_path)
		# else: SX command line is be empty because an error happens in generate_cmd_line. Let's do nothing
	
	def save_cmd_line(self):
		cmd_line = self.generate_cmd_line()
		if cmd_line:
			file_name_out = QtGui.QFileDialog.getSaveFileName(self, "Generate Command Line", options = QtGui.QFileDialog.DontUseNativeDialog)
			if file_name_out != '':
				file_out = open(file_name_out,'w')
				file_out.write(cmd_line + '\n')
				file_out.close()
				print 'Saved the following command to %s:' % file_name_out
				print cmd_line
		# else: Do nothing
	
	def write_params(self, file_name_out):
		file_out = open(file_name_out,'w')
		
		# Write script name for consistency check upon loading
		file_out.write('@@@@@ %s gui setting - ' % (self.sxcmd.name))
		file_out.write(EMANVERSION + ' (CVS' + CVSDATESTAMP[6:-2] +')')
		file_out.write(' @@@@@ \n')
		
		# Define list of (tab) groups
		group_main = "main"
		group_advanced = "advanced"
		
		# Loop through all states
		for group in [group_main, group_advanced]:
			# Loop through all command tokens
			for cmd_token in self.sxcmd.token_list:
				if cmd_token.group == group:
					val_str = ''
					if cmd_token.type == 'bool':
						if cmd_token.widget.checkState() == Qt.Checked:
							val_str = 'YES'
						else:
							val_str = 'NO'
					else:
						# For now, use line edit box for the other type
						val_str = str(cmd_token.widget.text())
						# if cmd_token.type == "output":
						# elif cmd_token.type == "directory":
						# elif cmd_token.type == "image":
						# elif cmd_token.type == "parameters":
						# elif cmd_token.type == "pdb":
						# elif cmd_token.type == "function":
						# else:
						#	if cmd_token.type not in ["int", "float", "string"]: ERROR("Logical Error: Encountered unsupported type (%s). Consult with the developer."  % line_wiki, "%s in %s" % (__name__, os.path.basename(__file__)))
					if cmd_token.is_required == False:
						file_out.write('<%s> %s (default %s) == %s \n' % (cmd_token.key_base, cmd_token.label, cmd_token.default, val_str))
					else:
						file_out.write('<%s> %s (default required %s) == %s \n' % (cmd_token.key_base, cmd_token.label, cmd_token.type, val_str))
				# else: do nothig
		
		# At the end of parameter file...
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
			
	def read_params(self, file_name_in):
		file_in = open(file_name_in,'r')
	
		# Check if this parameter file is for this sx script
		line_in = file_in.readline()
		if line_in.find('@@@@@ %s gui setting' % (self.sxcmd.name)) != -1:
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
					# Extract key_base of this command token
					target_operator = "<"
					item_tail = label_in.find(target_operator)
					if item_tail != 0: ERROR("Paramter File Format Error: Command token entry should start from \"%s\" for key base name in line (%s)" % (target_operator, line_in), "%s in %s" % (__name__, os.path.basename(__file__)))
					label_in = label_in[item_tail + len(target_operator):].strip() # Get the rest of line
					target_operator = ">"
					item_tail = label_in.find(target_operator)
					if item_tail == -1: ERROR("Paramter File Format Error: Command token entry should have \"%s\" closing key base name in line (%s)" % (target_operator, line_in), "%s in %s" % (__name__, os.path.basename(__file__)))
					key_base = label_in[0:item_tail]
					# Get corresponding cmd_token
					if key_base not in self.sxcmd.token_dict.keys(): ERROR("Paramter File Format Error: Command token entry should start from \"%s\" for key base name in line %s" % (target_operator, line_in), "%s in %s" % (__name__, os.path.basename(__file__)))
					cmd_token = self.sxcmd.token_dict[key_base]
					
					if cmd_token.type == "bool":
						# construct new widget(s) for this command token
						if val_str_in == 'YES':
							cmd_token.widget.setChecked(True)
						else: # val_str_in == 'NO'
							cmd_token.widget.setChecked(False)
					else:
						# For now, use line edit box for the other type
						cmd_token.widget.setText(val_str_in)
						# if cmd_token.type == "output":
						# elif cmd_token.type == "directory":
						# elif cmd_token.type == "image":
						# elif cmd_token.type == "parameters":
						# elif cmd_token.type == "pdb":
						# elif cmd_token.type == "function":
						# else:
						#	if cmd_token.type not in ["int", "float", "string"]: ERROR("Logical Error: Encountered unsupported type (%s). Consult with the developer."  % line_wiki, "%s in %s" % (__name__, os.path.basename(__file__)))
		else:
			QMessageBox.warning(self, 'Fail to load paramters', 'The specified file is not paramter file for %s.' % self.name)
		
		file_in.close()
	
	def save_params(self):
		file_name_out = QtGui.QFileDialog.getSaveFileName(self, "Save Parameters", options = QtGui.QFileDialog.DontUseNativeDialog)
		if file_name_out != '':
			self.write_params(file_name_out)
	
	def load_params(self):
		file_name_in = QtGui.QFileDialog.getOpenFileName(self, "Load parameters", options = QtGui.QFileDialog.DontUseNativeDialog)
		if file_name_in != '':
			self.read_params(file_name_in)
	
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

	def __init__(self, parent):
		QWidget.__init__(self, parent)
		
		# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		# class variables
		self.sxpopup = parent
		
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
		temp_label = QtGui.QLabel('<b>%s</b>' % (self.sxpopup.sxcmd.name), self)
		temp_label.move(self.x1,self.y1)
		# NOTE: 2015/11/17 Toshio Moriya
		# Necessary to separate '<b>%s</b>' from the information for avoiding to invoke the tag interpretations of string
		# e.g. < becomes the escape character
		temp_label = QtGui.QLabel('%s' % (self.sxpopup.sxcmd.short_info), self)
		temp_label.setWordWrap(True)
		temp_label.setFixedWidth(600)
		temp_label.move(self.x1 + 100, self.y1)
		self.y1 += 50
		
		# Add load paramater button 
		self.load_params_btn = QPushButton("Load parameters", self)
		self.load_params_btn.move(self.x1-5,self.y1)
		self.load_params_btn.setToolTip('Load gui parameter settings to retrieve a previously-saved one')
		self.connect(self.load_params_btn, SIGNAL("clicked()"), self.sxpopup.load_params)
		
		# Add widget for editing command args and options
		for cmd_token in self.sxpopup.sxcmd.token_list:
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
					#	if cmd_token.type not in ["int", "float", "string"]: ERROR("Logical Error: Encountered unsupported type (%s). Consult with the developer."  % line_wiki, "%s in %s" % (__name__, os.path.basename(__file__)))
				cmd_token_widget.move(self.x2,self.y2 - 7)
				cmd_token_widget.setToolTip(cmd_token.help)		
				
				self.y2 = self.y2+25
				
				# Register this widget
				cmd_token.widget = cmd_token_widget
				
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
		self.set_widget_enable_state(self.mpi_nproc_edit, self.sxpopup.sxcmd.mpi_support)
		self.set_widget_enable_state(self.mpi_cmd_line_edit, self.sxpopup.sxcmd.mpi_support)

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
		self.qsub_job_name_edit.setText(self.sxpopup.sxcmd.name)
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
		self.connect(self.save_params_btn, SIGNAL("clicked()"), self.sxpopup.save_params)
		
		# self.y4 = self.y4+30
		self.y2 = self.y2+30

		self.cmd_line_btn = QPushButton("Generate command line", self)
		# self.cmd_line_btn.move(self.x1-5,  self.y4)
		self.cmd_line_btn.move(self.x1-5,  self.y2)
		self.cmd_line_btn.setToolTip('Generate command line from gui parameter settings')
		self.connect(self.cmd_line_btn, SIGNAL("clicked()"), self.sxpopup.save_cmd_line)
		
		self.y2 = self.y2+30

		# Add a run button
		self.execute_btn = QtGui.QPushButton('Run %s' % self.sxpopup.sxcmd.name, self)
		# make 3D textured push button look
		s = "QPushButton {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #8D0);min-width:90px;margin:5px} QPushButton:pressed {font: bold; color: #000;border: 1px solid #333;border-radius: 11px;padding: 2px;background: qradialgradient(cx: 0, cy: 0,fx: 0.5, fy:0.5,radius: 1, stop: 0 #fff, stop: 1 #084);min-width:90px;margin:5px}"
		self.execute_btn.setStyleSheet(s)
		# self.execute_btn.move(self.x5,  self.y5)
		self.execute_btn.move(self.x5,  self.y2)
		self.connect(self.execute_btn, SIGNAL("clicked()"), self.sxpopup.execute_cmd_line)

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
		if self.sxpopup.sxcmd.mpi_support:
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
					#	if cmd_token.type not in ["int", "float", "string"]: ERROR("Logical Error: Encountered unsupported type (%s). Consult with the developer."  % line_wiki, "%s in %s" % (__name__, os.path.basename(__file__)))
				cmd_token_widget.move(self.x2,self.y1)
				cmd_token_widget.setToolTip(cmd_token.help)		
				
				self.y1 = self.y1+25
				
				# Register this widget
				cmd_token.widget = cmd_token_widget

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
		# SX Commands (sx*.py)
		# --------------------------------------------------------------------------------
		self.y2 = 95
		
		# Construct list of sxscript objects (extracted from associated wiki documents)
		sxcmd_list = construct_sxcmd_list()
		
		for sxcmd in sxcmd_list:
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

