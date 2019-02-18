#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
#
# Author: Markus Stabrin 2018/12/05 (markus.stabrin@mpi-dortmund.mpg.de)
# Copyright (c) 2018 MPI-Dortmund
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

import argparse
import subprocess
import glob
import re
import global_def
import os

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('submission_command', type=str, help='Submission command, e.g., qsub, qsub -V, sbatch, bash')
parser.add_argument('input_run_dir', type=str, help='Directory containin the pipeline submission files')
parser.add_argument('--hold_flag', type=str, default=None, help='Hold flag for the submission command, e.g. -hold_jid')
parser.add_argument('--first_hold_number', type=str, default=None, help='Wait number of an already running job')
args = parser.parse_args()

global_def.BATCH = True
if not os.path.exists(args.input_run_dir):
    global_def.ERROR('Input directory does not exist!', 'sxbatch.py', 1)

qsub_dict = {
    'qsub': re.compile('Your job (\w+)'),
    'sbatch': re.compile('Submitted batch job (\w+)'),
    }
if args.submission_command.split()[0] not in qsub_dict and args.hold_flag:
    global_def.ERROR('Qsub return output not known! Please contact the SPHIRE authors!', 'sxbatch.py', 1)

if args.first_hold_number:
    prev_hold = args.first_hold_number
else:
    prev_hold = 'aaa'
for idx, file_name in enumerate(sorted(glob.glob('{0}/*'.format(args.input_run_dir)))):
    command = args.submission_command.split()
    if args.hold_flag and (idx != 0 or args.first_hold_number):
        command.append(args.hold_flag)
        command.append(prev_hold)
    else:
        pass
    command.append(file_name)
    if args.hold_flag:
        print(' '.join(command))
        stdout = subprocess.check_output(command)
        print(stdout)
        prev_hold = qsub_dict[command[0]].match(stdout).group(1)
    else:
        subprocess.Popen(command).wait()

global_def.BATCH = False
