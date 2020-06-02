#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
import os
import glob
import re

bin_files = glob.glob('../sphire/bin/*.py')
lib_files = glob.glob('../sphire/libpy/*.py')


sequence = '(?:[^\w]|^)((?:{0})\.[\w]+)'.format('|'.join([
    os.path.splitext(os.path.basename(entry))[0].replace('sparx_', '')
    for entry in lib_files
    if 'sparx.py' not in entry
    ]))
SPARX_FUNC_RE = re.compile(sequence+'\s*\(')
SPARX_FUNC_RE_2 = re.compile(sequence)
sequence = '^import\s+({0})'.format('|'.join([
    os.path.splitext(os.path.basename(entry))[0].replace('sparx_', '')
    for entry in lib_files
    if 'sparx.py' not in entry
    ]))
SPARX_FUNC_RE_3 = re.compile(sequence)

ok_set = (
    'global_def.MPI',
    'global_def.BATCH',
    'global_def.interpolation_method_2D',
    'global_def.LOGFILE_HANDLE',
    'global_def.IS_LOGFILE_OPEN',
    'global_def.LOGFILE',
    'global_def.CACHE_DISABLE',
    'user_functions.factory',
    'alignment.ali_vol_func',
    'global_def.SPARXVERSION',
    'global_def.SPARX_MPI_TAG_UNIVERSAL',
    )

#for file_path in sorted(['../sphire/libpy/sparx_multi_shc.py']):
for file_path in sorted(bin_files + lib_files):
    print(file_path)
    with open(file_path, 'r') as read:
        lines = read.readlines()

    output_lines = []
    did_change = True
    while did_change:
        did_change = False
        for idx, line in enumerate(lines[:]):
            match_1 = SPARX_FUNC_RE.findall(line)
            match_2 = SPARX_FUNC_RE_2.findall(line)
            match_3 = SPARX_FUNC_RE_3.findall(line)
            if match_3:
                for entry in set(match_3):
                    lines[idx] = line.replace(entry, 'sparx_{0}'.format(entry))
            elif match_1 == match_2 and match_1:
                did_change = True
                for entry in set(match_1):
                    lines[idx] = line.replace(entry, 'sparx_{0}'.format(entry))
            elif match_1:
                did_change = True
                for entry in set(match_1):
                    lines[idx] = line.replace(entry, 'sparx_{0}'.format(entry))
            elif match_2:
                for entry in set(match_2):
                    if entry in ok_set:
                        lines[idx] = line.replace(entry, 'sparx_{0}'.format(entry))
                    elif entry.endswith('.py'):
                        continue
                    else:
                        print('NOT SURE ABOUT THIS!', idx, entry, line)
            else:
                lines[idx] = line

    with open(file_path, 'w') as write:
        write.write(''.join(lines))
