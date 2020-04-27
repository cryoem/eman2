#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

# Author: Markus Stabrin, 2018/11/20 (markus.stabrin@mpi-dortmund.mpg.de
# Copyright (c) 2018 Max Plank of molecular physiology Dortmund
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA

import os
import shutil
import sys
import glob
import re
import importlib
import argparse

from pyflakes import reporter as modReporter
from pyflakes.checker import Checker
import pyflakes.messages as pym
import pyflakes.api as pyfl

parser = argparse.ArgumentParser()
parser.add_argument('--silent', action='store_true', help='Do not write any output to disc')
parser.add_argument('--leave_imports', action='store_true', help='Do not clear tmp import lines to not touch the original structure.')
options = parser.parse_args()


def my_report(self, messageClass, *args, **kwargs):
    if pym.UndefinedName == messageClass:
        message = messageClass(self.filename, *args, **kwargs)
        self.undefined_names.append(message.to_list())
    elif pym.UnusedImport == messageClass:
        if 'sparx.py' not in self.filename:
            message = messageClass(self.filename, *args, **kwargs)
            self.unused_imports.append(message.to_list())


def my_to_list(self):
    lineno = int(self.lineno) - 1
    col = int(self.col)
    name = re.match(".*'(.*)'", self.message % self.message_args).group(1)
    return (lineno, col, name)


def my_exception(self, filename, msg, lineno, offset, text):
    if msg == 'expected an indented block':
        print(filename, msg, lineno, offset, text.strip())
        self.undefined_names.append([int(lineno)-1, text])


def index_search(lines, index, no_index):
    if lines[index].strip().endswith('\\'):
        no_index.append(index+1)
        index_search(lines, index+1, no_index)


folders = ['tmp', 'no_import', 'new']
for entry in folders:
    try:
        shutil.rmtree(entry)
    except Exception as e:
        print(e)
    try:
        os.mkdir(entry)
    except Exception as e:
        print(e)

IMPORT_DEF_RE = re.compile(r'^(?:def|class) ([^(]*).*')
IMPORT_IMPORTS_RE = re.compile(r'^(\s*).*from\s+([\w.]*)\s+import\s+([\w.]*)\s*(?:as\s+([\w.]*).*|)')
IMPORT_SINGLE_IMPORT_RE = re.compile(r'^(\s*)(?:from\s+[\w.]\s+|)import\s+([\w.,\s]*)\s*(?:as\s*([\w.]*).*|)')
IMPORT_LEADING_RE = re.compile("^(\s*)[^\s]*.*")
IMPORT_COMMENT_RE = re.compile("^(\s*)#[^\s]*.*")
IMPORT_MATPLOTLIB_RE = re.compile("^(\s*)matplotlib.use")
IMPORT_FIND_SYNTAX_RE = re.compile("^\s*(if|else|try|except).*:")


lib_files = sorted(glob.glob('../sphire/libpy/*.py'))
lib_files_2 = sorted(glob.glob('../sphire/templates/sxgui_template.py'))
lib_eman2_files = sorted(glob.glob('../libpyEM/*.py'))
lib_eman2_files_2 = sorted(glob.glob('../libpyEM/qtgui/*.py'))
lib_eman2_files_3 = sorted(glob.glob('../libpyEM/*.py.in'))
bin_files = sorted(glob.glob('../sphire/bin/*.py'))
bin_files_2 = sorted(glob.glob('../sphire/templates/*.py'))

qtgui_files = [os.path.splitext(os.path.basename(entry))[0] for entry in lib_eman2_files_2]

transform_dict = {}

lib_modules = {}
lib_modules_ext = {}
for lib_file in lib_files + lib_files_2 + lib_eman2_files + lib_eman2_files_2 + lib_eman2_files_3:
    name = os.path.splitext(os.path.basename(lib_file))[0].replace('.py', '')
    with open(lib_file) as read:
        lib_modules[name] = [
            IMPORT_DEF_RE.match(entry).group(1)
            for entry in read.readlines()
            if IMPORT_DEF_RE.match(entry)
            ]
    if 'sparx_global_def' == name:
        lib_modules[name].append('SPARXVERSION')
        lib_modules[name].append('CACHE_DISABLE')
        lib_modules[name].append('SPARX_MPI_TAG_UNIVERSAL')
        lib_modules[name].append('interpolation_method_2D')
        lib_modules[name].append('Eulerian_Angles')
        lib_modules[name].append('BATCH')
        lib_modules[name].append('MPI')
        lib_modules[name].append('LOGFILE')
        lib_modules[name].append('SPARX_DOCUMENTATION_WEBSITE')
    elif 'global_def' == name:
        lib_modules[name].append('SPARXVERSION')
        lib_modules[name].append('CACHE_DISABLE')
        lib_modules[name].append('SPARX_MPI_TAG_UNIVERSAL')
        lib_modules[name].append('interpolation_method_2D')
        lib_modules[name].append('Eulerian_Angles')
        lib_modules[name].append('BATCH')
        lib_modules[name].append('MPI')
        lib_modules[name].append('LOGFILE')
        lib_modules[name].append('SPARX_DOCUMENTATION_WEBSITE')
    elif 'EMAN2_meta' == name:
        lib_modules[name].append('EMANVERSION')
        lib_modules[name].append('DATESTAMP')
    elif 'emapplication' in name:
        lib_modules[name].append('get_application')

    if name.startswith('sparx_'):
        transform_dict[name.replace('sparx_', '')] = name


external_modules = [
    'subprocess',
    'random',
    'mpi',
    'sys',
    'traceback',
    'time',
    'numpy',
    'numpy.random',
    'math',
    'operator',
    'optparse',
    'EMAN2_cppwrap',
    'sets',
    'copy',
    'inspect',
    'scipy',
    'scipy.optimize',
    'scipy.stats',
    'datetime',
    'heapq',
    'matplotlib',
    'os',
    'string',
    'builtins',
    'shutil',
    'glob',
    'types',
    'pickle',
    'zlib',
    'struct',
    'fractions',
    'socket',
    ]
for entry in external_modules:
    importlib.import_module(entry)
    lib_modules_ext[entry] = dir(sys.modules[entry])

python_files = bin_files + lib_files + bin_files_2

modReporter.Reporter.syntaxError = my_exception
reporter = modReporter._makeDefaultReporter()
reporter.undefined_names = []

Checker.report = my_report
Checker.undefined_names = []
Checker.unused_imports = []
pym.Message.to_list = my_to_list


#python_files = glob.glob('../sphire/bin/sxchains.py')
#python_files = glob.glob('../sphire/libpy/applications.py')
rounds = 0
while True:
    rounds += 1
    ok = 0
    replace = 0
    fatal = [0, []]
    confusion = [0, []]
    syntax = [0, []]
    for file_name in python_files:
        if 'sparx.py' in file_name:
            continue
        print('\n\n\n######################################\n\n\n')
        print(file_name)

        with open(file_name, 'r') as read:
            lines = read.readlines()

        bad_index = []
        no_index = []
        file_modules = []
        local_func_import = {}
        for index, entry in enumerate(lines):
            match = IMPORT_IMPORTS_RE.match(entry)
            if IMPORT_COMMENT_RE.match(entry):
                continue

            elif match and \
                    'future' not in entry and \
                    'qt' not in entry.lower() and \
                    not 'builtins' in entry:
                string = match.group(1)
                local_func_import.setdefault(match.group(2), []).append(match.group(3))
                file_modules.append([index, string, [match.group(2)]])
                if match.group(4):
                    try:
                        lib_modules[match.group(2)].append(match.group(4))
                    except KeyError:
                        lib_modules_ext[match.group(2)].append(match.group(4))
                if match.group(2) in external_modules:
                    try:
                        lib_modules[match.group(2)].append(match.group(3))
                    except KeyError:
                        lib_modules_ext[match.group(2)].append(match.group(3))
                bad_index.append([index, string])
                if lines[index].strip().endswith('\\'):
                    no_index.append(index+1)
                    index_search(lines, index+1, no_index)

            else:
                match = IMPORT_SINGLE_IMPORT_RE.match(entry)
                if match and \
                        'future' not in entry and \
                        'qt' not in entry.lower() and \
                        not 'builtins' in entry:
                    string = match.group(1)
                    name = [loc_entry.strip() for loc_entry in match.group(2).split(',')]
                    file_modules.append([index, string, name])
                    if lines[index].strip().endswith('\\'):
                        no_index.append(index+1)
                        index_search(lines, index+1, no_index)

            match = IMPORT_MATPLOTLIB_RE.match(entry)
            if match:
                string = match.group(1)
                no_index.append(index)


        no_from_import_lines = lines[:]
        no_import_lines = lines[:]

        for idx, string in bad_index:
            if '#IMPORTIMPORTIMPORT' in lines[idx]:
                continue
            elif lines[idx].strip().startswith('#'):
                find_multiline = '\n{0}'.format(string).join([entry.strip() for entry in lines[idx].split(';')])
                no_from_import_lines[idx] = '#{0}pass#IMPORTIMPORTIMPORT {1}\n'.format(string, find_multiline)
                no_import_lines[idx] = '#{0}pass#IMPORTIMPORTIMPORT {1}\n'.format(string, find_multiline)
            else:
                find_multiline = '\n{0}'.format(string).join([entry.strip() for entry in lines[idx].split(';')])
                no_from_import_lines[idx] = '{0}pass#IMPORTIMPORTIMPORT {1}\n'.format(string, find_multiline)
                no_import_lines[idx] = '{0}pass#IMPORTIMPORTIMPORT {1}\n'.format(string, find_multiline)

        for idx in no_index:
            if '#IMPORTIMPORTIMPORT' in lines[idx]:
                continue
            else:
                string = IMPORT_LEADING_RE.match(lines[idx]).group(1)
                find_multiline = '\n{0}'.format(string).join([entry.strip() for entry in lines[idx].split(';')])
                no_from_import_lines[idx] = '#{0}pass#IMPORTIMPORTIMPORT {1}\n'.format(string, find_multiline)
                no_import_lines[idx] = '#{0}pass#IMPORTIMPORTIMPORT {1}\n'.format(string, find_multiline)

        correct_imports = []
        for idx, string, module in file_modules:
            correct_imports.extend([entry for entry in module if entry not in transform_dict])
            if '#IMPORTIMPORTIMPORT' in lines[idx]:
                continue
            elif lines[idx].strip().startswith('#'):
                find_multiline = '\n{0}'.format(string).join([entry.strip() for entry in lines[idx].split(';')])
                no_import_lines[idx] = '#{0}pass#IMPORTIMPORTIMPORT {1}\n'.format(string, find_multiline)
            else:
                find_multiline = '\n{0}'.format(string).join([entry.strip() for entry in lines[idx].split(';')])
                no_import_lines[idx] = '{0}pass#IMPORTIMPORTIMPORT {1}\n'.format(string, find_multiline)
        correct_imports = list(set(correct_imports))

        while True:
            stop = False
            need_intervention = False
            Checker.undefined_names = []
            Checker.unused_imports = []
            reporter.undefined_names = []
            file_content = ''.join(no_from_import_lines)
            pyfl.check(file_content, file_name, reporter)

            if not options.silent:
                with open(os.path.join('tmp', os.path.basename(file_name)), 'w') as write:
                    write.write(file_content)

            if not options.silent:
                file_content = ''.join(no_import_lines)
                with open(os.path.join('no_import', os.path.basename(file_name)), 'w') as write:
                    write.write(file_content)

            if not reporter.undefined_names:
                stop = True
            else:
                maxi = 30
                for start_num, text in reporter.undefined_names:
                    stop_idx = 0
                    for i in range(1, maxi):
                        line = no_from_import_lines[start_num-i]
                        match = IMPORT_FIND_SYNTAX_RE.match(line)
                        if match:
                            stop_idx = i
                            break
                    if match.group(1) == 'if':
                        no_import_lines[start_num - stop_idx] = '\n'
                        no_from_import_lines[start_num - stop_idx] = '\n'
                        if 'else:' in text:
                            no_import_lines[start_num] = '\n'
                            no_from_import_lines[start_num] = '\n'
                    else:
                        stop = True
                        need_intervention = True
            if stop:
                if need_intervention:
                    syntax[0] += 1
                    syntax[1].append(file_name)
                    print('Needs manual intervention!')
                break
            else:
                print('Resolved!')

        fatal_list = []
        ok_list = []
        confusion_list = []
        replace_list = []
        for line_number, column, name in sorted(Checker.undefined_names):
            mod_list = []
            for key, values in lib_modules.items():
                for val in values:
                    if name == val:
                        mod_list.append(key)
            mod_list = list(set(mod_list))
            if len(mod_list) == 1:
                out_list = ok_list
            else:
                for key, values in lib_modules_ext.items():
                    for val in values:
                        if name == val:
                            mod_list.append(key)
                mod_list = list(set(mod_list))

                if name == 'os':
                    mod_list = ['os']
                    out_list = replace_list
                elif name == 'random':
                    mod_list = ['random', 'numpy']
                    out_list = confusion_list
                elif name == 'copy':
                    mod_list = ['copy']
                    out_list = replace_list
                elif name in transform_dict:
                    mod_list = [transform_dict[name]]
                    out_list = replace_list
                elif not mod_list:
                    out_list = fatal_list
                elif len(mod_list) == 1:
                    out_list = ok_list
                #elif 'numpy' in mod_list:
                #    mod_list = ['numpy']
                #    out_list = ok_list
                elif ['math', 'numpy', 'scipy'] == sorted(mod_list):
                    mod_list = ['numpy']
                    out_list = ok_list
                elif ['numpy', 'scipy'] == sorted(mod_list):
                    mod_list = ['numpy']
                    out_list = ok_list
                elif ['numpy', 'scipy', 'scipy.optimize'] == sorted(mod_list):
                    mod_list = ['numpy']
                    out_list = ok_list
                else:
                    local_imports = []
                    for key in local_func_import:
                        for module_name in local_func_import[key]:
                            if module_name == name:
                                local_imports.append(key)
                    local_imports = list(set(local_imports))

                    if len(local_imports) == 1:
                        mod_list = local_imports
                        out_list = ok_list
                    else:
                        local_imports = []
                        for entry in file_modules:
                            for module_name in entry[2]:
                                if module_name == name:
                                    local_imports.append(module_name)
                        local_imports = list(set(local_imports))

                        if len(local_imports) == 1:
                            mod_list = local_imports
                            out_list = ok_list
                        else:
                            #print(name, local_imports)
                            out_list = confusion_list
            out_list.append([line_number, column, name, mod_list])

        print('Typos that needs to be resolved:')
        template = 'name: {2:>25s}, line: {0: 6d}, column: {1: 6d}, module(s): {3}'
        fatal[0] += len(fatal_list)
        if len(fatal_list):
            fatal[1].append(file_name)
        for entry in fatal_list:
            print(template.format(*entry))

        print('')
        print('Replace list:')
        used_modules = []
        replace += len(replace_list)
        idx_line = 0
        idx_column = 1
        idx_name = 2
        idx_mod = 3
        for entry in replace_list:
            print(template.format(*entry))
            used_modules.extend(entry[idx_mod])
        print('')

        print('')
        print('Confusion list:')
        confusion[0] += len(confusion_list)
        if len(confusion_list):
            confusion[1].append(file_name)
        for entry in confusion_list:
            print(template.format(*entry))
        print('')

        print('RESOLVED THINGS:')
        ok += len(ok_list)
        len_dict = {}
        for entry in ok_list:
            print(template.format(*entry))
            used_modules.extend(entry[idx_mod])
            out = []

            new_line = []
            current_line = no_import_lines[entry[idx_line]]
            len_name = len(entry[idx_name])
            len_mod = len(entry[idx_mod][0])+1
            len_adjust = len_dict.setdefault(entry[idx_line], 0)

            new_line.append(current_line[:entry[idx_column] + len_adjust])
            new_line.append('{0}.{1}'.format(entry[idx_mod][0], entry[idx_name]))
            new_line.append(current_line[entry[idx_column]+len_name+len_adjust:])
            no_import_lines[entry[idx_line]] = ''.join(new_line)
            try:
                len_dict[entry[idx_line]] += len_mod
            except KeyError:
                len_dict[entry[idx_line]] = len_mod

        correct_imports_clean = []
        for entry in correct_imports:
            try:
                importlib.import_module(entry.split()[0])
            except ImportError:
                if 'eman2_gui' in entry:
                    pass
                elif 'sparx_' in entry:
                    pass
                elif 'sxgui_template' in entry:
                    pass
                else:
                    continue
            try:
                lib = transform_dict[entry]
            except KeyError:
                lib = entry
            correct_imports_clean.append(lib)
        correct_imports_clean = list(set(correct_imports_clean))

        used_modules_clean = []
        for entry in used_modules:
            try:
                lib = transform_dict[entry]
            except KeyError:
                lib = entry
            used_modules_clean.append(lib)
        used_modules = list(set(used_modules_clean))


        if not replace_list and not ok_list:
            print('Remove imports', Checker.unused_imports)
            imports_to_remove = []
            for _, _, name in Checker.unused_imports:
                if name not in ('matplotlib'):
                    imports_to_remove.append(name)
            imports = ['import {0}\n'.format(entry) if entry not in qtgui_files else 'import eman2_gui.{0} as {0}\n'.format(entry) for entry in list(set(used_modules)) if entry not in imports_to_remove]
            imports.extend(['import {0}\n'.format(entry) if entry.split('.')[-1] not in qtgui_files else 'import eman2_gui.{0} as {0}\n'.format(entry.split('.')[-1]) for entry in correct_imports_clean if entry not in imports_to_remove])
        else:
            imports = ['import {0}\n'.format(entry) if entry not in qtgui_files else 'import eman2_gui.{0} as {0}\n'.format(entry) for entry in list(set(used_modules))]
            imports.extend(['import {0}\n'.format(entry) if entry.split('.')[-1] not in qtgui_files else 'import eman2_gui.{0} as {0}\n'.format(entry.split('.')[-1]) for entry in correct_imports_clean])
        imports = sorted(list(set(imports)))
        inserted = False
        for idx, entry in enumerate(imports[:]):
            if entry == 'import matplotlib\n' and not inserted:
                imports.insert(idx+1, 'matplotlib.use("Agg")\n')
                inserted = True
                break
            elif 'matplotlib' in entry and not inserted:
                imports.insert(idx, 'matplotlib.use("Agg")\n')
                imports.insert(idx, 'import matplotlib\n')
                inserted = True
                break

        imports = ''.join(imports)
        first_1 = False
        first_2 = False
        first_3 = False
        index = 3
        for idx, line in enumerate(no_import_lines[2:], 2):
            if line.startswith("class") or line.startswith('def'):
                    if not first_1 and not first_2 and not first_3:
                        index = idx
                        break
            elif line.startswith("'''"):
                if first_1:
                    index = idx+1
                    break
                else:
                    first_1 = True
            elif line.startswith('"""'):
                if first_2:
                    index = idx+1
                    break
                else:
                    first_2 = True
            elif line.startswith('#') and not first_1 and not first_2:
                first_3 = True
            elif first_3:
                index = idx+1
                break
        no_import_lines.insert(index, imports)

        remove_indices = []
        first_1 = False
        first_2 = False
        for idx, line in enumerate(no_import_lines[2:], 2):
            if line.startswith("'''"):
                if first_1:
                    first_1 = False
                else:
                    first_1 = True
            elif line.startswith('"""'):
                if first_2:
                    first_2 = False
                else:
                    first_2 = True
            elif line.startswith('#'):
                pass
            elif line.startswith("class") or line.startswith('def'):
                if not first_1 and not first_2:
                    pass
            if '#IMPORTIMPORTIMPORT' in line and not options.leave_imports:
                remove_indices.append(idx)


        output_lines = []
        for idx, line in enumerate(no_import_lines):
            if idx in remove_indices:
                pass
            else:
                output_lines.append(line)

        if not options.silent:
            file_content = ''.join(output_lines)
            with open(os.path.join('new', os.path.basename(file_name)), 'w') as write:
                write.write(file_content)
            with open(file_name, 'w') as write:
                write.write(file_content)
        print('')

        if not ok_list  and replace_list:
            print('REPLACE')
            len_dict = {}
            for entry in replace_list:
                print(template.format(*entry))
                out = []

                new_line = []
                current_line = no_import_lines[entry[idx_line]+1]
                len_name = len(entry[idx_name])
                len_mod = len(entry[idx_mod][0])
                len_adjust = len_dict.setdefault(entry[idx_line]+1, 0)

                new_line.append(current_line[:entry[idx_column] + len_adjust])
                new_line.append(entry[idx_mod][0])
                new_line.append(current_line[entry[idx_column] + len_name + len_adjust:])
                no_import_lines[entry[idx_line]+1] = ''.join(new_line)
                try:
                    len_dict[entry[idx_line]+1] += len_mod - len_name
                except KeyError:
                    len_dict[entry[idx_line]+1] = len_mod - len_name
            output_lines = []
            for idx, line in enumerate(no_import_lines):
                if idx in remove_indices:
                    pass
                else:
                    output_lines.append(line)

            file_content = ''.join(output_lines)
            with open(os.path.join('new', os.path.basename(file_name)), 'w') as write:
                write.write(file_content)
            with open(file_name, 'w') as write:
                write.write(file_content)

    print('FATAL:', fatal)
    print('CONFUSION:', confusion)
    print('RESOLVED:', ok)
    print('REPLACE:', replace)
    print('SYNTAX:', syntax)
    if options.silent:
        print('Resolved after', rounds, 'rounds')
        break
    elif ok == 0 and replace == 0:
        print('Resolved after', rounds, 'rounds')
        break
