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
import importlib
import sys
import shutil
import re
import glob

from pyflakes import reporter as modReporter
from pyflakes.checker import Checker
import pyflakes.messages as pym
import pyflakes.api as pyfl

USED_FOLDER = ('bin', 'libpy', 'templates')
USED_FOLDER_EMAN2 = ('libpyEM', 'libpyEM/qtgui')
BIN_FOLDER = ('bin', 'templates')
LIB_FOLDER = ('libpy', 'templates', 'libpyEM', 'libpyEM/qtgui')

EMAN2_GUI_DICT = {}

IGNORE_FILES = (
    'sp_sparx.py',
    'sp_development.py',
    )
IGNORE_MODULES = (
    '__future__',
    'builtins',
    'future',
    'matplotlib.backends',
    'ipython',
    )
EXTERNAL_LIBS = (
    'EMAN2_cppwrap',
    'PyQt5.QtCore',
    'PyQt5.QtGui',
    )

IGNORE_LIST = (
    'sys',
    'random',
    'os',
    'traceback',
    'mpi',
    'collections',
    're',
    'six',
    'json',
    'EMAN2db',
    'numpy',
    'shutil',
    'logging',
    'argparse',
    'configparser',
    'user_functions',
    'sp_user_functions',
    'pickle',
    'time',
    'string',
    'subprocess',
    'fundamentals',
    'sp_fundamentals',
    'types',
    'EMAN2',
    'EMAN2_cppwrap',
    'operator',
    'errno',
    'global_def',
    'sp_global_def',
    )
IGNORE_RE_DICT = {}
for entry in IGNORE_LIST:
    IGNORE_RE_DICT[entry] = '(?:[^\w]|^){0}\.'.format(entry)

REPLACE_DICT = {
    'np': ['numpy', 'numpy'],
    'interpolate': ['scipy.interpolate', 'scipy.interpolate'],
    'stats': ['scipy.stats', 'scipy.stats'],
    'plt': ['matplotlib.pyplot', 'matplotlib.pyplot'],
    'scipy_spatial': ['scipy.spatial', 'scipy.spatial'],
    'pylab': ['matplotlib.pyplot', 'matplotlib.pyplot'],
    }

RE_DICT = {
    'strftime': [['(?:[^\w]|^)strftime\(', 'time']],
    'localtime': [['(?:[^\w]|^)localtime\(', 'time']],
    'random': [['(?:[^\w]|^)random\(', 'random']],
    'randint': [['(?:[^\w]|^)randint\(', 'random']],
    'seed': [['(?:[^\w]|^)seed\(', 'random']],
    'time': [['(?:[^\w]|^)time\(', 'time']],
    'Qt': [['(?:[^\w]|^)Qt\.', 'PyQt5.QtCore']],
    'QtCore': [['(?:[^\w]|^)QtCore\.', 'PyQt5']],
    'QtGui': [['(?:[^\w]|^)QtGui\.', 'PyQt5']],
    'fft': [['(?:[^\w]|^)fft\(', 'sp_fundamentals']],
    'pad': [['(?:[^\w]|^)pad\(', 'sp_utilities']],
    'path': [['(?:[^\w]|^)path\.(?:join|exists|isfile|realpath|basename)\(', 'os']],
    'Set': [['(?:[^\w]|^)Set\(', 'sets']],
    'copy': [['(?:[^\w]|^)copy\.(?:deepcopy|copy)\(', 'copy']],
    'array': [['(?:[^\w]|^)array\(', 'numpy']],
    'sin': [['(?:[^\w]|^)sin\(', 'numpy']],
    'tanh': [['(?:[^\w]|^)tanh\(', 'numpy']],
    'cos': [['(?:[^\w]|^)cos\(', 'numpy']],
    'log': [['(?:[^\w]|^)log\(', 'numpy']],
    'pi': [['(?:[^\w]|^)pi(?:[^\w]|^)', 'numpy']],
    'sqrt': [['(?:[^\w]|^)sqrt\(', 'numpy']],
    'radians': [['(?:[^\w]|^)radians\(', 'numpy']],
    'exp': [['(?:[^\w]|^)exp\(', 'numpy']],
    'copysign': [['(?:[^\w]|^)copysign\(', 'numpy']],
    'log': [['(?:[^\w]|^)log\(', 'numpy']],
    'ceil': [['(?:[^\w]|^)ceil\(', 'numpy']],
    'degrees': [['(?:[^\w]|^)degrees\(', 'numpy']],
    'floor': [['(?:[^\w]|^)floor\(', 'numpy']],
    'tan': [['(?:[^\w]|^)tan\(', 'numpy']],
    'log10': [['(?:[^\w]|^)log10\(', 'numpy']],
    'fmod': [['(?:[^\w]|^)fmod\(', 'numpy']],
    'argsort': [['(?:[^\w]|^)argsort\(', 'numpy']],
    'log2': [['(?:[^\w]|^)log2\(', 'sp_alignment']],
    'linalg': [['(?:[^\w]|^)linalg\.', 'numpy']],
    'zeros': [['(?:[^\w]|^)zeros\(', 'numpy']],
    'fmin': [['(?:[^\w]|^)fmin\(', 'scipy.optimize']],
    'split': [['(?:[^\w]|^)split\([\s\w]+(?:,\s*(?:\'|").*?|\s*)\)', 'string'], ['(?:[^\w]|^)split\((?:\'|").*?\)', 're']],
    'power': [['(?:[^\w]|^)power\(verr.*?\)', 'numpy'], ['(?:[^\w]|^)power\(periodogram.*?\)', 'sp_morphology'], ['(?:[^\w]|^)power\(EMAN2_cppwrap\.periodogram.*?\)', 'sp_morphology']],
    'square': [['(?:[^\w]|^)square\(', 'sp_morphology']],
    'loads': [['(?:[^\w]|^)loads\(', 'pickle']],
    'dumps': [['(?:[^\w]|^)dumps\(', 'pickle']],
    'compress': [['(?:[^\w]|^)compress\(', 'zlib']],
    }

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
NO_COMMENTS_DIR = os.path.join(CURRENT_DIR, '00_NO_COMMENTS')
NO_IMPORTS_DIR = os.path.join(CURRENT_DIR, '01_NO_IMPORTS')
try:
    shutil.rmtree(NO_IMPORTS_DIR) 
except OSError:
    pass
NO_WILDCARDS_DIR = os.path.join(CURRENT_DIR, '01_NO_WILDCARDS')
try:
    shutil.rmtree(NO_WILDCARDS_DIR) 
except OSError:
    pass

FUNCDEF_RE = re.compile('^(?:def|class)\s+([^\(]+)')
IMPORT_RE = re.compile('^(\s*)(import\s+([\w.]+)\s*.*|from\s+([\w.]+)\s+import.*)')
COMMENT_RE = re.compile("^\s*#")
BLOCK_STRING_SINGLE_RE = re.compile('\'\'\'')
BLOCK_STRING_DOUBLE_RE = re.compile('"""')
MULTILINE_RE = re.compile('.*\\\s*$')

FILE_NAME = None
PRINT_LINE = None
def my_to_list(self):
    lineno = int(self.lineno)
    col = int(self.col)
    name = re.match(".*'(.*)'", self.message % self.message_args).group(1)
    return (lineno, col, name)


def my_exception(self, filename, msg, lineno, offset, text):
    if msg == 'expected an indented block':
        print('ERROR:', filename, msg, lineno, offset, text.strip())
        self.indent_error.append([int(lineno), text])
    else:
        print('ERROR:', filename, msg, lineno, offset, text.strip())
        self.general_error.append([int(lineno), text])


ERRORS = {}
def my_report(self, messageClass, *args, **kwargs):
    if pym.UndefinedName == messageClass:
        message = messageClass(self.filename, *args, **kwargs)
        self.undefined_names.append(message.to_list())
    elif pym.UnusedImport == messageClass:
        if self.filename not in IGNORE_FILES:
            message = messageClass(self.filename, *args, **kwargs)
            self.unused_imports.append(message.to_list())
    elif pym.UnusedVariable == messageClass:
        message = messageClass(self.filename, *args, **kwargs)
        self.unused_var.append(message.to_list())
    elif pym.RedefinedInListComp == messageClass:
        message = messageClass(self.filename, *args, **kwargs)
        self.shadowed_var.append(message.to_list())
    elif pym.RedefinedWhileUnused == messageClass:
        message = messageClass(self.filename, *args, **kwargs)
        self.dublicated_funcs.append(message.to_list())
    else:
        if self.filename not in IGNORE_FILES:
            print(messageClass(self.filename, *args, **kwargs))
            ERRORS[str(messageClass)] = messageClass(self.filename, *args, **kwargs).to_list()


def reset_lists():
    GLOBAL_CHECKER.undefined_names = []
    GLOBAL_CHECKER.unused_imports = []
    GLOBAL_CHECKER.unused_var = []
    GLOBAL_CHECKER.shadowed_var = []
    GLOBAL_CHECKER.dublicated_funcs = []

modReporter.Reporter.syntaxError = my_exception
GLOBAL_REPORTER = modReporter._makeDefaultReporter()
GLOBAL_REPORTER.indent_error = []
GLOBAL_REPORTER.general_error = []

GLOBAL_CHECKER = Checker
GLOBAL_CHECKER.report = my_report
pym.Message.to_list = my_to_list


def print_all_info(**kwargs):
    if PRINT_LINE:
        if kwargs['line_idx'] in PRINT_LINE:
            for key in sorted(kwargs.keys()):
                print(key, ':', kwargs[key])
            print('')


def get_file_dict():
    file_dict = {}
    for folder_name in USED_FOLDER:
        files_dir = os.path.join(CURRENT_DIR, NO_COMMENTS_DIR, folder_name, '*.py')
        file_dict[folder_name] = sorted(glob.glob(files_dir))
    for folder_name in USED_FOLDER_EMAN2:
        files_dir = os.path.join(CURRENT_DIR, '../..', folder_name, '*.py*')
        file_dict[folder_name] = [entry for entry in sorted(glob.glob(files_dir)) if not entry.endswith('.pyc')]
        if 'qtgui' in folder_name:
            for file_path in file_dict[folder_name]:
                mod_name = os.path.splitext(os.path.basename(file_path))[0]
                EMAN2_GUI_DICT[mod_name] = 'eman2_gui.{0}'.format(mod_name)
    return file_dict


def get_external_libs(external_libs, lib_modules=None):
    if lib_modules is None:
        lib_modules = []
    lib_modules_ext = {}
    for entry in external_libs:
        if entry in lib_modules:
            continue
        try:
            importlib.import_module(entry)
        except:
            if 'pyqt' in entry.lower():
                continue
            else:
                raise
        lib_modules_ext[entry] = dir(sys.modules[entry])
    return lib_modules_ext


def get_library_funcs(file_dict):
    lib_modules = {}
    for key in file_dict:
        if key not in LIB_FOLDER:
            continue
        for file_name in file_dict[key]:
            name = os.path.splitext(os.path.basename(file_name))[0].replace('.py', '')
            with open(file_name) as read:
                lines = read.readlines()
            lib_modules[name] = [
                FUNCDEF_RE.match(entry).group(1)
                for entry in lines
                if FUNCDEF_RE.match(entry)
                ]
            if 'global_def' in name:
                lib_modules[name].append('IS_LOGFILE_OPEN')
                lib_modules[name].append('LOGFILE_HANDLE')
                lib_modules[name].append('SPARXVERSION')
                lib_modules[name].append('CACHE_DISABLE')
                lib_modules[name].append('SPARX_MPI_TAG_UNIVERSAL')
                lib_modules[name].append('interpolation_method_2D')
                lib_modules[name].append('Eulerian_Angles')
                lib_modules[name].append('BATCH')
                lib_modules[name].append('MPI')
                lib_modules[name].append('LOGFILE')
                lib_modules[name].append('SPARX_DOCUMENTATION_WEBSITE')
                lib_modules[name].append('ERROR')
                lib_modules[name].append('sxprint')
            elif 'EMAN2_meta' == name:
                lib_modules[name].append('EMANVERSION')
                lib_modules[name].append('DATESTAMP')
            elif 'EMAN2' == name:
                lib_modules[name].append('HOMEDB')
                lib_modules[name].append('T')
                lib_modules[name].append('bispec_invar_parm')
                lib_modules[name].append('outplaceprocs')
                lib_modules[name].append('this_file_dirname')
                lib_modules[name].append('GUIMode')
                lib_modules[name].append('app')
                lib_modules[name].append('GUIbeingdragged')
                lib_modules[name].append('originalstdout')
                lib_modules[name].append('to_numpy')
                lib_modules[name].append('from_numpy')
                lib_modules[name].append('file_mode_imap')
                lib_modules[name].append('file_mode_intmap')
                lib_modules[name].append('file_mode_range')
                lib_modules[name].append('good_box_sizes')
                lib_modules[name].append('parseparmobj1')
                lib_modules[name].append('parseparmobj2')
                lib_modules[name].append('parseparmobj3')
                lib_modules[name].append('parseparmobj4')
                lib_modules[name].append('parseparmobj_ob')
                lib_modules[name].append('parseparmobj_logical')
                lib_modules[name].append('parseparmobj_op_words')
                lib_modules[name].append('parseparmobj_logical_words')
                lib_modules[name].append('glut_inited')
            elif 'emapplication' in name:
                lib_modules[name].append('get_application')

    return lib_modules


def remove_imports(file_dict, lib_modules):
    local_imports = {}
    for key, file_names in file_dict.items():
        if key in USED_FOLDER_EMAN2:
            continue
        output_dir = os.path.join(NO_IMPORTS_DIR, key)
        try:
            os.makedirs(output_dir)
        except OSError:
            pass
        for file_path in file_names:
            basename = os.path.basename(file_path)
            output_file_name = os.path.join(output_dir, basename)
            with open(file_path) as read:
                lines = read.readlines()

            local_modules = []
            if basename not in IGNORE_FILES:
                comment1 = False
                comment2 = False
                multiline_idx = []
                for line_idx, line in enumerate(lines[:]):
                    if COMMENT_RE.match(line):
                        continue
                    match = IMPORT_RE.match(line)
                    match_block_1 = BLOCK_STRING_SINGLE_RE.findall(line)
                    match_block_2 = BLOCK_STRING_DOUBLE_RE.findall(line)
                    do_final = False
                    if match_block_1 and not comment1 and not comment2:
                        if len(match_block_1) % 2 == 1:
                            comment1 = True
                    elif match_block_2 and not comment1 and not comment2:
                        if len(match_block_2) % 2 == 1:
                            comment2 = True
                    elif match_block_1:
                        comment1 = False
                        do_final = True
                    elif match_block_2:
                        comment2 = False
                        do_final = True

                    if match and not comment1 and not comment2 and not do_final:
                        indent = match.group(1)
                        content = match.group(2)
                        module = None
                        if match.group(3) is not None:
                            module = match.group(3)
                        elif match.group(4) is not None:
                            module = match.group(4)
                        assert module is not None
                        do_continue = False
                        for ignored_module in IGNORE_MODULES:
                            if ignored_module.lower() in module.lower():
                                do_continue = True
                                break
                        if do_continue:
                            continue
                        if ';' in content:
                            split_content = content.split(';')
                            new_line = '{0}{1}; pass#IMPORTIMPORTIMPORT {2}\n'.format(
                                indent,
                                '; '.join([entry.strip() for entry in split_content[1:]]),
                                split_content[0].strip(),
                                )
                        else:
                            new_line = '{0}pass#IMPORTIMPORTIMPORT {1}\n'.format(
                                indent,
                                content.strip(),
                                )
                        lines[line_idx] = new_line
                        idx = 0
                        while MULTILINE_RE.match(lines[line_idx+idx]):
                            idx += 1
                            multiline_idx.append(line_idx+idx)
                        local_modules.append(module)
                for idx in multiline_idx:
                    lines[idx] = '\n'
                try:
                    local_imports[basename] = get_external_libs(set(local_modules), lib_modules)
                except:
                    print(file_path)
                    raise
            else:
                local_imports[basename] = {}

            with open(output_file_name, 'w') as write:
                write.write(''.join(lines))
    return local_imports


def identify_missing(file_dict):
    local_imports = {}
    for key, file_names in file_dict.items():
        if key in USED_FOLDER_EMAN2:
            local_imports[key] = []
            continue
        input_dir = os.path.join(NO_IMPORTS_DIR, key)
        for file_path in file_names:
            basename = os.path.basename(file_path)
            input_file_path = os.path.join(input_dir, basename)
            reset_lists()
            with open(input_file_path) as read:
                file_content = read.read()
            pyfl.check(file_content, basename, GLOBAL_REPORTER)
            local_imports[basename] = GLOBAL_CHECKER.undefined_names

    return local_imports


def index_search(lines, index, no_index):
    if lines[index].strip().endswith('\\'):
        no_index.append(index+1)
        index_search(lines, index+1, no_index)


def fix_missing(file_dict, missing_modules_local, lib_modules, lib_modules_ext, lib_modules_local):
    for key, file_names in file_dict.items():
        if key in USED_FOLDER_EMAN2:
            continue
        output_dir = os.path.join(NO_WILDCARDS_DIR, key)
        input_dir = os.path.join(NO_IMPORTS_DIR, key)
        try:
            os.makedirs(output_dir)
        except OSError:
            pass
        for file_path in file_names:
            local_imports = []
            basename = os.path.basename(file_path)
            output_file_path = os.path.join(output_dir, basename)
            input_file_path = os.path.join(input_dir, basename)
            with open(input_file_path) as read:
                lines = read.readlines()
            for missing in reversed(sorted(missing_modules_local[basename])):
                matches = []
                dicts = [lib_modules, lib_modules_ext, lib_modules_local[basename]]
                for dictionary in dicts:
                    for key, values in dictionary.items():
                        for val in values:
                            if val == missing[2]:
                                matches.append(key)
                ignore = False
                replace = False
                if 'sp_development' in matches:
                    matches.remove('sp_development')
                if 'PyQt5.Qt' in matches:
                    matches.remove('PyQt5.Qt')
                for entry in matches[:]:
                    if 'eman2_gui.' in entry:
                        matches.remove(entry)
                    elif 'sparx' in entry:
                        matches.remove(entry)

                search_line = lines[missing[0]-1][max(missing[1]-1, 0):]
                if len(matches) == 1 or len(set(matches)) == 1:
                    used_module = matches[0]
                else:
                    try:
                        match_ignore = re.match(IGNORE_RE_DICT[missing[2]], search_line)
                    except KeyError:
                        match_ignore = None
                    try:
                        re_idx = None
                        for idx, value in enumerate(RE_DICT[missing[2]]):
                            match_re = re.match(value[0], search_line)
                            if match_re:
                                re_idx = idx
                                break
                    except KeyError:
                        match_re = None

                    if match_ignore:
                        used_module = missing[2]
                        ignore = True
                    elif 'EMAN2_cppwrap' in matches:
                        used_module = 'EMAN2_cppwrap'
                    elif match_re:
                        used_module = RE_DICT[missing[2]][re_idx][1]
                    elif missing[2] in REPLACE_DICT:
                        used_module = REPLACE_DICT[missing[2]][1]
                        replace = REPLACE_DICT[missing[2]][0]
                    elif matches:
                        print('CONFUSION', basename, missing, matches)
                        continue
                    else:
                        #print('TYPO', basename, missing, matches)
                        continue

                local_imports.append(used_module)
                if used_module in EMAN2_GUI_DICT:
                    used_module = EMAN2_GUI_DICT[used_module]
                if not ignore and not replace:
                    lines[missing[0]-1] = '{0}{1}.{2}'.format(
                        lines[missing[0]-1][:missing[1]],
                        used_module,
                        lines[missing[0]-1][missing[1]:]
                        )
                elif replace:
                    lines[missing[0]-1] = '{0}{1}.{2}'.format(
                        lines[missing[0]-1][:missing[1]],
                        used_module,
                        lines[missing[0]-1][missing[1] + len(missing[2]) + 1:]
                        )

            with open(output_file_path, 'w') as write:
                write.write(''.join(lines))
            with open(output_file_path + '_imports', 'w') as write:
                write.write(''.join(['{0}\n'.format(entry) for entry in sorted(list(set(local_imports)))]))


def main():
    file_dict = get_file_dict()

    #Extract all function names from library files
    lib_modules = get_library_funcs(file_dict)

    # Look for content of external_modules
    lib_modules_ext = get_external_libs(EXTERNAL_LIBS)

    # Remove all imports from files
    lib_modules_local = remove_imports(file_dict, lib_modules)

    missing_modules_local = identify_missing(file_dict)

    fix_missing(file_dict, missing_modules_local, lib_modules, lib_modules_ext, lib_modules_local)

if __name__ == '__main__':
    main()
