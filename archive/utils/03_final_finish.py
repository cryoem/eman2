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

IGNORE_FILES = (
    'sp_sparx.py',
    )

USED_FOLDER = ('bin', 'libpy', 'templates')
FUNCDEF_RE = re.compile('^(?:def|class)\s+([^\(]+)')
EMPTY_RE = re.compile("^\s*$")
BAD_RE = re.compile("(?:pass#IMPORTIMPORTIMPORT|#MULTILINEMULTILINEMULTILINE)")
IF_RE = re.compile("^\s*if .*:")
ELSE_RE = re.compile("^\s*else\s*:")
LINE_RE = re.compile("^\s*\w")

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))

NO_USED = os.path.join(CURRENT_DIR, '02_NO_USED')
NO_UNUSED = os.path.join(CURRENT_DIR, '02_NO_UNUSED')
NO_WILDCARDS_DIR = os.path.join(CURRENT_DIR, '01_NO_WILDCARDS')

IMPORTS_SPHIRE_DIR = os.path.join(CURRENT_DIR, '03_SPHIRE_IMPORTS')
SPHIRE_DIR = os.path.join(CURRENT_DIR, 'SPHIRE')
try:
    shutil.rmtree(IMPORTS_SPHIRE_DIR)
except:
    pass
try:
    shutil.rmtree(SPHIRE_DIR)
except:
    pass

from pyflakes import reporter as modReporter
from pyflakes.checker import Checker
import pyflakes.messages as pym
import pyflakes.api as pyfl

def my_to_list(self):
    lineno = int(self.lineno)
    col = int(self.col)
    name = re.match(".*'(.*)'", self.message % self.message_args).group(1)
    return (lineno, col, name)


def my_exception(self, filename, msg, lineno, offset, text):
    if msg == 'expected an indented block':
        #print(filename, msg, lineno, offset, text.strip())
        self.indent_error.append([int(lineno), text])
    else:
        print('WOHOO')
        print(filename, msg, lineno, offset, text.strip())


def my_report(self, messageClass, *args, **kwargs):
    message = messageClass(self.filename, *args, **kwargs)
    if self.filename not in IGNORE_FILES:
        self.my_error_dict.setdefault(str(messageClass), []).append(message.to_list())


def reset_lists():
    GLOBAL_CHECKER.my_error_dict = {}
    GLOBAL_REPORTER.indent_error = []

modReporter.Reporter.syntaxError = my_exception
GLOBAL_REPORTER = modReporter._makeDefaultReporter()

GLOBAL_CHECKER = Checker
GLOBAL_CHECKER.report = my_report
GLOBAL_CHECKER.my_error_dict = {}
pym.Message.to_list = my_to_list


USED_FOLDER_EMAN2 = ('libpyEM/qtgui', )
EMAN2_GUI_DICT = {}
def get_file_dict():
    file_dict = {}
    for folder_name in USED_FOLDER:
        #files_dir = os.path.join(CURRENT_DIR, NO_USED, folder_name, '*.py')
        #file_dict.setdefault('no_used', {})[folder_name] = sorted(glob.glob(files_dir))
        files_dir = os.path.join(CURRENT_DIR, NO_UNUSED, folder_name, '*.py')
        file_dict.setdefault('no_unused', {})[folder_name] = sorted(glob.glob(files_dir))
    for folder_name in USED_FOLDER_EMAN2:
        files_dir = os.path.join(CURRENT_DIR, '../..', folder_name, '*.py')
        for file_path in glob.glob(files_dir):
            mod_name = os.path.splitext(os.path.basename(file_path))[0]
            EMAN2_GUI_DICT[mod_name] = 'eman2_gui.{0}'.format(mod_name)
    return file_dict


def put_imports(file_dict):
    for key in file_dict:
        #if key == 'no_used':
        #    output_base = IMPORTS_SPARX_DIR
        #    second = True
        #else:
        output_base = IMPORTS_SPHIRE_DIR
        second = False
        for dir_name in file_dict[key]:
            output_dir = os.path.join(output_base, dir_name)
            try:
                os.makedirs(output_dir)
            except OSError:
                pass
            for file_path in file_dict[key][dir_name]:
                basename = os.path.basename(file_path)
                imports_file = '{0}_imports'.format(file_path.replace(NO_USED, NO_WILDCARDS_DIR).replace(NO_UNUSED, NO_WILDCARDS_DIR))
                with open(imports_file, 'r') as read:
                    imports = ''.join(['import {0}'.format(entry) if entry.strip() not in EMAN2_GUI_DICT else 'import {0}\n'.format(EMAN2_GUI_DICT[entry.strip()]) for entry in read.readlines()])
                with open(file_path, 'r') as read:
                    lines = read.readlines()
                if second:
                    lines.insert(1, imports)
                else:
                    for idx, line in enumerate(lines[:]):
                        if ('import' in line and '__future__' not in line) or FUNCDEF_RE.match(line):
                            lines.insert(idx, imports)
                            break
                with open(os.path.join(output_dir, basename), 'w') as write:
                    write.write(''.join(lines))


def remove_unused(file_dict):
    for key in file_dict:
        #if key == 'no_used':
        #    input_base = IMPORTS_SPARX_DIR
        #    output_base = SPARX_DIR
        #else:
        input_base = IMPORTS_SPHIRE_DIR
        output_base = SPHIRE_DIR
        for dir_name in file_dict[key]:
            input_dir = os.path.join(input_base, dir_name)
            output_dir = os.path.join(output_base, dir_name)
            try:
                os.makedirs(output_dir)
            except OSError:
                pass
            for file_path in file_dict[key][dir_name]:
                basename = os.path.basename(file_path)
                input_file_path = os.path.join(input_dir, basename)
                output_file_path = os.path.join(output_dir, basename)

                with open(input_file_path, 'r') as read:
                    file_content = read.read()
                reset_lists()
                pyfl.check(file_content, input_file_path, GLOBAL_REPORTER)
                bad_idx = []
                try:
                    unused_imports = GLOBAL_CHECKER.my_error_dict[str(pym.UnusedImport)]
                except KeyError:
                    unused_imports = []
                for lineno, _, _ in unused_imports:
                    bad_idx.append(lineno-1)
                lines = file_content.splitlines(1)
                for idx, line in enumerate(lines):
                    if FUNCDEF_RE.match(line):
                        curr_idx = idx
                        while EMPTY_RE.match(lines[curr_idx-1]):
                            bad_idx.append(curr_idx-1)
                            curr_idx -= 1

                for idx in reversed(sorted(bad_idx)):
                    del lines[idx]

                for idx, line in reversed(list(enumerate(lines[:]))):
                    if FUNCDEF_RE.match(line):
                        lines.insert(idx, '\n')
                        lines.insert(idx, '\n')
                    elif BAD_RE.search(line):
                        del lines[idx]

                reset_lists()
                pyfl.check(''.join(lines), output_file_path, GLOBAL_REPORTER)
                bad_idx = []
                while GLOBAL_REPORTER.indent_error:
                    do_break = False
                    for idx, _ in GLOBAL_REPORTER.indent_error:
                        curr_idx = idx-1
                        while not LINE_RE.match(lines[curr_idx-1]) and curr_idx >= 0:
                            curr_idx -= 1
                        final_idx = curr_idx - 1
                        if idx - final_idx == 2 and IF_RE.match(lines[final_idx]) and ELSE_RE.match(lines[idx-1]):
                            print(basename, lines[idx-1], lines[final_idx])
                            bad_idx.append(idx-1)
                            bad_idx.append(final_idx)
                        elif IF_RE.match(lines[final_idx]):
                            bad_idx.append(final_idx)
                        else:
                            do_break = True
                            break
                    if do_break:
                        print('ERROR', output_file_path, GLOBAL_REPORTER.indent_error)
                        break
                    for idx in reversed(sorted(bad_idx)):
                        del lines[idx]
                    reset_lists()
                    pyfl.check(''.join(lines), output_file_path, GLOBAL_REPORTER)

                with open(output_file_path, 'w') as write:
                    write.write(''.join(lines))


def main():
    file_dict = get_file_dict()

    put_imports(file_dict)

    remove_unused(file_dict)

if __name__ == '__main__':
    main()
