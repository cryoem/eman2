#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
import collections
import re
import glob
import shutil
import os

used_files = (
    'sx.py',
    'sx_real.py',
    'sx3dvariability.py',
    'sxcompute_isac_avg.py',
    'sxcter.py',
    'sxchains.py',
    'e2bdb.py',
    'e2boxer_old.py',
    'e2boxer.py',
    'e2display.py',
    'e2proc3d.py',
    'sxfilterlocal.py',
    'sxgui_cter.py',
    'sxgui_meridien.py',
    'sxgui_unblur.py',
    'sxheader.py',
    'sxisac2.py',
    'sxlocres.py',
    'sxmeridien.py',
    'sxpdb2em.py',
    'sxpipe.py',
    'sxgui.py',
    'sxcpy.py',
    'sxprocess.py',
    'sxproj_compare.py',
    'sxrelion2sphire.py',
    'sxsphire2relion.py',
    'sxrewindow.py',
    'sxrviper.py',
    'sxsort3d_depth.py',
    'sxsort3d.py',
    'sxsummovie.py',
    'sxunblur.py',
    'sxviper.py',
    'sxwindow.py',
    )

template_files = ['../sphire/templates/sxgui_template.py', '../sphire/templates/wikiparser.py']

bin_dir = '../sphire/bin'
lib_dir = '../sphire/libpy'

tmp_bin_dir = '../sphire/sphire_tmp_bin'
tmp_lib_dir = '../sphire/sphire_tmp_libpy'

unused_bin_dir = '../sphire/sphire_bin_unused'
unused_lib_dir = '../sphire/sphire_libpy_unused'

old_bin_dir = '../sphire/sparx_bin'
new_bin_dir = '../sphire/sphire_bin'
old_lib_dir = '../sphire/sparx_libpy'
new_lib_dir = '../sphire/sphire_libpy'

#tmp_lib_dir = '../sphire/sphire_libpy_tmp'

IMPORT_RE = re.compile('^(?:def|class)\s+([^\(]+)')
END_RE = re.compile('^(?:\w|"|\')')
CONSTANT_RE = re.compile('^(\w+)')
START_RE = re.compile('^(?:def)')
SPARX_FUNC_RE = re.compile('(?:(?<=[^\w])|^)({0})\.([\w]+)'.format('|'.join([
    os.path.splitext(os.path.basename(entry))[0]
    for entry in glob.glob('{0}/*.py'.format(lib_dir))
    if 'sparx.py' not in entry
    ])))
COMMENT_RE = re.compile('^\s*#')
MULTILINE1_RE = re.compile('\'\'\'')
MULTILINE2_RE = re.compile('"""')
FUNCDEF_RE = re.compile('^.*(?:def|class)\s+([^\(]+)')

for dirname in [new_bin_dir, new_lib_dir, unused_bin_dir, unused_lib_dir, old_bin_dir, old_lib_dir, tmp_bin_dir, tmp_lib_dir]:
    try:
        os.mkdir(dirname)
    except OSError as e:
        pass

bin_files = sorted(glob.glob('{0}/*.py'.format(bin_dir)))
lib_files = sorted(glob.glob('{0}/*.py'.format(lib_dir)))

# Remove all comments from the file
#bin_files = sorted(glob.glob('{0}/sxprocess.py'.format(bin_dir)))
#lib_files = []
for idx, out_dir in enumerate([[tmp_bin_dir, old_bin_dir, unused_bin_dir], [tmp_lib_dir, old_lib_dir, unused_lib_dir]]):
    if idx == 0:
        used_list = bin_files
    else:
        used_list = lib_files

    for file_path in used_list:
        basename = os.path.basename(file_path)
        # Do not touch unused files.
        if basename not in used_files and idx == 0:
            shutil.copy2(file_path, out_dir[1])
            continue

        # Import lines
        with open(file_path, 'r') as read:
            lines = read.readlines()
        comment1 = False
        comment2 = False
        usage = False
        wrote_unused = False
        docstring = False
        begin = False
        comment = 0
        with open(os.path.join(out_dir[0], basename), 'w') as write_0, \
                open(os.path.join(out_dir[2], basename), 'w') as write_1:
            for line_idx, line in enumerate(lines):
                ok = False
                match_1 = MULTILINE1_RE.findall(line)
                match_2 = MULTILINE2_RE.findall(line)
                is_comment = COMMENT_RE.match(line)
                is_usage = bool('usage =' in line or 'usage=' in line)
                is_doc = False
                is_first = False
                for i in range(6):
                    if line_idx - i < 0:
                        break
                    elif re.match('^\s*return', lines[line_idx-i]) and (match_1 or match_2):
                        break
                    elif FUNCDEF_RE.match(lines[line_idx-i]) and (match_1 or match_2):
                        is_doc = True
                        break

                if IMPORT_RE.match(line) and not begin and not comment2 and not comment1:
                    begin = True
                elif IMPORT_RE.match(line) and not begin and (comment2 or comment1):
                    print("WOHOOO CHECK THIS OUT", basename, line_idx, line)
                elif not begin:
                    pass

                if is_comment:
                    pass
                elif match_1 and not comment2 and not comment1:
                    if len(match_1) % 2 == 1:
                        is_first = True
                        if usage:
                            usage = False
                        elif is_usage:
                            usage = True
                        elif docstring:
                            docstring = False
                        elif is_doc:
                            docstring = True
                        else:
                            comment1 = True
                elif match_2 and not comment1 and not comment2:
                    if len(match_2) % 2 == 1:
                        is_first = True
                        if usage:
                            usage = False
                        elif is_usage:
                            usage = True
                        elif docstring:
                            docstring = False
                        elif is_doc:
                            docstring = True
                        else:
                            comment2 = True
                elif match_1:
                    if len(match_1) % 2 == 1:
                        comment1 = False
                        ok = True
                elif match_2:
                    if len(match_2) % 2 == 1:
                        comment2 = False
                        ok = True

                if (comment1 or comment2 or ok) and begin:
                    wrote_unused = True
                    if is_first:
                        string = re.match('^(\s*)', line).group(1)
                        write_0.write('{0}"""Multiline Comment{1}"""\n'.format(string, comment))
                        comment += 1
                        write_1.write('{0}{1}\n'.format(line.strip(), comment))
                    else:
                        write_1.write(line)
                else:
                    write_0.write(line)

        if not wrote_unused:
            os.remove(os.path.join(out_dir[2], basename))

bin_files = sorted(glob.glob('{0}/*.py'.format(tmp_bin_dir)))
lib_files = sorted(glob.glob('{0}/*.py'.format(tmp_lib_dir)))

# Extract all the functions from each file
used_functions = {'bin': {}, 'lib': {}}
print('Extract all functions')
for name, files in zip(['bin', 'lib'], [bin_files, lib_files]):
    for file_path in files:
        mod_name = os.path.splitext(os.path.basename(file_path))[0]
        with open(file_path) as read:
            lines = read.readlines()
        start = False
        for idx, line in enumerate(lines):

            match = IMPORT_RE.match(line)
            if match:
                used_functions[name].setdefault(mod_name, {'used': []}).setdefault('functions', []).append(match.group(1))

            if COMMENT_RE.match(line) or '#IMPORTIMPORTIMPORT' in line:
                pass
            elif match and not start:
                func_name = match.group(1)
                start = True
            elif match and start:
                func_name = match.group(1)
                start = True
            elif START_RE.match(line) and start:
                start = False

            if start:
                used_functions[name].setdefault(mod_name, {'used': []}).setdefault(func_name, {}).setdefault('lines', []).append(line)
                used_functions[name].setdefault(mod_name, {'used': []}).setdefault(func_name, {}).setdefault('indices', []).append(idx)

# Search for doubled entries
print('Search for doubles in bin files')
for key, items in used_functions.items():
    for mod_name, values in items.items():
        nr_funcs = [(item, count) for item, count in collections.Counter(values['functions']).items() if count > 1]
        if nr_funcs:
            print(key, mod_name, nr_funcs)

# Do a recursive library search
IGNORED = []
def recursive_check(function_dict, module, key):
    try:
        lines = function_dict[module][key]['lines']
    except KeyError:
        global IGNORED
        if (module, key) not in IGNORED:
            print('IGNORE:', module, key)
            IGNORED.append((module,key))
        return None
    else:
        function_dict[module]['used'].append(key)
        function_dict[module]['used'] = list(set(function_dict[module]['used']))

    current_func_re = re.compile('(?:(?<=[^\w])|^)({0})(?=[^\w]|$)'.format('|'.join([
        entry
        for entry in function_dict[module]
        if entry not in ('used', 'functions')
        ])))
    for line in lines:
        if COMMENT_RE.match(line) or '#IMPORTIMPORTIMPORT' in line:
            continue

        match = SPARX_FUNC_RE.findall(line)
        match_2 = current_func_re.findall(line)
        for mod, name in match:
            if name in function_dict[mod]['used']:
                continue
            else:
                function_dict[mod]['used'].append(name)
                function_dict[mod]['used'] = list(set(function_dict[mod]['used']))
                recursive_check(function_dict, mod, name)

        for name in match_2:
            if name in function_dict[module]['used']:
                continue
            else:
                function_dict[module]['used'].append(name)
                function_dict[module]['used'] = list(set(function_dict[module]['used']))
                recursive_check(function_dict, module, name)

# Extract all the internally used functions from the bin files
print('Find unused functions in bin files')
for name, values in used_functions['bin'].items():
    if name in ('sxgui'):
        print('IGNORE', name)
        continue
    with open(os.path.join(tmp_bin_dir, '{0}.py'.format(name))) as read:
        lines = read.readlines()

    unused_functions = []
    len_unused = -999
    bad_idx = []
    while len(unused_functions) != len_unused:
        len_unused = len(unused_functions)
        tmp_used_functions = []
        tmp_unused_functions = []
        for entry in values['functions']:
            re_comp = re.compile('[^\w^\d^_]{0}[^\w^%]'.format(entry))
            function_re = re.compile('^(?:def|class)\s+{0}[\(]'.format(entry))
            for idx, line in enumerate(lines):
                if idx in set(bad_idx):
                    continue
                if re_comp.search(line) and not function_re.match(line):
                    tmp_used_functions.append(entry)
            tmp_unused_functions = [entry for entry in values['functions'] if entry not in set(tmp_used_functions)]
        for func_name in tmp_unused_functions:
            bad_idx.extend(values[func_name]['indices'])

        unused_functions.extend(tmp_unused_functions)
        unused_functions = list(set(unused_functions))
    print('UNUSED:', name, unused_functions)

    used_lines = []
    unused_lines = []
    for idx, line in enumerate(lines):
        is_used = True
        for entry in unused_functions:
            if idx in values[entry]['indices']:
                is_used = False
                break
        if is_used:
            used_lines.append(line)
        else:
            unused_lines.append(line)

    with open(os.path.join(new_bin_dir, '{0}.py'.format(name)), 'w') as write:
        write.write(''.join(used_lines))

    if unused_lines:
        with open(os.path.join(unused_bin_dir, '{0}.py'.format(name)), 'a') as write:
            write.write(''.join(unused_lines))

    for line in used_lines:
        match = SPARX_FUNC_RE.findall(line)
        for mod, name in match:
            recursive_check(used_functions['lib'], mod, name)
    print('')

for key, vals in used_functions['lib'].items():
    vals['unused'] = [entry for entry in vals['functions'] if entry not in vals['used']]
    print('UNUSED FUNCTIONS', key)
    print(sorted(vals['unused']))

# Removed unused functions
print('Remove unused lib functions')
shutil.copy2(os.path.join(tmp_lib_dir, '{0}.py'.format('sparx')), os.path.join(new_lib_dir, '{0}.py'.format('sparx')))
for name in used_functions['lib']:
    with open(os.path.join(tmp_lib_dir, '{0}.py'.format(name))) as read:
        lines = read.readlines()

    if not used_functions['lib'][name]['used']:
        if name in ('sparx_user_functions', 'user_functions'):
            with open(os.path.join(new_lib_dir, '{0}.py'.format(name)), 'w') as write:
                write.write(''.join(lines))
        else:
            print('UNUSED', name)
        continue

    used_lines = []
    unused_lines = []
    for idx, line in enumerate(lines):
        is_used = True
        for entry in used_functions['lib'][name]['unused']:
            try:
                indices = set(used_functions['lib'][name][entry]['indices'])
            except KeyError:
                indices = []
            if idx in indices:
                is_used = False
                break
        if is_used:
            used_lines.append(line)
        else:
            unused_lines.append(line)

    with open(os.path.join(new_lib_dir, '{0}.py'.format(name)), 'w') as write:
        write.write(''.join(used_lines))

    if unused_lines:
        with open(os.path.join(unused_lib_dir, '{0}.py'.format(name)), 'a') as write:
            write.write(''.join(unused_lines))
