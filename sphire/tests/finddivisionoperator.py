#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
import re
from past.utils import old_div

import numpy


# IMPORT_RE = re.compile('^(?:def|class)\s+([^\(]+)')
# END_RE = re.compile('^(?:\w|"|\')')
# CONSTANT_RE = re.compile('^(\w+)')
# START_RE = re.compile('^(?:def)')
COMMENT_RE = re.compile('^\s*#')
MULTILINE1_RE = re.compile('\'\'\'')
MULTILINE2_RE = re.compile('"""')
FUNCDEF_RE = re.compile('^.*(?:def|class)\s+([^\(]+)')



DIVLINE_RE =  re.compile('(/)(?:(?=.*)|$)')
DIVEQUALSIGN_RE = re.compile('(/=)(?:(?=.*)|$)')

with open('sphire/tests/testdivisionchange.py') as read:
    lines = read.readlines()
    comment1 = False
    comment2 = False
    usage = False
    wrote_unused = False
    docstring = False
    begin = False
    comment = 0

    for idx, line in enumerate(lines[:]):

        ok = False
        match_1 = MULTILINE1_RE.findall(line)
        match_2 = MULTILINE2_RE.findall(line)
        is_division = DIVLINE_RE.findall(line)
        is_comment = COMMENT_RE.match(line)
        is_divequalsign = DIVEQUALSIGN_RE.findall(line)
        is_usage = bool('usage =' in line or 'usage=' in line)
        is_doc = False
        is_first = False

        if is_comment:
            continue
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
            else:
                ok = True
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
            else:
                ok = True
        elif match_1:
            if len(match_1) % 2 == 1:
                comment1 = False
                ok = True
            else:
                ok = True
        elif match_2:
            if len(match_2) % 2 == 1:
                comment2 = False
                ok = True
            else:
                ok = True
        if not comment1 and not comment2 and not ok:
            for i in is_division:
                if(len(line) > 50 ):   # if length of lines is greater than 50 thats means I have to handle those lines manually
                    continue
                if is_divequalsign:
                    # print(line)
                    string = re.match('^(\s*)', line).group(1)
                    print('old syntex', line)
                    print('\n')
                    print('new syntex', string + line.split()[0] + ' = ' + "old_div" + '(' + line.split()[0] + ', ' + line.split()[2] + ')')
                    print('\n')
                    lines[idx] = string + line.split()[0] + ' = ' + "old_div" + '(' + line.split()[0] + ', ' + line.split()[2] + ')' + '\n'

                elif line.count('/') == 1:
                    # print(line )
                    # print(line.split())
                    string = re.match('^(\s*)', line).group(1)
                    print('old syntex', line)
                    print('\n')
                    print('new syntex', string + line.split()[0] + ' ' + line.split()[1] + ' ' + "old_div" + '(' + line.split()[2] + ', ' + line.split()[4] + ')')
                    print('\n')



                # else:
                #     # continue
                #     print(i, idx+1, len(line),line)
                #     print(line.count('/'))
                    # for m in re.finditer('/', line):
                    #     print(' / found', m.start(), m.end())
                # print(idx, len(line), line)



with open('sphire/tests/testdivisionchangenew.py','w') as newfile:
    for idx, line in enumerate(lines[:]):
        # print(lines[idx])
        newfile.write(str(line))
    # RE1 = re.compile('(/)(?:(?=.*)|$)')
    # RE2 = re.compile('^\s*#')
    # RE3 = re.compile('"""')


    # match1 = RE1.findall(line)
    # match2 = RE2.match(line)
    # match3 = RE3.match(line)
    # if (match2 or match3):
    #     continue
    # for i in match1:
    #     print(i, idx+1, len(line),line)