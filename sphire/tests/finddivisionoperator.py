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
GREATEREQUAL_RE = re.compile('(>=)')
BRACKETS_RE = re.compile('\(')
STARSIGN_RE = re.compile('\*')
DOUBLESLASH_RE = re.compile('//')
ROUNDINT_RE = re.compile('int')
SELFBRACKETS_RE = re.compile('self.brackets')
EQUALSIGN_RE = re.compile('=')

with open('sphire/libpy/sparx_fundamentals.py') as read:
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
        is_greaterqual = GREATEREQUAL_RE.findall(line)
        is_bracket = BRACKETS_RE.findall(line)
        is_starsign = STARSIGN_RE.findall(line)
        is_doubleslash = DOUBLESLASH_RE.findall(line)
        is_roundint = ROUNDINT_RE.findall(line)
        is_selfbrackets = SELFBRACKETS_RE.findall(line)
        is_equalsign = EQUALSIGN_RE.findall(line)


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
                if len(line) < 1:   # if length of lines is greater than 50 thats means I have to handle those lines manually
                    continue
                if is_divequalsign:
                    string = re.match('^(\s*)', line).group(1)
                    print('old syntex', line)
                    print('\n')
                    print('new syntex', string + line.split()[0] + ' = ' + "old_div" + '(' + line.split()[0] + ', ' + line.split()[2] + ')')
                    print('\n')
                    lines[idx] = string + line.split()[0] + ' = ' + "old_div" + '(' + line.split()[0] + ', ' + line.rstrip().split()[2] + ')' + '\n'


                #  #### For one / operator #######
                elif line.count('/') == 1 and not is_bracket and not is_starsign and not is_roundint \
                        and not line.strip().split()[0] == 'wedgeFactor' and not line.strip().split()[0] == 'theta':
                    string = re.match('^(\s*)', line).group(1)
                    print('old syntex', line)
                    # print(line.partition('#')[0],'\n')
                    print('new syntex', string + line.split()[0] + ' ' + line.split()[1] + ' ' + "old_div" + '(' + line.split()[2] \
                          + ', ' + line.partition('#')[0].rstrip().split('/')[1] + ')')
                    print('\n')
                    lines[idx] = string + line.split()[0] + ' ' + line.split()[1] + ' ' + "old_div" + '(' + line.split()[2] \
                          + ', ' + line.partition('#')[0].rstrip().split('/')[1] + ')' + '\n'

                elif line.count('/') == 1 and not is_bracket and not is_starsign and not is_roundint \
                         and not line.strip().split()[0] == 'wedgeFactor' and not line.strip().split()[0] == 'theta' and line.split()[4] == 'self.nsym':
                    string = re.match('^(\s*)', line).group(1)
                    print('old syntex', line)
                    print('new syntex',
                          string + line.split()[0] + ' ' + line.split()[1] + ' ' + "old_div" + '(' + line.split()[2] \
                          + ', ' +  line.split()[4] +  ')' + line.partition('#')[0].rstrip().split('self.nsym')[1] )
                    print('\n')
                    lines[idx] = string + line.split()[0] + ' ' + line.split()[1] + ' ' + "old_div" + '(' + line.split()[2] \
                          + ', ' +  line.split()[4] +  ')' + line.partition('#')[0].rstrip().split('self.nsym')[1] + '\n'

                elif line.count('/') == 1 and not is_bracket and is_starsign and not is_roundint \
                        and not is_selfbrackets and not line.strip().split('.append')[0] == 'self.symangles' \
                        and not line.strip().split()[0] == 'wedgeFactor' and not line.strip().split()[0] == 'theta':
                    string = re.match('^(\s*)', line).group(1)
                    print('old syntex', line)
                    print('\n')
                    print('new syntex', line.split('*')[0] + '* ' + "old_div" + '(' + line.split('*')[1].split('/')[0] \
                          + ', ' + line.partition('#')[0].rstrip().split('/')[1]+ ')')
                    print('\n')
                    lines[idx] = line.split('*')[0] + '* ' + "old_div" + '(' + line.split('*')[1].split('/')[0] + ', ' + \
                                 line.partition('#')[0].rstrip().split('/')[1] + ')' + '\n'

                elif line.count('/') == 1 and is_bracket and is_starsign and not is_roundint \
                        and not is_selfbrackets and not line.strip().split('.append')[0] == 'self.symangles' \
                        and not line.strip().split()[0] == 'wedgeFactor' and not line.strip().split()[0] == 'theta' \
                        and is_equalsign:
                    string = re.match('^(\s*)', line).group(1)
                    print('old syntex', line)
                    print('\n')
                    print(line.rstrip().split())
                    print(line.partition('='))

                    print(line.split('=')[1].split('/')[0])
                    print(line.partition('=')[2].partition('/')[2].split('*')[0])


                    if len(line.rstrip().split()) > 5 :
                        print('new syntex', string + line.split()[0] + ' = ' + "old_div" + '(' + line.split('=')[1].split('/')[0] \
                              + ',' + line.split()[6] + ')' + line.split()[7] + line.split()[8] )
                        print('\n')
                        lines[idx] = string + line.split()[0] + ' = ' + "old_div" + '(' + line.split('=')[1].split('/')[0] \
                              + ',' + line.split()[6] + ')' + line.split()[7] + line.split()[8] + '\n'


                    if len(line.rstrip().split()) <= 5 :
                        print('new syntex', string + line.partition('/')[0]  + "old_div" + '(' + line.split('=')[1].split('/')[0] \
                              + ',' + line.partition('=')[2].partition('/')[2].split('*')[0]  + + ')'  + line.partition('=')[2].partition('*')[2]     )
                        print('\n')
                        # lines[idx] = string + line.split()[0] + ' = ' + "old_div" + '(' + line.split('=')[1].split('/')[0] \
                        #       + ',' + line.split()[6] + ')' + line.split()[7] + line.split()[8] + '\n'


                #### For two / operator #######


                elif line.count('/') == 2 and not is_bracket and not is_starsign and not is_doubleslash \
                        and not is_roundint and not is_selfbrackets  and not line.strip().split('.append')[0] == 'self.symangles' \
                        and not line.rstrip().split()[-1] == 'self.nsym' and is_equalsign:
                    string = re.match('^(\s*)', line).group(1)
                    print('old syntex', line)
                    print('\n')
                    # print(line.split())
                    print('new syntex',
                          string + line.split()[0] + ' ' + line.split()[1] + ' ' + "old_div" + '(' + \
                          "old_div" + '(' + line.split()[2] + ', ' + line.split('/')[1] + ')' \
                          + ',' + line.strip().split('/')[-1] + ')')
                    print('\n')
                    lines[idx] = string + line.split()[0] + ' ' + line.split()[1] + ' ' + "old_div" + '(' + \
                          "old_div" + '(' + line.split()[2] + ', ' + line.split('/')[1] + ')' \
                          + ',' + line.strip().split('/')[-1] + ')' + '\n'

                elif line.count('/') == 2    and not is_bracket and not is_starsign and not is_doubleslash \
                        and not is_roundint and not is_selfbrackets  and not line.strip().split('.append')[0] == 'self.symangles' \
                        and line.rstrip().split()[-1] == 'self.nsym' and is_equalsign:
                    string = re.match('^(\s*)', line).group(1)
                    print('old syntex', line)
                    print('\n')
                    print('new syntex',
                          string + line.split()[0] + ' ' + line.split()[1] + ' ' + \
                          "old_div" + '(' + line.split()[2] + ', ' + line.split()[4] + ')' \
                          + ' ' + line.split()[5]  + ' ' + line.split()[6] + ' ' + line.split()[7] + ' ' + "old_div" + '(' + line.split()[8] + ', ' + line.split()[10] + ')' + ')')
                    print('\n')
                    lines[idx] = string + line.split()[0] + ' ' + line.split()[1] + ' ' + \
                          "old_div" + '(' + line.split()[2] + ', ' + line.split()[4] + ')' \
                          + ' ' + line.split()[5]  + ' ' + line.split()[6] + ' ' + line.split()[7] + \
                                 ' ' + "old_div" + '(' + line.split()[8] + ', ' + line.split()[10] + ')' + ')' + '\n'

                elif line.count('/') == 2 and not is_bracket and not is_starsign and is_doubleslash and not is_roundint and is_equalsign:
                    string = re.match('^(\s*)', line).group(1)
                    print('old syntex', line)
                    print('\n')
                    print('new syntex',
                          string + line.split()[0] + ' = ' +
                          "old_div" + '(' + line.split('=')[1].split('//')[0] + ', ' + line.strip().split('=')[1].split('//')[1] + ')')
                    print('\n')
                    lines[idx] = string + line.split()[0] + ' = ' + \
                          "old_div" + '(' + line.split('=')[1].split('//')[0] + ', ' + line.strip().split('=')[1].split('//')[1] + ')' + '\n'


                elif line.count('/') == 2 and not is_bracket and is_starsign and not is_doubleslash and not is_roundint and is_equalsign:
                    string = re.match('^(\s*)', line).group(1)
                    print('old syntex', line)
                    print('\n')
                    print('new syntex',
                          string + line.split()[0] + ' = ' +
                          "old_div" + '(' + line.split('=')[1].split('//')[0] + ', ' +
                          line.strip().split('=')[1].split('//')[1] + ')')
                    print('\n')
                    lines[idx] = string + line.split()[0] + ' = ' + \
                          "old_div" + '(' + line.split('=')[1].split('//')[0] + ', ' + \
                          line.strip().split('=')[1].split('//')[1] + ')' + '\n'


                elif line.count('/') == 2 and not is_bracket and is_starsign and is_doubleslash and not is_roundint and is_equalsign:
                    string = re.match('^(\s*)', line).group(1)
                    print('old syntex', line)
                    print('\n')
                    print('new syntex',
                          string + line.split()[0] + ' = ' +
                          "old_div" + '(' + line.split('=')[1].split('//')[0] + ', ' +
                          line.strip().split('=')[1].split('//')[1] + ')')
                    print('\n')
                    lines[idx] = string + line.split()[0] + ' = ' + \
                          "old_div" + '(' + line.split('=')[1].split('//')[0] + ', ' + \
                          line.strip().split('=')[1].split('//')[1] + ')' + '\n'

                elif line.count('/') == 2 and is_bracket and not is_starsign and not is_doubleslash and not is_roundint \
                        and not line.strip().split()[0] == 'elif' and not line.strip().split()[0] == 'if' and is_equalsign:
                    string = re.match('^(\s*)', line).group(1)
                    print('old syntex', line)
                    print('\n')
                    print('new syntex',
                          string + line.split()[0] + ' = ' +
                          "old_div" + '(' + line.split('=')[1].split('//')[0] + ', ' +
                          line.strip().split('=')[1].split('//')[1] + ')')
                    print('\n')
                    lines[idx] = string + line.split()[0] + ' = ' + \
                          "old_div" + '(' + line.split('=')[1].split('//')[0] + ', ' + \
                          line.strip().split('=')[1].split('//')[1] + ')' + '\n'

                elif line.count('/') == 2 and is_bracket and not is_starsign and is_doubleslash and not is_roundint \
                      and not line.strip().split()[0] == 'elif' and not line.strip().split()[0] == 'if' and is_equalsign:
                    string = re.match('^(\s*)', line).group(1)
                    print('old syntex', line)
                    print('\n')
                    print('new syntex',
                          string + line.split()[0] + ' = ' +
                          "old_div" + '(' + line.split('=')[1].split('//')[0] + ', ' +
                          line.strip().split('=')[1].split('//')[1] + ')')
                    print('\n')
                    lines[idx] = string + line.split()[0] + ' = ' + \
                          "old_div" + '(' + line.split('=')[1].split('//')[0] + ', ' + \
                          line.strip().split('=')[1].split('//')[1] + ')' + '\n'

                elif line.count('/') == 2 and is_bracket and is_starsign and not is_doubleslash and not is_roundint \
                        and not line.strip().split()[0] == 'elif' and not line.strip().split()[0] == 'if' and is_equalsign:
                    string = re.match('^(\s*)', line).group(1)
                    print('old syntex', line)
                    print('\n')
                    print('new syntex',
                          string + line.split()[0] + ' = ' +
                          "old_div" + '(' + line.split('=')[1].split('//')[0] + ', ' +
                          line.strip().split('=')[1].split('//')[1] + ')')
                    print('\n')
                    lines[idx] = string + line.split()[0] + ' = ' + \
                          "old_div" + '(' + line.split('=')[1].split('//')[0] + ', ' + \
                          line.strip().split('=')[1].split('//')[1] + ')' + '\n'

                elif line.count('/') == 2 and is_bracket and is_starsign and is_doubleslash and not is_roundint and is_equalsign:
                    string = re.match('^(\s*)', line).group(1)
                    print('old syntex', line)
                    print('\n')
                    print('new syntex',
                          string + line.split()[0] + ' = ' +
                          "old_div" + '(' + line.split('=')[1].split('//')[0] + ', ' +
                          line.strip().split('=')[1].split('//')[1] + ')')
                    print('\n')
                    lines[idx] = string + line.split()[0] + ' = ' + \
                          "old_div" + '(' + line.split('=')[1].split('//')[0] + ', ' + \
                          line.strip().split('=')[1].split('//')[1] + ')' + '\n'

                elif line.strip().split()[0] == 'if' and is_doubleslash and is_bracket:
                    string = re.match('^(\s*)', line).group(1)
                    # print('with if statement on line', idx + 1)
                    print('old_syntex', line )
                    print('\n')
                    print('new_syntex',  string + line.split()[0]  + ' (' + "old_div" + '(' + line.split()[1].split('(')[-1] + ',' + \
                          line.split('//')[1].split(')')[0]  + ')' + line.rstrip().split(')')[1] + '):')
                    print('\n')
                    lines[idx] = string + line.split()[0]  + ' (' + "old_div" + '(' + line.split()[1].split('(')[-1] + ',' + \
                          line.split('//')[1].split(')')[0]  + ')' + line.rstrip().split(')')[1] + '):' + '\n'

                elif line.count('/') == 1 and is_bracket and is_greaterqual and line.strip().split()[0] == 'if' :
                    string = re.match('^(\s*)', line).group(1)
                    print('old syntex', line)
                    print('\n')
                    print('new syntex', string + line.split()[0]  + ' ' + line.strip().split()[1] + ' ' + line.strip().split()[2]  + \
                          ' (' + "old_div" + '(' + line.strip().split()[3] +',' + line.rstrip().partition('if')[2].split('/')[1].replace(':',')):') )
                    print('\n')
                    lines[idx] = string + line.split()[0]  + ' ' + line.strip().split()[1] + ' ' + line.strip().split()[2]  + \
                          ' (' + "old_div" + '(' + line.strip().split()[3] +',' + line.rstrip().partition('if')[2].split('/')[1].replace(':',')):') + '\n'

                elif line.count('/') > 2  and is_bracket and is_greaterqual and line.strip().split()[0] == 'if':
                    print('\n')
                    print('old_syntex', line)
                    print('new syntex', line.rstrip().split('/')[0].partition('360.0')[0]  + ' (' + "old_div" + '('  + "old_div" + '(' + \
                          line.rstrip().split('/')[0].partition('360.0')[1] + ',' +  line.strip().split('/')[1] + ')' + ',' + \
                          line.strip().split('/')[2].split()[0] + ')' + line.strip().split('/')[2].partition('360')[0].rstrip().replace('2', '') + \
                          ' (' + "old_div" + '(' + line.strip().split('/')[2].split()[-1] + ',' + line.strip().split('/')[3].partition(')')[0]  + ')' + '):' )
                    print('\n')
                    lines[idx] = line.rstrip().split('/')[0].partition('360.0')[0]  + ' (' + "old_div" + '('  + "old_div" + '(' + \
                          line.rstrip().split('/')[0].partition('360.0')[1] + ',' +  line.strip().split('/')[1] + ')' + ',' + \
                          line.strip().split('/')[2].split()[0] + ')' + line.strip().split('/')[2].partition('360')[0].rstrip().replace('2', '') + \
                          ' (' + "old_div" + '(' + line.strip().split('/')[2].split()[-1] + ',' + line.strip().split('/')[3].partition(')')[0]  + ')' + '):' + '\n'


                elif line.strip().split()[0] == 'elif' and is_doubleslash and is_bracket:
                    print('old_syntex', line )
                    print('\n')
                    print('new_syntex', line.partition('and')[0] + line.partition('and')[1]  ,
                          ' (' + "old_div" + '(' + line.rstrip().split()[5].split('(')[1] + ',' +  line.split()[7] + ')' + line.partition('2)')[2])
                    print('\n')
                    lines[idx] = line.partition('and')[0] + line.partition('and')[1]  , \
                          ' (' + "old_div" + '(' + line.rstrip().split()[5].split('(')[1] + ',' +  line.split()[7] + ')' + line.partition('2)')[2] + '\n'

                elif line.count('/') == 1 and is_bracket and is_greaterqual and line.strip().split()[0] == 'elif':
                    # print('with elif statement on line', idx + 1)
                    print('old syntex', line)
                    print('\n')
                    print('new syntex', line.split()[0] + ' ' + line.strip().split()[1] + ' ' + line.strip().split()[2] + \
                          ' (' + "old_div" + '(' + line.strip().split()[3] + ',' + \
                          line.rstrip().partition('if')[2].split('/')[1].replace(':', ')):'))
                    print('\n')
                    lines[idx] =  line.split()[0] + ' ' + line.strip().split()[1] + ' ' + line.strip().split()[2] + \
                          ' (' + "old_div" + '(' + line.strip().split()[3] + ',' + \
                          line.rstrip().partition('if')[2].split('/')[1].replace(':', ')):') + '\n'


                elif is_selfbrackets:
                    print('ignore selfbrackets on line', idx + 1)


                elif line.strip().split('.append')[0] == 'self.symangles':
                    print('ignore self.symangles on line', idx + 1)


                elif line.strip().split()[0] == 'wedgeFactor':
                    print('ignore wedgeFactor command on line', idx+1)

                elif line.strip().split()[0] == 'theta':
                    print('ignore theta command on line', idx+1)

                elif not is_equalsign:
                    print ('equal sign is missing on line', idx+1)


# with open('sphire/tests/testdivisionchangenew.py','w') as newfile:
#     for idx, line in enumerate(lines[:]):
#         # print(lines[idx])
#         newfile.write(str(line))