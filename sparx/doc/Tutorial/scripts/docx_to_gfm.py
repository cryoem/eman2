# -*- coding: utf-8 -*-
import os
import re
import argparse


def sanity_check(lines, sanity_list):
    inlines = set(['|||SI|||', '|||URL|||', '|||SIEND|||', '|||URLEND|||'])
    for end_name in sanity_list:
        start = False
        last_occurance = 0
        start_name = end_name.replace(r'END|||', r'|||')
        for idx, line in enumerate(lines):

            if start_name in line and not start:
                last_occurance = idx
                start = True
            elif start_name in line and start:
                start = True
                if start_name in inlines:
                    index = last_occurance
                else:
                    index = last_occurance + 2
                print('line:', last_occurance, ' Missing END for:', start_name)
                print(lines[index])
            else:
                pass

            if end_name in line and start:
                start = False
            elif end_name in line and not start:
                start = False
                if end_name in inlines:
                    index = idx
                else:
                    index = idx - 2
                print('line:', idx, ' Missing START for:', end_name)
                print(lines[index])
            else:
                pass


    figure_re = re.compile(r'.*|||FIGURE[A-Z]*|||(.*?)|||FIGURE[A-Z]*|||.*'.replace('|', '\\|'))
    for idx, line in enumerate(lines):
        if '|||FIGURE' in line:
            figure_names = figure_re.match(line).groups()
            for figure_name in figure_names:
                if '.png' in figure_name:
                    figure_name = 'latex/media/{0}'.format(figure_name.replace(r'\_', '_'))
                    if not os.path.exists(figure_name):
                        print('line:', idx, ' Missing figure:', figure_name)


def main():

    ### PARSE ARGUMENTS

    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='input file in docx format')
    parser.add_argument('output', help='output file in gfm format')
    parser.add_argument('--tags_file', help='tags file')
    args = parser.parse_args()

    ### CONVERT FROM DOCX TO GFM

    command = "pandoc -f docx+empty_paragraphs -t gfm -o {0} --wrap=none {1}".format(
        args.output,
        args.input
        )
    os.system(command)

    ### FIX URL TAGS

    with open(args.output, 'r') as read:
        lines = read.read()
    lines = lines.replace(r'|||URL|||<', r'|||URL|||')
    lines = lines.replace(r'|||URLEND|||>', r'|||URLEND|||')
    lines = lines.replace(r'>|||URLEND|||', r'|||URLEND|||')

    new_lines = []
    match_re = re.compile(r'.*(<https.*?>).*')
    for line in lines.split('\n'):
        match = match_re.match(line)
        if match is not None:
            for entry in match.groups():
                line = line.replace(entry, entry[1:-1])
        new_lines.append(line)

    with open(args.output, 'w') as write:
        write.write('\n'.join(new_lines))

    ### Do sanity check

    if args.tags_file is not None:
        sanity_list = []
        typo_list = []
        with open(args.tags_file, 'r') as read:
            lines = read.readlines()

        for line in lines:
            line = line.strip()
            if not line:
                continue
            else:
                key = [value for value in line.split('\t') if value][0]
                if 'END|||' in key:
                    sanity_list.append(key)
                typo_list.append(key)

        for name in typo_list:
            name_raw = name.replace('|||', '')
            re_check_list = []
            re_check_list.append(re.compile('(?<!\|)(?:\|{{1,2}}|\|{{4,5}}|\|{{7,}}){0}(?:\|{{3}}|\|{{6}})(?!\|)'.format(name_raw)))
            re_check_list.append(re.compile('(?<!\|)(?:\|{{3}}|\|{{6}}){0}(?:\|{{1,2}}|\|{{4,5}}|\|{{7,}})(?!\|)'.format(name_raw)))
            re_check_list.append(re.compile('(?<!\|)(?:\|{{1,2}}|\|{{4,5}}|\|{{7,}}){0}(?:\|{{1,2}}|\|{{4,5}}|\|{{7,}})(?!\|)'.format(name_raw)))
            for idx, line in enumerate(new_lines):
                for entry in re_check_list:
                    match_list = entry.findall(line)
                    for match in match_list:
                        print('line:', idx, ' Tag typo for:', name)
                        print(line)

        sanity_check(new_lines, sanity_list)
    else:
        print('No tag file present for sanity checks')
        pass

if __name__ == '__main__':
    main()
