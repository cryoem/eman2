# -*- coding: utf-8 -*-
from __future__ import print_function
from builtins import range
import os
import re
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='input in gfm format')
    parser.add_argument('output', help='output in latex format')
    parser.add_argument('--acronyms_file', help='acronyms file used by lates')
    parser.add_argument('--tags_file', help='tags file')
    return parser.parse_args()


def load_tag_dict(tags_file):
    dictionary = {}
    if tags_file is not None:
        with open(tags_file, 'r') as read:
            for line in read.readlines():
                line = line.strip()
                if not line:
                    continue
                else:
                    key, value = [value for value in line.split('\t') if value]
                    dictionary[key.strip()] = value.strip()
    else:
        pass

    return dictionary


def load_lines(input_file):
    with open(input_file, 'r') as read:
        lines = read.read()
    return lines


def replace_tags(lines, replace_dict):
    sanity_list = []
    typo_list = []
    for key, value in list(replace_dict.items()):
        if 'END|||' in key:
            sanity_list.append(key)
        typo_list.append(key)
    sanity_check(lines, sanity_list)
    for name in typo_list:
        name_raw = name.replace('|||', '')
        re_check_list = []
        re_check_list.append(re.compile('(?<!\|)(?:\|{{1,2}}|\|{{4,5}}|\|{{7,}}){0}(?:\|{{3}}|\|{{6}})(?!\|)'.format(name_raw)))
        re_check_list.append(re.compile('(?<!\|)(?:\|{{3}}|\|{{6}}){0}(?:\|{{1,2}}|\|{{4,5}}|\|{{7,}})(?!\|)'.format(name_raw)))
        re_check_list.append(re.compile('(?<!\|)(?:\|{{1,2}}|\|{{4,5}}|\|{{7,}}){0}(?:\|{{1,2}}|\|{{4,5}}|\|{{7,}})(?!\|)'.format(name_raw)))
        for idx, line in enumerate(lines.split('\n')):
            for entry in re_check_list:
                match_list = entry.findall(line)
                for match in match_list:
                    print('line:', idx, ' Tag typo for:', name)
                    print(line)

    for key, value in list(replace_dict.items()):
        if not key in lines:
            print(key, 'not found')
        lines = lines.replace(key, value)

    return lines


def sanity_check(lines, sanity_list):
    for end_name in sanity_list:
        start = False
        last_occurance = 0
        start_name = end_name.replace(r'END|||', r'|||')
        for idx, line in enumerate(lines.split('\n')):
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
    for idx, line in enumerate(lines.split('\n')):
        if '|||FIGURE' in line:
            figure_names = figure_re.match(line).groups()
            for figure_name in figure_names:
                if '.png' in figure_name:
                    figure_name = 'latex/media/{0}'.format(figure_name.replace(r'\_', '_'))
                    if not os.path.exists(figure_name):
                        print('line:', idx, ' Missing figure:', figure_name)


def fix_whitespaces(lines):
    new_lines = []
    for line in lines.split('\n'):
        if line:
            new_lines.append(line)
        else:
            pass
    return '\n'.join(new_lines).replace('|||NEWLINE|||', '\n')


def replaces(lines):
    in_quote = "“"
    out_quote = "”"
    count = 0
    new_lines = []
    for char in lines:
        if char == '"':
            if count % 2 == 0:
                new_lines.append(r'\enquote{')
            else:
                new_lines.append(r'}')
            count += 1
        else:
            new_lines.append(char)
    lines = ''.join(new_lines)
    lines = lines.replace(in_quote, r'\enquote{')
    lines = lines.replace(out_quote, r'}')
    lines = lines.replace('\\\\', '\\')
    lines = lines.replace(r'%', r'\%')
    lines = lines.replace(r'$', r'\$')
    lines = lines.replace(r'\$\Rightarrow\$', r'$\Rightarrow$')
    lines = lines.replace(r'\$\cdot\$', r'$\cdot$')
    lines = lines.replace(r'\$\alpha\$', r'$\alpha$')
    lines = lines.replace(r'\$\beta\$', r'$\beta$')
    lines = lines.replace(r'\$\sigma\$', r'$\sigma$')
    lines = lines.replace(r'\[', r'[')
    lines = lines.replace(r'\-', r'-')
    lines = lines.replace(r'\>', r'>')
    lines = lines.replace(r'>>>', r'\textgreater\hbox{}\textgreater\hbox{}\textgreater')
    lines = lines.replace(r'\<\<\<', r'\textless\hbox{}\textless\hbox{}\textless')
    lines = lines.replace(r'=>', r'=\textgreater')
    lines = lines.replace(r'<=', r'\textless=')
    lines = lines.replace(r'\]', r']')
    lines = lines.replace(r'\]', r']')
    lines = lines.replace(r'&', r'\&')
    lines = lines.replace(r'\*', r'*')
    lines = lines.replace(r'A^2', r'A\textasciicircum2')
    for i in reversed(list(range(2, 30))):
        dashes = ['-' for i in range(i)]
        lines = lines.replace(r''.join(dashes), r'\texttt{{{0}}}'.format('{}'.join(dashes)))
    lines = lines.replace(r'\labelitemi{\texttt{-{}-}}', r'\labelitemi{--}')
    lines = lines.replace(r'\caption{', r'\caption*{')
    lines = lines.replace(r'\begin{equation}', r'\begin{equation*}')
    lines = lines.replace(r'\end{equation}', r'\end{equation*}')
    return lines


def regex_italic(lines):
    new_lines = []
    template_out = '\\textit{{{0}}}'
    template_in = '*{0}*'
    italic_re = re.compile(r".*\*([^(\*)]*)\*")
    for line in lines.split('\n'):
        while True:
            match = italic_re.match(line)
            if match is None:
                break
            else:
                phrase = match.group(1)
                if phrase.endswith(r'\enquote{'):
                    break
                else:
                    line = line.replace(template_in.format(phrase), template_out.format(phrase))
                pass
        new_lines.append(line)
    return '\n'.join(new_lines)


def regex_bolt(lines):
    new_lines = []
    template_out = '\\textbf{{{0}}}'
    template_in = '**{0}**'
    bolt_re = re.compile(r".*\*\*(.*)\*\*")
    for line in lines.split('\n'):
        while True:
            match = bolt_re.match(line)
            if match is None:
                break
            else:
                phrase = match.group(1)
                line = line.replace(template_in.format(phrase), template_out.format(phrase))
                pass
        new_lines.append(line)
    return '\n'.join(new_lines)


def regex_si(lines):
    new_lines = []
    si_re = re.compile(r"\\SI{([^}]*)}")
    for line in lines.split('\n'):
        match = si_re.findall(line)
        for phrase in match:
            phrase_split = phrase.split(' ')
            line = line.replace(phrase, '{0}}}{{\\{1}'.format(phrase_split[0].replace('−', '-'), '\\'.join(phrase_split[1:])))
            line = line.replace(r'}{\1', r'}{1')
        new_lines.append(line)
    return '\n'.join(new_lines)


def regex_num(lines):
    new_lines = []
    num_re = re.compile("[ =][0-9.]*[0-9]+[^a-zA-Z}|]+")
    match_re = re.compile('[0-9.]*[0-9]+')
    template = '\\num{{{0}}}'
    for line in lines.split('\n'):
        match = num_re.findall(line)
        for phrase in match:
            num_match = match_re.findall(phrase)
            new_phrase = phrase.replace(num_match[0], template.format(num_match[0]))
            line = line.replace(phrase, new_phrase)
        new_lines.append(line)
    return '\n'.join(new_lines)


def regex_sirange(lines):
    new_lines = []
    sirange_re = re.compile(r"\\SIrange{([^}]*)}")
    for line in lines.split('\n'):
        match = sirange_re.findall(line)
        for phrase in match:
            phrase_split = phrase.split(' ')
            number_split = phrase_split[0].split('|')
            line = line.replace(phrase, '{0}}}{{{1}}}{{\\{2}'.format(number_split[0].replace('−', '-'), number_split[1].replace('−', '-'), '\\'.join(phrase_split[1:])))
        new_lines.append(line)
    return '\n'.join(new_lines)


def regex_media(lines):
    new_lines = []
    media_re = re.compile(r"media\/([^}]*)}")
    for line in lines.split('\n'):
        match = media_re.findall(line)
        for phrase in match:
            new_phrase = phrase.replace('\\', '')
            line = line.replace(phrase, new_phrase)
        new_lines.append(line)
    return '\n'.join(new_lines)


def regex_equation(lines):
    equation_in = []
    equation_out = []
    lines_split = lines.split('\n')
    for idx, phrase in enumerate(lines_split):
        if phrase == '\\begin{equation*}':
            equation_in.append(idx)
        elif phrase == '\\end{equation*}':
            equation_out.append(idx)
        else:
            pass
    assert len(equation_in) == len(equation_out)
    for idx_in, idx_out in zip(equation_in, equation_out):
        new_line = []
        for line in lines_split[idx_in+1:idx_out]:
            if line.strip():
                new_line.append(line)
            else:
                pass
        lines_in = '\n'.join(lines_split[idx_in+1:idx_out])
        lines_out = '{0}'.format('\n'.join(new_line))
        lines = lines.replace(lines_in, lines_out)
    return lines


def regex_numrange(lines):
    new_lines = []
    numrange_re = re.compile(r"\\numrange{([^}]*)}")
    for line in lines.split('\n'):
        match = numrange_re.findall(line)
        for phrase in match:
            phrase_split = phrase.split(' ')
            line = line.replace(phrase, '{0}}}{{{1}'.format(phrase_split[0], phrase_split[1]))
        new_lines.append(line)
    return '\n'.join(new_lines)


def regex_acro(lines, acronyms_file):
    if acronyms_file is not None:
        with open(acronyms_file, 'r') as read:
            acro_lines = read.readlines()
        acro_re = re.compile(r"\\acro{([^}]*)}")
        template_out = '\\ac{{{0}}}'
        template_in = '|||{0}|||'
        for line in acro_lines:
            if line.strip().startswith(r'\acro{'):
                match = acro_re.findall(line)
                for phrase in match:
                    lines = lines.replace(template_in.format(phrase), template_out.format(phrase))
    else:
        pass
    return lines


def main():
    args = parse_args()
    replace_dict = load_tag_dict(tags_file=args.tags_file)
    lines = load_lines(args.input)
    #lines = regex_num(lines)
    lines = replace_tags(lines, replace_dict)
    lines = regex_sirange(lines)
    lines = regex_bolt(lines)
    lines = regex_italic(lines)
    lines = regex_si(lines)
    lines = regex_acro(lines, args.acronyms_file)
    lines = regex_numrange(lines)
    lines = regex_media(lines)
    lines = regex_equation(lines)
    lines = replaces(lines)
    lines = fix_whitespaces(lines)

    if '|||' in lines:
        print('Pattern ||| found in output file: Most probably failed')
    else:
        print('Most probably worked')

    with open(args.output, 'w') as wrt:
        wrt.write(lines)

main()

