#!/usr/bin/env python
#
# Author: Markus Stabrin 2020 (markus.stabrin@mpi-dortmund.mpg.de)
#
# Copyright (c) 2020 Max Planck Institute of Molecular Physiology
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete SPHIRE and EMAN2 software packages have some GPL dependencies,
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
# ========================================================================================
import sys
import pandas
import re
import argparse



class StarFile:

    def __init__(self, star_file):
        self.star_file = star_file
        self.imported_content = {}
        self.line_dict = {}
        self.analyse_star_file()

    def analyse_star_file(self):
        with open(self.star_file) as read:
            content = read.read()

        # https://regex101.com/r/D7O06N/1
        for tag_match in re.finditer('^data_([^\s]*)\s*$', content, re.M):
            tag = tag_match.group(1)
            self.line_dict[tag] = {
                'block': [None, None],
                'header': [None, None],
                'content': [None, None],
                'is_loop': None,
                }

            current_flag = 0
            prev_content = content[:tag_match.start() + current_flag]
            current_content = content[tag_match.start() + current_flag:]
            current_flag += tag_match.start()

            # https://regex101.com/r/4o3dNy/1/
            self.line_dict[tag]['block'][0] = \
                len(re.findall('\n', prev_content)) + 1

            # https://regex101.com/r/h7Wm8y/2
            header_match = re.search(
                '((?:(?:loop_\s*)?^_.*$\r?\n?)+)',
                current_content,
                re.M
                )
            prev_content = content[:header_match.start() + current_flag]
            current_content = content[header_match.start() + current_flag:]
            current_flag += header_match.start()

            self.line_dict[tag]['is_loop'] = header_match.group(1).startswith('loop_')
            # https://regex101.com/r/4o3dNy/1/
            self.line_dict[tag]['header'][0] = \
                len(re.findall('\n', prev_content)) + 1 + self.line_dict[tag]['is_loop']

            prev_content = content[:header_match.end() + current_flag - header_match.start()]
            current_content = content[header_match.end() + current_flag - header_match.start():]
            current_flag += header_match.end() - header_match.start()
            # https://regex101.com/r/4o3dNy/1/
            self.line_dict[tag]['header'][1] = \
                len(re.findall('\n', prev_content))

            if not self.line_dict[tag]['is_loop']:
                self.line_dict[tag]['content'] = self.line_dict[tag]['header']
            else:
                self.line_dict[tag]['content'][0] = self.line_dict[tag]['header'][1] + 1
                # https://regex101.com/r/HYnKMl/1
                newline_match = re.search('^\s*$', current_content, re.M)

                prev_content = content[:newline_match.start() + current_flag]
                current_content = content[newline_match.start() + current_flag:]
                current_flag += newline_match.start()

                # https://regex101.com/r/4o3dNy/1/
                self.line_dict[tag]['content'][1] = \
                    len(re.findall('\n', prev_content))

            self.line_dict[tag]['block'][1] = self.line_dict[tag]['content'][1]

            self.read_tag(tag, self.line_dict[tag])

    def read_tag(self, tag, line_dict):
        if not line_dict['is_loop']:
            data = self.read_without_loop(line_dict)
        else:
            data = self.read_with_loop(line_dict)
        self.imported_content[tag] = data

    def read_with_loop(self, line_dict):
        header_names = pandas.read_csv(
            self.star_file,
            usecols=[0],
            skiprows=line_dict['header'][0] - 1,
            nrows=line_dict['header'][1] - line_dict['header'][0] + 1,
            skip_blank_lines=False,
            header=None,
            delim_whitespace=True,
            squeeze=True,
            )

        return pandas.read_csv(
            self.star_file,
            index_col=None,
            names=header_names,
            skiprows=line_dict['content'][0]-1,
            nrows=line_dict['content'][1] - line_dict['content'][0] + 1,
            skip_blank_lines=False,
            header=None,
            delim_whitespace=True,
            )

    def read_without_loop(self, line_dict):
        return pandas.read_csv(
            self.star_file,
            index_col=0,
            names=['', '0'],
            skiprows=line_dict['content'][0]-1,
            nrows=line_dict['content'][1] - line_dict['content'][0] + 1,
            skip_blank_lines=False,
            header=None,
            delim_whitespace=True,
            ).transpose()

    def __getitem__(self, tag):
        return self.imported_content[tag]

    def __setitem__(self, tag, data):
        self.imported_content[tag] = data

    def write_star(self, star_file, tags):
        for idx, tag in enumerate(tags):
            if idx == 0:
                mode = 'w'
            else:
                mode = 'a'
            df = self.imported_content[tag]
            is_loop = self.line_dict[tag]['is_loop']

            if is_loop:
                export_header = '\ndata_\n\nloop_\n' + '\n'.join([
                    '{} #{}'.format(entry, idx)
                    for idx, entry
                    in enumerate(df, 1)
                    ])

                with open(star_file, mode) as write:
                    write.write(f'{export_header}\n')
                df.to_csv(star_file, sep='\t', header=False, index=False, mode='a')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('output')
    return parser.parse_args()


def main():
    args = parse_args()
    star_file = StarFile(args.input)
    tag_key = 'particles'
    try:
        star_file[tag_key]
    except KeyError:
        tag_key = 'micrographs'
        try:
            star_file[tag_key]
        except KeyError:
            print('Could not read particles or micrographs entry!')
            return 2

    try:
        optics_keys = star_file['optics'].groupby('_rlnOpticsGroup')
    except KeyError:
        return 1

    try:
        particle_keys = star_file[tag_key].groupby('_rlnOpticsGroup')
    except KeyError:
        return 1

    for val_particles, df_particles in particle_keys:
        for val_optics, df_optics in optics_keys:
            if val_particles == val_optics:
                for key in df_optics:
                    if key == '_rlnImagePixelSize':
                        new_key = '_rlnDetectorPixelSize'
                        star_file[tag_key].loc[df_particles.index, '_rlnMagnification'] = 10000
                    elif key == '_rlnHelicalTrackLengthAngst':
                        new_key = '_rlnHelicalTrackLength'
                    else:
                        new_key = key

                    star_file[tag_key].loc[df_particles.index, new_key] = star_file['optics'].loc[df_optics.index, key].iloc[0]

                try:
                    star_file[tag_key].loc[df_particles.index, '_rlnOriginX'] = \
                        star_file[tag_key].loc[df_particles.index, '_rlnOriginXAngst'] / \
                        star_file['optics'].loc[df_optics.index, '_rlnImagePixelSize'].iloc[0]
                except KeyError:
                    pass
                try:
                    star_file[tag_key].loc[df_particles.index, '_rlnOriginY'] = \
                        star_file[tag_key].loc[df_particles.index, '_rlnOriginYAngst'] / \
                        star_file['optics'].loc[df_optics.index, '_rlnImagePixelSize'].iloc[0]
                except KeyError:
                    pass

    star_file.write_star(args.output, [tag_key])
    return 0

if __name__ == '__main__':
    sys.exit(main())
