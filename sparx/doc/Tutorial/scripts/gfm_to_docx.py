# -*- coding: utf-8 -*-
import os
import argparse


def main():

    ### PARSE ARGUMENTS

    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='input file in gfm format')
    parser.add_argument('output', help='output file in docx format')
    args = parser.parse_args()

    ### CONVERT FROM DOCX TO GFM

    command = "pandoc -f gfm -t docx -o {0} --wrap=none {1}".format(
        args.output,
        args.input
        )
    os.system(command)


if __name__ == '__main__':
    main()
