#!/usr/bin/env python

import glob
import os
import shutil
import sys
from pathlib import Path
from subprocess import run
from multiprocessing import Pool


THIS_DIR = Path(__file__).resolve().parent


def get_program_files():
	source_path = THIS_DIR.parent.parent / 'programs'

	progs_exclude = {
		"e2.py",
		"e2projectmanager.py",
		"e2unwrap3d.py",
		"e2version.py",
		"e2fhstat.py",
		"e2_real.py",
		"e2proc3d.py",  # uses OptionParser
		"e2seq2pdb.py",  # no help provided
	}

	progs = set(Path(p).name for p in source_path.glob('e2*.py')) - progs_exclude

	return sorted(progs)


def write_parser_options_of_prog(prog):
	output = run([prog, '--help-to-html'], capture_output=True)

	if output.stderr:
		print(f'Error while trying to run {prog} ...')
		print(output.stderr)
		sys.exit(1)

	prog = Path(prog).stem
	txt = output.stdout.decode()
	ftxt = THIS_DIR / str(prog + '.html')

	with open(ftxt, 'w') as f:
		print("Writing {} ...".format(f.name))
		f.write(txt)


def main():
	progs = get_program_files()

	with Pool() as pool:
		pool.map(write_parser_options_of_prog, progs)

	outdir = 'programs_help_html'
	os.makedirs(outdir, exist_ok=True)

	shutil.copy2('styles.css', outdir)
	for f in glob.glob('*.html'):
		fnew = os.path.join(outdir, f)
		print(f"Moving: {f} -> {fnew}")
		shutil.move(f, fnew)


if __name__ == "__main__":
	main()
