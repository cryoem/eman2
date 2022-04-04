#!/usr/bin/env python

import glob
import os
import shutil
import sys
from pathlib import Path
from subprocess import run
from multiprocessing import Pool


THIS_DIR = Path(__file__).resolve().parent


def parser_options_html(prog, txt, has_detailed_info_page):
	detailed_info_txt = f'''
<p>
For more information go to <a href="https://blake.bcm.edu/emanwiki/EMAN2/Programs/{prog}">emanwiki/EMAN2/Programs/{prog}</a>.
</p>
'''

	if not has_detailed_info_page:
		detailed_info_txt = ''

	return f'''<!DOCTYPE html>
<html lang="en-us">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <link rel="stylesheet" href="styles.css">
  <title> {prog} - help</title>
</head>

<body >
<h1> {prog} </h1>

{txt}

{detailed_info_txt}

</body>
</html>
'''


def html_link(prog_html_path, prog):
	return f'''    <li><a href="{prog_html_path}">{prog}</a></li>\n'''

with open('programs_with_more_info.txt') as fin:
	PROGS_WITH_DETAILED_INFO = fin.read().splitlines()


def get_program_files():
	source_path = THIS_DIR.parent.parent / 'programs'

	progs_exclude = {
		"e2.py",
		"e2projectmanager.py",
		"e2unwrap3d.py",
		"e2version.py",
		"e2fhstat.py",
		"e2_real.py",
		"e2proc3d.py",              # uses OptionParser
		"e2seq2pdb.py",             # no help provided
		"e2spt_autoboxer.py",       # broken
		"e2spt_classaverage.py",    # broken
		"e2spt_resolutionplot.py",  # broken
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
		f.write(parser_options_html(prog, txt,
		                            prog + '.py' in PROGS_WITH_DETAILED_INFO))


def main():
	progs = get_program_files()

	with open(THIS_DIR / 'index.html.tmpl', 'r') as fin:
		template_lines = fin.read().partition('================')

	with open(THIS_DIR / 'index.html', 'w') as fout:
		fout.write(template_lines[0])

		for prog in progs:
			fout.write(html_link(prog_html_path=Path(prog).stem + '.html', prog=prog))

		fout.write(template_lines[2])

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
