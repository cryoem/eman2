#!/usr/bin/env python

from pathlib import Path


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


def main():
	progs = get_program_files()


if __name__ == "__main__":
	main()
