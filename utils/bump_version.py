#!/usr/bin/env python

import argparse
from collections import namedtuple

FileKeyVerPair = namedtuple('FileKeyVerPair', ['fname', 'keyword', 'version_pair'])
OldNew  = namedtuple('OldNew', ['old', 'new'])


def update_versions_in_files(fname, keyword, version_pair):
	print(f"  ...Updating: {fname}, keyword: {keyword}, version pair: {version_pair.old} -> {version_pair.new}")

	with open(fname, 'r') as fin:
		lines = fin.readlines()

	for i, k in enumerate(lines):
		if keyword in k and version_pair.old in k:
			print(f"  ......Found    line: {k}", end='')
			lines[i] = k.replace(version_pair.old, version_pair.new)
			print(f"  ......Replaced line: {lines[i]}")

	with open(fname, 'w') as fout:
		fout.writelines(lines)


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('current_version')
	parser.add_argument('new_version')

	options = parser.parse_args()
	print(f"Updating version: {options.current_version} -> {options.new_version} ...")

	versions_pair = OldNew(options.current_version, options.new_version)

	eman_versions = (
		FileKeyVerPair('CMakeLists.txt',                             'VERSION',      versions_pair),
		FileKeyVerPair('ci_support/set_env_vars.sh',                 'EMAN_VERSION', versions_pair),
		FileKeyVerPair('ci_support/construct.yaml', 'eman2',        versions_pair),
	)

	for f in eman_versions:
		update_versions_in_files(f.fname, f.keyword, f.version_pair)


if __name__ == "__main__":
	main()
