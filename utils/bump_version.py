#!/usr/bin/env python

import argparse
from collections import namedtuple

FileKeyVerPair = namedtuple('FileKeyVerPair', ['fname', 'keyword', 'version_pair'])
OldNew  = namedtuple('OldNew', ['old', 'new'])


def update_versions_in_files(fname, keyword, version_pair):
	print(f"  ...Updating: {fname}, keyword: {keyword}, version pair: {version_pair.old} -> {version_pair.new}")


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
		FileKeyVerPair('ci_support/constructor-mini/construct.yaml', 'eman2',        versions_pair),
		FileKeyVerPair('ci_support/constructor-huge/construct.yaml', 'eman2',        versions_pair),
	)

	for f in eman_versions:
		fname = f.fname
		keyword = f.keyword
		version_pair = f.version_pair

		update_versions_in_files(fname, keyword, version_pair)


if __name__ == "__main__":
	main()
