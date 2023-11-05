#!/usr/bin/env python

import argparse
import os

import requests
import subprocess as sp


bump_choices = ('major', 'minor', 'patch')
bump_choices_idx = {v:i for i,v in enumerate(bump_choices)}

def get_latest_version_from_gh():
	tags = [t['name'] for t in requests.get('https://api.github.com/repos/cryoem/eman2/tags').json()]
	print(f"Received GitHub tags:\n{tags}")

	# tags = ['v2.99', 'v2.99.33', 'release', 'v2.44a-windoiws'] ->
	# tags = [['2', '99'], ['2', '99', '33'], ['2', '44a-windoiws']]
	tags = [t[1:].split('.') for t in tags if t.startswith('v')]

	# tags = [['2', '99'], ['2', '99', '33'], ['2', '44a-windoiws']] ->
	# tags = [[2, 99], [2, 99, 33]]
	tags = [[int(word) for word in t] for t in tags if all([word.isnumeric() for word in t])]

	# tags = [[2, 99], [2, 99, 33]] ->
	# tags = [[2, 99, 0], [2, 99, 33]]
	for t in tags:
		t.extend((len(bump_choices)-len(t))*[0])  # pad with 0s

	# tags = [[2, 99, 0], [2, 99, 33]] ->
	# tags = [[2, 99, 33]], [2, 99, 0]
	tags = sorted(tags, reverse=True)
	print(f"Version tags (sorted: latest to oldest):\n{tags}")

	version = tags[0]

	print(f"Current version:  {version}")

	return version

def bump_version_tuple(bump_bit, version):
	version = version.copy()
	i_bump = bump_choices_idx[bump_bit]

	# Bump requested version bit
	version[i_bump] += 1

	# Reset lower version bits
	for i in range(i_bump + 1, len(version)):
		version[i] = 0

	print(f"Bumped version bits: {version}")

	return version

def push_tag(version_cur, version_new):
	branch = os.getenv('GITHUB_REF_NAME')

	version_cur = '.'.join([str(i) for i in version_cur])
	version_new = '.'.join([str(i) for i in version_new])

	for cmd in (
			f'python utils/bump_version.py {version_cur} {version_new}',
			f'git commit -a -m v{version_new}',
			f'git tag -f v{version_new}',
			f'git push origin {branch} v{version_new}',
			):
		print(f"> {cmd}")
		output = sp.run(cmd.split(), capture_output=True, text=True).stdout
		print(output)

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--bump', choices=bump_choices, default='patch')

	options = parser.parse_args()

	bump_bit = options.bump.lower()

	print(f"Bumping '{bump_bit}' version...")

	version_cur = get_latest_version_from_gh()
	version_new = bump_version_tuple(bump_bit, version_cur)
	push_tag(version_cur, version_new)

	with open(os.getenv('GITHUB_ENV'), "a") as fout:
		fout.write(f"new_version={'.'.join([str(i) for i in version_new])}")


if __name__ == "__main__":
	main()
