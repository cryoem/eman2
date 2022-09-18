#!/usr/bin/env python

import argparse


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('current_version')
	parser.add_argument('new_version')

	options = parser.parse_args()
	print(f"Updating version: {options.current_version} -> {options.new_version} ...")


if __name__ == "__main__":
	main()
