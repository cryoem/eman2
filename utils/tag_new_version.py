#!/usr/bin/env python

import argparse


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--bump', choices=['major', 'minor', 'patch'], default='patch')

	options = parser.parse_args()
	print(f"Bumping '{options.bump.upper()}' version...")


if __name__ == "__main__":
	main()
