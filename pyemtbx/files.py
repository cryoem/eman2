import os
import re

# remove files under 'absdir' using pattern
def remove_files(pattern, absdir=None):
	pattern = pattern.replace(".", "\\.")
	pattern = pattern.replace("*", ".*")
	if not absdir:
		absdir = os.getcwd()
		
	files = os.listdir(absdir)

	for file in files:
		if re.match(pattern, file):
			os.remove(file)
			
