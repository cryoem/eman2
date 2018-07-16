from __future__ import print_function


with open('py3_divs.txt') as f:
	lines = f.readlines()
	
for l in lines:
	cols = l.split()
	fl = cols[0].split(':')
	file = fl[0]
	num = int(fl[1])
	
	strings = ' '.join(cols[1:]).split('|')
	
	print(file, num, strings)
	
	with open(file, 'r') as f:
		lines = f.readlines()
	
	lines[num - 1] = lines[num - 1].replace(strings[0], strings[1])
	
	with open(file, 'w') as f:
		f.writelines(lines)
