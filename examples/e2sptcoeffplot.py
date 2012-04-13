#!/usr/bin/env python

from EMAN2 import *
from sys import argv
from operator import itemgetter	
import matplotlib.pyplot as plt
import sys

data=argv[1]
cutoff=None
if len(argv) > 2:
	cutoff=argv[2]

scores=[]
dataset=[]
x=[]

n=EMUtil.get_image_count(data)
n2=round(n/2.0)

for i in range(n):
	a=EMData(data,i)
	hdr=a.get_attr_dict()
	
	if 'spt_score' in hdr:
		score=-1*a['spt_score']
	elif 'spt_coefficient' in hdr:
		score=a['spt_coefficient']
	else:
		print "No score found in header. Terminating!"
		sys.exit()
	scores.append(score)
	dataset.append({'score':float(score), 'particle':a})
	x.append(i)

dataset = sorted(dataset, key=itemgetter('score'), reverse=True)

halfptclnum=int(round(n/2.0))
halfptcl=dataset[keyptclnum]
print "keyptcl is", keyptcl
halfscore=keyptcl['score']

for s in scores:
	if s < halfscore:
		scores.remove(s)
		print "I have removed this aberrantly low score", s

scores.sort()
scores.reverse()

thresh=cutoff

if cutoff=='half':
	thresh="%0.2f" % (halfscore)

thresh=str(thresh).replace('.','p')

group1=[]
g1name = data.replace('.hdf','_th' + thresh + 'G1.hdf')
group2=[]
g2name = data.replace('.hdf','_th' + thresh + 'G2.hdf')

if cutoff == 'half':
	
	k1=0
	k2=0
	for i in dataset:
		if i['score'] >= keyscore:
			group1.append(i['particle'])
			i['particle'].write_image(g1name,k1)
			k1+=1
			print "I have appended a particle to group 1"
		else:
			if ['score'] > halfscore/2
				group2.append(i['particle'])
				i['particle'].write_image(g2name,k2)
				k2+=1	
				print "I have appended a particle to group 2"
			else:
				print "\n\n@@@@ Found a garbage particle!\n\n"


if cutoff and cutoff != 'half':
	cutoff=float(cutoff)
	k1=0
	k2=0
	for i in dataset:
		if i['score'] >= cutoff:
			group1.append(i['particle'])
			i['particle'].write_image(g1name,k1)
			k1+=1
		else:
			if ['score'] > halfscore/2
				group2.append(i['particle'])
				i['particle'].write_image(g2name,k2)
				k2+=1
			else:
				print "\n\n@@@@ Found a garbage particle!\n\n"	
	
if group1:
	g1avg = sum(group1)/len(group1)
	g1avg.write_image(g1name.replace('.hdf','_avg.hdf'),0)

if group2:
	g2avg = sum(group2)/len(group2)
	g2avg.write_image(g2name.replace('.hdf','_avg.hdf'),0)
			
plot_name = data.replace('.hdf','_scores.png')
plt.plot(x, scores, linewidth=1, marker='o', linestyle='--', color='r')
plt.title(plot_name)
plt.ylabel('score')
plt.xlabel('n ptcl')
a = plt.gca()

#a.set_xlim(1,mults[-1])
#a.set_ylim(0,max(x))

plt.savefig(plot_name)
plt.clf()
