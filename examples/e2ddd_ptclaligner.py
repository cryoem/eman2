#!/usr/bin/env python

#
# Author: James Michael Bell 08/06/2015
# Copyright (c) 2015 Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these
# packages if you opt to use BSD licensing. The warranty disclaimer below
# holds in either instance.
#
# This complete copyright notice must be included in any revised version
# of the source code. Additional authorship citations may be added, but
# existing author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA	2111-1307 USA
#

from EMAN2 import *

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """e2ddd_ptclaligner.py [options]

	Determines the optimal per-particle alignment of boxed particles from DDD movie frames.
	
	Example: e2ddd_ptclaligner.py --dddmovie 20140926_43522_raw.hdf --average 20140926_43522_raw_avg.hdf --project .
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--dddmovie",type=str,default=None,help="The DDD movie whose particles you wish to align",required=True)
	parser.add_argument("--average",type=str,default=None,help="The averaged DDD movie frames you used for boxing",required=True)
	parser.add_argument("--project",type=str,default=None,help="Location of eman2 'project' containing e2boxercache and info directories", required=True)
	parser.add_argument("--maxiters",type=int,default=5,help="How many times you wish to iteratively improve the alignment")
	
	(options, args) = parser.parse_args()
	
	pid=E2init(sys.argv)
	
	nfs = EMUtil.get_image_count(options.dddmovie)
	coords = js_open_dict(options.project + '/info/' + options.average[:-4] + '_info.json')
	boxes = [b[0:2] for b in coords['boxes']]
	base = js_open_dict(options.project + '/e2boxercache/base.json')
	boxsize = base['box_size']
	ptcls = options.project+'/particles/'+options.average[:-4] + '_ptcls.hdf'
	prepost = []
	for b,box in enumerate(boxes[0:5]):
		print('box {}/{}'.format(b+1,len(boxes)))
		ptcl = EMData(ptcls,b)
		ptcl.process_inplace('normalize.edgemean')
		# iteratively align particle frames to the avg of all of the particle's frames
		for iter in xrange(options.maxiters):
			if iter == 0: bavg=Averagers.get('mean')
			aavg=Averagers.get('mean')
			for i in xrange(nfs):
				r = Region(box[0]-boxsize/2,box[1]-boxsize/2,boxsize,boxsize)
				d = EMData(options.dddmovie,i,False,r)
				d.process_inplace('normalize.edgemean')
				if iter == 0: bavg.add_image(d)
				d2 = d.align('translational',ptcl,{'intonly':0, 'masked':0, 'maxshift':5, 'nozero':0, 'useflcf':1})
				t = d2.get_attr('xform.align2d')
				d2.transform(t)
				aavg.add_image(d2)
			if iter == 0: before = bavg.finish()
			ptcl = aavg.finish()
		
		prepost.append(before)
		prepost.append(ptcl)
	
	display(prepost)
	
	# optionally(?) use interpolation to shift pixels in regions not containing particles (weighted by particles whos shifts have been calculated)
	
	E2end(pid)

if __name__ == "__main__":
	main()
