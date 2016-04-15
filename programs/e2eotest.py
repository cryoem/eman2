#!/usr/bin/env python


# Author: David Woolford, 12/9/2008 (woolford@bcm.edu)
# Copyright (c) 2000-2007 Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#




def main():
	print """Resolution testing in EMAN2.1 is performed using 'gold standard' procedures integrated into 
e2refine-easy. The old-style eotest result could be subject to noise bias and other problems and is
no longer considered a valid measure of resolution. This statement also extends to EMAN1. That is,
eotest results from EMAN1 may be rejected by reviewers as not following modern practices. Note that
this is a good thing in many different respects. The new method:

1. no longer requires a separate time-consuming program to run
2. actually helps the refinement converge faster
3. gives a true resolution curve for every iteration
4. provides automatic optimal filtering of the final output map (if CTF correction is performed)

After completing a run of e2refine_easy, please look at the refine_xx/report/index.html file, which
will show a variety of resolution and convergence plots for your convenience. The raw FSC curves are
also available as refine_xx/fsc*, and can be opened by double-clicking on them in the browser"""

if __name__ == "__main__":
    main()
