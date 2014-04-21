#!/usr/bin/env python

#
# Author: Vernon Williams (vewillia@bcm.edu)
# Copyright (c) 2014- Baylor College of Medicine
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
#
#-------------------------------------------------------------------------
#
# Functions:
#
# have_match = matches (string, pattern)
# have_match = matches_pats (string, patterns)
#
# matches returns a boolean value telling whether the string in the
# first argument matches a wildcard pattern string in the second argument.
#
# matches_pat returns a boolean value telling whether the string in the
# first argument matches any of the blank separated wildcard patterns
# in the string in the second argument.
#
# A wildcard pattern string may contain any of the 5 following wildcard
# characters:
#
#   * - match zero or more of any characters
#   ? - match zero or one of any characters
#   % - match one of any character
#   # - match one decimal digit 0 to 9
#   ! - match one upper or lower case letter (a to z or A to Z)
#
# The patterns argument to matches_pats may contain a special value
# 'not', meaning that a match will occur only if the string matches
# all of the wildcard patterns before the 'not' but not any after the
# 'not'.
#
# Examples of pattern arguments:
#
# 1.  "%??.txt" - match any file names ending with ".txt" and
#     having a 1, 2, or 3 character base name
#
# 2.  "!*.jpg not *#.jpg" - match any file names ending in ".jpg"
#     and beginning with a letter but the base name not ending in
#     a decimal digit
#
# Written  by Vernon Williams, May      26, 1994.
# Modified by Vernon Williams, November  3, 2009,
#    to make new version with "w_" prefixed to the name
#    to match my current file naming conventions,
#    to add documentation, to add local variable matches,
#    to return false if string or pattern is a null pointer,
#    to add include file string.h.
# Modified by Vernon Williams, April    18, 2014,
#    to convert to Python.
#
#--------------------------------------------------------------------------

# For internal use:

def pat_match (string, pattern, str_len, pat_len) :
	# Constants:

	any_string = '*'  # Match any chars
	any_char   = '?'  # Match 0 or 1 chars
	one_char   = '%'  # Match 1 char
	alp_char   = '!'  # Match 1 letter
	num_char   = '#'  # Match 1 digit

	done = False

	if pat_len < str_len :
		nchars = pat_len
	else :
		nchars = str_len

	for ic in xrange (0, nchars) :
		str_char = string  [ic]
		pat_char = pattern [ic]

		if pat_char == any_string  or  pat_char == any_char :
			if pat_char == any_char :
				ic_last = ic + 1
			else :
				ic_last = str_len

			ic_pat = ic + 1

			match  = False

			for ic_str in xrange (ic, ic_last + 1) :
				match = pat_match (string [ic_str : ], pattern [ic_pat : ],
							str_len - ic_str, pat_len - ic_pat)
				if match :
					break

			done = True
		else :
			if pat_char == one_char :
				match = True
			elif pat_char == alp_char :
				match = str_char.isalpha()
			elif pat_char == num_char :
				match = str_char.isdigit()
			else :
				match = (pat_char == str_char)

			if not match :
				done = True

		if done :
			break

	if not done :
		if max (str_len, 0) == max (pat_len, 0) :
			match = True
		elif str_len > pat_len :
			match = False
		else : # str_len < pat_len
			match = True

			for ic in xrange (max (str_len, 0), pat_len) :
				pat_char = pattern [ic]

				if pat_char != any_string  and  pat_char != any_char :
					match = False
					break
	return match

#--------------------------------------------------------------------------

# For export:

def matches (string, pattern) :
	if isinstance (string, str)  and  isinstance (pattern, str) :
		str_len = len (string)
		pat_len = len (pattern)
		matches = pat_match (string, pattern, str_len, pat_len)
	else :
		matches = false

	return matches

#--------------------------------------------------------------------------

# For export:

def matches_pats (string, patterns) :
#	constants:

	exclude_flag = "not"

	if isinstance (string, str)  and  isinstance (patterns, str) :
		pats = patterns.split ( )

		not_in  = False

		if len (pats) > 0  and  pats [0] == exclude_flag :
			matched = True
		else :
			matched = False

		for pat in pats :
			if pat == exclude_flag  and  not not_in :
				not_in = True
			else :
				match = matches (string, pat)

				if not_in :
					if match :
						matched = False
						break
				else :
					if match :
						matched = True
	else :
		matched = False

	return matched
