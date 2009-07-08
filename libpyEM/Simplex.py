#!/usr/bin/env python
# 
# Copyright (c) 2001 Vivake Gupta (v@omniscia.org).  All rights reserved.
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
#
# This software is maintained by Vivake (v@omniscia.org) and is available at:
#     http://www.omniscia.org/~vivake/python/Simplex.py

# Modified (debugged?) 7/16/2004 Michele Vallisneri (vallis@vallis.org)

""" Simplex - a regression method for arbitrary nonlinear function minimization

Simplex minimizes an arbitrary nonlinear function of N variables by the
Nelder-Mead Simplex method as described in:

Nelder, J.A. and Mead, R., "A Simplex Method for Function Minimization",
   Computer Journal, Vol. 7, 1965, pp. 308-313.

It makes no assumptions about the smoothness of the function being minimized.
It converges to a local minimum which may or may not be the global minimum
depending on the initial guess used as a starting point.
"""

import math
import copy
import sys

class Simplex:
    def __init__(self, testfunc, guess, increments, kR = -1, kE = 2, kC = 0.5, data=None):
        """Initializes the simplex.
        INPUTS
        ------
        testfunc      the function to minimize
        guess[]       an list containing initial guesses
        increments[]  an list containing increments, perturbation size
        kR            reflection constant  (alpha =-1.0)
        kE            expansion constant   (gamma = 2.0)
        kC            contraction constant (beta  = 0.5)
        """
        self.testfunc = testfunc
        self.guess = guess
        self.increments = increments
        self.data=data
        self.kR = kR
        self.kE = kE
        self.kC = kC
        self.numvars = len(self.guess)
        self.simplex = []

        self.lowest = -1
        self.highest = -1
        self.secondhighest = -1

        self.errors = []
        self.currenterror = 0

        # Initialize vertices
        # MV: the first vertex is just the initial guess
        #     the other N vertices are the initial guess plus the individual increments
        #     the last two vertices will store the centroid and the reflected point
        #     the compute errors at the ... vertices
        
        for vertex in range(0, self.numvars + 3):
            self.simplex.append(copy.copy(self.guess))

        for vertex in range(0, self.numvars + 1):
            for x in range(0, self.numvars):
                if x == (vertex - 1):
                    self.simplex[vertex][x] = self.guess[x] + self.increments[x]
            self.errors.append(0)
        self.calculate_errors_at_vertices()

    def minimize(self, epsilon = 0.0001, maxiters = 250, monitor = 1):
        """Walks to the simplex down to a local minima.
        INPUTS
        ------
        epsilon       convergence requirement
        maxiters      maximum number of iterations
        monitor       if non-zero, progress info is output to stdout  

        OUTPUTS
        -------
        an array containing the final values
        lowest value of the error function
        number of iterations taken to get here
        """
        
        iter = 0
        
        for iter in range(0, maxiters):
            # Identify highest and lowest vertices
            
            self.highest = 0
            self.lowest = 0
            for vertex in range(1, self.numvars + 1):
                if self.errors[vertex] > self.errors[self.highest]:
                    self.highest = vertex
                if self.errors[vertex] < self.errors[self.lowest]:
                    self.lowest = vertex

            # Identify second-highest vertex

            self.secondhighest = self.lowest
            for vertex in range(0, self.numvars + 1):
                if vertex == self.highest:
                    continue
                elif vertex == self.secondhighest:
                    continue
                elif self.errors[vertex] > self.errors[self.secondhighest]:
                    self.secondhighest = vertex

            # Test for convergence:
            #   compute the average merit figure (ugly)

            S = 0.0
            for vertex in range(0, self.numvars + 1):
                S = S + self.errors[vertex]
            F2 = S / (self.numvars + 1)

            #   compute the std deviation of the merit figures (ugly)

            S1 = 0.0
            for vertex in range(0, self.numvars + 1):
                S1 = S1 + (self.errors[vertex] - F2)**2
            T = math.sqrt(S1 / self.numvars)
            
            # Optionally, print progress information

            if monitor:
                print '\r' + 72 * ' ',
                print '\rIteration = %d   Best = %f   Worst = %f' % \
                      (iter,self.errors[self.lowest],self.errors[self.highest]),
                sys.stdout.flush()
                
            if T <= epsilon:
                # We converged!  Break out of loop!
                
                break;
            else:
                # Didn't converge.  Keep crunching.
                
                # Calculate centroid of simplex, excluding highest vertex
                # store centroid in element N+1

                # loop over coordinates
                for x in range(0, self.numvars):
                    S = 0.0
                    for vertex in range(0, self.numvars + 1):
                        if vertex == self.highest:
                            continue
                        S = S + self.simplex[vertex][x]
                    self.simplex[self.numvars + 1][x] = S / self.numvars

                # reflect the simplex across the centroid
                # store reflected point in elem. N + 2 (and self.guess)
                
                self.reflect_simplex()
                self.currenterror = self.testfunc(self.guess,self.data)

                if self.currenterror < self.errors[self.highest]:
                    self.accept_reflected_point()

                if self.currenterror <= self.errors[self.lowest]:
                    self.expand_simplex()
                    self.currenterror = self.testfunc(self.guess,self.data)

                    # at this point we can assume that the highest
                    # value has already been replaced once
                    if self.currenterror < self.errors[self.highest]:
                        self.accept_expanded_point()
                elif self.currenterror >= self.errors[self.secondhighest]:
                    # worse than the second-highest, so look for
                    # intermediate lower point

                    self.contract_simplex()
                    self.currenterror = self.testfunc(self.guess,self.data)

                    if self.currenterror < self.errors[self.highest]:
                        self.accept_contracted_point()
                    else:
                        self.multiple_contract_simplex()
                        
        # Either converged or reached the maximum number of iterations.
        # Return the lowest vertex and the currenterror.

        for x in range(0, self.numvars):
            self.guess[x] = self.simplex[self.lowest][x]
        self.currenterror = self.errors[self.lowest]
        return self.guess, self.currenterror, iter

    # same as expand, but with alpha < 1; kC = 0.5 fine with NR

    def contract_simplex(self):
        for x in range(0, self.numvars):
            self.guess[x] = self.kC * self.simplex[self.highest][x] + (1 - self.kC) * self.simplex[self.numvars + 1][x]
        return

    # expand: if P is vertex and Q is centroid, alpha-expansion is Q + alpha*(P-Q),
    #         or (1 - alpha)*Q + alpha*P; default alpha is 2.0; agrees with NR
    def expand_simplex(self):
        for x in range(0, self.numvars):
            self.guess[x] = self.kE * self.guess[x]                 + (1 - self.kE) * self.simplex[self.numvars + 1][x]
        return

    # reflect: if P is vertex and Q is centroid, reflection is Q + (Q-P) = 2Q - P,
    #          which is achieved for kR = -1 (default value); agrees with NR
    def reflect_simplex(self):
        # loop over variables
        for x in range(0, self.numvars):
            self.guess[x] = self.kR * self.simplex[self.highest][x] + (1 - self.kR) * self.simplex[self.numvars + 1][x]
            # store reflected point in elem. N + 2
            self.simplex[self.numvars + 2][x] = self.guess[x]
        return

    # multiple contraction: around the lowest point; agrees with NR

    def multiple_contract_simplex(self):
        for vertex in range(0, self.numvars + 1):
            if vertex == self.lowest:
                continue
            for x in range(0, self.numvars):
                self.simplex[vertex][x] = 0.5 * (self.simplex[vertex][x] + self.simplex[self.lowest][x])
        self.calculate_errors_at_vertices()
        return

    def accept_contracted_point(self):
        self.errors[self.highest] = self.currenterror
        for x in range(0, self.numvars):
            self.simplex[self.highest][x] = self.guess[x]
        return

    def accept_expanded_point(self):
        self.errors[self.highest] = self.currenterror
        for x in range(0, self.numvars):
            self.simplex[self.highest][x] = self.guess[x]
        return

    def accept_reflected_point(self):
        self.errors[self.highest] = self.currenterror
        for x in range(0, self.numvars):
            self.simplex[self.highest][x] = self.simplex[self.numvars + 2][x]
        return

    def calculate_errors_at_vertices(self):
        for vertex in range(0, self.numvars + 1):
            if vertex == self.lowest:
                # compute the error unless we're the lowest vertex
                continue
            for x in range(0, self.numvars):
                self.guess[x] = self.simplex[vertex][x]
            self.currenterror = self.testfunc(self.guess,self.data)
            self.errors[vertex] = self.currenterror
        return

def myfunc(args):
    return abs(args[0] * args[0] * args[0] * 5 - args[1] * args[1] * 7 + math.sqrt(abs(args[0])) - 118)

def main():
    s = Simplex(myfunc, [1, 1, 1], [2, 4, 6])
    values, err, iter = s.minimize()
    print 'args = ', values
    print 'error = ', err
    print 'iterations = ', iter

if __name__ == '__main__':
    main()
