#!/usr/bin/env python
#
# Author: James Michael Bell, 07/16/2015 (jmbell@bcm.edu)
# Copyright (c) 2015 Baylor College of Medicine
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

from EMAN2 import *
import time
import datetime
import sys
import os
import signal
import copy
import pickle
import abc
import numpy as np

def main():

    progname = os.path.basename(sys.argv[0])
    usage = """Anneal.py [options]

    Optimize a function via simulated annealing.
    """

    parser = EMArgumentParser(usage=usage,version=EMANVERSION)
    parser.add_argument("--tstep",type=float,default=0.9,help="Amount to decrease T by each temperature step")
    parser.add_argument("--nt",type=int,default=10000,help="Number of T steps to cycle through")
    parser.add_argument("--nsteps",type=int,default=1000,help="Number of points to sample at each T step")
    parser.add_argument("--nlimit",type=int,default=100000,help="Maximum successful changes before continuing")
    parser.add_argument("--guess",type=float,default=(-1.0,1.0),nargs="+",help="Initial guess for each dimension. Default: -1.0 1.0. NOTE: Separate these two parameters with spaces.")
    parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
    parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

    (options, args) = parser.parse_args()

    if len(options.guess) > 2:
        print("You should only specify a lower and upper bound with the --guess parameter")
        sys.exit(1)
    elif len(options.guess) == 1:
        high = np.abs(options.guess[0])
        low = -1*high
    else:
        low,high = options.guess

    x0 = [0 for i in range(10)]
    fx0 = testfunc(x0)
    print("ANALYTIC:")
    print("argmin:\t{}\nminval:\t{}\n".format(x0,fx0))

    print("ANNEALED:")
    bounds = [(low,high) for i in x0]
    sa = SimpleAnnealer(testfunc,bounds)
    params = {'t_step':options.tstep,'nt':options.nt,'n_steps':options.nsteps,'n_limit':options.nlimit}
    x1,fx1,err1,iters1=sa.anneal(**params)
    print("argmin:\t{}\nminval:\t{}\n".format(x1,fx1))
    error1 = np.abs(fx1 - fx0)
    print("Iterations:\t\t{}\nEstimated Error:\t{}\nActual Error:\t\t{}\n\n".format(iters1,err1*0.01,error1))

    print("Traveling Sales Person Problem:")
    example_tsp_problem()

def testfunc(x,*kwargs): return np.sum([xi*xi-np.cos(10*xi) for xi in x])

def geodesic_distance(a, b):
    """Calculates distance between two latitude-longitude coordinates."""
    R = 3963  # radius of Earth (miles)
    lat1, lon1 = np.radians(a[0]), np.radians(a[1])
    lat2, lon2 = np.radians(b[0]), np.radians(b[1])
    return np.arccos(np.sin(lat1) * np.sin(lat2) + np.cos(lat1) * np.cos(lat2) * np.cos(lon1 - lon2)) * R

def example_tsp_problem():
    # latitude and longitude for the twenty largest U.S. cities
    cities = {'New York City': (40.72, 74.00),'Los Angeles': (34.05, 118.25),'Chicago': (41.88, 87.63),'Houston': (29.77, 95.38),'Phoenix': (33.45, 112.07),'Philadelphia': (39.95, 75.17),'San Antonio': (29.53, 98.47),'Dallas': (32.78, 96.80),'San Diego': (32.78, 117.15),'San Jose': (37.30, 121.87),'Detroit': (42.33, 83.05),'San Francisco': (37.78, 122.42), 'Jacksonville': (30.32, 81.70),'Indianapolis': (39.78, 86.15),'Austin': (30.27, 97.77),'Columbus': (39.98, 82.98),'Fort Worth': (32.75, 97.33),'Charlotte': (35.23, 80.85),'Memphis': (35.12, 89.97),'Baltimore': (39.28, 76.62)}
    # initial state, a randomly-ordered itinerary
    init_state = list(cities.keys())
    np.random.shuffle(init_state)
    # create a distance matrix
    distance_matrix = {}
    for ka, va in cities.items():
        distance_matrix[ka] = {}
        for kb, vb in cities.items():
            if kb == ka: distance_matrix[ka][kb] = 0.0
            else: distance_matrix[ka][kb] = geodesic_distance(va, vb)
    tsps = TSPSolver(init_state, distance_matrix)
    # since our state is just a list, slice is the fastest way to copy
    tsps.copy_strategy = "slice"
    state, e = tsps.anneal()
    while state[0] != 'New York City':
        state = state[1:] + state[:1]  # rotate NYC to start
    print("%i mile route:" % e)
    for city in state:
        print("\t"+city)

class SimpleAnnealer:

    """
    Multi-dimensional Simulated Annealing Optimization

    objective :    A function of (x1,x2,...) , bounded by bounds
    bounds    :    List of tuples [(min1,max1),(min2,max2), ...]

    Modified from Anneal.py - Simple Simulated Annealing in Python

    Copyright (c) 2003 Richard P. Muller (rmuller@sandia.gov). All rights
    reserved. See the LICENSE file for licensing details.
    """

    def __init__(self,objective,bounds):
        self.objective = objective
        self.bounds = bounds

    def variator(self,x):
        xnew = []
        for i in range(len(x)):
            lower,upper = self.bounds[i]
            var = 0.1*(upper-lower)
            xnew.append(x[i] + var * (2*np.random.random()-1))
        return xnew

    def anneal(self,t_step=0.9,nt=10000,n_steps=1000,n_limit=10000):
        x = self.make_initial_guess(self.bounds) # random point between xmin and xmax
        T = self.initialize(x)
        iters = 0
        for i in range(nt):
            n_succ = 0
            for j in range(n_steps):
                xnew = self.variator(x)
                dE = self.objective(xnew)-self.objective(x)
                if dE < 0 or np.random.random() < np.exp(-dE/T):
                    x = xnew
                    n_succ += 1
                iters += 1
                if n_succ >= n_limit:
                    break
            if n_succ == 0:
                break # no successes -> converged
            T *= t_step
        return (x,self.objective(x),dE,iters)

    @staticmethod
    def make_initial_guess(bounds):
        """Get random values between the bounds for the functions"""
        return [lower+np.random.random()*(upper-lower) for (lower,upper) in bounds]

    def initialize(self,x):
        """
        Run 10 random variations from the starting guess to obtain
        an initial temperature
        """
        E0 = self.objective(x)
        Es = []
        for i in range(10):
            xnew = self.variator(x)
            Es.append(self.objective(xnew))
        T = 2*max(map(lambda E,E0=E0: abs(E-E0),Es))
        return T

class BaseAnnealer(object):

    # This software is maintained by perrygeo and is available at:
    #     https://github.com/perrygeo/simanneal/blob/master/simanneal/anneal.py
    # Last Modified: kim0 on Apr 9, 2015

    """
    Performs simulated annealing by calling functions to calculate
    energy and make moves on a state.  The temperature schedule for
    annealing may be provided manually or estimated automatically.
    """

    __metaclass__ = abc.ABCMeta
    Tmax = 25000.0
    Tmin = 2.5
    steps = 50000
    updates = 100
    copy_strategy = 'slice'
    user_exit = False
    save_state_on_exit = False

    @staticmethod
    def round_figures(x, n):
        """Returns x rounded to n significant figures."""
        return round(x, int(n - np.ceil(np.log10(abs(x)))))

    @staticmethod
    def time_string(seconds):
        """Returns time in seconds as a string formatted HHHH:MM:SS."""
        s = int(round(seconds))  # round to nearest second
        h, s = divmod(s, 3600)   # get hours and remainder
        m, s = divmod(s, 60)     # split remainder into minutes and seconds
        return '%4i:%02i:%02i' % (h, m, s)

    def __init__(self, initial_state=None, load_state=None):
        if initial_state:
            self.state = self.copy_state(initial_state)
        elif load_state:
            with open(load_state, 'rb') as fh:
                self.state = pickle.load(fh)
        else:
            raise ValueError('No valid values supplied for neither \
            initial_state nor load_state')

        signal.signal(signal.SIGINT, self.set_user_exit)

    def save_state(self, fname=None):
        """Saves state"""
        if not fname:
            date = datetime.datetime.now().isoformat().split(".")[0]
            fname = date + "_energy_" + str(self.energy()) + ".state"
        print("Saving state to: %s" % fname)
        with open(fname, "w") as fh:
            pickle.dump(self.state, fh)

    @abc.abstractmethod
    def move(self):
        """Create a state change"""
        pass

    @abc.abstractmethod
    def energy(self):
        """Calculate state's energy"""
        pass

    def set_user_exit(self, signum, frame):
        """
        Raises the user_exit flag, further iterations are stopped
        """
        self.user_exit = True

    def set_schedule(self, schedule):
        """
        Takes the output from `auto` and sets the attributes
        """
        self.Tmax = schedule['tmax']
        self.Tmin = schedule['tmin']
        self.steps = int(schedule['steps'])
        self.updates = schedule['updates']

    def copy_state(self, state):
        """
        Returns an exact copy of the provided state
        Implemented according to self.copy_strategy, one of
        * deepcopy : use copy.deepcopy (slow but reliable)
        * slice: use list slices (faster but only works if state is list-like)
        * method: use the state's copy() method
        """
        if self.copy_strategy == 'deepcopy':
            return copy.deepcopy(state)
        elif self.copy_strategy == 'slice':
            return state[:]
        elif self.copy_strategy == 'method':
            return state.copy()

    def update(self, step, T, E, acceptance, improvement):
        """Prints the current temperature, energy, acceptance rate,
        improvement rate, elapsed time, and remaining time.
        The acceptance rate indicates the percentage of moves since the last
        update that were accepted by the Metropolis algorithm.  It includes
        moves that decreased the energy, moves that left the energy
        unchanged, and moves that increased the energy yet were reached by
        thermal excitation.
        The improvement rate indicates the percentage of moves since the
        last update that strictly decreased the energy.  At high
        temperatures it will include both moves that improved the overall
        state and moves that simply undid previously accepted moves that
        increased the energy by thermal excititation.  At low temperatures
        it will tend toward zero as the moves that can decrease the energy
        are exhausted and moves that would increase the energy are no longer
        thermally accessible."""

        elapsed = time.time() - self.start
        if step == 0:
            print(' Temperature        Energy    Accept   Improve     Elapsed   Remaining')
            sys.stdout.write('\r%12.2f  %12.4f                      %s            ' % \
                (T, E, self.time_string(elapsed)))
            sys.stdout.flush()
        else:
            remain = (self.steps - step) * (elapsed / step)
            sys.stdout.write('\r%12.2f  %12.4f  %7.2f%%  %7.2f%%  %s  %s' % \
            (T, E, 100.0 * acceptance, 100.0 * improvement,\
            self.time_string(elapsed), self.time_string(remain))),
            sys.stdout.flush()

    def anneal(self):
        """Minimizes the energy of a system by simulated annealing.
        Parameters
        state : an initial arrangement of the system
        Returns
        (state, energy): the best state and energy found.
        """
        niters = 0
        step = 0
        self.start = time.time()

        # Precompute factor for exponential cooling from Tmax to Tmin
        if self.Tmin <= 0.0:
            raise Exception('Exponential cooling requires a minimum "\
                "temperature greater than zero.')
        Tfactor = -np.log(self.Tmax / self.Tmin)

        # Note initial state
        T = self.Tmax
        E = self.energy()
        prevState = self.copy_state(self.state)
        prevEnergy = E
        bestState = self.copy_state(self.state)
        bestEnergy = E
        trials, accepts, improves = 0, 0, 0
        if self.updates > 0:
            updateWavelength = self.steps / self.updates
            self.update(step, T, E, None, None)

        # Attempt moves to new states
        while step < self.steps and not self.user_exit:
            step += 1
            T = self.Tmax * np.exp(Tfactor * step / self.steps)
            self.move()
            E = self.energy()
            dE = E - prevEnergy
            trials += 1
            if dE > 0.0 and np.exp(-dE / T) < np.random.random():
                # Restore previous state
                self.state = self.copy_state(prevState)
                E = prevEnergy
            else:
                # Accept new state and compare to best state
                accepts += 1
                if dE < 0.0:
                    improves += 1
                prevState = self.copy_state(self.state)
                prevEnergy = E
                if E < bestEnergy:
                    bestState = self.copy_state(self.state)
                    bestEnergy = E
            if self.updates > 1:
                if step // updateWavelength > (step - 1) // updateWavelength:
                    self.update(
                        step, T, E, accepts / trials, improves / trials)
                    trials, accepts, improves = 0, 0, 0
            niters += 1

        # line break after progress output
        print('')

        self.state = self.copy_state(bestState)
        if self.save_state_on_exit:
            self.save_state()
        # Return best state and energy
        return (bestState, bestEnergy, niters)

    def auto(self, minutes, steps=2000):
        """Minimizes the energy of a system by simulated annealing with
        automatic selection of the temperature schedule.
        Keyword arguments:
        state -- an initial arrangement of the system
        minutes -- time to spend annealing (after exploring temperatures)
        steps -- number of steps to spend on each stage of exploration
        Returns the best state and energy found."""

        def run(T, steps):
            """
            Anneals a system at constant temperature and returns the state,
            energy, rate of acceptance, and rate of improvement.
            """
            E = self.energy()
            prevState = self.copy_state(self.state)
            prevEnergy = E
            accepts, improves = 0, 0
            for step in range(steps):
                self.move()
                E = self.energy()
                dE = E - prevEnergy
                if dE > 0.0 and np.exp(-dE / T) < np.random.random():
                    self.state = self.copy_state(prevState)
                    E = prevEnergy
                else:
                    accepts += 1
                    if dE < 0.0:
                        improves += 1
                    prevState = self.copy_state(self.state)
                    prevEnergy = E
            return E, float(accepts) / steps, float(improves) / steps

        step = 0
        self.start = time.time()

        # Attempting automatic simulated anneal...
        # Find an initial guess for temperature
        T = 0.0
        E = self.energy()
        self.update(step, T, E, None, None)
        while T == 0.0:
            step += 1
            self.move()
            T = abs(self.energy() - E)

        # Search for Tmax - a temperature that gives 98% acceptance
        E, acceptance, improvement = run(T, steps)

        step += steps
        while acceptance > 0.98:
            T = self.round_figures(T / 1.5, 2)
            E, acceptance, improvement = run(T, steps)
            step += steps
            self.update(step, T, E, acceptance, improvement)
        while acceptance < 0.98:
            T = self.round_figures(T * 1.5, 2)
            E, acceptance, improvement = run(T, steps)
            step += steps
            self.update(step, T, E, acceptance, improvement)
        Tmax = T

        # Search for Tmin - a temperature that gives 0% improvement
        while improvement > 0.0:
            T = self.round_figures(T / 1.5, 2)
            E, acceptance, improvement = run(T, steps)
            step += steps
            self.update(step, T, E, acceptance, improvement)
        Tmin = T

        # Calculate anneal duration
        elapsed = time.time() - self.start
        duration = self.round_figures(int(60.0 * minutes * step / elapsed), 2)

        print('') # New line after auto() output
        # Don't perform anneal, just return params
        return {'tmax': Tmax, 'tmin': Tmin, 'steps': duration}

class TSPSolver(BaseAnnealer):

    """
    Test annealer with a travelling salesman problem.
    """

    # pass extra data (the distance matrix) into the constructor
    def __init__(self, state, distance_matrix):
        self.distance_matrix = distance_matrix
        super(TSPSolver, self).__init__(state)  # important!

    def move(self):
        """Swaps two cities in the route."""
        a = np.random.randint(0, len(self.state) - 1)
        b = np.random.randint(0, len(self.state) - 1)
        self.state[a], self.state[b] = self.state[b], self.state[a]

    def energy(self):
        """Calculates the length of the route."""
        e = 0
        for i in range(len(self.state)):
            e += self.distance_matrix[self.state[i-1]][self.state[i]]
        return e

if __name__ == '__main__':
    main()
