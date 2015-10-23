# Attempt to write a simple time series generator
# Should sample uniformly at each time step within range
# Range ultimately determined by empirical data sets
# Value at step n+1 not constrained by value at step n (!!!)

import random
from sys import argv

script, start, end, ns, r = argv

time_series = []
time_steps = []
start, end, ns, r = int(argv[1]), int(argv[2]), int(argv[3]), int(argv[4])


# produces a list, the items of which are n time steps in the interval [x, y)
def steps(x, y, n):
    inc = (y - x) / float(n)
    while x < y:
        time_steps.append(x)
        x = x + inc


# generates values for the time series from time_steps
# sample is drawn from uniform distribution from -t to t
# produces time_series list (vector of length = number of steps)
def val(list_name):
    for t in list_name:
        time_series.append(random.uniform(-r*t, r*t))


steps(start, end, ns)

val(time_steps)

# outputs time series as a line-by-line display
# (could be a list/vector but this is a dummy script for ease of reading)
print "Your sample time series is:"

for value in time_series:
    print "\t", value