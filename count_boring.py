#!/usr/bin/env python

#counts those to which C=I search strategy applies and returns that C=I indeed in these graphs 
import pyximport
pyximport.install() # automatic compliation of cython modules
import enumdags

# Settings
MAXSIZE = 7
NUMPROCESSES = 7 # number of parallel processes

# Count total and boring DAGs of size at most MAXSIZE
for n in range(1, MAXSIZE+1):
    enumdags.countdags(n, NUMPROCESSES)
