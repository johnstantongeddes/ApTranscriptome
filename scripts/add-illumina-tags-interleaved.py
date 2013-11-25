#! /usr/bin/env python
import screed, sys, itertools

srfile = sys.argv[1]

n = 1

for read in screed.open(srfile):
    name = read.name 
    if n%2!=0:
	    name += '/1'
    else:
	    name += '/2'
    n += 1
    print '@%s\n%s\n+\n%s' % (name,
                              read.sequence,
                              read.accuracy,
                             )
