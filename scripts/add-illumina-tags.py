#! /usr/bin/env python
import screed, sys, itertools

srfile = sys.argv[1]
rnum = sys.argv[2]

n = 1

for read in screed.open(srfile):
    name = read.name 
	
    name += '/' + str(rnum)
    n += 1
    print '@%s\n%s\n+\n%s' % (name,
                              read.sequence,
                              read.accuracy,
                             )
