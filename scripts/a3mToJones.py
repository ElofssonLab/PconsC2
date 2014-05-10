#!/usr/bin/env python

import sys, os

infilef = sys.argv[1]
infile = open(infilef)

count = 0
for l in infile:
    if '>' in l:
        if count != 0:
            sys.stdout.write('\n')
        count += 1
        continue
    l = l.strip()
    upperseq = ''.join([c for c in l if not c.islower()])
    upperseq = upperseq.replace('X', '-')
    sys.stdout.write(upperseq)
