#!/usr/bin/python

import sys, subprocess, os

pdbid = sys.argv[1]
pdb = pdbid.split(':')[0]

if os.path.exists(pdbid + '/netSurf.local'):
	sys.exit(0)

subprocess.call('/home/mjs/sw/netsurfp-1.0/netsurfp -a -i ' + pdbid + '/sequence.fa > ' + pdbid + '/netSurf.local', shell=True) 
