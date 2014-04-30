#!/usr/bin/env python

"""Script to predict contacts, given the 16 input files"""

import numpy as np
import pickle, sys
import os

if len(sys.argv) != 21:
	print 'Usage: ' + sys.argv[0] + ' <input files> <output file>'
	print 'Output files need to come in *order*!'
	print 'That is:'
	print ' JackHMMER 1e-4 Psicov'
	print ' JackHMMER 1e-4 plmDCA'
	print ' JackHMMER 1e-0 Psicov'
	print ' JackHMMER 1e-0 plmDCA'
	print ' JackHMMER 1e-10 Psicov'
	print ' JackHMMER 1e-10 plmDCA'
	print ' JackHMMER 1e-40 Psicov'
	print ' JackHMMER 1e-40 plmDCA'
	print ' HHblits 1e-4 Psicov'
	print ' HHblits 1e-4 plmDCA'
	print ' HHblits 1e-0 Psicov'
	print ' HHblits 1e-0 plmDCA'
	print ' HHblits 1e-10 Psicov'
	print ' HHblits 1e-10 plmDCA'
	print ' HHblits 1e-40 Psicov'
	print ' HHblits 1e-40 plmDCA'
	print ' NetSurf RSA'
	print ' PSIPRED'
	print ' Alignment for PSSM evaluation (e.g. HHblits at E=1)'
	print ' Output file'
	sys.exit(1)


def predict(X, forest):
	probability = []
	for t in range(len(forest)):
		tree = forest[t]
		while len(tree) > 2:
			if X[tree[0][0]] <= tree[0][1]:
				tree = tree[1]
			else:
				tree = tree[2]
		probability.append(tree[1]/float(tree[0] + tree[1]))
	return sum(probability)/len(probability)

def parsePSIPRED(f):
	SSdict = {}
	try:
		x = open(f).read().split('\n')
	except:
		return SSdict
	for l in x:
		y = l.split()
		if len(y) != 6:
			continue
		i = int(y[0])
		SSdict[i] = [float(y[3]), float(y[4]), float(y[5])]
	return SSdict
		
def parseNetSurfP(f):
	netSurfdict = {}
	for l in open(f).readlines():
		al = []
		x = l.split()
		if l.find('#') == 0:
			continue
		if l[0] not in ['B', 'E']:
			y = ['E']
			y.extend(x)
			x = y
		for y in [4,6,7, 8, 9]:
			al.append(float(x[y]) )
		netSurfdict[ int(x[3] )] = al
	return netSurfdict

def parsePSSM(alignment):
	pssm = {}
	one2number = 'ARNDCEQGHILKMFPSTWYV-'
	bi = [ 0.0825, 0.0553, 0.0406, 0.0545, 0.0137, 0.0393, 0.0675, 0.0707, 0.0227, 0.0595, 0.0966, 0.0584, 0.0242, 0.0386, 0.0470, 0.0657, 0.0534, 0.0108, 0.0292, 0.0687 ]
	b = {}
	for i in one2number[:-1]:
		b[i] = bi[one2number.find(i)] 

	freqs = {}
	seqcount = 0.
	gapcount = 0
	for l in open(alignment):
		if l.find('>') > -1:
			continue
		x = l.strip()
		seqcount += 1
		for i in range(len(x)):
			try:
				freqs[i][x[i]] += 1 
			except:
				try:
					freqs[i][x[i]] = 1 
				except:
					freqs[i] = {} 
					freqs[i][x[i]] = 1 
			if x[i] == '-':
				gapcount += 1

	b['-'] = gapcount/(seqcount * len(freqs.keys()))

	for i in sorted(freqs.keys()):
		q = []
		for l in one2number:
			try:
				q.append(np.log( freqs[i][l]/(b[l] * seqcount)))
			except Exception as e:
				q.append(np.log( 0.1/(b[l] * seqcount)))
		pssm[i+1] = q
	return pssm

files = sys.argv[1:]
X = []
Y = []
maxres = -1
acceptable = []
sys.stderr.write('Parsing NetSurfP\n')
accessibility = parseNetSurfP(files[16])
sys.stderr.write('Parsing PSIPRED\n')
SSdict = parsePSIPRED(files[17])
sys.stderr.write('Parsing PSSM\n')
pssm = parsePSSM(files[18])
contacts = {}
pairs = set()

for index in range(16):
	contacts[index] = {}
	d = files[index]
	sys.stderr.write('Parsing ' + d + '\n')
	if not os.path.exists(d):
		sys.stderr.write(d + ' does not exist!\n')
		sys.exit(1)
	infile = open(d).readlines()
	for m in infile:
		if d.find('psicov') > -1:
			x = m.split()
			if len(x) != 5:
				print d + ' has wrong format!'
				sys.exit(1)
			c = 4
		elif d.find('plmdca') > -1:
			x = m.split(',')
			if len(x) != 3:
				print d + ' has wrong format!'
				sys.exit(1)
			c = 2
		aa1 = int(x[0])
		aa2 = int(x[1])
		if aa1 > maxres:
			maxres = aa1
		if aa2 > maxres:
			maxres = aa2	
		if x[c].find('nan') > -1:
			score = 0
		else:
			score = float(x[c])
		contacts[index][(aa1, aa2)] = score
		pairs.add((aa1,aa2))

sys.stderr.write('Coalescing data')
pairlist = sorted(list(pairs))
count = 0.
for s in pairlist:
	count += 1
	q = []
	q.append(abs(s[0]-s[1]))
	for index in range(16):
		try:
			q.append(contacts[index][s])
		except:
			q.append(0)

        for i in range(-4, 5):
                try:
                        q.extend(SSdict[s[0]] )
                except:
                        q.extend([0,0,0])

        for i in range(-4, 5):
                try:
                        q.extend(SSdict[s[1]] )
                except:
                        q.extend([0,0,0])
	
	for i in range(-4, 5):
		try:
			q.extend(accessibility[s[0]] )
		except:
			q.extend([0,0,0,0,0])
	
	for i in range(-4, 5):
		try:
			q.extend(accessibility[s[1]] )
		except:
			q.extend([0,0,0,0,0])

	q.extend(pssm[s[0]])
	q.extend(pssm[s[1]])

	X.append(q)
	Y.append(s)

sys.stderr.write('\nLoading RF 0\n')
forest = pickle.load(open(os.path.dirname(os.path.abspath(sys.argv[0])) + '/layer0.dat'))
previouslayer = {}
sys.stderr.write('\nPredicting layer 0\n')
for l in range(len(Y)):
	(aa1, aa2) = (Y[l][0], Y[l][1])
	previouslayer[(aa1, aa2)] = predict(X[l], forest)
#	print 'L0', aa1, aa2, '{:6.4f}'.format(previouslayer[(aa1, aa2)])

sys.stderr.write('Coalescing data\n')
for layer in range(1,6):
	X = []
	Y = []
	count = 0
	for s in pairlist:
		count += 1
		q = []
		q.append(abs(s[0]-s[1]))
		for index in range(16):
			try:
				q.append(contacts[index][s])
			except:
				q.append(0)

        	for i in range(-4, 5):
        	        try:
        	                q.extend(SSdict[s[0]] )
        	        except:
        	                q.extend([0,0,0])

	        for i in range(-4, 5):
	                try:
	                        q.extend(SSdict[s[1]] )
	                except:
        	                q.extend([0,0,0])
	
		for i in range(-4, 5):
			try:
				q.extend(accessibility[s[0]] )
			except:
				q.extend([0,0,0,0,0])
	
		for i in range(-4, 5):
			try:
				q.extend(accessibility[s[1]] )
			except:
				q.extend([0,0,0,0,0])

		q.extend(pssm[s[0]])
		q.extend(pssm[s[1]])
                for i in range(-5,6):
                        for j in range(-5, 6):
                                try:
                                     	q.append(previouslayer[(s[0] + i, s[1] + j)])
                                except Exception as e:
#					sys.stderr.write(str(e) + '\n')
                                        q.append(-3)
		X.append(q)
		Y.append(s)

	sys.stderr.write('Loading RF {:d}\n'.format(layer))
	forest = pickle.load(open(os.path.dirname(os.path.abspath(sys.argv[0])) + '/layer{:d}.dat'.format(layer)))
	previouslayer = {}
	sys.stderr.write('Predicting layer {:d}\n'.format(layer))
	for l in range(len(Y)):
		(aa1, aa2) = (Y[l][0], Y[l][1])
		previouslayer[(aa1, aa2)] = predict(X[l], forest)
#		print 'L{:d}'.format(layer), aa1, aa2, '{:6.4f}'.format(previouslayer[(aa1, aa2)])

f = open(files[19], 'w')
for l in range(len(Y)):
	(aa1, aa2) = (Y[l][0], Y[l][1])
	f.write('{:d} {:d} {:6.4f}\n'.format(aa1, aa2, previouslayer[(aa1, aa2)]))

