#!/usr/bin/env python
import sys

sys.path.append("/scratch/mirco_local/bioinfo-toolbox")
from parsing import parse_fasta
from parsing import parse_contacts

if len(sys.argv) != 5:
    sys.stderr.write('Incorrect number of command line arguments.\n')
    sys.stderr.write('Usage: ' + sys.argv[0] + ' <sequence file> <contact file> <CASP target ID> <output filename>\n\n')
    sys.exit(0)


sfile = sys.argv[1]
cfile = sys.argv[2]
target = sys.argv[3]

seq = parse_fasta.read_fasta(open(sfile)).items()[0][1][0]

contacts = parse_contacts.parse(open(cfile), min_dist=0)

print len(contacts)
print contacts[0]
print seq

ofile = open(sys.argv[4], 'w')

ofile.write("PFRMAT RR\nTARGET %s\nAUTHOR 6685-2065-9124\nMETHOD Pcons-net\nREMARK PconsC2\nMETHOD Improved contact predictions using the\nMETHOD recognition of protein like contact\nMETHOD patterns.\nMODEL  1\n" % target)

tmp_i = 1
for aa in seq:
    ofile.write(aa)
    if tmp_i % 50 == 0:
        ofile.write("\n")
    tmp_i += 1

if tmp_i % 50 != 0:
    ofile.write("\n")

count = 0
for c in contacts:
    ofile.write("%d %d 0 8 %s\n" % (c[1], c[2], c[0]))
    count += 1
    if count > 31998:
        break


ofile.write("END")
