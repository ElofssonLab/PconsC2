#!/cvmfs/fgi.csc.fi/apps/sl6/python/2.7.3/bin/python
#!/pdc/vol/python/2.7.3/bin/python

import sys, subprocess, os
import string as s

def check_output(command):
	return subprocess.Popen(command, stdout=subprocess.PIPE).communicate()[0]
	
#root = '/triton/ics/work/mjs/pconse/'
root = '/triton/ics/work/mjs/pconse/'
hhblits = '/triton/ics/work/mjs/hhsuite-2.0.16-linux-x86_64/bin/hhblits'
trim = root + 'scripts/trim.py'
trim2 = root + 'scripts/trimToFasta.py'
trim3 = root + 'scripts/trimToFastaDCA.py'
reformat = root + 'scripts/reformat.pl'
psicov= root + 'bin/psicov'
# uniref100 = '/home/mjs/db/uniref/MessyUniref100'
uniref100 = '/triton/ics/work/mjs/uniref100/uniref100.fasta'
matlab = '/share/apps/matlab/R2012b/bin/matlab'
dca = root + 'bin/dca.m'

if len(sys.argv) < 2:
	print sys.argv[0], '<target>'
	sys.exit(0)

target = '9xxxx'
seqfile = sys.argv[1]
infile = seqfile

rundir = seqfile.rfind('/')
if rundir < 0:
	rundir = '.'
else:
	rundir = seqfile[:rundir]

print seqfile, sys.argv[0]


if not os.path.exists(seqfile):
	print seqfile, 'does not exist'
	sys.exit(0)

f = open(seqfile).read()

if os.path.exists(seqfile + '.fa'):
	subprocess.call(['mv', seqfile + '.fa', seqfile +'.bak'])

f2 = open(seqfile +'.fa', 'w')
if f[0] != '>':
        f2.write('>target\n' + f +'\n')
else:  
        x = f.split('\n')
        if len(x[0]) > 6:
                target = x[0][1:5] + x[0][6]
        f2.write('>target\n' + "".join(x[1:]) + '\n')
f2.close()

print 'HHblits'
if not os.path.exists(seqfile + '.hhE10.jones'):
	t = check_output([hhblits, '-n', '5', '-all', '-maxfilt', '999999', '-realign_max', '999999', '-oa3m', seqfile +'.hhE10.a3m', '-e', '1e-10', '-i', seqfile + '.fa', '-d', '/triton/ics/work/mjs/uniprot20_2013_03/uniprot20_2013_03'])
	check_output([reformat, 'a3m', 'fas', seqfile + '.hhE10.a3m', seqfile + '.hhE10.fas'])

if not os.path.exists(seqfile + '.hhE10.psicov'):
        t = check_output([trim, seqfile + '.hhE10.fas'])
	f = open(seqfile + '.hhE10.jones', 'w')
	f.write(t)
	f.close()
	print "Running psicov"
	if not os.path.exists(seqfile + '.hhE10.psicov'):
		t = check_output([psicov, '-j', '1', seqfile + '.hhE10.jones'])
		f = open(seqfile + '.hhE10.psicov', 'w')
		f.write(t)
		f.close()

t = check_output([trim2, seqfile + '.hhE10.fas'])
f = open(seqfile + '.hhE10.trimmed', 'w')
f.write(t)
f.close()

os.chdir(os.path.abspath(infile)[:os.path.abspath(infile).rfind('/')])
infilestem = infile.split('/')[-1] + '.hhE10'

if not os.path.exists(infilestem + ".plmdca"):
        print "Running plmDCA"
	t = check_output([matlab, '-nodesktop', '-r', "path(path, '/triton/ics/work/mjs/pconse/scripts/plmDCA_asymmetric'); path(path, '/triton/ics/work/mjs/pconse/scripts/plmDCA_asymmetric/functions'); path(path, '/triton/ics/work/mjs/pconse/scripts/plmDCA_asymmetric/3rd_party_code/minFunc/'); plmDCA_asymmetric ( '" + infilestem + ".trimmed', '" + infilestem + ".plmdca', 0.01, 0.01, 0.1, 4); exit"])

