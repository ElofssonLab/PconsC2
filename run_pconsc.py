#!/usr/bin/env python
import sys, os
import shutil
import subprocess
import multiprocessing

sys.stderr.write("""
***********************************************************************
          PconsC2 : Improved contact predictions using the recognition 
                    of protein like contact patterns.
***********************************************************************

If you use PconsC2 for contact prediction, please cite:

Skwark M.J., Raimondi D., Michel M.. and Elofsson A.  
"Improved contact predictions using the recognition of protein
 like contact patterns."
-----------------------------------------------------------------------

""")

if len(sys.argv) < 4:
    sys.stderr.write('Incorrect number of command line arguments.\n')
    sys.stderr.write('Usage: ' + sys.argv[0] + ' [-c n_cores] [-p n_jobs] [--pconsc1] <hhblits db> <jackhmmer db> <sequence file>\n\n')
    sys.exit(0)

from localconfig import *
from src import predict_all
from src import predict_all_2

sys.stderr.write('\nTesting dependencies...\n')

### Check Jackhmmer ###
try:
    f = open(os.devnull, "w")
    x  = subprocess.call([jackhmmer, '-h'], stdout=f, stderr=f)
    f.close()
except Exception as e:
    sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
    sys.stderr.write('Chosen jackhmmer binary does not seem to work!\n')
    sys.exit(1)

### Check HHblits ###
try:
    f = open(os.devnull, "w")
    x  = subprocess.call([hhblits, '-h'], stderr=f, stdout=f)
    f.close()
    pass
except:
    sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
    sys.stderr.write('Chosen HHblits binary does not seem to work!\n')
    sys.exit(1)

### Check PSICOV ###
try:
    f = open(os.devnull, "w")
    x  = subprocess.call([psicov, root + '/extras/psicovtest.fas'], stdout=f, stderr=f)
    f.close()
except Exception as e:
    sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
    sys.stderr.write('Chosen PSICOV binary does not seem to work!\n')
    sys.exit(1)

if x == 255 and not psicovfail:
    sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
    sys.stderr.write('Your version of PSICOV refuses to handle low-complexity alignments.\n')
    sys.stderr.write('We recommend patching the PSICOV code to allow this. See 00README\n')
    sys.stderr.write('If you _really_ do not want to do that, please change psicovfail flag in \n')
    sys.stderr.write(os.path.abspath(sys.argv[0]) + ' to True.\n')
    sys.stderr.write('This will (most probably) affect the prediction performance.\n')
    sys.exit(1)


### Check plmDCA ###
if plmdca:
    try:
        f = open(os.devnull, "w")
        x  = subprocess.call([plmdca, '-h'], stdout=f, stderr=f)
        f.close()
    except Exception as e:
        sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
        sys.stderr.write('Chosen plmdca binary does not seem to work!\n')
        sys.exit(1)
elif matlab:
    try:
        f = open(os.devnull, "w")
        x  = subprocess.call([matlab, '-h'], stdout=f, stderr=f)
        f.close()
    except:
        sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
        sys.stderr.write('Chosen MATLAB binary does not seem to work!\n')
        sys.stderr.write('You can get MCR \n')
        sys.stderr.write('http://www.mathworks.se/products/compiler/mcr/\n')
        sys.exit(1)
else:
    sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
    sys.stderr.write('You must set one of plmdca or matlab in localconfig.py!\n')
    sys.exit(1)


sys.stderr.write('Dependencies OK.\n')


### parse parameters
pconsc1_flag = False
n_jobs_plm = min(2, n_cores)
n_jobs_psi = min(2, n_cores)

if '-c' in sys.argv:
    idx = sys.argv.index('-c')
    try:
        n_cores = int(sys.argv[idx+1])
    except:
        print 'Number of cores -c must be an integer, %r is not. Default is %s.' % (sys.argv[idx+1], n_cores)
        sys.exit(1)
    del sys.argv[idx]
    del sys.argv[idx]

if '--p_plm' in sys.argv:
    idx = sys.argv.index('--p_plm')
    try:
        n_jobs_plm = int(sys.argv[idx+1])
    except:
        print 'Number of parallel plmDCA jobs --p_plm must be an integer, %r is not. Default is %s.' % (sys.argv[idx+1], n_jobs_plm)
        sys.exit(1)
    del sys.argv[idx]
    del sys.argv[idx]

if '--p_psi' in sys.argv:
    idx = sys.argv.index('--p_psi')
    try:
        n_jobs_psi = int(sys.argv[idx+1])
    except:
        print 'Number of parallel PSICOV jobs -p_psi must be an integer, %r is not. Default is %s.' % (sys.argv[idx+1], n_jobs_psi)
        sys.exit(1)
    del sys.argv[idx]
    del sys.argv[idx]

if '--pconsc1' in sys.argv:
    idx = sys.argv.index('--pconsc1')
    pconsc1_flag = True
    del sys.argv[idx]

hhblitsdb = os.path.abspath(sys.argv[1])
jackhmmerdb = os.path.abspath(sys.argv[2])
seqfile = os.path.abspath(sys.argv[3])

shutil.copyfile(root + 'localconfig.py', root + 'src/localconfig.py')

predict_all.main(hhblitsdb, jackhmmerdb, seqfile, n_cores=n_cores, n_jobs_plm=n_jobs_plm, n_jobs_psi=n_jobs_psi, pconsc1_flag=pconsc1_flag)
#predict_all_2.main(hhblitsdb, jackhmmerdb, seqfile, n_cores=n_cores)
