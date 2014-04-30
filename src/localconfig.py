#!/usr/bin/env python
import os
import sys
import multiprocessing
import subprocess

if __name__ == '__main__':
    print 'Please do not run me! Use run_pconsc.py'
    print '\n\tYours sincerely,\n\n\t', sys.argv[0]
    sys.exit(0)

# Directory where PconsC in the distributable package is located
root = os.path.dirname(os.path.abspath(sys.argv[0])) + '/'


########################################################
### Please adjust the following paths to your system ###
########################################################

### Jackhmmer executable ###
jackhmmer = root + 'dependencies/hmmer-3.0/src/jackhmmer'

### HHblits executable ###
hhblits = root + 'dependencies/hhsuite-2.0.16/bin/hhblits'

### PSICOV executable ###
psicov = root + 'dependencies/psicov-1.11/psicov'

### NetSurfP executable ###
netsurf = root + 'dependencies/netsurfp-1.0/netsurfp'

### PSIPRED executable ###
psipred = root + 'dependencies/psipred/runpsipred'

### MATLAB executable ###
# Please set this variable to None if you don't have access to matlab. 
# PconsFold will then try to use the compiled version. 
#matlab = '/sw/apps/matlab/x86_64/8.1/bin/matlab'
matlab = None

### Path to MATLAB compiler ###
# Only needed if matlab is not available.
matlabdir = '' 



########################################################
###  Please do not change anything below this line   ###
########################################################

# Paths to included scripts
trim2jones = root + 'scripts/a3mToJones.py'
trim2trimmed = root + 'scripts/a3mToTrimmed.py'
#trim = root + 'scripts/trim.py'
#trim2 = root + 'scripts/trimToFasta.py'

# Reformat script scavenged from HHsuite. Please cite the HHblits paper!
reformat = root + 'scripts/reformat.pl'

# Maximum amount of cores to use per default
n_cores = multiprocessing.cpu_count()

# Enable work-around for PSICOV not handling low complexity alignments
psicovfail = True

# Adjust plmdca path to either standalone or compiled, 
# depending on presence of matlab.
if matlab:
    plmdca = None # matlab licence present: do not use compiled version
    plmdcapath = root + 'dependencies/plmDCA_symmetric_v2'
else:
    plmdca = root + 'dependencies/plmdca_standalone/2012/build01/bin/plmdca'
    plmdcapath = None

