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

### Working directory ###
rundir = os.getcwd() + '/'  ## Manual addition

if not 'pcons2fold' in root:
    root += 'pcons2fold/'

# Look if PconsC or Rosetta dependencies need to be checked
rosetta_flag = False
if 'rosetta' in root:
    rosetta_flag = True
    root = '/'.join(root.split('/')[:-3]) + '/pconsc/'

########################################################
### Please adjust the following paths to your system ###
########################################################

### Path to root folder of your Rosetta installation ###
# REQUIRES: Rosetta 3.5 or weekly build
rosettadir = '/home/felix/Software/rosetta_2014.34.57213_bundle'

### Jackhmmer executable ###
jackhmmer = '/usr/bin/jackhmmer'

### HHblits executable ###
hhblits = '/usr/bin/hhblits'

### PSICOV executable ###
psicov = '/usr/local/bin/psicov'

### NetSurfP executable ###
netsurf = '/home/felix/Software/netsurfp-1.0/netsurfp'

### PSIPRED executable ###
psipred = '/home/felix/Software/psipred/bin/psipred'

### Path to TM-score ###
# only needed if result should be compared to native structure
tmscore_binary = '/usr/bin/TMscore'

### MATLAB executable ###
# Please set this variable to None if you don't have access to matlab.
# PconsFold will then try to use the compiled version.
#matlab = '/sw/apps/matlab/x86_64/8.1/bin/matlab'
matlab = '/home/felix/Software/MATLAB/R2014b/bin/matlab'

### Path to MATLAB compiler ###
# Only needed if matlab is not available.
matlabdir = ''

########################################################
###  Please do not change anything below this line   ###
########################################################

# Internal Rosetta paths
rosetta_db_dir = rosettadir + '/main/database'
rosetta_binary_dir = rosettadir + '/main/source/bin'
rosetta_make_fragments = rosettadir + '/tools/fragment_tools/make_fragments.pl'
rosetta_abinitiorelax = rosetta_binary_dir + '/AbinitioRelax.linuxgccrelease'
rosetta_extract = rosetta_binary_dir + '/extract_pdbs.linuxgccrelease'
rosetta_relax = rosetta_binary_dir + '/relax.linuxgccrelease'

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
    plmdcapath = '/home/felix/Software/plmDCA_symmetric_v3'
    #plmdcapath = '/home/felix/Software/plmDCA_asymmetric_v2'
else:
    plmdca = root + 'dependencies/plmdca_standalone/2012/build01/bin/plmdca'
    plmdcapath = None
