#!/usr/bin/env python
from localconfig import *
import prepare_input as prep
from datetime import datetime
import sys, subprocess, os


plot_flag = False
try:
    from plotting.plot_contact_map import plot_map
    plot_flag = True
except:
    sys.stderr.write('\nWARNING:\nBiopython or Matplotlib not available, skip contact map plotting.\n')
    pass


def check_output(command):
    return subprocess.Popen(command, stdout=subprocess.PIPE).communicate()[0]


def main(hhblitsdb, jackhmmerdb, seqfile, n_cores=1):
   
    prep.run_alignments(hhblitsdb, jackhmmerdb, seqfile, n_cores=1)
   
    jhpredictionnames = []
    hhpredictionnames = []
    jhpredictionnames, hhpredictionnames = prep.run_psicov(seqfile, jhpredictionnames, hhpredictionnames, n_cores=n_cores)
    jhpredictionnames, hhpredictionnames = prep.run_plmdca(seqfile, jhpredictionnames, hhpredictionnames, n_cores=n_cores)

    sys.stderr.write("Predicting...\n")
    #l = [os.path.dirname(os.path.abspath(sys.argv[0])) + '/predict.py']
    l = [root + '/src/predict.py']
    l.extend(jhpredictionnames)
    l.extend(hhpredictionnames)
    results = check_output(l)

    f = open(seqfile + '.pconsc.out', 'w')
    f.write(results)
    f.close()


    # plot the top L*1 contacts in a contact map
    # those contacts are later used during protein folding
    if plot_flag:
        if os.path.exists('native.pdb') and os.path.exists(seqfile + '.horiz'):
            plot_map(seqfile, seqfile + '.pconsc.out', 1.0, pdb_filename='native.pdb', psipred_filename=seqfile + '.horiz')
        elif os.path.exists('native.pdb'):
            plot_map(seqfile, seqfile + '.pconsc.out', 1.0, pdb_filename='native.pdb')
        elif os.path.exists(seqfile + '.horiz'):
            plot_map(seqfile, seqfile + '.pconsc.out', 1.0, psipred_filename=seqfile + '.horiz')
        else:
            plot_map(seqfile, seqfile + '.pconsc.out', 1.0)


if __name__ == "__main__":

    ### parse parameters

    if len(sys.argv) < 4:
        print sys.argv[0], '[-c n_cores] <hhblits db> <jackhmmer db> <sequence file>'
        sys.exit(0)

    if '-c' in sys.argv:
        idx = sys.argv.index('-c')
        try:
            n_cores = int(sys.argv[idx+1])
        except:
            print 'Number of cores -c must be an integer, %r is not. Default is %s.' % (sys.argv[idx+1], n_cores)
            sys.exit(1)
        del sys.argv[idx +1]
        del sys.argv[idx]

    hhblitsdb = sys.argv[1]
    jackhmmerdb = sys.argv[2]
    seqfile = sys.argv[3]

    main(hhblitsdb, jackhmmerdb, seqfile, n_cores=n_cores)
