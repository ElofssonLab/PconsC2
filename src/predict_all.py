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


def main(hhblitsdb, jackhmmerdb, seqfile, n_cores=1, layers=5, pconsc1_flag=False):

    prep.run_alignments(hhblitsdb, jackhmmerdb, seqfile, n_cores=n_cores)
   
    psicovnames = prep.run_contact_pred(seqfile, 'psicov', n_cores=n_cores)
    plmdcanames = prep.run_contact_pred(seqfile, 'plmdca', n_cores=n_cores)
    predictionnames = dict(psicovnames.items() + plmdcanames.items())

    l = [root + '/src/predict2.py']
    names = ['jhE4', 'jhE0', 'jhE10', 'jhE40', 'hhE4', 'hhE0', 'hhE10', 'hhE40']
    for key in names:
        l.append(predictionnames[key + 'psicov'])
        l.append(predictionnames[key + 'plmdca'])

    if pconsc1_flag:
        sys.stderr.write("Running PconsC1...\n")
        sys.stderr.write("Predicting...\n")
        result_name = seqfile + '.pconsc.out'
        l[0] = root + '/src/predict.py'
        results = check_output(l)

        f = open(result_name, 'w')
        f.write(results)
        f.close()
    else:
        sys.stderr.write("Running PconsC2...\n")
        netsurfpredictionname, sspredictionname, pssmaliname = prep.run_pconsc2_dependencies(hhblitsdb, seqfile, n_cores=1)
        sys.stderr.write("Predicting...\n")
        result_name = seqfile + '.pconsc2.out'
        l.extend([netsurfpredictionname, sspredictionname, pssmaliname, result_name])
        check_output(l)

    # plot the top L*1 contacts in a contact map
    # those contacts are later used during protein folding
    if plot_flag:
        if os.path.exists('native.pdb') and os.path.exists(seqfile + '.horiz'):
            plot_map(seqfile, result_name, 1., pdb_filename='native.pdb', psipred_filename=seqfile + '.horiz')
        elif os.path.exists('native.pdb'):
            plot_map(seqfile, result_name, 1., pdb_filename='native.pdb')
        elif os.path.exists(seqfile + '.horiz'):
            plot_map(seqfile, result_name, 1., psipred_filename=seqfile + '.horiz')
        else:
            plot_map(seqfile, result_name, 1.)


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

    # number of layers to compute (deep learning)
    # must be between 0 and 4 as there are no further trained forests given (layerX.forest)
    layers = int(sys.argv[4])

    main(hhblitsdb, jackhmmerdb, seqfile, n_cores=n_cores, layers=layers)



