#!/usr/bin/env python
from datetime import datetime
import sys, subprocess, os, math
import multiprocessing

sys.path.append("../")
from localconfig import *

names = ['jhE4', 'jhE0', 'jhE10', 'jhE40', 'hhE4', 'hhE0', 'hhE10', 'hhE40']
cutoffs = ['1e-4', '1', '1e-10', '1e-40','1e-4', '1', '1e-10', '1e-40']


def check_output(command):
    #print ' '.join(command)
    return subprocess.Popen(command, stdout=subprocess.PIPE).communicate()[0]


def run_alignments(hhblitsdb, jackhmmerdb, seqfile, n_cores=1):

    if hhblitsdb.endswith('_a3m_db'):
        hhblitsdb = hhblitsdb[:-7]
    if not os.path.exists(hhblitsdb + '_a3m_db'):
        sys.stderr.write('\n' + hhblitsdb + '_a3m_db' + 'does not exist\n')
        sys.exit(1)
    if not os.path.exists(jackhmmerdb):
        sys.stderr.write('\n' + jackhmmerdb + 'does not exist\n')
        sys.exit(1)
    if not os.path.exists(seqfile):
        sys.stderr.write('\n' + seqfile + 'does not exist\n')
        sys.exit(0)

    f = open(seqfile).read()

    if os.path.exists(seqfile + '.fasta'):
        subprocess.call(['mv', seqfile + '.fasta', seqfile +'.bak'])

    f2 = open(seqfile +'.fasta', 'w')
    if f[0] != '>':
        f2.write('>target\n' + f +'\n')
    else:
        x = f.split('\n')
        if len(x[0]) > 6:
            target = x[0][1:5] + x[0][6]
        f2.write('>target\n' + "".join(x[1:]) + '\n')
    f2.close()

    for i in range(len(names)):

        exists_a3m = os.path.exists(seqfile + '.' + names[i] + '.a3m')
        exists_psicov = os.path.exists(seqfile + '.' + names[i] + '.psicov')
        exists_plmdca = os.path.exists(seqfile + '.' + names[i] + '.plmdca')

        # only create alignment file if at least one of the contact maps is missing
        if not exists_a3m and (not exists_psicov or not exists_plmdca):
            if 'jh' in names[i]:
                sys.stderr.write(str(datetime.now()) + ' ' + names[i] + ': generating Jackhmmer alignment\nThis may take quite a few minutes!\n ')
                t = check_output([jackhmmer, '--cpu', str(n_cores), '-N', '5', '-E', cutoffs[i], '-A', seqfile +'.' + names[i] + '.ali', seqfile + '.fasta', jackhmmerdb])
                check_output([reformat, 'sto', 'a3m', seqfile + '.' + names[i] + '.ali', seqfile + '.' + names[i] + '.a3m'])
                check_output(['rm', seqfile + '.' + names[i] + '.ali'])

            elif 'hh' in names[i]:
                sys.stderr.write(str(datetime.now()) + ' ' + names[i] + ': generating HHblits alignment\nThis may take quite a few minutes!\n ')
                t = check_output([hhblits, '-all', '-oa3m', seqfile + '.' + names[i] + '.a3m', '-e', cutoffs[i], '-cpu', str(n_cores), '-i', seqfile + '.fasta', '-d', hhblitsdb])



def run_cpred_job(i, seqfile, method, n_cores, n_jobs):

        c_fname = seqfile + '.' + names[i] + '.' + method
        exists_cfile = os.path.exists(c_fname)

        plmdca_cores = int(math.floor(n_cores/n_jobs))
        #plmdca_cores = int(math.floor(6/n_jobs))

        if not exists_cfile:
            if method == 'plmdca':
                input_fname = seqfile + '.' + names[i] + '.trimmed'
                t = check_output([trim2trimmed, seqfile + '.' + names[i] + '.a3m'])
                f = open(input_fname, 'w')
                f.write(t)
                f.close()

                sys.stderr.write(str(datetime.now()) + ' ' + names[i] + ': running plmDCA\nThis may take more than an hour.\n')
                if plmdca:
                    #t = check_output([plmdca, matlabdir, seqfile + '.jh' + names[i] + ".trimmed", seqfile + '.jh' + names[i] + ".plmdca", "0.01", "0.01", "0.1", str(n_cores)])
                    #t = check_output([plmdca, input_fname, c_fname, "0.01", "0.01", "0.1", str(n_cores)])

                    # plmDCA_symmetric
                    #t = check_output([plmdca, input_fname, c_fname, "0.01", "0.01", "0.1", str(plmdca_cores)])
                    # plmDCA_asymmetric
                    t = check_output([plmdca, input_fname, c_fname, "0.1", str(plmdca_cores)])
                else:
                    #t = check_output([matlab, '-nodesktop', '-nosplash', '-r', "path(path, '" + plmdcapath + "'); path(path, '" + plmdcapath + "/functions'); path(path, '" + plmdcapath + "/3rd_party_code/minFunc/'); plmDCA_symmetric ( '" + seqfile + '.' + names[i] + ".trimmed', '" + seqfile + '.' + names[i] + ".plmdca', 0.01, 0.01, 0.1, " + str(n_cores) + "); exit"])
                    #t = check_output([matlab, '-nodesktop', '-nosplash', '-r', "path(path, '" + plmdcapath + "'); path(path, '" + plmdcapath + "/functions'); path(path, '" + plmdcapath + "/3rd_party_code/minFunc/'); plmDCA_symmetric ( '" + seqfile + '.hh' + names[i] + ".trimmed', '" + seqfile + '.hh' + names[i] + ".plmdca', 0.01, 0.01, 0.1, " + str(n_cores) + "); exit"])
                    # plmDCA_symmetric
                    #t = check_output([matlab, '-nodesktop', '-nosplash', '-r', "path(path, '" + plmdcapath + "'); path(path, '" + plmdcapath + "/functions'); path(path, '" + plmdcapath + "/3rd_party_code/minFunc/'); plmDCA_symmetric ( '" + input_fname  + "', '" + c_fname + "', 0.01, 0.01, 0.1, " + str(plmdca_cores) + "); exit"])
                    # plmDCA_asymmetric
                    t = check_output([matlab, '-nodesktop', '-nosplash', '-r', "path(path, '" + plmdcapath + "'); path(path, '" + plmdcapath + "/functions'); path(path, '" + plmdcapath + "/3rd_party_code/minFunc/'); plmDCA_asymmetric ( '" + input_fname  + "', '" + c_fname + "', 0.1, " + str(plmdca_cores) + "); exit"])

            elif method == 'psicov':
                input_fname = seqfile + '.' + names[i] + '.jones'
                t = check_output([trim2jones, seqfile + '.' + names[i] + '.a3m'])
                f = open(input_fname, 'w')
                f.write(t)
                f.close()

                t = ''
                sys.stderr.write(str(datetime.now()) + ' ' + names[i] + ': running PSICOV\nThis may take more than an hour.\n')
                try:
                    # Joel @ NSC: Added -o flag, in case the psicov binary has not
                    # been compiled with MINEFSEQS=0.
                    t = check_output([psicov, '-o', input_fname])
                except:
                    t = ''
                f = open(c_fname, 'w')
                f.write(t)
                f.close()

            else:
                sys.stderr.write('\nWrong method string!\n Must be one of ["psicov", "plmdca"]\n')
                sys.exit(0)



def run_contact_pred(seqfile, method, n_cores=1, n_jobs=1):

    predictionnames = {}

    if not os.path.exists(seqfile):
        sys.stderr.write('\n' + seqfile + 'does not exist\n')
        sys.exit(0)

    f = open(seqfile).read()

    if os.path.exists(seqfile + '.fasta'):
        subprocess.call(['mv', seqfile + '.fasta', seqfile +'.bak'])

    f2 = open(seqfile +'.fasta', 'w')
    if f[0] != '>':
        f2.write('>target\n' + f +'\n')
    else:
        x = f.split('\n')
        if len(x[0]) > 6:
            target = x[0][1:5] + x[0][6]
        f2.write('>target\n' + "".join(x[1:]) + '\n')
    f2.close()

    pool = multiprocessing.Pool(n_jobs)
    for i in range(len(names)):
        c_fname = seqfile + '.' + names[i] + '.' + method
        predictionnames[names[i] + method] = c_fname
        pool.apply_async(run_cpred_job, [i, seqfile, method, n_cores, n_jobs])

    pool.close()
    pool.join()

    return predictionnames



def run_pconsc2_dependencies(hhblitsdb, seqfile, n_cores=1):

    if not os.path.exists(seqfile):
        sys.stderr.write('\n' + seqfile + 'does not exist\n')
        sys.exit(0)

    if not os.path.exists(seqfile + '.rsa'):
        sys.stderr.write(str(datetime.now()) + ': running NetSurfP\nThis may take quite a few minutes!\n')
        t = check_output([netsurf, '-i', seqfile, '-a'])
        f = open(seqfile + '.rsa', 'w')
        f.write(t)
        f.close()

    netsurfpredictionname = seqfile + '.rsa'

    if not os.path.exists(seqfile + '.ss2'):
        sys.stderr.write(str(datetime.now()) + ': running Psipred\nThis may take quite a few minutes!\n')
        check_output([psipred, seqfile + '.fasta'])

    sspredictionname = seqfile + '.ss2'

    pssmaliname =  seqfile + '.hhE0.a3m'
    if not os.path.exists(pssmaliname):
        sys.stderr.write(str(datetime.now()) + ' HHblits E0: generating alignment for PSSM feature\nThis may take quite a few minutes!\n ')
        check_output([hhblits, '-all', '-oa3m', pssmaliname, '-e', '1', '-cpu', str(n_cores), '-i', seqfile + '.fasta', '-d', hhblitsdb])

    return netsurfpredictionname, sspredictionname, pssmaliname
