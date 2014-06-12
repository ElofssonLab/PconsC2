#!/usr/bin/env python
from localconfig import *
from datetime import datetime
import sys, subprocess, os


def check_output(command):
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

    names = ['E4', 'E0', 'E10', 'E40']
    cutoffs = ['1e-4', '1', '1e-10', '1e-40']

    failed = []

    for i in range(4):
        
        exists_jh = os.path.exists(seqfile + '.jh' + names[i] + '.a3m')
        exists_jh_psicov = os.path.exists(seqfile + '.jh' + names[i] + '.psicov')
        exists_jh_plmdca = os.path.exists(seqfile + '.jh' + names[i] + '.plmdca')
        exists_hh = os.path.exists(seqfile + '.hh' + names[i] + '.a3m')
        exists_hh_psicov = os.path.exists(seqfile + '.hh' + names[i] + '.psicov')
        exists_hh_plmdca = os.path.exists(seqfile + '.hh' + names[i] + '.plmdca')

        # only create alignment file if at least one of the contact maps is missing
        if not exists_jh and (not exists_jh_psicov or not exists_jh_plmdca):
            sys.stderr.write(str(datetime.now()) + ' jackhmmer ' + names[i] + ': generating alignment\nThis may take quite a few minutes!\n ')
            t = check_output([jackhmmer, '--cpu', str(n_cores), '-N', '5', '-E', cutoffs[i], '-A', seqfile +'.jh' + names[i] + '.ali', seqfile + '.fasta', jackhmmerdb])
            check_output([reformat, 'sto', 'a3m', seqfile + '.jh' + names[i] + '.ali', seqfile + '.jh' + names[i] + '.a3m'])
            check_output(['rm', seqfile + '.jh' + names[i] + '.ali'])

        # only create alignment file if at least one of the contact maps is missing
        if not exists_hh and (not exists_hh_psicov or not exists_hh_plmdca):
            sys.stderr.write(str(datetime.now()) + ' HHblits' + names[i] + ': generating alignment\nThis may take quite a few minutes!\n ')
            t = check_output([hhblits, '-all', '-oa3m', seqfile + '.hh' + names[i] + '.a3m', '-e', cutoffs[i], '-cpu', str(n_cores), '-i', seqfile + '.fasta', '-d', hhblitsdb])
            #check_output([reformat, 'a3m', 'fas', seqfile + '.hh' + names[i] + '.a3m', seqfile + '.hh' + names[i] + '.fas'])
        


def run_psicov(seqfile, jhpredictionnames, hhpredictionnames, n_cores=1):

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

    names = ['E4', 'E0', 'E10', 'E40']
    cutoffs = ['1e-4', '1', '1e-10', '1e-40']

    failed = []

    for i in range(4):
        
        exists_jh_psicov = os.path.exists(seqfile + '.jh' + names[i] + '.psicov')
        exists_hh_psicov = os.path.exists(seqfile + '.hh' + names[i] + '.psicov')

        if not exists_jh_psicov:
            #t = check_output([trim, seqfile + '.jh' + names[i] + '.fas'])
            t = check_output([trim2jones, seqfile + '.jh' + names[i] + '.a3m'])
            f = open(seqfile + '.jh' + names[i] + '.jones', 'w')
            f.write(t)
            f.close()

            t = ''
            sys.stderr.write(str(datetime.now()) + ' jackhmmer ' + names[i] + ': running PSICOV\nThis may take more than an hour.\n')
            try:
                # Joel @ NSC: Added -o flag, in case the psicov binary has not
                # been compiled with MINEFSEQS=0.
                t = check_output([psicov, '-o', seqfile + '.jh' + names[i] + '.jones'])
            except:
                t = ''
            f = open(seqfile + '.jh' + names[i] + '.psicov', 'w')
            f.write(t)
            f.close()

        jhpredictionnames[names[i] + 'psicov'] = seqfile + '.jh' + names[i] + '.psicov'
        
        if not exists_hh_psicov:
            #t = check_output([trim, seqfile + '.hh' + names[i] + '.fas'])
            t = check_output([trim2jones, seqfile + '.hh' + names[i] + '.a3m'])
            f = open(seqfile + '.hh' + names[i] + '.jones', 'w')
            f.write(t)
            f.close()
            
            sys.stderr.write(str(datetime.now()) + ' HHblits ' + names[i] + ': running PSICOV\nThis may take more than an hour.\n')
            t = ''
            try:
                # Joel @ NSC: Added -o flag, in case the psicov binary has not
                # been compiled with MINEFSEQS=0.
                t = check_output([psicov, '-o', seqfile + '.hh' + names[i] + '.jones'])
            except:
                t = ''
            f = open(seqfile + '.hh' + names[i] + '.psicov', 'w')
            f.write(t)
            f.close()

        hhpredictionnames[names[i] + 'psicov'] = seqfile + '.hh' + names[i] + '.psicov'

    return jhpredictionnames, hhpredictionnames



def run_plmdca(seqfile, jhpredictionnames, hhpredictionnames, n_cores=1):

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

    names = ['E4', 'E0', 'E10', 'E40']
    cutoffs = ['1e-4', '1', '1e-10', '1e-40']

    failed = []

    for i in range(4):
        
        exists_jh_plmdca = os.path.exists(seqfile + '.jh' + names[i] + '.plmdca')
        exists_hh_plmdca = os.path.exists(seqfile + '.hh' + names[i] + '.plmdca')

        if not exists_jh_plmdca:
            t = check_output([trim2trimmed, seqfile + '.jh' + names[i] + '.a3m'])
            f = open(seqfile + '.jh' + names[i] + '.trimmed', 'w')
            f.write(t)
            f.close()

            sys.stderr.write(str(datetime.now()) + ' jackhmmer ' + names[i] + ': running plmDCA\nThis may take more than an hour.\n')
            if plmdca:
                #t = check_output([plmdca, matlabdir, seqfile + '.jh' + names[i] + ".trimmed", seqfile + '.jh' + names[i] + ".plmdca", "0.01", "0.01", "0.1", str(n_cores)])
                t = check_output([plmdca, seqfile + '.jh' + names[i] + ".trimmed", seqfile + '.jh' + names[i] + ".plmdca", "0.01", "0.01", "0.1", str(n_cores)])
            else:
                t = check_output([matlab, '-nodesktop', '-nosplash', '-r', "path(path, '" + plmdcapath + "'); path(path, '" + plmdcapath + "/functions'); path(path, '" + plmdcapath + "/3rd_party_code/minFunc/'); plmDCA_symmetric ( '" + seqfile + '.jh' + names[i] + ".trimmed', '" + seqfile + '.jh' + names[i] + ".plmdca', 0.01, 0.01, 0.1, " + str(n_cores) + "); exit"])

        jhpredictionnames[names[i] + 'plmdca'] = seqfile + '.jh' + names[i] + '.plmdca'

        if not exists_hh_plmdca:
            #t = check_output([trim2, seqfile + '.hh' + names[i] + '.fas'])
            t = check_output([trim2trimmed, seqfile + '.hh' + names[i] + '.a3m'])
            f = open(seqfile + '.hh' + names[i] + '.trimmed', 'w')
            f.write(t)
            f.close()

            sys.stderr.write(str(datetime.now()) + ' HHblits ' + names[i] + ': running plmDCA\nThis may take more than an hour.\n')
            if plmdca:
                #t = check_output([plmdca, matlabdir, seqfile + '.hh' + names[i] + ".trimmed", seqfile + '.hh' + names[i] + ".plmdca", "0.01", "0.01", "0.1", str(n_cores)])
                t = check_output([plmdca, seqfile + '.hh' + names[i] + ".trimmed", seqfile + '.hh' + names[i] + ".plmdca", "0.01", "0.01", "0.1", str(n_cores)])
            else:
                t = check_output([matlab, '-nodesktop', '-nosplash', '-r', "path(path, '" + plmdcapath + "'); path(path, '" + plmdcapath + "/functions'); path(path, '" + plmdcapath + "/3rd_party_code/minFunc/'); plmDCA_symmetric ( '" + seqfile + '.hh' + names[i] + ".trimmed', '" + seqfile + '.hh' + names[i] + ".plmdca', 0.01, 0.01, 0.1, " + str(n_cores) + "); exit"])

        hhpredictionnames[names[i] + 'plmdca'] = seqfile + '.hh' + names[i] + '.plmdca'

    return jhpredictionnames, hhpredictionnames



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

    if not os.path.exists(seqfile + '.horiz'):
        sys.stderr.write(str(datetime.now()) + ': running Psipred\nThis may take quite a few minutes!\n')
        check_output([psipred, seqfile + '.fasta'])

    sspredictionname = seqfile + '.horiz'

    pssmaliname =  seqfile + '.hhE0.a3m'
    if not os.path.exists(pssmaliname):
        sys.stderr.write(str(datetime.now()) + ' HHblits E0: generating alignment for PSSM feature\nThis may take quite a few minutes!\n ')
        t = check_output([hhblits, '-all', '-oa3m', seqfile + '.hhE0.a3m', '-e', '1', '-cpu', str(n_cores), '-i', seqfile + '.fasta', '-d', hhblitsdb])

    return netsurfpredictionname, sspredictionname
