PconsC2
=======

Improved contact predictions using the recognition of protein like contact patterns.


## Dependencies:

- [NetSurfP 1.1](http://www.cbs.dtu.dk/services/NetSurfP/)
- [PSIPRED v3.5](http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/)
- [Jackhmmer from HMMER v3.0 or higher](http://hmmer.janelia.org/)
- [HHblits from HHsuite v2.0.16](http://toolkit.tuebingen.mpg.de/hhblits)
- [PSICOV v1.11](http://bioinfadmin.cs.ucl.ac.uk/downloads/PSICOV/)
- [plmDCA asymmetric](http://plmdca.csc.kth.se/)
- either MATLAB v8.1 or higher
- or [MATLAB Compiler Runtime (MCR) v8.1](http://www.mathworks.se/products/compiler/mcr/)

MATLAB is needed to run plmDCA. However, if MATLAB is not available you can also use a compiled version of plmDCA. For the compiled version to run you need to provide a path to MCR.


## How to run it:

__Make sure__ all dependencies are working correctly and adjust the paths in `localconfig.py`.

To run PconsC2 use:
```
./run_pconsc.py [-c <n_cores>] [--p_plm <n_plmdca_jobs>] 
             [--p_psi <n_psicov_jobs>] [--pconsc1]
             <hhblits_database> <jackhmmer_database> <sequence_file>
```
- Required:
  - `hhblits_database` and `jackhmmer_database` are paths to the databases used by HHblits and Jackhmmer
  - `sequence_file` is the path to the input protein sequence in FASTA format (only single sequences). 
- Optional:
  - `n_cores` specifies the number of cores to use during computation (default: number of available cores). 
  - `n_plmdca_jobs` specifies the number of plmDCA instances run in parallel (default: min(2, n_cores)).
  - `n_psicov_jobs` specifies the number of PSICOV instances run in parallel (default: min(2, n_cores)).
  - `--pconsc1` flag runs PconsC1 instead of PconsC2

---

