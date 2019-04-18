# RotatedKleeMintyProblem
A linear constrained optimization benchmark for probabilistic search algorithms

Considering the domain of constrained optimization, the number of available benchmark environments bears no relation to the amount of distinct problem features. The rotated Klee-Minty problem represents the proposal of a scalable linear constrained optimization problem which is suitable for benchmarking Evolutionary Algorithms.

This repository provides a Matlab implementation of the benchmarking environment introduced in the paper:
"A Linear Constrained Optimization Benchmark For Probabilistic Search Algorithms: The Rotated Klee-Minty Problem"
by M. Hellwig and H.-G. Beyer, preprint available at \url{avialable soon}.

The reporting of implementation issues as well as suggestions for improvement (or modification) are welcome. 
Please contact 
hemi_at_fhv_dot_at

# Benchmarking
The Matlab code has been extended for the GECCO 2019
workshop paper "Comparison of Contemporary Evolutionary Algorithms on
the Rotated Klee-Minty Problem" by M. Hellwig, P. Spettel, and H.-G. Beyer.

Instructions for running the code accompanying that paper
follow below. Clone the repository and checkout the relevant branch
using the commands

    $ cd /to/some/path
    $ git clone https://github.com/patsp/RotatedKleeMintyProblem.git
    $ cd RotatedKleeMintyProblem
    $ git checkout ea_comparison

Adjust the general settings in `experimentRotatedKleeMintyProblem.m`
to your needs (e.g., the dimensions to run). Then, change
`input.strategy` in `experimentRotatedKleeMintyProblem.m`
to 'RandS', 'runlcCMSAESOnRotatedKleeMintyProblem', and
'runepsMAgESOnRotatedKleeMintyProblem' one after the other.
For every setting of `input.strategy`,
issue the following command to perform the experiments:

    $ matlab -nodesktop -r 'experimentRotatedKleeMintyProblem;exit;'

Note that these calls create folders with the experimental data
in Matlab `.mat` files.

After that, run

    $ matlab -nodesktop -r 'pprocRotatedKleeMintyProblem;exit;'

for the post-processing. It reads the experimental data and
creates folders containing the generated plots.

