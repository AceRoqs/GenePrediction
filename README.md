This is basic code for doing biological sequence analysis.

* _smith_waterman_ implements local alignment using the
[Smith-Waterman](http://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)
dynamic programming algorithm.
* _viterbi_ is an implementation of the
[Viterbi algorithm](http://en.wikipedia.org/wiki/Viterbi_algorithm)
for solving Hidden Markov Models.
* _protein_coding_ scans genomes for ORFs \(Open reference frames\) using
a probability model.

The code is largely based off of ideas in the book
[Biological Sequence Analysis](http://amzn.to/odfdWC), and uses genome data
that can be obtained from the
[NCBI Genome bank](http://www.ncbi.nlm.nih.gov/genome) and elsewhere.

The code was written in 2006, but I have modernized it to use some nicer
C++ constructs, so it requires a C++0x compiler due to the requirement
of rvalue references and lambdas. The code was tested on Visual C++ 2010,
but should work with GCC or any compiler regardless of platform without
too much trouble.

A word of caution that the code is extremely numeric heavy, and for longer
gene sequences, can take several minutes (in 2010) to compute the results,
even with the dynamic programming optimizations.

Toby Jones \([www.turbohex.com](http://www.turbohex.com), [ace.roqs.net](http://ace.roqs.net)\)

