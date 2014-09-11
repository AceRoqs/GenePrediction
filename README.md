This is basic code for doing biological sequence analysis.

* SmithWaterman implements local alignment using the
[Smith-Waterman](http://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)
dynamic programming algorithm.
* Viterbi is an implementation of the
[Viterbi algorithm](http://en.wikipedia.org/wiki/Viterbi_algorithm)
for solving Hidden Markov Models.
* ProteinCoding scans genomes for ORFs \(open reference frames\) using
a probability model.

The code is largely based off of ideas in the book
[Biological Sequence Analysis](http://amzn.to/odfdWC), and uses genome data
that can be obtained from the
[NCBI Genome bank](http://www.ncbi.nlm.nih.gov/genome) and elsewhere.

The code was written in 2006, but I have modernized it to use some nicer
C++ constructs, so it requires a C++11 compiler due to the requirement
of rvalue references and lambdas. The code was tested on Visual C++ 2013,
but should work with GCC or any compiler regardless of platform without
too much trouble.

A word of caution that the code is extremely numeric heavy, and for longer
gene sequences, can take several seconds to several minutes (in 2013) to
compute the results, even with the dynamic programming optimizations.

Some of these algorithms can be useful outside the computational biology
space. For example, Smith-Waterman can be useful in doing text correction,
and Viterbi can be used to do signal analysis for speech-to-text.

Toby Jones \([www.turbohex.com](http://www.turbohex.com), [ace.roqs.net](http://ace.roqs.net)\)

