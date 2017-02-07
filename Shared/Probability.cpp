#include "PreCompile.h"
#include "Probability.h"

// Based on Numerically Stable Hidden Markov Model Implementation
// by Tobias Mann
// http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf
//
// Given two probabilities x and y, represented by their logs lx, ly,
// return the log of their sum log(x+y) = log(exp(lx) + exp(ly)).
//
// Assume log(0) is represented by NaN.
//
// The "lx > ly" trick is some protection from underflow:
//   log(a+b) = log(a(1+b/a)) = log(a)+log(1+b/a),
// which will be most accurate when b/a < 1.
double log_of_sum_of_logs(double lx, double ly)
{
    if(_isnan(lx))
    {
        return ly;
    }

    if(_isnan(ly))
    {
        return lx;
    }

    if(lx > ly)
    {
        return lx + log(1 + exp(ly - lx));
    }
    else
    {
        return ly + log(1 + exp(lx - ly));
    }
}

