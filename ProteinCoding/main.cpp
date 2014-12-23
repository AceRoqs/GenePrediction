// Prediction of Protein Coding Genes.
// ORF, as used in this program means Open reading frame.
// http://en.wikipedia.org/wiki/Open_reading_frame

#include "PreCompile.h"
#include "fna.h"
#include "ProteinCoding.h"

//---------------------------------------------------------------------------
template<typename Predicate>
void process(
    const std::string& sample_data,
    const std::multimap<size_t, size_t>& ORFs,
    std::vector<size_t>& count,
    std::vector<double>& odds,
    Predicate predicate)
{
    // For all ORFs that satisfy the predicate, calculate probabilities for each k-tuple.
    size_t singles_count = 0;
    size_t doubles_count = 0;
    size_t triples_count = 0;
    size_t quad_count = 0;

    // Build up the probabilities against each ORF that satisfy the predicate.
    auto iter = std::find_if(std::begin(ORFs), std::end(ORFs), predicate);

    while(iter != ORFs.end())
    {
        size_t length = iter->first;

        // Save the combined lengths to use in division later.
        singles_count += length;
        doubles_count += length - 1;
        triples_count += length - 2;
        quad_count    += length - 3;

        // Add up singles.  Get the letter, convert unknown bases to T's,
        // and map that to the proper index.  Increment that index.
        for(size_t ix = 0; ix < length; ++ix)
        {
            count[index_probability(normalize_nucleotide(sample_data[iter->second + ix]))]++;
        }

        // Add up doubles.  Same as with the singles, except that pairs of
        // nucleotides are encoded and that index incremented.
        for(size_t ix = 0; ix < length - 1; ++ix)
        {
            unsigned int term = 0;
            term |= normalize_nucleotide(sample_data[iter->second + ix]) << 8;
            term |= normalize_nucleotide(sample_data[iter->second + ix + 1]);
            count[index_probability(term)]++;
        }

        // Add up triples.
        for(size_t ix = 0; ix < length - 2; ++ix)
        {
            unsigned int term = 0;
            term |= normalize_nucleotide(sample_data[iter->second + ix]) << 16;
            term |= normalize_nucleotide(sample_data[iter->second + ix + 1]) << 8;
            term |= normalize_nucleotide(sample_data[iter->second + ix + 2]);
            count[index_probability(term)]++;
        }

        // Add up quads.
        for(size_t ix = 0; ix < length - 3; ++ix)
        {
            unsigned int term = 0;
            term |= normalize_nucleotide(sample_data[iter->second + ix]) << 24;
            term |= normalize_nucleotide(sample_data[iter->second + ix + 1]) << 16;
            term |= normalize_nucleotide(sample_data[iter->second + ix + 2]) << 8;
            term |= normalize_nucleotide(sample_data[iter->second + ix + 3]);
            count[index_probability(term)]++;
        }

        // Find the next ORF, starting where the iteration left off.
        iter = std::find_if(++iter, ORFs.end(), predicate);
    }

    // All of the sequences have been counted, so get the frequency and save the log.
    for(size_t ix = 0; ix < 4; ++ix)
    {
        odds[ix] = log(static_cast<double>(count[ix]) / singles_count);
    }
    for(size_t ix = 4; ix < 20; ++ix)
    {
        odds[ix] = log(static_cast<double>(count[ix]) / doubles_count);
    }
    for(size_t ix = 20; ix < 84; ++ix)
    {
        odds[ix] = log(static_cast<double>(count[ix])/ triples_count);
    }
    const size_t the_rest = count.size();
    for(size_t ix = 84; ix < the_rest; ++ix)
    {
        odds[ix] = log(static_cast<double>(count[ix]) / quad_count);
    }
}

//---------------------------------------------------------------------------
int main()
{
    std::cout << "Program start." << std::endl;

    // Read in the gbk file and get a map of gene start nucleotides keyed to their stop nucleotide.
    // NC_000909 is M. jannaschii.
    std::cout << "Reading NC_000909.gbk..." << std::endl;
    std::map<size_t, size_t> genes;
    read_gbk("NC_000909.gbk", genes);

    // Read in fna file and extract the sequence information.
    std::cout << "Reading NC_000909.fna..." << std::endl;
    std::string sample_data;
    read_fna("NC_000909.fna", sample_data);

    // Scan through the sequence in one pass, recording the ORFs.
    std::cout << "Recording ORFs..." << std::endl;
    std::multimap<size_t, size_t> ORFs;

    size_t max_ORF = record_ORFs(sample_data, ORFs);

    // Build the log probability tables for each k-tuple.  Use a 3rd order
    // Markov model, which requires space proportional to the sum of the
    // square of each order.  This is more-or-less hard coded for 3rd
    // order, but it can be generalized.
    std::cout << "Building probability tables..." << std::endl;

    // Tables for the log probabilities.
    std::vector<double> odds_ORF(4 + (4 * 4) + (4 * 4 * 4) + (4 * 4 * 4 * 4));
    std::vector<double> odds_background(4 + (4 * 4) + (4 * 4 * 4) + (4 * 4 * 4 * 4));

    // Save the count of each term.  We'll later divide each term by the count
    // associated with that order, in order to get the probabilities for the above
    // tables.
    std::vector<size_t> count_ORF(4 + (4 * 4) + (4 * 4 * 4) + (4 * 4 * 4 * 4));
    std::vector<size_t> count_background(4 + (4 * 4) + (4 * 4 * 4) + (4 * 4 * 4 * 4));

    // For all ORFs larger than 1400 nucleotides, calculate probabilities for each k-tuple.
    process(sample_data, ORFs, count_ORF, odds_ORF, [](const entry_type& entry)
    {
        return 1400 <= entry.first;
    });

    // Do the exact same thing, but do it for sequences less than 50 nucleotides,
    // so the background frequencies are obtained.
    process(sample_data, ORFs, count_background, odds_background, [](const entry_type& entry)
    {
        return 50 >= entry.first;
    });

    // Calculate Markov model for ORFs and print histogram.
    print_histogram(std::cout, sample_data, genes, ORFs, max_ORF, odds_ORF, odds_background);

    std::cout << "Program done." << std::endl;
    return 0;
}

