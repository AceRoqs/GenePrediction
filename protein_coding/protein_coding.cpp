/*
Copyright (C) 2006-2011 by Toby Jones.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

// Prediction of Protein Coding Genes.

#include <cassert>
#include <fstream>
#include <algorithm>    // for find_if
#include <map>          // for map/multimap
#include <vector>
#include "protein_coding.h"
#include "probability.h"

//---------------------------------------------------------------------------
// Index into array of probabilities for each of the k-tuples.
size_t index_probability(unsigned int nucleotide)
{
    assert(nucleotide != 0);

    if('G' == nucleotide)
    {
        return 0;
    }
    if('C' == nucleotide)
    {
        return 1;
    }
    if('A' == nucleotide)
    {
        return 2;
    }
    if('T' == nucleotide)
    {
        return 3;
    }

    // Fancy math to return 'G' - 'T'     => indices 0-3
    //                      'GG' - 'TT'   => indices 4-19
    //                      'GGG' - 'TTT' => indices 20-83
    //                      etc...
    return ((index_probability(nucleotide >> 8) + 1) * 4) + index_probability(nucleotide & 0xff);
}

//---------------------------------------------------------------------------
// Substitute any character non G, C or A with T.
char normalize_nucleotide(char nucleotide)
{
    if((nucleotide == 'G') ||
       (nucleotide == 'C') ||
       (nucleotide == 'A'))
    {
        return nucleotide;
    }

    return 'T';
}

//---------------------------------------------------------------------------
// Take a codon which has been packed into an unsigned int, and check
// if it's a stop codon (TAA, TAG, or TGA).
bool is_stop_codon(unsigned int codon)
{
    bool stop_codon = false;
    stop_codon |= (codon == (('T' << 16) | ('A' << 8) | 'A'));
    stop_codon |= (codon == (('T' << 16) | ('A' << 8) | 'G'));
    stop_codon |= (codon == (('T' << 16) | ('G' << 8) | 'A'));

    return stop_codon;
}

//---------------------------------------------------------------------------
// Parse the gbk file looking for coding sequences.
// gbk files are GenBank files that describe genomes.
// NCBI has genomes available for download.
// NOTE: The file is assumed to come from a trusted source and be well-formed.
void read_gbk(_In_ const char* filename, std::map<size_t, size_t>& genes)
{
    std::ifstream input_file(filename);
    char line[256];

    while(input_file.good())
    {
        // Ignore everything until the CDS tag comes along.
        input_file >> line;
        if(strcmp(line, "CDS") == 0)
        {
            // Differentiate between complement and non-complement tags.
            input_file >> line;
            if(strncmp(line, "complement", sizeof "complement" - 1) != 0)
            {
                // This is a non-complement tag.  Replace all '.' with nulls,
                // and then extract the two values.
                char* start_nucleotide = line;
                char* stop_nucleotide = "";
                char* scan_character = line;
                while('\0' != scan_character[0])
                {
                    if('.' == scan_character[0])
                    {
                        scan_character[0] = '\0';
                        stop_nucleotide = scan_character + 1;
                    }

                    ++scan_character;
                }

                // Get the start and stop nucleotides, and make them 0-based instead of 1-based.
                int start = atoi(start_nucleotide) - 1;
                int stop = atoi(stop_nucleotide) - 1;

                // Insert into the gene map, keyed to the stop position for searching.
                genes.insert(std::pair<size_t, size_t>(stop, start));
            }
        }
    }

    input_file.close();
}

//---------------------------------------------------------------------------
// Do a one pass scan through the sequence data, recording the ORFs.
// Return the length of the longest ORF.
size_t record_ORFs(const std::string& sample_data, std::multimap<size_t, size_t>& ORFs)
{
    size_t max_ORF = 0;

    // Encode the first two characters of the sequence data to prime the loop.
    unsigned int codon = (sample_data[0] << 8) | sample_data[1];

    // Save the start position of each codon reading frame.
    size_t start_nucleotide[3] = { 0, 1, 2 };

    const size_t nucleotide_count = sample_data.size();
    for(size_t index = 2; index < nucleotide_count; ++index)
    {
        // Drop the first character of the codon, and OR in the next nucleotide.
        codon = codon << 8;
        codon |= normalize_nucleotide(sample_data[index]);
        codon &= 0x00FFFFFF;

        // If this is a stop codon, record the reading frame.
        if(is_stop_codon(codon))
        {
            size_t start = start_nucleotide[(index - 2) % 3];
            size_t length = index - start_nucleotide[(index - 2) % 3] + 1;

            // Assert that it's a multiple of the codon length.
            // Also record the max length.
            assert((length % 3) == 0);
            max_ORF = std::max(max_ORF, length);

            // Save the start nucleotide keyed against the length.
            ORFs.insert(std::pair<size_t, size_t>(length, start));

            // Set a new start nucleotide for this reading frame.
            start_nucleotide[(index - 2) % 3] = index + 1;
        }
    }

    return max_ORF;
}

//---------------------------------------------------------------------------
// Calculate Markov model for ORF and print a histogram.
void print_histogram(
    std::ostream& output_stream,
    const std::string& sample_data,
    const std::map<size_t, size_t>& genes,
    const std::multimap<size_t, size_t>& ORFs,
    size_t max_ORF,
    const std::vector<double>& odds_ORF,
    const std::vector<double>& odds_background)
{
    output_stream << "Printing matches..." << std::endl;

    // Print out stats for each ORF length.
    for(size_t index = 3; index < max_ORF; index += 3)
    {
        size_t gene_count = 0;
        size_t ORF_count = 0;
        double avg_log_odds = 0.0;
        size_t positive = 0;
        size_t pos_and_real = 0;

        auto iter = std::find_if(std::begin(ORFs), std::end(ORFs), [index](const entry_type& entry)
        {
            return index == entry.first;
        });

        while(iter != ORFs.end())
        {
            // Accumulate log probabilities for the Markov model.
            double P_prob = 0;
            double Q_prob = 0;
            unsigned int term = 0;

            for(size_t ii = 0; ii < iter->first; ++ii)
            {
                term = term << 8;
                term &= 0xFFFFFFFF;
                term |= normalize_nucleotide(sample_data[iter->second + ii]);

                size_t index_prob = index_probability(term);

                P_prob = log_of_sum_of_logs(P_prob, odds_ORF[index_prob]);
                Q_prob = log_of_sum_of_logs(Q_prob, odds_background[index_prob]);
            }

            // Accumulate log odds for average later.
            double log_odds = log(P_prob / Q_prob);
            avg_log_odds += log_odds;
            ++ORF_count;

            // Record if the log odds are positive.
            if(log_odds > 0)
            {
                ++positive;
            }

            // Is there a gene from gbk file that has the same stop codon?
            // Since the genes std::map (from the gbk file) is considered correct,
            // this validates that the genes that were predicted are valid.
            // (second = start, first = length)
            size_t stop = iter->second + iter->first - 1;
            auto genes_iterator = genes.find(stop);
            if(genes_iterator != genes.end())
            {
                // Record the match, and also mark if there are positive log_odds for this match.
                ++gene_count;

                if(log_odds > 0)
                {
                    ++pos_and_real;
                }
            }

            // Continue with next ORF, starting the search where the iteration left off.
            iter = std::find_if(++iter, ORFs.end(), [index](const entry_type& entry)
            {
                return index == entry.first;
            });
        }

        // If at least one gene for this length was matched, print out the stats.
        if(gene_count > 0)
        {
            avg_log_odds = avg_log_odds / static_cast<double>(ORF_count);

            output_stream << "ORF Len: " << index
                          << ", gene/ORF: " << gene_count << "/" << ORF_count
                          << ", real/pos: " << pos_and_real << "/" << positive
                          << ", avg log odds: " << avg_log_odds
                          << std::endl;
        }
    }
}

