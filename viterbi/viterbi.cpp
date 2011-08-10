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

#include <cassert>
#include <iostream>
#include <vector>
#include "viterbi.h"
#include "probability.h"

//---------------------------------------------------------------------------
// 300 roll dice example taken from page 57 in Durbin, et al.
// http://amzn.to/odfdWC
const char durbin_dice[] = "315116246446644245311321631164152133625144543631656626566666"
                           "651166453132651245636664631636663162326455236266666625151631"
                           "222555441666566563564324364131513465146353411126414626253356"
                           "366163666466232534413661661163252562462255265252266435353336"
                           "233121625364414432335163243633665562466662632666612355245242";

// Dice rolls taken from page 57 in Durbin, et al., to generate the durbin_dice output.
// (F) Fair die, (L) Loaded die.
const char die_type[] =    "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFLLLLLLLLLLLLLLL"
                           "LLLLLLFFFFFFFFFFFFLLLLLLLLLLLLLLLLFFFLLLLLLLLLLLLLLFFFFFFFFF"
                           "FFFFFFFFLLLLLLLLLLLLLFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFLL"
                           "LLLLLLLLFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
                           "FFFFFFFFFFFFFFFFFFFFFFFFFFFLLLLLLLLLLLLLLLLLLLLLLFFFFFFFFFFF";

//---------------------------------------------------------------------------
double Probability_table::log_prob_at(size_t row, size_t column)
{
    return m_log_prob_matrix[row * m_columns + column];
}

//---------------------------------------------------------------------------
void Probability_table::set_log_prob_at(double log_prob, size_t row, size_t column)
{
    m_log_prob_matrix[row * m_columns + column] = log_prob;
}

//---------------------------------------------------------------------------
// Print out the list of probabilities and log probabilities.
void Probability_table::print_parameters(std::ostream& output_stream)
{
    output_stream << "Emission log prob:\n";
    for(auto prob = m_emission_probabilities.cbegin(); prob != m_emission_probabilities.cend(); ++prob)
    {
        output_stream << log(*prob) << " ("<< *prob << ") \n";
    }

    output_stream << "Initial log prob:\n";
    for(auto prob = m_initial_probabilities.cbegin(); prob != m_initial_probabilities.cend(); ++prob)
    {
        output_stream << log(*prob) << " ("<< *prob << ") \n";
    }

    output_stream << "transition log prob:\n";
    for(auto prob = m_edges.cbegin(); prob != m_edges.cend(); ++prob)
    {
        output_stream << log(*prob) << " ("<< *prob << ") \n";
    }
}

//---------------------------------------------------------------------------
// build_table() does the work for unrolling the HMM to a probability table used for dynamic programming.
void Probability_table::build_table()
{
    // Initialize the first column with the (log) probability of choosing the node,
    // multiplied by (added to) the (log) probability of emitting what the node emitted.
    for(size_t ii = 0; ii < m_rows; ++ii)
    {
        double prob = log_of_sum_of_logs(log(m_initial_probabilities[ii]),
                                         log(m_emission_probabilities[ii * m_emission_count + m_emission_index(m_sample_data[0])]));

        set_log_prob_at(prob, ii, 0);
    }

    // Visit each entry in the table (besides the base cases) and score each.
    // Viterbi needs to be done column-by-column (as opposed to row-by-row).
    for(size_t jj = 1; jj < m_columns; ++jj)
    {
        for(size_t ii = 0; ii < m_rows; ++ii)
        {
            // Take the (log) probability of the previous node, multiplied by (added to) the (log) probability of
            // taking the path/edge of that node to the current node, multiplied by (added to) the (log) probability
            // of emitting what the current node emitted.
            double prob = log_of_sum_of_logs(log_prob_at(0, jj - 1), log(m_edges[0 * m_rows + ii]));
            prob = log_of_sum_of_logs(prob, log(m_emission_probabilities[ii * m_emission_count + m_emission_index(m_sample_data[jj])]));

            // For this node, walk each of the edges that points at this node.
            // Do the probability calculation for that edge, and take the max of all of the calculations.
            for(size_t kk = 1; kk < m_rows; ++kk)
            {
                double new_prob = log_of_sum_of_logs(log_prob_at(kk, jj - 1), log(m_edges[kk * m_rows + ii]));
                new_prob = log_of_sum_of_logs(new_prob, log(m_emission_probabilities[ii * m_emission_count + m_emission_index(m_sample_data[jj])]));

                prob = std::max(new_prob, prob);
            }

            // Save the calculated probability.
            set_log_prob_at(prob, ii, jj);
        }
    }
}

//---------------------------------------------------------------------------
Probability_table::Probability_table(
    std::string&& sample_data,                      // Sample data.
    std::vector<double>&& initial_probabilities,    // Probabilities of transition from begin state.
    std::vector<double>&& edges,                    // Probabilities for each edge (transitions between states).
    std::vector<double>&& emission_probabilities,   // Probabilities for each emission for each model.
    size_t (*emission_index)(char))                 // Function to map emission to an index in the probability vector.
    : m_edges(std::move(edges))
    , m_emission_probabilities(std::move(emission_probabilities))
    , m_emission_count(m_emission_probabilities.size() / initial_probabilities.size()) // number of emission probabilities (i.e. dice=12, 6 emissions * 2 types of dice)
    , m_emission_index(emission_index)
    , m_sample_data(std::move(sample_data))
    , m_initial_probabilities(std::move(initial_probabilities))
    , m_columns(m_sample_data.length())             // Number of samples in the sample data.
    , m_rows(m_initial_probabilities.size())        // Number of Markov models being combined.
{
    m_log_prob_matrix.resize(m_columns * m_rows);
    build_table();
}

//---------------------------------------------------------------------------
// Do trace back of highest probability path (Viterbi traceback).
// Search the final column for the highest score, then save
// the backtrace (m_probable_path) from that entry.
void Probability_table::trace_back_and_save(std::ostream& output_stream)
{
    // Search the final column for the highest log probability.
    // Begin the trace back from that score.
    double max_score = log_prob_at(0, m_columns - 1);
    size_t high_row = 0;
    for(size_t ii = 1; ii < m_rows; ++ii)
    {
        if(log_prob_at(ii, m_columns - 1) > max_score)
        {
            max_score = log_prob_at(ii, m_columns - 1);
            high_row = ii;
        }
    }

    output_stream << "Viterbi path log probability: " << max_score << "\n";

    // The probable path is a list of the followed rows.
    m_probable_path.resize(m_columns);
    m_probable_path[m_columns - 1] = high_row;

    // Walk the columns in reverse order for the traceback.  Calculate the previous nodes'
    // (log) probabilities, and follow the path with the max score.
    for(size_t jj = m_columns - 1; jj > 0; --jj)
    {
        size_t new_high_row = 0;  // Assume initial row has the new highest probability.

        double prob = log_of_sum_of_logs(log_prob_at(0, jj - 1), log(m_edges[0 * m_rows + high_row]));
        prob = log_of_sum_of_logs(prob, log(m_emission_probabilities[high_row * m_emission_count + m_emission_index(m_sample_data[jj])]));

        for(size_t kk = 1; kk < m_rows; ++kk)
        {
            double new_prob = log_of_sum_of_logs(log_prob_at(kk, jj - 1), log(m_edges[kk * m_rows + high_row]));
            new_prob = log_of_sum_of_logs(new_prob, log(m_emission_probabilities[high_row * m_emission_count + m_emission_index(m_sample_data[jj])]));

            // If this row scored a higher probability than the previous max,
            // take this row as the new max.
            if(new_prob > prob)
            {
                new_high_row = kk;
                prob = new_prob;
            }
        }

        // After walking all the rows in this column, save the highest scoring one and use
        // that as the basis for the next column's score.
        high_row = new_high_row;
        m_probable_path[jj - 1] = high_row;
    }
}

//---------------------------------------------------------------------------
// This function makes the assumption of two states (unlike the rest of the HMM code).
void Probability_table::print_found_sequences(std::ostream& output_stream, size_t max_hits, size_t min_nucleotide_count)
{
    print_parameters(output_stream);

    output_stream << "Printing hits:\n";

    bool in_sequence = false;
    size_t start_index = 0;
    size_t hit_count = 0;

    // Walk the columns.  When start of a hit is found, save the index, and
    // when hits stop, then print out the sequence from the saved index
    // to the current index.
    for(size_t jj = 0; jj < m_columns; ++jj)
    {
        if(0 == m_probable_path[jj])
        {
            if(in_sequence)
            {
                in_sequence = false;

                // Only print significant hits.
                if(jj - start_index < min_nucleotide_count)
                {
                    continue;
                }

                // Print sequence from start_index to jj - 1.
                output_stream << "Hit " << ++hit_count << ": location: " << start_index << ".." << jj - 1 << " length: " << jj - start_index << "\n";
                while(start_index < jj)
                {
                    output_stream << m_sample_data[start_index];
                    ++start_index;
                }

                output_stream << "\n\n";

                // Exit early once the max number of hits has been reached.
                if((0 != max_hits) && (hit_count == max_hits))
                {
                    break;
                }
            }
        }
        else
        {
            // Save the start index if a sequence is not in progress.
            if(!in_sequence)
            {
                start_index = jj;
                in_sequence = true;
            }
        }
    }
}

//---------------------------------------------------------------------------
// This behaves the same as print_found_sequences(), but just returns a count
// instead of printing the sequences.
size_t Probability_table::count_hits()
{
    bool in_sequence = false;
    size_t start_index = 0;
    size_t hit_count = 0;

    for(size_t jj = 0; jj < m_columns; ++jj)
    {
        if(0 == m_probable_path[jj])
        {
            // Count the if a sequence was in progress.
            if(in_sequence)
            {
                in_sequence = false;
                ++hit_count;

                while(start_index < jj)
                {
                    ++start_index;
                }
            }
        }
        else
        {
            // Save the start index if a sequence is not in progress.
            if(!in_sequence)
            {
                start_index = jj;
                in_sequence = true;
            }
        }
    }

    return hit_count;
}

//---------------------------------------------------------------------------
// Implement Viterbi training across edges (Ak,l) and pEmissionProbabilities (Ek(b))
void Probability_table::train_and_print(std::ostream& output_stream)
{
    // Calculate new estimate of transition probabilities.
    std::vector<size_t> edges(m_rows);

    // TODO: A O(n) algorithm would be better here.
    for(size_t ii = 0; ii < m_rows; ++ii)
    {
        edges.assign(m_rows, 0);

        // Count the number of each transition taken.
        // Walk m_columns - 1 because there are only edges between the nodes.
        size_t total_edges = 0;
        for(size_t jj = 0; jj < m_columns - 1; ++jj)
        {
            // If this node was in the trace, count its attached edge.
            // Also increase the count of the specific edge that was followed.
            if(m_probable_path[jj] == ii)
            {
                edges[m_probable_path[jj + 1]] += 1;
                ++total_edges;
            }
        }

        // Calculate new estimate of this edge probability.
        for(size_t jj = 0; jj < m_rows; ++jj)
        {
            m_edges[ii * m_rows + jj] = edges[jj] / static_cast<double>(total_edges);
        }
    }

    // Calculate new estimate of emission probabilities.
    std::vector<size_t> emissions(m_emission_count);

    // TODO: A O(n) algorithm would be better here.
    for(size_t ii = 0; ii < m_rows; ++ii)
    {
        emissions.assign(m_emission_count, 0);

        // Count the number of emissions, keeping a separate count for each row.
        size_t total_emissions = 0;
        for(size_t jj = 0; jj < m_columns; ++jj)
        {
            // If this node was in the trace, then count its emission towards the total,
            // and also up the count of the specific letter that was emitted.
            // Keep a separate total for each row.
            if(m_probable_path[jj] == ii)
            {
                emissions[m_emission_index(m_sample_data[jj])] += 1;
                ++total_emissions;
            }
        }

        // Calculate estimate of new probability of emitting this character.
        for(size_t jj = 0; jj < m_emission_count; ++jj)
        {
            m_emission_probabilities[ii * m_emission_count + jj] = emissions[jj] / static_cast<double>(total_emissions);
        }
    }

    // Print out the newly estimated parameters.
    build_table();
    trace_back_and_save(output_stream);
    print_parameters(output_stream);
    output_stream << "Hits: " << count_hits() << "\n" << std::endl;
}

//---------------------------------------------------------------------------
#ifdef _DEBUG
void Probability_table::print_dice_rolls(std::ostream& output_stream)
{
    // Test code for dice example.
    const size_t display_length = 60;
    const size_t total_rows = (m_columns + display_length - 1) / display_length;

    // For ease of testing, ensure that the sample data is a multiple of the display length.
    assert((m_columns % display_length) == 0);

    for(size_t blocks = 0; blocks < total_rows; ++blocks)
    {
        output_stream << "Rolls:   ";
        for(size_t length = 0; length < m_columns / total_rows; ++length)
        {
            output_stream << m_sample_data[display_length * blocks + length];
        }
        output_stream << "\nDie:     ";
        for(size_t length = 0; length < m_columns / total_rows; ++length)
        {
            output_stream << die_type[display_length * blocks + length];
        }
        output_stream << "\nViterbi: ";
        for(size_t length = 0; length < m_columns / total_rows; ++length)
        {
            if(m_probable_path[display_length * blocks + length])
            {
                output_stream << "L";
            }
            else
            {
                output_stream << "F";
            }
        }
        output_stream << "\n" << std::endl;
    }
}
#endif

