#include "PreCompile.h"
#include "SmithWaterman.h"
#include "ScorePolicy.h"

//---------------------------------------------------------------------------
int Alignment_table::score_at(size_t row, size_t column) const
{
    return m_score_table[row * m_columns + column];
}

//---------------------------------------------------------------------------
void Alignment_table::set_score_at(int score, size_t row, size_t column)
{
    m_score_table[row * m_columns + column] = score;
}

//---------------------------------------------------------------------------
// Helper function for printing trace backs.
void Alignment_table::print_trace_back(std::ostream& output_stream, size_t row, size_t column, std::vector<residue_pair> optimal_alignment) const
{
    bool print_when_done = true;

    // Don't continue the back trace if the base case is hit.
    if((row > 0) && (column > 0))
    {
        // Calculate potential score paths.  Calculate the scores as done originally, but only
        // follow scores that match what was actually used.
        int current_score =  score_at(row, column);
        int diagonal_score = score_at(row - 1, column - 1) + m_score_policy(m_sequence2[row - 1], m_sequence1[column - 1]);
        int above_score =    score_at(row - 1, column)     + m_score_policy(m_sequence2[row - 1], gap_character);
        int left_score =     score_at(row, column - 1)     + m_score_policy(gap_character,        m_sequence1[column - 1]);

        // Matching scores are not expected to be less than 0, but don't follow a trace of 0's.
        //
        // The way this works is that a residue pair is pushed into the optimial_alignment list,
        // and a _copy_ of the list is sent to each recursion.  This way, each trace can be followed independently.
        //
        // TODO: There is an inefficiency here in that the list is copied for each pair instead of
        // for each branch taken.  This is fine since the running cost is mostly in building the table,
        // rather than the trace back.  This could be an area for improvement.
        if((above_score == current_score) && (above_score > 0))
        {
            residue_pair pair;
            pair.residue1 = gap_character;
            pair.residue2 = m_sequence2[row - 1];

            optimal_alignment.push_back(pair);

            print_trace_back(output_stream, row - 1, column, optimal_alignment);

            optimal_alignment.pop_back();

            print_when_done = false;
        }

        if((left_score == current_score) && (left_score > 0))
        {
            residue_pair pair;
            pair.residue1 = m_sequence1[column - 1];
            pair.residue2 = gap_character;

            optimal_alignment.push_back(pair);

            print_trace_back(output_stream, row, column - 1, optimal_alignment);

            optimal_alignment.pop_back();

            print_when_done = false;
        }

        if((diagonal_score == current_score) && (diagonal_score > 0))
        {
            residue_pair pair;
            pair.residue1 = m_sequence1[column - 1];
            pair.residue2 = m_sequence2[row - 1];

            optimal_alignment.push_back(pair);

            print_trace_back(output_stream, row - 1, column - 1, optimal_alignment);

            optimal_alignment.pop_back();

            print_when_done = false;
        }
    }

    // If no work was done, then the end of the local alignment has been reached.  Print it out.
    if(print_when_done)
    {
        print_alignment(output_stream, optimal_alignment);
    }
}

//---------------------------------------------------------------------------
void Alignment_table::print_alignment(std::ostream& output_stream, const std::vector<residue_pair>& alignment) const
{
    // The vector contains a reverse list of residue pair alignments.  Walk it backwards
    // and print each character.
    for(auto residue = alignment.rbegin(); residue != alignment.rend(); ++residue)
    {
        output_stream << residue->residue1;
    }

    output_stream << "\n";

    for(auto residue = alignment.rbegin(); residue != alignment.rend(); ++residue)
    {
        output_stream << residue->residue2;
    }

    output_stream << "\n";
}

//---------------------------------------------------------------------------
// Method to permute a sequence.
// TODO: The choice of rand() may have an effect here.  Something to test would be
// a mersenne twister rand() which has equidistribution properties, and see if the p-value
// is significantly changed.
void permute_sequence(std::string& sequence)
{
    for(size_t ix = sequence.size(); ix > 0; --ix)
    {
        std::swap(sequence[ix - 1], sequence[rand() % ix]);
    }
}

//---------------------------------------------------------------------------
// Constructor that takes a pair of sequences to align.
Alignment_table::Alignment_table(const std::string& sequence1, const std::string& sequence2, int (*score_policy)(char char1, char char2))
    : m_columns(sequence1.length() + 1)
    , m_rows(sequence2.length() + 1)
    , m_sequence1(sequence1)
    , m_sequence2(sequence2)
    , m_score_policy(score_policy)
{
    // Create and init entries to 0 score.
    m_score_table.resize(m_columns * m_rows);

    // Visit each entry in the table (besides the base cases) and score each.
    for(size_t row = 1; row < m_rows; ++row)
    {
        for(size_t column = 1; column < m_columns; ++column)
        {
            int diagonal_score = score_at(row - 1, column - 1) + m_score_policy(m_sequence2[row - 1], m_sequence1[column - 1]);
            int above_score =    score_at(row - 1, column)     + m_score_policy(m_sequence2[row - 1], gap_character);
            int left_score =     score_at(row, column - 1)     + m_score_policy(gap_character,        m_sequence1[column - 1]);

            // Take the max score of 0 and the three potential scores and save it.
            int score = std::max(0, diagonal_score);
            score = std::max(score, left_score);
            score = std::max(score, above_score);

            m_max_score = std::max(m_max_score, score);
            set_score_at(score, row, column);
        }
    }
}

//---------------------------------------------------------------------------
// Print all of the trace backs.  Search the score table
// for scores that match the maximum and then call a helper
// method to print all of the traces from that entry.
void Alignment_table::print_trace_back(std::ostream& output_stream) const
{
    output_stream << "Optimal score: " << m_max_score << "\nTrace back sequences:\n";

    std::vector<residue_pair> optimal_alignment;

    for(size_t row = 0; row < m_rows; ++row)
    {
        for(size_t column = 0; column < m_columns; ++column)
        {
            if(score_at(row, column) == m_max_score)
            {
                print_trace_back(output_stream, row, column, optimal_alignment);
            }
        }
    }
}

//---------------------------------------------------------------------------
// Print the score table by walking the matrix score by score and printing each number.
void Alignment_table::print_table(std::ostream& output_stream) const
{
    for(size_t row = 0; row < m_rows; ++row)
    {
        for(size_t column = 0; column < m_columns; ++column)
        {
            output_stream << std::setw(4) << score_at(row, column);
        }

        output_stream << "\n";
    }

    output_stream << "\n";
}

//---------------------------------------------------------------------------
// Calculate p-value for current sequence pair.
// A sequence is chosen, permuted, and scored against the other sequence.
// The p-value is the score k/N:
//      k=number of scores higher than original alignment
//      N=number of permutations total
//
// There are multiple alignments resulting from multiple trace backs, but
// only the score the mostly recently cached alignment.
void Alignment_table::calc_pvalue(std::ostream& output_stream, unsigned int num_permutations) const
{
    unsigned int num_better_scores = 0;

    std::string permuted_sequence(m_sequence2);
    for(unsigned int ix = 0; ix < num_permutations; ++ix)
    {
        permute_sequence(permuted_sequence);
        Alignment_table test_table(m_sequence1, permuted_sequence, m_score_policy);

        if(test_table.m_max_score > m_max_score)
        {
            ++num_better_scores;
        }
    }

    output_stream << "p-value: " << static_cast<float>(num_better_scores) / num_permutations
                  << " (" << num_better_scores << " / " << num_permutations << ")\n\n";
}

