#ifndef VITERBI_H
#define VITERBI_H

//---------------------------------------------------------------------------
// Definition of a probability table for dynamic programming.
class Probability_table
{
private:
    // Because multiplication of successive probabilities produces extremely
    // small numbers which are numerically unstable, use log probabilities
    // instead, which can be added without the numerical issues.
    std::vector<double> m_log_prob_matrix;              // 2D matrix of log probabilities.
    std::vector<double> m_edges;                        // [m_rows x m_rows] matrix of edge probabilities.
    std::vector<double> m_emission_probabilities;       // Probability of each emission.
    std::vector<size_t> m_probable_path;                // List of rows indicating the probable path.

    size_t m_emission_count;                            // Number of potential emissions.

    size_t (*m_emission_index)(char);                   // Function to map emission to an index in the probability vector.

    const std::string m_sample_data;                    // String of sample data to model (on j axis).
    const std::vector<double> m_initial_probabilities;  // Vector of initial probabilities (number of Markov models being combined).
    const size_t m_columns;                             // Width of matrix (m_sample_data.length()).
    const size_t m_rows;                                // Height of matrix (m_initial_probabilities.size()).

    // Not implemented to prevent accidental copying.
    Probability_table(const Probability_table&);
    Probability_table& operator=(const Probability_table&);

protected:
    double log_prob_at(size_t row, size_t column);
    void set_log_prob_at(double log_prob, size_t row, size_t column);
    void print_parameters(std::ostream& output_stream);
    void build_table();

public:
    Probability_table(
        std::string&& sample_data,                      // Sample data.
        std::vector<double>&& initial_probabilities,    // Probabilities of transition from begin state.
        std::vector<double>&& edges,                    // Probabilities for each edge.
        std::vector<double>&& emission_probabilities,   // Probabilities for each emission for each model.
        size_t (*emission_index)(char));                // Function to map emission to an index in the probability vector.

    void trace_back_and_save(std::ostream& output_stream);
    void print_found_sequences(std::ostream& output_stream, size_t max_hits, size_t min_nucleotide_count);
    size_t count_hits();
    void train_and_print(std::ostream& output_stream);
#ifdef _DEBUG
    void print_dice_rolls(std::ostream& output_stream);
#endif
};

#ifdef _DEBUG
extern const char durbin_dice[];
#endif

#endif

