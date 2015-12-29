// Viterbi algorithm for Hidden Markov Models.

#include "PreCompile.h"
#include "Viterbi.h"
#include <Shared/fna.h>

//---------------------------------------------------------------------------
// Remap dice emissions (1-6) to indexes (0-5).
// ASCII/UTF-8 is assumed.
size_t dice_emission_index(char roll)
{
    return roll - '1';
}

//---------------------------------------------------------------------------
// Remap nucleotide emissions (ACGT) to indices (0-3).
// ASCII/UTF-8 is assumed.
size_t nucleotide_emission_index(char nucleotide)
{
    if('A' == nucleotide)
    {
        return 0;
    }
    if('C' == nucleotide)
    {
        return 1;
    }
    if('G' == nucleotide)
    {
        return 2;
    }

    // Treat everything else as a 'T'.
    return 3;
}

//---------------------------------------------------------------------------
int main()
{
#ifdef _DEBUG
    {
        // Exercise the Viterbi algorithm on the dice example in Durbin.
        std::cout << "HMM of Durbin Dice:\n";

        // Two models: fair die/loaded die.
        const size_t model_count = 2;
        std::vector<double> initial_probabilities(model_count);
        initial_probabilities[0] = 0.95;
        initial_probabilities[1] = 0.05;

        // Transition probabilities.
        std::vector<double> edges(model_count * model_count);
        edges[0 * model_count + 0] = 0.95;
        edges[0 * model_count + 1] = 0.05;
        edges[1 * model_count + 0] = 0.1;
        edges[1 * model_count + 1] = 0.9;

        // Emission probabilities.
        const size_t emission_count = 6;
        std::vector<double> emission_probabilities(model_count * emission_count);
        emission_probabilities[0 * emission_count + 0] = (1.0 / 6.0);   // Probabilities of fair dice.
        emission_probabilities[0 * emission_count + 1] = (1.0 / 6.0);
        emission_probabilities[0 * emission_count + 2] = (1.0 / 6.0);
        emission_probabilities[0 * emission_count + 3] = (1.0 / 6.0);
        emission_probabilities[0 * emission_count + 4] = (1.0 / 6.0);
        emission_probabilities[0 * emission_count + 5] = (1.0 / 6.0);
        emission_probabilities[1 * emission_count + 0] = (1.0 / 10.0);  // Probabilities of loaded dice.
        emission_probabilities[1 * emission_count + 1] = (1.0 / 10.0);
        emission_probabilities[1 * emission_count + 2] = (1.0 / 10.0);
        emission_probabilities[1 * emission_count + 3] = (1.0 / 10.0);
        emission_probabilities[1 * emission_count + 4] = (1.0 / 10.0);
        emission_probabilities[1 * emission_count + 5] = (1.0 / 2.0);

        {
            std::string dice(durbin_dice);
            Probability_table table(std::move(dice),
                                    std::move(initial_probabilities),
                                    std::move(edges),
                                    std::move(emission_probabilities),
                                    dice_emission_index);
            table.trace_back_and_save(std::cout);
            table.print_dice_rolls(std::cout);
        }
    }
#endif

    {
        // Exercise the Viterbi algorithm on M. jannaschii.
        std::cout << "HMM Viterbi of M. jannaschii:\n";

        // Initial probabilities.
        // Two models: low G-C base pair content/high G-C base pair content.
        const size_t model_count = 2;
        std::vector<double> initial_probabilities(model_count);
        initial_probabilities[0] = 0.9999;
        initial_probabilities[1] = 0.0001;

        // Transition probabilities.
        std::vector<double> edges(model_count * model_count);
        edges[0 * model_count + 0] = 0.9999;
        edges[0 * model_count + 1] = 0.0001;
        edges[1 * model_count + 0] = 0.01;
        edges[1 * model_count + 1] = 0.99;

        // Emission probabilities.
        const size_t emission_count = 4;
        std::vector<double> emission_probabilities(model_count * emission_count);
        emission_probabilities[0 * emission_count + nucleotide_emission_index('A')] = 0.25; // Probabilities of low GC genomic background.
        emission_probabilities[0 * emission_count + nucleotide_emission_index('C')] = 0.25;
        emission_probabilities[0 * emission_count + nucleotide_emission_index('G')] = 0.25;
        emission_probabilities[0 * emission_count + nucleotide_emission_index('T')] = 0.25;
        emission_probabilities[1 * emission_count + nucleotide_emission_index('A')] = 0.20; // Probabilities of high GC genomic background.
        emission_probabilities[1 * emission_count + nucleotide_emission_index('C')] = 0.30;
        emission_probabilities[1 * emission_count + nucleotide_emission_index('G')] = 0.30;
        emission_probabilities[1 * emission_count + nucleotide_emission_index('T')] = 0.20;

        // Read in the sequence data.
        std::cout << "Reading NC_000909.fna..." << std::endl;

        std::string sample_data;
        read_fna("NC_000909.fna", sample_data);

        std::cout << "Beginning analysis..." << std::endl;
        {
            Probability_table table(std::move(sample_data),
                                    std::move(initial_probabilities),
                                    std::move(edges),
                                    std::move(emission_probabilities),
                                    nucleotide_emission_index);
            table.trace_back_and_save(std::cout);
            table.print_found_sequences(std::cout, 0, 0);

            // Do Viterbi training 10 times.
            for(unsigned int ii = 0; ii < 10; ++ii)
            {
                table.train_and_print(std::cout);
            }

            // Print first 10 sequences of at least 50 nucleotides.
            table.print_found_sequences(std::cout, 10, 50);
        }
    }

    std::cout << "Program done." << std::endl;
    return 0;
}

