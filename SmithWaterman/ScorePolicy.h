#pragma once

const char gap_character = '-';

// Basic scoring policy, with no gap penalty.
// This is simply for testing.
int basic_calc_score(char char1, char char2);

// Scoring policy against the BLOSUM-62 matrix, with a gap penalty.
int BLOSUM62_calc_score_with_penalty(char char1, char char2, int gap_penalty);

template<int GAP_PENALTY>
int BLOSUM62_calc_score(char char1, char char2)
{
    // Redirect to a function that has access to the BLOSUM matrix data.
    return BLOSUM62_calc_score_with_penalty(char1, char2, GAP_PENALTY);
}

