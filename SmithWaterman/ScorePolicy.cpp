/*
Copyright (C) 2006-2014 by Toby Jones.

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

#include "PreCompile.h"
#include "ScorePolicy.h"    // Pick up forward declarations to ensure correctness.

//---------------------------------------------------------------------------
// BLOSUM 62 score matrix.
// http://en.wikipedia.org/wiki/BLOSUM
static const int BLOSUM62_matrix[] = {
//   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
//   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V

     4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0,     // A
    -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3,     // R
    -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,     // N
    -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,     // D
     0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1,     // C
    -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,     // Q
    -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,     // E
     0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3,     // G
    -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,     // H
    -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3,     // I
    -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1,     // L
    -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,     // K
    -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1,     // M
    -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1,     // F
    -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2,     // P
     1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,     // S
     0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0,     // T
    -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3,     // W
    -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1,     // Y
     0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4,     // V
};

static const size_t BLOSUM62_width = 20;

// Character mapping from amino acid letter to index in BLOSUM62_matrix to maintain O(1) lookup.
static const int BLOSUM62_index[] = { 0, -1, 4, 3, 6, 13, 7, 8, 9, -1, 11, 10, 12, 2, -1, 14, 5, 1, 15, 16, -1, 19, 17, -1, 18 };

//---------------------------------------------------------------------------
// Score +2 for a match, -1 for a mismatch.
int basic_calc_score(char char1, char char2)
{
    return char1 == char2 ? 2 : -1;
}

//---------------------------------------------------------------------------
// Scoring policy against the BLOSUM-62 matrix, with a gap penalty.
int BLOSUM62_calc_score_with_penalty(char char1, char char2, int gap_penalty)
{
    static_assert((sizeof(BLOSUM62_matrix) / sizeof(BLOSUM62_matrix[0])) == (BLOSUM62_width * BLOSUM62_width),
                  "BLOSUM62 matrix is an incorrect size.");

    // Check for linear gap cost.
    if((gap_character == char1) || (gap_character == char2))
    {
        return gap_penalty;
    }

    // Lookup score in BLOSUM-62 matrix.
    // ASCII/UTF-8 is assumed.
    size_t row    = BLOSUM62_index[toupper(char1) - 'A'];
    size_t column = BLOSUM62_index[toupper(char2) - 'A'];

    // If the input strings ever come from an untrusted source, the
    // invalid character error case must be handled instead of asserting
    // (including bounds check on BLOSUM62_index).
    assert(row != -1);
    assert(column != -1);

    return BLOSUM62_matrix[row * BLOSUM62_width + column];
}

