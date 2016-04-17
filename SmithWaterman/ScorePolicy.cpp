#include "PreCompile.h"
#include "ScorePolicy.h"    // Pick up forward declarations to ensure correctness.

//---------------------------------------------------------------------------
// BLOSUM 62 score matrix.
// http://en.wikipedia.org/wiki/BLOSUM
constexpr int BLOSUM62_matrix[] = {
//   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22
//   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X

     4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0,     // A
    -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1,     // R
    -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1,     // N
    -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1,     // D
     0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2,     // C
    -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1,     // Q
    -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1,     // E
     0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1,     // G
    -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1,     // H
    -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1,     // I
    -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1,     // L
    -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1,     // K
    -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1,     // M
    -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1,     // F
    -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2,     // P
     1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0,     // S
     0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0,     // T
    -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2,     // W
    -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1,     // Y
     0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1,     // V
    -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1,     // B
    -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1,     // Z
     0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1,     // X
};

constexpr size_t BLOSUM62_width = 23;

// Character mapping from amino acid letter to index in BLOSUM62_matrix to maintain O(1) lookup.
//                                 A   B   C   D   E   F   G   H   I   J   K   L   M   N   O   P   Q   R   S   T   U   V   W   X   Y   Z
constexpr int BLOSUM62_index[] = { 0, 20,  4,  3,  6, 13,  7,  8,  9, -1, 11, 10, 12,  2, -1, 14,  5,  1, 15, 16, -1, 19, 17, 22, 18, 21 };

//---------------------------------------------------------------------------
// Score +2 for a match, -1 for a mismatch.
int basic_calc_score(char char1, char char2)
{
    return char1 == char2 ? 2 : -1;
}

//---------------------------------------------------------------------------
static size_t index_from_char(char char1)
{
    const size_t index = toupper(char1) - 'A';
    assert(index < ('Z' - 'A' + 1));
    return index;
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
    size_t row    = BLOSUM62_index[index_from_char(char1)];
    size_t column = BLOSUM62_index[index_from_char(char2)];

    // If the input strings ever come from an untrusted source, the
    // invalid character error case must be handled instead of asserting
    // (including bounds check on BLOSUM62_index).
    assert(row != -1);
    assert(column != -1);

    return BLOSUM62_matrix[row * BLOSUM62_width + column];
}

