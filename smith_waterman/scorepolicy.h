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

#ifndef SCOREPOLICY_H
#define SCOREPOLICY_H

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

#endif

