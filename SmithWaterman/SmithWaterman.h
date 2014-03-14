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

#ifndef SMITHWATERMAN_H
#define SMITHWATERMAN_H

//---------------------------------------------------------------------------
// Definition of a sequence alignment table.
class Alignment_table
{
private:
    std::vector<int> m_score_table; // 2D matrix of scores
    const size_t m_columns;         // width of matrix
    const size_t m_rows;            // height of matrix
    int m_max_score;                // maximum score in this matrix
    std::string m_sequence1;        // represents sequence on j axis
    std::string m_sequence2;        // represents sequence on i axis

    // This is a function that represents the scoring policy (BLOSUM-62 or otherwise)
    int (*m_score_policy)(char char1, char char2);

    // Not implemented to prevent accidental copying.
    Alignment_table(const Alignment_table&);
    Alignment_table& operator=(const Alignment_table&);

protected:
    struct residue_pair
    {
        char residue1;    // residue on j axis
        char residue2;    // residue on i axis
    };

    int score_at(size_t row, size_t column) const;
    void set_score_at(int score, size_t row, size_t column);
    void print_trace_back(std::ostream& output_stream, size_t row, size_t column, std::vector<residue_pair> optimal_alignment) const;
    void print_alignment(std::ostream& output_stream, const std::vector<residue_pair>& alignment) const;

public:
    Alignment_table(const std::string& sequence1, const std::string& sequence2, int (score_policy)(char char1, char char2));
    void print_trace_back(std::ostream& output_stream) const;
    void print_table(std::ostream& output_stream) const;
    void calc_pvalue(std::ostream& output_stream, unsigned int cSequences) const;
};

#endif

