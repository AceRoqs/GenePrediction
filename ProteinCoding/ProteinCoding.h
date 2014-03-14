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

#ifndef PROTEINCODING_H
#define PROTEINCODING_H

//---------------------------------------------------------------------------
// Typedef for use in find_if.
typedef std::multimap<size_t, size_t>::value_type entry_type;

//---------------------------------------------------------------------------
size_t index_probability(unsigned int nucleotide);
char normalize_nucleotide(char nucleotide);
void read_gbk(_In_ const char* filename, std::map<size_t, size_t>& genes);
size_t record_ORFs(const std::string& sample_data, std::multimap<size_t, size_t>& ORFs);
void print_histogram(
    std::ostream& output_stream,
    const std::string& sample_data,
    const std::map<size_t, size_t>& genes,
    const std::multimap<size_t, size_t>& ORFs,
    size_t max_ORF,
    const std::vector<double>& odds_ORF,
    const std::vector<double>& odds_background);

#endif

