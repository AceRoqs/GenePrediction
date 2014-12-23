#pragma once

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

