#pragma once

//---------------------------------------------------------------------------
// Typedef for use in find_if or lower_bound.
typedef std::vector<std::pair<size_t, size_t>>::value_type entry_type;

//---------------------------------------------------------------------------
size_t index_probability(unsigned int nucleotide);
char normalize_nucleotide(char nucleotide);
std::vector<size_t> read_gbk(_In_ const char* filename);
std::tuple<size_t, std::vector<std::pair<size_t, size_t>>> record_ORFs(const std::string& sample_data);
void print_histogram(
    std::ostream& output_stream,
    const std::string& sample_data,
    const std::vector<size_t>& stop_codons,
    const std::vector<std::pair<size_t, size_t>>& ORFs,
    size_t max_ORF,
    const std::vector<double>& odds_ORF,
    const std::vector<double>& odds_background);

