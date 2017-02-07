#include "PreCompile.h"
#include "fasta.h"    // Pick up forward declarations to ensure correctness.

static std::string read_fasta_stream(std::istream& is, size_t size)
{
    std::string sample_data;
    sample_data.reserve(size);

    std::string line;
    while(!std::getline(is, line).eof())
    {
        if(line[0] != '>')
        {
            sample_data.append(line);
        }
    }

    return sample_data;
}

// Reads a FASTA format file.
// https://en.wikipedia.org/wiki/FASTA_format
// This implementation will generally be used for FNA (fasta
// nucleic acid) files.
std::string read_fasta_file(_In_ const char* filename)
{
    const auto size = [filename]() -> size_t
    {
        std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
        const auto size = in.tellg();
        return size <= SIZE_MAX ? static_cast<size_t>(size) : SIZE_MAX;
    }();

    // Read in the sequence data.
    std::ifstream input_file(filename);
    return read_fasta_stream(input_file, size);
}

