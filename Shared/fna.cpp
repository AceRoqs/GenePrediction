#include "PreCompile.h"
#include "fna.h"    // Pick up forward declarations to ensure correctness.

//---------------------------------------------------------------------------
// Read in fna file and extract the sequence information.
// NOTE: The file is assumed to come from a trusted source and be well-formed.
void read_fna(_In_ const char* filename, std::string& sample_data)
{
    sample_data.reserve(2 * 1024 * 1024);

    // Read in the sequence data.
    std::ifstream input_file(filename);

    std::string line;
    while(!std::getline(input_file, line).eof())
    {
        if(line[0] != '>')
        {
            sample_data.append(line);
        }
    }

    input_file.close();
}

