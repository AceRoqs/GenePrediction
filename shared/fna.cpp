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

#include <string>
#include <fstream>
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

