// This file is part of OpenMOOR, an Open-source simulation program for MOORing
// systems in offshore renewable energy applications.
//
// Created by Lin Chen on Mar 6, 2017.
//
// Copyright 2018 Lin Chen <l.chen.tj@gmail.com> & Biswajit Basu <basub@tcd.ie>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
////////////////////////////////////////////////////////////////////////////////


#ifndef reader_h
#define reader_h

#include "numeric.h"
#include "utility.h"
#include "setting.h"
#include "input.h"
#include "cable.h"
#include "mooring.h"
#include "moorerror.h"
#include "rapidxml-1.13/rapidxml.hpp"
#include "rapidxml-1.13/rapidxml_print.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

namespace moor {
/// \brief Reader reads input files.
///
/// A group of functions for handling file reading.
class Reader
{
    
public:
    
    /// Read input data from main input file.
    Reader(void) {};
    
    /// Read input data from the XML file. The current, initial cable state and
    /// excitation files are to be found in the same folder. And the outputs are
    /// to be saved in subfolders.
    ErrorCode read_to(InputData& input_data, const std::string input_file_name);
    
    /// Read current velocity and corresponding coordinate.
    int read_to(MatrixXd& data_mat, const std::string data_file_name,
                const int expected_cols, const int skip_rows);
    
    /// Read setting.
    ErrorCode read_to(Setting& setting);
    
private:
    
    // Used when reading main input data file.
    int check_file_existence(const std::string file_name);
    int read_number(const rapidxml::xml_node<>* node);
    int check_availability(const rapidxml::xml_node<>* node,
                           const std::vector<std::string>& names);
    int check_is_number(const rapidxml::xml_node<>* node,
                        const std::vector<std::string>& names);
    int extract_multi_number(const std::string token,
                             std::vector<std::string>& number_string);
};

} // End of namespace moor.

#endif // reader_h
