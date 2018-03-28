// This file is part of OpenMOOR, an Open-source simulation program for MOORing
// systems in offshore renewable energy applications.
//
// Created by Lin Chen on Aug 25, 2017.
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


#ifndef writer_h
#define writer_h

#include "cable.h"
#include "platform.h"
#include "moorerror.h"
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip> // setprecision
#include <sstream> // stringstream

namespace moor {
    
/// \brief Writer writes to files.
///
/// A group functions used to output analysis results of the platform and
/// cables into .dat files.
class Writer
{
public:
    Writer(void) {};
    
    // Write platform file header.
    void write(const std::string target_file_name, const Platform& platform);
    
    // Write platform state to file in dynamic analysis.
    void write(const std::string target_file_name, const Platform& platform,
               const double time);
    
    // Write the static analysis result using dynamic relaxation.
    void write(const std::string target_file_name, const Platform& platform,
               const int iteration_step);
     
    // Write cable state result.
    void write(const std::string target_file_name, const Cable& cable);
    
    // Write log file to store the error message.
    int write_error_message(std::string file_name, ErrorCode err,
                            ErrorCategory err_cat);
};

} // End of namespace moor.

#endif // writer_h
