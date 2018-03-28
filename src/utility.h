// This file is part of OpenMOOR, an Open-source simulation program for MOORing
// systems in offshore renewable energy applications.
//
// Created by Lin Chen on May 1, 2017.
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


#ifndef utility_h
#define utility_h

#include "Eigen/QR"
#include <math.h>
#include <vector>
#include <iostream>
#include <string>
#include <regex>

namespace moor {
    
/// \brief Unitility functions.

/// Polynomial curve fitting for the given vector y as a function of vector x.
/// This function is used to fit discrete Current data.
void polyfit(const Eigen::VectorXd &x, const Eigen::VectorXd &y,
             Eigen::VectorXd &p, const int order);

/// Linear interpolation for considering non-uniform current profile:
/// interpolating at a vector of points.
double linear_interpolate(const double x, const std::vector<double> &x_grid,
                          const std::vector<double> &y_known);

/// Linear interpolation for considering non-uniform current profile for
/// interpolating at one point.
void linear_interpolate(const std::vector<double> &x, std::vector<double> &y,
                        const std::vector<double> &x_grid,
                        const std::vector<double> &y_known);

/// Check if the string can be converted to a number.
bool is_number(const std::string &str);

} // End of namespace moor.

#endif // utility_h
