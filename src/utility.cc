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


#include "utility.h"

namespace moor {
    
////////////////////////////////////////////////////////////////////////////////
// The following code is downloaded from the website as below
//   http://svn.clifford.at/handicraft/2014/polyfit/ which is referred herein
//void polyfit(const std::vector<double> &xv, const std::vector<double> &yv,
//             std::vector<double> &coeff, int order)
//{
//    Eigen::MatrixXd A(xv.size(), order+1);
//    Eigen::VectorXd yv_mapped = Eigen::VectorXd::Map(&yv.front(), yv.size());
//    Eigen::VectorXd result;
//
//    assert(xv.size() == yv.size());
//    assert(xv.size() >= order+1);
//
//    // create matrix
//    for (size_t i = 0; i < xv.size(); i++)
//        for (size_t j = 0; j < order+1; j++)
//            A(i, j) = pow(xv.at(i), j);
//
//    // solve for linear least squares fit
//    result = A.householderQr().solve(yv_mapped);
//
//    coeff.resize(order+1);
//    for (size_t i = 0; i < order+1; i++)
//        coeff[i] = result[i];
//}
////////////////////////////////////////////////////////////////////////////////
void polyfit(const Eigen::VectorXd &x, const Eigen::VectorXd &y,
             Eigen::VectorXd &p, const int order)
{
    Eigen::MatrixXd A(x.rows(), order+1);
    
    assert(x.rows() == y.rows());
    assert(x.rows() >= order+1);
    
    // Create coefficient matrix.
    for (size_t i = 0; i < x.rows(); i++)
        for (size_t j = 0; j < order+1; j++)
            A(i, j) = pow(x(i), j);
    
    // Solve using linear least squares fitting.
    p.resize(order+1);
    p = A.householderQr().solve(y);
}

////////////////////////////////////////////////////////////////////////////////
// Linear interpolation for one point.
////////////////////////////////////////////////////////////////////////////////
double linear_interpolate(const double x, const std::vector<double> &x_grid,
                          const std::vector<double> &y_known)
{
    int i = 0;
    while (x > x_grid[++i]) { if (i > (x_grid.size()-1)) break; }
    
    double y = (y_known[i-1] + (x - x_grid[i-1]) / (x_grid[i] - x_grid[i-1])
                * (y_known[i] - y_known[i-1]));

    return y;
}


////////////////////////////////////////////////////////////////////////////////
// Linear interpolation for a set of points.
////////////////////////////////////////////////////////////////////////////////
void linear_interpolate(const std::vector<double> &x, std::vector<double> &y,
                        const std::vector<double> &x_grid,
                        const std::vector<double> &y_known)
{
    assert(x_grid.size() == y_known.size());
    
    for (int k=0; k<x.size(); k++)
    {
        int i=0;
        
        while (x[k] > x_grid[++i] ) { if (i > (x_grid.size()-1)) break; };
        
        y[k] = (y_known[i-1] + (x[k] - x_grid[i-1]) / (x_grid[i] - x_grid[i-1])
                * (y_known[i] - y_known[i-1]));
    }
}


////////////////////////////////////////////////////////////////////////////////
// For check input data.
////////////////////////////////////////////////////////////////////////////////
bool is_number(const std::string &token)
{
    return std::regex_match(token, std::regex("[+\\-]?(?:0|[1-9]\\d*)(?:\\.\\d*)?(?:[eE][+\\-]?\\d+)?$"));
}
    
} // End of namespace moor.
