// This file is part of OpenMOOR, an Open-source simulation program for MOORing
// systems in offshore renewable energy applications.
//
// Created by Lin Chen on May 17, 2017.
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


#ifndef meshoptimizer_h
#define meshoptimizer_h

#include "utility.h"
#include <vector>
#include <iostream>
#include <iostream>
#include <boost/numeric/odeint.hpp>

namespace moor {
    
/// \brief MeshOptimizer class.
///
class MeshOptimizer
{
    
public:
    
    MeshOptimizer(std::vector<double>& present_mesh,
                  std::vector<double>& reference_variable);
    
    /// Optimize the mesh.
    void optimize(const int n_node, const double weight_factor);
    /// Optimized mesh point coordinates.
    std::vector<double> optimized_mesh;
    
private:
    
    /// Definition of the differential equation for optimization.
    void diff_fun(const double &s, double &dsdq, const double q,
                  const double beta);
    
private:
    
    /// Input mesh point coordinates.
    const std::vector<double> present_mesh;
    const std::vector<double> reference_variable;
    std::vector<double> reference_variable_scaled;
    
    /// Optimization parameter.
    double weight_factor;
    
    /// Observer for saving the intermediate integration.
    struct Observer
    {
        std::vector<double>& states;
        
        Observer(std::vector<double>& state) : states(state)
        {
        };
        void operator() (const double& x, double t)
        {
            states.push_back(x);
        }
    };
};

} // End of namespace moor.

#endif // meshoptimizer_h
