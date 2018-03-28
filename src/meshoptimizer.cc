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


#include "meshoptimizer.h"
using namespace boost::numeric::odeint;
namespace moor {

////////////////////////////////////////////////////////////////////////////////
/// Optimize the node distribution. The segment length between neighboring
/// points is optimized using the method proposed in \cite gobat2006time. A
/// uniform discretization is used as the reference.
////////////////////////////////////////////////////////////////////////////////
MeshOptimizer::MeshOptimizer(std::vector<double>& present_mesh_,
                             std::vector<double>& reference_variable_) :
    present_mesh(present_mesh_),
    reference_variable(reference_variable_)
{
    reference_variable_scaled.resize(reference_variable.size()-1);
    
    for (int i=0; i<reference_variable_scaled.size(); i++)
        reference_variable_scaled[i] = abs(reference_variable[i+1] -
                                           reference_variable[i]);
    
    double ref_var_max = *std::max_element(reference_variable_scaled.begin(),
                                           reference_variable_scaled.end());
    
    for (int i=0; i<reference_variable_scaled.size(); i++)
        reference_variable_scaled[i] = reference_variable_scaled[i] /
        ref_var_max;
    
    weight_factor = 1; // Default value.
};


////////////////////////////////////////////////////////////////////////////////
/// Carry out the optimization.
////////////////////////////////////////////////////////////////////////////////
void MeshOptimizer::optimize(const int n_node, const double weight_factor_)
{
    weight_factor = weight_factor_;
    
    double beta0 = 1.0, beta1 = 1.0 + weight_factor;
    
    auto fun0 = [&] (double &s, double &dsdq, const double q) -> void {
        diff_fun(s, dsdq, q, beta0);};
    
    auto fun1 = [&] (double &s, double &dsdq, const double q) -> void {
        diff_fun(s, dsdq, q, beta1);};
    
    double s0 = 0, s1 = 0;
    
    double dq = present_mesh.back()/(n_node-1.0);
    
    integrate_const(runge_kutta4<double>(), fun0, s0, 0.0, present_mesh.back(), dq);
    
    integrate_const(runge_kutta4<double>(), fun1, s1, 0.0, present_mesh.back(), dq);
    
    double beta2 = beta1, s2, length = present_mesh.back();
    
    while (abs(s1 - length) / length > 1E-10) 
    {
        beta2 = beta1 - (s1 - length) * (beta1 - beta0) / (s1 - s0);
        auto fun2 = [&] (double &s, double &dsdq, const double q) -> void {
            diff_fun(s, dsdq, q, beta2);};
        
        s2 = 0;
        integrate_const(runge_kutta4<double>(), fun2, s2, 0.0,
                        present_mesh.back(), dq);
        beta0 = beta1;
        beta1 = beta2;
        s0    = s1;
        s1    = s2;
    };

    auto fun = [&] (double &s, double &dsdq, const double q) -> void {
        diff_fun(s, dsdq, q, beta2);};
    double s = 0;
    
    integrate_const(runge_kutta4<double>(),
                    fun,
                    s,
                    0.0,
                    present_mesh.back(),
                    dq,
                    Observer(optimized_mesh));
    
    for (int i=0; i<n_node; i++)
        optimized_mesh[i] = optimized_mesh[i] / optimized_mesh.back() * length;
}

////////////////////////////////////////////////////////////////////////////////
/// Differential equation of the optimization problem.
////////////////////////////////////////////////////////////////////////////////
void MeshOptimizer::diff_fun(const double &s, double &dsdq,
                             const double q,
                             const double beta)
{
    double ref_var_interp = linear_interpolate(s, present_mesh,
                                               reference_variable_scaled);
    
    dsdq = beta / (weight_factor * ref_var_interp + 1.0);
};
    
} // End of namespace moor.
