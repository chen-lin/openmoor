// This file is part of OpenMOOR, an Open-source simulation program for MOORing
// systems in offshore renewable energy applications.
//
// Created by Lin Chen on May 5, 2017.
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


#ifndef solver_h
#define solver_h

#include "Eigen/Dense"
#include <iostream>
#include <vector>

namespace moor {

/// \brief Solver uses Newton method with relaxation.
///
/// Solver has been implemented based on the function \a solvde provided in
/// Numerical Recipes \cite press2007numerical. For the parameters defined for
/// using Generalized Alpha method, refer to the papers \cite gobat2000dynamics
/// \cite gobat2001application \cite gobat2002generalized \cite gobat2006time
/// for details. Those parameters are used in Cable class.
class Solver
{
    
public:
    
    Solver(void){};
    void initialize(const int n, const int n_b);
    
    // Adaptive relaxation.
    void adjust_relaxation(double present_error, double prev_error);
    // Manipulate the augmented matrix to find the solution.
    int solve(std::vector< Eigen::MatrixXd >& aug_mat);
    bool check_convergence(const double error) const
    {
        if (error <= convergence_tolerance)
            return true;
        else
            return false;
    }
    /// Reset relaxation factor at the beginning of solving procedure at each
    /// time step.
    void reset_relaxation_factor(void)
    {
        relaxation_factor = initial_relaxation_factor;
    }
    /// Update number of interations taken to solving.
    void update_n_iteration(int it) { n_iteration = it; }
    /// @name Getters.
    ///@{
    int get_iteration_limit(void) const { return n_iteration_limit; }
    double get_relaxation_factor(void) const { return relaxation_factor; }
    int get_n_iteration(void) const { return n_iteration; }
    ///@}
    
    /// Weight factor for stiffness matrix and force vector: \f$\alpha_m\f$.
    double alpha_k;
    /// Weight factor for mass matrix: \f$\alpha_k\f$
    double alpha_m;
    /// Parameter in the generalized trapzoidal rule:
    /// \f$\gamma=\frac{1}{2}-\alpha_m + \alpha_k\f$
    double gamma;
    
    /// @name Constants
    /// They are used in the Generalized-\f$\alpha\f$ method and are calculated
    /// from \f$\alpha_k\f$, \f$\alpha_m\f$, and \f$\gamma\f$.
    ///@{
    double alpha_k1;        ///< \f$1-\alpha_k\f$
    double alpha_k_square;  ///< \f$\alpha_k^2\f$
    double alpha_k1_square; ///< \f$(1-\alpha_k)^2\f$
    double alpha_k_cross;   ///< \f$\alpha_k(1-\alpha_k)\f$
    double alpha_m1;        ///< \f$1-\alpha_m\f$
    double alpha_m_square;  ///< \f$\alpha_m^2\f$
    double alpha_m1_square; ///< \f$(1-\alpha_m)^2\f$
    double alpha_m_cross;   ///< \f$\alpha_k(1-\alpha_k)\f$
    double gamma1;          ///< \f$1-\gamma\f$
    ///@}
    
    int n_iteration_limit; ///< Allowed maximal iteration number.
    double initial_relaxation_factor;
    ///< Initial relaxation factor: normally less than 1.
    double increment_factor;
    ///< Incremental factor in dynamic relaxation.
    double decrement_factor;
    ///< Decremental factor in dynamic relaxation.
    double convergence_tolerance; ///< Convergence criterion.
    /// For using generalized alpha method.
    double lambda_infinity;
    
private:
    int n_nodal_state;        ///< Number of nodal state.
    int n_bound_constraint;   ///< Number of left boundary equations
    int n_iteration;          ///< Save interation number after convergence.
    double relaxation_factor; ///< For dynamic relaxation.
    
    int  gauss_jordan_eliminate(Eigen::MatrixXd & aug_mat,
                                int i_start, int i_end, int j_start);
    void reduce(std::vector< Eigen::MatrixXd >& s, int i);
    void back_substitute(std::vector<Eigen::MatrixXd >& s);
};


} // End of namespace moor.

#endif // solver_h
