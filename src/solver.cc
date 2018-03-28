// This file is part of OpenMOOR, an Open-source simulation program for MOORing
// systems in offshore renewable energy applications.
//
// Created by Lin Chen on Sep 15, 2017.
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


#include "solver.h"

namespace moor {

////////////////////////////////////////////////////////////////////////////////
/// When basic parameters are obatined from the input data. Initialize the other
/// parameters.
////////////////////////////////////////////////////////////////////////////////
void Solver::initialize(const int n, const int n_b)
{
    n_nodal_state = n;
    n_bound_constraint = n_b;
    
    n_iteration = 0;
    relaxation_factor = initial_relaxation_factor;
    alpha_k = lambda_infinity/(lambda_infinity - 1.0);
    alpha_m = (3*lambda_infinity+1)/2/(lambda_infinity - 1.0);
    gamma = 0.5 - alpha_m + alpha_k;
    
    // Initialization constants for use.
    alpha_k1 = 1-alpha_k;
    alpha_k_square = pow(alpha_k, 2.0);
    alpha_k1_square = pow(alpha_k1, 2.0);
    alpha_k_cross = alpha_k * alpha_k1;
    alpha_m1 = 1 - alpha_m;
    alpha_m_square = pow(alpha_m, 2.0);
    alpha_m1_square = pow(alpha_m1, 2.0);
    alpha_m_cross = alpha_m * alpha_m1;
    gamma1 = 1 - gamma;
}
    

////////////////////////////////////////////////////////////////////////////////
/// Solve the augmented matrix for state increment, by three steps:
///   - Reduction;
///   - Gauss-Jordan elimination;
///   - Back substitution.
////////////////////////////////////////////////////////////////////////////////
int Solver::solve(std::vector< Eigen::MatrixXd >& aug_mat)
{
    int fail = 0, n_node = aug_mat.size() - 1;
    
    fail = gauss_jordan_eliminate(aug_mat[0],
                                  n_bound_constraint,
                                  n_nodal_state-1,
                                  n_nodal_state);
    
    int k = 1;
    for (k=1; k<n_node; k++)
    {
        reduce(aug_mat, k);
        fail = gauss_jordan_eliminate(aug_mat[k], 0,
                                      n_nodal_state-1,
                                      n_bound_constraint);
    }
    
    reduce(aug_mat,n_node);
    fail = gauss_jordan_eliminate(aug_mat[n_node],
                                  0, (n_nodal_state-n_bound_constraint)-1,
                                  n_bound_constraint);
    back_substitute(aug_mat);

    // Sort the last column of the augmented matrix.
    for (k=0; k<n_node; k++)
    {
        aug_mat[k].block(0,2*n_nodal_state,n_bound_constraint,1)
        = aug_mat[k].block(n_nodal_state-n_bound_constraint,2*n_nodal_state,n_bound_constraint,1);
        aug_mat[k].block(n_bound_constraint,2*n_nodal_state,n_nodal_state-n_bound_constraint,1)
        = aug_mat[k+1].block(0,2*n_nodal_state,n_nodal_state-n_bound_constraint,1);
    }
    return fail;
}


////////////////////////////////////////////////////////////////////////////////
/// Adjust the relaxation factor according to the error change trend.
////////////////////////////////////////////////////////////////////////////////
void Solver::adjust_relaxation(double present_error, double prev_error)
{
    if (present_error > prev_error || prev_error == 0.0)
        relaxation_factor = relaxation_factor / decrement_factor;
    else
        relaxation_factor = relaxation_factor * increment_factor;
    
    relaxation_factor = (relaxation_factor < initial_relaxation_factor ?
                         relaxation_factor : initial_relaxation_factor);
    relaxation_factor = relaxation_factor > 1E-5 ? relaxation_factor : 1E-5;
}

////////////////////////////////////////////////////////////////////////////////
/// Diagonalize the square block of augmented matrix by Gauss Jordan Elimination
/// using pivoting.
/// <pre>
///         0 0 0 X X X X X X X B           0 0 0 1 0 0 0 0 S S C
///         0 0 0 X X X X X X X B           0 0 0 0 1 0 0 0 S S C
///         0 0 0 X X X X X X X B     =>    0 0 0 0 0 1 0 0 S S C
///         0 0 0 X X X X X X X B           0 0 0 0 0 0 1 0 S S C
///         0 0 0 X X X X X X X B           0 0 0 0 0 0 0 1 S S C
/// </pre>
/// The rows are swapped to form a block diagonal matrix.
////////////////////////////////////////////////////////////////////////////////
int  Solver::gauss_jordan_eliminate(Eigen::MatrixXd& aug_mat,
                                    int i_start, int i_end,
                                    int j_start)
{
    int n_dim = i_end - i_start + 1;
    int i_pivot, j_pivot; double pivot;
    Eigen::RowVectorXd temp_row;
    
    for (int k=0; k<n_dim; k++) 
    {
        pivot = aug_mat.block(i_start+k,j_start+k,n_dim-k,1).array().
        abs().maxCoeff(&i_pivot,&j_pivot);
        
        if (pivot <= 1E-20)
            return 1; // Singularity. // assert(pivot > 1E-20);
        
        // Swap rows.
        if (i_pivot!=0) 
        {
            temp_row = aug_mat.row(i_start+k+i_pivot);
            aug_mat.row(i_start+k+i_pivot) = aug_mat.row(i_start+k);
            aug_mat.row(i_start+k) = temp_row;
        }
        
        aug_mat.row(i_start+k) = aug_mat.row(i_start+k) *
        (1.0/aug_mat(i_start+k,j_start+k));
        aug_mat(i_start+k,j_start+k) = 1.0;
        
        // Elimination.
        for (int i=k+1; i<n_dim; i++) 
        {
            aug_mat.row(i_start+i) -= aug_mat.row(i_start+k) *
            aug_mat(i_start+i,j_start+k);
        }
        
        // Set zeros.
        if (n_dim > 1)
            aug_mat.block(i_start+k+1,j_start+k,n_dim-k-1,1).setZero();
    }
    
    // Back substitution.
    for (int j=n_dim-1; j>0; j--) 
    {
        for (int i=j-1; i>=0; i--) 
        {
            aug_mat.row(i_start+i) -= aug_mat.row(i_start+j) *
            aug_mat(i_start+i,j_start+j);
        }
        aug_mat.block(i_start,j_start+j,j,1).setZero();
    }
    return 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Reduce columns jz1 .. jz2-1 of the s matrix, using previous results as
/// stored in the c matrix. Only colums jm1 .. jm2-1 and jmf are affected by
/// prior results.
/// <pre>
///        X X X X X X X X X X B           0 0 0 S S S S S S S C
///        X X X X X X X X X X B           0 0 0 S S S S S S S C
///        X X X X X X X X X X B     =>    0 0 0 S S S S S S S C
///        X X X X X X X X X X B           0 0 0 S S S S S S S C
///        X X X X X X X X X X B           0 0 0 S S S S S S S C
/// </pre>
////////////////////////////////////////////////////////////////////////////////
void Solver::reduce(std::vector< Eigen::MatrixXd >& s,int i)
{
    // Alter the columns of the coefficient matrix.
    s[i].block(0,n_bound_constraint,n_nodal_state,n_nodal_state-n_bound_constraint)
    -= (s[i].block(0,0,n_nodal_state,n_bound_constraint)
        * s[i-1].block(n_nodal_state-n_bound_constraint,
                       n_nodal_state+n_bound_constraint,
                       n_bound_constraint,n_nodal_state-n_bound_constraint));
    
    // Alter the b column.
    s[i].col(2*n_nodal_state)
    -= (s[i].block(0,0,n_nodal_state,n_bound_constraint)
        * s[i-1].block(n_nodal_state-n_bound_constraint,
                       2*n_nodal_state,n_bound_constraint,1));
    
    // For testing.
    s[i].block(0,0,n_nodal_state,n_bound_constraint)
    -= (s[i].block(0,0,n_nodal_state,n_bound_constraint)
        * s[i-1].block(n_nodal_state-n_bound_constraint,
                       n_nodal_state,n_bound_constraint,n_bound_constraint));
}


////////////////////////////////////////////////////////////////////////////////
/// Back substitute to dealing with the following structure
/// <pre>
///    1 X X            V  B
///      1         X X  V  B
///        1       X X  V  B
///          1     X X  V  B
///            1   X X  V  B
///              1 X X  V  B
///                1    V  B
///                  1  V  B
/// </pre>
/// Note: Values of B after back substitution are the solution.
////////////////////////////////////////////////////////////////////////////////
void Solver::back_substitute(std::vector< Eigen::MatrixXd >& s)
{
    int n = s.size();
    for (int i=n-2; i>=0; i--)
    {
        s[i].col(2*n_nodal_state)
        -= (s[i].block(0,n_nodal_state+n_bound_constraint,
                       n_nodal_state, n_nodal_state-n_bound_constraint) *
            s[i+1].block(0,2*n_nodal_state,n_nodal_state-n_bound_constraint,1));
    }
}

} // End of namespace moor.
