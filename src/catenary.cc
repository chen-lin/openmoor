// This file is part of OpenMOOR, an Open-source simulation program for MOORing
// systems in offshore renewable energy applications.
//
// Created by Lin Chen on May 4, 2017.
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



#include "catenary.h"
namespace moor {

////////////////////////////////////////////////////////////////////////////////
/// Construct a catenary from the chord vector, undeformed cable length, weight
/// per unit length (wet weight) and axial stiffness.
////////////////////////////////////////////////////////////////////////////////
Catenary::Catenary(const Vector2d& chord_vector_, const double length_,
                   const double unit_length_weight_,
                   const double axial_stiffness_) :
        chord_vector(chord_vector_),
        axial_stiffness(axial_stiffness_),
        length(length_),
        unit_length_weight(unit_length_weight_)
{
    equation.setZero();
    tangent_fairlead_stiffness.setZero();
    
    total_weight = unit_length_weight * length;
    
    /// -# Solve boundary forces and suspended length assuming no seabed contact.
    // Initializing guess on the support forces.
    bound_force.setConstant(total_weight / 2.0);
    suspended_length = length;
    
    double solution_error = 1.0;
    while (solution_error > 1.0E-9)
    {
        form_nonlinear_equation();
        Vector2d increment;
        increment = (-equation.block(0,0,2,2).colPivHouseholderQr()
                     .solve(equation.col(2)));
        bound_force.col(1) += increment;
        solution_error = (bound_force.col(1).array().inverse() *
                         increment.array()).abs().sum() / 2.0;
    };
    
    bound_force(1,0) = bound_force(1,1) - total_weight;
    form_nonlinear_equation();
    tangent_fairlead_stiffness = equation.block(0,0,2,2).inverse();
    
    const int is_grounded = 1;
    
    /// -# Check if part of the catenary is below seabed. If yes resolve the
    /// boundary forces using a different set of equation.
    if ((bound_force(1,1) < total_weight) && is_grounded)
    {
        bound_force.col(1).setConstant(total_weight / 2.0);
        solution_error = 1.0;
        
        while (solution_error > 1.0E-9)
        {
            form_nonlinear_equation(is_grounded);
            Vector2d increment;
            increment = (-equation.block(0,0,2,2).colPivHouseholderQr()
                         .solve(equation.col(2)));
            bound_force.col(1) += increment;
            solution_error = (bound_force.col(1).array().inverse() *
                              increment.array()).abs().sum() / 2.0;
            
        };
        suspended_length = bound_force(1,1) / unit_length_weight;
        bound_force(1,0) = 0.0;
        
        // Tangential stiffness is estimated at the mean time.
        form_nonlinear_equation(is_grounded);
        tangent_fairlead_stiffness = equation.block(0,0,2,2).inverse();
    };
    
    bound_force(0,0) = -bound_force(0,1);
    suspended_weight =  suspended_length * unit_length_weight;
};


////////////////////////////////////////////////////////////////////////////////
/// Calculate discrete state based on catenary solution according to coordinate
/// vector.
////////////////////////////////////////////////////////////////////////////////
void Catenary::discretize(const int n_node, const VectorXd &coordinate,
                            const double bending_stiffness)
{
    assert(coordinate.size()==n_node);
    discrete_state.resize(n_node, 5);
    discrete_state.setZero();
    discrete_state.col(0) = coordinate;
    
    double s, T;
    for (int k=0; k<n_node; k++) 
    {
        s = discrete_state(k,0);
        
        if (s < (length - suspended_length))
        {
            T = bound_force(0,1);
            discrete_state(k,2) = 0.0;
            discrete_state(k,4) = 0.0;
        }
        else
        {
            T = sqrt(pow(bound_force(0,1),2.0)
                     + pow((bound_force(1,1) - (length - s)
                            * unit_length_weight),2.0));
            discrete_state(k,4) = (-bound_force(0,1) * unit_length_weight
                                   / pow(T, 2.0));
            discrete_state(k,2)
            = (2.0 * bound_force(0,1) * pow(unit_length_weight,2.0)
               * (bound_force(1,1) - (length - s)*unit_length_weight)
               / pow(T,4.0) * (-bending_stiffness)
               / pow(1 + discrete_state(k,1),3.0));
        }
        // Get strain.
        discrete_state(k,1) = T / axial_stiffness;
        
        // Get Euler Angle.
        discrete_state(k,3) = asin(bound_force(0,1) / T);
    }
}


////////////////////////////////////////////////////////////////////////////////
/// Calculate discrete state based on catenary solution using uniform segment
/// length.
////////////////////////////////////////////////////////////////////////////////
void Catenary::discretize(const int n_node, const double bending_stiffness)
{
    VectorXd coordinate;
    coordinate.resize(n_node);
    coordinate.setZero();
    for (int k=0; k<n_node; k++)
        coordinate(k) = length / (n_node-1) * k;
    
    discretize(n_node, coordinate, bending_stiffness);
}

////////////////////////////////////////////////////////////////////////////////
/// This method is based on Jain (1980) which is for inextensive cables.
////////////////////////////////////////////////////////////////////////////////
void Catenary::estimate_fairlead_stiffness(void)
{
    double contact_angle =
    atan((bound_force(1,1) - suspended_weight) / bound_force(0,1));
    
    double H = bound_force(0,1); // Of the same along the catenary.
    
    // Subscript _a represents the variables at the lower anchor and subscript
    // _f represents the top anchor.
    double T_a = H / cos(contact_angle);
    
    // Get the hypothetical anchor where the contact angle is zero.
    double V_a = T_a * sin(contact_angle);
    
    double added_length = V_a / unit_length_weight;
    
    // Note vertical force at the anchor is balanced by additional cable.
    Vector2d imaginary_anchor_location;
    double c1 = H / unit_length_weight;
    
    imaginary_anchor_location(0) = c1 * asinh(1/c1 * added_length);
    imaginary_anchor_location(1) = c1 * (sqrt(1+pow(1/c1 * added_length,2.0))
                                         - 1.0);
    
    double total_depth  = chord_vector(1) + imaginary_anchor_location(1);
    double total_length = sqrt(pow(total_depth / c1+1.0,2.0)-1) * c1;
    double total_radius = c1 * asinh(1/c1 * total_length);
    double radius       = total_radius - imaginary_anchor_location(0);
    
    double T_f = sqrt( H*H + pow(unit_length_weight * total_length,2.0));
    
    // Get the fairlead stiffnesses.
    tangent_fairlead_stiffness.setZero();
    
    double c2 = T_a * total_length - T_f * added_length;
    
    tangent_fairlead_stiffness(0,1) = (H/( unit_length_weight * c2 / (T_f - T_a)
                                          * (radius/ H - c2/(T_f * T_a))
                                          - pow(H,2.0) * chord_vector(1)
                                          / (T_f * T_a) ));
    
    tangent_fairlead_stiffness(1,0) = tangent_fairlead_stiffness(0,1);
    
    tangent_fairlead_stiffness(0,0) = (tangent_fairlead_stiffness(0,1)/c1
                                       * c2 / (T_f - T_a));
    
    tangent_fairlead_stiffness(1,1) = (tangent_fairlead_stiffness(0,1)
                                       /c1 * (T_f * T_a) / (T_f - T_a)
                                       * (radius/H - c2/(T_f * T_a)));
}


////////////////////////////////////////////////////////////////////////////////
// Nonlinear characteristic equation for case when no part is grounded.
////////////////////////////////////////////////////////////////////////////////
void Catenary::form_nonlinear_equation(void)
{
    double W  = total_weight;
    double H  = bound_force(0,1);
    double V2 = bound_force(1,1);
    double r  = V2/H, r1 = (V2- W)/H;
    
    equation(0,2) = (H * length / axial_stiffness
                     + H / unit_length_weight
                     * (asinh(r) - asinh(r1)) - chord_vector(0));
    
    equation(1,2) = (W * length / axial_stiffness * (V2/W - 0.5)
                     + H / unit_length_weight * (sqrt(1+r*r) - sqrt(1+r1*r1))
                     - chord_vector(1));
    
    // Flexibility matrix.
    equation(0,0) = (length / axial_stiffness + length / W
                     * ( asinh(r) - asinh(r1) - r/sqrt(1+r*r) + r1/(1+r1*r1)));
    
    equation(0,1) = length / W * (1/sqrt(1+r*r) - 1/sqrt(1+r1*r1));
    
    equation(1,0) = equation(0,1);

    equation(1,1) = (length / axial_stiffness + length/ W
                     * ( r/ sqrt(1+r*r) - r1/sqrt(1+r1*r1) ));
}


////////////////////////////////////////////////////////////////////////////////
// Nonlinear characteristic equation when part of the cable is grounded.
////////////////////////////////////////////////////////////////////////////////
void Catenary::form_nonlinear_equation(const int is_grounded)
{
    if (is_grounded==0) return;
    
    double W  = total_weight;
    double H  = bound_force(0,1);
    double V2 = bound_force(1,1);
    double r  = V2/H;
    
    equation(0,2) = ((1 - V2/W) * (H * length / axial_stiffness + length)
                     + H * V2/ axial_stiffness /unit_length_weight
                     + H / unit_length_weight * asinh(r) - chord_vector(0));
    
    equation(1,2) = (V2 * V2 /2.0 / axial_stiffness /unit_length_weight
                     + H / unit_length_weight * (sqrt(1+r*r) - 1.0)
                     - chord_vector(1));
    
    // Flexibility matrix.
    equation(0,0) = (length / axial_stiffness + 1.0/unit_length_weight
                     * asinh(r) - 1.0/unit_length_weight * r/ sqrt(1 + r*r));
    
    equation(0,1) = 1.0/unit_length_weight*(1.0 / sqrt(1 + r*r) - 1.0);
    
    equation(1,0) = equation(0,1);
    
    equation(1,1) = (V2 / axial_stiffness / unit_length_weight
                     + 1.0/ unit_length_weight * r / sqrt(1+r*r));
}
    
} // End of namespace moor.
