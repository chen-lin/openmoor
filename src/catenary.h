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


#ifndef catenary_h
#define catenary_h

#include "numeric.h"
#include <iostream>

namespace moor {
    
/// \brief Two-dimensional elastic Catenary.
///
/// A two dimensional elastic catenary is defined by the chord vector, the
/// undeformed cable length (m), cable weight (wet weight) per unit length (N/m),
/// and the axial stiffness (N*m). The forces at the supports are obtained by
/// solving the characteristic nonlinear equation using Newton's iteration. The
/// strain, Euler angle, and curvature along the catenary can then be obtained
/// analytically. The equation in \cite irvine1981cable are used. The
///
/// The catenary solution is not accurate because the cable bending stiffness
/// and the fluid current effect cannot be considered. But the catenary solution
/// is required for initializing the nonlinear cable analysis:
///
///  - The approximate fairlead forces are used for initializing the shooting
///    procedure for placing the fairlead position to the desired spot in the
///    nonlinear analysis considering bending stiffness and current effect;
///  - The approximate fairlead stiffness matrix (2*2) \cite al2016stiffness is
///    also used in the shooting procedure for updating the fairlead tension.
///    Besides, an assumed value of the fairlead tension in the out-of-plane
///    direction is used when there exists current velocity in that direction;
///  - The approximate fairlead stiffness matrix from each cable is assembled to
///    form the approximate mooring system stiffness for updating the floater
///    static location when considering the current effect on the whole mooring
///    system \cite al2016stiffness;
///  - The strain, Euler angle and curvature at discrete locations (node
///    locations) are used for initializing the nonlinear analysis. Besides, the
///    shear force is approximately estimated from the Euler angle and cable
///    bending stiffness.
class Catenary
{
    
public:
    
    Catenary(const Vector2d& chord_vector,
             const double    length,
             const double    unit_length_weight,
             const double    axial_stiffness);
    
    /// Estimate 2*2 stiffness matrix at the cable top end.
    void estimate_fairlead_stiffness(void);
    /// Discretize the cable uniformly using the given number of nodes.
    void discretize(const int    n_node,
                    const double bending_stiffness);
    /// Discretize the cable using the given node number and nodal coordinate.
    void discretize(const int       n_node,
                    const VectorXd& coordinate,
                    const double    bending_stiffness);
    
private:
    
    /// Characteristic equation for solving the bound forces of the catenary.
    void form_nonlinear_equation(void);
    
    /// Characteristic equation for solving fairlead forces when part of the
    /// cable is grounded.
    void form_nonlinear_equation(const int is_grounded);
    
    /// Augmented matrix representing nonlinear equation for solving the higher
    /// support forces.
    Eigen::Matrix<double,2,3> equation;

    const Vector2d chord_vector;       ///< [horizontal vertical]' distances.
    const double   axial_stiffness;    ///< EA in N/m.
    const double   unit_length_weight; ///< Uniform along the cable.
    const double   length;             ///< Undeformed total length.
    
    double total_weight;     ///< Total weight of the cable.
    double suspended_weight; ///< Total weight excluding the grounded part.
    double suspended_length; ///< Suspended length if part grounded
    
public:
    
    /// The discrete states obtained from catenary solution are arranged as
    /// a matrix of size n_node * 5. The five columns are the nodal coordiante,
    /// strain, in-plane shear force, in-plane Euler Angle and in-plane
    /// curvature respectively. All other states needed for initializing the
    /// subsequent analysis are assumed zeros.
    MatrixXd discrete_state;
    Matrix2d bound_force; ///< Forces at two supports: [lower upper].
    /// Fairlead stiffness based on elastic catenary solutions.
    Matrix2d tangent_fairlead_stiffness;
};

} // End of namespace moor.

#endif // catenary_h
