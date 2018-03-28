// This file is part of OpenMOOR, an Open-source simulation program for MOORing
// systems in offshore renewable energy applications.
//
// Created by Lin Chen on Mar 5, 2017.
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


#ifndef node_h
#define node_h

#include "numeric.h"
#include "mooring.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

namespace moor {

class Cable;

/// \brief NodeEA finite difference node.
///
/// Nodes used for approximating the continuous cable. When the local coordinate
/// is parametrized using Euler Angles.
class NodeEA
{
    friend class Cable;
    
public:
    
    // Construct a cable from given node coordinate, cable structural property
    // and the environment it is subjected to.
    NodeEA(double s,
           StructProperty* struct_prop,
           HydroProperty*  hydro_prop,
           SeabedProperty* seabed_prop,
           Constant* constant);
    
    /// Update node state including internal forces, velocity, Euler Angle and
    /// curvature from the state vector after each iteration. The trignometric
    /// functions of Euler Angle is updated as well for further use.
    inline void update_from_state(void)
    {
        // For parameterization using Euler Angles.
        strain            = state(index(0));
        apply_constitutive_relation();
        
        internal_force(1) = state(index(1));
        internal_force(2) = state(index(2));
        velocity(0)       = state(index(3));
        velocity(1)       = state(index(4));
        velocity(2)       = state(index(5));
        euler_angle(0)    = state(index(6));
        euler_angle(1)    = state(index(7));
        curvature(1)      = state(index(8));
        curvature(2)      = state(index(9));
        
        // State dependent.
        sin_angle = euler_angle.array().sin();
        cos_angle = euler_angle.array().cos();
        tan_angle = euler_angle.array().tan();
    };
    
    /// @name Getters.
    ///@{
    double   get_coordinate(void) const { return coordinate; };
    Vector3d get_internal_force(void) const { return internal_force; };
    Vector3d get_velocity(void) const { return velocity; };
    Vector3d get_euler_angle(void) const { return euler_angle; };
    Vector3d get_position(void) const { return position; };
    Vector3d get_global_position(void) const { return global_position; };
    Vector3d get_curvature(void) const { return curvature; };
    double get_total_force(void) const { return sqrt(internal_force.array().square().sum());};
    ///@}

private:
    
    /// Set constant components in mass and stiffness matrices and constant
    /// components in force vector.
    void set_state_independent(void);
    /// Set components in mass and stiffness matrices and in force vector which
    /// depend on the state.
    void update_state_dependent(void);
    /// Consider the seabed effect by modifying the effective cable weight when
    /// grounded on seabed.
    void update_effective_weight(const int is_static);
    /// Reset hydrodynamic property if the cable node is submerged again. The
    /// following parameters need to be modified:
    /// unit_length_weight (buoyancy effect), unit_length_added_mass,
    /// inertia_force_constant, drag_force_constant.
    void update_hydrodynamics(void);
    /// Evaluate the current velocity based on cuurent nodal position and
    /// the current profile function.
    void update_current_velocity(MatrixXd& current_velocity_poly);
    
    inline void apply_constitutive_relation(void)
    {
        internal_force(0) = (strain > 0.0 ?
                             (struct_property->axial_stiffness * strain) : 0);
		axial_stiffness = strain >= 0.0 ? struct_property->axial_stiffness : 0.0;
    };
    /// Set boundary equations for the fixed anchor.
    void add_constraint(void);
    /// Set boundary equations for the fairlead subject to static force.
    void add_constraint(Vector3d& external_force);
    /// Transform the internal force in Lagrangian coordinate into Cartesian.
    Vector3d transform_internal_force(void);
    /// Set boundary equations for the moving fairlead.
    void add_constraint(Vector3d& forced_velocity, double time);
    
private:
    Constant *constant;      ///< Constants.
    double coordinate;       ///< Nodal coordinate.
    double strain;           ///< Epsilon.
    /// The axial stiffness is defined as the tension derivative with respect to
    /// which scan be nonlinear.
    double axial_stiffness;
    
    /// Nodal tangential tension, in-plane and out-of-plane shear forces, i.e.
    /// \f$T, S_n, S_b\f$, the tension strain relationship can be nonlinear.
    Vector3d internal_force;
    Vector3d velocity;       ///< Nodal velocity in moving frame
    Vector3d euler_angle;    ///< Nodal Euler Angle
    Vector3d curvature;      ///< Curvatures and torsional deformation
    Vector3d position;       ///< Nodal position in cable Cartesian coordinate.
    Vector3d global_position;///< Nodal position in global coordinate.
    
    /// Inertial force constant calculated from added mass coefficient and
    /// cable structural property.
    double inertial_force_constant;
    /// Stiffness related to bending and torsional stiffnesses.
    double stiffness_constant;
    /// Drag coefficient multiplied by cable diameter.
    Vector3d drag_force_constant;
    /// The added mass coefficient multiplied by unit length cable displaced
    /// water mass.
    Vector3d unit_length_added_mass;
    Vector3d relative_velocity; ///< Nodal velocity relative to fluid.
    
    const HydroProperty*  hydro_property;  ///< Nodal hydro-property.
    const SeabedProperty* seabed_property; ///< Nodal seabed property.
    const StructProperty* struct_property; ///< Nodal structural property.
    
    /// Trignometric functions of nodal Euler Angles, with three columns
    /// respectively corresponding to sine, cosine and tangent.
    Vector3d sin_angle;
    Vector3d cos_angle;
    Vector3d tan_angle;
    
    Eigen::Matrix<int,   10,1> index;    ///< Index for sorting state vector.
    Eigen::Matrix<double,10,1> state;    ///< Nodal state vector.
    
    /// Current velocity at the node in cable Cartesian coordinate
    Vector3d current_velocity;
    Eigen::Matrix<double,10,10> mass;      ///< Nodal mass matrix.
    Eigen::Matrix<double,10,10> stiffness; ///< Nodal stiffness matrix.
    Eigen::Matrix<double,10,1 > force;     ///< Nodal force vector.
    
    /// Jacobian matrix of the nodal force.
    Eigen::Matrix<double,10,10> force_jacobian;
    
    /// Effective unit length weight of cable considering seabed contact which
    /// is updated in each iteration.
    double unit_length_weight; ///< Wet weight per unit length.
    
    /// 3D matrix for storing the Jacobian of the nodal mass matrix.
    std::vector< Eigen::Matrix<double,10,10> > mass_jacobian;
    /// 3D matrix for storing the Jacobian of the nodal stiffness matrix.
    std::vector< Eigen::Matrix<double,10,10> > stiffness_jacobian;
    /// Nodal equation represented as an augmented matrix for a node at bound.
    Eigen::Matrix<double,10,11> equation;
};

} // End of namespace moor.

#endif // node_h
