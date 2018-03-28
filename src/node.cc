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


#include "node.h"
namespace moor {

////////////////////////////////////////////////////////////////////////////////
/// Node is defined by its arc-length coordiante \f$s\f$ and the associated
/// structural property, hydrodynamic property, seabed contact property. The
/// constants to be used in evaluating nodal mass and stiffness matrices and
/// nodal force are also defined.
////////////////////////////////////////////////////////////////////////////////
NodeEA::NodeEA(double s,
               StructProperty* struct_prop,
               HydroProperty*  hydro_prop,
               SeabedProperty* seabed_prop,
               Constant *c) :
    coordinate(s),
    struct_property(struct_prop),
    hydro_property(hydro_prop),
    seabed_property(seabed_prop),
    constant(c)
{
    // Froude-Krylov component.
    
    inertial_force_constant = ((1 + hydro_property->added_mass_coefficient(1))
                               * M_PI * pow(struct_prop->diameter, 2.0)/4.0
                               * constant->WATER_DENSITY);
    
    // Ignoring the Froude-Krylov component for comparison purpose.
    axial_stiffness = struct_prop->axial_stiffness;
    
    drag_force_constant = (constant->WATER_DENSITY * struct_prop->diameter
                           * hydro_property->drag_coefficient.array());
    
    unit_length_added_mass = (constant->WATER_DENSITY * M_PI
                              * pow(struct_prop->diameter,2.0) / 4.0
                              * hydro_property->added_mass_coefficient.array());
    
    unit_length_weight = struct_property->unit_length_weight;
    
    stiffness_constant = (struct_property->bending_stiffness -
                          struct_property->torsional_stiffness);

    // Set sizes of vectors, matrices and tensors.
    internal_force.setZero();
    velocity.setZero();
    euler_angle.setZero();
    sin_angle.setZero();
    cos_angle.setZero();
    tan_angle.setZero();
    
    curvature.setZero();
    position.setZero();
    global_position.setZero();
    
    current_velocity.setZero();
    relative_velocity.setZero();
    equation.setZero();
    
    // For the case of using Euler Angles for parametrization to be changed
    // when Quaternions are used.
    index << 5,6,7,2,3,4,8,9,0,1;
    state.setZero();
    mass.setZero();
    stiffness.setZero();
    force.setZero();
    force_jacobian.setZero();
    mass_jacobian.resize(10);
    stiffness_jacobian.resize(10);
    
    for (int i=0; i < 10; i++)
    {
        mass_jacobian[i].setZero();
        stiffness_jacobian[i].setZero();
    }
}


////////////////////////////////////////////////////////////////////////////////
/// Set the time independents in mass and stiffness matrices.
////////////////////////////////////////////////////////////////////////////////
void NodeEA::set_state_independent()
{
    /// -# Set mass matrix components.
    mass(0,index(0)) =  struct_property->damping_coefficient;
    mass(3,index(0)) = -1.0;
    mass(0,index(3)) = -((struct_property->unit_length_mass
                          + unit_length_added_mass(0)));
    mass(1,index(4)) = -(struct_property->unit_length_mass
                         + unit_length_added_mass(1));
    mass(2,index(5)) = -(struct_property->unit_length_mass
                         + unit_length_added_mass(2));
    
    /// -# Set stiffness matrix components.
    stiffness(0,index(0)) = axial_stiffness;
    for (int i=1; i<8; i++) stiffness(i,index(i)) = 1.0;
    stiffness(8,index(8)) = struct_property->bending_stiffness;
    stiffness(9,index(9)) = struct_property->bending_stiffness;
}


////////////////////////////////////////////////////////////////////////////////
/// Update state dependent components of nodal matrices and tensors.
////////////////////////////////////////////////////////////////////////////////
void NodeEA::update_state_dependent()
{
    /// -# Update relative velocity based on the updated current velocity.
    relative_velocity = (velocity-AngleAxisd(-euler_angle(1), Vector3d::UnitY())
                         * AngleAxisd(-euler_angle(0),Vector3d::UnitZ())
                         * current_velocity);
    
    // Define derived variabels.
    double sqrt_strain = sqrt(1.0 + strain);
    double sqrt_vel = sqrt(pow(relative_velocity(1),2.0)
                           + pow(relative_velocity(2),2.0));
    
    /// -# Update mass matrix.
    mass(0,index(6)) = (struct_property->unit_length_mass * velocity(1)
                        * cos_angle(1));
    
    mass(1,index(6)) = -(struct_property->unit_length_mass *
                         (velocity(2) * sin_angle(1) + velocity(0) * cos_angle(1))
                         + inertial_force_constant *
                         (current_velocity(0) * cos_angle(0)
                          + current_velocity(1) * sin_angle(0)));
    
    mass(2,index(6)) = (struct_property->unit_length_mass * velocity(1)
                        * sin_angle(1) + inertial_force_constant
                        * (-current_velocity(0) * sin_angle(0) * sin_angle(1)
                           + current_velocity(1) * cos_angle(0) * sin_angle(1)));
    
    mass(4,index(6)) = -(1.0+strain) * cos_angle(1);
    
    mass(0,index(7)) = (-struct_property->unit_length_mass * velocity(2));
    
    mass(2,index(7)) =  (struct_property->unit_length_mass * velocity(0)
                         + inertial_force_constant
                         * (current_velocity(0) * cos_angle(0) * cos_angle(1)
                            + current_velocity(1) * sin_angle(0) * cos_angle(1)
                            - current_velocity(2) * sin_angle(1)));
    
    mass(5,index(7)) = 1.0 + strain;
    
    /// -# Update stiffness matrix.
    stiffness(0,index(6)) = -internal_force(1) * cos_angle(1);
    
    stiffness(1,index(6)) =  (internal_force(0) * cos_angle(1)
                              + internal_force(2) * sin_angle(1));
    
    stiffness(2,index(6)) = -internal_force(1) * sin_angle(1);
    
    stiffness(3,index(6)) = -velocity(1) * cos_angle(1);
    
    stiffness(4,index(6)) =  (velocity(0) * cos_angle(1) +
                              velocity(2) * sin_angle(1));
    
    stiffness(5,index(6)) = -velocity(1) * sin_angle(1);
    
    stiffness(6,index(6)) =  cos_angle(1);
    
    stiffness(0,index(7)) =  internal_force(2);
    
    stiffness(2,index(7)) = -internal_force(0);
    
    stiffness(3,index(7)) =  velocity(2);
    
    stiffness(5,index(7)) = -velocity(0);
    
    /// -# Update nodal force. Not that here the weight per unit length is the
    /// effective weight considering seabed effect.
    force(0) = (-unit_length_weight * cos_angle(0) * cos_angle(1)
                - 0.5 * drag_force_constant(0)
                * M_PI * relative_velocity(0) * fabs(relative_velocity(0))
                * sqrt_strain );
    
    force(1) = (unit_length_weight * sin_angle(0)
                - 0.5 * drag_force_constant(1) * relative_velocity(1)
                * sqrt_vel * sqrt_strain );
    
    force(2) = (-unit_length_weight * cos_angle(0) * sin_angle(1)
                - 0.5 * drag_force_constant(2) * relative_velocity(2)
                * sqrt_vel * sqrt_strain);
    
    force(6) = -curvature(2);
    
    force(7) = -curvature(1);
    
    force(8) = (stiffness_constant * tan_angle(1)
                * pow(curvature(2),2.0) - internal_force(2)
                * pow((1.0 + strain),3.0));
    
    force(9) = (-stiffness_constant * tan_angle(1)
                * curvature(2) * curvature (1)
                + internal_force(1) * pow((1.0 + strain),3.0));
    
    /// -# Update the tensor matrix related to mass matrix.
    // Each slice stores the Jacobian matrix of the row vector of the mass
    // matrix with respect to the state.
    mass_jacobian[0](index(6),index(4)) = (struct_property->unit_length_mass
                                           * cos_angle(1));
    
    mass_jacobian[0](index(6),index(7)) = (-struct_property->unit_length_mass
                                           * sin_angle(1) * velocity(1));
    
    mass_jacobian[0](index(7),index(5)) = -struct_property->unit_length_mass;
    
    mass_jacobian[1](index(6),index(3)) = (-struct_property->unit_length_mass
                                           * cos_angle(1));
    
    mass_jacobian[1](index(6),index(5)) = (-struct_property->unit_length_mass
                                           * sin_angle(1));
    
    mass_jacobian[1](index(6),index(6)) = (inertial_force_constant
                                           * (current_velocity(0) * sin_angle(0)
                                              - current_velocity(1)
                                              * cos_angle(0)) );
    
    mass_jacobian[1](index(6),index(7)) = (-struct_property->unit_length_mass
                                           * velocity(2) * cos_angle(1)
                                           + struct_property->unit_length_mass
                                           * velocity(0) * sin_angle(1));
    
    mass_jacobian[2](index(6),index(4)) = (struct_property->unit_length_mass
                                           * sin_angle(1));
    
    mass_jacobian[2](index(6),index(6)) = (-inertial_force_constant
                                           * current_velocity(0)
                                           * cos_angle(0) * sin_angle(1)
                                           - inertial_force_constant
                                           * current_velocity(1)
                                           * sin_angle(0) * sin_angle(1));
    
    mass_jacobian[2](index(6),index(7)) = (struct_property->unit_length_mass
                                           * velocity(1) * cos_angle(1)
                                           - inertial_force_constant
                                           * current_velocity(0) * sin_angle(0)
                                           * cos_angle(1)
                                           + inertial_force_constant
                                           * current_velocity(1) * cos_angle(0)
                                           * cos_angle(1));
    
    mass_jacobian[2](index(7),index(3)) =  struct_property->unit_length_mass;
    
    mass_jacobian[2](index(7),index(6)) = (-inertial_force_constant
                                           * current_velocity(0)
                                           * sin_angle(0) * cos_angle(1)
                                           + inertial_force_constant
                                           * current_velocity(1)
                                           * cos_angle(0) * cos_angle(1));
    
    mass_jacobian[2](index(7),index(7)) = (-inertial_force_constant
                                           * current_velocity(0)
                                           * cos_angle(0) * sin_angle(1)
                                           - inertial_force_constant
                                           * current_velocity(1)
                                           * sin_angle(0) * sin_angle(1)
                                           - inertial_force_constant
                                           * current_velocity(2)
                                           * cos_angle(1));
    
    mass_jacobian[4](index(6),index(0)) = -cos_angle(1);
    
    mass_jacobian[4](index(6),index(7)) = (1.0+strain)* sin_angle(1);
    
    mass_jacobian[5](index(7),index(0)) = 1.0;
    
    /// -# Update the tensor related to Jacobain of stiffness term.
    // Each slice stores the Jacobian matrix of the row vector of the stiffness
    // matrix with respect to the state.
    stiffness_jacobian[0](index(6),index(1)) = -cos_angle(1);
    
    stiffness_jacobian[0](index(6),index(7)) = internal_force(1) * sin_angle(1);
    
    stiffness_jacobian[0](index(7),index(2)) = 1.0;
    
    stiffness_jacobian[1](index(6),index(0)) = axial_stiffness * cos_angle(1);
    
    stiffness_jacobian[1](index(6),index(2)) = sin_angle(1);
    
    stiffness_jacobian[1](index(6),index(7)) = (-internal_force(0) * sin_angle(1)
                                                + internal_force(2) * cos_angle(1));
    
    stiffness_jacobian[2](index(6),index(1)) = -sin_angle(1);
    
    stiffness_jacobian[2](index(6),index(7)) = (-internal_force(1) * cos_angle(1));
    
    stiffness_jacobian[2](index(7),index(0)) = -axial_stiffness;
    
    stiffness_jacobian[3](index(6),index(4)) = -cos_angle(1);
    
    stiffness_jacobian[3](index(6),index(7)) =  (velocity(1) * sin_angle(1));
    
    stiffness_jacobian[3](index(7),index(5)) =  1.0;
    
    stiffness_jacobian[4](index(6),index(3)) =  cos_angle(1);
    
    stiffness_jacobian[4](index(6),index(5)) =  sin_angle(1);
    
    stiffness_jacobian[4](index(6),index(7)) = (-velocity(0) * sin_angle(1)
                                                + velocity(2) * cos_angle(1));
    
    stiffness_jacobian[5](index(6),index(4)) = -sin_angle(1);
    
    stiffness_jacobian[5](index(6),index(7)) = (-velocity(1) * cos_angle(1));
    
    stiffness_jacobian[5](index(7),index(3)) = -1.0;
    
    stiffness_jacobian[6](index(6),index(7)) = -sin_angle(1);
    
    
    /// -# Update the Jacobian matrix related to nodal force.
    double dur2dfi = -(-current_velocity(0) * sin_angle(0) * cos_angle(1)
                       + current_velocity(1) * cos_angle(0) * cos_angle(1));
    
    double dur2dst = -(-current_velocity(0) * cos_angle(0) * sin_angle(1)
                       - current_velocity(1) * sin_angle(0) * sin_angle(1)
                       - current_velocity(2) * cos_angle(1));
    
    double dvr2dfi = -(-current_velocity(0) * cos_angle(0)
                       -current_velocity(1) * sin_angle(0));
    
    double dvr2dst =  0.0;
    
    double dwr2dfi = -(-current_velocity(0) * sin_angle(0) * sin_angle(1)
                       + current_velocity(1) * cos_angle(0) * sin_angle(1));
    
    double dwr2dst = -( current_velocity(0) * cos_angle(0) * cos_angle(1)
                       + current_velocity(1) * sin_angle(0) * cos_angle(1)
                       - current_velocity(2) * sin_angle(1));
    
    force_jacobian(0,index(0)) = (-0.25 * M_PI * drag_force_constant(0)
                                  * relative_velocity(0)
                                  * fabs(relative_velocity(0)) / sqrt_strain);
    
    force_jacobian(0,index(3)) = (-M_PI * drag_force_constant(0)
                                  * sqrt_strain * fabs(relative_velocity(0)));
    
    force_jacobian(0,index(6)) = (unit_length_weight * sin_angle(0) * cos_angle(1)
                                  - M_PI * drag_force_constant(0)
                                  * sqrt_strain * fabs(relative_velocity(0))
                                  * dur2dfi);
    
    force_jacobian(0,index(7)) = (unit_length_weight * cos_angle(0)
                                  * sin_angle(1) - M_PI * drag_force_constant(0)
                                  * sqrt_strain * fabs(relative_velocity(0))
                                  * dur2dst);
    
    force_jacobian(1,index(0)) = (-0.25 * drag_force_constant(1)
                                  * relative_velocity(1) * sqrt_vel
                                  / sqrt_strain);
    
    force_jacobian(2,index(0)) = (-0.25 * drag_force_constant(2)
                                  * relative_velocity(2) * sqrt_vel
                                  / sqrt_strain);
    
    force_jacobian(6,index(9)) = -1.0;
    
    force_jacobian(7,index(8)) = -1.0;
    
    force_jacobian(8,index(0)) = (-internal_force(2) * 3.0
                                  * pow(1.0+strain, 2.0));
    
    force_jacobian(8,index(2)) = -pow((1.0+strain),3.0);
    
    force_jacobian(8,index(7)) = (stiffness_constant * pow(curvature(2),2.0)
                                  / pow(cos_angle(1), 2.0));
    
    force_jacobian(8,index(9)) = (stiffness_constant * tan_angle(1)
                                  * 2.0 * curvature(2));
    
    force_jacobian(9,index(0)) = (internal_force(1) * 3.0
                                  * pow(1.0+strain, 2.0));
    
    force_jacobian(9,index(1)) =  pow((1 + strain),3.0);
    
    force_jacobian(9,index(7)) = (-stiffness_constant * curvature(1)
                                  * curvature(2) / pow(cos_angle(1), 2.0));
    
    force_jacobian(9,index(8)) = (-stiffness_constant * tan_angle(1)
                                  * curvature(2));
    
    force_jacobian(9,index(9)) = (-stiffness_constant * tan_angle(1)
                                  * curvature(1));
    
    if (std::isinf(1.0/sqrt_vel))
    {
        force_jacobian(1,index(4)) =  0.0;
        force_jacobian(1,index(6)) =  (unit_length_weight * cos_angle(0));
        force_jacobian(1,index(7)) =  0.0;
        force_jacobian(2,index(5)) =  0.0;
        force_jacobian(2,index(6)) = (unit_length_weight * sin_angle(0)
                                      * sin_angle(1));
        force_jacobian(2,index(7)) = (-unit_length_weight * cos_angle(0)
                                      * cos_angle(1));
    }
    else
    {
        force_jacobian(1,index(4)) = (-0.5 * drag_force_constant(1) * sqrt_strain
                                      * (sqrt_vel + pow(relative_velocity(1),2.0)
                                         / sqrt_vel));
        
        force_jacobian(1,index(6)) = (unit_length_weight * cos_angle(0)
                                      - 0.5 * drag_force_constant(1) * sqrt_strain
                                      * (sqrt_vel * dvr2dfi
                                         + relative_velocity(1)/sqrt_vel
                                         * (relative_velocity(1) * dvr2dfi
                                            + relative_velocity(2) * dwr2dfi)));
        
        force_jacobian(1,index(7)) = (-0.5 * drag_force_constant(1) * sqrt_strain
                                      * ( sqrt_vel * dvr2dst
                                         + relative_velocity(1) / sqrt_vel
                                         * (relative_velocity(1) * dvr2dst
                                            + relative_velocity(2) * dwr2dst)));
        
        force_jacobian(2,index(5)) = (-0.5 * drag_force_constant(2) * sqrt_strain
                                      * (sqrt_vel + pow(relative_velocity(2),2.0)
                                         / sqrt_vel));
        
        force_jacobian(2,index(6)) = (unit_length_weight * sin_angle(0)
                                      * sin_angle(1)
                                      - 0.5 * drag_force_constant(2) * sqrt_strain
                                      * (sqrt_vel * dwr2dfi
                                         + relative_velocity(2) / sqrt_vel
                                         * (relative_velocity(1) * dvr2dfi
                                            + relative_velocity(2) * dwr2dfi)));
        
        force_jacobian(2,index(7)) = (-unit_length_weight * cos_angle(0)
                                      * cos_angle(1)
                                      - 0.5 * drag_force_constant(2) * sqrt_strain
                                      * (sqrt_vel * dwr2dst
                                         + relative_velocity(2) / sqrt_vel
                                         * (relative_velocity(1) * dvr2dst
                                            + relative_velocity(2) * dwr2dst)));
    }
}


////////////////////////////////////////////////////////////////////////////////
/// \brief Static analysis.
///  Add constraint for node which is the fairlead to be excited.
////////////////////////////////////////////////////////////////////////////////
void NodeEA::add_constraint(Vector3d& external_force)
{
    // Residual.
    equation(0,10) = (internal_force(0) * cos_angle(0) * cos_angle(1)
                      - internal_force(1) * sin_angle(0)
                      + internal_force(2) * cos_angle(0) * sin_angle(1)
                      - external_force(0));

    equation(1,10) = (internal_force(0) * sin_angle(0) * cos_angle(1)
                      + internal_force(1) * cos_angle(0)
                      + internal_force(2) * sin_angle(0) * sin_angle(1)
                      - external_force(1));

    equation(2,10) = (- internal_force(0) * sin_angle(1)
                      + internal_force(2) * cos_angle(1)
                      - external_force(2));

    equation(3,10) =  curvature(1);
    equation(4,10) =  curvature(2);

    // Jacobian.
    equation(0,index(0)) = axial_stiffness * cos_angle(0) * cos_angle(1);

    equation(0,index(1)) = -sin_angle(0);
    equation(0,index(2)) = (cos_angle(0) * sin_angle(1));

    equation(0,index(6)) = (- internal_force(0) * sin_angle(0) * cos_angle(1)
                            - internal_force(1) * cos_angle(0)
                            - internal_force(2) * sin_angle(0) * sin_angle(1));

    equation(0,index(7)) = (- internal_force(0) * cos_angle(0) * sin_angle(1)
                            + internal_force(2) * cos_angle(0) * cos_angle(1));

    equation(1,index(0)) = axial_stiffness * sin_angle(0) * cos_angle(1);

    equation(1,index(1)) =  cos_angle(0);
    equation(1,index(2)) =  sin_angle(0) * sin_angle(1);

    equation(1,index(6)) = (internal_force(0) * cos_angle(0) * cos_angle(1)
                            - internal_force(1) * sin_angle(0)
                            + internal_force(2) * cos_angle(0) * sin_angle(1));

    equation(1,index(7)) = (-internal_force(0) * sin_angle(0) * sin_angle(1)
                            + internal_force(2) * sin_angle(0) * cos_angle(1));

    equation(2,index(0)) = -axial_stiffness * sin_angle(1);

    equation(2,index(2)) =  cos_angle(1);

    equation(2,index(7)) = (- internal_force(0) * cos_angle(1)
                            - internal_force(2) * sin_angle(1));

    equation(3,index(8)) =  1.0;
    equation(4,index(9)) =  1.0;
}


////////////////////////////////////////////////////////////////////////////////
/// Update the effective weight to consider seabed effect.
////////////////////////////////////////////////////////////////////////////////
void NodeEA::update_effective_weight(const int is_static)
{
    // Here should use the WATER_DEPTH to determine the seabed elevation.
    if (position(0) < 0.0)
    {
        unit_length_weight = (struct_property->unit_length_weight
                              - seabed_property->stiffness_coefficient 
							  * struct_property->diameter
                              * fabs(position(0)) 
							  - seabed_property->damping_coefficient 
							  * struct_property->diameter
							  * velocity(1));
        
        // For static problem.
        if (is_static && unit_length_weight < 0)
            unit_length_weight = 0;
    }
    else
        unit_length_weight = struct_property->unit_length_weight;
}


////////////////////////////////////////////////////////////////////////////////
/// Add constraint for node which is fixed.
////////////////////////////////////////////////////////////////////////////////
void NodeEA::add_constraint()
{
    // Residual.
    equation(5,10) = velocity(0);
    equation(6,10) = velocity(1);
    equation(7,10) = velocity(2);
    equation(8,10) = curvature(1);
    equation(9,10) = curvature(2);
    
    // Jacobian.
    equation(5,index(3)) = 1.0;
    equation(6,index(4)) = 1.0;
    equation(7,index(5)) = 1.0;
    equation(8,index(8)) = 1.0;
    equation(9,index(9)) = 1.0;
}


////////////////////////////////////////////////////////////////////////////////
/// For problems when the boundary condition is given as the fairlead velocity.
////////////////////////////////////////////////////////////////////////////////
void NodeEA::add_constraint(Vector3d& forced_velocity, double time)
{
    // Residual.
    equation(0,10) = (velocity(0) * cos_angle(0) * cos_angle(1)
                      - velocity(1) * sin_angle(0)
                      + velocity(2) * cos_angle(0) * sin_angle(1)
                      - forced_velocity(0));

    equation(1,10) = (velocity(0) * sin_angle(0) * cos_angle(1)
                      + velocity(1) * cos_angle(0)
                      + velocity(2) * sin_angle(0) * sin_angle(1)
                      - forced_velocity(1));

    equation(2,10) = ((-velocity(0) * sin_angle(1)
                       + velocity(2) * cos_angle(1) - forced_velocity(2)));

    equation(3,10) =  curvature(1);
    equation(4,10) =  curvature(2);

    // Jacobian.
    equation(0,index(3)) = (cos_angle(0) * cos_angle(1));
    equation(0,index(4)) = -sin_angle(0);
    equation(0,index(5)) = (cos_angle(0) * sin_angle(1));

    equation(0,index(6)) = (-velocity(0) * sin_angle(0) * cos_angle(1)
                            - velocity(1) * cos_angle(0)
                            - velocity(2) * sin_angle(0) * sin_angle(1));

    equation(0,index(7)) = (-velocity(0) * cos_angle(0) * sin_angle(1)
                            + velocity(2) * cos_angle(0) * cos_angle(1));

    equation(1,index(3)) = (sin_angle(0) * cos_angle(1));
    equation(1,index(4)) =  cos_angle(0);
    equation(1,index(5)) = (sin_angle(0) * sin_angle(1));
    equation(1,index(6)) = (velocity(0) * cos_angle(0) * cos_angle(1)
                            - velocity(1) * sin_angle(0)
                            + velocity(2) * cos_angle(0) * sin_angle(1));

    equation(1,index(7)) = (-velocity(0) * sin_angle(0) * sin_angle(1)
                            + velocity(2) * sin_angle(0) * cos_angle(1));

    equation(2,index(3)) = -sin_angle(1);
    equation(2,index(5)) =  cos_angle(1);
    equation(2,index(7)) = (-velocity(0) * cos_angle(1)
                            -velocity(2) * sin_angle(1));

    equation(3,index(8)) =  1.0;
    equation(4,index(9)) =  1.0;
}


////////////////////////////////////////////////////////////////////////////////
/// Transform the internal force in Lagrangian coordinate to cable Cartesian
/// coordinate.
////////////////////////////////////////////////////////////////////////////////
Vector3d NodeEA::transform_internal_force(void)
{
    Vector3d transformed_force;
    transformed_force = (AngleAxisd(euler_angle(0), Vector3d::UnitZ())
                         * AngleAxisd(euler_angle(1), Vector3d::UnitY())
                         * internal_force);
    return transformed_force;
}


////////////////////////////////////////////////////////////////////////////////
/// Zero hydrostatic and hydrodynamic effects if the node is no longer
/// submerged.
////////////////////////////////////////////////////////////////////////////////
void NodeEA::update_hydrodynamics(void)
{
    if (global_position(2)<=0)
    {
        // Froude-Krylov component.
        inertial_force_constant = ((1 + hydro_property->added_mass_coefficient(1))
                                   * M_PI * pow(struct_property->diameter, 2.0)/4.0
                                   * constant->WATER_DENSITY);
        
        drag_force_constant = (constant->WATER_DENSITY * struct_property->diameter
                               * hydro_property->drag_coefficient.array());
        
        unit_length_added_mass = (constant->WATER_DENSITY * M_PI
                                  * pow(struct_property->diameter,2.0) / 4.0
                                  * hydro_property->added_mass_coefficient.array());
        
        unit_length_weight = struct_property->unit_length_weight;
    }
    else
    {
        inertial_force_constant = 0.0;
        drag_force_constant.setZero();
        unit_length_added_mass.setZero();
        unit_length_weight = (struct_property->unit_length_weight
                              + constant->WATER_DENSITY * M_PI/4.0
                              * pow(struct_property->diameter,2.0)
                              * constant->GRAV_ACC);
    }
}

////////////////////////////////////////////////////////////////////////////////
/// Current velocity update.
////////////////////////////////////////////////////////////////////////////////
void NodeEA::update_current_velocity(MatrixXd& current_velocity_poly)
{
    for (int j=0; j<3; j++)
        current_velocity(j) = poly_eval(current_velocity_poly.col(j),
                                        global_position(2));
}

} // End of namespace moor.
