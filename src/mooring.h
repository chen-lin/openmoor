// This file is part of OpenMOOR, an Open-source simulation program for MOORing
// systems in offshore renewable energy applications.
//
// Created by Lin Chen on May 24, 2017.
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


#ifndef mooring_h
#define mooring_h

#include "numeric.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

namespace moor {

////////////////////////////////////////////////////////////////////////////////
/// Constants used in mooring system analysis.
////////////////////////////////////////////////////////////////////////////////
struct Constant
{
    double WATER_DENSITY;
    double GRAV_ACC;
    double WATER_DEPTH;
};
    
    
////////////////////////////////////////////////////////////////////////////////
/// \brief Rigid body movement.
///
/// Rigid body motion described by a reference point of 6 degrees of freedom.
////////////////////////////////////////////////////////////////////////////////
struct RigidBody
{
    Vector6d position;
};


////////////////////////////////////////////////////////////////////////////////
/// \brief Cable generator.
///
/// Cable length and mapping properties for cable with varied properties and
/// also options in solving cable. In case of both n_node and initial_state_file
/// are provided. The number of nodes as given in initial_state_file should be
/// consistent with n_node.
////////////////////////////////////////////////////////////////////////////////
struct CableMap
{
    int         i_current;
    int         i_solver;
    Vector2i    i_connection;
    // Cable divided into segments according to varied properties.
    VectorXd    segment_length;
    VectorXi    i_struct_property;
    VectorXi    i_hydro_property;
    VectorXi    i_seabed_property;
    int         n_node;
    int         is_saving;
    std::string initial_state_file;
    MatrixXd    initial_state;
};

    
////////////////////////////////////////////////////////////////////////////////
/// \brief Structural property.
///
/// Structural property is to be associated to each node to model cable with
/// varied properties along the cable coordinate.
////////////////////////////////////////////////////////////////////////////////
struct StructProperty
{
    double diameter;           ///< %Cable outer diameter, in meter.
    double unit_length_mass;   ///< %Cable mass per unit length, in kg/m.
    /// %Cable weight per unit length: wet weight in N/m.
    double unit_length_weight;
    double axial_stiffness;    ///< Axial stiffness in N*m.
    double bending_stiffness;  ///< Bending stiffness.
    double torsional_stiffness;///< Torsional stiffness.
    double damping_coefficient;///< %Cable internal damping.
};

    
////////////////////////////////////////////////////////////////////////////////
/// Types of connections.
////////////////////////////////////////////////////////////////////////////////
enum ConnectionType
{
    FAIRLEAD,        ///< Fairlead attached to a platform.
    ANCHOR,          ///< Fixed connection.
    INTERCONNECTION  ///< Interconnection which is presently used.
};


////////////////////////////////////////////////////////////////////////////////
/// Connection with a particular type and position. For fairlead the position is
/// initial position.
////////////////////////////////////////////////////////////////////////////////
struct Connection
{
    ConnectionType type;
    Vector3d position;
};


////////////////////////////////////////////////////////////////////////////////
/// \brief Hydrodynamic property.
///
/// Hydro-property is to be associated with particular cable nodes.
////////////////////////////////////////////////////////////////////////////////
struct HydroProperty
{
    /// The added mass coefficients in tangential, normal and binormal direction
    /// of the cable. The tangential coefficient is included here to consider
    /// the case when a mooring chain is equivalently modeled as a circular
    /// cable.
    Vector3d added_mass_coefficient;
    
    /// Drag coefficients in tangent, normal and binormal directions for using
    /// Morison's formula, i.e. \f$C_{dt},~~C_{dn},~~C_{db}\f$.
    Vector3d drag_coefficient;
};

////////////////////////////////////////////////////////////////////////////////
/// \brief Seabed property.
///
/// Seabed property is associated with particular cable nodes. Different cables
/// and nodes can be associated with different seabed conditions.
////////////////////////////////////////////////////////////////////////////////
struct SeabedProperty
{
    double damping_coefficient;  //!< Not used currently.
    double stiffness_coefficient;//!< Nondimensional parameter.
};


////////////////////////////////////////////////////////////////////////////////
/// \brief Current class.
///
/// Fluid current velocity profile defined using several discrete data points.
/// The data points will be fitted using polynomials and used for cable
/// solution.
////////////////////////////////////////////////////////////////////////////////
struct Current
{
    /// For fitting the discrete data points. Currently the velocity is assumed
    /// the function of the vertical position only.
    int polyfit_order;

    /// File saving data points with the first line the header and the data
    /// starts from the second line.
    std::string profile_file;
    
    /// Current profile data including the current velocity and measurement
    /// points in global reference system. Coordinate in global referential
    /// frame, i.e. \f$x, y, z\f$. However, here the current profile is only
    /// considered to be varying in vertical direction, i.e. along \f$z\f$
    /// direction. The other two coordinates are actually not used presently.
    MatrixXd profile_data;
};
    
} // End of namespace moor.

#endif // mooring_h
