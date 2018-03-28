// This file is part of OpenMOOR, an Open-source simulation program for MOORing
// systems in offshore renewable energy applications.
//
// Created by Lin Chen on Aug 27, 2017.
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


#ifndef input_h
#define input_h

#include "numeric.h"
#include "moorerror.h"
#include "solver.h"
#include "mooring.h"
#include <vector>
#include <iostream>
#include <string>

namespace moor {

/// \brief Input class storing input data.
///
/// The input data is extracted from the input file and saved as an object for
/// initializing the platform, cables, structural/hydro/seabd/current and other
/// objects for subsequent analysis. The parameters are described as follows.
///  - Platform: reference point position and the number of fairleads.
///  - Fairleads: initial positions.
///  - Anchors: position.
///  - Note that for the current time being the number of the fairleads is equal
///    to the number of the anchors.
///  - Parameters of cables: number of cables; node number for each cable;
///    length of cables; cable solver index; initial cable state file name.
///  - Structural, hydro and seabed properties.
///  - Current properties: number of data points and the current file name.
///  - Solver parameters.
class InputData
{
    
public:
    /// Input.
    InputData(Constant *cons,
              RigidBody *ref_p,
              std::vector<CableMap> *c_map,
              std::vector<Connection> *conn,
              std::vector<StructProperty> *stp,
              std::vector<HydroProperty> *hyp,
              std::vector<SeabedProperty> *sbp,
              std::vector<Current> *crrt,
              std::vector<Solver> *slvr) :
        cable_maps(c_map),
        reference_point(ref_p),
        constant(cons),
        connections(conn),
        struct_properties(stp),
        hydro_properties(hyp),
        seabed_properties(sbp),
        solvers(slvr),
        currents(crrt)
    {
    };
    
    InputData(void)
    {
        cable_maps=NULL;
        reference_point=NULL;
        constant=NULL;
        connections=NULL;
        struct_properties=NULL;
        hydro_properties=NULL;
        seabed_properties=NULL;
        solvers=NULL;
        currents=NULL;
    };
    
    /// Check the input data and return an error code.
    ErrorCode validate(void);
    
public:
    
    /// Constants.
    Constant *constant;
    
    /// Reference point of the mooring system for calculating mooring force and
    /// computing fairlead velocity and displacement due to platform motion.
    RigidBody *reference_point;
    
    /// The number of structural/hydro/seabed properties, current types and the
    /// types of solvers used for solving cable states. Each cable with one
    /// solver and one reference current but could have varied structural/hydro/
    /// seabed properties for varied nodes.
    int n_struct_property;
    int n_hydro_property;
    int n_seabed_property;
    int n_current;
    int n_solver;
    
    /// Component numbers: one platform; multiple cables and the corresponding
    /// connections.
    int n_connection;
    int n_cable;

    /// Cable geometric parameters and also a map of cable parameters with
    /// corresponding properties, current and solver.
    std::vector<CableMap>       *cable_maps;
    
    /// Properties and options for solving cable motion.
    std::vector<Current>        *currents;
    std::vector<StructProperty> *struct_properties;
    std::vector<HydroProperty>  *hydro_properties;
    std::vector<SeabedProperty> *seabed_properties;
    std::vector<Solver>         *solvers;
    std::vector<Connection>     *connections;
};

} // End of namespace moor.

#endif // input_h
