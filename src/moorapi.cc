// This file is part of OpenMOOR, an Open-source simulation program for MOORing
// systems in offshore renewable energy applications.
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


#include "moorapi.h"
#include "numeric.h"
#include "input.h"
#include "reader.h"
#include "writer.h"
#include "solver.h"
#include "utility.h"
#include "platform.h"
#include "mooring.h"
#include "moorerror.h"
#include "cable.h"
#include "node.h"
#include <string>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;
using namespace moor;

////////////////////////////////////////////////////////////////////////////////
/// Create all global objects to be used.
////////////////////////////////////////////////////////////////////////////////
ErrorCode err_code = ErrorCode::SUCCESS;
ErrorCategory error_domain;
Reader reader;
Writer writer;
Platform platform;
Constant constant;
vector<StructProperty> struct_properties;
vector<HydroProperty> hydro_properties;
vector<SeabedProperty> seabed_properties;
vector<CableMap> cable_maps;
vector<Solver> solvers;
vector<Cable> cables;
string work_folder;


////////////////////////////////////////////////////////////////////////////////
/// Initialization.
////////////////////////////////////////////////////////////////////////////////
int DECLDIR initialize(char file_name[])
{
    string input_file_name(file_name);
    size_t work_folder_index = input_file_name.find_last_of("/\\");
    work_folder = input_file_name.substr(0, work_folder_index + 1);
    
    if (struct_properties.size() > 0) struct_properties.clear();
    if (hydro_properties.size() > 0) hydro_properties.clear();
    if (seabed_properties.size() > 0) seabed_properties.clear();
    if (cable_maps.size() > 0) cable_maps.clear();
    if (solvers.size() > 0) solvers.clear();
    if (cables.size() > 0) cables.clear();
    
    /// The following variables are used only for initialization.
    RigidBody reference_point;
    vector<Connection> connections;
    vector<Current> currents;
    
    InputData input_data(&constant,
                         &reference_point,
                         &cable_maps,
                         &connections,
                         &struct_properties,
                         &hydro_properties,
                         &seabed_properties,
                         &currents,
                         &solvers);
    
    err_code = reader.read_to(input_data, input_file_name);
    
    if (err_code != ErrorCode::SUCCESS)
    {
        writer.write_error_message(work_folder+"openmoor.log",
                                   err_code, error_domain);
        return 1;
    }
    
    err_code = input_data.validate();
    
    if (err_code != ErrorCode::SUCCESS)
    {
        writer.write_error_message(work_folder+"openmoor.log",
                                   err_code, error_domain);
        return 1;
    }
    
    for (int i=0; i < cable_maps.size(); i++)
    {
        Cable cable(cable_maps[i].segment_length.sum(),
                    connections[cable_maps[i].i_connection(0)].position,
                    connections[cable_maps[i].i_connection(1)].position);
        
        cables.push_back(cable);
    }
    
    int n_fairlead = cables.size(); // Current requirement.
    Matrix3Xd fairlead_position;
    fairlead_position.resize(3, n_fairlead);
    for (int i=0; i<n_fairlead; i++)
    {
        fairlead_position.col(i)
        = connections[cable_maps[i].i_connection(1)].position;
    }
    
    for (int i=0; i<cables.size(); i++)
    {
        cables[i].initialize(cable_maps[i],
                             currents[cable_maps[i].i_current],
                             solvers[cable_maps[i].i_solver],
                             struct_properties,
                             hydro_properties,
                             seabed_properties,
                             constant);
    }
    
    /// Initialize platform position and fairlead number and positions.
    platform.initialize(reference_point.position,
                        n_fairlead,
                        fairlead_position,
                        cables);
    
    string cable_file_name;
    // Output cable states if desired.
    for (int i=0; i<cables.size(); i++)
    {
        if (cables[i].is_saving)
        {
            cable_file_name = (work_folder + "cable_" + to_string(i)
                               + "_time_" + to_string(0.0) + ".dat");
            writer.write(cable_file_name, cables[i]);
        }
    }

    // Initialize solvers before assign the solver pointer to cables.
    for (int i=0; i<solvers.size(); i++)
        solvers[i].initialize(10, 5); // For using Euler Angle formulation.
    
    err_code = ErrorCode::SUCCESS;
    return 0;
}
    

////////////////////////////////////////////////////////////////////////////////
/// Time stepping.
////////////////////////////////////////////////////////////////////////////////
int DECLDIR update(double* mooring_load, const double displacement[],
                   const double velocity[], const double time,
                   const double dt, const double save_dt)
{
    Vector6d forced_displacement, forced_velocity;
    
    for (int i=0; i<6; i++)
    {
        forced_displacement(i) = displacement[i];
        forced_velocity(i) = velocity[i];
    }
    
    if (dt < 1E-6)
    {
        err_code = ErrorCode::TIME_STEP_TOO_SMALL;
    }
    else
    {
        err_code =
        platform.solve_mooring_load(forced_displacement,
                                    forced_velocity,
                                    time,
                                    dt,
                                    cables);
    }
    
    for (int i=0; i<6; i++) mooring_load[i] = platform.get_mooring_load(i);
    
    if (err_code != ErrorCode::SUCCESS
        && err_code != ErrorCode::TIME_STEP_TOO_SMALL)
    {
        writer.write_error_message(work_folder+"openmoor.log",
                                   err_code, error_domain);
        return 1;
    }
    
    // Output cable states if desired.
    double s_dt = save_dt > dt ? save_dt : dt;
    string cable_file_name;
    for (int i=0; i < cables.size(); i++)
    {
        if (cables[i].is_saving && fabs(remainder(time, s_dt)) < 1E-4)
        {
            
            cable_file_name = (work_folder + "cable_" + to_string(i)
                               + "_time_" + to_string(time) + ".dat");
            writer.write(cable_file_name, cables[i]);
        }
    }
    
    return 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Finish and clear all objects.
////////////////////////////////////////////////////////////////////////////////
int DECLDIR finish(void)
{
    if (struct_properties.size() > 0) struct_properties.clear();
    if (hydro_properties.size() > 0) hydro_properties.clear();
    if (seabed_properties.size() > 0) seabed_properties.clear();
    if (cable_maps.size() > 0) cable_maps.clear();
    if (solvers.size() > 0) solvers.clear();
    if (cables.size() > 0) cables.clear();
    return 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Get mooring load component.
////////////////////////////////////////////////////////////////////////////////
double DECLDIR get_mooring_load(const int i_component)
{
    if (i_component >=0 && i_component <6)
        return platform.get_mooring_load(i_component);
    else
        return NAN;
}


////////////////////////////////////////////////////////////////////////////////
/// Check reference point velocity.
////////////////////////////////////////////////////////////////////////////////
double DECLDIR get_reference_point_displacement(const int i_component)
{
    if (i_component >=0 && i_component <6)
        return platform.get_displacement(i_component);
    else
        return NAN;
}


////////////////////////////////////////////////////////////////////////////////
/// Check reference point velocity.
////////////////////////////////////////////////////////////////////////////////
double DECLDIR get_reference_point_velocity(const int i_component)
{
    if (i_component >=0 && i_component <6)
        return platform.get_velocity(i_component);
    else
        return NAN;
}


////////////////////////////////////////////////////////////////////////////////
/// Get cable total force at fairlead.
////////////////////////////////////////////////////////////////////////////////
double DECLDIR get_cable_fairlead_force(const int i_cable)
{
	if (i_cable >= 0 && i_cable < cables.size())
		return cables[i_cable].get_total_fairlead_force();
	else
		return NAN;
}


////////////////////////////////////////////////////////////////////////////////
/// Get cable total force at anchor.
////////////////////////////////////////////////////////////////////////////////
double DECLDIR get_cable_anchor_force(const int i_cable)
{
    if (i_cable >= 0 && i_cable < cables.size())
        return cables[i_cable].get_total_anchor_force();
    else
        return NAN;
}


////////////////////////////////////////////////////////////////////////////////
/// Get total nodal force at a node of a cable.
////////////////////////////////////////////////////////////////////////////////
double DECLDIR get_cable_nodal_force(const int i_cable, const int i_node)
{
    if (i_cable >= 0 && i_cable < cables.size())
        return cables[i_cable].get_total_nodal_force(i_node);
    else
        return NAN;
}
