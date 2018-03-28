// This file is part of OpenMOOR, an Open-source simulation program for MOORing
// systems in offshore renewable energy applications.
//
// Created by Lin Chen on Aug 28, 2017.
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


#ifndef Setting_h
#define Setting_h

#include "numeric.h"
#include "moorerror.h"
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

namespace moor {

/// Simulation type.
enum SimulationType
{
    FORCED_MOTION,    ///< Dynamic time history analysis.
    SHOOTING,         ///< Static analysis using shooting method.
    RELAXATION        ///< Static analysis using dynamic relaxation.
};


/// Required for shooting static analysis.
struct Shooting
{
    Vector3d fairlead_position_tolerance;
    double   fairlead_force_relaxation_factor;
    int      fairlead_position_iteration_limit;
    int      platform_position_iteration_limit;
    double   platform_displacement_tolerance;
    double   platform_displacement_relaxation_factor;
    Matrix6d platform_hydrostatic_stiffness;
    Vector6d platform_other_load;
    double   cable_out_of_plane_stiffness;
};
    

/// Required for dynamic relaxation.
struct Relaxation
{
    Matrix6d platform_mass;
    Matrix6d platform_hydrostatic_stiffness;
    Matrix6d platform_damping;
    Vector6d platform_other_load;
    double   platform_velocity_tolerance;
    double   cable_velocity_tolerance;
    double   stop_time;
    double   time_step;
};


/// Required for dynamic time history analysis.
struct ForcedMotion
{
    std::string data_file;
    int n_time;
    MatrixXd time_history;  // Velocity time history.
};
    
/// \brief Setting sets up Simulation.
///
/// Setting deals with file arrangement and analysis options.
class Setting
{
    
public:
    
    Setting(const std::string path,
            ForcedMotion *f_motion,
            Shooting *shoot,
            Relaxation *relax);
    
    ErrorCode validate(void);
    
    /// The full path of the compiled application OpenMOOR.exe or OpenMOOR.app.
    const std::string moor_path;
    
    /// Important: The Setting.xml should be in the same folder as the OpenMOOR
    /// program.
    const std::string setting_file;
    
    /// The folder where the Setting.xml is located and the other input files
    /// will be find in the subfolders and the output will be written into
    /// files located in the subfolders as well.
    std::string moor_folder;
    
    /// The main input file including cable informantion, reference point of the
    /// platform, structural/hydro/seabed/current properties and solvers. A
    /// relative path can be provided and the work_folder will be extraced.
    std::string mooring_file;
    
    /// Extracted from full path for the mooring input file.
    std::string work_folder;
    
    SimulationType simulation_type;
    
    /// Is saving platform states as each step.
    int is_saving_platform;
    
    Shooting *shooting;
    Relaxation *relaxation;
    ForcedMotion *forced_motion;
};

} // End of namespace moor.

#endif // setting_h
