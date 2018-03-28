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


#include "setting.h"

using namespace std;

namespace moor {
    
Setting::Setting(const std::string path,
                 ForcedMotion *f_motion,
                 Shooting *shoot,
                 Relaxation *relax) :
    moor_path(path),
    setting_file("Setting.xml")
{
    size_t folder_index = moor_path.find_last_of("/\\");
    moor_folder = moor_path.substr(0, folder_index + 1);
    
    shooting = shoot;
    relaxation = relax;
    forced_motion = f_motion;
    
    // The only default value.
    is_saving_platform = 1;
}


ErrorCode Setting::validate(void)
{
    switch(simulation_type)
    {
        case SimulationType::SHOOTING:
        {
            if ( shooting->fairlead_position_tolerance(0) <= 0
                || shooting->fairlead_position_tolerance(1) <= 0
                || shooting->fairlead_position_tolerance(2) <= 0
                || shooting->fairlead_force_relaxation_factor <= 0
                || shooting->fairlead_position_iteration_limit < 0
                || shooting->platform_position_iteration_limit < 0)
                return ErrorCode::SHOOTING_BAD_PARA;
        }
            break;
            
        case SimulationType::RELAXATION:
        {
            if ( relaxation->platform_velocity_tolerance <= 0
                || relaxation->cable_velocity_tolerance <= 0
                || relaxation->stop_time < 0
                || relaxation->time_step <= 0 )
                return ErrorCode::RELAXATION_BAD_PARA;
        }
            break;
            
        case SimulationType::FORCED_MOTION:
        {
            if (forced_motion->time_history(0,0) < 0)
                return ErrorCode::FORCED_MOTION_BAD_TIME;
            else
            {
                for (int i=1; i<forced_motion->n_time; i++)
                    if (forced_motion->time_history(i,0)
                        < forced_motion->time_history(i-1,0))
                        return ErrorCode::FORCED_MOTION_BAD_TIME;
            }
        }
            break;
    }
    return ErrorCode::SUCCESS;
}

} // End of namespace moor.
