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


#include "input.h"

namespace moor {
    
ErrorCode InputData::validate(void)
{
    // Numbers of components have also been validated to ensure that the values
    // are positive values. The consistence of number and specific definition is
    // is also validated in reading.
    
    // Validate constants data.
    if (constant->GRAV_ACC <=0 || constant->WATER_DENSITY <=0
        || constant->WATER_DEPTH <=0)
        return ErrorCode::MOORING_INVALID_CONSTANT;
    
    ////////////////////////////////////////////////////////////////////////////
    // Check cable data.
    ////////////////////////////////////////////////////////////////////////////
    int valid_current = 1, valid_solver = 1, valid_n_node = 1;
    int valid_connection = 1, valid_segment = 1;
    int valid_end = 1, valid_increase = 1, valid_state = 1;
    for (int i=0; i<cable_maps->size(); i++)
    {
        valid_current
        = (valid_current
           && (cable_maps->at(i).i_current >= 0
               && cable_maps->at(i).i_current < currents->size()));
        
        valid_solver
        = (valid_solver
           && (cable_maps->at(i).i_solver >= 0
               && cable_maps->at(i).i_solver < solvers->size()));
        
        valid_n_node = (valid_n_node && cable_maps->at(i).n_node > 0);
        
        valid_segment
        = (valid_segment
           && cable_maps->at(i).segment_length.rows() > 0
           && cable_maps->at(i).segment_length.rows() == cable_maps->at(i).i_struct_property.rows()
           && cable_maps->at(i).segment_length.rows() == cable_maps->at(i).i_hydro_property.rows()
           && cable_maps->at(i).segment_length.rows() == cable_maps->at(i).i_seabed_property.rows());
        
        valid_end
        = (valid_end
           && cable_maps->at(i).i_connection(0) != cable_maps->at(i).i_connection(1)
           && cable_maps->at(i).i_connection(0) >= 0
           && cable_maps->at(i).i_connection(0) < connections->size()
           && cable_maps->at(i).i_connection(1) >= 0
           && cable_maps->at(i).i_connection(1) < connections->size()
           && (connections->at(cable_maps->at(i).i_connection(0)).type == ConnectionType::ANCHOR)
           && (connections->at(cable_maps->at(i).i_connection(1)).type == ConnectionType::FAIRLEAD));
        
        for (int j=0; j < cable_maps->at(i).segment_length.rows(); j++)
            valid_increase = (valid_increase && cable_maps->at(i).segment_length(j) > 0);
        
        if (cable_maps->at(i).initial_state.rows() > 2)
        {
            valid_state
            = (valid_state
               && cable_maps->at(i).initial_state(0,0) < 1E-3
               && cable_maps->at(i).initial_state(0,0) >= 0
               && abs(cable_maps->at(i).initial_state.bottomRows<1>()(0,0)
                      - cable_maps->at(i).segment_length.sum()) < 1e-3);
            
            int valid = 1;
            for (int j=1; j < cable_maps->at(i).initial_state.rows(); j++)
            {
                valid
                = (valid && cable_maps->at(i).initial_state(j,0) > 0
                   && (cable_maps->at(i).initial_state(j,0) - cable_maps->at(i).initial_state(j-1,0)) > 0 );
            }
            valid_state = valid_state && valid;
        }
    }
    
    if (!valid_current) return ErrorCode::MOORING_INVALID_CABLE_CURRENT_INDEX;
    if (!valid_solver) return ErrorCode::MOORING_INVALID_CABLE_SOLVER_INDEX;
    if (!valid_n_node) return ErrorCode::MOORING_INVALID_CABLE_NODE_NUM;
    if (!valid_segment) return ErrorCode::MOORING_INVALID_CABLE_PROP_ASSOCIATION;
    if (!valid_end) return ErrorCode::MOORING_INVALID_CABLE_CONNECTION_INDEX;
    if (!valid_increase) return ErrorCode::MOORING_INVALID_CABLE_SEGMENT_LENGTH;
    if (!valid_state) return ErrorCode::MOORING_INVALID_CABLE_STATE;
    

    ////////////////////////////////////////////////////////////////////////////
    // Check structproperty data.
    ////////////////////////////////////////////////////////////////////////////
    int valid_positive = 1, valid_nonnegative = 1;
    for (int i=0; i<struct_properties->size(); i++)
    {
        valid_positive = (valid_positive
                         && struct_properties->at(i).diameter > 0
                         && struct_properties->at(i).unit_length_mass > 0
                         && struct_properties->at(i).unit_length_weight > 0
                         && struct_properties->at(i).axial_stiffness > 0);
        valid_nonnegative = (valid_nonnegative
                            && struct_properties->at(i).bending_stiffness >=0
                            && struct_properties->at(i).torsional_stiffness >=0
                            && struct_properties->at(i).damping_coefficient >=0);
    }
    
    if (!valid_positive) return ErrorCode::MOORING_INVALID_STRUCT_PROP;
    if (!valid_nonnegative) return ErrorCode::MOORING_INVALID_STRUCT_PROP;
    
    ////////////////////////////////////////////////////////////////////////////
    // Check hydroproperty data.
    ////////////////////////////////////////////////////////////////////////////
    valid_nonnegative = 1;
    for (int i=0; i<hydro_properties->size(); i++)
    {
        for (int j=0; j<3; j++)
            valid_nonnegative = (valid_nonnegative
                                && hydro_properties->at(i).added_mass_coefficient(j) >=0
                                && hydro_properties->at(i).drag_coefficient(j) >=0);
    }
    
    if (!valid_nonnegative) return ErrorCode::MOORING_INVALID_HYDRO_COEFFICIENT;
    
    ////////////////////////////////////////////////////////////////////////////
    // Check seabedproperty data.
    ////////////////////////////////////////////////////////////////////////////
    valid_nonnegative = 1;
    for (int i=0; i<seabed_properties->size(); i++)
    {
        valid_nonnegative = (valid_nonnegative
                            && seabed_properties->at(i).damping_coefficient >=0
                            && seabed_properties->at(i).stiffness_coefficient >=0);
    }
    
    if (!valid_nonnegative) return ErrorCode::MOORING_INVALID_SEABED;
    
    ////////////////////////////////////////////////////////////////////////////
    // Check current data.
    ////////////////////////////////////////////////////////////////////////////
    int valid= 1;
    for (int i=0; i<currents->size(); i++)
    {
        valid = (valid && currents->at(i).polyfit_order>=0
                 && currents->at(i).profile_data.rows() >= 1
                 && currents->at(i).polyfit_order < currents->at(i).profile_data.rows());
        
        for (int j=0; j<currents->at(i).profile_data.rows(); j++)
        {
            valid = (valid && currents->at(i).profile_data(j, 2) <= 0);
            if (j > 0)
            {
                valid = valid && (currents->at(i).profile_data(j, 2)
                                  < currents->at(i).profile_data(j-1,2));
            }
        }
    }
    
    if (!valid) return ErrorCode::MOORING_INVALID_CURRENT;
    
    ////////////////////////////////////////////////////////////////////////////
    // Validate solver data.
    ////////////////////////////////////////////////////////////////////////////
    for (int i=0; i<solvers->size(); i++)
    {
        valid = (valid && solvers->at(i).n_iteration_limit>0
                 && solvers->at(i).convergence_tolerance >0
                 && solvers->at(i).initial_relaxation_factor >0
                 && solvers->at(i).increment_factor >=1.0
                 && solvers->at(i).decrement_factor >=1.0
                 && abs(solvers->at(i).lambda_infinity - 1) > 1E-6);
    }
    
    if (!valid) return ErrorCode::MOORING_INVALID_SOLVER;
    
    return ErrorCode::SUCCESS;
}
    
} // End of namespace moor.
