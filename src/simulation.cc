// This file is part of OpenMOOR, an Open-source simulation program for MOORing
// systems in offshore renewable energy applications.
//
// Created by Lin Chen on Jul 26, 2017.
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


#include "simulation.h"
using namespace std;

namespace moor {

////////////////////////////////////////////////////////////////////////////////
/// Run simulation according to the simulation type.
////////////////////////////////////////////////////////////////////////////////
ErrorCode Simulation::run(void)
{
    // Perform analysis.
    cout << endl;
    
    ErrorCode err_code;
    switch(setting->simulation_type)
    {
        case SimulationType::SHOOTING:
        {
            cout << "--- STATIC ANSLYSIS USING SHOOTING METHOD" << endl;
            err_code = shoot(*setting->shooting);
        }
            break;

        case SimulationType::FORCED_MOTION:
        {
            cout << "--- FORCED MOTION ANALYSIS" << endl;
            err_code = analyze_forced_motion(*setting->forced_motion);
        }
            break;

        case SimulationType::RELAXATION:
        {
            cout << "--- STATIC ANSLYSIS USING DYNAMIC RELAXATION" << endl;
            err_code = relax(*setting->relaxation);
        }
            break;
    }
    return ErrorCode::SUCCESS;
}
    

////////////////////////////////////////////////////////////////////////////////
/// A two-level shooting method for static analysis.
////////////////////////////////////////////////////////////////////////////////
ErrorCode Simulation::shoot(Shooting &shooting)
{
    // Set parameters for static analysis using shooting method.
    Vector6d displacement;
    displacement.setZero();

    int i_shoot = 1, crit_shoot  = 1;
    string platform_file_name = setting->work_folder+"s_platform_state.dat";

    // Write the header of the log file.
    if (setting->is_saving_platform)
    {
        writer->write(platform_file_name, *platform);
        writer->write(platform_file_name, *platform, i_shoot);
    }

    Matrix6d total_stiffness_matrix;
    Vector6d force_difference, displacement_increment;

    // Modify the cable out-of-plane stiffness.
    for (int i=0; i < platform->get_n_fairlead(); i++)
        cables->at(i).linearized_stiffness(2,2)
        = shooting.cable_out_of_plane_stiffness;

    ErrorCode err_code;
    while (crit_shoot &&
           i_shoot < shooting.platform_position_iteration_limit)
    {
        cout << endl << "  -- platform balance shooting step: "
        << i_shoot << endl;
        
        err_code =
        platform->solve_mooring_load(displacement,
                                     shooting.fairlead_force_relaxation_factor,
                                     shooting.fairlead_position_tolerance,
                                     shooting.fairlead_position_iteration_limit,
                                     *cables);
        
//        if (err_code != ErrorCode::SUCCESS) return err_code;

        // This needs to be performed after displacement movement.
        platform->approximate_mooring_stiffness(*cables);

        for (int i=0; i<platform->get_n_fairlead(); i++)
        {
            
            // Output computation results.
            if (cables->at(i).is_saving)
            {
                string cable_file_name = (setting->work_folder+"s_cable_"
                                          + std::to_string(i)+"_iteration_"
                                          + std::to_string(i_shoot) + ".dat");

                // Write cable file.
                writer->write(cable_file_name, cables->at(i));
            }
        }

        force_difference = (platform->get_mooring_load()
                            + shooting.platform_other_load
                            - shooting.platform_hydrostatic_stiffness
                            * displacement);

        total_stiffness_matrix
        = (platform->get_approximated_mooring_stiffness()
           + shooting.platform_hydrostatic_stiffness);
        
        displacement_increment
        = (shooting.platform_displacement_relaxation_factor
           * total_stiffness_matrix.colPivHouseholderQr().solve(force_difference));

        displacement  += displacement_increment;

        crit_shoot
        = (displacement_increment.array().abs().maxCoeff()
           > shooting.platform_displacement_tolerance);

        writer->write(platform_file_name, *platform, i_shoot);

        i_shoot++;
    }
    return ErrorCode::SUCCESS;
}


////////////////////////////////////////////////////////////////////////////////
/// The dynamic analysis of the mooring system for given excitation time
/// history, following the steps:
////////////////////////////////////////////////////////////////////////////////
ErrorCode Simulation::analyze_forced_motion(ForcedMotion &forced_motion)
{
    VectorXd times;
    times.resize(forced_motion.n_time);
    times = forced_motion.time_history.col(0);

    // - Perform dynamic analysis now.
    int i_time = 0;
    double time_step;

    // Write the initial solution to file.
    string platform_file_name = (setting->work_folder + "f_platform_state.dat");

    if (setting->is_saving_platform)
    {
        writer->write(platform_file_name, *platform); // Header.
        writer->write(platform_file_name, *platform, times(i_time));
    }

    // Print on screen.
    cout << endl;
    cout << "  -- t = " << setw(8) << times(i_time) << " s : " << endl;
    cout << "     | initialization data saved " << endl << endl;

    for (int i=0; i < platform->get_n_fairlead(); i++)
    {
        if (cables->at(i).is_saving)
        {
            // Output cable computation results
            string cable_file_name = (setting->work_folder
                                      + "f_cable_" + std::to_string(i)
                                      + "_time_" +
                                      std::to_string(times(i_time))
                                      + ".dat");

            // Write cable state to file.
            writer->write(cable_file_name, cables->at(i));
        }
    }

    ErrorCode err_code;
    
    Vector6d displacement, velocity;
    displacement.setZero();
    /// - Perform dynamic analysis.
    for (int i_time = 1; i_time < forced_motion.n_time; i_time++)
    {
        cout << "  -- t = " << setw(8) << times(i_time) << " s: " << endl;

        velocity = forced_motion.time_history.block(i_time, 1, 1,6).transpose();
        time_step = times(i_time) - times(i_time - 1);

        displacement = displacement + velocity * time_step;

        err_code =
        platform->solve_mooring_load(displacement, velocity,
                                     times(i_time),
                                     time_step,
                                     *cables);
        
        if (err_code != ErrorCode::SUCCESS) return err_code;
        
        for (int i=0; i < platform->get_n_fairlead(); i++)
        {
            // Screen
            cout << "     | cable " << setw(3) << i << " solved with " << setw(3)
            << cables->at(i).solver->get_n_iteration() << " iterations" << endl;
            
            if (cables->at(i).is_saving)
            {
                // Output cable states if desired.
                string cable_file_name = setting->work_folder +
                "f_cable_" + std::to_string(i) + "_time_" +
                std::to_string(times(i_time)) + ".dat";
                
                writer->write(cable_file_name, cables->at(i));
            }
        }
        cout << endl;
        
        /// - Write platform file.
        if (setting->is_saving_platform)
            writer->write(platform_file_name, *platform, times(i_time));
    }
    return ErrorCode::SUCCESS;
}


////////////////////////////////////////////////////////////////////////////////
/// Carrying out dynamic relaxation analysis using OpenMOOR needs to supply
/// the platform structural property.
////////////////////////////////////////////////////////////////////////////////
ErrorCode Simulation::relax(Relaxation &relaxation)
{
    /// - Set platform matrix and force vector.
    platform->set_structural_matrix(relaxation.platform_mass,
                                    relaxation.platform_damping,
                                    relaxation.platform_hydrostatic_stiffness,
                                    relaxation.platform_other_load);
    
    double t = 0, dt = relaxation.time_step;

    cout << endl;

    string platform_file_name = (setting->work_folder + "r_platform_state.dat");

    // Write header.
    if (setting->is_saving_platform)
    {
        writer->write(platform_file_name, *platform);
        writer->write(platform_file_name, *platform, t);
    }
    
    Vector6d displacement, velocity;

    ErrorCode err_code;
    
    /// - Perform dynamic analysis of the coupled system.
    while (t < relaxation.stop_time)
    {
        platform->update_state(t, dt);

        cout << "  -- t = " << setw(8) << t << " s: " << endl;

        displacement = platform->get_displacement();
        velocity = platform->get_velocity();

        err_code =
        platform->solve_mooring_load(displacement,
                                     velocity,
                                     t,
                                     dt,
                                     *cables);
        if (err_code != ErrorCode::SUCCESS) return err_code;
        
        for (int i=0; i < platform->get_n_fairlead(); i++)
        {
            // Screen print.
            cout << "     | cable " << setw(3) << i << " solved with " << setw(3)
            << cables->at(i).solver->get_n_iteration() << " iterations" << endl;
            
            if (cables->at(i).is_saving)
            {
                // Output cable states if desired.
                string cable_file_name = setting->work_folder +
                "r_cable_" + std::to_string(i) + "_time_" +
                std::to_string(t) + ".dat";
                
                writer->write(cable_file_name, cables->at(i));
            }
        }
        cout << endl;
        t = t + dt;
    }
    return ErrorCode::SUCCESS;
}

} // End of namespace moor.
