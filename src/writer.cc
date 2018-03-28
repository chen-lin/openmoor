// This file is part of OpenMOOR, an Open-source simulation program for MOORing
// systems in offshore renewable energy applications.
//
// Created by Lin Chen on Aug 25, 2017.
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


#include "writer.h"
using namespace std;

namespace moor {

////////////////////////////////////////////////////////////////////////////////
/// Write platform file header.
////////////////////////////////////////////////////////////////////////////////
void Writer::write(const string target_file_name, const Platform& platform)
{
    // Write header into file and check the existence of the target file.
    ofstream target_file;
    target_file.open(target_file_name, std::ofstream::out);
    
    // Write the header.
    target_file <<
    setw(12) << "step/time(s)" <<
    setw(14) << "x (m)" <<
    setw(14) << "y (m)" <<
    setw(14) << "z (m)" <<
    setw(14) << "roll (rad)" <<
    setw(14) << "pitch (rad)" <<
    setw(14) << "yaw (rad)" <<
    setw(14) << "v_x (m/s)" <<
    setw(14) << "v_y (m/s)" <<
    setw(14) << "v_z (m/s)" <<
    setw(14) << "v_rll (rad/s)" <<
    setw(14) << "v_pth (rad/s)" <<
    setw(14) << "v_yaw (rad/s)" <<
    setw(14) << "F_x (N)" <<
    setw(14) << "F_y (N)" <<
    setw(14) << "F_z (N)" <<
    setw(14) << "M_x (N*m)" <<
    setw(14) << "M_y (N*m)" <<
    setw(14) << "M_z (N*m)" << endl;
    
    target_file.close();
}


////////////////////////////////////////////////////////////////////////////////
/// Write platform in shooting for static analysis.
////////////////////////////////////////////////////////////////////////////////
void Writer::write(const string target_file_name, const Platform& platform,
                   const int iter_step)
{
    ofstream target_file;
    target_file.open(target_file_name, std::ofstream::app);
    
    target_file <<
    setw(12) << iter_step <<
    setw(14) << scientific << platform.get_displacement(0) <<
    setw(14) << scientific << platform.get_displacement(1) <<
    setw(14) << scientific << platform.get_displacement(2) <<
    setw(14) << scientific << platform.get_displacement(3) <<
    setw(14) << scientific << platform.get_displacement(4) <<
    setw(14) << scientific << platform.get_displacement(5) <<
    setw(14) << scientific << platform.get_velocity(0) <<
    setw(14) << scientific << platform.get_velocity(1) <<
    setw(14) << scientific << platform.get_velocity(2) <<
    setw(14) << scientific << platform.get_velocity(3) <<
    setw(14) << scientific << platform.get_velocity(4) <<
    setw(14) << scientific << platform.get_velocity(5) <<
    setw(14) << scientific << platform.get_mooring_load(0) <<
    setw(14) << scientific << platform.get_mooring_load(1) <<
    setw(14) << scientific << platform.get_mooring_load(2) <<
    setw(14) << scientific << platform.get_mooring_load(3) <<
    setw(14) << scientific << platform.get_mooring_load(4) <<
    setw(14) << scientific << platform.get_mooring_load(5) << endl;
    
    target_file.close();
}


////////////////////////////////////////////////////////////////////////////////
// Writer platform.
////////////////////////////////////////////////////////////////////////////////
void Writer::write(const string target_file_name, const Platform& platform,
                   const double time)
{
    ofstream target_file;
    target_file.open(target_file_name, std::ofstream::app);
    
    target_file <<
    setw(12) << time <<
    setw(14) << scientific << platform.get_displacement(0) <<
    setw(14) << scientific << platform.get_displacement(1) <<
    setw(14) << scientific << platform.get_displacement(2) <<
    setw(14) << scientific << platform.get_displacement(3) <<
    setw(14) << scientific << platform.get_displacement(4) <<
    setw(14) << scientific << platform.get_displacement(5) <<
    setw(14) << scientific << platform.get_velocity(0) <<
    setw(14) << scientific << platform.get_velocity(1) <<
    setw(14) << scientific << platform.get_velocity(2) <<
    setw(14) << scientific << platform.get_velocity(3) <<
    setw(14) << scientific << platform.get_velocity(4) <<
    setw(14) << scientific << platform.get_velocity(5) <<
    setw(14) << scientific << platform.get_mooring_load(0) <<
    setw(14) << scientific << platform.get_mooring_load(1) <<
    setw(14) << scientific << platform.get_mooring_load(2) <<
    setw(14) << scientific << platform.get_mooring_load(3) <<
    setw(14) << scientific << platform.get_mooring_load(4) <<
    setw(14) << scientific << platform.get_mooring_load(5) << endl;
    
    target_file.close();
}


////////////////////////////////////////////////////////////////////////////////
// Write cable state into the target file.
////////////////////////////////////////////////////////////////////////////////
void Writer::write(const string target_file_name, const Cable& cable)
{
    // The FILE is used instead of ofstream because using ofstream causes
    // problems in dynamic linking library which we have not fixed.
    FILE *target_file;
    target_file = fopen(target_file_name.c_str(), "w");
    fprintf(target_file,  "        s (m)         x (m)         y (m) "
            "        z (m)         T (N)       S_n (N)       S_b (N) "
            "      u (m/s)       v (m/s)       w (m/s)     phi (rad) "
            "  theta (rad)   kappa_2 (-)   kappa_3 (-) \n");
    
    // Write cable state.
    for (int k=0; k <cable.n_node; k++)
    {
        fprintf(target_file, "% .6E % .6E % .6E % .6E % .6E % .6E % .6E"
                " % .6E % .6E % .6E % .6E % .6E % .6E % .6E\n",
                cable.nodes[k].get_coordinate(),
                cable.nodes[k].get_global_position()(0),
                cable.nodes[k].get_global_position()(1),
                cable.nodes[k].get_global_position()(2),
                cable.nodes[k].get_internal_force()(0),
                cable.nodes[k].get_internal_force()(1),
                cable.nodes[k].get_internal_force()(2),
                cable.nodes[k].get_velocity()(0),
                cable.nodes[k].get_velocity()(1),
                cable.nodes[k].get_velocity()(2),
                cable.nodes[k].get_euler_angle()(0),
                cable.nodes[k].get_euler_angle()(1),
                cable.nodes[k].get_curvature()(1),
                cable.nodes[k].get_curvature()(2));
    }
    fclose(target_file);
}

    
////////////////////////////////////////////////////////////////////////////////
/// Write log file for error check.
////////////////////////////////////////////////////////////////////////////////
int Writer::write_error_message(string file_name,
                                ErrorCode err,
                                ErrorCategory err_cat)
{
    std::time_t result = std::time(nullptr);
    ofstream log_file;
    log_file.open(file_name, std::ofstream::app);
    log_file << std::asctime(std::localtime(&result)) << endl;
    log_file << "  " << err_cat.message(err) << endl << endl;
    log_file.close();
    return 0;
}
} // End of namespace moor.
