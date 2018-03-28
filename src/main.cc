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


#include "input.h"
#include "moorconfig.h"
#include "reader.h"
#include "cable.h"
#include "solver.h"
#include "setting.h"
#include "moorerror.h"
#include "mooring.h"
#include "platform.h"
#include "simulation.h"
#include "writer.h"
#include <string>
#include <iostream>

using namespace std;
using namespace moor;

////////////////////////////////////////////////////////////////////////////////
/// Print description on the screen.
////////////////////////////////////////////////////////////////////////////////
void print_version(void)
{
    cout << endl << endl;
    cout << "------------------------------------------------------------------"
    "--------------" << endl << endl;
    cout << "      Open Source Simulation Program for MOORing Systems in "
    "Renewable Energy\n";
    fprintf(stdout, " %s Version %d.%d\n", "                          OpenMOOR",
            OpenMOOR_VERSION_MAJOR, OpenMOOR_VERSION_MINOR);
    cout << "                  Copyright (c) 2018 Lin Chen & Biswajit Basu\n\n";
    cout << "------------------------------------------------------------------"
    "--------------" << endl << endl;
}


////////////////////////////////////////////////////////////////////////////////
/// Check the status and quit if an error is found.
////////////////////////////////////////////////////////////////////////////////
void check_status(ErrorCode err_code, ErrorCategory error_domain)
{
    if (err_code != ErrorCode::SUCCESS)
    {
        cout << "  ... Unseccessful: " + error_domain.message(err_code) << endl
        << endl;

		char answer[2];
		printf("\n\n Close the window? (y/n) \n\n");
		scanf("%s", answer);

		if (strcmp(answer, "y") == 0) exit(1);
    }
    else
        cout << "    ... " + error_domain.message(err_code) << endl << endl;
}


////////////////////////////////////////////////////////////////////////////////
/// Driver of the simulation.
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
    ////////////////////////////////////////////////////////////////////////////
    // Preparation.
    ////////////////////////////////////////////////////////////////////////////
    print_version();
    
    cout << "--- PREPARATION" << endl << endl;
    
    ErrorCode err_code;
    ErrorCategory error_domain;
    
    // Get reader ready for reading input files.
    Reader reader;
    // For handling the output.
    Writer writer;
    // Simulation types with corresponding parameters.
    Shooting shooting;
    Relaxation relaxation;
    ForcedMotion forced_motion;

    // Set work folders and input files.
    Setting setting(argv[0],
                    &forced_motion,
                    &shooting,
                    &relaxation);
    
    // Initialization of setting from Setting.xml.
    cout << "  Reading Setting.xml ..." << endl;
    err_code = reader.read_to(setting);
    check_status(err_code, error_domain);
    cout << "  Validating setting data ..." << endl;
    err_code = setting.validate();
    check_status(err_code, error_domain);
    
    // For input data including cable and floater information.
    Constant constant;
    RigidBody reference_point;
    std::vector<Connection> connections;
    std::vector<StructProperty> struct_properties;
    std::vector<HydroProperty> hydro_properties;
    std::vector<SeabedProperty> seabed_properties;
    std::vector<Current> currents;
    std::vector<Solver> solvers;
    std::vector<CableMap> cable_maps;
    
    ////////////////////////////////////////////////////////////////////////////
    // Initialize data from mooring input file and validate.
    ////////////////////////////////////////////////////////////////////////////
    InputData input_data(&constant,
                         &reference_point,
                         &cable_maps,
                         &connections,
                         &struct_properties,
                         &hydro_properties,
                         &seabed_properties,
                         &currents,
                         &solvers);
    
    cout << "  Reading mooring input xml file ..." << endl;
    err_code = reader.read_to(input_data, setting.mooring_file);
    check_status(err_code, error_domain);

    cout << "  Validating mooring input data ..." << endl;
    err_code = input_data.validate();
    check_status(err_code, error_domain);
    
    ////////////////////////////////////////////////////////////////////////////
    // Create main objects: cables, platform, and simulation.
    ////////////////////////////////////////////////////////////////////////////
    vector<Cable> cables;
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
        cables[i].initialize(cable_maps[i],
                             currents[cable_maps[i].i_current],
                             solvers[cable_maps[i].i_solver],
                             struct_properties,
                             hydro_properties,
                             seabed_properties,
                             constant);

    /// Create platform with reference position and fairlead position.
    Platform platform;
    /// Initialize platform position and fairlead number and positions.
    platform.initialize(reference_point.position, n_fairlead,
                        fairlead_position, cables);
    
    // Initialize solvers before assign the solver pointer to cables.
    for (int i=0; i<solvers.size(); i++)
        solvers[i].initialize(10, 5); // For using Euler Angle formulation.

    ////////////////////////////////////////////////////////////////////////////
    // Run simulation.
    ////////////////////////////////////////////////////////////////////////////
    Simulation simulation(&platform, &cables, &setting, &writer);
    err_code = simulation.run();
    check_status(err_code, error_domain);
	char answer[2];
	printf("\n\n Close the window? (y/n) \n\n");
	scanf("%s", answer);

	if (strcmp(answer, "y") == 0) exit(0);
    return 0;
}
