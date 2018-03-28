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


#ifndef moorapi_h
#define moorapi_h

#ifndef __APPLE__
	#define DECLDIR __declspec(dllexport)
#else
	#define DECLDIR
#endif

#ifdef __cplusplus
extern "C"
{
#endif

// Initialization.
int DECLDIR initialize(char input_file_name[]);

// Time stepping.
int DECLDIR update(double* mooring_load, const double displacement[],
                   const double velocity[], const double time,
                   const double dt, const double save_dt);
// Clear all global variables and close the program.
int DECLDIR finish(void);

// Get mooring load component if the index is not in the range, a NAN will be
// returned.
double DECLDIR get_mooring_load(const int i_component);

// Get cable total force at the fairlead if the index is not in the range, a NAN
// will be returned.
double DECLDIR get_cable_fairlead_force(const int i_cable);

// Get cable total force at the fairlead if the index is not in the range, a NAN
// will be returned.
double DECLDIR get_cable_anchor_force(const int i_cable);
    
// Get cable total force at the fairlead if the index is not in the range, a NAN
// will be returned.
double DECLDIR get_cable_nodal_force(const int i_cable, const int i_node);
    
// Get platform reference point displacement component and if the index is not
// in the range, a NAN will be returned.
double DECLDIR get_reference_point_displacement(const int i_component);

// Get flaoter reference point velocity component and if the index is not in
// the range, a NAN will be returned.
double DECLDIR get_reference_point_velocity(const int i_component);

#ifdef __cplusplus
}
#endif

#endif // moorapi_h
