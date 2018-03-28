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


#ifndef simulation_h
#define simulation_h

#include "numeric.h"
#include "cable.h"
#include "writer.h"
#include "setting.h"
#include "moorerror.h"
#include "platform.h"
#include <vector>
#include <iostream>
#include <string>

namespace moor {
    
/// \brief Simulation handles varied types of analyses.
///
/// Simulation carries out static analysis or dynamic analysis for a given
/// mooring system with multiple cables. For the static problem, dynamic
/// relaxation and shooting method can be used.
class Simulation
{
    
public:
    
    /// Create a project for analysis.
    Simulation(Platform *p, std::vector<Cable> *c, Setting *s, Writer *w) :
    platform(p), cables(c), setting(s), writer(w) {};
    
    /// Run simulation.
    ErrorCode run(void);
    
    /// Conduct static analysis using Newton-like method.
    ErrorCode shoot(Shooting &shooting);

    /// Solve responses for given excitation time history.
    ErrorCode analyze_forced_motion(ForcedMotion &forced_motion);

    /// Static analysis using dynamic relaxation.
    ErrorCode relax(Relaxation &relaxation);
    
private:
    
    Writer   *writer;
    Setting  *setting;
    Platform *platform;
    
    /// Cable list.
    std::vector< Cable> *cables;
};

} // End of namespace moor.

#endif // simulation_h
