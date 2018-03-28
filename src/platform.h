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


#ifndef platform_h
#define platform_h

#include "numeric.h"
#include "node.h"
#include "cable.h"
#include "moorerror.h"
#include "catenary.h"
#include "solver.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <iostream>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <omp.h>
typedef boost::array<double,12> PlatformState;

namespace moor {
    
/// \brief Platform class.
///
/// A platform of multiple fairleads is considered and each of them is connected
/// to a mooring cable. At present, only one platform can be considered.
class Platform
{
    
public:
    
    Platform(void);
    
    /// Initialize platform position and fairlead number and positions.
    void initialize(const Vector6d&  initial_reference_point_position,
                    const int        n_fairlead,
                    const Matrix3Xd& initial_fairlead_position,
                    std::vector<Cable>& cables);
    
    /// Update fairlead position and velocity for given platform motion.
    void move(const  Vector6d& forced_displacement,
			  const  Vector6d& forced_velocity,
              double time);
    
    // Estimate mooring stiffness based on catenary solution.
    void approximate_mooring_stiffness(std::vector<Cable>& cables);
    
    /// Solve nonlinear mooring load.
    ErrorCode solve_mooring_load(Vector6d& displacement,
                                 Vector6d& velocity,
                                 double    time,
                                 double    dt,
                                 std::vector<Cable>& cables);
    
    // Shooting for solving mooring load for given static displacement.
    ErrorCode solve_mooring_load(Vector6d& displacement,
                                 const double relaxation_factor,
                                 const Vector3d& tolerance,
                                 int shoot_n_limit,
                                 std::vector<Cable>& cables);
    
    /// Time stepping for solving platform state.
    void update_state(const double t, const double dt);
    
    /// Assemble mooring load from cables: nonlinear solution.
    void assemble_mooring_load(std::vector<Cable>& cables);
    
    /// Set platform property for carrying out dynamic relaxation to obtain
    /// the static solution of the coupled system. For dynamic analysis and
    /// using OpenMOOR as a dynamic linking library, this is not required.
    void set_structural_matrix(const Matrix6d& mass,
                               const Matrix6d& damping,
                               const Matrix6d& stiffness,
                               const Vector6d& load);
    
    /// @name Getters.
    ///@{
    Matrix6d get_approximated_mooring_stiffness(void)
    {
        return approximated_mooring_stiffness;
    }
    
    int get_n_fairlead(void) const { return n_fairlead; };
    
    Vector6d get_displacement(void) const { return reference_point_displacement; }
    
    double get_displacement(int i_component) const
    {
        return reference_point_displacement(i_component);
    }
    
    Vector6d get_velocity(void) const { return reference_point_velocity; }
    
    double get_velocity(int i_component) const
    {
        return reference_point_velocity(i_component);
    }
    
    Vector6d get_mooring_load(void) const { return mooring_load; }
    double get_mooring_load(int i_component) const
    {
        return mooring_load(i_component);
    }
    ///@}

private:
    
    /// Update reference point displacement/velocity from the state vector which
    /// is used in time integration.
    void update_from_state(void)
    {
        for (int i=0; i<6; i++)
        {
            reference_point_displacement(i) = state[i];
            reference_point_velocity(i) = state[i+6];
        }
    }
    
    /// Ordinary differential equation of platform motion.
    void diff_fun(const PlatformState& x, PlatformState &dxdt,
                  double t, Vector6d& force);
    
    
private:
    
    /// Displacement of the platform represented by a reference point.
    Vector6d  reference_point_displacement;
    
    /// Velocity of the platform represented by a reference point.
    Vector6d  reference_point_velocity;
    
    /// Platform state observer for time integartion.
    struct Observer
    {
        std::vector<PlatformState>& states;
        Observer (std::vector<PlatformState> &state) : states(state)
        {
        };
        void operator() (const PlatformState &x, double t)
        {
            states.push_back(x);
        }
    };
    
    /// Stiffness due to the buoyancy load is required for static analysis and
    /// dynamic relaxation. Mass and damping matrices are additionally required
    /// for dynamic relaxation.
    Matrix6d mass_matrix;
    
    /// Platform restoring stiffness.
    Matrix6d stiffness_matrix;
    
    /// Platform damping excluding damping effect induced by mooring load.
    Matrix6d damping_matrix;
    
    /// External force considered here is the summation of gravitational and
    /// buoyancy load, for static analysis.
    Vector6d external_load;
    
    /// Mass matrix inversed.
    Matrix6d mass_matrix_inverse;
    
    /// Sub-block of the state matrix in state-spaced model related to platform
    /// mass.
    Matrix6d state_matrix_damping;
    
    /// Sub-block of the state matrix in state-spaced model related to platform
    /// restoring stiffness.
    Matrix6d state_matrix_stiffness;

    /// Number of fairleads.
    /// -- Feb 28, 2017: presently n_cable = n_fairlead = n_anchor.
    int n_fairlead;
    
    /// Platform state for integration.
    PlatformState state;
    
    Vector6d  mooring_load;  ///< Mooring load of six components.
    
    /// Mooring stiffness approximated based on catenary solutions
    Matrix6d approximated_mooring_stiffness;
    
    /// Initialization of the platform position.
    Vector6d initial_reference_point_position;
    
    /// Dependent of initial reference point and platform configuration.
    Matrix3Xd initial_fairlead_position;
    
    /// Derived parameters - in global coordinate.
    Matrix3Xd initial_reference_point_to_fairlead;
    
    /// Present values in global coordinate when in motion.
    Matrix3Xd reference_point_to_fairlead;
    
    /// State of platform in motion.
    Vector6d  reference_point_position;

    /// State of platform in motion.
    Matrix3Xd fairlead_position;
    
    /// State of platform in motion.
    Matrix3Xd fairlead_displacement;
    
    /// State of platform in motion.
    Matrix3Xd fairlead_velocity;
    
    /// Here the fairlead force is collected from the cables and transformed to
    /// the platform reference frame.
    Matrix3Xd fairlead_force;
    
    /// The error between the fairlead position resulting from the platform and
    /// the fairlead position solved from the cables: in cable Cartesian
    /// coordinate system.
    Matrix3Xd fairlead_position_error;
};

} // End of namespace moor.

#endif // platform_h
