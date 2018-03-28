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


#ifndef cable_h
#define cable_h

#include "numeric.h"
#include "node.h"
#include "mooring.h"
#include "moorerror.h"
#include "meshoptimizer.h"
#include "solver.h"
#include "catenary.h"
#include "utility.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "stdlib.h"

namespace moor {

/// \brief Cable submerged in fluid.
///
/// Cable is the main class for simulating mooring systems. Currently it is
/// spatially discretized using finite difference method. Refer to
/// \cite tjavaras1996dynamics \cite tjavaras1998mechanics and
/// \cite gobat2000dynamics
/// for the formulation and equations.
class Cable
{

public:
    
    Cable(const double   length_,
          const Vector3d anchor_position_,
          const Vector3d fairlead_position_);
    
    // Initialization.
    void initialize(CableMap &cable_map,
                    Current current,
                    Solver& solver,
                    std::vector< StructProperty >& struct_prop,
                    std::vector< HydroProperty  >& hydro_prop,
                    std::vector< SeabedProperty >& seabed_prop,
                    Constant &constant);
    
    /// Computation of fairlead position error: shooting for static solution.
    Vector3d match_fairlead_position(Vector3d target_fairlead_position);
    
    /// For static problems.
    void form_equations(void);
    
    /// For dynamic problems.
    void form_equations(double time, double dt);
    
    /// Calculate cable nodal position in Cartesian coordinate system from the
    /// first node (seabed anchor: 0,0,0).
    void update_nodal_position(void);
    
    /// Transform the solved fairead force to the floater coordiante system for
    /// assembling the mooring load from cables.
    Vector3d transform_fairlead_force_to_global(void);
    
    /// Transform the fairlead velocity in floater reference frame to the cable
    /// reference frame for solving cable state.
    Vector3d transform_fairlead_velocity(Vector3d fairlead_velocity_global);

    /// Transfrom nodal internal force in Lagrangian reference frame at the
    /// fairlead node to floater referential frame.
    void update_fairlead_force(void);
    
    /// Set current profile function.
    void set_current_profile_function(Current& current);
    
    /// For static problem, solve cable nodal state for given fairlead load.
    ErrorCode update_nodal_state(Vector3d external_load);
    
    /// For dynamic problem when the fairlead velocity is given at a time.
    ErrorCode update_nodal_state(Vector3d forced_velocity,double time,double dt);
    
    /// @name Getters.
    ///
    ///@{
    Vector3d get_fairlead_force(void) const { return fairlead_force; }
    double get_total_fairlead_force(void);
    double get_total_anchor_force(void);
    double get_total_nodal_force(const int i_node);
    double get_averaged_unit_length_weight(void) const
    {
        return averaged_unit_length_weight;
    }
    double get_averaged_axial_stiffness(void) const
    {
        return averaged_axial_stiffness;
    }
    ///@}
    
public:
    
    /// File for initializing cable state. Association matrix will also be read
    /// from this file. When providing initial cable state file, the number of
    /// node in this file should be consistent with n_node. If the initial state
    /// is provided as a space, a catenary solution will be used for
    /// initialization.
    std::string initial_state_file_name;
    
    /// Matrix of [structProp hydroProp seabed] associated with cable nodes.
    MatrixXi association_matrix;
    
    int is_saving;
    
    /// Unstretched length.
    const double length;
    
    /// Number of nodes including the anchor and fairlead.
    int n_node;
    
    /// Cable composed of a number of nodes: at least 3 nodes with two-point
    /// boundary conditions: by default the first and the last nodes are
    /// boundary nodes.
    std::vector< NodeEA > nodes;
    
    /// State variable matrix is assembled from nodes. When using Euler Angle
    /// formulation, the state vector for each node is given as
    /// <br>
    /// \f$[\epsilon~~S_n~~S_b~~u~~v~~w~~\phi~~\theta~~\kappa_2~~\kappa_3]^\top\f$
    /// <br> where <br>
    /// \f$\epsilon\f$ -- strain; <br>
    /// \f$S_n\f$      -- in-plane shear force;<br>
    /// \f$S_b\f$      -- out-of-plane shear force;<br>
    /// \f$u,v,w\f$    -- tangential, normal, and binormal velocities;<br>
    /// \f$\kappa_2\f$  -- in-plane curvature;<br>
    /// \f$\kappa_3\f$  -- out-of-plane curvature.<br>
    MatrixXd state;
    
    VectorXd coordinate;     ///< Vector of all nodal coordinates.
    
    /// Cable stiffness at the fairlead approximated at the initial catenary
    /// configuration.
    Matrix3d linearized_stiffness;
    
    /// Pointing to a solver.
    Solver* solver;
    
    /// Anchor to fairlead in floater coordinate system.
    Vector3d anchor_to_fairlead_global;
    
    /// Chord vector.
    Vector3d anchor_to_fairlead;
    
private:
    
    /// Cable anchor positon.
    const Vector3d anchor_position;
    
    /// Polynomials fitting the transformed current velocity in three direction
    /// as a function of the vertical location.
    MatrixXd current_velocity_polynomials;
    
    /// Define Local coordiante for all cables: the Cartesian coordiante is
    /// located at each anchor and with x-y in the vertical plane determined
    /// by the anchor and the initial fairlead position and x is pointing
    /// upwards. The axis referred to here is defined based on the initial
    /// position of the relevant fairlead.
    Matrix3d axis_rotation_matrix;
    
    /// Averaged structural property along the cable.
    void average_structural_property(CableMap &c_map,
                                     std::vector<StructProperty> &s_props);
    
    double averaged_unit_length_weight;
    double averaged_axial_stiffness;
    double averaged_bending_stiffness;
    
    /// Fairlead force in platform coordinate system.
    Vector3d fairlead_force;
    
    VectorXd segment_length; ///< Distance between nodes.
    
    // The augmented matrix for the coupled neighboring nodes in combination of
    // constraint equations for the boundary condition.
    std::vector< Eigen::MatrixXd > equations;
    
    // For dynamic analysis - to save the states, stiffness and mass matrices,
    // and force vector related to last time step for use at the present step.
    std::vector< Eigen::MatrixXd > last_stiffness;
    
    // Mass matrix of the last time step for all nodes.
    std::vector< Eigen::MatrixXd > last_mass;
    
    // Force matrix of the last time step.
    MatrixXd last_force;
    // State matrix of the last time step.
    MatrixXd last_state;
    // Store the state time derivative for next step.
    MatrixXd state_derivative;
    // Store the state time derivative for next step.
	MatrixXd last_state_derivative;
};

} // End of namespace moor.

#endif // cable_h
