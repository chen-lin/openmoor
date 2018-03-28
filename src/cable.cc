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


#include "cable.h"
using namespace std;
namespace moor {

////////////////////////////////////////////////////////////////////////////////
/// Define a cable by its initial end positions and undeformed length. Here both
/// anchor and fairlead positions are given in global coordiante system. The
/// input variables are:
/// - \a length_ is the undeformed cable length hanging between to points; and
/// - Coordinates of the two global position of the two end points.
////////////////////////////////////////////////////////////////////////////////
Cable::Cable(const double length_,
             const Vector3d anchor_position_,
             const Vector3d fairlead_position):
    length(length_),
    anchor_position(anchor_position_)
{
    // Initial anchor to fairlead vector.
    anchor_to_fairlead_global = fairlead_position - anchor_position;
    
    // The Cartesian coordinate system is defined by the anchor position and
    // the initial fairlead position.
    axis_rotation_matrix.setZero();
    Vector3d axis_rotation_angle;
    
    // Initializing the axis rotation angle. Assuming that the X-axis is the
    // in vertical direction.
    axis_rotation_angle.fill(-M_PI/2);
    axis_rotation_angle(2) = atan(anchor_to_fairlead_global(1)
                                  / anchor_to_fairlead_global(0));
    
    if (anchor_to_fairlead_global(0) < 0.0) axis_rotation_angle(2) += M_PI;

    axis_rotation_matrix
    = (AngleAxisd(-axis_rotation_angle(0), Vector3d::UnitZ())
       * AngleAxisd(-axis_rotation_angle(1), Vector3d::UnitX())
       * AngleAxisd(-axis_rotation_angle(2), Vector3d::UnitZ()));
    
    // Get the chord vector for catenary initialization.
    anchor_to_fairlead = axis_rotation_matrix * anchor_to_fairlead_global;
    
    // For saving;
    is_saving = 0;
};


////////////////////////////////////////////////////////////////////////////////
/// Initialize cable according to cable map and associated properties.
////////////////////////////////////////////////////////////////////////////////
void Cable::initialize(CableMap &cable_map,
                       Current current,
                       Solver& solver_,
                       std::vector< StructProperty >& struct_props,
                       std::vector< HydroProperty  >& hydro_props,
                       std::vector< SeabedProperty >& seabed_props,
                       Constant &constant)
{
    // Assignment.
    is_saving = cable_map.is_saving;
    solver = &solver_;
    
    // Use read input state or catenary solution.
    n_node = (cable_map.initial_state.rows() > 1 ?
              cable_map.initial_state.rows() : cable_map.n_node);
    
    // Resizing.
    coordinate.resize(n_node);
    segment_length.resize(n_node-1);
    state.resize(10,n_node);
    association_matrix.resize(3,n_node);
    equations.resize(n_node+1);
    
    // Initialization.
    coordinate.setZero();
    fairlead_force.setZero();
    state.setZero();
    association_matrix.setZero();
    
    // For dynamic analysis.
    state_derivative.resize(10,n_node);
    state_derivative.setZero();
    last_state.resize(10,n_node);
    last_state.setZero();
    last_force.resize(10,n_node);
    last_force.setZero();
    last_mass.resize(n_node+1);
    last_stiffness.resize(n_node+1);
    last_state_derivative.resize(10,n_node);
    last_state_derivative.setZero();
    
    for (int i=0; i<n_node+1; i++)
    {
        equations[i].resize(10,21);
        equations[i].setZero();
        last_mass[i].resize(10,10);
        last_mass[i].setZero();
        last_stiffness[i].resize(10,10);
        last_stiffness[i].setZero();
    }
    
    average_structural_property(cable_map, struct_props);
    
    /// Copy initial state data
    if (cable_map.initial_state.rows() > 1)
    {
        coordinate = cable_map.initial_state.col(0);
        state = cable_map.initial_state.block(0,1,n_node,10).transpose();
    }
    else
    {
        /// - Initialization using catenary approximation. Note that when using
        /// the catenary solution the association matrix is assumed all zeros
        /// and hence uniform structural property along the cable is considered.
        Vector2d chord;
        chord(0) = pow((pow(anchor_to_fairlead(1),2.0)
                        + pow(anchor_to_fairlead(2),2.0)), 0.5);
        chord(1) = anchor_to_fairlead(0);
        
        // Use the averaged structural property as the property along the cable
        // for calculating the catenary solution.
        Catenary cat(chord, length, averaged_unit_length_weight,
                     averaged_axial_stiffness);
        
//        cat.estimate_fairlead_stiffness();
        
        linearized_stiffness.setZero();
        
        // Noting the index of the catenary fairlead stiffness if different
        // from the cable linearized stiffness matrix.
        linearized_stiffness(0,0) = cat.tangent_fairlead_stiffness(1,1);
        linearized_stiffness(1,1) = cat.tangent_fairlead_stiffness(0,0);
        linearized_stiffness(0,1) = cat.tangent_fairlead_stiffness(1,0);
        linearized_stiffness(1,0) = cat.tangent_fairlead_stiffness(0,1);
//        linearized_stiffness.block(0,0,2,2) = cat.tangent_fairlead_stiffness;

        // SET A PROPER VALUE AS THE OUT-OF-PLANE STIFFNESS.
        linearized_stiffness(2,2) = 1000;

        // Initialize and approximate fairlead load by catenary solution.
        fairlead_force(0) = cat.bound_force(1, 1);
        fairlead_force(1) = cat.bound_force(0, 1);

        cat.discretize(n_node, averaged_bending_stiffness);

        vector<double> present_mesh, reference_variable;

        for (int i=0; i<n_node; i++)
        {
            present_mesh.push_back(cat.discrete_state(i,0));
            reference_variable.push_back(cat.discrete_state(i,4));
        }

        MeshOptimizer mesh_optimizer(present_mesh, reference_variable);
        mesh_optimizer.optimize(n_node, 0.1);

        VectorXd new_coordinate;
        new_coordinate.resize(mesh_optimizer.optimized_mesh.size());

        for (int i=0; i<mesh_optimizer.optimized_mesh.size(); i++)
            new_coordinate(i) = mesh_optimizer.optimized_mesh[i];

        // Get state estimate according to the new mesh.
        cat.discretize(n_node, new_coordinate, averaged_bending_stiffness);

        coordinate   = cat.discrete_state.col(0);
        state.row(0) = cat.discrete_state.col(1).transpose();
        state.row(1) = cat.discrete_state.col(2).transpose();
        state.row(6) = cat.discrete_state.col(3).transpose();
        state.row(9) = cat.discrete_state.col(4).transpose();
    }

    // Get segment length based on coordinate.
    segment_length = (coordinate.block(1,0,n_node-1,1)
                      - coordinate.block(0,0,n_node-1,1));
    
    // Get current velocity polynomials for current velocity interpolation.
    set_current_profile_function(current);

    // Calculate association matrix based on the segment length;
    VectorXi ind;
    ind.resize(coordinate.rows());
    ind.setZero();

    VectorXd segment_end_length;
    segment_end_length.resize(cable_map.segment_length.rows()+1);
    segment_end_length.setZero();
    
    for (int i=0; i<cable_map.segment_length.rows();i++)
        segment_end_length(i+1) = cable_map.segment_length.block(0,0,i+1,1).sum();

    for (int i=0; i<n_node; i++)
    {
        for (int j=0; j<cable_map.segment_length.rows(); j++)
            if (coordinate(i) <= segment_end_length(j+1)
                && coordinate(i) >= segment_end_length(j))
                ind(i) = j;
    }
    
    // - Create nodes: requiring structural, hydrodynamic, seabed properties and
    // association matrix.
    for(int i=0; i<n_node; i++)
    {
        NodeEA node(coordinate(i),
                    &struct_props[cable_map.i_struct_property(ind(i))],
                    &hydro_props[cable_map.i_hydro_property(ind(i))],
                    &seabed_props[cable_map.i_seabed_property(ind(i))],
                    &constant);

        // Reorder the variables, very important.
        for (int j=0; j<10; j++)
            node.state(node.index(j)) = state(j,i);

        // Update the cable state as well, important!
        state.col(i) = node.state;
        nodes.push_back(node);
    }
    
    /// - Add fixed boundary conditions: the seabed anchor as the first node.
    nodes[0].add_constraint();

    /// - Set state-independent for each node: mass and stiffness matrices and
    /// force vector. For dynamic problem it is also important to initialize
    /// the cable state, nodal mass and stiffness matrice and force vector at
    /// the last step.
    for (int k=0; k<n_node; k++)
    {
        nodes[k].set_state_independent();
        nodes[k].update_from_state();
    }
    update_nodal_position();
    update_fairlead_force();
}
    

// Static problem: form augmented matrix for solving cable states.
void Cable::form_equations(void)
{
    equations[0].block(0,10,10,11) = nodes[0].equation;
    equations[n_node].block(0,0,10,10) = nodes[n_node-1].equation.block(0,0,10,10);
    equations[n_node].col(20) = nodes[n_node-1].equation.col(10);
    
    Matrix10d jk0, jk1;
    // For static problem the central difference is used.
    for (int i=0; i<n_node-1; i++) 
    {
        jk0.setZero();
        jk1.setZero();
        for (int j=0; j<10; j++) 
        {
            jk0.row(j)
            = ((nodes[i+1].state - nodes[i].state).transpose()
               * nodes[i].stiffness_jacobian[j] / segment_length(i));
            
            jk1.row(j)
            = ((nodes[i+1].state - nodes[i].state).transpose()
               * nodes[i+1].stiffness_jacobian[j] / segment_length(i));
        }
        
        equations[i+1].block(0,0,10,10)
        = ((-nodes[i+1].stiffness - nodes[i].stiffness)/segment_length(i)
           + nodes[i].force_jacobian + jk0);
        
        equations[i+1].block(0,10,10,10)
        = (( nodes[i+1].stiffness + nodes[i].stiffness)/segment_length(i)
           + nodes[i+1].force_jacobian + jk1);
        
        equations[i+1].col(20)
        = (( nodes[i+1].stiffness + nodes[i].stiffness)
           * ( nodes[i+1].state - nodes[i].state) / segment_length(i)
           + nodes[i+1].force + nodes[i].force);
    }
}


////////////////////////////////////////////////////////////////////////////////
/// Dynamic problem: form augmented matrix for solving cable states.
////////////////////////////////////////////////////////////////////////////////
void Cable::form_equations(double time, double dt)
{
    // Set matrices for boundary nodes.
    equations[0].block(0,10,10,11) = nodes[0].equation;
    equations[n_node].block(0,0,10,10) = nodes[n_node-1].equation.block(0,0,10,10);
    equations[n_node].col(20) = nodes[n_node-1].equation.col(10);
    
    Matrix10d k0, k1, jm0, jm1, jk0, jk1;
    Vector10d y0ds, y1ds;
    k0.setZero(); k1.setZero(); jm0.setZero(); jm1.setZero(); jk0.setZero(); 
	jk1.setZero(); y0ds.setZero(); y1ds.setZero(); 
    
    for (int i=0; i<n_node-1; i++)
    {
        y0ds = (last_state.col(i+1) - last_state.col(i)) / segment_length(i);
        y1ds = (nodes[i+1].state - nodes[i].state) / segment_length(i);
        
        k1 = nodes[i].stiffness + nodes[i+1].stiffness;
        k0 = last_stiffness[i] + last_stiffness[i+1];
        
        Vector10d ys, yt0, yt1;
        yt0.setZero(); ys.setZero(); yt1.setZero();
        ys  = solver->alpha_k1_square * y1ds + solver->alpha_k_cross * y0ds ;
        yt0 = (solver->alpha_m1_square * state_derivative.col(i)
               + solver->alpha_m_cross * last_state_derivative.col(i));
        yt1 = (solver->alpha_m1_square * state_derivative.col(i+1)
               + solver->alpha_m_cross * last_state_derivative.col(i+1));
        
        // Jacobian matrix.
        for (int j=0; j<10; j++)
        {
            jk0.row(j) = ys.transpose() * nodes[i].stiffness_jacobian[j];
            jk1.row(j) = ys.transpose() * nodes[i+1].stiffness_jacobian[j];
            jm0.row(j) = yt0.transpose() * nodes[i].mass_jacobian[j];
            jm1.row(j) = yt1.transpose() * nodes[i+1].mass_jacobian[j];
        }
        
        equations[i+1].block(0,0,10,10)
        = (jk0 + jm0
           - solver->alpha_k1_square / segment_length(i) * k1
           - solver->alpha_k_cross / segment_length(i) * k0
		   + (solver->alpha_m1_square * nodes[i].mass 
			   + solver->alpha_m_cross * last_mass[i]) / dt / solver->gamma
           + solver->alpha_k1 * nodes[i].force_jacobian);

        equations[i+1].block(0,10,10,10)
        = (jk1 + jm1
           + solver->alpha_k1_square / segment_length(i) * k1
           + solver->alpha_k_cross / segment_length(i) * k0
		   + (solver->alpha_m1_square * nodes[i+1].mass
			   + solver->alpha_m_cross * last_mass[i+1]) / dt / solver->gamma
           + solver->alpha_k1 * nodes[i+1].force_jacobian );
        
		equations[i+1].col(20)
			= ((solver->alpha_m1_square * nodes[i].mass
				+ solver->alpha_m_cross*last_mass[i]) * state_derivative.col(i)
				+ (solver->alpha_m1_square * nodes[i + 1].mass
					+ solver->alpha_m_cross * last_mass[i + 1]) * state_derivative.col(i + 1)
				+ (solver->alpha_m_cross * nodes[i].mass
					+ solver->alpha_m_square * last_mass[i]) * last_state_derivative.col(i)
				+ (solver->alpha_m_cross * nodes[i + 1].mass
					+ solver->alpha_m_square * last_mass[i + 1]) * last_state_derivative.col(i + 1)
				+ (solver->alpha_k1_square * k1 + solver->alpha_k_cross * k0) * y1ds
				+ (solver->alpha_k_cross * k1 + solver->alpha_k_square * k0) * y0ds
				+ solver->alpha_k1 * (nodes[i].force + nodes[i + 1].force)
				+ solver->alpha_k * (last_force.col(i) + last_force.col(i + 1)));
    }
}


////////////////////////////////////////////////////////////////////////////////
/// Update nodal position according to the nodal state.
////////////////////////////////////////////////////////////////////////////////
void Cable::update_nodal_position(void)
{
    nodes[0].position.setZero();
    for (int k=1; k<n_node; k++)
    {
        nodes[k].position(0) = (nodes[k-1].position(0)
                                + ((1+nodes[k].strain)
                                   * nodes[k].cos_angle(0)
                                   * nodes[k].cos_angle(1)
                                   + (1+nodes[k-1].strain)
                                   * nodes[k-1].cos_angle(0)
                                   * nodes[k-1].cos_angle(1) )
                                * 0.5 * segment_length(k-1));
        
        nodes[k].position(1) = (nodes[k-1].position(1)
                                + ((1+nodes[k].strain)
                                   * nodes[k].sin_angle(0)
                                   * nodes[k].cos_angle(1)
                                   + (1+nodes[k-1].strain)
                                   * nodes[k-1].sin_angle(0)
                                   * nodes[k-1].cos_angle(1) )
                                * 0.5 *segment_length(k-1));
        
        nodes[k].position(2) = (nodes[k-1].position(2)
                                + (-(1+nodes[k].strain)
                                   * nodes[k].sin_angle(1)
                                   - (1+nodes[k-1].strain)
                                   * nodes[k-1].sin_angle(1))
                                * 0.5 *segment_length(k-1));
    }
    
    // Also update the global position of nodes.
    for (int k=0; k<n_node; k++)
    {
        nodes[k].global_position = (axis_rotation_matrix.transpose()
                                    * nodes[k].position + anchor_position);
    }
}

    
////////////////////////////////////////////////////////////////////////////////
/// Solve and update the state of each node.
////////////////////////////////////////////////////////////////////////////////
ErrorCode Cable::update_nodal_state(Vector3d fairlead_force)
{
    for (int k=0; k<n_node; k++)
        nodes[k].update_from_state();
    
    // Update nodal position.
    update_nodal_position();
    
    for (int k=0; k<n_node; k++)
    {
        nodes[k].update_effective_weight(1);
        
        // Update current velocity according to nodal vertical position.
        nodes[k].update_current_velocity(current_velocity_polynomials);
        
        nodes[k].update_state_dependent();
    }
    
    // Add fairlead constraint and the corresponding matrices.
    nodes[n_node-1].add_constraint(fairlead_force);
    
    // Compute nodal solution.
    double previous_error = 0.0, added_error, state_error;
    solver->reset_relaxation_factor();
    
    int fail;
    for (int it=1; it<solver->get_iteration_limit(); it++)
    {
        form_equations();
        
        fail = solver->solve(equations);
        if (fail) return ErrorCode::SINGULAR_MATRIX_SOLVER;
            
        added_error = 0.0;
        for (int i=0; i<state.rows(); i++) 
        {
            state_error = 0.0;
            
            for (int k=0; k<n_node; k++)
                state_error += fabs(equations[k](i,20));
            
            added_error += state_error;
        }
        added_error = added_error / (state.rows() * n_node);
        
        // Adjust relaxation factor for fast convergence.
        solver->adjust_relaxation(added_error, previous_error);
        previous_error = added_error;
        
        // Get nodal solution from the manipulated augmented matrix.
        for (int k=0; k<n_node; k++) 
        {
            state.col(k) -= (solver->get_relaxation_factor()
                             * equations[k].col(20));
            
            // Update each node state.
            nodes[k].state = state.col(k);
            nodes[k].update_from_state();
        }
        if (state.hasNaN()) return ErrorCode::NAN_CABLE_SOLUTION;
        
        update_nodal_position();
        for (int k=0; k < n_node; k++) 
        {
            // Update hydrodynamic parameters.
            nodes[k].update_hydrodynamics();
            nodes[k].update_effective_weight(1);
            // Update current velocity according to nodal vertical position.
            nodes[k].update_current_velocity(current_velocity_polynomials);
            nodes[k].update_state_dependent();
        }

        // The augmented matrix of the last node is altered with state.
        nodes[n_node-1].add_constraint(fairlead_force);
        
        if (solver->check_convergence(added_error))
        {
			solver->update_n_iteration(it);
            break;
        }
    }
    
    return ErrorCode::SUCCESS;
}


////////////////////////////////////////////////////////////////////////////////
/// Solve and update the state of each node.
////////////////////////////////////////////////////////////////////////////////
ErrorCode Cable::update_nodal_state(Vector3d forced_velocity,
                                    double time,
                                    double dt)
{
    for (int k=0; k<n_node; k++)
        nodes[k].update_from_state();
    
    // Update nodal position.
    update_nodal_position();
    
    for (int k=0; k<n_node; k++)
    {
        nodes[k].update_effective_weight(0);
        // Update current velocity according to nodal vertical position.
        nodes[k].update_current_velocity(current_velocity_polynomials);
        nodes[k].update_state_dependent();
    }
    
    for (int k=0; k<n_node; k++)
    {
        last_state.col(k) = nodes[k].state;
        last_force.col(k) = nodes[k].force;
        last_stiffness[k] = nodes[k].stiffness;
        last_mass[k] = nodes[k].mass;
    }
    
	last_state_derivative = state_derivative;
	state_derivative = -solver->gamma1/solver->gamma * last_state_derivative;

    // Add fairlead constraint and the corresponding matrices.
    nodes[n_node-1].add_constraint(forced_velocity, time);

    // Get nodal solution.
    double previous_error = 0.0, added_error, state_error;

    solver->reset_relaxation_factor();
    int fail;
    for (int it=1; it < solver->get_iteration_limit(); it++)
    {
        form_equations(time, dt);
        
        fail = solver->solve(equations);
        if (fail) return ErrorCode::SINGULAR_MATRIX_SOLVER;
            
        added_error = 0.0;
        
        for (int i=0; i<state.rows(); i++)
        {
            state_error = 0.0;
            for (int k=0; k<n_node; k++)
                state_error += fabs(equations[k](i,20));
            
            added_error += state_error;
        }
        added_error = added_error / (state.rows() * n_node);

        // Adjust relaxation for convergence.
        solver->adjust_relaxation(added_error, previous_error);
        previous_error = added_error;

        // Get nodal solution from the manipulated augmented matrix.
        for (int k=0; k<n_node; k++) 
        {
            state.col(k) -= (solver->get_relaxation_factor()
                             * equations[k].col(20));
            
            // Update each node state.
            nodes[k].state = state.col(k);
            nodes[k].update_from_state();
        }
        
        if (state.hasNaN()) return ErrorCode::NAN_CABLE_SOLUTION;

		// Update state time derivative.
		state_derivative = (1/dt/solver->gamma * (state - last_state)
                            - (1-solver->gamma)/solver->gamma
                            * last_state_derivative);
		
        update_nodal_position();
        
        for (int k=0; k < n_node; k++)
        {
            // Update hydrodynamic parameters.
            nodes[k].update_hydrodynamics();
            nodes[k].update_effective_weight(0);
            // Update current velocity according to nodal vertical position.
            nodes[k].update_current_velocity(current_velocity_polynomials);
            nodes[k].update_state_dependent();
        }
        
        // The augmented matrix of the last node is altered along with state.
        nodes[n_node-1].add_constraint(forced_velocity, time);
        
        if (solver->check_convergence(added_error))
        {
            solver->update_n_iteration(it);
            break;
        }
    }
    return ErrorCode::SUCCESS;
}


////////////////////////////////////////////////////////////////////////////////
/// Obtain the fairlead velocity in cable Cartesian coordinate system.
////////////////////////////////////////////////////////////////////////////////
Vector3d Cable::transform_fairlead_velocity(Vector3d fairlead_velocity_global)
{
    Vector3d fairlead_velocity = (axis_rotation_matrix
                                  * fairlead_velocity_global);
    return fairlead_velocity;
}


////////////////////////////////////////////////////////////////////////////////
/// Update the fairlead force of the cable in the Cartesian coordinate.
////////////////////////////////////////////////////////////////////////////////
void Cable::update_fairlead_force()
{
    fairlead_force = nodes[n_node-1].transform_internal_force();
}


////////////////////////////////////////////////////////////////////////////////
/// Tranformatin of fairlead force to global coordinate system.
////////////////////////////////////////////////////////////////////////////////
Vector3d Cable::transform_fairlead_force_to_global(void)
{
    Vector3d fairlead_force_global = (axis_rotation_matrix.transpose()
                                      * fairlead_force);
    return fairlead_force_global;
}


////////////////////////////////////////////////////////////////////////////////
/// For static analysis.
////////////////////////////////////////////////////////////////////////////////
Vector3d Cable::match_fairlead_position(Vector3d target_fairlead_position)
{
    // Fairlead position in cable Cartesian reference frame.
    anchor_to_fairlead_global = target_fairlead_position - anchor_position;
    anchor_to_fairlead = axis_rotation_matrix * anchor_to_fairlead_global;
    
    Vector3d fairlead_position_error = (nodes.back().get_position()
                                        - anchor_to_fairlead);
    return fairlead_position_error;
}


////////////////////////////////////////////////////////////////////////////////
/// Average the structural property over the cable for estimating the mooring
/// stiffness and also to use catenary solution for initialization.
////////////////////////////////////////////////////////////////////////////////
void Cable::average_structural_property(CableMap &c_map,
                                        std::vector<StructProperty> &s_props)
{
    double w = 0, EA = 0, EI = 0;
    for (int i=0; i<c_map.segment_length.rows(); i++)
    {
        w += (c_map.segment_length(i)
              * s_props[c_map.i_struct_property(i)].unit_length_weight);
        EA += (c_map.segment_length(i)
               * s_props[c_map.i_struct_property(i)].axial_stiffness);
        EI += (c_map.segment_length(i)
               * s_props[c_map.i_struct_property(i)].bending_stiffness);
    }
    // Consider that the full length of the cable is submerged.
    averaged_unit_length_weight = w/length;
    averaged_axial_stiffness = EA/length;
    averaged_bending_stiffness = EI/length;
}

    
////////////////////////////////////////////////////////////////////////////////
/// Calculate total fairlead force.
////////////////////////////////////////////////////////////////////////////////
double Cable::get_total_fairlead_force(void)
{
    return nodes[n_node-1].get_total_force();
};
    
    
////////////////////////////////////////////////////////////////////////////////
/// Calculate total anchor force.
////////////////////////////////////////////////////////////////////////////////
double Cable::get_total_anchor_force(void)
{
    return nodes[0].get_total_force();
};


////////////////////////////////////////////////////////////////////////////////
/// Calculate total force at a specific node.
////////////////////////////////////////////////////////////////////////////////
double Cable::get_total_nodal_force(const int i_node)
{
    if (i_node >=0 && i_node < n_node)
        return nodes[i_node].get_total_force();
    else
        return NAN;
};
    
////////////////////////////////////////////////////////////////////////////////
// Setup the current profile function during initialization.
////////////////////////////////////////////////////////////////////////////////
void Cable::set_current_profile_function(Current& current)
{
    // Get current velocity polynomials for current velocity interpolation.
    MatrixXd coordinate_transformed(3, current.profile_data.rows());
    MatrixXd velocity_transformed(3, current.profile_data.rows());
    
    coordinate_transformed
    = (axis_rotation_matrix
       * current.profile_data.block(0,0,current.profile_data.rows(),3).transpose());
    
    velocity_transformed
    = (axis_rotation_matrix
       * current.profile_data.block(0,3,current.profile_data.rows(),3).transpose());
    
    current_velocity_polynomials.resize(current.polyfit_order+1,3);
    
    for (int i=0; i<3; i++)
    {
        VectorXd x = coordinate_transformed.transpose().col(0);
        VectorXd y = velocity_transformed.transpose().col(i);
        VectorXd z;
        polyfit(x,y,z,current.polyfit_order);
        current_velocity_polynomials.col(i) = z;
    }
}

} // End of namespace moor.
