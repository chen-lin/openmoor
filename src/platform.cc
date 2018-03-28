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


#include "platform.h"
using namespace std;
using namespace boost::numeric::odeint;

namespace moor {
    
Platform::Platform(void)
{
    initial_reference_point_position.Zero();
    reference_point_position.Zero();
    reference_point_displacement.Zero();
    reference_point_velocity.Zero();
    mooring_load.Zero();
    state = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	
	mass_matrix.Zero();
	stiffness_matrix.Zero();
	damping_matrix.Zero();
	external_load.Zero();
	mass_matrix_inverse.Zero();
}


////////////////////////////////////////////////////////////////////////////////
// Initialization.
////////////////////////////////////////////////////////////////////////////////
void Platform::initialize(const Vector6d&  initial_reference_point_position_,
                          const int        n_fairlead_,
                          const Matrix3Xd& initial_fairlead_position_,
                          std::vector<Cable>& cables)
{
    // Assignment.
    initial_reference_point_position = initial_reference_point_position_;
    n_fairlead = n_fairlead_;
    
    // Resize the fairlead position matrix according to the fairlead number.
    initial_fairlead_position.resize(3,n_fairlead);
    initial_fairlead_position = initial_fairlead_position_;

    // Get derived parameters.
    initial_reference_point_to_fairlead = (initial_fairlead_position
                                           - initial_reference_point_position
                                           .block(0,0,3,1).rowwise()
                                           .replicate(n_fairlead));
    
    reference_point_to_fairlead = initial_reference_point_to_fairlead;
    reference_point_position = initial_reference_point_position;
    
    // Initialization of matrices.
    fairlead_position.setZero(3,n_fairlead);
    fairlead_displacement.setZero(3,n_fairlead);
    fairlead_velocity.setZero(3,n_fairlead);

    fairlead_force.setZero(3,n_fairlead);
    fairlead_position_error.setZero(3,n_fairlead);
    
    // Initialization of the mooring load;
    reference_point_displacement.setZero();
    reference_point_velocity.setZero();
    assemble_mooring_load(cables);
}


////////////////////////////////////////////////////////////////////////////////
// Rigid body motion.
////////////////////////////////////////////////////////////////////////////////
void Platform::move(const  Vector6d& forced_displacement,
                    const  Vector6d& forced_velocity,
                    double time)
{
    // Update platform velocity and displacement by assignment.
    reference_point_displacement = forced_displacement;
    reference_point_velocity = forced_velocity;
    reference_point_position = (initial_reference_point_position
                                + reference_point_displacement);
    
    // Update relative fairlead position using updated rotation matrix.
    Matrix3d rotate_matrix;
    rotate_matrix
    = (AngleAxisd(reference_point_displacement(5), Vector3d::UnitZ())
       * AngleAxisd(reference_point_displacement(4), Vector3d::UnitY())
       * AngleAxisd(reference_point_displacement(3), Vector3d::UnitX()));
    
    reference_point_to_fairlead = (rotate_matrix
                                   * initial_reference_point_to_fairlead);
    
    Vector3d rotate_velocity = reference_point_velocity.block(3,0,3,1);
    
    // Get the fairlead velocity in the platform (global) coordinate.
    for (int i_fairlead = 0; i_fairlead<n_fairlead; i_fairlead++)
        fairlead_velocity.col(i_fairlead)
        = (rotate_velocity.transpose()
           .cross(reference_point_to_fairlead.col(i_fairlead)).transpose()
           + reference_point_velocity.block(0,0,3,1));
    
    // Get the fairlead displacement in the floater coordinate.
    fairlead_position = (reference_point_to_fairlead
                         + reference_point_position.block(0,0,3,1)
                         .rowwise().replicate(n_fairlead));
    
    fairlead_displacement = fairlead_position - initial_fairlead_position;
}

    
////////////////////////////////////////////////////////////////////////////////
// Set the structural matrices for the platform for dynamic relaxation
// or coupled analysis.
////////////////////////////////////////////////////////////////////////////////
void Platform::set_structural_matrix(const Matrix6d& mass_matrix_,
                                     const Matrix6d& damping_matrix_,
                                     const Matrix6d& stiffness_matrix_,
                                     const Vector6d& external_load_)
{
    mass_matrix = mass_matrix_;
    damping_matrix = damping_matrix_;
    stiffness_matrix = stiffness_matrix_;
    external_load = external_load_;
    mass_matrix_inverse = mass_matrix.inverse();
    state_matrix_damping = - mass_matrix_inverse * damping_matrix;
    state_matrix_stiffness = - mass_matrix_inverse * stiffness_matrix;
}
 

////////////////////////////////////////////////////////////////////////////////
/// First-order differential equation describing the floater motion.
////////////////////////////////////////////////////////////////////////////////
void Platform::diff_fun(const PlatformState& x,
                        PlatformState& dx,
                        double t,
                        Vector6d& load)
{
    for (int k=0; k<6; k++)
        dx[k] = x[k+6];
    
    for (int k=6; k<12; k++)
    {
        double temp = 0.0;
        for (int j=0; j<6; j++)
        {
            temp = (temp + state_matrix_damping(k-6,j) * x[k]
                    + state_matrix_stiffness(k-6,j) * x[k-6]
                    + mass_matrix_inverse(k-6,j) * load(k-6));
        }
        dx[k] = temp;
    }
}


////////////////////////////////////////////////////////////////////////////////
// Assembling mooring load from the fairlead load in cable reference frame.
////////////////////////////////////////////////////////////////////////////////
void Platform::assemble_mooring_load(std::vector<Cable>& cables)
{
    // Cable fairlead force transformed to floater reference frame.
    mooring_load.setZero();
    
    for (int i=0; i<n_fairlead; i++)
    {
        // Force transformation and reaction force.
        fairlead_force.col(i) = -cables[i].transform_fairlead_force_to_global();
        mooring_load.block(0,0,3,1) += fairlead_force.col(i);
        mooring_load.block(3,0,3,1) += (reference_point_to_fairlead.col(i)
                                        .cross(fairlead_force.col(i)));
    }
}


////////////////////////////////////////////////////////////////////////////////
// Solve mooring load for dynamic problem.
////////////////////////////////////////////////////////////////////////////////
ErrorCode Platform::solve_mooring_load(Vector6d& displacement,
                                       Vector6d& velocity,
                                       double time,
                                       double dt,
                                       std::vector<Cable>& cables)
{
    move(displacement,velocity,time);
    vector<ErrorCode> err_code;
    err_code.resize(n_fairlead);
    
    //=========================== Start parallel computing. ====================
    int i, thread_id;
#pragma omp parallel private(thread_id, i)
    {
        i = 0; // cable index
#pragma omp for
        
        for (i=0; i<n_fairlead; i++)
        {
            Vector3d fairlead_velocity_transformed
            = cables[i].transform_fairlead_velocity(fairlead_velocity.col(i));
            
            err_code[i] = cables[i].update_nodal_state(fairlead_velocity_transformed,
                                                    time, dt);
            cables[i].update_fairlead_force();
        }
    }
    //======================== End of parallel computing. ======================
    
    for (i=0; i<n_fairlead; i++)
        if (err_code[i] != ErrorCode::SUCCESS) return err_code[i];
    
    assemble_mooring_load(cables);
    return ErrorCode::SUCCESS;
}
    

////////////////////////////////////////////////////////////////////////////////
/// Use shooting method for solving mooring load for given static displacement
/// following steps as below:
////////////////////////////////////////////////////////////////////////////////
ErrorCode Platform::solve_mooring_load(Vector6d& displacement,
                                       const double relaxation_factor,
                                       const Vector3d& tolerance,
                                       int shoot_n_limit,
                                       std::vector<Cable>& cables)
{
    Vector6d velocity;
    velocity.setZero();
    
    /// - Force platform to move.
    move(displacement, velocity, 0.0);
    vector<ErrorCode> err_code; err_code.resize(n_fairlead);
    
    /// - Initialization of trial fairlead forces.
    for (int i=0; i < n_fairlead; i++)
        fairlead_force.col(i) = cables[i].get_fairlead_force();
    
    /// - Update cable state using parallel computation.
    // ======================= Start of parallel computing.================== //
    int i, thread_id;
    
#pragma omp parallel private(thread_id, i)
    {
        i = 0; // Cable index.
        
#pragma omp for
        
        for (i=0; i<n_fairlead; i++)
        {
            int criterion = 1, i_shoot_fairlead = 1;
            
            while (criterion)
            {
                err_code[i] = cables[i].update_nodal_state(fairlead_force.col(i));
                
                fairlead_position_error.col(i)
                = cables[i].match_fairlead_position(fairlead_position.col(i));
                
                fairlead_force.col(i) -= (relaxation_factor
                                          * cables[i].linearized_stiffness
                                          * fairlead_position_error.col(i));
                
                criterion
                = ((fabs(fairlead_position_error(0,i)) > tolerance(0) ||
                    fabs(fairlead_position_error(1,i)) > tolerance(1) ||
                    fabs(fairlead_position_error(2,i)) > tolerance(2) ) &&
                   (i_shoot_fairlead < shoot_n_limit) );
                
                std::cout << "     | cable " << i  << "  fairlead shooting: "
                << i_shoot_fairlead << " solved with " << setw(3)
                << cables[i].solver->get_n_iteration() << std::endl <<
                "         position err. (xc,yc,zc):" << scientific <<
                fairlead_position_error.col(i).transpose() << std::endl;
                
                i_shoot_fairlead++;
            }
            cables[i].update_fairlead_force();
        }
    }
    
    // ======================= End of parallel computing. =================== //
    for (i=0; i<n_fairlead; i++)
        if (err_code[i] != ErrorCode::SUCCESS) return err_code[i];
    
    /// - Assemble mooring load from cables.
    assemble_mooring_load(cables);
    return ErrorCode::SUCCESS;
}
    

////////////////////////////////////////////////////////////////////////////////
// Time integration to update platform state.
////////////////////////////////////////////////////////////////////////////////
void Platform::update_state(const double t, const double dt)
{
    Vector6d total_load = external_load + mooring_load;
    std::vector<PlatformState> states;
    
    auto fun = [&] (PlatformState &x, PlatformState &dx, double t) ->
    void { diff_fun(x,dx,t,total_load); };
    
    states.clear();
    integrate_const(runge_kutta4<PlatformState>(), fun, state, t, t+dt, dt,
                    Observer(states));
    
    state = states.back();
    update_from_state();
}


////////////////////////////////////////////////////////////////////////////////
/// Estimate mooring stiffness from catenary solutions using the formulation
/// described in \cite al2016stiffness. Note that two-dimensional catenary
/// solutions are used as the basis and hence current and inertia effects are
/// not included.
////////////////////////////////////////////////////////////////////////////////
void Platform::approximate_mooring_stiffness(std::vector<Cable>& cables)
{
    double phi, theta, psi, s1, c1, s2, c2, s3, c3;
    
    phi   = reference_point_displacement(3);
    theta = reference_point_displacement(4);
    psi   = reference_point_displacement(5);
    
    s1 = sin(phi);   c1 = cos(phi);
    s2 = sin(theta); c2 = cos(theta);
    s3 = sin(psi);   c3 = cos(psi);
    
    Vector3d dXp_phi, dXp_theta, dXp_psi;
    
    Matrix3d mat1, mat2, mat3;
    
    mat1 <<  0.0  ,  c3*s2*c1+s3*s1, -c3*s2*s1+s3*c1,
             0.0  ,  s3*s2*c1-c3*s1, -s3*s2*s1-c3*c1,
             0.0  ,  c2*c1         , -c2*s1         ;
    
    mat2 << -c3*s2,  c3*c2*s1      ,  c3*c2*c1      ,
            -s3*s2,  s3*c2*s1      ,  s3*c2*c1      ,
            -c2   , -s2*s1         , -s2*c1         ;
    
    mat3 << -s3*c2, -s3*s2*s1-c3*c1, -s3*s2*c1+c3*s1,
             c3*c2,  c3*s2*s1-s3*c1,  c3*s2*c1+s3*s1,
             0.0  ,  0.0           ,  0.0;

    std::vector<Matrix6d> K(n_fairlead);
    
    for (int i = 0; i < n_fairlead; i++)
    {
        dXp_phi   = mat1 * initial_reference_point_to_fairlead.col(i);
        dXp_theta = mat2 * initial_reference_point_to_fairlead.col(i);
        dXp_psi   = mat3 * initial_reference_point_to_fairlead.col(i);
        
        double l = sqrt(pow(cables[i].anchor_to_fairlead(1),2.0) +
                        pow(cables[i].anchor_to_fairlead(2),2.0) );
        
        double h = cables[i].anchor_to_fairlead(0);
        
        Vector2d chord;
        chord << l,h;
        
        Catenary cat(chord, cables[i].length,
                     cables[i].get_averaged_unit_length_weight(),
                     cables[i].get_averaged_axial_stiffness());
        
//        cat.estimate_fairlead_stiffness();
        // Also update cable linearized stiffness;
        cables[i].linearized_stiffness(0,0)=cat.tangent_fairlead_stiffness(1,1);
        cables[i].linearized_stiffness(1,1)=cat.tangent_fairlead_stiffness(0,0);
        cables[i].linearized_stiffness(0,1)=cat.tangent_fairlead_stiffness(1,0);
        cables[i].linearized_stiffness(1,0)=cat.tangent_fairlead_stiffness(0,1);
        
        double H = cat.bound_force(0,1), V = cat.bound_force(1,1);
        
        Matrix2d Kc;
        Kc = cat.tangent_fairlead_stiffness.block(0,0,2,2);
        
        double beta;
        beta = atan(cables[i].anchor_to_fairlead_global(1)
                    / cables[i].anchor_to_fairlead_global(0));
        
        if (cables[i].anchor_to_fairlead_global(0) < 0.0)
            beta += M_PI;
        
        double sin_beta = sin(beta), cos_beta = cos(beta);
        double dl_phi, dl_theta, dl_psi, dh_phi, dh_theta, dh_psi, dbeta_phi,
        dbeta_theta, dbeta_psi;
        
        dl_phi   = cos_beta * dXp_phi(0) + sin_beta * dXp_phi(1);
        dl_theta = cos_beta * dXp_theta(0) + sin_beta * dXp_theta(1);
        dl_psi   = cos_beta * dXp_psi(0) + sin_beta * dXp_psi(1);
        
        dh_phi   = dXp_phi(2);
        dh_theta = dXp_theta(2);
        dh_psi   = dXp_psi(2);
        
        dbeta_phi = 1/l*(cos_beta * dXp_phi(1) - sin_beta * dXp_phi(0));
        dbeta_theta = 1/l*(cos_beta * dXp_theta(1) - sin_beta * dXp_theta(0));
        dbeta_psi = 1/l*(cos_beta * dXp_psi(1) - sin_beta * dXp_psi(0));
        
        K[i](0,0) = Kc(0,0) * pow(cos_beta,2.0) + H/l * pow(sin_beta,2.0);
        K[i](0,1) = sin_beta * cos_beta * (Kc(0,0) - H/l);
        K[i](0,2) = cos_beta * Kc(0,1);
        K[i](0,3) = (cos_beta * (Kc(0,0) * dl_phi + Kc(0,1) * dh_phi)
                     - H * sin_beta  * dbeta_phi);
        K[i](0,4) = (cos_beta * (Kc(0,0) * dl_theta + Kc(0,1) * dh_theta)
                     - H * sin_beta * dbeta_theta);
        K[i](0,5) = (cos_beta * (Kc(0,0) * dl_psi + Kc(0,1) * dh_psi)
                     - H * sin_beta * dbeta_psi);
        
        K[i](1,0) = K[i](0,1);
        K[i](1,1) = pow(sin_beta,2.0) * Kc(0,0) + pow(cos_beta,2.0) * H/l;
        K[i](1,2) = sin_beta * Kc(0,1);
        K[i](1,3) = (sin_beta * (Kc(0,0) * dl_phi + Kc(0,1) * dh_phi)
                     + H * cos_beta  * dbeta_phi);
        K[i](1,4) = (sin_beta * (Kc(0,0) * dl_theta + Kc(0,1) * dh_theta)
                     + H * cos_beta  * dbeta_theta);
        K[i](1,5) = (sin_beta * (Kc(0,0) * dl_psi + Kc(0,1) * dh_psi)
                     + H * cos_beta  * dbeta_psi);
        
        K[i](2,0) = K[i](0,2);
        K[i](2,1) = K[i](1,2);
        K[i](2,2) = Kc(1,1);
        K[i](2,3) = Kc(1,0) * dl_phi + Kc(1,1) * dh_phi;
        K[i](2,4) = Kc(1,0) * dl_theta + Kc(1,1) * dh_theta;
        K[i](2,5) = Kc(1,0) * dl_psi + Kc(1,1) * dh_psi;
        
        double dX, dY, dZ;
        dX = fairlead_position(0,i) - reference_point_position(0);
        dY = fairlead_position(1,i) - reference_point_position(1);
        dZ = fairlead_position(2,i) - reference_point_position(2);
        
        K[i](3,0) = dY * K[i](2,0) - dZ * K[i](1,0);
        K[i](3,1) = dY * K[i](2,1) - dZ * K[i](1,1);
        K[i](3,2) = dY * K[i](2,2) - dZ * K[i](1,2);
        K[i](3,3) = (dY * K[i](2,3) - dZ * K[i](1,3) + V * dXp_phi(1)
                     - H * sin_beta * dXp_phi(2));
        K[i](3,4) = (dY * K[i](2,4) - dZ * K[i](1,4) + V * dXp_theta(1)
                     - H * sin_beta * dXp_theta(2));
        K[i](3,5) = (dY * K[i](2,5) - dZ * K[i](1,5) + V * dXp_psi(1)
                     - H * sin_beta * dXp_psi(2));
        
        K[i](4,0) = dZ * K[i](0,0) - dX * K[i](2,0);
        K[i](4,1) = dZ * K[i](0,1) - dX * K[i](2,1);
        K[i](4,2) = dZ * K[i](0,2) - dX * K[i](2,2);
        K[i](4,3) = (dZ * K[i](0,3) - dX * K[i](2,3) + H *cos_beta * dXp_phi(2)
                     - V * dXp_phi(0));
        K[i](4,4) = (dZ * K[i](0,4) - dX * K[i](2,4) + H *cos_beta * dXp_theta(2)
                     - V * dXp_theta(0));
        K[i](4,5) = (dZ * K[i](0,5) - dX * K[i](2,5) + H *cos_beta * dXp_psi(2)
                     - V * dXp_psi(0));
        
        K[i](5,0) = K[i](0,5);
        K[i](5,1) = K[i](1,5);
        K[i](5,2) = K[i](2,5);
        K[i](5,3) = (dX * K[i](1,3) - dY * K[i](0,3) + H *sin_beta *dXp_phi(0)
                     - H * cos_beta * dXp_phi(1));
        K[i](5,4) = (dX * K[i](1,3) - dY * K[i](0,3) + H *sin_beta *dXp_theta(0)
                     - H * cos_beta * dXp_theta(1));
        K[i](5,5) = (dX * K[i](1,3) - dY * K[i](0,3) + H *sin_beta *dXp_psi(0)
                     - H * cos_beta * dXp_psi(1));
    }
    
    approximated_mooring_stiffness.setZero();
    for (int i=0; i<n_fairlead; i ++)
        approximated_mooring_stiffness += K[i];
}

} // End of namespace moor.
