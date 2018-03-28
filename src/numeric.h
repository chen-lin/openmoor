// This file is part of OpenMOOR, an Open-source simulation program for MOORing
// systems in offshore renewable energy applications.
//
// Created by Lin Chen on Aug 29, 2017.
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


#ifndef moor_h
#define moor_h

#include "Eigen/Dense"
#include "Eigen/Geometry"
#include "Eigen/Polynomials"
#include <string>

namespace moor {

/// Types imported and defined from Eigen for matrix manipulation.
using Eigen::Array33d;
using Eigen::VectorXd;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector4d;

using Eigen::VectorXi;
using Eigen::Vector2i;

using Eigen::MatrixXd;
using Eigen::Matrix2d;
using Eigen::Matrix3d;
using Eigen::Matrix3Xd;
using Eigen::Matrix3Xi;
using Eigen::MatrixXi;

using Eigen::AngleAxisd;
    
typedef Eigen::Matrix<double, 6, 1 > Vector6d;
typedef Eigen::Matrix<double,10,10 > Matrix10d;
typedef Eigen::Matrix<double,10, 1 > Vector10d;
typedef Eigen::Matrix<double, 6, 6 > Matrix6d;
typedef Eigen::Matrix<double, 6, Eigen::Dynamic> Matrix6Xd;
typedef Eigen::Matrix<std::string, 1, Eigen::Dynamic> VectorXs;

} // End of namespace moor.

#endif // moor_h
