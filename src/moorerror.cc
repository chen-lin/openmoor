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


#include "moorerror.h"

namespace moor {

std::string ErrorCategory::message(ErrorCode err)
{
    switch (err)
    {
        case ErrorCode::SUCCESS:
			return "Successful.";
            
        case ErrorCode::MOORING_FILE_NONEXISTENT:
			return "Failed to find main mooring file. It should be a xml file.";

        case ErrorCode::MOORING_FILE_ERROR_PARSE:
            return "Parse error found in main input xml file. Check xml format.";

        case ErrorCode::MOORING_FILE_NO_CASE:
            return "Root node of main input file must be 'case'.";

        case ErrorCode::MOORING_FILE_INCOMPLETE:
            return ("Incomplete components found in mooring input file. "
                    "Required are: 'constants', 'platform', 'connections', "
                    "'cables', 'structuralproperties', 'hydroproperties', "
                    "'seabedproperties', 'currents' and 'solvers'.");

        case ErrorCode::MOORING_FILE_INCOMPLETE_CONSTANTS:
            return ("Incomplete constants definition. Required are "
                    "'gravitationalacceleration', 'waterdensity' and "
                    "'waterdepth'.");

        case ErrorCode::MOORING_FILE_NAN_CONSTANT:
            return "NaN found in constant definitions.";

        case ErrorCode::MOORING_FILE_INCOMPLETE_PLATFORM:
            return ("Failed to find platform position definition or incomplete "
                    "position components provided. Position of the reference "
                    "point is required to have six component: 'x', 'y', 'z' "
                    "'roll', 'pitch' and 'yaw'.");

        case ErrorCode::MOORING_FILE_NAN_PLATFORM_POSITION:
            return "NaN found in platform position definition.";

        case ErrorCode::MOORING_FILE_BAD_CONNECTIONS_NUM:
            return ("Unacceptable 'number' of 'connections': "
					"NaN or nonpositive.");

        case ErrorCode::MOORING_FILE_NO_CONNECTION:
            return "Failed to find at least one connection definition.";

        case ErrorCode::MOORING_FILE_INCOMPLETE_CONNECTION_DEF:
            return ("Incomplete connection definition. Required are 'id', "
                    "'type', 'x', 'y' and 'z'.");

        case ErrorCode::MOORING_FILE_NAN_CONNECTION_DEF:
            return "NaN found in connection 'id' or position components.";

        case ErrorCode::MOORING_FILE_OUTOFRANGE_CONNECTION_ID:
            return ("Found connection 'id' out of range of the 'number' of "
                    "'connections'. ");

        case ErrorCode::MOORING_FILE_UNKNOWN_CONNECTION_TYPE:
            return ("Found unknown connection 'type'. Currently supported "
                    "types are 'anchor' and 'fairlead'.");

        case ErrorCode::MOORING_FILE_MISSED_REPEATED_CONNECTION_ID:
            return ("Found connection missed or repeatedly defined. Check "
                    "connection 'id' and the 'number' of 'connections'.");

        // Cable error.
        case ErrorCode::MOORING_FILE_BAD_CABLES_NUM:
            return "Unacceptable 'number' of 'cables': NaN or nonpositive.";

        case ErrorCode::MOORING_FILE_NO_CABLE_DEF:
            return "Failed to find at least one cable.";

        case ErrorCode::MOORING_FILE_INCOMPLETE_CABLE_DEF:
            return ("Incomplete cable definition. Required are: 'id', "
                    "'initialstatefile (can be empty)', 'icurrent', 'isolver', "
                    "'segmentlength', 'istructproperty', "
                    "'ihydroproperty', 'iseabedproperty', 'iconnection' and "
                    "'saveflag'.");

        case ErrorCode::MOORING_FILE_NAN_CABLE_DEF:
            return ("NaN found in cable 'id', 'icurrent', 'isolver', "
                    "'nodenumber', or 'saveflag'.");

        case ErrorCode::MOORING_FILE_OUTOFRANGE_CABLE_ID:
            return "Found cable 'id' out of range of the 'number' of 'cables'.";

        case ErrorCode::MOORING_FILE_NAN_CABLE_LENGTH:
            return "NaN found in cable 'segmentlength'.";

        case ErrorCode::MOORING_FILE_NAN_CABLE_STRUCTPROP_INDEX:
            return "NaN found in cable 'istructproperty'.";

        case ErrorCode::MOORING_FILE_NAN_CABLE_HYDROPROP_INDEX:
            return "NaN found in cable 'ihydroproperty'.";

        case ErrorCode::MOORING_FILE_NAN_CABLE_SEABEDPROP_INDEX:
            return "NaN found in cable 'iseabedproperty'.";

        case ErrorCode::MOORING_FILE_NAN_CABLE_CONNECTION_INDEX:
            return "NaN found in cable 'iconnection'.";

        case ErrorCode::MOORING_FILE_BAD_CABLE_CONNECTION_INDEX:
            return "Failed to find two connection indexes for a cable.";

        case ErrorCode::MOORING_FILE_MISSED_REPEATED_CABLE_ID:
            return ("Found cable missed or repeatedly defined. Check cable 'id'"
                    " and the 'number' of cables.");

        // Structural property error.
        case ErrorCode::MOORING_FILE_BAD_STRUCTPROPS_NUM:
            return "Unacceptable number of 'structuralproperties'.";

        case ErrorCode::MOORING_FILE_NO_STRUCTPROP:
            return "Failed to find at least one 'structuralproperty'.";

        case ErrorCode::MOORING_FILE_INCOMPLETE_STRUCTPROP_DEF:
            return ("Incomplete 'structuralproperty' data. Required are: 'id', "
                    "'diameter', 'unitlengthmass', 'unitlengthweight', "
                    "'bendingstiffness', 'torsionalstiffness', and "
                    "'dampingcoefficient'");

        case ErrorCode::MOORING_FILE_NAN_STRUCTPROP:
            return "NaN found in 'structuralproperty'.";

        case ErrorCode::MOORING_FILE_OUTOFRANGE_STRUCTPROP_ID:
            return ("Found 'structuralproperty' out of range of the number of "
                    "'structuralproperties'.");

        case ErrorCode::MOORING_FILE_MISSED_REPEATED_STRUCTPROP_ID:
            return ("Found 'structuralproperty' missed or repeatedly defined. "
                    "Check 'structuralproperty' 'id' and the 'number' of "
                    "'structuralproperties'.");

        // Hydro-property error.
        case ErrorCode::MOORING_FILE_BAD_HYDROPROPS_NUM:
            return ("Unacceptable number of 'hydroperoperties': "
					"NaN or nonpositive.");

        case ErrorCode::MOORING_FILE_NO_HYDROPROP:
            return "Unable to find at least one 'hydroproperty'.";

        case ErrorCode::MOORING_FILE_INCOMPLETE_HYDROPROP_DEF:
            return ("Incomplete hydroproperty data. Required are 'id', "
                    "'addedmasscoefficient' and 'dragcoefficient'");

        case ErrorCode::MOORING_FILE_NAN_HYDROPROP:
            return "NaN found in 'hydroproperty' 'id'.";

        case ErrorCode::MOORING_FILE_OUTOFRANGE_HYDROPROP_ID:
            return "Found 'hydroproperty' 'id' out of range.";

        case ErrorCode::MOORING_FIEL_INCOMPLETE_HYDRO_COEFFICIENTS:
            return ("Incomplete addedmasscoefficient or dragcoefficient data."
					"Required are 'tangential', 'normal' and 'binormal' "
					"components");
            
        case ErrorCode::MOORING_FIEL_NAN_HYDRO_COEFFICIENTS:
            return ("NaN found in 'addedmasscoefficient' or 'dragcoefficient' "
					"data");
            
        case ErrorCode::MOORING_FILE_MISSED_REPEATED_HYDROPROP_ID:
            return ("Found 'hydroproperty' missed or repeatedly defined. "
                    "Check 'hydroproperty' 'id' and the 'number' of "
                    "'hydroproperties'.");

        // Seabed property error.
        case ErrorCode::MOORING_FILE_BAD_SEABEDPROPS_NUM:
            return "Unacceptable number of 'seabedproperties'.";

        case ErrorCode::MOORING_FILE_NO_SEABEDPROP:
            return "Unable to find at least one 'seabedproperty' definition.";

        case ErrorCode::MOORING_FILE_INCOMPLETE_SEABEDPROP_DEF:
            return ("Incomplete seabedproperty data. Required are 'id', "
                    "'dampingcoefficient' and 'stiffnesscoefficient'.");

        case ErrorCode::MOORING_FILE_NAN_SEABEDPROP:
            return ("NaN found in 'seabedproperty' data.");

        case ErrorCode::MOORING_FILE_OUTOFRANGE_SEABEDPROP_ID:
            return ("Seabedproperty 'id' is out of range!\n");

        case ErrorCode::MOORING_FILE_MISSED_REPEATED_SEABEDPROP_ID:
            return ("Found 'seabedproperty' missed or repeatedly defined. "
                    "Check 'seabedproperty' 'id' and the number of "
                    "'seabedproperties'.");

        // Current error.
        case ErrorCode::MOORING_FILE_BAD_CURRENTS_NUM:
            return "Unacceptable of number of 'currents'.";

        case ErrorCode::MOORING_FILE_NO_CURRENT:
            return "Unable find at least one current definition";

        case ErrorCode::MOORING_FILE_INCOMPLETE_CURRENT_DEF:
            return ("Incomplete current data. Required are 'id', "
                    "'polyorder' and 'profilefile'.");

        case ErrorCode::MOORING_FILE_NAN_CURRENT:
            return "NaN found in current 'id' or 'polyorder'.";

        case ErrorCode::MOORING_FILE_OUTOFRANGE_CURRENT_ID:
            return "Found out of range current 'id'.";

        case ErrorCode::MOORING_FILE_MISSED_REPEATED_CURRENT_ID:
            return ("Found 'current' missed or repeatedly defined. "
                    "Check 'current' 'id' and the number of 'currents'.");

        case ErrorCode::MOORING_FILE_BAD_CURRENT_DATA:
            return ("Unacceptable current profile data or NaN found. Check "
					"current data file. Current data file should have one line "
					"of header and data matrix with 6 columns containing "
					"three coordinates in global reference system and three "
					"current velocities.");
            
        // Solver error.
        case ErrorCode::MOORING_FILE_BAD_SOLVERS_NUM:
            return "Unacceptable number of 'solvers'.";

        case ErrorCode::MOORING_FILE_NO_SOLVER:
            return "Unable to find at least one 'solver' definition.";

        case ErrorCode::MOORING_FILE_INCOMPLETE_SOLVER_DEF:
            return ("Incomplete solver definition. Required are 'id', "
                    "'iterationnumberlimit', 'convergencetolerance', "
                    "'initialrelaxationfactor', 'relaxationincreasefactor', "
                    "'relaxationdecreasefactor' and 'lambdainfinity'");

        case ErrorCode::MOORING_FILE_NAN_SOLVER:
            return "NaN found in 'solver' definition.";

        case ErrorCode::MOORING_FILE_OUTOFRANGE_SOLVER_ID:
            return ("Found solver 'id' out of range of the 'number' of "
					"'solvers'.");

        case ErrorCode::MOORING_FILE_MISSED_REPEATED_SOLVER_ID:
            return ("Found 'solver' missed or repeatedly defined. "
                    "Check 'solver' 'id' and the number of 'solvers'.");

        // Setting error.
        case ErrorCode::SETTING_FILE_NONEXISTENT:
            return "Failed to find Setting.xml file.";

        case ErrorCode::SETTING_FILE_ERROR_PARSE:
            return "Found xml parse error in Setting.xml.";

        case ErrorCode::SETTING_FILE_NO_SETTING_NODE:
            return "Root node of Setting.xml must be setting.";

        case ErrorCode::SETTING_FILE_NO_SIMULATION_TYPE:
            return "Failed to find 'simulationtype' definition in Setting.xml.";
            
        case ErrorCode::SETTING_FILE_UNKNOWN_SIMULATION_TYPE:
            return ("Found unknown 'simulationtype'. Currently supported are "
                    "'shooting', 'relaxation', and 'forcedmotion'.");

        case ErrorCode::SETTING_FILE_NO_MOORING_FILE_DEF:
            return "Failed to find 'mooringinputfile' definition in Setting.xml.";

        case ErrorCode::SETTING_FILE_NO_SHOOTING_PARA:
            return "Unable to find 'shooting' parameter definition.";

        case ErrorCode::SETTING_FILE_INCOMPLETE_SHOOTING_PARA:
            return ("Incomplete parameters for shooting. Required are "
                    "'fairleadpositiontolerance', "
                    "'fairleadforcerelaxationfactor', "
                    "'fairleadpositioniterationlimit', "
                    "'platformpositioniterationlimit', "
                    "'platformdisplacementtolerance', "
                    "'platformdisplacementrelaxationfactor', "
                    "'cableoutofplanestiffness', "
                    "'platformhydrostaticstiffness' and "
                    "'platformotherload'.");

        case ErrorCode::SETTING_FILE_NAN_SHOOTING_PARA:
            return ("NaN found in shooting parameters: "
                    "'fairleadpositiontolerance', "
                    "'fairleadforcerelaxationfactor', "
                    "'fairleadpositioniterationlimit', "
                    "'platformpositioniterationlimit', "
                    "'platformdisplacementtolerance' or "
                    "'platformdisplacementrelaxationfactor'.");

        case ErrorCode::SETTING_FILE_INCOMPLETE_SHOOTING_FAIRLEAD_TOLERANCE:
            return ("Incomplete definition of 'shooting' parameter "
                    "fairleadpositiontolerance.");

        case ErrorCode::SETTING_FILE_NAN_SHOOTING_FAIRLEAD_TOLERANCE:
            return ("NaN found in 'shooting' parameter "
					"'fairleadpositiontolerance'.");

        case ErrorCode::SETTING_FILE_BAD_SHOOTING_STIFFNESS:
            return ("NaN or not 36 elements found in "
                    "'platformhydrostaticstiffness' for "
                    "'shooting'.");

        case ErrorCode::SETTING_FILE_BAD_SHOOTING_LOAD:
            return ("NaN or not 6 elements found in 'platformotherload' for "
                    "'shooting'.");

        case ErrorCode::SETTING_FILE_NO_RELAXATION_PARA:
            return "Unable to find 'relaxation' parameter definition.";

        case ErrorCode::SETTING_FILE_INCOMPLETE_RELAXATION:
            return ("Incomplete parameters for 'relaxation'. Required are "
                    "'platformvelocitytolerance', 'cablevelocitytolerance', "
                    "'stoptime', 'timestep', 'platformmass', "
                    "'platformdamping', 'platformhydrostaticstiffness' and "
                    "'platformotherload'.");

        case ErrorCode::SETTING_FILE_NAN_RELAXATION_PARA:
            return ("NaN found in 'relaxation' parameters: "
                    "'platformvelocitytolerance', 'cablevelocitytolerance'"
                    "'stoptime', or 'timestep'.");

        case ErrorCode::SETTING_FILE_BAD_RELAXATION_MASS:
            return ("NaN or not 36 elements found in "
                    "'platformmass' for 'relaxtion'.");

        case ErrorCode::SETTING_FILE_BAD_RELAXATION_DAMPING:
            return ("NaN or not 36 elements found in "
                    "'platformdamping' for 'relaxtion'.");

        case ErrorCode::SETTING_FILE_BAD_RELAXATION_STIFFNESS:
            return ("NaN or not 36 elements found in "
                    "'platformhydrostaticstiffness' for 'relaxtion'.");

        case ErrorCode::SETTING_FILE_BAD_RELAXATION_LOAD:
            return ("NaN or not 6 elements found in 'platformotherload' for "
                   "'relaxtion'.");

        case ErrorCode::SETTING_FILE_NO_MOTION:
            return "Unable to find 'forcedmotion' node in Setting.xml.";

        case ErrorCode::SETTING_FILE_NO_TIME_HISTORY:
            return "Unable to find 'timehistory' file definition.";

        // Validation of mooring input data.
        case ErrorCode::MOORING_INVALID_CONSTANT:
            return "Found nonpositive constants.";
        
        case ErrorCode::MOORING_INVALID_CABLE_CURRENT_INDEX:
            return "Cable 'icurrent' is out of range.";
        
        case ErrorCode::MOORING_INVALID_CABLE_SOLVER_INDEX:
            return "Cable 'isolver' is out of range.";
            
        case ErrorCode::MOORING_INVALID_CABLE_NODE_NUM:
            return "Cable 'nodenumber' is negative.";
            
        case ErrorCode::MOORING_INVALID_CABLE_PROP_ASSOCIATION:
            return ("Failed to find at least one group of cable length and "
                    "associated property or found inconsistent cable segments "
                    "and properties definitions.");
            
        case ErrorCode::MOORING_INVALID_CABLE_CONNECTION_INDEX:
            return ("Two same ends found for a cable or cable connection index "
                    "out of range. In addtion, currently, the first connection "
                    "must be an 'anchor' and the second must be a 'fairlead'.");
            
        case ErrorCode::MOORING_INVALID_CABLE_SEGMENT_LENGTH:
            return "Cable segmentlength should be positive.";
            
        case ErrorCode::MOORING_INVALID_STRUCT_PROP:
            return ("Positive values required for diameter, unitlengthmass, "
                    " unitlengthweight, and axialstiffness.");
            
        case ErrorCode::MOORING_INVALID_HYDRO_COEFFICIENT:
            return ("Nonnegative values required for bending stiffness, "
                    "torsional stiffness and damping coefficient.");
            
        case ErrorCode::MOORING_INVALID_SEABED:
            return ("Nonnegative values required hydrodynamic coefficients.");
            
        case ErrorCode::MOORING_INVALID_CURRENT:
            return ("Invalid current data. Required are that polyorder is "
                    "nonnegative and less then the number of profile data "
                    "points and at least one data point provided. In addition, "
                    "the Z coordinate should be negative and decreasing "
                    "monotonically");
            
        case ErrorCode::MOORING_INVALID_SOLVER:
            return ("Invalid solver found: iterationnumberlimit, "
                    "convergencetolerance, "
                    "initialrelaxationfactor, relaxationincreasefactor, "
                    "and relaxationdecreasefactor should be positive numbers; "
                    "relaxationincreasefactor and relaxationincreasefactor "
                    "should be no less than 1; lambdainfinity should not be equal "
                    "to 1 and often between -1 and 0.");
            
        case ErrorCode::MOORING_INVALID_CABLE_STATE:
            return ("Initial cable state matrix read successfully, however, "
                    "the cable coordinate is not valid. The cable coordinate "
                    "should be positive, begining with zero, increasing "
                    "monotonically and terminating with full cable length "
                    "consistent with the sum of segmentlength and at "
                    "least two nodes are required. Check the first column of "
                    "initial state file and the 'segmentlength'.");
            
        case ErrorCode::SHOOTING_BAD_PARA:
            return ("Found nonpositive tolerance or relaxation factor or "
                    "negative iteration limit in shooting control.");

        case ErrorCode::RELAXATION_BAD_PARA:
            return ("Found nonpositive tolerance, steptime, or timestep.");
            
        case ErrorCode::CURRENT_FILE_NONEXISTENT:
            return "Unable to find current profile data file.";

        case ErrorCode::CURRENT_FILE_BAD_DATA:
            return ("Found NaN or wrong column number (expect 6) of current "
                    "profile data.");
            
        case ErrorCode::MOTION_FILE_NONEXISTENT:
            return "Unable to find forcedmotion time history data.";
            
        case ErrorCode::MOTION_FILE_BAD_DATA:
            return ("Found NaN or wrong column number (expect 7) of forced "
                    "motion time history.");
            
        case ErrorCode::FORCED_MOTION_BAD_TIME:
            return "Found negative or decreasing time in the forced velocity.";
            
        case ErrorCode::TIME_STEP_TOO_SMALL:
            return "Time step too small, no need to update.";
            
        case ErrorCode::SINGULAR_MATRIX_SOLVER:
            return ("Singularity found in solving the cable equation. "
                    "Check input.");
            
        case ErrorCode::NAN_CABLE_SOLUTION:
            return ("NaN found in cable state solution.");
    }
}
    
} // End of namespace moor.
