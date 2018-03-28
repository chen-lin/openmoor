// This file is part of OpenMOOR, an Open-source simulation program for MOORing
// systems in offshore renewable energy applications.
//
// Created by Lin Chen on Mar 6, 2017.
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


#include "reader.h"

using namespace std;
using namespace rapidxml;

namespace moor {

////////////////////////////////////////////////////////////////////////////////
/// Get data from the input file.
////////////////////////////////////////////////////////////////////////////////
ErrorCode Reader::read_to(InputData& data_in, const std::string input_file_name)
{
    int flag = check_file_existence(input_file_name);
    if (flag)
        return ErrorCode::MOORING_FILE_NONEXISTENT;
    
    // Get the work folder.
    size_t work_folder_index = input_file_name.find_last_of("/\\");
    std::string work_folder = input_file_name.substr(0, work_folder_index + 1);
    
    // Load main input file providing cable information.
    ifstream file(input_file_name);
    
    // Check status of the input file.
    stringstream buffer;
    buffer << file.rdbuf();
    file.close();
    
    // Decode and extract the data in buffer.
    try
    {
        xml_document<> doc;
        std::string content(buffer.str());
        doc.parse<0> (&content[0]);
    }
    catch (const rapidxml::parse_error& e)
    {
        return ErrorCode::MOORING_FILE_ERROR_PARSE;
    }
    
    xml_document<> doc;
    std::string content(buffer.str());
    doc.parse<0> (&content[0]);
    
    xml_node<> *root_node;
    root_node = doc.first_node("case");
    if (root_node == 0)
        return ErrorCode::MOORING_FILE_NO_CASE;
    
    vector<string> component_name;
    component_name.push_back("constants");
    component_name.push_back("platform");
    component_name.push_back("connections");
    component_name.push_back("cables");
    component_name.push_back("structuralproperties");
    component_name.push_back("hydroproperties");
    component_name.push_back("seabedproperties");
    component_name.push_back("currents");
    component_name.push_back("solvers");
    
    if (!check_availability(root_node, component_name))
        return ErrorCode::MOORING_FILE_INCOMPLETE;
    
    int index = 0;
    xml_node<>* subnode;
    xml_node<>* sub2node;
    VectorXi index_count;
    
    for (xml_node<> *child_node = root_node->first_node();
         child_node; child_node = child_node->next_sibling())
    {
        ////////////////////////////////////////////////////////////////////////
        // Constants.
        ////////////////////////////////////////////////////////////////////////
        if (strcmp(child_node->name(), "constants") == 0)
        {
            vector<string> constant_name;
            constant_name.push_back("gravitationalacceleration");
            constant_name.push_back("waterdensity");
            constant_name.push_back("waterdepth");
            
            if (!check_availability(child_node, constant_name))
                return ErrorCode::MOORING_FILE_INCOMPLETE_CONSTANTS;

            if (!check_is_number(child_node, constant_name))
                return ErrorCode::MOORING_FILE_NAN_CONSTANT;

            // Define constant of gravity.
            subnode = child_node->first_node("gravitationalacceleration");
            data_in.constant->GRAV_ACC = stod(subnode->value());

            // Define water density.
            subnode = child_node->first_node("waterdensity");
            data_in.constant->WATER_DENSITY = stod(subnode->value());

            // Define water depth.
            subnode = child_node->first_node("waterdepth");
            data_in.constant->WATER_DEPTH = stod(subnode->value());
        }
        ////////////////////////////////////////////////////////////////////////
        // Platform input data.
        ////////////////////////////////////////////////////////////////////////
        else if (strcmp(child_node->name(), "platform") == 0)
        {
            subnode = child_node->first_node("position");
            if (subnode==0)
                return ErrorCode::MOORING_FILE_INCOMPLETE_PLATFORM;
            else
            {
                vector<string> names;
                names.push_back("x");
                names.push_back("y");
                names.push_back("z");
                names.push_back("roll");
                names.push_back("pitch");
                names.push_back("yaw");
                
                if (!check_availability(subnode, names))
                    return ErrorCode::MOORING_FILE_INCOMPLETE_PLATFORM;

                if (!check_is_number(subnode, names))
                    return ErrorCode::MOORING_FILE_NAN_PLATFORM_POSITION;

                data_in.reference_point->position
                <<  stod(subnode->first_attribute("x")->value()),
                    stod(subnode->first_attribute("y")->value()),
                    stod(subnode->first_attribute("z")->value()),
                    stod(subnode->first_attribute("roll")->value()),
                    stod(subnode->first_attribute("pitch")->value()),
                    stod(subnode->first_attribute("yaw")->value());
            }
        }
        ////////////////////////////////////////////////////////////////////////
        // Connection input data.
        ////////////////////////////////////////////////////////////////////////
        else if (strcmp(child_node->name(), "connections") == 0)
        {
            data_in.n_connection = read_number(child_node);
            if (data_in.n_connection == -1)
                return ErrorCode::MOORING_FILE_BAD_CONNECTIONS_NUM;

            // Resize the connection container.
            data_in.connections->resize(data_in.n_connection);
            if (child_node->first_node("connection")==0)
                return ErrorCode::MOORING_FILE_NO_CONNECTION;

            index_count.resize(data_in.n_connection);
            index_count.setConstant(0);

            for (subnode = child_node->first_node("connection"); subnode;
                 subnode = subnode->next_sibling("connection"))
            {
                vector<string> names;
                names.push_back("id");
                names.push_back("x");
                names.push_back("y");
                names.push_back("z");
                names.push_back("type");

                if (!check_availability(subnode, names))
                    return ErrorCode::MOORING_FILE_INCOMPLETE_CONNECTION_DEF;

                names.erase(names.end()-1);
                if (!check_is_number(subnode, names))
                    return ErrorCode::MOORING_FILE_NAN_CONNECTION_DEF;

                index = stoi(subnode->first_attribute("id")->value());

                if (index < 0 || index >= data_in.n_connection)
                    return ErrorCode::MOORING_FILE_OUTOFRANGE_CONNECTION_ID;

                index_count(index) = index_count(index) + 1;

                if (strcmp(subnode->first_attribute("type")->value(), "anchor") == 0)
                    data_in.connections->at(index).type = ConnectionType::ANCHOR;
                else if (strcmp(subnode->first_attribute("type")->value(), "fairlead") == 0)
                    data_in.connections->at(index).type = ConnectionType::FAIRLEAD;
                else
                    return ErrorCode::MOORING_FILE_UNKNOWN_CONNECTION_TYPE;

                data_in.connections->at(index).position
                <<  stod(subnode->first_attribute("x")->value()),
                    stod(subnode->first_attribute("y")->value()),
                    stod(subnode->first_attribute("z")->value());
            }

            int crit = 1;
            for (int i=0; i<index_count.rows(); i++)
                crit = (crit && (index_count(i) == 1));
            if (!crit)
                return ErrorCode::MOORING_FILE_MISSED_REPEATED_CONNECTION_ID;
        }
        ////////////////////////////////////////////////////////////////////////
        // Cable input data.
        ////////////////////////////////////////////////////////////////////////
        else if (strcmp(child_node->name(), "cables") == 0)
        {
            data_in.n_cable = read_number(child_node);
            
            if (data_in.n_cable == -1)
                return ErrorCode::MOORING_FILE_BAD_CABLES_NUM;
            
            data_in.cable_maps->resize(data_in.n_cable);

            if (child_node->first_node("cable")==0)
                return ErrorCode::MOORING_FILE_NO_CABLE_DEF;

            index_count.resize(data_in.n_cable);
            index_count.setConstant(0);

            for (subnode = child_node->first_node("cable"); subnode;
                 subnode = subnode->next_sibling("cable"))
            {
                vector<string> names;
                names.push_back("id");
                names.push_back("initialstatefile");
                names.push_back("icurrent");
                names.push_back("isolver");
                names.push_back("nodenumber");
                names.push_back("segmentlength");
                names.push_back("istructproperty");
                names.push_back("ihydroproperty");
                names.push_back("iseabedproperty");
                names.push_back("iconnection");
                names.push_back("saveflag");
                if (!check_availability(subnode, names))
                    return ErrorCode::MOORING_FILE_INCOMPLETE_CABLE_DEF;
                
                names.clear();
                names.push_back("id");
                names.push_back("icurrent");
                names.push_back("isolver");
                names.push_back("nodenumber");
                names.push_back("saveflag");

                if (!check_is_number(subnode, names))
                    return ErrorCode::MOORING_FILE_NAN_CABLE_DEF;

                index = stoi(subnode->first_attribute("id")->value());

                if (index < 0 || index >= data_in.n_cable)
                    return ErrorCode::MOORING_FILE_OUTOFRANGE_CABLE_ID;
                
                index_count(index) = index_count(index) + 1;

                data_in.cable_maps->at(index).n_node
                = stoi(subnode->first_node("nodenumber")->value());

                data_in.cable_maps->at(index).i_current
                = stoi(subnode->first_node("icurrent")->value());

                data_in.cable_maps->at(index).i_solver
                = stoi(subnode->first_node("isolver")->value());

                data_in.cable_maps->at(index).is_saving
                = stoi(subnode->first_node("saveflag")->value());

                data_in.cable_maps->at(index).initial_state_file
                = work_folder + subnode->first_node("initialstatefile")->value();

                int flag
                = check_file_existence(data_in.cable_maps->at(index).initial_state_file);
                if (!flag)
                    int failure = read_to(data_in.cable_maps->at(index).initial_state,
                                          data_in.cable_maps->at(index).initial_state_file,
                                          11, 1);
                
                vector<string> number_string;

                // Cable segment indexes.
                flag = extract_multi_number(subnode->first_node("segmentlength")->value(), number_string);

                if (!flag)
                    return ErrorCode::MOORING_FILE_NAN_CABLE_LENGTH;

                data_in.cable_maps->at(index).segment_length.resize(number_string.size());
                for (int i=0; i<number_string.size(); i++)
                    data_in.cable_maps->at(index).segment_length(i)
                    = stod(number_string[i]);

                // Structural property indexes.
                number_string.clear();
                flag = extract_multi_number(subnode->first_node("istructproperty")->value(), number_string);

                if (!flag)
                    return ErrorCode::MOORING_FILE_NAN_CABLE_STRUCTPROP_INDEX;

                data_in.cable_maps->at(index).i_struct_property.resize(number_string.size());
                for (int i=0; i<number_string.size(); i++)
                    data_in.cable_maps->at(index).i_struct_property(i)
                    = stoi(number_string[i]);

                // Hydro-property indexes.
                number_string.clear();
                flag = extract_multi_number(subnode->first_node("ihydroproperty")->value(), number_string);

                if (!flag)
                    return ErrorCode::MOORING_FILE_NAN_CABLE_HYDROPROP_INDEX;

                data_in.cable_maps->at(index).i_hydro_property.resize(number_string.size());
                for (int i=0; i<number_string.size(); i++)
                    data_in.cable_maps->at(index).i_hydro_property(i)
                    = stoi(number_string[i]);

                // Seabed property indexes.
                number_string.clear();
                flag = extract_multi_number(subnode->first_node("iseabedproperty")->value(), number_string);

                if (!flag)
                    return ErrorCode::MOORING_FILE_NAN_CABLE_SEABEDPROP_INDEX;

                data_in.cable_maps->at(index).i_seabed_property.resize(number_string.size());
                for (int i=0; i<number_string.size(); i++)
                    data_in.cable_maps->at(index).i_seabed_property(i)
                    = stoi(number_string[i]);

                // Connection indexes.
                number_string.clear();
                flag = extract_multi_number(subnode->first_node("iconnection")->value(), number_string);

                if (!flag)
                    return ErrorCode::MOORING_FILE_NAN_CABLE_CONNECTION_INDEX;

                if (number_string.size() != 2)
                    return ErrorCode::MOORING_FILE_BAD_CABLE_CONNECTION_INDEX;

                for (int i=0; i<number_string.size(); i++)
                    data_in.cable_maps->at(index).i_connection(i)
                    = stoi(number_string[i]);
            }

            int crit = 1;
            for (int i=0; i<index_count.rows(); i++)
                crit = (crit && (index_count(i) == 1));
            if (!crit)
                return ErrorCode::MOORING_FILE_MISSED_REPEATED_CABLE_ID;
        }
        ////////////////////////////////////////////////////////////////////////
        // Structural properties input data.
        ////////////////////////////////////////////////////////////////////////
        else if (strcmp(child_node->name(), "structuralproperties") == 0)
        {
            data_in.n_struct_property = read_number(child_node);
            if (data_in.n_struct_property == -1)
                return ErrorCode::MOORING_FILE_BAD_STRUCTPROPS_NUM;
            
            data_in.struct_properties->resize(data_in.n_struct_property);
            if (child_node->first_node("structuralproperty")==0)
                return ErrorCode::MOORING_FILE_NO_STRUCTPROP;

            index_count.resize(data_in.n_struct_property);
            index_count.setConstant(0);

            for (subnode = child_node->first_node("structuralproperty"); subnode;
                 subnode = subnode->next_sibling("structuralproperty"))
            {
                vector<string> names;
                names.push_back("id");
                names.push_back("diameter");
                names.push_back("unitlengthmass");
                names.push_back("unitlengthweight");
                names.push_back("axialstiffness");
                names.push_back("bendingstiffness");
                names.push_back("torsionalstiffness");
                names.push_back("dampingcoefficient");

                if (!check_availability(subnode, names))
                    return ErrorCode::MOORING_FILE_INCOMPLETE_STRUCTPROP_DEF;

                if (!check_is_number(subnode, names))
                    return ErrorCode::MOORING_FILE_NAN_STRUCTPROP;

                index = stoi(subnode->first_attribute("id")->value());

                if (index < 0 || index >= data_in.n_struct_property)
                    return ErrorCode::MOORING_FILE_OUTOFRANGE_STRUCTPROP_ID;

                index_count(index) = index_count(index) + 1;

                data_in.struct_properties->at(index).diameter
                = stod(subnode->first_node("diameter")->value());

                data_in.struct_properties->at(index).unit_length_mass
                = stod(subnode->first_node("unitlengthmass")->value());

                data_in.struct_properties->at(index).unit_length_weight
                = stod(subnode->first_node("unitlengthweight")->value());

                data_in.struct_properties->at(index).axial_stiffness
                = stod(subnode->first_node("axialstiffness")->value());

                data_in.struct_properties->at(index).bending_stiffness
                = stod(subnode->first_node("bendingstiffness")->value());

                data_in.struct_properties->at(index).torsional_stiffness
                = stod(subnode->first_node("torsionalstiffness")->value());

                data_in.struct_properties->at(index).damping_coefficient
                = stod(subnode->first_node("dampingcoefficient")->value());
            }

            int crit = 1;
            for (int i=0; i<index_count.rows(); i++)
                crit = (crit && (index_count(i) == 1));
            if (!crit)
            {
                printf("Missed or repearted structuralproperty definition.");
                return ErrorCode::MOORING_FILE_MISSED_REPEATED_STRUCTPROP_ID;
            }
        }
        ////////////////////////////////////////////////////////////////////////
        // Hydro properties input data.
        ////////////////////////////////////////////////////////////////////////
        else if (strcmp(child_node->name(), "hydroproperties") == 0)
        {
            data_in.n_hydro_property = read_number(child_node);
            if (data_in.n_hydro_property == -1)
                return ErrorCode::MOORING_FILE_BAD_HYDROPROPS_NUM;
            
            data_in.hydro_properties->resize(data_in.n_hydro_property);

            index_count.resize(data_in.n_hydro_property);
            index_count.setConstant(0);

            if (child_node->first_node("hydroproperty")==0)
                return ErrorCode::MOORING_FILE_NO_HYDROPROP;

            for (subnode = child_node->first_node("hydroproperty"); subnode;
                 subnode = subnode->next_sibling("hydroproperty"))
            {
                vector<string> names;
                names.push_back("id");
                names.push_back("addedmasscoefficient");
                names.push_back("dragcoefficient");

                if (!check_availability(subnode, names))
                    return ErrorCode::MOORING_FILE_INCOMPLETE_HYDROPROP_DEF;

                if (!is_number(subnode->first_attribute("id")->value()))
                    return ErrorCode::MOORING_FILE_NAN_HYDROPROP;

                index = stoi(subnode->first_attribute("id")->value());

                if (index < 0 || index >= data_in.n_hydro_property)
                    return ErrorCode::MOORING_FILE_OUTOFRANGE_HYDROPROP_ID;

                index_count(index) = index_count(index) + 1;

                sub2node = subnode->first_node("addedmasscoefficient");

                names.clear();
                names.push_back("tangential");
                names.push_back("normal");
                names.push_back("binormal");

                if (!check_availability(sub2node, names))
                    return ErrorCode::MOORING_FIEL_INCOMPLETE_HYDRO_COEFFICIENTS;
                
                if (!check_is_number(sub2node, names))
                    return ErrorCode::MOORING_FIEL_NAN_HYDRO_COEFFICIENTS;

                data_in.hydro_properties->at(index).added_mass_coefficient
                <<  stod(sub2node->first_attribute("tangential")->value()),
                    stod(sub2node->first_attribute("normal")->value()),
                    stod(sub2node->first_attribute("binormal")->value());

                sub2node = subnode->first_node("dragcoefficient");

                if (!check_availability(sub2node, names))
                    return ErrorCode::MOORING_FIEL_INCOMPLETE_HYDRO_COEFFICIENTS;

                if (!check_is_number(sub2node, names))
                    return ErrorCode::MOORING_FIEL_NAN_HYDRO_COEFFICIENTS;
                
                sub2node = subnode->first_node("dragcoefficient");
                data_in.hydro_properties->at(index).drag_coefficient
                <<  stod(sub2node->first_attribute("tangential")->value()),
                    stod(sub2node->first_attribute("normal")->value()),
                    stod(sub2node->first_attribute("binormal")->value());
            }

            int crit = 1;
            for (int i=0; i<index_count.rows(); i++)
                crit = (crit && (index_count(i) == 1));
            if (!crit)
                return ErrorCode::MOORING_FILE_MISSED_REPEATED_HYDROPROP_ID;
        }
        ////////////////////////////////////////////////////////////////////////
        // Seabed properties input data.
        ////////////////////////////////////////////////////////////////////////
        else if (strcmp(child_node->name(), "seabedproperties") == 0)
        {
            data_in.n_seabed_property = read_number(child_node);
            if (data_in.n_seabed_property == -1)
                return ErrorCode::MOORING_FILE_BAD_SEABEDPROPS_NUM;

            data_in.seabed_properties->resize(data_in.n_seabed_property);
            if (child_node->first_node("seabedproperty")==0)
                return ErrorCode::MOORING_FILE_NO_SEABEDPROP;

            index_count.resize(data_in.n_seabed_property);
            index_count.setConstant(0);

            for (subnode = child_node->first_node("seabedproperty"); subnode;
                 subnode = subnode->next_sibling("seabedproperty"))
            {
                vector<string> names;
                names.push_back("id");
                names.push_back("dampingcoefficient");
                names.push_back("stiffnesscoefficient");

                if (!check_availability(subnode, names))
                    return ErrorCode::MOORING_FILE_INCOMPLETE_SEABEDPROP_DEF;

                if (!check_is_number(subnode, names))
                    return ErrorCode::MOORING_FILE_NAN_SEABEDPROP;

                index = stoi(subnode->first_attribute("id")->value());

                if (index < 0 || index >= data_in.n_seabed_property)
                    return ErrorCode::MOORING_FILE_OUTOFRANGE_SEABEDPROP_ID;

                index_count(index) = index_count(index) + 1;

                data_in.seabed_properties->at(index).damping_coefficient
                = stod(subnode->first_node("dampingcoefficient")->value());
                
                data_in.seabed_properties->at(index).stiffness_coefficient
                = stod(subnode->first_node("stiffnesscoefficient")->value());
            }

            int crit = 1;
            for (int i=0; i<index_count.rows(); i++)
                crit = (crit && (index_count(i) == 1));
            if (!crit)
                return ErrorCode::MOORING_FILE_MISSED_REPEATED_SEABEDPROP_ID;
        }
        ////////////////////////////////////////////////////////////////////////
        // Currents input data.
        ////////////////////////////////////////////////////////////////////////
        else if (strcmp(child_node->name(), "currents") == 0)
        {
            data_in.n_current = read_number(child_node);
            if (data_in.n_current == -1)
                return ErrorCode::MOORING_FILE_BAD_CURRENTS_NUM;

            data_in.currents->resize(data_in.n_current);
            if (child_node->first_node("current")==0)
                return ErrorCode::MOORING_FILE_NO_CURRENT;

            index_count.resize(data_in.n_current);
            index_count.setConstant(0);

            for (subnode = child_node->first_node("current"); subnode;
                 subnode = subnode->next_sibling("current"))
            {
                vector<string> names;
                names.push_back("id");
                names.push_back("polyorder");
                names.push_back("profilefile");


                if (!check_availability(subnode, names))
                    return ErrorCode::MOORING_FILE_INCOMPLETE_CURRENT_DEF;

                names.erase(names.end()-1);

                if (!check_is_number(subnode, names))
                    return ErrorCode::MOORING_FILE_NAN_CURRENT;

                index = stoi(subnode->first_attribute("id")->value());

                if (index < 0 || index >= data_in.n_current)
                    return ErrorCode::MOORING_FILE_OUTOFRANGE_CURRENT_ID;

                index_count(index) = index_count(index) + 1;

                data_in.currents->at(index).polyfit_order =
                stoi(subnode->first_node("polyorder")->value());

                data_in.currents->at(index).profile_file
                = work_folder + subnode->first_node("profilefile")->value();
            }

            int crit = 1;
            for (int i=0; i<index_count.rows(); i++)
                crit = (crit && (index_count(i) == 1));
            if (!crit)
                return ErrorCode::MOORING_FILE_MISSED_REPEATED_CURRENT_ID;

            for (int i=0; i<data_in.currents->size(); i++)
            {
                int flag = check_file_existence(data_in.currents->at(i).profile_file);
                if (flag)
                    return ErrorCode::CURRENT_FILE_NONEXISTENT;
                
                int failure = read_to(data_in.currents->at(i).profile_data,
                                      data_in.currents->at(i).profile_file,
                                      6, 1);
                if (failure)
                    return ErrorCode::MOORING_FILE_BAD_CURRENT_DATA;
            }
        }
        ////////////////////////////////////////////////////////////////////////
        // Solver input data.
        ////////////////////////////////////////////////////////////////////////
        else if (strcmp(child_node->name(), "solvers") == 0)
        {
            data_in.n_solver = read_number(child_node);
            if (data_in.n_solver == -1)
                return ErrorCode::MOORING_FILE_BAD_SOLVERS_NUM;
            
            data_in.solvers->resize(data_in.n_solver);

            if (child_node->first_node("solver")==0)
                return ErrorCode::MOORING_FILE_NO_SOLVER;

            index_count.resize(data_in.n_solver);
            index_count.setConstant(0);

            for (subnode = child_node->first_node("solver"); subnode;
                 subnode = subnode->next_sibling("solver"))
            {
                vector<string> names;
                names.push_back("id");
                names.push_back("iterationnumberlimit");
                names.push_back("convergencetolerance");
                names.push_back("initialrelaxationfactor");
                names.push_back("relaxationincreasefactor");
                names.push_back("relaxationdecreasefactor");
                names.push_back("lambdainfinity");

                if (!check_availability(subnode, names))
                    return ErrorCode::MOORING_FILE_INCOMPLETE_SOLVER_DEF;

                if (!check_is_number(subnode, names))
                    return ErrorCode::MOORING_FILE_NAN_SOLVER;

                index = stoi(subnode->first_attribute("id")->value());

                if (index < 0 || index >= data_in.n_solver)
                    return ErrorCode::MOORING_FILE_OUTOFRANGE_SOLVER_ID;

                index_count(index) = index_count(index) + 1;

                data_in.solvers->at(index).n_iteration_limit =
                stoi(subnode->first_node("iterationnumberlimit")->value());

                data_in.solvers->at(index).convergence_tolerance =
                stod(subnode->first_node("convergencetolerance")->value());

                data_in.solvers->at(index).initial_relaxation_factor =
                stod(subnode->first_node("initialrelaxationfactor")->value());

                data_in.solvers->at(index).increment_factor =
                stod(subnode->first_node("relaxationincreasefactor")->value());

                data_in.solvers->at(index).decrement_factor =
                stod(subnode->first_node("relaxationdecreasefactor")->value());

                data_in.solvers->at(index).lambda_infinity =
                stod(subnode->first_node("lambdainfinity")->value());
            }

            int crit = 1;
            for (int i=0; i<index_count.rows(); i++)
                crit = (crit && (index_count(i) == 1));
            if (!crit)
                return ErrorCode::MOORING_FILE_MISSED_REPEATED_SOLVER_ID;
        }
    }
    return ErrorCode::SUCCESS;
}


////////////////////////////////////////////////////////////////////////////////
/// Initialization of the setting.
////////////////////////////////////////////////////////////////////////////////
ErrorCode Reader::read_to(Setting& setting)
{
    // Get the absolute path to Setting.xml file.
    std::string file_name = setting.moor_folder + setting.setting_file;
    
    ////////////////////////////////////////////////////////////////////////////
    // Setting file.
    ////////////////////////////////////////////////////////////////////////////
    int flag = check_file_existence(file_name);
    if (flag)
        return ErrorCode::SETTING_FILE_NONEXISTENT;
    
    ifstream file(file_name);
    std::stringstream buffer;
    buffer << file.rdbuf();
    file.close();
    
    try
    {
        xml_document<> set_doc;
        std::string content(buffer.str());
        set_doc.parse<0> (&content[0]);
    }
    catch (const rapidxml::parse_error& e)
    {
        return ErrorCode::SETTING_FILE_ERROR_PARSE;
    }
    
    xml_document<> set_doc;
    std::string content(buffer.str());
    set_doc.parse<0> (&content[0]);
    
    xml_node<> *root_node;
    root_node = set_doc.first_node("setting");
    if (root_node == 0)
        return ErrorCode::SETTING_FILE_NO_SETTING_NODE;
    
    ////////////////////////////////////////////////////////////////////////////
    // Simulation type.
    ////////////////////////////////////////////////////////////////////////////
    xml_node<> *child_node = root_node->first_node("simulationtype");
    if (child_node==0)
        return ErrorCode::SETTING_FILE_NO_SIMULATION_TYPE;
    
    std::string type_string = child_node -> value();

    if (type_string.compare("shooting") == 0)
        setting.simulation_type = SimulationType::SHOOTING;
    else if (type_string.compare("forcedmotion") == 0)
        setting.simulation_type = SimulationType::FORCED_MOTION;
    else if (type_string.compare("relaxation") ==0 )
        setting.simulation_type = SimulationType::RELAXATION;
    else
        return ErrorCode::SETTING_FILE_UNKNOWN_SIMULATION_TYPE;
    
    ////////////////////////////////////////////////////////////////////////////
    // Mooring input file.
    ////////////////////////////////////////////////////////////////////////////
    // Get the mooringinputfile.
    child_node = root_node->first_node("mooringinputfile");
    if (child_node==0)
        return ErrorCode::SETTING_FILE_NO_MOORING_FILE_DEF;
    
    // Whether the file exists has not been checked yet.
    std::string relative_path = child_node -> value();
    setting.mooring_file = setting.moor_folder + relative_path;
    
    flag = check_file_existence(setting.mooring_file);
    if (flag) return ErrorCode::MOORING_FILE_NONEXISTENT;
    
    size_t work_folder_index = setting.mooring_file.find_last_of("/\\");
    setting.work_folder = setting.mooring_file.substr(0, work_folder_index+1);
    
    // If this flag is not set properly, the default value is used.
    child_node = root_node->first_node("platformsaveflag");
    if (child_node != 0 && is_number(child_node->value()))
        setting.is_saving_platform = stoi(child_node->value());
    
    xml_node<>* subnode;
    std::string cell;

    ////////////////////////////////////////////////////////////////////////////
    // Parameters for specific simulation type.
    ////////////////////////////////////////////////////////////////////////////
    switch(setting.simulation_type)
    {
        ////////////////////////////////////////////////////////////////////////
        // Shooting parameters.
        ////////////////////////////////////////////////////////////////////////
        case SimulationType::SHOOTING:
        {
            child_node = root_node->first_node("shooting");
            if (child_node == 0)
                return ErrorCode::SETTING_FILE_NO_SHOOTING_PARA;
            
            std::vector<std::string> para_names;
            para_names.push_back("fairleadpositiontolerance");
            para_names.push_back("fairleadforcerelaxationfactor");
            para_names.push_back("fairleadpositioniterationlimit");
            para_names.push_back("platformpositioniterationlimit");
            para_names.push_back("platformdisplacementtolerance");
            para_names.push_back("platformdisplacementrelaxationfactor");
            para_names.push_back("cableoutofplanestiffness");
            para_names.push_back("platformhydrostaticstiffness");
            para_names.push_back("platformotherload");
            
            if (!check_availability(child_node, para_names))
                return ErrorCode::SETTING_FILE_INCOMPLETE_SHOOTING_PARA;
            else
            {
                para_names.erase(para_names.end()-2, para_names.end());
                para_names.erase(para_names.begin());
                if (!check_is_number(child_node, para_names))
                    return ErrorCode::SETTING_FILE_NAN_SHOOTING_PARA;
            }
            
            subnode = child_node->first_node("fairleadpositiontolerance");
            para_names.clear();
            para_names.push_back("xc");
            para_names.push_back("yc");
            para_names.push_back("zc");
            if (!check_availability(subnode, para_names))
                return ErrorCode::SETTING_FILE_INCOMPLETE_SHOOTING_FAIRLEAD_TOLERANCE;
            else if (!check_is_number(subnode, para_names))
                return ErrorCode::SETTING_FILE_NAN_SHOOTING_FAIRLEAD_TOLERANCE;
            
            setting.shooting->fairlead_position_tolerance
            << stod(subnode->first_attribute("xc")->value()),
               stod(subnode->first_attribute("yc")->value()),
               stod(subnode->first_attribute("zc")->value());

            setting.shooting->fairlead_force_relaxation_factor =
            stod(child_node->first_node("fairleadforcerelaxationfactor")->value());

            setting.shooting->fairlead_position_iteration_limit =
            stoi(child_node->first_node("fairleadpositioniterationlimit")->value());

            setting.shooting->platform_position_iteration_limit =
            stoi(child_node->first_node("platformpositioniterationlimit")->value());

            setting.shooting->platform_displacement_tolerance =
            stod(child_node->first_node("platformdisplacementtolerance")->value());

            setting.shooting->platform_displacement_relaxation_factor =
            stod(child_node->first_node("platformdisplacementrelaxationfactor")->value());

            setting.shooting->cable_out_of_plane_stiffness =
            stod(child_node->first_node("cableoutofplanestiffness")->value());

            
            // Connection indexes.
            std::vector<std::string> number_string;
            int flag;
            flag = extract_multi_number(child_node->first_node("platformhydrostaticstiffness")->value(), number_string);
            
            if (!flag || number_string.size() != 36)
                return ErrorCode::SETTING_FILE_BAD_SHOOTING_STIFFNESS;
            else
            {
                for (int i=0; i<6; i++)
                    for (int j=0; j<6; j++)
                        setting.shooting->platform_hydrostatic_stiffness(i,j)
                        = stod(number_string[6*i+j]);
            }
            
            number_string.clear();
            flag = extract_multi_number(child_node->first_node("platformotherload")->value(), number_string);
            if (!flag || number_string.size() != 6)
                return ErrorCode::SETTING_FILE_BAD_SHOOTING_LOAD;
            else
            {
                for (int i=0; i<6; i++)
                    setting.shooting->platform_other_load(i)= stod(number_string[i]);
            }
        }
            break;
        
        ////////////////////////////////////////////////////////////////////////
        // Relaxation parameters.
        ////////////////////////////////////////////////////////////////////////
        case SimulationType::RELAXATION:
        {
            child_node = root_node->first_node("relaxation");
            if (child_node==0)
                return ErrorCode::SETTING_FILE_NO_RELAXATION_PARA;
            
            std::vector<std::string> para_names;
            para_names.push_back("platformvelocitytolerance");
            para_names.push_back("cablevelocitytolerance");
            para_names.push_back("stoptime");
            para_names.push_back("timestep");
            para_names.push_back("platformmass");
            para_names.push_back("platformdamping");
            para_names.push_back("platformhydrostaticstiffness");
            para_names.push_back("platformotherload");
            
            if (!check_availability(child_node, para_names))
                return ErrorCode::SETTING_FILE_INCOMPLETE_RELAXATION;
            else
            {
                para_names.erase(para_names.end()-4, para_names.end());
                if (!check_is_number(child_node, para_names))
                    return ErrorCode::SETTING_FILE_NAN_RELAXATION_PARA;
            }
            
            setting.relaxation->platform_velocity_tolerance =
            stod(child_node->first_node("platformvelocitytolerance")->value());
            
            setting.relaxation->cable_velocity_tolerance =
            stod(child_node->first_node("cablevelocitytolerance")->value());
            
            setting.relaxation->stop_time =
            stod(child_node->first_node("stoptime")->value());
            
            setting.relaxation->time_step =
            stod(child_node->first_node("timestep")->value());
            
            std::vector<std::string> number_string;
            int flag;
            flag = extract_multi_number(child_node->first_node("platformmass")->value(), number_string);
            
            if (!flag || number_string.size() != 36)
                return ErrorCode::SETTING_FILE_BAD_RELAXATION_MASS;
            else
            {
                for (int i=0; i<6; i++)
                    for (int j=0; j<6; j++)
                        setting.relaxation->platform_mass(i,j) = stod(number_string[6*i+j]);
            }
            
            number_string.clear();
            flag = extract_multi_number(child_node->first_node("platformdamping")->value(), number_string);
            
            if (!flag || number_string.size() != 36)
                return ErrorCode::SETTING_FILE_BAD_RELAXATION_DAMPING;
            else
            {
                for (int i=0; i<6; i++)
                    for (int j=0; j<6; j++)
                        setting.relaxation->platform_damping(i,j)
                        = stod(number_string[6*i+j]);
            }
            
            number_string.clear();
            flag = extract_multi_number(child_node->first_node("platformhydrostaticstiffness")->value(), number_string);
            
            if (!flag || number_string.size() != 36)
                return ErrorCode::SETTING_FILE_BAD_RELAXATION_STIFFNESS;
            else
            {
                for (int i=0; i<6; i++)
                    for (int j=0; j<6; j++)
                        setting.relaxation->platform_hydrostatic_stiffness(i,j)
                        = stod(number_string[6*i+j]);
            }
            
            number_string.clear();
            flag = extract_multi_number(child_node->first_node("platformotherload")->value(), number_string);
            if (!flag || number_string.size() != 6)
                return ErrorCode::SETTING_FILE_BAD_RELAXATION_LOAD;
            else
            {
                for (int i=0; i<6; i++)
                    setting.relaxation->platform_other_load(i)= stod(number_string[i]);
            }
            
        }
            break;

        ////////////////////////////////////////////////////////////////////////
        // Forcedmotion parameters.
        ////////////////////////////////////////////////////////////////////////
        case SimulationType::FORCED_MOTION:
        {
            child_node = root_node->first_node("forcedmotion");
            if (child_node==0)
                return ErrorCode::SETTING_FILE_NO_MOTION;

            if (child_node->first_node("timehistory")==0)
                return ErrorCode::SETTING_FILE_NO_TIME_HISTORY;
            else
            {
                setting.forced_motion->data_file
                = setting.work_folder + child_node->first_node("timehistory")->value();
            }
            
            int err_code = check_file_existence(setting.forced_motion->data_file);
            if (err_code)
                return ErrorCode::MOTION_FILE_NONEXISTENT;
            
            err_code = read_to(setting.forced_motion->time_history,
                               setting.forced_motion->data_file, 7, 1);
            
            if (err_code)
                return ErrorCode::MOTION_FILE_BAD_DATA;
            
            setting.forced_motion->n_time
            = setting.forced_motion->time_history.rows();
        }
            break;
    }
    return ErrorCode::SUCCESS;
}


////////////////////////////////////////////////////////////////////////////////
/// Check whether file exist.
////////////////////////////////////////////////////////////////////////////////
int Reader::check_file_existence(const std::string file_name)
{
    ifstream file(file_name);
    if (!file.good() || (file_name.back() == '/' || file_name.back() == '\\'
                         || file_name.back() == '.'))
        return 1; // Open file failed.
    else
    {
        file.close();
        return 0;
    }
}

////////////////////////////////////////////////////////////////////////////////
/// Read data matrix with header lines which applies to current profile data,
/// cable, initial state data, and excitation data.
////////////////////////////////////////////////////////////////////////////////
int Reader::read_to(MatrixXd& data_mat, const std::string data_file_name,
                    const int expected_cols, const int skip_rows)
{
    string line, cell;
    int i_line = 0, i_cell;
    ifstream data_file(data_file_name);
    
    if (!data_file.good() || (data_file_name.back() == '/'
                              || data_file_name.back() == '\\'
                              || data_file_name.back() == '.'))
    {
        data_file.close();
        return 1; // Open file failed.
    }
    else
    {
        // Get number of lines.
        int i_line = 0;
        while (std::getline(data_file, line))
            i_line++;
        
        int n_points = i_line - 1;
        data_mat.resize(n_points, expected_cols);
        
        data_file.close();
        data_file.open(data_file_name);
        
        i_line = 0;
        while (std::getline(data_file, line))
        {
            std::stringstream line_stream(line);
            if (i_line >= skip_rows)
            {
                i_cell = 0;
                while (line_stream >> cell)
                {
                    if (!(i_cell < expected_cols))
                        return 3; // More than expected columns found.
                    else if (!is_number(cell))
                        return 4; // NaN found in data.
                    else
                    {
                        data_mat(i_line - skip_rows, i_cell) = stod(cell);
                    }
                    i_cell++;
                }
                if (i_cell != expected_cols)
                    return 2;
            }
            i_line++;
        }
        data_file.close();
        return 0; // Success.
    }
}


////////////////////////////////////////////////////////////////////////////////
/// Read the number of each component in the input file.
////////////////////////////////////////////////////////////////////////////////
int Reader::read_number(const xml_node<>* node)
{
    if (node->first_attribute("number")==0 ||
        !is_number(node->first_attribute("number")->value()))
        return -1;
    else
    {
        return (stoi(node->first_attribute("number")->value()) > 0 ?
                stoi(node->first_attribute("number")->value()) : -1);
    }
}


////////////////////////////////////////////////////////////////////////////////
/// Check whether expected subnodes or attributes are available.
////////////////////////////////////////////////////////////////////////////////
int Reader::check_availability(const rapidxml::xml_node<>* node,
                               const std::vector<std::string>& names)
{
    int success = 1;
    for (int i=0; i<names.size(); i++)
    {
        success = (success && (node->first_node(names[i].c_str()) !=0
                               || node->first_attribute(names[i].c_str()) !=0));
    }
    return success;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if values of a group of nodes or attributes are numbers.
////////////////////////////////////////////////////////////////////////////////
int Reader::check_is_number(const rapidxml::xml_node<>* node,
                            const std::vector<std::string>& names)
{
    int success = 1;
    for (int i=0; i<names.size(); i++)
    {
        if (node->first_node(names[i].c_str()))
            success = (success && (is_number(node->first_node(names[i].c_str())->value())));
        else if (node->first_attribute(names[i].c_str()))
            success = success && (is_number(node->first_attribute(names[i].c_str())->value()));
    }
    return success;
}
    
////////////////////////////////////////////////////////////////////////////////
/// Check whether each of the string in a whitespace separated text is a number.
////////////////////////////////////////////////////////////////////////////////
int Reader::extract_multi_number(const std::string token,
                                 std::vector<std::string>& number_string)
{
    std::stringstream line_stream;
    line_stream << token;
    number_string.clear();
    string cell;
    int success = 1;
    while (line_stream >> cell)
    {
        number_string.push_back(cell);
        success = success && is_number(cell);
    }
    return success;
}

} // End of namespace moor.
