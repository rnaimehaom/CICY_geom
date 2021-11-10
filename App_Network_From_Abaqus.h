//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Calculate electrical resistivity of a network read from an Abaqus output file (.odb)
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#ifndef APP_NETWORK_FROM_ABAQUS_H
#define APP_NETWORK_FROM_ABAQUS_H

#include "time.h"
#include "Hns.h"
#include "Backbone_Network.h"
#include "Shells.h"
#include "Contact_grid.h"
#include "Cutoff_Wins.h"
#include "Direct_Electrifying.h"
#include "Electrical_analysis.h"
#include "Hoshen_Kopelman.h"
#include "Input_Reader.h"
#include "Read_Network.h"
#include "GNP_reconstruction.h"
//Include for Abaqus
#include <odb_API.h>

using namespace hns;

//---------------------------------------------------------------------------
class App_Network_From_Abaqus
{
public:

    //Constructor
    App_Network_From_Abaqus() {};

    //Member Functions
    int Nanoparticle_resistor_network_from_odb(Input* Init)const;
    int Get_gnps_partially_outside_sample(const Geom_sample& geom_sample, const vector<GNP>& gnps, vector<int> &gnps_outside, vector<vector<int> >& vertices_gnps_out, vector<vector<bool> > &vertex_flags)const;
    //Generate a network of nanoparticles from the data in the Abaqus database
    int Apply_displacements_from_Abaqus(const string& particle_type, const int& n_cnts, const vector<vector<long int> >& structure, const vector<int>& gnps_outside, const vector<vector<int> >& vertices_gnps_out, const vector<vector<bool> >& vertex_flags, odb_Assembly& root_assy, odb_Frame& previous_frame, odb_Frame& current_frame, Geom_sample& geom_sample, vector<Point_3D>& points_cnt, vector<GNP>& gnps)const;
    int Apply_displacements_to_sample(odb_Assembly& root_assy, odb_FieldOutput& previous_fieldU, odb_FieldOutput& current_fieldU, Geom_sample& geom_sample)const;
    int Get_displacement_change_from_single_node_set(const string& setname, odb_Assembly& root_assy, odb_FieldOutput& previous_fieldU, odb_FieldOutput& current_fieldU, Point_3D& disp)const;
    int Get_displacement_from_single_node_set(const string& setname, odb_Assembly& root_assy, odb_FieldOutput& fieldU, Point_3D& disp)const;
    int Apply_displacements_to_cnts(const vector<vector<long int> >& structure, odb_Assembly& root_assy, odb_FieldOutput& previous_fieldU, odb_FieldOutput& current_fieldU, const int& n_cnts, vector<Point_3D>& points_cnt)const;
    string Get_cnt_set_name(const int& cnt_i)const;
    int Apply_displacements_to_gnps(const vector<int>& gnps_outside, const vector<vector<int> >& vertices_gnps_out, const vector<vector<bool> >& vertex_flags, odb_Assembly& root_assy, odb_FieldOutput& previous_fieldU, odb_FieldOutput& current_fieldU, vector<GNP>& gnps)const;
    vector<int> All_gnp_vertices()const;
    int Apply_displacements_to_gnp_vertices(const vector<int>& vertices, odb_Assembly& root_assy, odb_FieldOutput& previous_fieldU, odb_FieldOutput& current_fieldU, GNP& gnp_i)const;
    string Get_gnp_set_name(const int& gnp_i, const int& vertex)const;
};
//---------------------------------------------------------------------------
#endif
//===========================================================================