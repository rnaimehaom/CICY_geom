//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Calculate electrical resistivity of a network read from dissplacements saved in binary files (obtained from an Abaqus output odb file)
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#ifndef APP_NETWORK_FROM_DISPLACEMENTS_H
#define APP_NETWORK_FROM_DISPLACEMENTS_H

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

using namespace hns;

//---------------------------------------------------------------------------
class App_Network_From_Displacements
{
public:

    //Constructor
    class App_Network_From_Displacements () {};

    //Member Functions
    int Nanoparticle_resistor_network_from_displacements(Input* Init)const;
    int Get_gnps_partially_outside_sample(const Geom_sample& geom_sample, const vector<GNP>& gnps, vector<vector<int> >& vertices_gnps_in)const;
    int Apply_displacements_from_files(const int& frame, const string& particle_type, const vector<vector<long int> >& structure_cnt, const vector<vector<int> >& vertices_gnps_in, Geom_sample& geom_sample, vector<Point_3D>& points_cnt, vector<GNP>& gnps)const;
    int Apply_displacements_to_sample(const int& frame, Geom_sample& geom_sample)const;
    int Apply_incremental_displacements_to_sample(const int& frame, Geom_sample& geom_sample)const;
    double Get_difference(const streamsize& double_size, ifstream& previous, ifstream& current)const; 
    int Apply_displacements_to_cnts(const int& frame, const vector<vector<long int> >& structure_cnt, vector<Point_3D>& points_cnt)const;
    int Apply_displacements_to_gnps(const int& frame, const vector<vector<int> >& vertices_gnps_in, vector<GNP>& gnps)const;
    vector<int> All_gnp_vertices()const;
    int Apply_displacements_to_vertices(const int& frame, const streamsize& double_size, const vector<int>& vertices, GNP& gnp_i, ifstream& previous, ifstream& current)const;
};
//---------------------------------------------------------------------------
#endif
//===========================================================================
