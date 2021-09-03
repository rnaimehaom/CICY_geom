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
#include "Network_From_Abaqus.h"
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
    //Generate a network of nanoparticles from the data in the Abaqus database
    int Generate_nanoparticle_network_from_file(const Simu_para& simu_para, const Visualization_flags& vis_flags, Geom_sample& geom_sample, vector<Point_3D>& points_cnt, vector<double>& radii, vector<vector<long int> >& structure, vector<GNP>& gnps)const;
    int Read_cnt_data_from_csv(vector<Point_3D>& points_cnt, vector<double>& radii, vector<vector<long int> >& structure)const;
    int Read_gnp_data_from_csv(vector<GNP>& gnps)const;
    int Read_sample_geometry(Geom_sample& geom_sample)const;
};
//---------------------------------------------------------------------------
#endif
//===========================================================================