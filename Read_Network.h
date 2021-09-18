//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Read 3D nanoparticle network from file
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#ifndef READ_NETWORK_H
#define READ_NETWORK_H

#include<iomanip>
#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>
#include "Input_Reader.h"
#include "Geometry_3D.h"
#include "Hns.h"
#include "Generate_Network.h"

using namespace hns;

class Read_Network
{

public:

    //Constructor
    Read_Network() {};

    //Member Functions
    int Generate_nanoparticle_network_from_file(const Simu_para& simu_para, const Visualization_flags& vis_flags, const Output_data_flags& out_flags, Geom_sample& geom_sample, vector<Point_3D>& points_cnt, vector<double>& radii, vector<vector<long int> >& structure, vector<GNP>& gnps)const;
    int Read_sample_geometry(Geom_sample& geom_sample)const;
    int Read_cnt_data_from_csv(vector<Point_3D>& points_cnt, vector<double>& radii, vector<vector<long int> >& structure)const;
    int Read_cnt_data_from_dat(vector<Point_3D>& points_cnt, vector<double>& radii, vector<vector<long int> >& structure)const;
    int Read_gnp_data_from_csv(const cuboid& sample_geom, vector<GNP>& gnps)const;
};
//---------------------------------------------------------------------------
#endif
//===========================================================================

