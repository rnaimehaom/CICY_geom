//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Calculate the volume or weight fraction for each observation window
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================


#ifndef APP_CONTENT_DIST_H
#define APP_CONTENT_DIST_H

#include "time.h"
#include "Hns.h"
using namespace hns;

#include "Shells.h"
#include "Cutoff_Wins.h"
#include "Generate_Network.h"
#include "Input_Reader.h"

//---------------------------------------------------------------------------
class App_Content_Dist
{
public:
    
    //Constructor
    App_Content_Dist(){};
    
    //Member Functions
    int Calculate_content_on_each_window(Input *Init)const;
    int Calculate_content_in_window(const cuboid &window_geo, const string &particle_type, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<vector<long int> > &structure, const vector<int> &cnts_inside, const vector<GNP> &gnps, const vector<int> &gnps_inside)const;
};

#endif
