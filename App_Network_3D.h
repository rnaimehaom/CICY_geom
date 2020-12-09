//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Create 3D nanoparticle network. Find backbone, dead branches and isolated clusters. Calculate electrical conductivity of sample (and observation windows)
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#ifndef APPNETWORK3D_H
#define APPNETWORK3D_H

#include "time.h"
#include "Hns.h"
using namespace hns;

#include "Backbone_Network.h"
#include "Shells.h"
#include "Contact_grid.h"
#include "Cutoff_Wins.h"
#include "Direct_Electrifying.h"
#include "Electrical_analysis.h"
#include "Generate_Network.h"
#include "Hoshen_Kopelman.h"
#include "Input_Reader.h"

//---------------------------------------------------------------------------
class App_Network_3D
{
public:
    
    //Constructor
    App_Network_3D(){};
    
    //Member Functions
    int Generate_nanoparticle_resistor_network(Input *Init)const;
    int Update_obseravtion_window_geometry(const int &window, const int &window_domain, const Geom_sample &sample_geo, cuboid &window_geo)const;
};
//---------------------------------------------------------------------------
#endif
//===========================================================================

