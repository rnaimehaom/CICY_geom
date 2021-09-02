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
using namespace hns;

#include "Backbone_Network.h"
#include "Shells.h"
#include "Contact_grid.h"
#include "Cutoff_Wins.h"
#include "Direct_Electrifying.h"
#include "Electrical_analysis.h"
#include "Hoshen_Kopelman.h"
#include "Input_Reader.h"

//---------------------------------------------------------------------------
class App_Network_From_Abaqus
{
public:

    //Constructor
    App_Network_From_Abaqus() {};

    //Member Functions
    int Nanoparticle_resistor_network_from_odb(Input* Init)const;
};
//---------------------------------------------------------------------------
#endif
//===========================================================================