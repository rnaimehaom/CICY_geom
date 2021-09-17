//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Calculate electrical resistivity of a network read from a file (csv or binary)
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#ifndef APP_NETWORK_FROM_FILE_H
#define APP_NETWORK_FROM_FILE_H

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


using namespace hns;

//---------------------------------------------------------------------------
class App_Network_From_Abaqus
{
public:

	int Generate_nanoparticle_resistor_network_from_file(Input* Init)const;
};
//---------------------------------------------------------------------------
#endif
//===========================================================================