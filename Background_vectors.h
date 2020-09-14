//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GS networks
//OBJECTIVE:    Generate shells (background vectors) to map each CNT into an observation window and facilitate their trimming
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#ifndef BACKGROUND_VECTORS_H
#define BACKGROUND_VECTORS_H

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
#include "GenNetwork.h"

//-------------------------------------------------------
class Background_vectors
{
public:
    //Data Member
    
    //Constructor
    Background_vectors(){};
    
    //Member Functions
    int Generate_shells(const struct Geom_sample &sample, const struct Nanotube_Geo &cnts, const vector<Point_3D> &points_in, const vector<GCH> &hybrid_particles, vector<vector<int> > &shells_cnt, vector<vector<int> > &shells_gnps);
    int Add_to_shell(const double midpoints[], const double boundary_layer[], const double core[], const double half_step[], const Point_3D &point, const int &n_shells, vector<vector<int> > &shells_cnt);
    int Find_minimum_shell(const double midpoints[], const double boundary_layer[], const double core[], const double half_step[], const Point_3D &point, const int &n_shells);
    int Find_shell(const double &x_in, const double &x_m, const double &x_layer, const double &x_core, const double &dx_half, const int &n_shells);
    int Add_to_shells(const double midpoints[], const double boundary_layer[], const double core[], const double half_step[], const GCH &hybrid, const int &n_shells, vector<vector<int> > &shells_gnp);
    
private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================
