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
    int Add_to_shell(const struct Geom_sample &sample, const Point_3D &point, vector<vector<int> > &shells_cnt);
    int Find_minimum_shell(const struct Geom_sample &sample, const Point_3D &point, const int &num_shells);
    int Find_shell(const double &x_in, const double &x_0, const double &len_x, const double &dx, const double &win_min_x, const double &win_max_x, const int &num_shells);
    int Add_to_shells(const struct Geom_sample &sample, const GCH &hybrid, vector<vector<int> > &shells_gnp);
    
private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================
