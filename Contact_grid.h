//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GS networks
//OBJECTIVE:    Generate a grid to facilitate findic contact points
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#ifndef CONTACT_GRID_H
#define CONTACT_GRID_H

#include <iostream>

#include "Input_Reader.h"
#include "Generate_Network.h"

//-------------------------------------------------------
class Contact_grid
{
public:
    
    //Sectioned domain for CNTs
    vector<vector< long int> > sectioned_domain_cnts;
    //Sectioned domain for GNP numbers
    vector<vector<int> > sectioned_domain_gnps;
    
    
    //Constructor
    Contact_grid(){};
    
    //Member Functions
    int Generate_contact_grid(const int &window, const string &particle_type, const Geom_sample &sample_geom, const cuboid &window_geom, const vector<int> &cnts_inside, vector<Point_3D> &points_cnt, const vector<vector<long int> > &structure, const vector<int> &gnps_inside, const vector<GNP> &gnps);
    int Adjust_regions_if_needed(const double &overlapping, const Geom_sample &sample_geom, const cuboid &window_geom, int n_regions[], double l_regions[]);
    int Fill_sectioned_domain_cnts(const cuboid &window_geom, const vector<int> &cnts_inside, const vector<vector<long int> > &structure, const vector<Point_3D> &points_cnt, const double &overlapping, const int n_regions[], const double l_regions[]);
    int Calculate_region_coordinates(const cuboid &window_geom, const Point_3D &point, const int n_regions[], const double l_regions[], int &a, int &b, int &c);
    int Calculate_overlapping_flags(const cuboid &window_geom, const Point_3D &point, const double &overlapping, const int &a, const int &b, const int &c, const int n_regions[], const double l_regions[], int f_regions[]);
    int Assign_point_to_regions_cnts(const int &a, const int &b, const int &c, const int f_regions[], const int n_regions[], const long int &P);
    int Calculate_t(const int &a, const int &b, const int &c, const int &sx, const int &sy);
    int Fill_sectioned_domain_gnps(const cuboid &window_geom, const vector<GNP> &gnps, const vector<int> &gnps_inside, const double &overlapping, const int n_regions[], const double l_regions[]);
    int Fill_sectioned_domain_single_gnp(const cuboid &window_geom, const GNP &gnp, const double &overlapping, const int n_regions[], const double l_regions[]);
    int Assign_point_to_regions_gnps(const int &a, const int &b, const int &c, const int f_regions[], const int n_regions[], const int &GNP);
    
private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================
