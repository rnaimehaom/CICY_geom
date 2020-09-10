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

//-------------------------------------------------------
class Contact_grid
{
public:
    //Data Member
    vector<vector< long int> > sectioned_domain; //Sectioned domain for CNTs
    vector<vector< long int> > sectioned_domain_gnps; //Sectioned domain for GNPs
    vector<vector<int> > sectioned_domain_hyb; //Sectioned domain for GNP numbers
    
    //Constructor
    Contact_grid(){};
    
    //Member Functions
    int Generate_contact_grid(const int &window, const string &particle_type, const struct Geom_sample &sample, const struct Cutoff_dist &cutoffs, const struct Nanotube_Geo &cnts, const vector<int> &cnts_inside, vector<Point_3D> &points_in, const vector<vector<long int> > &structure, const vector<int> &gnps_inside, const vector<Point_3D> &points_gnp, const vector<vector<long int> > &structure_gnp);
    int Generate_window_geometry(const int &window, const struct Geom_sample &sample, struct Geom_sample &window_geom);
    int Adjust_regions_if_needed(const double &cutoff, struct Geom_sample &window_geom, int &sx, int &sy, int &sz);
    int Fill_sectioned_domain(const struct Geom_sample &window_geom, const vector<int> &particles_inside, const vector<vector<long int> > &structure, const vector<Point_3D> &points_in, const double &cutoff, const int &sx, const int &sy, const int &sz, vector<vector< long int> > &sectioned_domain);
    int Calculate_region_coordinates(const struct Geom_sample &window_geom, const Point_3D &point, const int &sx, const int &sy, const int &sz, int &a, int &b, int &c);
    int Calculate_postion_flags(const struct Geom_sample &window_geom, const Point_3D &point, const double &cutoff, const int &a, const int &b, const int &c, const int &sx, const int &sy, const int &sz, int &fx, int &fy, int &fz);
    int Assign_point_to_region(const int &a, const int &b, const int &c, const int &fx, const int &fy, const int &fz, const int &sx, const int &sy, const long int &P, vector<vector< long int> > &sectioned_domain);
    int Calculate_t(const int &a, const int &b, const int &c, const int &sx, const int &sy);
    int Fill_sectioned_domain(const struct Geom_sample &window_geom, const vector<int> &hybs_inside, const vector<vector<long int> > &structure, const vector<Point_3D> &points_in, const double &cutoff, const int &sx, const int &sy, const int &sz, vector<vector<int> > &sectioned_domain_hyb);
    
    
private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================
