//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GS networks
//OBJECTIVE:    Generate 3D nanoparticle network
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#ifndef GENERATE_NETWORK_H
#define GENERATE_NETWORK_H

#include<iomanip>
#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>
#include <set>
#include <random> 
#include "Input_Reader.h"
#include "Geometry_3D.h"
#include "MathMatrix.h"
#include "Hns.h"
#include "Collision_detection.h"
#include "VTK_Export.h"

using namespace hns;

//-------------------------------------------------------
class Generate_Network
{
public:
    //Data Member
    
    //Constructor
    Generate_Network(){};
    
    //Generate a network of nanoparticles
    int Generate_nanoparticle_network(const Simu_para &simu_para, const Geom_sample &geom_sample, const Nanotube_Geo &nanotube_geo, const GNP_Geo &gnp_geo, const Cutoff_dist &cutoffs, const Tecplot_flags &tec360_flags, vector<Point_3D> &cpoints, vector<double> &cnts_radius_out, vector<vector<long int> > &cstructures, vector<GNP> &gnps)const;
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    //CNTs
    //Generate a network defined by points and connections
    int Generate_cnt_network_threads_mt(const Simu_para &simu_para, const Geom_sample &geom_sample, const Nanotube_Geo &nanotube_geo, const Cutoff_dist &cutoffs, vector<vector<Point_3D> > &cnts_points, vector<double> &cnts_radius)const;
    int CNT_seeds(const vector<unsigned int> &CNT_seeds, unsigned int net_seeds[])const;
    //This functions initializes the vector sub-regions
    int Initialize_cnt_subregions(const Geom_sample &sample_geom, int n_subregion[], vector<vector<long int> > &sectioned_domain)const;
    int Get_length_and_radius(const Nanotube_Geo &nanotube_geo, mt19937 &engine, uniform_real_distribution<double> &dist, double &cnt_length, double &cnt_rad)const;
    //Check if the current CNT is penetrating another CNT, i.e. is the new point is overlapping other point
    int Check_penetration(const Geom_sample &geom_sample, const Nanotube_Geo &nanotube_geo, const vector<vector<Point_3D> > &cnts, const vector<vector<int> > &global_coordinates, const vector<vector<long int> > &sectioned_domain, const vector<double> &radii, const vector<Point_3D> &cnt_new, const int n_subregions[], const double &cnt_rad, const double &d_vdW, int &point_overlap_count, int &point_overlap_count_unique, Point_3D &point)const;
    int Get_cnt_point_subregion(const Geom_sample &geom_sample, const int n_subregions[], const Point_3D &point)const;
    void Get_penetrating_points(const vector<vector<Point_3D> > &cnts, const vector<vector<int> > &global_coordinates, const vector<long int> &subregion_vec, const vector<double> &radii, const double &rad_plus_dvdw, Point_3D &point, vector<Point_3D> &affected_points, vector<double> &cutoffs_p, vector<double> &distances)const;
    void Move_point(const vector<Point_3D> &cnt_new, Point_3D &point, vector<double> &cutoffs, vector<double> &distances, vector<Point_3D> &affected_points)const;
    void One_overlapping_point(const vector<double> &cutoffs, const vector<double> &distances, const vector<Point_3D> &affected_points, Point_3D &point)const;
    void Two_overlapping_points(const vector<double> &cutoffs, const vector<Point_3D> &affected_points, Point_3D &point)const;
    void Three_or_more_overlapping_points(const vector<double> &cutoffs, const vector<double> &distances, const vector<Point_3D> &affected_points, Point_3D &point)const;
    int Check_segment_orientation(const Point_3D &point, const vector<Point_3D> &cnt_new)const;
    //Calculate the effective portion (length) which falls into the given region defined by a cuboid
    double Length_inside_sample(const Geom_sample &geom_sample, const Point_3D &prev_point, const Point_3D &new_point, const bool &is_prev_inside_sample, bool &is_new_inside_sample)const;
    void Add_cnt_point_to_overlapping_regions(const Geom_sample &geom_sample, Point_3D point, long int global_num, const int n_subregions[], vector<vector<long int> > &sectioned_domain)const;
    int Calculate_t(int a, int b, int c, int sx, int sy)const;
    int Transform_points_cnts(const Geom_sample &geom_sample, const Nanotube_Geo &nano_geo, vector<vector<Point_3D> > &cnts_points, vector<Point_3D> &cpoints, vector<double> &radii_in, vector<double> &radii_out, vector<vector<long int> > &cstructures)const;
    int Add_cnts_inside_sample(const Geom_sample &geom_sample, const Nanotube_Geo &nano_geo, const int &CNT_old, vector<Point_3D> &cnt, vector<Point_3D> &cpoints, vector<double> &radii_in, vector<double> &radii_out, vector<vector<long int> > &cstructures, long int &point_count, int &cnt_count)const;
    int Add_cnt_segment(const Geom_sample &geom_sample, const bool &is_first_inside_sample, const int &start, const int &end, const int &min_points, const int &CNT_old, vector<Point_3D> &cnt, vector<Point_3D> &cpoints, vector<double> &radii_in, vector<double> &radii_out, vector<vector<long int> > &cstructures, long int &point_count, int &cnt_count)const;
    int Add_boundary_point(const Geom_sample &geom_sample, const Point_3D &p_outside, const Point_3D &p_inside, const int &cnt_count, vector<Point_3D> &cpoints, vector<long int> &struct_temp, long int &point_count)const;
    Point_3D Find_intersection_at_boundary(const Geom_sample &geom_sample, const Point_3D &p_outside, const Point_3D &p_inside)const;
    int Recalculate_vol_fraction_cnts(const Geom_sample &geom_sample, const Simu_para &simu_para, const Nanotube_Geo &nano_geo, const vector<Point_3D> &cpoints, const vector<double> &radii, const vector<vector<long int> > &cstructures)const;
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    //Common operations
    int Get_random_value_mt(const string &dist_type, mt19937 &engine, uniform_real_distribution<double> &dist, const double &min, const double &max, double &value)const;
    int Get_point_in_cuboid_mt(const cuboid &cub, Point_3D &point, mt19937 &engine_x, mt19937 &engine_y, mt19937 &engine_z, uniform_real_distribution<double> &dist)const;
    int Get_initial_direction_mt(const string &dir_distrib_type, const double &ini_theta, const double &ini_phi, mt19937 &engine_inital_direction, uniform_real_distribution<double> &dist_initial, MathMatrix &rotation)const;
    int Get_normal_direction_mt(const double &omega, double &cnt_theta, double &cnt_phi, mt19937 &engine_theta, mt19937 &engine_phi, uniform_real_distribution<double> &dist)const;
    //Transform angles into matrix
    MathMatrix Get_transformation_matrix(const double &theta, const double &phi)const;
    //Calculate the coordinates of the new CNT point (transformation of coordinates)
    Point_3D Get_new_point(MathMatrix &Matrix, const double &Rad)const;
    //This function checks if a point is inside a cuboid
    int Is_point_inside_cuboid(const cuboid &cub, const Point_3D &point)const;
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    //GNPs
    int Generate_gnp_network_mt(const Simu_para &simu_para, const GNP_Geo &gnp_geo, const Geom_sample &geom_sample, const Cutoff_dist &cutoffs, vector<GNP> &gnps, vector<vector<int> > &sectioned_domain, double &gnp_vol_tot, double &gnp_wt_tot)const;
    int Initialize_gnp_subregions(const Geom_sample &sample_geom, int n_subregion[], vector<vector<int> > &sectioned_domain)const;
    int GNP_seeds(const vector<unsigned int> &GNP_seeds, unsigned int net_seeds[])const;
    int Generate_gnp(const GNP_Geo &gnp_geo, GNP &gnp, mt19937 &engine_l, mt19937 &engine_t, uniform_real_distribution<double> &dist)const;
    int Update_gnp_plane_equations(GNP &gnp)const;
    int Deal_with_gnp_interpenetrations(const Geom_sample &geom_sample, const Cutoff_dist &cutoffs, const vector<GNP> &gnps, const int n_subregions[], GNP &gnp_new, vector<vector<int> > &sectioned_domain, bool &rejected)const;
    int Get_gnp_subregions(const Geom_sample &geom_sample, const GNP &gnp_new, const int n_subregions[], set<int> &subregions)const;
    int Add_gnp_subregions_to_set_for_gnp_point(const Geom_sample &geom_sample, const Point_3D &new_point, const int n_subregions[], set<int> &subregions)const;
    int Get_gnps_in_subregions(const vector<vector<int> > &sectioned_domain, const set<int> &subregions, set<int> &gnp_set)const;
    int Move_gnps_if_needed(const Cutoff_dist &cutoffs, const vector<GNP> &gnps, set<int> &gnp_set, GNP &gnp_new, bool &displaced)const;
    int Find_direction_of_touching_gnps(Collision_detection &GJK_EPA, const GNP &gnpA, const GNP &gnpB, Point_3D &N)const;
    int Move_gnp(const Point_3D &displacement, GNP &gnp)const;
    int Add_valid_gnp_to_subregions(const int &gnp_new_idx, const set<int> &subregions, vector<vector<int> > &sectioned_domain)const;
    int Calculate_generated_gnp_vol_and_update_total_vol(const GNP_Geo gnp_geom, const Geom_sample &sample_geom, const GNP &gnp, double &gnp_vol_tot, bool &is_all_outside)const;
    int Approximate_gnp_volume_inside_sample(const Geom_sample &sample_geom, const GNP &gnp, double &gnp_vol, bool &is_all_outside)const;
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    //Mixed Particles
    int Generate_cnt_network_threads_among_gnps_mt(const Simu_para &simu_para, const GNP_Geo &gnp_geo, const Nanotube_Geo &nanotube_geo, const Geom_sample &geom_sample, const Cutoff_dist &cutoffs, const vector<GNP> &gnps, const vector<vector<int> > &sectioned_domain_gnp, const double &gnp_vol_tot, const double &gnp_wt_tot, vector<vector<Point_3D> > &cnts_points, vector<double> &cnts_radius)const;
    int Generate_initial_point_of_cnt(const Geom_sample &geom_sample, const Simu_para &simu_para, const vector<vector<Point_3D> > &cnts, const vector<double> &radii, vector<Point_3D> &new_cnt, const double &rad_plus_dvdw, const vector<vector<int> > &global_coordinates, const vector<vector<long int> > &sectioned_domain_cnt, const vector<GNP> &gnps, const vector<vector<int> > &sectioned_domain_gnp, const int n_subregions[], Point_3D &new_point, mt19937 &engine_x, mt19937 &engine_y, mt19937 &engine_z, uniform_real_distribution<double> &dist)const;
    int Check_mixed_interpenetration(const Geom_sample &geom_sample, const vector<vector<Point_3D> > &cnts, const vector<double> &radii, vector<Point_3D> &new_cnt, const double &rad_plus_dvdw, const vector<vector<int> > &global_coordinates, const vector<vector<long int> > &sectioned_domain_cnt, const vector<GNP> &gnps, const vector<vector<int> > &sectioned_domain_gnp, const int n_subregions[], Point_3D &new_point)const;
    int Get_gnp_penetrating_points(const vector<GNP> &gnps, const vector<int> &subregion_gnp, const double &cutoff, const Point_3D &new_point, vector<Point_3D> &affected_points, vector<double> &cutoffs_p, vector<double> &distances, GNP &affected_gnp)const;
    Point_3D Get_gnp_point_closest_to_point(const GNP &gnp, const Point_3D &P, double &distance) const;
    Point_3D Get_point_closest_to_large_gnp_face(const GNP &gnp, const int &V0, const int &V1, const int &V2, const int &V3, const int &F, const Point_3D &P, double &distance)const;
    Point_3D Distance_from_point_to_edge(const Point_3D &P, const Point_3D &V1, const Point_3D &V2, double &distance)const;
    int Deal_with_point_inside_gnp(const Geom_sample &geom_sample, const GNP &gnp, const double &cutoff, vector<Point_3D> &new_cnt, Point_3D &new_point)const;
    Point_3D Find_closest_face_and_relocate(const GNP &gnp, const Point_3D &new_point, const double &cutoff)const;
    bool Is_point_inside_gnp(const GNP &gnp, const Point_3D &P)const;
    Point_3D Find_new_position_one_point_inside_gnp(const GNP &gnp, const Point_3D &P1, const Point_3D &P2, const double &cutoff)const;
    int Find_new_position_two_points_inside_gnp(const Geom_sample &geom_sample, const GNP &gnp, const double &cutoff, vector<Point_3D> &new_cnt, Point_3D &new_point)const;
    vector<int> Find_gnp_faces_intersecting_boundary(const Geom_sample &geom_sample, const GNP &gnp)const;
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    
};
//-------------------------------------------------------
#endif
//===========================================================================
