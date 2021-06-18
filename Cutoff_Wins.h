//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Cut out an observation window
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#ifndef CUTOFFWINS_H
#define CUTOFFWINS_H

#include "Input_Reader.h"
#include "Shells.h"
#include "Printer.h"

//-------------------------------------------------------
class Cutoff_Wins
{
public:
    
    //List of CNTs inside the observation window
    vector<int> cnts_inside;
    //List of GNPs inside the observation window
    vector<int> gnps_inside;
    //CNTs located at a window boundary
    //boundary_cnt[i] contains CNTs in contact with boundary i
    vector<vector<int> > boundary_cnt;
    //CNT points located at a window boundary
    //boundary_cnt_pts[i] contains CNT points at boundary i
    vector<vector<long int> > boundary_cnt_pts;
    //GNPs located at a window boundary
    //boundary_gnp[i] contains GNPs in contact with boundary i
    vector<vector<int> > boundary_gnp;
    //GNP points located at a window boundary
    //boundary_gnp_pts[i] contains GNP points at boundary i
    vector<vector<long int> > boundary_gnp_pts;
    

    
    //Constructor
    Cutoff_Wins(){};
    
    //Member Functions
    int Extract_observation_window(const int &window, const string &particle_type, const Geom_sample &sample_geo, const cuboid &window_geo, const Nanotube_Geo &cnts_geo, vector<GNP> &gnps, vector<vector<long int> > &structure_cnt, vector<double> &radii, vector<Point_3D> &points_cnt, vector<vector<int> > &shells_cnt, vector<Shell> &shells_gnp, vector<vector<long int> > &structure_gnp, vector<Point_3D> &points_gnp);
    int Trim_boundary_cnts(const int &window, const string &particle_type, const Geom_sample &sample_geo, const cuboid &window_geo, const Nanotube_Geo &cnts, vector<Point_3D> &points_cnt, vector<vector<long int> > &structure_cnt, vector<vector<int> > &shells_cnt, vector<double> &radii);
    int Get_percolation_layer_cuboid(const double &cnt_rad, const cuboid &window_geo, cuboid &layer_geom);
    int Add_cnt_segment_to_structure(const string &particle_type, const cuboid &window_geo, const cuboid &layer_geom, const double var_shells[][3], const int &start, const int &end, const int &min_points, const int &CNT, const string &last_point_loc, vector<Point_3D> &points_cnt, vector<vector<long int> > &structure_cnt, vector<vector<int> > &shells_cnt, vector<double> &radii, int &segments, int &first_idx, int &last_idx);
    int Substitute_boundary_point(const cuboid &window_geo, const Point_3D &p_inside, Point_3D &p_outside);
    string Where_is(const Point_3D &point, const cuboid &window_geo);
    string Where_is_with_layer(const Point_3D &point, const cuboid &window_geo, const cuboid &layer_geom, const int &flag = 0, const double &cnt_rad = 0, const Point_3D &V = Point_3D(0.0,0.0,0.0));
    int Check_criteria_for_cnt_at_boundary(const double &p_coord, const double &min_boundary, const double &max_boundary, const double &min_layer, const double &max_layer, const double &V_cood, const double &cnt_rad);
    string Where_is_with_boundary(const Point_3D &point, const cuboid &window_geo, int &boundary);
    int Add_cnt_point_to_boundary_vectors(const cuboid &layer_geom, const Point_3D &P, const long int &P_num, const int &CNT_num);
    int In_which_boundary(const Point_3D &P, const cuboid &layer_geom);
    int Fill_cnts_inside(const vector<vector<long int> > &structure);
    int Fill_gnps_inside(const int &window, const cuboid &window_geo, const vector<Shell> &shells_gnp, vector<GNP> &gnps, vector<vector<long int> > &structure_gnp, vector<Point_3D> &points_gnp);
    int Find_gnp_boundary_points(const cuboid &window_geo, GNP &gnp, vector<vector<long int> > &structure_gnp, vector<Point_3D> &points_gnp);
    int Accumulate_boundary_points_due_to_intersections(const cuboid &window_geo, const GNP &gnp, const vector<string> &locations, vector<vector<Point_3D> > &points_acc);
    int Find_intersections_of_gnp_edges_with_window_boundaries(const cuboid &window_geo, const GNP &gnp, const vector<string> &locations, vector<vector<Point_3D> > &points_acc);
    int Find_two_intersections_of_gnp_edges_with_window(const cuboid &window_geo, const GNP &gnp, const Point_3D &V1, const Point_3D &V2, vector<Point_3D> &Pts);
    int Find_intersections_of_window_edges_with_gnp_faces(const cuboid &window_geo, const GNP &gnp, vector<vector<Point_3D> > &points_acc);
    int Does_edge_intersect_gnp(const Edge &edg, const Point_3D window_vertices[], const GNP &gnp, Point_3D &P);
    int Find_window_vertex_inside_gnp(const cuboid &window_geo, const GNP &gnp, vector<vector<Point_3D> > &points_acc);
    
private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================
