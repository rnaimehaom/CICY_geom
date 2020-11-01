//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GS networks
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
    //GNPs located at a window boundary
    //boundary_gnp[i] contains GNPs in contact with boundary i
    vector<vector<int> > boundary_gnp;
    //GNP points located at a window boundary
    //boundary_gnp_pts[i] contains GNP pointss at boundary i
    vector<vector<int> > boundary_gnp_pts;
    //GNP points that are located at the boundary
    vector<Point_3D> gnp_boundary_pts;
    
    //Deprecated:
    //This vector will help find points on the boundary. This is used in the direct electrifying algorithm
    vector<vector<short int> > boundary_flags_cnt;
    //This vector will help find points on the boundary. This is used in the direct electrifying algorithm
    vector<vector<short int> > boundary_flags_gnp;
    

    
    //Constructor
    Cutoff_Wins(){};
    
    //Member Functions
    int Extract_observation_window(const int &window, const string &particle_type, const Geom_sample &sample_geo, const cuboid &window_geo, const Nanotube_Geo &cnts_geo, vector<GNP> &gnps, vector<vector<long int> > &structure, vector<double> &radii, vector<Point_3D> &points_in, vector<vector<int> > &shells_cnt, vector<Shell> &shells_gnp);
    int Trim_boundary_cnts(const int &window, const Geom_sample &sample, const cuboid &window_geo, const Nanotube_Geo &cnts, vector<Point_3D> &points_in, vector<vector<long int> > &structure, vector<vector<int> > &shells_cnt, vector<double> &radii);
    int Add_cnt_segment_to_structure(const Geom_sample &sample_geo, const cuboid &window_geo, const double var_shells[][3], const int &start, const int &end, const int &min_points, const int &CNT, const string &last_point_loc, vector<Point_3D> &points_in, vector<vector<long int> > &structure, vector<vector<int> > &shells_cnt, vector<double> &radii, int &segments, int &first_idx, int &last_idx);
    int Substitute_boundary_point(const cuboid &window_geo, const Point_3D &p_inside, Point_3D &p_outside);
    string Where_is(const Point_3D &point, const cuboid &window_geo);
    string Where_is_with_boundary(const Point_3D &point, const cuboid &window_geo, int &boundary);
    void Add_to_boundary_vectors(const cuboid &window_geo, const Point_3D &point3d, const long int &point, const int &new_CNT);
    void Add_CNT_to_boundary(vector<int> &boundary, const int &CNT, const long int &point, const short int &flag1, const short int &flag2);
    int Fill_cnts_inside(const vector<vector<long int> > &structure);
    int Fill_gnps_inside(const int &window, const cuboid &window_geo, const vector<GNP> &gnps, const vector<Shell> &shells_gnp);
    int Find_gnp_boundary_points(const cuboid &window_geo, const GNP &gnp);
    int Accumulate_boundary_points_due_to_intersections(const cuboid &window_geo, const GNP &gnp, const vector<string> &locations, vector<vector<Point_3D> > &points_acc);
    int Find_intersections_of_gnp_edges_with_window_boundaries(const cuboid &window_geo, const GNP &gnp, const vector<string> &locations, vector<vector<Point_3D> > &points_acc);
    int Find_two_intersections_of_gnp_edges_with_window(const cuboid &window_geo, const GNP &gnp, const Point_3D &V1, const Point_3D &V2, vector<Point_3D> &Pts);
    int Find_intersections_of_window_edges_with_gnp_faces(const cuboid &window_geo, const GNP &gnp, vector<vector<Point_3D> > &points_acc);
    int Does_edge_intersect_gnp(const Point_3D &V1, const Point_3D &V2, const GNP &gnp, Point_3D &P);
    int Find_window_vertex_inside_gnp(const cuboid &window_geo, const GNP &gnp, vector<vector<Point_3D> > &points_acc);
    
    
    
    int Find_boundary_points_partially_inside_case(const cuboid &window_geo, const GNP &gnp, const vector<string> &locations, const vector<int> &inside_v);
    int Find_boundary_points_in_multiple_boundaries(const GNP &gnp, const vector<vector<int> > &boundary_acc);
    Point_3D Find_corner_inside_gnp(const cuboid &window_geo, const GNP &gnp);
    
    
    //Deprecated:
    int Add_GNPs_to_boundary(const vector<long int> &gnp_discrete, const vector<vector<int> > &boundaries, vector<Point_3D> &points_gnp);
    int Find_average_boundary_point(const vector<long int> &gnp_discrete, const vector<int> &boundary, vector<Point_3D> &points_gnp);
    int Change_repeated_seed(const int &CNT_original, const int &CNT_previous, int &index2_previous, int &index1_current, vector<vector<long int> > &structure, vector<Point_3D> &points_in);
    int Save_seeds(const vector<GCH> &hybrid_particles, const vector<vector<long int> > &structure, vector<long int> &seeds);
    int Compare_seeds(vector<GCH> &hybrid_particles, const vector<vector<long int> > &structure, const vector<long int> &seeds);
    
private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================
