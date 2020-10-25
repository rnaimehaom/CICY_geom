//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GS networks
//OBJECTIVE:    Functions to export Tecplot visialization files
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#ifndef  TECPLOTEXPORT_H
#define TECPLOTEXPORT_H

#include "Input_Reader.h"
#include "Generate_Network.h"

//-------------------------------------------------------
class Tecplot_Export
{
public:
    //Data Member
    
    //Constructor
    Tecplot_Export(){};
    
    //Member Functions
    //The geometric structure of CNT network (by threads in Tecplot) in a cuboid
    int Export_network_threads(const struct cuboid &cub, const vector<vector<Point_3D> > &cnts_points)const;
    int Export_network_threads(const struct cuboid &cub, const vector<vector<Point_3D> > &cnts_points, const string &filename)const;
    int Export_triangulation_network_3dlines(const GCH &hybrid, const vector<Point_3D> &cnts_points, const vector<Point_3D> &gnps_points, const vector<vector<long int> > &structure, string filename)const;
    int Export_network_threads(const struct cuboid &cub, const int &n_cluster, const vector<vector<int> > &gnp_clusters, const vector<vector<int> > &cnt_clusters, const vector<vector<long int> > &structure, const vector<Point_3D> &cnts_points, const vector<GCH> &hybrid_particles, string &filename, string &family)const;
    int Export_network_3dlines(const struct cuboid &cub, const int &n_cluster, const vector<vector<int> > &gnp_clusters, const vector<vector<int> > &cnt_clusters, const vector<vector<long int> > &structure, const vector<Point_3D> &cnts_points, const vector<GCH> &hybrid_particles, string &filename, string &family)const;
    int Export_nano_3dlines(ofstream &otec, const vector<Point_3D> &cnts_points, const vector<vector<long int> > &structure)const;
    //The geometric structure of CNT network (by quadrilaterial elements in Tecplot)
    int Export_cnt_network_meshes(const int &zone_flag, const struct cuboid &cub, const vector<vector<Point_3D> > &cnts_points, const vector<double> &cnts_radius)const;
    //The geometric structure of CNT network (by quadrilaterial elements in Tecplot) with a specific filename. This function uses a 1D point vector, a 2D structure vector that references the point vector and a vector of hybrid particles
    int Export_network_meshes(const struct cuboid &cub, const int &n_cluster, const vector<vector<int> > &gnp_clusters, const vector<vector<int> > &cnt_clusters, const vector<vector<long int> > &structure, const vector<Point_3D> &cnts_points, const vector<double> &cnts_radius, const vector<GCH> &hybrid_particles, string &filename, string &family)const;
    
    //Export a 3D cuboid
    int Export_cuboid(ofstream &otec, const struct cuboid &cub)const;
    //Export 3D nanotube threads
    int Export_nano_threads(ofstream &otec, const vector<vector<Point_3D> > &cnts_points)const;
    int Export_nano_threads(ofstream &otec, const vector<Point_3D> &cnts_points, const vector<vector<long int> > &structure)const;
    //Export nanotube network by tetrahedron elements (Multiple zones in tecplot: each nanotube by one zone)
    int Export_cnts_meshes_multizones(const struct cuboid &cub, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles)const;
    //Export nanotube network by tetrahedron elements (Single zones in tecplot: all nanotubes by one zone)
    int Export_cnts_meshes_singlezone(const struct cuboid &cub, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles)const;
    //Export hybrid particle network by tetrahedron elements (Single zones in tecplot: all nanotubes by one zone) with a specific filename
    int Export_cnt_meshes_singlezone(ofstream &otec, const struct cuboid &cub, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles, string &filename, string &family)const;
    //Export a 3D cuboid with a random orientation
    int Generate_cnts_nodes_elements(vector<vector<Node> > &nodes, vector<vector<Element> > &eles, const vector<vector<Point_3D> > &cnts_points, const vector<double> &cnts_radius)const;
    int Get_angles_vector_in_spherial_coordinates(const Point_3D &normal, double &sita, double &pha)const;
    int Get_points_circle_in_plane(const Point_3D &center, const double &trans_sita, const double &trans_pha, const double &radius, const int &num_sec, vector<Node> &nod_temp)const;
    int Get_projected_points_in_plane(const Point_3D &center, const Point_3D &normal, const Point_3D &line, const int &num_sec, vector<Node> &nod_temp)const;
    MathMatrix Get_transformation_matrix(const double &theta, const double &phi)const;
    int Generate_cnts_nodes_elements(vector<vector<Node> > &nodes, vector<vector<Element> > &eles, const vector<Point_3D> &cnts_points, const vector<double> &cnts_radius, const vector<vector<long int> > &structure)const;
    int Export_randomly_oriented_gnps(ofstream &otec, const vector<GCH> &hybrid_particles, const vector<int> &gnp_cluster, string &family)const;
    int Export_randomly_oriented_gnps(ofstream &otec, const vector<GNP> &gnps)const;
    int Export_singlegnp(const GNP &gnps, const string &filename);
};
//-------------------------------------------------------
#endif
//===========================================================================
