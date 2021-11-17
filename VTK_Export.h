//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Export visualization files in ASCII VTK legacy file format
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#ifndef VTK_Export_h
#define VTK_Export_h

#include<vector>
#include<stdlib.h>
#include<stdio.h>
#include<vector>
#include<string>
#include "Geometry_3D.h"

//-------------------------------------------------------
class VTK_Export
{
public:
    
    //Constructor
    VTK_Export(){};
    
    //---------------------------------------------------------------------------
    //CNTs
    int Export_cnts_1D_vector(const vector<Point_3D> &points, const vector<vector<long int> > &structure, const string &filename)const;
    int Add_header(ofstream &otec)const;
    int Add_points_from_vector(const vector<Point_3D> &points, ofstream &otec)const;
    int Add_offsets_for_structure(const vector<vector<long int> > &structure, ofstream &otec)const;
    int Add_connectivity_for_structure(const vector<vector<long int> > &structure, ofstream &otec)const;
    int Export_cnts_2D_vector(const vector<vector<Point_3D> > &points, const string &filename)const;
    long int Count_number_of_points(const vector<vector<Point_3D> > &points)const;
    int Add_points_from_2D_vector(const vector<vector<Point_3D> > &points, ofstream &otec)const;
    int Add_offsets_for_2D_vector(const vector<vector<Point_3D> > &points, ofstream &otec)const;
    int Add_connectivity_for_2D_vector(const vector<vector<Point_3D> > &structure, ofstream &otec)const;
    //---------------------------------------------------------------------------
    int Export_from_cnt_indices(const vector<Point_3D> &points, const vector<vector<long int> > &indices, const string &filename)const;
    int Count_points_and_lines_from_indices(const vector<vector<long int> > &indices, long int &n_points, int &n_lines)const;
    int Add_points_from_indices(const vector<Point_3D> &points, const vector<vector<long int> > &indices, ofstream &otec)const;
    int Add_offsets_from_indices(const vector<vector<long int> > &indices, ofstream &otec)const;
    int Add_connectivity_from_indices(const vector<vector<long int> > &indices, ofstream &otec)const;
    //---------------------------------------------------------------------------
    int Export_cnts_in_cluster(const vector<Point_3D> &points, const vector<vector<long int> > &structure, const vector<int> &cluster, const string &filename)const;
    int Count_points_in_cluster(const vector<int> &cluster, const vector<vector<long int> > &structure, long int &n_points)const;
    int Add_points_in_cluster(const vector<Point_3D> &points, const vector<vector<long int> > &structure, const vector<int> &cluster, ofstream &otec)const;
    int Add_offsets_from_cluster(const vector<vector<long int> > &structure, const vector<int> &cluster, ofstream &otec)const;
    int Add_connectivity_from_cluster(const vector<vector<long int> > &structure, const vector<int> &cluster, ofstream &otec)const;
    //---------------------------------------------------------------------------
    int Export_single_cnt(const vector<Point_3D> &points, const string &filename)const;
    //---------------------------------------------------------------------------
    int Export_from_cnt_structure(const vector<Point_3D> &points, const vector<vector<long int> > &structure, const string &filename)const;
    int Count_points_in_structure(const vector<vector<long int> > &structure, long int &n_points)const;
    int Add_points_from_structure(const vector<Point_3D> &points, const vector<vector<long int> > &structure, ofstream &otec)const;
    int Add_offsets_from_structure(const vector<vector<long int> > &structure, ofstream &otec)const;
    int Add_connectivity_from_structure(const vector<vector<long int> > &structure, ofstream &otec)const;
    //---------------------------------------------------------------------------
    //GNPs
    int Export_gnps(const vector<GNP> &gnps, const string &filename)const;
    int Add_all_gnp_vertices(const vector<GNP> &gnps, ofstream &otec)const;
    int Add_points_from_array(const Point_3D points[], const int &arr_size, ofstream &otec)const;
    int Add_ofsets_for_gnps(const int &n_gnps, ofstream &otec)const;
    int Add_connectivity_for_gnps(const long int &n_gnps, ofstream &otec, const long int &cnt_points_offset = 0)const;
    //---------------------------------------------------------------------------
    int Export_gnps_in_cluster(const vector<GNP> &gnps, const vector<int> &cluster, const string &filename);
    int Add_all_gnp_vertices_from_cluster(const vector<GNP> &gnps, const vector<int> &cluster, ofstream &otec)const;
    //---------------------------------------------------------------------------
    int Export_gnps_single_files(const vector<GNP> &gnps, const string &base_filename)const;
    int Export_single_gnp(const GNP& gnp_i, const string& filename)const;
    //---------------------------------------------------------------------------
    //Cuboid
    int Export_cuboid(const cuboid &cub, const string &filename);
    //---------------------------------------------------------------------------
    //Mixed
    int Export_hybrid_material(const vector<Point_3D> &points, const vector<vector<long int> > &structure, const vector<GNP> &gnps, const string &filename)const;
    //---------------------------------------------------------------------------
    //triangulations
    int Export_triangulation(const vector<Point_3D> &points, const vector<EdgeL> &triangulation, const string &filename)const;
    int Add_points_from_triangulation_edges(const vector<Point_3D> &points, const vector<EdgeL> &triangulation, ofstream &otec)const;
    int Add_offsets_for_trinagulation(const vector<EdgeL> &triangulation, ofstream &otec)const;
    int Add_connectivity_for_trinagulation(const vector<EdgeL> &triangulation, ofstream &otec)const;
    //---------------------------------------------------------------------------
    //Tetrahedron/Triangles
    int Export_triangles(const vector<Point_3D>& vertices, const vector<TrFaceL>& triangles, const string& filename)const;
    int Add_ofsets_for_triangles(const int& n_tri, ofstream& otec)const;
    int Add_connectivity_for_triangles(const vector<TrFaceL>& triangles, ofstream& otec)const;
    int Export_supertetrahedron(const vector<Point_3D>& vertices, const vector<TrFaceL>& triangles, const string& filename)const;
    //---------------------------------------------------------------------------
    //Points
    int Export_point_array(const Point_3D points[], const int& size, const string& filename)const;
    int Add_point_coordinates_from_array(const Point_3D points[], const int& size, ofstream& otec)const;
    int Add_vertices_offsets_connectivity_for_n_points(const int& n, ofstream& otec)const;
    int Add_consecutive_numbers(const int& n, ofstream& otec)const;
    int Export_point_vector(const vector<Point_3D>& points, const string& filename)const;
    int Add_point_coordinates_from_vector(const vector<Point_3D>& points, ofstream& otec)const;
    int Export_selected_points_in_array(const vector<int>& vertices, const Point_3D points[], const string& filename)const;
    int Add_point_coordinates_from_vertex_vector(const vector<int>& vertices, const Point_3D points[], ofstream& otec)const;
};

#endif /* VTK_Export_hpp */
