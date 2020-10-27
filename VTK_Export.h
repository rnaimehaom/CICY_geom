//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GS networks
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
    int Export_gnps(const vector<GNP> &gnps, const string &filename)const;
    int Add_all_gnp_vertices(const vector<GNP> &gnps, ofstream &otec)const;
    int Add_points_from_array(const Point_3D points[], const int &arr_size, ofstream &otec)const;
    int Add_ofsets_for_gnps(const vector<GNP> &gnps, ofstream &otec)const;
    int Add_connectivity_for_gnps(const vector<GNP> &gnps, ofstream &otec, const long int &cnt_points_offset = 0)const;
    int Export_hybrid_material(const vector<Point_3D> &points, const vector<vector<long int> > &structure, const vector<GNP> &gnps, const string &filename)const;
};

#endif /* VTK_Export_hpp */