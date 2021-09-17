//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GS networks
//OBJECTIVE:    Functions to print (output into a file) different types of data structures
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#ifndef PRINTER_H
#define PRINTER_H

#include "Input_Reader.h"
#include "Geometry_3D.h"

//---------------------------------------------------------------------------
class Printer
{
public:
    //Data Member
    
    //Constructor
    Printer(){};
    void Print_1d_vec(const vector<Point_3D> &list, const string &filename);
    void Print_1d_vec(const vector<char> &list, const string &filename);
    void Print_1d_vec(const vector<int> &list, const string &filename);
    void Print_1d_vec(const vector<double> &list, const string &filename);
    void Append(const double &value, const string &filename);
    void Append_1d_vec(const vector<double> &list, const string &filename);
    void Print_1d_vec(const vector<long int> &list, const string &filename);
    void Print_2d_vec(const vector<vector<int> > &num_mat, const string &filename);
    void Print_2d_vec(const vector<vector<long int> > &num_mat, const string &filename);
    void Print_2d_vec(const vector<vector<double> > &num_mat, const string &filename);
    void Print_4_vertices_gnps(const vector<GNP> &gnps, const int &prec, const string &filename);
    void Print_gnp_data(const vector<GNP> &gnps, const int &prec, const string &filename);
    double Recover_angle(const double &cos_alpha, const double &sin_alpha, const double &two_PI);
    void Print_cnt_points_and_structure(const cuboid &geom_sample, const vector<vector<long int> > &structure, const vector<Point_3D> &points_cnt, const vector<double> &radii, const int &prec, const string &filename_points, const string &filename_struct);
    void Check_if_close_enough_to_boundary(const cuboid &geom_sample, const Point_3D &P, const int &prec, int &cnt_points, ofstream &otec_points);
    void Print_cnt_points_and_structure_binary(const cuboid& geom_sample, const vector<vector<long int> >& structure, const vector<Point_3D>& points_cnt, const vector<double>& radii, const string& filename_points, const string& filename_struct);
};
#endif
