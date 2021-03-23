//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GS networks
//OBJECTIVE:    Functions to print (output into a file) different types of data structures
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#ifndef PRINTER_H
#define PRINTER_H

#include "Input_Reader.h"

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
    void Append_1d_vec(const vector<double> &list, const string &filename);
    void Print_1d_vec(const vector<long int> &list, const string &filename);
    void Print_2d_vec(const vector<vector<int> > &num_mat, const string &filename);
    void Print_2d_vec(const vector<vector<long int> > &num_mat, const string &filename);
    void Print_2d_vec(const vector<vector<double> > &num_mat, const string &filename);
    void Print_4_vertices_gnps(const vector<GNP> &gnps, const int &prec, const string &filename);
};
#endif
