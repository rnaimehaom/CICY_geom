//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Dlaunay triangulation for a given set of points
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#ifndef TRIANGULATION_H
#define TRIANGULATION_H

#include "Input_Reader.h"
#include "Generate_Network.h"

//-------------------------------------------------------
class Triangulation
{
public:
    //Data Member
    
    //Constructor
    Triangulation(){};
    
    //3D triangulation
    int Generate_3d_trangulation(const vector<Point_3D> &points_gnp, const vector<long int> &structure_i, GNP &gnp_i);
    int Generate_trivial_triangulation(const vector<long int> &structure_i, GNP &gnp_i);
    int Bowyer_watson(const vector<Point_3D> &points_gnp, const vector<long int> &structure_i, GNP &gnp_i);
    int Generate_supertetrahedron(const GNP &gnp_i, vector<Point_3D> &vertices, vector<TrFaceL> &triangles);
    int Find_bad_triangles(const Point_3D &vertex, const vector<Point_3D> &points_gnp, vector<Point_3D> &vertices, vector<TrFaceL> &triangles, vector<EdgeL> &bad_triangles_edges);
    bool Is_in_circumcircle(const Point_3D &vertex_i, const vector<Point_3D> &points_gnp, vector<Point_3D> &vertices, const TrFaceL &triangle);
    Point_3D Calculate_circumcenter(const vector<Point_3D> &points_gnp, const vector<Point_3D> &vertices, const TrFaceL &triangle);
    Point_3D Get_point(const long int &index, const vector<Point_3D> &points_gnp, const vector<Point_3D> &vertices);
    int Add_triangle_as_edges(const TrFaceL &triangle, vector<EdgeL> &edges);
    int Add_new_triangles(const long int &vertex, vector<TrFaceL> &triangles, vector<EdgeL> &bad_triangles_edges);
    int Final_triangulation_edges(const vector<TrFaceL> &triangles, GNP &gnp_i);
    
    
};
//-------------------------------------------------------
#endif
//===========================================================================
