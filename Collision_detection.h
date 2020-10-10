//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GS networks
//OBJECTIVE:    Detect interpenetrations (collisions) of GS (polyhedra)
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#ifndef Collision_detection_h
#define Collision_detection_h

#include "Geometry_3D.h"


//-------------------------------------------------------
class Collision_detection
{
public:
    
    //Constructor
    Collision_detection(){};
    
    //---------------------------------------------------------------------------
    //GJK
    int GJK(const GNP &gnp1, const GNP &gnp_new, vector<Point_3D> &simplex, bool &p_flag, bool &t_flag);
    Point_3D Support_AB(const Point_3D &D, const Point_3D verticesA[], const Point_3D verticesB[]);
    Point_3D Support_map(const Point_3D &D, const Point_3D vertices[]);
    bool Is_origin_in_simplex(vector<Point_3D> &simplex, Point_3D &A, Point_3D &D, bool &t_flag);
    bool Update_simplex_case2(vector<Point_3D> &simplex, Point_3D &A, Point_3D &D, bool &t_flag);
    bool Update_simplex_case3(vector<Point_3D> &simplex, Point_3D &A, Point_3D &D, bool &t_flag);
    Point_3D Normal_to_face_ouside_tetrahedron(const Point_3D &A, const Point_3D &B, const Point_3D &C, const Point_3D &not_in_face);
    bool Update_simplex_case2_tetrahedron(vector<Point_3D> &simplex, Point_3D &A, Point_3D &D, Point_3D &E, Point_3D &ABC, Point_3D &AEB, Point_3D &ACE, bool &t_flag);
    bool Check_actual_touch(const Point_3D verticesA[], const Point_3D verticesB[], const Point_3D &D, vector<Point_3D> &simplex);
    bool Is_point_in_simplex(const vector<Point_3D> &simplex, const Point_3D &A);
    //---------------------------------------------------------------------------
    //Extended Polytope Algorithm (EPA)
    int EPA(const Point_3D verticesA[], const Point_3D verticesB[], vector<Point_3D> &simplex, Point_3D &normal, double &PD);
    void Normal_and_distance_to_origin(const vector<Point_3D> &simplex, const TrFace &f, Point_3D &normal, double &distance);
    void Find_closest_face(const vector<TrFace> &Faces, const vector<Point_3D> &Normals, const vector<double> &distances, Point_3D &normal, double &PD);
    void Reconstruct_simplex(const vector<Point_3D> &simplex, const int &S, vector<TrFace> &Faces, vector<Point_3D> &Normals, vector<double> &distances);
    int Is_edge_in_edges(const Edge &e, const vector<Edge> &edges);
    //---------------------------------------------------------------------------
    //Distance estimation
    int Distance_and_direction_from_simplex_to_origin(const vector<Point_3D> &simplex, Point_3D &N, double &dist);
    
};


#endif /* Collision_detection_h */
