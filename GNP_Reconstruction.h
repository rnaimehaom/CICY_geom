//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Given a GNP that is partially inside the sample and that has undergone 
//				some deformation, reconstruct a parallelepiped with a squared base that
//				fits that given GNP
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#ifndef GNP_RECONSTRUCTION_H
#define GNP_RECONSTRUCTION_H

#include "Geometry_3D.h"
#include "Generate_Network.h"

class GNP_Reconstruction
{
public:

    //Constructor
    GNP_Reconstruction() {};

    //Member functions
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    //Full GNP
    int Reconstruct_full_gnp(GNP& gnp_i)const;
    int Get_plane_from_top_squared_face(GNP& gnp_i, Plane_3D& Pl_top)const;
    int Get_plane_from_bottom_squared_face(const Plane_3D& Pl_top, GNP& gnp_i, Plane_3D& Pl_bot, double& d_planes)const;
    int Fit_squared_faces_on_parallel_planes(const Plane_3D& Pl_top, const Plane_3D& Pl_bot, double& d_planes, GNP& gnp_i)const;
    int Adjust_reference_edge(const GNP& gnp_i, const Point_3D& N_top, const double& d_planes, Point_3D& R1, Point_3D& R2, Point_3D& R1R2_hat, Point_3D& N_r1r2)const;
    int Get_gnp_side_length(const GNP& gnp_i, const Point_3D& M, const Point_3D& R1R2_hat, const Point_3D& N_r1r2, double& l_gnp)const;
    int Recalculate_gnp_vertices(const Point_3D& M, const Point_3D& R1R2_hat, const Point_3D& N_r1r2, const Point_3D& N_bot, const double& l_gnp, const double& d_planes, GNP& gnp_i)const;
    int Update_gnp_center(GNP& gnp_i)const;
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    //Partial GNP
    int Reconstruct_partial_gnp(const vector<bool>& vertex_flags, GNP& gnp_i)const;
    int Find_reconstruction_case(const vector<bool>& vertex_flags, int& gnp_case)const;
    int Three_short_edges(const vector<bool>& vertex_flags, GNP& gnp_i, Point_3D& N_top)const;
    int Get_reference_edge_case3(const vector<bool>& vertex_flags, int& R1, int& R2, int& LT, int& LB)const;
    int Adjust_vertices_along_thickness_planes_case3(const vector<bool>& vertex_flags, const int& R1, const int& R2, GNP& gnp_i, Point_3D& N)const;
    int Adjust_thin_faces_case3(const vector<bool>& vertex_flags, const int& R1, const int& R2, const int& LT, const int& LB, GNP& gnp_i)const;
    int Get_thin_face_case3(const int& R1, const int& R2, const int& LT, const int& LB, const int& RT, GNP& gnp_i, Plane_3D& Pl_thin)const;
    int Get_parallel_thin_face_case3(const vector<bool>& vertex_flags, const int other_v[], const int& R1, const int& R2, const int& LT, const int& LB, const int& RT, const Plane_3D& Pl_thin, GNP& gnp_i, double& l1)const;
    int Calculate_gnp_vertices_case3(const vector<bool>& vertex_flags, const int& R1, const int& R2, const int& LT, const int& LB, const int& RT, const Point_3D& N_top, GNP& gnp_i)const;
    int Two_consecutive_short_edges(const vector<bool>& vertex_flags, GNP& gnp_i, Point_3D& N_top)const;
    int Get_reference_vertices_case2(const vector<bool>& vertex_flags, int& R1, int& R2, int& R3, int& R4, int& O1, int& O2)const;
    int Get_reference_plane_case2(const int& R1, const int& R2, const int& R3, const int& R4, GNP& gnp_i)const;
    int Set_long_reference_edges_as_parallel_case2(const int& R1, const int& R2, const int& R3, const int& R4, const int& O1, const int& O2, GNP& gnp_i)const;
    int Find_gnp_length_and_calculate_vertices_case2(const int& R1, const int& R2, const int& R3, const int& R4, const int& O1, const int& O2, GNP& gnp_i, Point_3D& N_top)const;
};
//-------------------------------------------------------
#endif
//===========================================================================

