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
    int Get_reference_edge(const vector<bool>& vertex_flags, int& R1, int& R2, int& LT, int& LB)const;
    int Adjust_vertices_along_thickness_planes(const vector<bool>& vertex_flags, const int& R1, const int& R2, GNP& gnp_i, Point_3D& N)const;
    int Adjust_thin_faces(const vector<bool>& vertex_flags, const int& R1, const int& R2, const int& LT, const int& LB, GNP& gnp_i)const;

};
//-------------------------------------------------------
#endif
//===========================================================================

