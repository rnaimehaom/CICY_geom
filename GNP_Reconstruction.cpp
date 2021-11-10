//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Given a GNP that is partially inside the sample and that has undergone 
//				some deformation, reconstruct a parallelepiped with a squared base that
//				fits that given GNP
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "GNP_reconstruction.h"

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//This function reconstrucst a GNP after it has undergone a deformation
//Here, the paralellepiped that contains the eight vertices of the GNP on 
//its deformed position is found
int GNP_Reconstruction::Reconstruct_full_gnp(GNP& gnp_i)const
{
    //Plane that contains vertices v0 to v3
    Plane_3D Pl_top;

    //Get the equation of the plane that contains vertices v0 to v3
    if (!Get_plane_from_top_squared_face(gnp_i, Pl_top))
    {
        hout << "Error in Reconstruct_full_gnp when calling Get_plane_from_top_squared_face" << endl;
        return 0;
    }

    //Plane that contains vertices v4 to v7, initialized to be equal to top plane
    Plane_3D Pl_bot = Pl_top;

    //Variable to store the distance between planes
    double d_planes = 0.0;

    //Use the top plane (and the GNP vertices) to find the bottom plane
    if (!Get_plane_from_bottom_squared_face(Pl_top, gnp_i, Pl_bot, d_planes))
    {
        hout << "Error in Reconstruct_full_gnp when calling Get_plane_from_bottom_squared_face" << endl;
        return 0;
    }

    //Fit squared faces on the parallel planes
    if (!Fit_squared_faces_on_parallel_planes(Pl_top, Pl_bot, d_planes, gnp_i))
    {
        hout << "Error in Reconstruct_full_gnp when calling Fit_squared_faces_on_parallel_planes" << endl;
        return 0;
    }

    //Now that all vertices are calculated, update GNP center
    if (!Update_gnp_center(gnp_i))
    {
        hout << "Error in Reconstruct_full_gnp when calling Update_gnp_center" << endl;
        return 0;
    }

    //Reset rotation matrix of GNP
    gnp_i.rotation = MathMatrix(3, 3);

    //Object needed to recalculate GNP rotation matrix and plane equations
    Generate_Network GN;

    //Recalculate rotation matrix using the normal from one of the parallel planes
    if (!GN.Rotation_matrix_from_direction(Pl_top.N.x, Pl_top.N.y, Pl_top.N.z, gnp_i.rotation))
    {
        hout << "Error in Fit_squared_faces_on_parallel_planes when calling GN.Rotation_matrix_from_direction" << endl;
        return 0;
    }

    //Recalculate plane equations of GNP
    if (!GN.Update_gnp_plane_equations(gnp_i))
    {
        hout << "Error in Fit_squared_faces_on_parallel_planes when calling GN.Update_gnp_plane_equations" << endl;
        return 0;
    }

    return 1;
}
//This function "fits" the top vertices of a GNP (v0 to v3) onto a plane
int GNP_Reconstruction::Get_plane_from_top_squared_face(GNP& gnp_i, Plane_3D& Pl_top)const
{
    //Get plane equation for plane v0v1v3
    Pl_top = Plane_3D(gnp_i.vertices[0], gnp_i.vertices[1], gnp_i.vertices[3]);

    //Make sure normal vector goes outside the GNP
    if (Pl_top.N.dot(gnp_i.vertices[4] - gnp_i.vertices[0]) > 0.0) {
        Pl_top.N = Pl_top.N * (-1.0);
    }

    //Get the distance from v2 to the plane v0v1v3
    double d_v2 = Pl_top.distance_to(gnp_i.vertices[2]);

    //Check if v2 is on the plane v0v1v3
    if (d_v2 > Zero)
    {
        //v2 is not on the plane v0v1v3, check on which side of the plane v2 is located
        if ((gnp_i.vertices[2] - gnp_i.vertices[0]).dot(Pl_top.N) > 0.0)
        {
            //v2 is above plane, so set plane as v0v1v2
            Pl_top = Plane_3D(gnp_i.vertices[0], gnp_i.vertices[1], gnp_i.vertices[2]);

            //Make sure normal vector goes outside the GNP
            if (Pl_top.N.dot(gnp_i.vertices[4] - gnp_i.vertices[0]) > 0.0)
            {
                Pl_top.N = Pl_top.N * (-1.0);
            }

            //Get the distance from plane to v3
            double d_v3 = Pl_top.distance_to(gnp_i.vertices[3]);

            //Project v3 onto plane using the calculated distance
            //v3 is "below" the plane (inside the GNP), so move v3 towards the plane
            gnp_i.vertices[3] = gnp_i.vertices[3] + Pl_top.N * d_v3;
        }
        else
        {
            //v2 is below plane, so project v2 onto plane
            gnp_i.vertices[2] = gnp_i.vertices[2] + Pl_top.N * d_v2;
        }
    }

    return 1;
}
//From vertices v4 to v7, this function finds the vertex furthest from the top plane (Pl_top)
//and places the bottom plane (Pl_bot) parallel to the top plane at the furthest point found
int GNP_Reconstruction::Get_plane_from_bottom_squared_face(const Plane_3D& Pl_top, GNP& gnp_i, Plane_3D& Pl_bot, double& d_planes)const
{
    //Vector to store distances to top plane
    vector<double> dists(4, 0.0);

    //Variables to store the largest distance and index for the vertex with the largest distance
    //Initialize with vertex 4
    double d_max = Pl_top.distance_to(gnp_i.vertices[4]);
    int idx = 4;

    //Save distance into vector of distances
    dists[0] = d_max;

    //Iterate over vertices v5 to v7 and calculate distances from plane to those vertices
    //Detemine the largest distance and the index of the vertex
    for (int j = 5; j <= 7; j++)
    {
        //Calculate distance to vertex j
        double new_dist = Pl_top.distance_to(gnp_i.vertices[j]);

        //Check if vertex j is further away
        if (new_dist - d_max > Zero)
        {
            //A new maximum distance has been found so update variables
            d_max = new_dist;
            idx = j;

            //Save distance into vector of distances
            dists[j - 4] = d_max;
        }
    }

    //The maximum distance found is the distance between planes
    d_planes = d_max;

    //The bottom plane is initialized to be the same as the top plane, so the first three
    //coefficents are the same for both planes, the fourth coefficient is different
    //Update the fourth coefficient of the bottom plane so that it contains vertex idx
    Pl_bot.coef[3] = -Pl_top.N.dot(gnp_i.vertices[idx]);

    //Normal vector of bottom plane must go in the opposite direction of the top plane
    Pl_bot.N = Pl_top.N * (-1.0);

    //Iterate over all vertices again to check if they need to be projected or not
    for (int j = 4; j <= 7; j++)
    {
        //Ignore vertex idx that was used to calculate the parallel plane
        if (j != idx)
        {
            //Check if vetex j is too far from the plane
            if (dists[j] > Zero)
            {
                //Move vertex j towards the bottom plane
                gnp_i.vertices[j] = gnp_i.vertices[j] + Pl_bot.N * dists[j];
            }
        }
    }

    return 1;
}
//This function finds the squared faces of a GNP by "fiting" the vertices onto squares
int GNP_Reconstruction::Fit_squared_faces_on_parallel_planes(const Plane_3D& Pl_top, const Plane_3D& Pl_bot, double& d_planes, GNP& gnp_i)const
{
    //Initialize vertex R1 and R2 with v0 and v1, respectively
    Point_3D R1 = gnp_i.vertices[0], R2 = gnp_i.vertices[1];

    //Get unit vector along edge R1R2
    Point_3D R1R2_hat = (R2 - R1).unit();

    //Get vector perpendicular to edge R1R2
    Point_3D N_r1r2 = R1R2_hat.cross(Pl_top.N);

    //Check if the reference edge needs to be adjusted
    if (!Adjust_reference_edge(gnp_i, Pl_top.N, d_planes, R1, R2, R1R2_hat, N_r1r2))
    {
        hout << "Error in Fit_squared_faces_on_parallel_planes when calling Adjust_reference_edge" << endl;
        return 0;
    }

    //Get midpoint of edge R1R2
    Point_3D M = (R1 + R2) * 0.5;

    //Get the GNP length
    double l_gnp = 0.0;
    if (!Get_gnp_side_length(gnp_i, M, R1R2_hat, N_r1r2, l_gnp))
    {
        hout << "Error in Fit_squared_faces_on_parallel_planes when calling Get_gnp_side_length" << endl;
        return 0;
    }

    //Re-calculate the 8 GNP vertices
    if (!Recalculate_gnp_vertices(M, R1R2_hat, N_r1r2, Pl_bot.N, l_gnp, d_planes, gnp_i))
    {
        hout << "Error in Fit_squared_faces_on_parallel_planes when calling Recalculate_gnp_vertices" << endl;
        return 0;
    }

    //Update the GNP length and thckness in the GNP object
    gnp_i.l = l_gnp;
    gnp_i.t = d_planes;

    return 1;
}
//This function checks if the reference edge to determine the squared faces of a GNP
//needs to be adjusted
int GNP_Reconstruction::Adjust_reference_edge(const GNP& gnp_i, const Point_3D& N_top, const double& d_planes, Point_3D& R1, Point_3D& R2, Point_3D& R1R2_hat, Point_3D& N_r1r2)const
{

    //Check if v4 is inside the GNP, taking edge R1R2 as part of the GNP
    if (N_r1r2.dot(gnp_i.vertices[4] - R2) > Zero)
    {
        //Set the projection of v4 onto the top plane as the reference R1
        R1 = gnp_i.vertices[4] + N_top * d_planes;

        //Re-calculate unit vector along edge R1R2
        R1R2_hat = (R2 - R1).unit();

        //Re-calculate vector perpendicular to edge R1R2
        N_r1r2 = R1R2_hat.cross(N_top);
    }

    //Check if v5 is inside the GNP, taking edge R1R2 as part of the GNP
    if (N_r1r2.dot(gnp_i.vertices[5] - R1) > Zero)
    {
        //Set the projection of v5 onto the top plane as the reference R1
        R1 = gnp_i.vertices[5] + N_top * d_planes;

        //Re-calculate unit vector along edge R1R2
        R1R2_hat = (R2 - R1).unit();

        //Re-calculate vector perpendicular to edge R1R2
        N_r1r2 = R1R2_hat.cross(N_top);
    }

    return 1;
}
//This function calculated the GNP side lengt
//This function is used when the GNP is completely inside the sample
int GNP_Reconstruction::Get_gnp_side_length(const GNP& gnp_i, const Point_3D& M, const Point_3D& R1R2_hat, const Point_3D& N_r1r2, double& l_gnp)const
{
    //Get the distance from M to v0 along vector R1R2_hat
    double dx_half = abs(R1R2_hat.dot(gnp_i.vertices[0] - M));

    //Iterate over the remaining vertices and find the maximum distance from M
    //along vector R1R2_hat
    for (int i = 1; i < 8; i++)
    {
        //Calculate distance to vertex i along vector R1R2_hat
        double d_new = abs(R1R2_hat.dot(gnp_i.vertices[i] - M));

        //Check if new distance is larger than the current maximum
        if (d_new - dx_half > Zero)
        {
            //Update maximum distance
            dx_half = d_new;
        }
    }

    //Get the distance from R1R2 to v2
    double dy = abs(N_r1r2.dot(gnp_i.vertices[2] - M));

    //Remaining vertices on the "opposite side" of R1R2
    int v_set[] = { 3, 6, 7 };

    //Iterate over the remaining vertices on the "opposite side" of R1R2 and
    //find the maximum distance from R1R2
    for (int i = 0; i < 3; i++)
    {
        //Get vertex number
        int v_idx = v_set[i];

        //Calculate distance to vertex v_idx
        double d_new = abs(N_r1r2.dot(gnp_i.vertices[v_idx] - M));

        //Check if new distance is larger than the current maximum
        if (d_new - dy > Zero)
        {
            //Update maximum distance
            dy = d_new;
        }
    }

    //The new GNP side-length is the maximum distance between the two found
    l_gnp = max(2.0 * dx_half, dy);

    return 1;
}
//This function recalculate the GNP vertices for the case when the GNP is completely inside 
//the sample. Here, the midpoint of the reference edge, the distance between planes and
//the distance found as teh new l_gnp are used to calculate the new GNP vertices
int GNP_Reconstruction::Recalculate_gnp_vertices(const Point_3D& M, const Point_3D& R1R2_hat, const Point_3D& N_r1r2, const Point_3D& N_bot, const double& l_gnp, const double& d_planes, GNP& gnp_i)const
{
    //Half GNP side length is also used
    double l_gnp_half = l_gnp * 0.5;

    //This is used more than once
    Point_3D R1R2_half = R1R2_hat * l_gnp_half;

    //Calculate v0 from M
    gnp_i.vertices[0] = M - R1R2_half;

    //Calculate v1 from M
    gnp_i.vertices[1] = M + R1R2_half;

    //This is used more than once
    Point_3D N_length = N_r1r2 * l_gnp;

    //Calculate v2 from v1
    gnp_i.vertices[2] = gnp_i.vertices[1] - N_length;

    //Calculate v3 from v0
    gnp_i.vertices[3] = gnp_i.vertices[0] - N_length;

    //This is used more than once
    Point_3D N_thick = N_bot * d_planes;

    //Calculate v4 from v0
    gnp_i.vertices[4] = gnp_i.vertices[0] + N_thick;

    //Calculate v5 from v1
    gnp_i.vertices[5] = gnp_i.vertices[1] + N_thick;

    //Calculate v6 from v2
    gnp_i.vertices[6] = gnp_i.vertices[2] + N_thick;

    //Calculate v7 from v3
    gnp_i.vertices[7] = gnp_i.vertices[3] + N_thick;

    return 1;
}
//This function updates the GNP center, provided all vertices were already updated
int GNP_Reconstruction::Update_gnp_center(GNP& gnp_i)const
{
    //Initialize a point at the origin
    Point_3D C(0.0, 0.0, 0.0);

    //Accumulate all vertices in point C
    for (size_t i = 0; i < 8; i++)
    {
        C = C + gnp_i.vertices[i];
    }

    //The average of C is the centroid of the GNP center
    gnp_i.center = C / 8.0;

    return 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//This function reconstructs a partial GNP
//This is done in a case by case basis to facilitate reconstruction of the GNP
int GNP_Reconstruction::Reconstruct_partial_gnp(const vector<bool>& vertex_flags, GNP& gnp_i)const
{
    //Variable to store the case for GNP reconstruction
    int gnp_case = -1;

    //Find the case
    if (!Find_reconstruction_case(vertex_flags, gnp_case))
    {
        hout << "Error in Reconstruct_partial_gnp when calling Find_reconstruction_case" << endl;
        return 0;
    }

    //Reconstruct depending on the case
    switch (gnp_case)
    {
    case 3:
        //3 consecutive short edges
        break;
    case 2:
        //2 consecutive short edges
        break;
    case 1:
        //1 short edge
        break;
    case 0:
        //0 short edges
        break;
    case 4:
        //2 non-consecutive short edges
        break;
    default:
        break;
    }

    return 1;
}
//This function finds the case of the GNP reconstruction
int GNP_Reconstruction::Find_reconstruction_case(const vector<bool>& vertex_flags, int& gnp_case)const
{
    //Variable to count the number of edges present
    int n_edges = 0;

    //Variables for the case of two non-consecutive edges
    bool e1 = false, e2 = false;

    //Check which short edges are present and count them
    if (vertex_flags[0] && vertex_flags[4])
    {
        n_edges++;
        e1 = true;
    }
    if (vertex_flags[1] && vertex_flags[5])
    {
        n_edges++;
        e2 = true;
    }
    if (vertex_flags[2] && vertex_flags[6])
    {
        n_edges++;
    }
    if (vertex_flags[3] && vertex_flags[7])
    {
        n_edges++;
    }

    //Set the variable gnp_case with the corresponding case
    if (n_edges == 2 && e1 != e2)
    {
        //Case 4
        gnp_case = 4;
    }
    else
    {
        //Case is the same as n_edges
        gnp_case = n_edges;
    }

    return 1;
}
