//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Given a GNP that is partially inside the sample and that has undergone 
//				some deformation, reconstruct a parallelepiped with a squared base that
//				fits that given GNP
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "GNP_Reconstruction.h"

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
    
    //Variable to store the normal vector of the top surface of the GNP
    //This is neded to reconstruct the rotation matrix of the GNP
    Point_3D N_top;

    //Reconstruct depending on the case
    switch (gnp_case)
    {
    case 3:
        //3 consecutive short edges
        if (!Three_short_edges(vertex_flags, gnp_i, N_top))
        {
            hout << "Error in Reconstruct_partial_gnp when calling Three_short_edges" << endl;
            return 0;
        }
        break;
    case 2:
        //2 consecutive short edges
        if (!Two_consecutive_short_edges(vertex_flags, gnp_i, N_top))
        {
            hout << "Error in Reconstruct_partial_gnp when calling Three_short_edges" << endl;
            return 0;
        }
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
//This function reconstructs a GNP partially inside the sample for the case when
//three short edges are inside the sample
int GNP_Reconstruction::Three_short_edges(const vector<bool>& vertex_flags, GNP& gnp_i, Point_3D& N_top)const
{
    //Vertices for the reference edge
    int R1, R2;

    //Vertices for a reconstructed face (edge to the "left" of the reference edge)
    int LT, LB;

    //Get the vertices of the reference edge
    if (!Get_reference_edge_case3(vertex_flags, R1, R2, LT, LB))
    {
        hout << "Error in Three_short_edges when calling Get_reference_edge" << endl;
        return 0;
    }

    //Adjust the vertices to two paralel planes that define the thickness of the GNP
    //In the process the thickness of the GNP is obtained
    //The normal vector to the top plane is also obtained during the calculations
    if (!Adjust_vertices_along_thickness_planes_case3(vertex_flags, R1, R2, gnp_i, N_top))
    {
        hout << "Error in Three_short_edges when calling Adjust_vertices_along_thickness_planes" << endl;
        return 0;
    }

    //At this point all vertices are in two parallel planes that corespond to the
    //square face

    //Reconstruct face with vertices R1, R2, LT, LB and its parallel face
    //In the proces find a first approximation to the GNP side length
    if (!Adjust_thin_faces_case3(vertex_flags, R1, R2, LT, LB, gnp_i))
    {
        hout << "Error in Three_short_edges when calling Adjust_thin_faces_case3" << endl;
        return 0;
    }
    
    //With the information obtained, the GNP can be reconstructed
    //That is, calculate the coordiantes of all eight vertices of the GNP
    if (!Calculate_gnp_vertices_case3(vertex_flags, R1, R2, LT, LB, SIGABRT, N_top, gnp_i))
    {
        hout << "Error in Three_short_edges when calling Calculate_gnp_vertices_case3" << endl;
        return 0;
    }

    return 1;
}
//This function finds the reference edge for the case of a GNP that has only three short edges
//inside the sample
//The edge is given by the vertex numbers R1 and R2
//The vertices for the edge to the "left" of the edge are also obtained
int GNP_Reconstruction::Get_reference_edge_case3(const vector<bool>& vertex_flags, int& R1, int& R2, int& LT, int& LB)const
{
    //Find the top reference vertex
    if (vertex_flags[0] && vertex_flags[1] && vertex_flags[2] && !vertex_flags[3] && !vertex_flags[7])
    {
        //Reference edge vertices
        R1 = 1; R2 = 5;

        //"Left" edge vertices
        LT = 0; LB = 4;
    }
    else if (vertex_flags[1] && vertex_flags[2] && vertex_flags[3] && !vertex_flags[0] && !vertex_flags[4])
    {
        //Reference edge vertices
        R1 = 2; R2 = 6;

        //"Left" edge vertices
        LT = 1; LB = 5;
    }
    else if (vertex_flags[2] && vertex_flags[3] && vertex_flags[0] && !vertex_flags[1] && !vertex_flags[5])
    {
        //Reference edge vertices
        R1 = 3; R2 = 7;

        //"Left" edge vertices
        LT = 2; LB = 6;
    }
    else if (vertex_flags[3] && vertex_flags[0] && vertex_flags[1] && !vertex_flags[2] && !vertex_flags[6])
    {
        //Reference edge vertices
        R1 = 0; R2 = 4;

        //"Left" edge vertices
        LT = 3; LB = 7;
    }
    else
    {
        hout << "Error when finding the reference edge for the case of three short edges. Boolean vector of vertices indicates that the vertices that shpuld be present are not there." << endl;
        hout << "Boolean vector: " << endl;
        for (size_t i = 0; i < vertex_flags.size(); i++)
            hout << "vertex_flags[" << i << "]=" << vertex_flags[i] << endl;
        return 0;
    }

    return 1;
}
//This function adjusts all vertices on two paralel planes for the squared faces
int GNP_Reconstruction::Adjust_vertices_along_thickness_planes_case3(const vector<bool>& vertex_flags, const int& R1, const int& R2, GNP& gnp_i, Point_3D& N)const
{
    //Get midpoint along reference edge
    Point_3D M = (gnp_i.vertices[R1] + gnp_i.vertices[R2]) * 0.5;

    //Vector to store distances to all vertices
    vector<double> t_halfs(8, 0.0);

    //Get the distance from midpoint to a reference vertex
    //This will be the half of the GNP thickness
    double t_half = M.distance_to(gnp_i.vertices[R1]);
    t_halfs[0] = t_half;

    //Unit vector along refernce edge
    N = (gnp_i.vertices[R1] - gnp_i.vertices[R2]).unit();

    //Iterate over the GNP vertices to find the one furthest from M along normal vector N
    for (int i = 0; i < 8; i++)
    {
        //Avoid vertices outside the sample and reference edge
        if (vertex_flags[i] && i != R1 && i != R2)
        {
            //Calculate the distance from M to vertex i along normal vector N
            //Store it in the vector of distances
            t_halfs[i] = abs(N.dot(gnp_i.vertices[R1] - M));

            //Check if t_half needs updating
            if (t_halfs[i] - t_half > Zero)
            {
                //Update t_half
                t_half = t_halfs[i];
            }
        }
    }
    
    //Update the thickness of the GNP
    gnp_i.t = 2.0*t_half;

    //Adjust top vertices to be at distance t_half from M along normal N
    for (int i = 0; i < 4; i++)
    {
        //Avoid vertices outside the sample
        if (vertex_flags[i])
        {
            //Calculate new location of vertex
            gnp_i.vertices[i] = gnp_i.vertices[i] + N * (t_half - t_halfs[i]);
        }
    }

    //Reverse normal vector N to adjust bottom vertices
    N.set(-N.x, -N.y, -N.z);

    //Adjust bottom vertices to be at distance t_half from M along U
    for (int i = 4; i < 8; i++)
    {
        //Avoid vertices outside the sample
        if (vertex_flags[i])
        {
            //Calculate new location of vertex
            gnp_i.vertices[i] = gnp_i.vertices[i] + N * (t_half - t_halfs[i]);
        }
    }

    return 1;
}
//This function adjusts the this faces of a GNP
int GNP_Reconstruction::Adjust_thin_faces_case3(const vector<bool>& vertex_flags, const int& R1, const int& R2, const int& LT, const int& LB, GNP& gnp_i)const
{
    //Plane that will contain vertices R1R2LT
    Plane_3D Pl_thin;

    //Get the vertex to the right (RT = Right Top)
    int RT = (R1 + 1) % 4;

    //Reconstruct a plane that contains vertices R1, R2, LT, and LB
    if (!Get_thin_face_case3(R1, R2, LT, LB, RT, gnp_i, Pl_thin))
    {
        hout<<"Error in Adjust_thin_faces_case3 when calling Get_thin_face_case3"<<endl;
        return 0;
    }

    //At this point R1, R2, LT, and LB are on the same plane

    //Array with remaining vertices
    int other_v[] = { RT, (R1 + 2) % 4, (R2 + 1) % 4, (R2 + 2) % 4 };
    
    //Variable to store a first approximation to GNP side length
    double l1 = 0.0;

    //Set all remaining vertices on a plane parallel to Pl_thin
    //Also get an approximation to the GNP side length
    if (!Get_parallel_thin_face_case3(vertex_flags, other_v, R1, R2, LT, LB, RT, Pl_thin, gnp_i, l1))
    {
        hout<<"Error in Adjust_thin_faces_case3 when calling Get_parallel_thin_face_case3"<<endl;
        return 0;
    }
    
    //Get a first approximation of the GNP side length
    gnp_i.l = l1;

    //At this point, there are two parallel thin faces

    return 1;
}
//This function gets an initial thin face to reconstruct a GNP for case 3
int GNP_Reconstruction::Get_thin_face_case3(const int& R1, const int& R2, const int& LT, const int& LB, const int& RT, GNP& gnp_i, Plane_3D& Pl_thin)const
{
    //Create a plane with vertices R1R2LT
    Pl_thin = Plane_3D(gnp_i.vertices[R1], gnp_i.vertices[R2], gnp_i.vertices[LT]);

    //Make sure the normal vector goes outside the GNP
    if (Pl_thin.N.dot(gnp_i.vertices[RT] - gnp_i.vertices[R1]) > Zero )
    {
        //Normal goes inside GNP, so reverse it
        Pl_thin.N = Pl_thin.N*(-1.0);
    }

    //Check if LB is outside the GNP
    if (Pl_thin.N.dot(gnp_i.vertices[LB] - gnp_i.vertices[R1]) > Zero)
    {
        //LB is outside the GNP
        //Reset plane to be R1R2LB
        Pl_thin = Plane_3D(gnp_i.vertices[R1], gnp_i.vertices[R2], gnp_i.vertices[LB]);

        //Make sure the normal vector goes outside the GNP
        if (Pl_thin.N.dot(gnp_i.vertices[RT] - gnp_i.vertices[R1]) > Zero)
        {
            //Normal goes inside GNP, so reverse it
            Pl_thin.N = Pl_thin.N * (-1.0);
        }

        //Get distance from LT to plane R1R2LB
        double d_lt = Pl_thin.distance_to(gnp_i.vertices[LT]);

        //Move LT towards plane R1R2LB
        gnp_i.vertices[LT] = gnp_i.vertices[LT] + Pl_thin.N * d_lt;
    }
    else
    {
        //Get distance from LB to plane R1R2LT
        double d_lb = Pl_thin.distance_to(gnp_i.vertices[LB]);

        //Move LB towards plane R1R2LT
        gnp_i.vertices[LB] = gnp_i.vertices[LB] + Pl_thin.N * d_lb;
    }
    
    return 1;
}
//This function find a second thin face, which is parallel to the first one
int GNP_Reconstruction::Get_parallel_thin_face_case3(const vector<bool>& vertex_flags, const int other_v[], const int& R1, const int& R2, const int& LT, const int& LB, const int& RT, const Plane_3D& Pl_thin, GNP& gnp_i, double& l1)const
{
    //Vector to store distances from plane Pl_thin to remaining vertices
    vector<double> dists(4, 0.0);

    //Initialize maximum distance with distance to vertex RT
    l1 = Pl_thin.distance_to(gnp_i.vertices[RT]);
    dists[0] = l1;

    //Vertex with maximum distance
    int idx = 0;

    //Iterate over remaining vertices
    for (int i = 1; i < 4; i++)
    {
        //Get vertex number
        int v = other_v[i];

        //Check if vertex v is present
        if (vertex_flags[v])
        {
            //Get the distance from plane to vertex v
            double d_new = Pl_thin.distance_to(gnp_i.vertices[v]);

            //Store distance to vertex v
            dists[i] = d_new;

            //Check if a new maximum distance was found
            if (d_new - l1 > Zero)
            {
                //Update maximum distance and index
                l1 = d_new;
                idx = i;
            }
        }
    }

    //Set all remaining vertices at a distance l1
    for (int i = 1; i < 4; i++)
    {
        //Ignore the vertex at maximum distance
        if (i != idx)
        {
            //Get vertex number
            int v = other_v[i];

            //Check if vertex v is present
            if (vertex_flags[v])
            {
                //Se the vertex to its new position
                //Need the normal going to the opposite direction and move the vertex a
                //distance l1 - dists[i]
                //Instead of creating a new vector, I multiply the distance by -1 and
                //move the vertex a "negative" distance dists[i] - l1 in the
                //direction outside of the GNP
                //In this way the vertex is move in the direction inside the GNP
                gnp_i.vertices[v] = gnp_i.vertices[v] + Pl_thin.N * (dists[i] - l1);
            }
        }
    }
    
    return 1;
}
//This function calculates the eight vertices of a GNP for case 3 (three short
//vertices inside the sample)
int GNP_Reconstruction::Calculate_gnp_vertices_case3(const vector<bool>& vertex_flags, const int& R1, const int& R2, const int& LT, const int& LB, const int& RT, const Point_3D& N_top, GNP& gnp_i)const
{
    //Get the midpoint of edge R1LT
    Point_3D M = (gnp_i.vertices[R1] + gnp_i.vertices[LT]) * 0.5;

    //Get unit vector R1LT
    Point_3D V = (gnp_i.vertices[LT] - gnp_i.vertices[R1]).unit();

    //Variable to store second approximation of GNP side length
    double l2_half = 0.0;
    
    //Initialize maximum distance along V with distance to vertex R1
    l2_half = abs(V.dot(gnp_i.vertices[R1] - M));

    //Iterate over all vertices and find the furthest from M along vector U
    for (int i = 0; i < 8; i++)
    {
        //Check that the vertex is present and it is not R1
        if (vertex_flags[i] && i != R1)
        {
            //Calculate distance to vertex i along V
            double t_new = abs(V.dot(gnp_i.vertices[i] - M));
            
            //Check if t_new is larger than previous maximum
            if (t_new - l2_half > Zero)
            {
                //Update maximum and index
                l2_half = t_new;
            }
        }
    }
    
    //Get the GNP side length
    //Recall gnp_i.l cunrrently is the first approximation of the GNP side length
    gnp_i.l = max(gnp_i.l, 2.0*l2_half);
    
    //Half GNP side length
    double l_half = 0.5*gnp_i.l;
    
    //Using M and the GNP side length, reconstruct the GNP vertices
    
    //Reconstruct vertices along edge R1LT
    gnp_i.vertices[R1] = M - V*l_half;
    gnp_i.vertices[LT] = M + V*l_half;
    
    //Displacement to go from a top edge on a squared face to the edge below it
    Point_3D disp_down = N_top*gnp_i.t;
    
    //Reconstruct vertices on edge "below" edge R1LT
    gnp_i.vertices[R2] = gnp_i.vertices[R1] - disp_down;
    gnp_i.vertices[LB] = gnp_i.vertices[LT] - disp_down;
    
    //Get unit vector that is normal to both N_top and V
    Point_3D N_surf = V.cross(N_top);
    
    //Displacement for remaining vertices from the four already calculated
    Point_3D disp_surf = N_surf*gnp_i.l;
    
    //Reconstruct vertices to the right of reference edge
    gnp_i.vertices[RT] = gnp_i.vertices[R1] + disp_surf;
    gnp_i.vertices[RT + 4] = gnp_i.vertices[R2] + disp_surf;
    
    //Last vertices
    int last = (RT + 1) % 4;
    gnp_i.vertices[last] = gnp_i.vertices[LT] + disp_surf;
    gnp_i.vertices[last + 4] = gnp_i.vertices[LB] + disp_surf;
    
    return 1;
}
//This function reconstructs a GNP partially inside the sample for the case when
//two consecutive short edges are inside the sample
int GNP_Reconstruction::Two_consecutive_short_edges(const vector<bool>& vertex_flags, GNP& gnp_i, Point_3D& N_top)const
{
    //Variables to store the reference plane, i.e., these are the reference vertices
    int R1, R2, R3, R4;
    
    //Variables to store other vertices if present
    //Initialized to -1 to indicate they are not present
    int O1 = -1, O2 = -1;
    
    //Get the actual values of the refrence vertices
    if (!Get_reference_vertices_case2(vertex_flags, R1, R2, R3, R4, O1, O2))
    {
        hout<<"Error in Two_consecutive_short_edges when calling Get_reference_vertices_case2"<<endl;
        return 0;
    }
    
    //Get all vertices R1-R4 on a reference plane
    if (!Get_reference_plane_case2(R1, R2, R3, R4, gnp_i))
    {
        hout<<"Error in Two_consecutive_short_edges when calling Get_reference_plane_case2"<<endl;
        return 0;
    }
    
    //Set the edges R1R4 and R2R3 to be parallel
    //In the process, the GNP thickness is found
    if (!Set_long_reference_edges_as_parallel_case2(R1, R2, R3, R4, O1, O2, gnp_i))
    {
        hout<<"Error in Two_consecutive_short_edges when calling Set_long_reference_edges_as_parallel_case2"<<endl;
        return 0;
    }
    
    //Find the GNP side length and calculate GNP vertices
    if (!Find_gnp_length_and_calculate_vertices_case2(R1, R2, R3, R4, O1, O2, gnp_i, N_top)) {
        hout<<"Error in Two_consecutive_short_edges when calling Find_gnp_length_and_calculate_vertices_case2"<<endl;
        return 0;
    }
    
    return 1;
}
//This function finds the reference vertices for the case of
//two consecutive short edges inside the sample
int GNP_Reconstruction::Get_reference_vertices_case2(const vector<bool>& vertex_flags, int& R1, int& R2, int& R3, int& R4, int& O1, int& O2)const
{
    //Find the thin face that is actually present
    if (vertex_flags[0] && vertex_flags[4])
    {
        if (vertex_flags[1] && vertex_flags[5])
        {
            
            //Consecutive short edges are 04 and 15
            R1 = 0; R2 = 4; R3 = 5; R4 = 1;
            
            //Check if other vertices are present
            //Only one vertex from 3 and 7 may be present
            O1 = (vertex_flags[3])? 3: (vertex_flags[7])? 7: -1;
            //Only one vertex from 2 and 6 may be present
            O2 = (vertex_flags[2])? 2: (vertex_flags[6])? 6: -1;
        }
        else if (vertex_flags[3] && vertex_flags[7])
        {
            //Consecutive short edges are 37 and 04
            R1 = 3; R2 = 7; R3 = 4; R4 = 0;
            
            //Check if other vertices are present
            //Only one vertex from 2 and 6 may be present
            O1 = (vertex_flags[2])? 2: (vertex_flags[6])? 6: -1;
            //Only one vertex from 1 and 5 may be present
            O2 = (vertex_flags[1])? 1: (vertex_flags[5])? 5: -1;
        }
        else {
            hout<<"Error in Get_reference_vertices_case2. Short edge 04 is present, but neither 37 or 15 are present."<<endl;
            hout << "Boolean vector: " << endl;
            for (size_t i = 0; i < vertex_flags.size(); i++)
                hout << "vertex_flags[" << i << "]=" << vertex_flags[i] << endl;
            return 0;
        }
    }
    else if (vertex_flags[2] && vertex_flags[6])
    {
        if (vertex_flags[1] && vertex_flags[5])
        {
            //Consecutive short edges are 15 and 26
            R1 = 1; R2 = 5; R3 = 6; R4 = 2;
            
            //Check if other vertices are present
            //Only one vertex from 0 and 4 may be present
            O1 = (vertex_flags[0])? 0: (vertex_flags[4])? 4: -1;
            //Only one vertex from 3 and 7 may be present
            O2 = (vertex_flags[3])? 3: (vertex_flags[7])? 7: -1;
        }
        else if (vertex_flags[3] && vertex_flags[7])
        {
            //Consecutive short edges are 26 and 37
            R1 = 2; R2 = 6; R3 = 7; R4 = 3;
            
            //Check if other vertices are present
            //Only one vertex from 1 and 5 may be present
            O1 = (vertex_flags[1])? 1: (vertex_flags[5])? 5: -1;
            //Only one vertex from 0 and 4 may be present
            O2 = (vertex_flags[0])? 0: (vertex_flags[4])? 4: -1;
            
        }
        else {
            hout<<"Error in Get_reference_vertices_case2. Short edge 26 is present, but neither 37 or 15 are present."<<endl;
            hout << "Boolean vector: " << endl;
            for (size_t i = 0; i < vertex_flags.size(); i++)
                hout << "vertex_flags[" << i << "]=" << vertex_flags[i] << endl;
            return 0;
        }
    }
    else {
        hout<<"Error in Get_reference_vertices_case2. The conditions for the case of two consecutive short edges inside the sample are not met."<<endl;
        hout << "Boolean vector: " << endl;
        for (size_t i = 0; i < vertex_flags.size(); i++)
            hout << "vertex_flags[" << i << "]=" << vertex_flags[i] << endl;
        return 0;
    }
    
    return 1;
}
//This function finds the reference plane for the case of
//two consecutive short edges inside the sample
int GNP_Reconstruction::Get_reference_plane_case2(const int& R1, const int& R2, const int& R3, const int& R4, GNP& gnp_i)const
{
    //Vector R1R2
    Point_3D R1R2 = gnp_i.vertices[R2] - gnp_i.vertices[R1];
    
    //Normal vector going inside the GNP (considering the proper numbering of vertices Ri)
    //This vector is R1R4xR1R2
    Point_3D N_in = (gnp_i.vertices[R4] - gnp_i.vertices[R1]).cross(R1R2);
    N_in.make_unit();
    
    //Get signed distance of R4 to plane R1R2R4
    double d_sign = (gnp_i.vertices[R3] - gnp_i.vertices[R4]).dot(N_in);
    
    //Check if vertex R3 is outside the GNP
    if (d_sign < Zero)
    {
        //R3 is outside the GNP, update normal vector
        N_in = (gnp_i.vertices[R3] - gnp_i.vertices[R1]).cross(R1R2);
        N_in.make_unit();
        
        //Get distance from R4 to plane R1R2R3
        double d = abs((gnp_i.vertices[R3] - gnp_i.vertices[R4]).dot(N_in));
        
        //Project R4 onto plane
        //Since N_in goes inside, and need to move R4 towards the outside,
        //the negative sign is used
        gnp_i.vertices[R4] = gnp_i.vertices[R4] - N_in*d;
    }
    else
    {
        //R3 is inside the GNP, project it onto the plane R1R2R4
        //The signed distance of R3 is positive since the vector R4R3 also goes
        //towards the inside of the GNP
        //However, since I need to move R3 towards the outside, I use a negative sign
        gnp_i.vertices[R3] = gnp_i.vertices[R3] - N_in*d_sign;
    }
    
    return 1;
}
//This function sets the edges R1R4 and R2R3 as parallel and in the precess the GNP
//thickness is calculated
//If other vertices are present, they are also considered to calculate the distance that
//separates edges R1R4 and R2R3 and, thus, the GNP thickness
int GNP_Reconstruction::Set_long_reference_edges_as_parallel_case2(const int& R1, const int& R2, const int& R3, const int& R4, const int& O1, const int& O2, GNP& gnp_i)const
{
    //Get midpoint of edge R1R2
    Point_3D M = (gnp_i.vertices[R1] + gnp_i.vertices[R2])*0.5;
    
    //Get unit vector along R1R2
    Point_3D U = (gnp_i.vertices[R2] - gnp_i.vertices[R1]).unit();
    
    //Variable to store distances
    //It is initialized with the distance from M to R1 and R2, which is the same
    //since M is the midpoint of R1R2
    //Then I only update two elements of vector dists
    vector<double> dists(4, M.distance_to(gnp_i.vertices[R2]));
    
    //Get distance from M to R3 and R4 along U
    dists[2] = U.dot(gnp_i.vertices[R3] - M);
    dists[3] = U.dot(M - gnp_i.vertices[R4]);
    
    //Find maximum distance and index
    double d_max = dists[1];
    int idx = 1;
    
    if (dists[2] - d_max > Zero)
    {
        //Update maximum distance and index
        d_max = dists[2];
        idx = 2;
    }
    if (dists[3] - d_max > Zero)
    {
        //Update maximum distance and index
        d_max = dists[3];
        idx = 3;
    }
    
    //Variables for the distances to other vertices, if any
    double do1 = 0.0, do2 = 0.0;
    
    //Check if there are other vertices that might affect the calculation of
    //the GNP thickness
    if (O1 != -1)
    {
        //Calculate the distance to O1 along U
        do1 = abs(U.dot(M - gnp_i.vertices[O1]));
        
        //Check if maximum distance needs to be updated
        if (do1 - d_max > Zero)
        {
            d_max = do1;
            idx = -1;
        }
    }
    if (O2 != -1)
    {
        //Calculate the distance to O2 along U
        do2 = abs(U.dot(M - gnp_i.vertices[O2]));
        
        //Check if maximum distance needs to be updated
        if (do2 - d_max > Zero)
        {
            d_max = do2;
            idx = -2;
        }
    }
    
    //Now that the maximum distance has been determined, move the vertices
    if (idx == 1)
    {
        //Only adjust R3 and R4
        gnp_i.vertices[R3] = gnp_i.vertices[R3] + U*(d_max - dists[2]);
        gnp_i.vertices[R4] = gnp_i.vertices[R4] + U*(dists[3] - d_max);
    }
    else
    {
        //Adjust all except idx
        gnp_i.vertices[R1] = M - U*d_max;
        gnp_i.vertices[R2] = M + U*d_max;
        if (idx != 2)
        {
            //Adjust R3 since maximum is not at R3
            gnp_i.vertices[R3] = gnp_i.vertices[R3] + U*(d_max - dists[2]);
        }
        if (idx != 3)
        {
            //Adjust R4 since maximum is not at R4
            gnp_i.vertices[R4] = gnp_i.vertices[R4] + U*(dists[3] - d_max);
        }
    }
    
    //Adjust O1 and O2 if present
    if (O1 != -1) {
        //Get the displacement
        Point_3D S = (O1 <= 3)? U*(do1 - d_max): U*(d_max - do1);
        gnp_i.vertices[O1] = gnp_i.vertices[O1] + S;
    }
    if (O2 != -1) {
        //Get the displacement
        Point_3D S = (O2 <= 3)? U*(do2 - d_max): U*(d_max - do2);
        gnp_i.vertices[O2] = gnp_i.vertices[O2] + S;
    }
    
    //Set the GNP thickness
    gnp_i.t = 2.0*d_max;
    
    return 1;
}
//This function finds the GNP side length and calculates the GNP vertices
int GNP_Reconstruction::Find_gnp_length_and_calculate_vertices_case2(const int& R1, const int& R2, const int& R3, const int& R4, const int& O1, const int& O2, GNP& gnp_i, Point_3D& N_top)const
{
    //Get midpoint of edge R1R4
    Point_3D M = (gnp_i.vertices[R1] + gnp_i.vertices[R4])*0.5;
    
    //Get unit vector along R1R4
    Point_3D V = (gnp_i.vertices[R4] - gnp_i.vertices[R1]).unit();
    
    //Vector with all distances
    //Initialized with the distance from M to R1,
    //which is the also the distance from M to R2
    vector<double> dists(6, M.distance_to(gnp_i.vertices[R1]));
    
    //Initialize maximum distance 1 and index of maximum distance
    double d_max = dists[0];
    int idx = 0;
    
    //Calcualte distances to R2 and R3
    dists[1] = V.dot(M - gnp_i.vertices[R2]);
    dists[2] = V.dot(gnp_i.vertices[R3] - M);
    
    //Check if maximum distance needs updating
    if (dists[1] - d_max > Zero)
    {
        //Update maximum distance and index
        d_max = dists[1];
        idx = 1;
    }
    if (dists[2] - d_max > Zero)
    {
        //Update maximum distance and index
        d_max = dists[2];
        idx = 2;
    }
    
    //Variables for the distances to other vertices, if any
    double do1 = 0.0, do2 = 0.0;
    
    //Check if there are other vertices that might affect the calculation of
    //the GNP thickness
    if (O1 != -1)
    {
        //Calculate the distance to O1 along V
        do1 = abs(V.dot(M - gnp_i.vertices[O1]));
        
        //Check if maximum distance needs to be updated
        if (do1 - d_max > Zero)
        {
            d_max = do1;
            idx = -1;
        }
    }
    if (O2 != -1)
    {
        //Calculate the distance to O2 along V
        do2 = abs(V.dot(M - gnp_i.vertices[O2]));
        
        //Check if maximum distance needs to be updated
        if (do2 - d_max > Zero)
        {
            d_max = do2;
            idx = -2;
        }
    }
    
    //Now that the maximum distance has been determined, set the GNP side length
    gnp_i.l = 2.0*d_max;
    
    //Calculate the eight GNP vertices
    
    //Update R1 and R4
    gnp_i.vertices[R1] = M - V*d_max;
    gnp_i.vertices[R4] = M + V*d_max;
    
    //Get a vector going inside the GNP through the refrence face R1R2R3R4
    Point_3D N_in = V.cross(gnp_i.vertices[R3] - gnp_i.vertices[R1]);
    N_in.make_unit();
    
    //Displacement from reference face
    Point_3D disp_face = N_in*gnp_i.l;
    
    //Calculate remaining two vertices on top face
    idx = (R4 + 1) % 4;
    gnp_i.vertices[idx] = gnp_i.vertices[R4] + disp_face;
    idx = (R4 + 2) % 4;
    gnp_i.vertices[idx] = gnp_i.vertices[R1] + disp_face;
    
    //Calculate vector normal to top face
    //Both V and N_in are normal and unitary, thus its dot product is a unitary vector
    //Thus there is no need to normalize N_top
    N_top = V.cross(N_in);
    
    //Calculate displacemente though the thickness
    Point_3D disp_thick = N_top*gnp_i.t;
    
    //Calculate vertices in bottom face
    gnp_i.vertices[R2] = gnp_i.vertices[R1] - disp_thick;
    gnp_i.vertices[R3] = gnp_i.vertices[R4] - disp_thick;
    int idx_top = (R4 + 1) % 4;
    int idx_bot = (R3 + 1) % 4;
    gnp_i.vertices[idx_bot] = gnp_i.vertices[idx_top] - disp_thick;
    idx_top = (R4 + 2) % 4;
    idx_bot = (R3 + 2) % 4;
    gnp_i.vertices[idx_bot] = gnp_i.vertices[idx_top] - disp_thick;
    
    return 1;
}
