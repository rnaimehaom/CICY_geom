//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Given a GNP that is partially inside the sample and that has undergone 
//				some deformation, reconstruct a parallelepiped with a squared base that
//				fits that given GNP
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "GNP_Reconstruction.h"

//This function reconstructs a GNP
//This function determines whether a full GNP or a partial GNP is reconstructed
int GNP_Reconstruction::Reconstruct_gnp(const vector<int>& vertices, const vector<bool>& vertex_flags, GNP& gnp_i)const
{
    //Initial gnp
    GNP gnp0 = gnp_i;

    //Variable to store the normal vector of the top surface of the GNP
    //This is neded to reconstruct the rotation matrix of the GNP
    Point_3D N_top;

    //Determine whether a full GNP or a partial one is to bre reconstructed
    if (vertices.size() == 8)
    {
        //Reconstruct a full GNP
        if (!Reconstruct_full_gnp(gnp_i, N_top))
        {
            hout << "Error in Reconstruct_gnp when calling Reconstruct_full_gnp" << endl;
            return 0;
        }
    }
    else
    {
        //Reconstruct a partial GNP
        if (!Reconstruct_partial_gnp(vertex_flags, gnp_i, N_top))
        {
            hout << "Error in Reconstruct_gnp when calling Reconstruct_partial_gnp" << endl;
            return 0;
        }
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
    if (!GN.Rotation_matrix_from_direction(N_top.x, N_top.y, N_top.z, gnp_i.rotation))
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

    //Check if all initial points are actually inside the GNP
    bool terminate = false;
    stringstream ss;
    for (size_t i = 0; i < vertex_flags.size(); i++)
    {
        if (vertex_flags[i])
        {
            if (gnp_i.Is_point_inside_gnp(gnp0.vertices[i]))
            {
                ss << i << " inside GNP ";
            }
            else
            {
                ss << i << " outside GNP ";
                terminate = true;
            }
            //ss << gnp_i.flag << " O=" << gnp0.vertices[i].str(15) << " N=" << gnp_i.vertices[i].str(15) << endl;
            ss << endl;
        }
    }
    if (terminate)
    {
        hout << "Displaced point is outside reconstructed GNP" << endl;
        hout << ss.str() << endl;
        return 0;
    }


    return 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//This function reconstrucst a GNP after it has undergone a deformation
//Here, the paralellepiped that contains the eight vertices of the GNP on 
//its deformed position is found
int GNP_Reconstruction::Reconstruct_full_gnp(GNP& gnp_i, Point_3D& N_top)const
{
    //Plane that contains vertices v0 to v3
    Plane_3D Pl_top;

    //Get the equation of the plane that contains vertices v0 to v3
    //hout << "Get_plane_from_top_squared_face" << endl;
    if (!Get_plane_from_top_squared_face(gnp_i, Pl_top))
    {
        hout << "Error in Reconstruct_full_gnp when calling Get_plane_from_top_squared_face" << endl;
        return 0;
    }

    //Variable to store the distance between planes
    double d_planes = 0.0;

    //Use the top plane (and the GNP vertices) to find the bottom plane
    //hout << "Get_plane_from_bottom_squared_face" << endl;
    if (!Get_plane_from_bottom_squared_face(Pl_top, gnp_i, d_planes))
    {
        hout << "Error in Reconstruct_full_gnp when calling Get_plane_from_bottom_squared_face" << endl;
        return 0;
    }

    //Variable to get first approximation of GNP side length
    double lx = 0.0;

    //Set parallel planes v0v1v5v4 and v2v6v7v3
    //hout << "Set_parallel_planes_along_x" << endl;
    if (!Set_parallel_planes_along_x(gnp_i, lx))
    {
        hout << "Error in Reconstruct_full_gnp when calling Set_parallel_planes_along_x" << endl;
        return 0;
    }

    //Set parallel planes v0v4v7v3 and v1v5v6v2 and recalculate all vetex coordinates
    //hout << "Set_parallel_planes_along_y_and_recalculate_vertices" << endl;
    if (!Set_parallel_planes_along_y_and_recalculate_vertices(Pl_top, d_planes, lx, gnp_i))
    {
        hout << "Error in Reconstruct_full_gnp when calling Set_parallel_planes_along_y_and_recalculate_vertices" << endl;
        return 0;
    }

    //Set the normal vector to the top face
    N_top = Pl_top.N;

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
int GNP_Reconstruction::Get_plane_from_bottom_squared_face(const Plane_3D& Pl_top, GNP& gnp_i, double& d_planes)const
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

    //Iterate over all vertices again to check if they need to be projected or not
    for (int j = 4; j <= 7; j++)
    {
        //Ignore vertex idx that was used to calculate the parallel plane
        if (j != idx)
        {
            //Move vertex j-4 towards the bottom plane to the location of vertex j
            gnp_i.vertices[j] = gnp_i.vertices[j] - Pl_top.N * (d_max - dists[j - 4]);
        }
    }

    return 1;
}
//Set parallel planes v0v1v5v4 and v2v6v7v3
int GNP_Reconstruction::Set_parallel_planes_along_x(GNP& gnp_i, double& lx)const
{
    //Set v3v0 as reference edge, get the unit vector along that edge
    Point_3D U = (gnp_i.vertices[0] - gnp_i.vertices[3]).unit();

    //Get midpoint of reference edge
    Point_3D M = (gnp_i.vertices[0] + gnp_i.vertices[3]) * 0.5;

    //Vector to store all distances in the direction of vector U
    vector<double> dists(8, 0.0);

    //Initialize maximum distance with distance from midpoint to any vertex
    //of the reference edge
    double d_max = M.distance_to(gnp_i.vertices[0]);
    //Index of maximum distance
    int idx = 0;

    //Calculate all the remaining distances 
    //Vertex 0 is ignored since distance to that vertices is known
    for (size_t i = 1; i < dists.size(); i++)
    {
        //Ignore vertex 3 since distance to that vertices is known
        if (i != 3)
        {
            //Calculate distance along vector U
            double new_d = abs( U.dot(M - gnp_i.vertices[i]) );

            //Save distance
            dists[i] = new_d;

            //Check if new maximum was found
            if (new_d - d_max > Zero)
            {
                //Update maximum distance and index
                d_max = new_d;
                idx = (int)i;
            }
        }
    }

    //Set the approximated GNP length along x
    lx = 2.0 * d_max;

    //Set all point to maximum distance from M

    //Check if maximum was not found in reference edge
    if (idx != 0)
    {
        //Displacement on reference edge form midpoint
        Point_3D M_disp = U * d_max;
        
        //Update vertices 0 and 3
        gnp_i.vertices[0] = M + M_disp;
        gnp_i.vertices[3] = M - M_disp;
    }

    //Vertices for wich a positive sign is used, ignoring the verte in the reference edge
    int pos[] = {1, 4, 5};

    //Update vertices in array pos
    for (int i = 0; i < 3; i++)
    {
        //Get vertex number
        int v = pos[i];

        //Ignore vertex at maximum distance from M
        if (v != idx)
        {
            //Calculate new position
            gnp_i.vertices[v] = gnp_i.vertices[v] + U * (d_max - dists[v]);
        }
    }

    //Vertices for wich a negative sign is used, ignoring the verte in the reference edge
    int neg[] = {2, 6, 7};

    //Update vertices in array pos
    for (int i = 0; i < 3; i++)
    {
        //Get vertex number
        int v = neg[i];

        //Ignore vertex at maximum distance from M
        if (v != idx)
        {
            //Calculate new position
            gnp_i.vertices[v] = gnp_i.vertices[v] - U * (d_max - dists[v]);
        }
    }

    return 1;
}
//Set parallel planes v0v4v7v3 and v1v5v6v2
int GNP_Reconstruction::Set_parallel_planes_along_y_and_recalculate_vertices(const Plane_3D& Pl_top, const double& d_planes, const double& lx, GNP& gnp_i)const
{
    //Set v1v0 as reference edge, get the unit vector along that edge
    Point_3D V = (gnp_i.vertices[0] - gnp_i.vertices[1]).unit();

    //Get midpoint of reference edge
    Point_3D M = (gnp_i.vertices[0] + gnp_i.vertices[1]) * 0.5;

    //Initialize maximum distance with distance from midpoint to any vertex
    //of the reference edge
    double d_max = M.distance_to(gnp_i.vertices[0]);

    //Calculate all the remaining distances 
    //Vertices 0 and 1 are ignored since distance to those vertices is known
    for (size_t i = 2; i < 8; i++)
    {
        //Calculate distance along vector V
        double new_d = abs(V.dot(M - gnp_i.vertices[i]));

        //Check if new maximum was found
        if (new_d - d_max > Zero)
        {
            //Update maximum distance and index
            d_max = new_d;
        }
    }

    //Get the GNP side length
    gnp_i.l = max(lx, 2.0 * d_max);

    //Recalculate all vertex positions

    //Displacement on reference edge form midpoint
    Point_3D M_disp = V * gnp_i.l * 0.5;

    //Update vertices 0 and 1
    gnp_i.vertices[0] = M + M_disp;
    gnp_i.vertices[1] = M - M_disp;

    //Displacement though the thickness
    Point_3D t_disp = Pl_top.N * d_planes;

    //Update vertices 4 and 5
    gnp_i.vertices[4] = gnp_i.vertices[0] - t_disp;
    gnp_i.vertices[5] = gnp_i.vertices[1] - t_disp;

    //Calculate unit vector from face v0v1v5v4 towards face v3v2v6v7
    Point_3D N_in = (gnp_i.vertices[5] - gnp_i.vertices[1]).cross(V).unit();

    //Dsplacement along length
    Point_3D l_disp = N_in * gnp_i.l;

    //Calculate remaining vertices
    gnp_i.vertices[3] = gnp_i.vertices[0] + l_disp;
    gnp_i.vertices[2] = gnp_i.vertices[1] + l_disp;
    gnp_i.vertices[7] = gnp_i.vertices[4] + l_disp;
    gnp_i.vertices[6] = gnp_i.vertices[5] + l_disp;

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
int GNP_Reconstruction::Reconstruct_partial_gnp(const vector<bool>& vertex_flags, GNP& gnp_i, Point_3D& N_top)const
{
    //Variable to store the case for GNP reconstruction
    int gnp_case = -1;

    //Find the case
    if (!Find_reconstruction_case(vertex_flags, gnp_case))
    {
        hout << "Error in Reconstruct_partial_gnp when calling Find_reconstruction_case" << endl;
        return 0;
    }
    //hout << "case " << gnp_case << endl;

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
        if (!One_short_edge(vertex_flags, gnp_i, N_top)) {
            hout << "Error in Reconstruct_partial_gnp when calling One_short_edge" << endl;
            return 0;
        }
        break;
    case 0:
        //0 short edges
        if (!No_short_edges(vertex_flags, gnp_i, N_top))
        {
            hout << "Error in Reconstruct_partial_gnp when calling No_short_edges" << endl;
            return 0;
        }
        break;
    case 4:
        //2 non-consecutive short edges
        if (!Two_non_consecutive_short_edges(vertex_flags, gnp_i, N_top))
        {
            hout << "Error in Reconstruct_partial_gnp when calling Two_non_consecutive_short_edges" << endl;
            return 0;
        }
        break;
    default:
            hout << "Unkonwn case for reconstructing partial GNP. Case: " << gnp_case << endl;
            return 0;
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
    }
    if (vertex_flags[2] && vertex_flags[6])
    {
        n_edges++;
        e2 = true;
    }
    if (vertex_flags[3] && vertex_flags[7])
    {
        n_edges++;
    }

    //Set the variable gnp_case with the corresponding case
    if (n_edges == 2 && e1 == e2)
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

    //Get the vertices of the reference edge
    if (!Get_reference_edge_case3(vertex_flags, R1, R2))
    {
        hout << "Error in Three_short_edges when calling Get_reference_edge" << endl;
        return 0;
    }
    //hout << "R1=" << R1 << " R2=" << R2 << endl;

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

    //Variable to store first approximation of GNP side length
    double l1 = 0.0;

    //Set two thin faces as parallel
    if (!Set_first_pair_of_parallel_thin_faces_case3(vertex_flags, R1, gnp_i, l1))
    {
        hout << "Error in Three_short_edges when calling Set_pair_of_parallel_thin_faces" << endl;
        return 0;
    }

    //Find the GNP side length by setting the other two pair of thin faces parallel
    //and recalculate vertex coordinates
    if (!Find_gnp_side_length_and_recalculate_vertices_case3(vertex_flags, R1, l1, N_top, gnp_i))
    {
        hout << "Error in Three_short_edges when calling Find_gnp_side_length_and_recalculate_vertices_case3" << endl;
        return 0;
    }

    return 1;
}
//This function finds the reference edge for the case of a GNP that has only three short edges
//inside the sample
//The edge is given by the vertex numbers R1 and R2
//The vertices for the edge to the "left" of the edge are also obtained
int GNP_Reconstruction::Get_reference_edge_case3(const vector<bool>& vertex_flags, int& R1, int& R2)const
{
    //Find the top reference vertex by finding the vertex that is not completely inside the sample
    if (!vertex_flags[3] || !vertex_flags[7])
    {
        //Reference edge vertices
        R1 = 1; R2 = 5;
    }
    else if (!vertex_flags[0] || !vertex_flags[4])
    {
        //Reference edge vertices
        R1 = 2; R2 = 6;
    }
    else if (!vertex_flags[1] || !vertex_flags[5])
    {
        //Reference edge vertices
        R1 = 3; R2 = 7;
    }
    else if (!vertex_flags[2] || !vertex_flags[6])
    {
        //Reference edge vertices
        R1 = 0; R2 = 4;
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
    t_halfs[R1] = t_half;
    t_halfs[R2] = t_half;

    //Index of maximum distance
    int idx = R1;

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
            t_halfs[i] = abs(N.dot(gnp_i.vertices[i] - M));

            //Check if t_half needs updating
            if (t_halfs[i] - t_half > Zero)
            {
                //Update t_half and index
                t_half = t_halfs[i];
                idx = i;
            }
        }
    }
    
    //Update the thickness of the GNP
    gnp_i.t = 2.0*t_half;

    //Adjust top vertices to be at distance t_half from M along normal N
    for (int i = 0; i < 4; i++)
    {
        //Avoid vertices outside the sample
        if (vertex_flags[i] && i != idx)
        {
            //Calculate new location of vertex
            gnp_i.vertices[i] = gnp_i.vertices[i] + N * (t_half - t_halfs[i]);
        }
    }

    //Adjust bottom vertices to be at distance t_half from M along U
    for (int i = 4; i < 8; i++)
    {
        //Avoid vertices outside the sample
        if (vertex_flags[i] && i != idx)
        {
            //Calculate new location of vertex
            gnp_i.vertices[i] = gnp_i.vertices[i] - N * (t_half - t_halfs[i]);
        }
    }

    return 1;
}
//This function sets a pair of thin faces as parallel
//One of these thin faces is the one that contains the reference vertex R1 and the vertex to its "left"
int GNP_Reconstruction::Set_first_pair_of_parallel_thin_faces_case3(const vector<bool>& vertex_flags, const int& R1, GNP& gnp_i, double& l1)const
{
    //Calculate the reference vertex to the "left" of R1
    int LT = (R1 + 3) % 4;

    //Calculate the midpoint of edge LTR1
    Point_3D M = (gnp_i.vertices[LT] + gnp_i.vertices[R1]) * 0.5;

    //Calculate the unit vector from LT to R1
    Point_3D U = (gnp_i.vertices[R1] - gnp_i.vertices[LT]).unit();

    //Initialize a maximum distance using one of the reference vertices
    double d_max = M.distance_to(gnp_i.vertices[R1]);
    //Index of maximum distance
    int idx = R1;

    //Vector to store distances, initialized with initial guess for maximum distance
    vector<double> dists(8, d_max);

    //Iterate over all vertices
    for (int i = 0; i < (int)vertex_flags.size(); i++)
    {
        //Ignore reference vertices and those outside the sample
        if (vertex_flags[i] && i != R1 && i != LT)
        {
            //Calculate distance along vector U
            dists[i] = abs(U.dot(M - gnp_i.vertices[i]));

            //Check if maximum distance needs updating
            if (dists[i] - d_max > Zero)
            {
                //Update maximum distance and index
                d_max = dists[i];
                idx = i;
            }
        }
    }

    //Set the first approximation for GNP side length
    l1 = 2.0 * d_max;

    //Update locations of all vertices inside the sample

    //Calculate displacement for reference edges
    Point_3D M_disp = U * d_max;

    //Update reference vertices if not at maximum distance
    if (idx != R1)
    {
        gnp_i.vertices[R1] = M + M_disp;
        gnp_i.vertices[LT] = M - M_disp;
    }

    //Update vertex below reference vertex R1 if not at maximum distance
    if (idx != R1 + 4)
    {
        gnp_i.vertices[R1 + 4] = gnp_i.vertices[R1 + 4] + U * (d_max - dists[R1 + 4]);
    }

    //Update vertex below reference vertex LT if not at maximum distance
    if (idx != LT + 4)
    {
        gnp_i.vertices[LT + 4] = gnp_i.vertices[LT + 4] - U * (d_max - dists[LT + 4]);
    }

    //Update vertex to the right of reference vertex R1 if not at maximum distance
    int idx1 = (R1 + 1) % 4;
    if (idx != idx1)
    {
        gnp_i.vertices[idx1] = gnp_i.vertices[idx1] + U * (d_max - dists[idx1]);
    }

    //Update vertex below vertex idx1 if not at maximum distance
    if (idx != idx1 + 4)
    {
        gnp_i.vertices[idx1 + 4] = gnp_i.vertices[idx1 + 4] + U * (d_max - dists[idx1 + 4]);
    }

    //Check if top vertex of fourth edge is present
    //If so, update vertex if not at maximum distance
    int idxt = (R1 + 2) % 4;
    if (vertex_flags[idxt] && idx != idxt)
    {
        gnp_i.vertices[idxt] = gnp_i.vertices[idxt] - U * (d_max - dists[idxt]);
    }

    //Check if bottom vertex of fourth edge is present
    //If so, update vertex if not at maximum distance
    if (idx != idxt + 4)
    {
        gnp_i.vertices[idxt + 4] = gnp_i.vertices[idxt + 4] - U * (d_max - dists[idxt + 4]);
    }

    return 1;
}
//With the information obtained, the GNP can be reconstructed
//That is, calculate the coordiantes of all eight vertices of the GNP
int GNP_Reconstruction::Find_gnp_side_length_and_recalculate_vertices_case3(const vector<bool>& vertex_flags, const int& R1, const double& l1, const Point_3D& N_top, GNP& gnp_i)const
{
    //Calculate the reference vertex to the "right" of R1
    int RT = (R1 + 1) % 4;

    //Calculate the midpoint of edge RTR1
    Point_3D M = (gnp_i.vertices[RT] + gnp_i.vertices[R1]) * 0.5;

    //Calculate the unit vector from RT to R1
    Point_3D V = (gnp_i.vertices[R1] - gnp_i.vertices[RT]).unit();

    //Initialize a maximum distance using one of the reference vertices
    double d_max = M.distance_to(gnp_i.vertices[R1]);

    //Iterate over all vertices
    for (int i = 0; i < (int)vertex_flags.size(); i++)
    {
        //Ignore reference vertices and those outside the sample
        if (vertex_flags[i] && i != R1 && i != RT)
        {
            //Calculate distance along vector U
            double new_d = abs(V.dot(M - gnp_i.vertices[i]));

            //Check if maximum distance needs updating
            if (new_d - d_max > Zero)
            {
                //Update maximum distance and index
                d_max = new_d;
            }
        }
    }

    //Set the GNP side length
    gnp_i.l = max(l1, 2.0*d_max);

    //Now that the GNP side length is known, reconstruct the GNP
    //Recall the thickness is already known

    //Displacement from M along reference edge
    Point_3D M_disp = V * gnp_i.l * 0.5;

    //Calculate vertices along reference edge
    gnp_i.vertices[R1] = M + M_disp;
    gnp_i.vertices[RT] = M - M_disp;
    //hout << "R1=" << R1 << " RT=" << RT << endl;

    //Calculate displacement though the thickness
    Point_3D t_disp = N_top * gnp_i.t;

    //Calculate vertices below reference edge
    gnp_i.vertices[R1 + 4] = gnp_i.vertices[R1] - t_disp;
    gnp_i.vertices[RT + 4] = gnp_i.vertices[RT] - t_disp;
    //hout << "R1+4=" << R1+4 << " RT+4=" << RT + 4 << endl;

    //Calculate displacement towards opposite thin face
    //No need to make unit vector, since vector in cross produc are unitary
    Point_3D N_in = V.cross(N_top);

    //Displacement along the length
    Point_3D d_disp = N_in * gnp_i.l;

    //Calculate the top vertices on opposite thin face
    int idx1 = (RT + 1) % 4;
    gnp_i.vertices[idx1] = gnp_i.vertices[RT] + d_disp;
    int idx2 = (RT + 2) % 4;
    gnp_i.vertices[idx2] = gnp_i.vertices[R1] + d_disp;
    //hout << "idx1=" << idx1 << " idx2=" << idx2 << endl;

    //Calculate the bottom vertices on opposite thin face
    gnp_i.vertices[idx1 + 4] = gnp_i.vertices[idx1] - t_disp;
    gnp_i.vertices[idx2 + 4] = gnp_i.vertices[idx2] - t_disp;
    //hout << "idx1+4=" << idx1 + 4 << " idx2+4=" << idx2 + 4 << endl;

    return 1;
}
//This function reconstructs a GNP partially inside the sample for the case when
//two consecutive short edges are inside the sample
int GNP_Reconstruction::Two_consecutive_short_edges(const vector<bool>& vertex_flags, GNP& gnp_i, Point_3D& N_top)const
{
    //Initial gnp
    //GNP gnp0 = gnp_i;

    //Variables to store the reference plane, i.e., these are the reference vertices
    int R1, R2, R3, R4;
    
    //Variables to store other vertices if present
    //Initialized to -1 to indicate they are not present
    int O1 = -1, O2 = -1;
    
    //Get the actual values of the reference vertices
    if (!Get_reference_vertices_case2(vertex_flags, R1, R2, R3, R4, O1, O2))
    {
        hout<<"Error in Two_consecutive_short_edges when calling Get_reference_vertices_case2"<<endl;
        return 0;
    }
    //hout << "R1=" << R1 << " R2=" << R2 << " R3=" << R3 << " R4=" << R4 << " O1=" << O1 << " O2=" << O2 << endl;

    //Get all vertices R1 to R4 on a reference plane
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

    /* /Plane with final GNP
    Plane_3D Pl(gnp_i.vertices[R1], gnp_i.vertices[R2], gnp_i.vertices[R3]);
    //Distance from initial GNP to final GNP face
    for (size_t i = 0; i < vertex_flags.size(); i++)
    {
        if (vertex_flags[i])
        {
            hout << "d" << i << "=" << Pl.distance_to(gnp0.vertices[i]) << endl;
        }
    }*/
    
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
    
    //Get signed distance of R3 to plane R1R2R4
    double d_sign = (gnp_i.vertices[R3] - gnp_i.vertices[R4]).dot(N_in);
    //hout << "d_sign=" << d_sign << endl;
    
    //Check if vertex R3 is outside the GNP
    if (d_sign < 0.0)
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
    
    //Distance from M to R1 and R2, which is the same
    //since M is the midpoint of R1R2
    double d1 = M.distance_to(gnp_i.vertices[R2]);
    
    //Get distance from M to R3 and R4 along U
    double d3 = U.dot(gnp_i.vertices[R3] - M);
    double d4 = U.dot(M - gnp_i.vertices[R4]);
    
    //Find maximum distance and index
    double d_max = d1;
    int idx = 1;
    
    if (d3 - d_max > Zero)
    {
        //Update maximum distance and index
        d_max = d3;
        idx = 3;
    }
    if (d4 - d_max > Zero)
    {
        //Update maximum distance and index
        d_max = d4;
        idx = 4;
    }
    
    //Variables for the distances to other vertices, if any
    double do1 = 0.0, do2 = 0.0;
    
    //Check if there are other vertices that might affect the calculation of
    //the GNP thickness
    if (O1 != -1)
    {
        //Calculate the distance to O1 along U
        do1 = abs(U.dot(M - gnp_i.vertices[O1]));
        //hout << "d_max=" << d_max << " do1=" << do1 << endl;
        
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
    //hout << "d1=" << d1 << " d3=" << d3 << " d4=" << d4 << " do1=" << do1 << " do2=" << do2 << " d_max=" << d_max << endl;
    
    //Now that the maximum distance has been determined, move the vertices
    if (idx != 1)
    {
        //Adjust R1 (and R2) since maximum is not at R1 (nor R2)
        gnp_i.vertices[R1] = M - U*d_max;
        gnp_i.vertices[R2] = M + U*d_max;
    }
    if (idx != 3)
    {
        //Adjust R3 since maximum is not at R3
        gnp_i.vertices[R3] = gnp_i.vertices[R3] + U*(d_max - d3);
    }
    if (idx != 4)
    {
        //Adjust R4 since maximum is not at R4
        gnp_i.vertices[R4] = gnp_i.vertices[R4] - U*(d_max - d4);
    }
    //hout << "d12=" << gnp_i.vertices[R1].distance_to(gnp_i.vertices[R2]);
    //hout << " d34=" << gnp_i.vertices[R3].distance_to(gnp_i.vertices[R4]) << endl;
    
    //Adjust O1 and O2 if present and they were not the ones with maximum distance
    if (O1 != -1 && idx != -1) {
        //Get the displacement
        Point_3D S = (O1 <= 3)? U*(do1 - d_max): U*(d_max - do1);
        gnp_i.vertices[O1] = gnp_i.vertices[O1] + S;
    }
    if (O2 != -1 && idx != -2) {
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
    //which is the also the distance from M to R4
    double d1 = M.distance_to(gnp_i.vertices[R1]);
    
    //Initialize maximum distance 1 and index of maximum distance
    double d_max = d1;
    
    //Calcualte distances to R2 and R3
    double d2 = V.dot(M - gnp_i.vertices[R2]);
    double d3 = V.dot(gnp_i.vertices[R3] - M);
    
    //Check if maximum distance needs updating
    if (d2 - d_max > Zero)
    {
        //Update maximum distance and index
        d_max = d2;
    }
    if (d3 - d_max > Zero)
    {
        //Update maximum distance and index
        d_max = d3;
    }
    
    //Variables for the distances to other vertices, if any
    double do1 = 0.0, do2 = 0.0;

    //Variable for second approximation of l_GNP
    double l2 = 0.0;

    //Plane with reference face
    Plane_3D Pl;
    if (O1 != -1 || O2 != -1)
    {
        Pl = Plane_3D(gnp_i.vertices[R1], gnp_i.vertices[R2], gnp_i.vertices[R3]);
    }
    
    //Check if there are other vertices that might affect the calculation of
    //the GNP thickness
    if (O1 != -1)
    {
        //Calculate the distance to O1 along V
        do1 = abs(V.dot(M - gnp_i.vertices[O1]));

        //Calculate the distance from the reference face to O1
        l2 = Pl.distance_to(gnp_i.vertices[O1]);
        
        //Check if maximum distance needs to be updated
        if (do1 - d_max > Zero)
        {
            d_max = do1;
        }
    }
    if (O2 != -1)
    {
        //Calculate the distance to O2 along V
        do2 = abs(V.dot(M - gnp_i.vertices[O2]));

        //Calculate the distance from the reference face to O2
        double l3 = Pl.distance_to(gnp_i.vertices[O2]);
        if (l3 - l2 > Zero)
        {
            l2 = l3;
        }
        
        //Check if maximum distance needs to be updated
        if (do2 - d_max > Zero)
        {
            d_max = do2;
        }
    }
    //hout << "d14=" << gnp_i.vertices[R1].distance_to(gnp_i.vertices[R4]) << " d23=" << gnp_i.vertices[R3].distance_to(gnp_i.vertices[R2]) << endl;
    //hout << "d1=" << d1 << " d2=" << d2 << " d3=" << d3 << " do1=" << do1 << " do2=" << do2 << endl;
    //hout << "d_max=" << d_max << " l2=" << l2 << endl;
    
    //Now that the maximum distance has been determined, set the GNP side length
    gnp_i.l = max(2.0*d_max, l2);
    //hout << "l_GNP=" << gnp_i.l << endl;
    
    //Calculate the eight GNP vertices
    
    //Update R1 and R4
    gnp_i.vertices[R1] = M - V*d_max;
    //hout << "updated " << R1 << endl;
    gnp_i.vertices[R4] = M + V*d_max;
    //hout << "updated " << R4 << endl;
    
    //Get a vector going inside the GNP through the refrence face R1R2R3R4
    Point_3D N_in = V.cross(gnp_i.vertices[R3] - gnp_i.vertices[R1]);
    N_in.make_unit();
    
    //Displacement from reference face
    Point_3D disp_face = N_in*gnp_i.l;
    
    //Calculate remaining two vertices on top face
    int idx_top1 = (R4 + 1) % 4;
    gnp_i.vertices[idx_top1] = gnp_i.vertices[R4] + disp_face;
    //hout << "updated " << idx_top1 << endl;
    int idx_top2 = (R4 + 2) % 4;
    gnp_i.vertices[idx_top2] = gnp_i.vertices[R1] + disp_face;
    //hout << "updated " << idx_top2 << endl;
    
    //Calculate vector normal to top face
    //Both V and N_in are normal and unitary, thus its dot product is a unitary vector
    //Thus there is no need to normalize N_top
    N_top = V.cross(N_in);
    
    //Calculate displacemente through the thickness
    Point_3D disp_thick = N_top*gnp_i.t;
    
    //Calculate vertices in bottom face
    gnp_i.vertices[R2] = gnp_i.vertices[R1] - disp_thick;
    //hout << "updated " << R2 << endl;
    gnp_i.vertices[R3] = gnp_i.vertices[R4] - disp_thick;
    //hout << "updated " << R3 << endl;
    //hout << "d14=" << gnp_i.vertices[R1].distance_to(gnp_i.vertices[R4]) << " d23=" << gnp_i.vertices[R3].distance_to(gnp_i.vertices[R2]) << endl;
    int idx_bot = idx_top1 + 4;
    gnp_i.vertices[idx_bot] = gnp_i.vertices[idx_top1] - disp_thick;
    //hout << "updated " << idx_bot << endl;
    idx_bot = idx_top2 + 4;
    gnp_i.vertices[idx_bot] = gnp_i.vertices[idx_top2] - disp_thick;
    //hout << "updated " << idx_bot << endl;
    
    return 1;
}
//This function reconstructs a GNP partially inside the sample for the case when
//two non consecutive short edges are inside the sample
int GNP_Reconstruction::Two_non_consecutive_short_edges(const vector<bool>& vertex_flags, GNP& gnp_i, Point_3D& N_top)const
{
    //Variables to store the reference vertices
    int R1, R2, R3, R4;
    
    //Variables to store other vertices, if any
    int Ou = -1, Ov = -1;
    
    //Get the reference vertices
    if (!Get_reference_vertices_case4(vertex_flags, R1, R2, R3, R4, Ou, Ov))
    {
        hout<<"Error in Two_non_consecutive_short_edges when calling Get_reference_vertices_case4"<<endl;
        return 0;
    }
    
    //Get the unit vector along R1R2, which will serve as normal of the top surface
    //So in this case we actually have a unitary vector R2R1
    N_top = (gnp_i.vertices[R1] - gnp_i.vertices[R2]).unit();
    
    //Find the maximum distance from M to any vertex along N_top in order to
    //set the vertices on two parallel planes that correspond to the square faces
    //In the process the GNP thickness is calculated
    if (!Set_parallel_planes_for_square_faces_case4(R1, R2, R3, R4, Ou, Ov, gnp_i, N_top))
    {
        hout<<"Error in Two_non_consecutive_short_edges when calling Set_parallel_planes_for_square_faces_case4"<<endl;
        return 0;
    }
    
    //Get unit vector U and V, which go along the top squared face of the GNP
    Point_3D U, V;
    if (!Get_unit_vectors_on_square_surface_case4(R1, R4, gnp_i, N_top, U, V))
    {
        hout<<"Error in Two_non_consecutive_short_edges when calling Get_unit_vectors_on_square_surface"<<endl;
        return 0;
    }
    
    //Find the GNP side length
    if (!Calculate_gnp_side_length_case4(R1, R2, R3, R4, Ou, Ov, U, V, gnp_i))
    {
        hout<<"Error in Two_non_consecutive_short_edges when calling Calculate_gnp_side_length_case4"<<endl;
        return 0;
    }
    
    //Recalculate the GNP vertices
    if (!Calculate_gnp_vertices_case4(R1, R2, U, V, gnp_i))
    {
        hout<<"Error in Two_non_consecutive_short_edges when calling Calculate_gnp_vertices_case4"<<endl;
        return 0;
    }
    
    return 1;
}
//This function finds the reference vertices for the case of
//two non consecutive short edges inside the sample
int GNP_Reconstruction::Get_reference_vertices_case4(const vector<bool>& vertex_flags, int& R1, int& R2, int& R3, int& R4, int& Ou, int& Ov)const
{
    //There are only two options, check which of the them is
    if (vertex_flags[0] && vertex_flags[4] && vertex_flags[2] && vertex_flags[6])
    {
        //Set reference vertices
        R1 = 0; R2 = 4; R3 = 6; R4 = 2;
        
        //Set other vertices if present
        Ou = (vertex_flags[3])? 3 : (vertex_flags[7])? 7 : -1;
        Ov = (vertex_flags[1])? 1 : (vertex_flags[5])? 5 : -1;
    }
    else if (vertex_flags[1] && vertex_flags[5] && vertex_flags[3] && vertex_flags[7])
    {
        //Set reference vertices
        R1 = 1; R2 = 5; R3 = 7; R4 = 3;
        
        //Set other vertices if present
        Ou = (vertex_flags[0])? 0 : (vertex_flags[4])? 4 : -1;
        Ov = (vertex_flags[2])? 2 : (vertex_flags[6])? 6 : -1;
    }
    else
    {
        hout<<"Error in Get_reference_vertices_case4. There are no two edges across the diagonal inside the sample (case 4)."<<endl;
        hout << "Boolean vector: " << endl;
        for (size_t i = 0; i < vertex_flags.size(); i++)
            hout << "vertex_flags[" << i << "]=" << vertex_flags[i] << endl;
        return 0;
    }
    
    return 1;
}
//This function sets the vertices on two parallel planes that correspond to the square faces
//for the case of two non consecutive short edges are inside the sample
int GNP_Reconstruction::Set_parallel_planes_for_square_faces_case4(const int& R1, const int& R2, const int& R3, const int& R4, const int& Ou, const int& Ov, GNP& gnp_i, Point_3D& N_top)const
{
    //Get the midpoint of R1R2
    Point_3D M = (gnp_i.vertices[R1] + gnp_i.vertices[R2])*0.5;
    
    //Calculate the distance from M to R1 and R2
    double d12 = M.distance_to(gnp_i.vertices[R1]);
    
    //Initialize the maximum distance with the calculated distance
    double d_max = d12;
    //Initialize index of maximum distance
    double idx = 1;
    
    //Calculate the distance to the remaining reference vertices
    //At the same time look for the maximum distance
    double d3 = N_top.dot(M - gnp_i.vertices[R3]);
    if (d3 - d_max > Zero) {
        d_max = d3;
        idx = 3;
    }
    double d4 = N_top.dot(gnp_i.vertices[R4] - M);
    if (d4 - d_max > Zero) {
        d_max = d4;
        idx = 4;
    }
    
    //Calculate the distances to the other vertices if present
    //At the same time look for the maximum distance
    double dou = (Ou != -1)? abs(N_top.dot(gnp_i.vertices[Ou] - M)) : 0.0;
    if (dou - d_max > Zero) {
        d_max = dou;
        idx = -1;
    }
    double dov = (Ov != -1)? abs(N_top.dot(gnp_i.vertices[Ov] - M)) : 0.0;
    if (dov - d_max > Zero) {
        d_max = dov;
        idx = -2;
    }
    
    //At this point the maximum distance was found
    //Update the GNP thickness
    gnp_i.t = 2.0*d_max;
    
    //Displacement along the thickness
    Point_3D disp_t = N_top*d_max;
    
    //Move all vertices to the maximum distance from M
    if (idx != 1)
    {
        //Vertices R1 and R2 are not at the maximum distance, so update them
        gnp_i.vertices[R1] = M + disp_t;
        gnp_i.vertices[R2] = M - disp_t;
    }
    if (idx != 3)
    {
        //Vertex R3 is not at the maximum distance, so update it
        gnp_i.vertices[R3] = gnp_i.vertices[R3] - N_top*(d_max - d3);
    }
    if (idx != 4)
    {
        //Vertex R3 is not at the maximum distance, so update it
        gnp_i.vertices[R4] = gnp_i.vertices[R4] + N_top*(d_max - d4);
    }
    if (idx != -1)
    {
        //Vertex Ou is not at the maximum distance, so update it
        double mult = (Ou != -1 && Ou <= 3)? (d_max - dou) : (dou - d_max) ;
        gnp_i.vertices[Ou] = gnp_i.vertices[Ou] - N_top*mult;
    }
    if (idx != -2)
    {
        //Vertex Ov is not at the maximum distance, so update it
        double mult = (Ov != -1 && Ov <= 3)? (d_max - dov) : (dov - d_max) ;
        gnp_i.vertices[Ov] = gnp_i.vertices[Ov] - N_top*mult;
    }
    
    return 1;
}
//This function calculates the vectors that go along the top squared face of the GNP
int GNP_Reconstruction::Get_unit_vectors_on_square_surface_case4(const int& R1, const int& R4, const GNP& gnp_i, const Point_3D& N_top, Point_3D& U, Point_3D& V)const
{
    //Get unit vector R1R4
    Point_3D R1R4 = (gnp_i.vertices[R4] - gnp_i.vertices[R1]).unit();
    
    //Get vector notmal to both N_top and R1R4
    //This vector is unitary since both N_top and R1R4 are also unitary
    Point_3D U_45off = N_top.cross(R1R4);
    
    //Vector U comes from normalizing the addition of vectors R1R4 and U_45off
    U = (R1R4 + U_45off).unit();
    
    //Vector V is normal to N_top and U
    V = U.cross(N_top);
    
    return 1;
}
//This functions determines the GNP side length
int GNP_Reconstruction::Calculate_gnp_side_length_case4(const int& R1, const int& R2, const int& R3, const int& R4, const int& Ou, const int& Ov, const Point_3D& U, const Point_3D& V, GNP& gnp_i)const
{
    //Get distances from R1 along U while saving the maximum distance
    
    //Initialize maximum distance with distance to R3 along U
    double d_max = U.dot(gnp_i.vertices[R3] - gnp_i.vertices[R1]);
    
    //Calculate distance to R4 along U and update maximum if needed
    double d_new = U.dot(gnp_i.vertices[R4] - gnp_i.vertices[R1]);
    if (d_new - d_max > Zero)
    {
        d_max = d_new;
    }
    
    //Calculate distance to Ou along U if present and update maximum if needed
    d_new = (Ou != -1)? U.dot(gnp_i.vertices[Ou] - gnp_i.vertices[R1]) : 0.0;
    if (d_new - d_max > Zero)
    {
        d_max = d_new;
    }
    
    //Calculate distance to R3 along V and update maximum if needed
    d_new = V.dot(gnp_i.vertices[R3] - gnp_i.vertices[R1]);
    if (d_new - d_max > Zero)
    {
        d_max = d_new;
    }
    
    //Calculate distance to R4 along V and update maximum if needed
    d_new = V.dot(gnp_i.vertices[R4] - gnp_i.vertices[R1]);
    if (d_new - d_max > Zero)
    {
        d_max = d_new;
    }
    
    //Calculate distance to Ov along V if present and update maximum if needed
    d_new = (Ov != -1)? V.dot(gnp_i.vertices[Ov] - gnp_i.vertices[R1]) : 0.0;
    if (d_new - d_max > Zero)
    {
        d_max = d_new;
    }
    
    //Update GNP side length
    gnp_i.l = d_max;
    
    return 1;
}
//This function calculates the vertices of the GNP for the case when
//two non consecutive short edges are inside the sample
int GNP_Reconstruction::Calculate_gnp_vertices_case4(const int& R1, const int& R2, const Point_3D& U, const Point_3D& V, GNP& gnp_i)const
{
    //Reconstruction starts from edge R1R2
    
    //Displacement along V
    Point_3D V_disp = V*gnp_i.l;
    
    //Calculate vertices in edge to the right of edge R1R2
    int v = (R1 + 1) % 4;
    gnp_i.vertices[v] = gnp_i.vertices[R1] + V_disp;
    gnp_i.vertices[v + 4] = gnp_i.vertices[R2] + V_disp;
    
    //Displacement along U
    Point_3D U_disp = U*gnp_i.l;
    
    //Calculate vertices in edge to the left of edge R1R2
    v = (R1 + 3) % 4;
    gnp_i.vertices[v] = gnp_i.vertices[R1] + U_disp;
    gnp_i.vertices[v + 4] = gnp_i.vertices[R2] + U_disp;
    
    //Calculate vertices in the edge across edge R1R2
    int a = (R1 + 2) % 4;
    gnp_i.vertices[a] = gnp_i.vertices[v] + V_disp;
    gnp_i.vertices[a + 4] = gnp_i.vertices[v + 4] + V_disp;
    
    return 1;
}
//This function reconstructs a GNP partially inside the sample for the case when
//one short edge is inside the sample
int GNP_Reconstruction::One_short_edge(const vector<bool>& vertex_flags, GNP& gnp_i, Point_3D& N_top)const
{
    //Variable to store the vertex that will be used to form a second short edge
    int V;

    //Variable to determine the case in which case 1 eill be transformed
    int new_case = -1;
    
    //Find a vertex that converts this case into case 2 (two consecutive short edges)
    if (!Find_vertex_for_two_short_edges_case1(vertex_flags, gnp_i, V, new_case))
    {
        hout<<"Error in One_short_edge when calling Find_vertex_for_two_short_edges_case1"<<endl;
        return 0;
    }
    
    //Create a new vertex_flags vector with V set to true
    vector<bool> new_vertex_flags = vertex_flags;
    //hout << "new_vertex_flags=" << new_vertex_flags.size()<<" V="<<V << endl;
    new_vertex_flags[V] = true;
    
    //Go to the corresponding case
    if (new_case == 2)
    {
        //Call case 2 with new vertex_flags vector
        //hout << "Two_consecutive_short_edges" << endl;
        if (!Two_consecutive_short_edges(new_vertex_flags, gnp_i, N_top))
        {
            hout << "Error in One_short_edge when calling Two_consecutive_short_edges" << endl;
            return 0;
        }
    }
    else if (new_case == 4)
    {
        //Call case 4 with new vertex_flags vector
        //hout << "Two_non_consecutive_short_edges" << endl;
        if (!Two_non_consecutive_short_edges(new_vertex_flags, gnp_i, N_top))
        {
            hout << "Error in One_short_edge when calling Two_non_consecutive_short_edges" << endl;
            return 0;
        }
    }
    else
    {
        hout << "Error in One_short_edge. Case 2 or 4 could not be completed. new_case=" << new_case << endl;
        return 0;
    }
    
    return 1;
}
//This function finds a vertex that completes a second consecutive short edge
int GNP_Reconstruction::Find_vertex_for_two_short_edges_case1(const vector<bool>& vertex_flags, GNP& gnp_i, int& V, int& new_case)const
{
    //Find the short edge inside the sample
    if (vertex_flags[0] && vertex_flags[4])
    {
        //Find adjacent vertex
        if (!Find_adjacent_vertex_case1(vertex_flags, 0, 4, V, new_case, gnp_i))
        {
            hout<<"Error in Find_vertex_for_two_short_edges_case1 when calling Find_adjacent_vertex_case1 on edge 04."<<endl;
            return 0;
        }
    }
    else if (vertex_flags[1] && vertex_flags[5])
    {
        //Find adjacent vertex
        if (!Find_adjacent_vertex_case1(vertex_flags, 1, 5, V, new_case, gnp_i))
        {
            hout<<"Error in Find_vertex_for_two_short_edges_case1 when calling Find_adjacent_vertex_case1 on edge 15."<<endl;
            return 0;
        }
    }
    else if (vertex_flags[2] && vertex_flags[6])
    {
        //Find adjacent vertex
        if (!Find_adjacent_vertex_case1(vertex_flags, 2, 6, V, new_case, gnp_i))
        {
            hout<<"Error in Find_vertex_for_two_short_edges_case1 when calling Find_adjacent_vertex_case1 on edge 26."<<endl;
            return 0;
        }
    }
    else if (vertex_flags[3] && vertex_flags[7])
    {
        //Find adjacent vertex
        if (!Find_adjacent_vertex_case1(vertex_flags, 3, 7, V, new_case, gnp_i))
        {
            hout<<"Error in Find_vertex_for_two_short_edges_case1 when calling Find_adjacent_vertex_case1 on edge 37."<<endl;
            return 0;
        }
    }
    else
    {
        hout<<"Error in Find_vertex_for_two_short_edges_case1. No short edge inside the sample was found although there should be exactly one."<<endl;
        return 0;
    }
    
    return 1;
}
//Given a short edge, this function finds an adjacent edge that completes a second short edge
int GNP_Reconstruction::Find_adjacent_vertex_case1(const vector<bool>& vertex_flags, const int& R1, const int& R2, int& V, int& new_case, GNP& gnp_i)const
{
    //R1 is a top edge and R2 is a bottom edge
    
    //Get vertex to the right of R1
    int R = (R1 + 1) % 4;
    //Get vertex to the left of R1
    int L = (R1 + 3) % 4;
    //Get vertex across R1
    int A = (R1 + 2) % 4;
    //hout << "R=" << R << " L=" << L << endl;
    
    //Get displacement vector going from bottom square face towards top square face
    Point_3D disp_t = gnp_i.vertices[R1] - gnp_i.vertices[R2];
    
    //Check vertex to the right of R1
    if (vertex_flags[R])
    {
        //Set V as the vertex below R
        V = R + 4;
        //hout << "V = R + 4=" << V << endl;
        
        //Check that V is indeed not present, otherwise there is an error
        if (vertex_flags[V]) {
            hout<<"Error in Find_adjacent_vertex_case1. There is a second consecutive short edge ('to the right') inside the sample. Initial short edge is "<<R1<<"-"<<R2<<". Second consecutive edge is "<<R<<"-"<<V<<". The case identified was the one with only one short edge inside the sample."<<endl;
            return 0;
        }
        
        //Calculate the coordinates of vertex V, which is just going down from vertex R
        gnp_i.vertices[V] = gnp_i.vertices[R] - disp_t;
    }
    //Check vertex to the right of R2, i.e., the one below R
    else if (vertex_flags[R + 4])
    {
        //Set V as R
        V = R;
        //hout << "V = R =" << V << endl;
        
        //There is no need to check if R is present, as this was already done
        //resulting in false and thus reaching this estatement
        
        //Calculate the coordinates of vertex V, which is just going up from vertex R+4
        gnp_i.vertices[V] = gnp_i.vertices[R + 4] + disp_t;
    }
    //Check vertex to the left of R1
    else if (vertex_flags[L])
    {
        //Set V as the vertex below L
        V = L + 4;
        //hout << "V = L + 4=" << V << endl;
        
        //Check that V is indeed not present, otherwise there is an error
        if (vertex_flags[V]) {
            hout<<"Error in Find_adjacent_vertex_case1. There is a second consecutive short edge ('to the left') inside the sample. Initial short edge is "<<R1<<"-"<<R2<<". Second consecutive edge is "<<L<<"-"<<V<<". The case identified was the one with only one short edge inside the sample."<<endl;
            return 0;
        }
        
        //Calculate the coordinates of vertex V, which is just going down from vertex L
        gnp_i.vertices[V] = gnp_i.vertices[L] - disp_t;
    }
    //Check vertex to the left of R2, i.e., the one below L
    else if (vertex_flags[L + 4])
    {
        //Set V as L
        V = L;
        //hout << "V = L =" << V << endl;
        
        //There is no need to check if L is present, as this was already done
        //resulting in false and thus reaching this estatement
        
        //Calculate the coordinates of vertex V, which is just going up from vertex L+4
        gnp_i.vertices[V] = gnp_i.vertices[L + 4] + disp_t;
    }
    //Check vertex across R1
    else if (vertex_flags[A])
    {
        //Set V as A + 4
        V = A + 4;
        //hout << "V = A + 4 =" << V << endl;
        
        //Check that V is indeed not present, otherwise there is an error
        if (vertex_flags[V]) {
            hout << "Error in Find_adjacent_vertex_case1. There is a second non-consecutive short edge inside the sample. Initial short edge is " << R1 << "-" << R2 << ". Second consecutive edge is " << A << "-" << V << ". The case identified was the one with only one short edge inside the sample." << endl;
            return 0;
        }

        //Calculate the coordinates of vertex V, which is just going down from vertex A
        gnp_i.vertices[V] = gnp_i.vertices[A] - disp_t;

        //New case is 4
        new_case = 4;
    }
    //Check vertex across R2
    else if (vertex_flags[A + 4])
    {
        //Set V as A
        V = A;
        //hout << "V = A =" << V << endl;
        
        //There is no need to check if A is present, as this was already done
        //resulting in false and thus reaching this estatement
        
        //Calculate the coordinates of vertex V, which is just going up from vertex A+4
        gnp_i.vertices[V] = gnp_i.vertices[A + 4] + disp_t;

        //New case is 4
        new_case = 4;
    }
    else
    {
        hout<<"Error in Find_adjacent_vertex_case1. No adjacent vertex to edge "<<R1<<"-"<<R2<<" was found. GNP with one short edge inside the sample cannot be reconstructed."<<endl;
        hout << "Boolean vector: " << endl;
        for (size_t i = 0; i < vertex_flags.size(); i++)
            hout << "vertex_flags[" << i << "]=" << vertex_flags[i] << endl;
        return 0;
    }

    //If this part of the code is reached, there no errors were found
    //So to avoid repeating code, check if case 2 is the new case
    if (new_case != 4)
    {
        new_case = 2;
    }
    
    return 1;
}
//This function reconstructs a GNP partially inside the sample for the case when
//there are no shorts edges inside the sample
int GNP_Reconstruction::No_short_edges(const vector<bool>& vertex_flags, GNP& gnp_i, Point_3D& N_top)const
{
    //Variables to store reference vertices
    //For this case R4 may not be present
    int R1, R2, R3, R4;

    //Get the reference vertices
    if (!Get_reference_vertices_case0(vertex_flags, R1, R2, R3, R4))
    {
        hout << "Error in No_short_edges when calling Get_reference_vertices_case0" << endl;
        return 0;
    }

    //Unit vectors needed to calculate GNP vertices in their final position
    Point_3D U, V;

    //Variables to approximate the GNP side length
    double l1, l2;

    //Set edges parallel along vector R1R2
    if (!Set_parallel_edges_along_u_case0(R1, R2, R3, R4, gnp_i, l1))
    {
        hout << "Error in No_short_edges when calling Set_parallel_edges_along_u_case0" << endl;
        return 0;
    }

    //Set edges parallel along vector R2R3
    if (!Set_parallel_edges_along_v_case0(R1, R2, R3, R4, gnp_i, l2, V))
    {
        hout << "Error in No_short_edges when calling Set_parallel_edges_along_v_case0" << endl;
        return 0;
    }

    //Get the GNP side length
    gnp_i.l = 2.0 * max(l1, l2);

    //Calculate GNP vertices
    if (!Calculate_gnp_vertices_case0(R1, R2, R3, R4, V, gnp_i, N_top))
    {
        hout << "Error in No_short_edges when calling Calculate_gnp_vertices_case0" << endl;
        return 0;
    }

    return 1;
}
//This function finde the reference vertices for the case when
//there are no shorts edges inside the sample
int GNP_Reconstruction::Get_reference_vertices_case0(const vector<bool>& vertex_flags, int& R1, int& R2, int& R3, int& R4)const
{
    //Find the reference vertices by scanning pairs across the diagonal
    //First check vertices in top face
    if (vertex_flags[0] && vertex_flags[2])
    {
        //Check for third vertex
        if (vertex_flags[1])
        {
            //Set the three reference vertices that for sure are present
            R1 = 0; R2 = 1; R3 = 2;
            //Check if the fouth vertex is present
            R4 = (vertex_flags[3]) ? 3 : -1;
        }
        else if (vertex_flags[3])
        {
            //Set the three reference vertices that for sure are present
            R1 = 2; R2 = 3; R3 = 0;
            //The fourth vertex is not present since the previous if-statement was false
            //in orther to reach this part of the code
            R4 = -1;
        }
        else
        {
            hout << "Error in Get_reference_vertices_case0. There are not enough vertices in top square face to reconsturct GNP (1)." << endl;
            return 0;
        }
    }
    else if (vertex_flags[1] && vertex_flags[3])
    {
        //Check for third vertex
        if (vertex_flags[2])
        {
            //Set the three reference vertices that for sure are present
            R1 = 0; R2 = 1; R3 = 2;
        }
        else if (vertex_flags[0])
        {
            //Set the three reference vertices that for sure are present
            R1 = 3; R2 = 0; R3 = 1;
        }
        else
        {
            hout << "Error in Get_reference_vertices_case0. There are not enough vertices in top square face to reconsturct GNP (2)." << endl;
            return 0;
        }

        //The fourth vertex is not present since the previous if-statement was false
        //in orther to reach this part of the code. For the previous if-statement
        //to be false, at least one vertex was not present, thus at most three are present
        //Then for sure R4 is not present
        R4 = -1;
    }
    //Now check vertices in bottom face
    else if (vertex_flags[4] && vertex_flags[6])
    {
        //Check for third vertex
        if (vertex_flags[5])
        {
            //Set the three reference vertices that for sure are present
            R1 = 4; R2 = 5; R3 = 6;
            //Check if the fouth vertex is present
            R4 = (vertex_flags[7]) ? 7 : -1;
        }
        else if (vertex_flags[7])
        {
            //Set the three reference vertices that for sure are present
            R1 = 6; R2 = 7; R3 = 4;
            //The fourth vertex is not present since the previous if-statement was false
            //in orther to reach this part of the code
            R4 = -1;
        }
        else
        {
            hout << "Error in Get_reference_vertices_case0. There are not enough vertices in bottom square face to reconsturct GNP (1)." << endl;
            return 0;
        }
    }
    else if (vertex_flags[5] && vertex_flags[7])
    {
        //Check for third vertex
        if (vertex_flags[4])
        {
            //Set the three reference vertices that for sure are present
            R1 = 7; R2 = 4; R3 = 5;
        }
        else if (vertex_flags[6])
        {
            //Set the three reference vertices that for sure are present
            R1 = 5; R2 = 6; R3 = 7;
        }
        else
        {
            hout << "Error in Get_reference_vertices_case0. There are not enough vertices in bottom square face to reconsturct GNP (2)." << endl;
            return 0;
        }

        //The fourth vertex is not present since the previous if-statement was false
        //in orther to reach this part of the code. For the previous if-statement
        //to be false, at least one vertex was not present, thus at most three are present
        //Then for sure R4 is not present
        R4 = -1;
    }
    else
    {
        hout << "Error in Get_reference_vertices_case0. There are not enough vertices in GNP to reconsturct it." << endl;
        return 0;
    }

    return 1;
}
//This function sets edges R1R4 and R2R3 as parallel
int GNP_Reconstruction::Set_parallel_edges_along_u_case0(const int& R1, const int& R2, const int& R3, const int& R4, GNP& gnp_i, double& l1)const
{
    //Edge R1R2 will be a reference edge
    //Calculate midpoint of reference edge
    Point_3D M = (gnp_i.vertices[R1] + gnp_i.vertices[R2])*0.5;

    //Calculate vector U, which is the unit vectro along R1R2
    Point_3D U = (gnp_i.vertices[R2] - gnp_i.vertices[R1]).unit();

    //Find maximum distance from M to all vertices
    //Start with distance to R1 and R2
    double d_max = M.distance_to(gnp_i.vertices[R1]);
    //Index to determine where the maxium distance was found
    int idx = 1;

    //Find distance to R3 and update if needed
    double d3 = U.dot(gnp_i.vertices[R3] - M);
    if (d3 - d_max > Zero)
    {
        //Update distance and index
        d_max = d3;
        idx = 3;
    }

    //Find distance to R4 and update if needed
    double d4 = (R4 != -1)? abs(U.dot(gnp_i.vertices[R4] - M)) : 0.0;
    if (d4 - d_max > Zero)
    {
        //Update distance and index
        d_max = d4;
        idx = 4;
    }

    //Set l1 as the maximum distance found
    l1 = d_max;

    //Move vertices so that edges R1R4 and R2R3 are parallel
    if (idx != 1)
    {
        //Update R1 and R2 since they are not at maximum distance
        gnp_i.vertices[R1] = M - U * d_max;
        gnp_i.vertices[R2] = M + U * d_max;
    }
    if (idx != 3)
    {
        //Update R3 since it is not at maximum distance
        gnp_i.vertices[R3] = gnp_i.vertices[R3] + U * (d_max - d3);
    }
    if (R4 != -1 && idx != 4)
    {
        //Update R4 since it is not at maximum distance
        gnp_i.vertices[R4] = gnp_i.vertices[R4] - U * (d_max - d4);
    }

    return 1;
}
//This function sets edges R1R2 and R4R3 as parallel
int GNP_Reconstruction::Set_parallel_edges_along_v_case0(const int& R1, const int& R2, const int& R3, const int& R4, GNP& gnp_i, double& l2, Point_3D& V)const
{
    //Edge R2R3 will be a reference edge
    //Calculate midpoint of reference edge
    Point_3D M = (gnp_i.vertices[R2] + gnp_i.vertices[R3]) * 0.5;

    //Calculate vector V, which is the unit vectro along R2R3
    V = (gnp_i.vertices[R3] - gnp_i.vertices[R2]).unit();

    //Find maximum distance from M to all vertices
    //Start with distance to R2 and R3
    double d_max = M.distance_to(gnp_i.vertices[R2]);
    //Index to determine where the maxium distance was found
    int idx = 2;

    //Find distance to R1 and update if needed
    double d1 = V.dot(M - gnp_i.vertices[R1]);
    if (d1 - d_max > Zero)
    {
        //Update distance and index
        d_max = d1;
        idx = 1;
    }

    //Find distance to R4 and update if needed
    double d4 = (R4 != -1) ? abs(V.dot(gnp_i.vertices[R4] - M)) : 0.0;
    if (d4 - d_max > Zero)
    {
        //Update distance and index
        d_max = d4;
        idx = 4;
    }

    //Set l2 as the maximum distance found
    l2 = d_max;

    //Move vertices so that edges R1R2 and R4R3 are parallel
    if (idx != 1)
    {
        //Update R1 since it is not at maximum distance
        gnp_i.vertices[R1] = gnp_i.vertices[R1] - V * (d_max - d1);
    }
    if (idx != 2)
    {
        //Update R2 and R3 since they are not at maximum distance
        gnp_i.vertices[R2] = M - V * d_max;
        gnp_i.vertices[R3] = M + V * d_max;
    }
    if (R4 != -1 && idx != 4)
    {
        //Update R4 since it is not at maximum distance
        gnp_i.vertices[R4] = gnp_i.vertices[R4] + V * (d_max - d4);
    }

    return 1;
}
//This function calcualtes the GNP vertices for the case when
//there are no shorts edges inside the sample
int GNP_Reconstruction::Calculate_gnp_vertices_case0(const int& R1, const int& R2, const int& R3, const int& R4, const Point_3D& V, GNP& gnp_i, Point_3D& N_top)const
{
    //Get mid point of R1R2
    Point_3D M = (gnp_i.vertices[R1] + gnp_i.vertices[R2]) * 0.5;

    //Calculate vector U, which is the unit vectro along R1R2
    Point_3D U = (gnp_i.vertices[R2] - gnp_i.vertices[R1]).unit();

    //Calculate normal to top face
    N_top = V.cross(U);

    //Calculate displacement from M to R1, R2
    Point_3D disp_u = U * gnp_i.l * 0.5;

    //Calculate positions of R1, R2
    gnp_i.vertices[R1] = M - disp_u;
    gnp_i.vertices[R2] = M + disp_u;

    //Calculate displacement from R1 to R4, which is the same from R2 to R3
    Point_3D disp_v = V * gnp_i.l;

    //Calculate position of R3, R4
    gnp_i.vertices[R3] = gnp_i.vertices[R2] + disp_v;
    gnp_i.vertices[R4] = gnp_i.vertices[R1] + disp_v;

    //Boolean to determine tif vertices R1-R4 are on the top face
    bool top_face = R1 <= 3;

    //Calculate displacement for vertices in opposite face
    Point_3D disp_t = (top_face) ? N_top*(-gnp_i.t) : N_top * gnp_i.t;

    //Determine if should add or subtract 4 to obtain vertices on opposite face
    int N = (top_face) ? 4 : -4;

    //Calculate vetices on opposite face
    gnp_i.vertices[R1 + N] = gnp_i.vertices[R1] + disp_t;
    gnp_i.vertices[R2 + N] = gnp_i.vertices[R2] + disp_t;
    gnp_i.vertices[R3 + N] = gnp_i.vertices[R3] + disp_t;
    gnp_i.vertices[R4 + N] = gnp_i.vertices[R4] + disp_t;

    return 1;
}
