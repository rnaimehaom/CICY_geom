//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GS networks
//OBJECTIVE:    Dlaunay triangulation for a given set of points
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "Triangulation.h"


//Function to make the 3D triangulation
int Triangulation::Generate_3d_trangulation(const vector<Point_3D> &points_gnp, const vector<long int> &structure_i, GNP &gnp_i)
{
    //Check if the number of points to triangulate is 2, 3, or 4, which would be a trivial case
    if (structure_i.size() <= 4) {
        
        //Contruct the trivial case
        if (!Generate_trivial_triangulation(structure_i, gnp_i)) {
            hout<<"Error in Generate_3d_trangulation when calling Generate_trivial_triangulation"<<endl;
            return 0;
        }
    }
    else {
        
        //There are more than 4 points, make the triangulation using the Bowyer-Watson algorithm
        
        //Add points to the triangulation sequentially using the Bowyer-Watson algorithm
        if (!Bowyer_watson(points_gnp, structure_i, gnp_i)) {
            hout<<"Error in Generate_3d_trangulation when calling Bowyer_watson"<<endl;
            return 0;
        }
    }
    
    return 1;
}
//This function generates a trivial triangulation of 2, 3, or 4 vertices to triangulate
//In these trivial cases, I just need to pair one vertex with all others
int Triangulation::Generate_trivial_triangulation(const vector<long int> &structure_i, GNP &gnp_i)
{
    //Iterate over the points in structure_i
    for (int i = 0; i < (int)structure_i.size()-1; i++) {
        
        //First vertex
        long int v1 = structure_i[i];
        
        //Iterate over the remainig points in structure_i
        for (int j = i+1; j < (int)structure_i.size(); j++) {
            
            //Second vertex
            long int v2 = structure_i[j];
            
            //Add and edge with vertices v1 and v2 to the triangulation
            gnp_i.triangulation.push_back(EdgeL(v1, v2));
        }
    }
    
    return 1;
}
//Bowyer Watson algorithm, which generates a triangulation of some vertices by sequentially
//adding each vertex to triangulate
int Triangulation::Bowyer_watson(const vector<Point_3D> &points_gnp, const vector<long int> &structure_i, GNP &gnp_i)
{
    //Generate SUPERTETRAHEDRON
    vector<Point_3D> vertices(4, Point_3D(0,0,0));
    vector<TrFaceL> triangles(4);
    if (!Generate_supertetrahedron(gnp_i, vertices, triangles)) {
        hout<<"Error in Generate_3d_trangulation when calling Generate_supertetrahedron"<<endl;
        return 0;
    }
    
    //Add each point in the GNP (in structure_i) to the triangulation
    for (int i = 0; i < (int)structure_i.size(); i++) {
        
        //Vector for bad triangles
        vector<EdgeL> bad_triangles_edges;
        
        //Get the point number of the vertex added to the triangulation
        long int v = structure_i[i];
        
        //Find the triangles whose circumcircle contains the current point
        //and remove them from the triangulation
        if (!Find_bad_triangles(points_gnp[v], points_gnp, vertices, triangles, bad_triangles_edges)) {
            hout<<"Error in Bowyer_watson when calling Find_bad_triangles"<<endl;
            return 0;
        }
        
        //Find the edges that will form new triangles and add them to the triangulation
        if (!Add_new_triangles(v, triangles, bad_triangles_edges)) {
            hout<<"Error in Bowyer_watson when calling Add_new_triangles"<<endl;
            return 0;
        }
    }
    
    //Triangulation is finished, so find the edges that make up the final triangulation:
    // -Remove those edges that contain a vertex of the supertriangle
    // -Remove repeated edges
    if (!Final_triangulation_edges(triangles, gnp_i)) {
        hout<<"Error in Bowyer_watson when calling Final_triangulation_edges"<<endl;
        return 0;
    }
    
    
    return 1;
}
//This function generates a tetrahedron, which is the initial triangulation
int Triangulation::Generate_supertetrahedron(const GNP &gnp_i, vector<Point_3D> &vertices, vector<TrFaceL> &triangles)
{
    //Calculate the edge length of the regular tetrahedron
    //This edge lenth is such that the GNP is completely contained in the tetrahedron's insphere
    double a = 2.5*sqrt(gnp_i.l*gnp_i.l + gnp_i.t*gnp_i.t/2);
    
    //Some quantities
    double tmp1 = a/sqrt(3);
    
    //Calculate the vertices of the tetrahedron centered in the origin
    //Since vertices is initialized with all points as (0,0,0), just calculate and update
    //the non-zero coordinates
    vertices[0].x = tmp1*a;
    tmp1 = tmp1*0.5;
    //
    vertices[1].x = -tmp1*a;
    vertices[1].y = 0.5*a;
    //
    vertices[2].x = vertices[1].x;
    vertices[2].y = -vertices[1].y;
    //
    vertices[3].z = -vertices[1].x;
    
    //Iterate over all vertices to apply the GNP's translation and rotation
    for (long int i = 0; i < (long int)vertices.size(); i++) {
        
        //Translate and rotate vertex i towards the GNP's location
        vertices[i] = vertices[i].rotation(gnp_i.rotation, gnp_i.center);
    }
    
    //Generate the vector of triangles for the triangulation that consists
    //of only the supertetrahedron vertices
    triangles[0] = TrFaceL(-1, -2, -3);
    triangles[1] = TrFaceL(-1, -2, -4);
    triangles[2] = TrFaceL(-1, -3, -4);
    triangles[3] = TrFaceL(-2, -4, -3);
    
    return 1;
}
//This function finds the bad triangles in a triangulation, i.e., the triangles that contain
//a given point in their circumcircle
//The bad triangles are removed from the triangulation and added as edges
int Triangulation::Find_bad_triangles(const Point_3D &vertex, const vector<Point_3D> &points_gnp, vector<Point_3D> &vertices, vector<TrFaceL> &triangles, vector<EdgeL> &bad_triangles_edges)
{
    //Scan all triangles in reverse order so elements can be deleted from the vector
    for (int j = (int)triangles.size()-1; j >= 0 ; j--) {
        
        //If point is inside the circumcircle then add the triangle to bad_triangles
        if (Is_in_circumcircle(vertex, points_gnp, vertices, triangles[j])) {
            
            //Add current triangle to bad_triangles as edges
            if (!Add_triangle_as_edges(triangles[j], bad_triangles_edges)) {
                hout<<"Error in Find_bad_triangles when calling Add_triangle_as_edges"<<endl;
                return 0;
            }
            
            //Remove triangle from triangulation
            triangles.erase(triangles.begin()+j);
        }
    }
    
    return 1;
}
//This function checks if a point is inside the circumcircle of a triangle
bool Triangulation::Is_in_circumcircle(const Point_3D &vertex_i, const vector<Point_3D> &points_gnp, vector<Point_3D> &vertices, const TrFaceL &triangle)
{
    //Calculate circumcenter of triangle
    Point_3D C = Calculate_circumcenter(points_gnp, vertices, triangle);
    
    //Calculate squared radius of circumcircle (i.e. squared distance from center to any point)
    //Comparing squared distances will save time by avoiding calculating squared roots
    Point_3D P = Get_point(triangle.v1, points_gnp, vertices);
    double rad2 = C.squared_distance_to(P);
    
    //Squared distance from circumcenter to vertex_i
    double dist = C.squared_distance_to(vertex_i);
    
    //If the squared distance between the point and C is smaller or equal than the squared radius,
    //then the point is inside the circumcircle so return 1
    return ( dist - rad2 <= Zero);
}
//Given a triangle, this function calculates its circumcenter
Point_3D Triangulation::Calculate_circumcenter(const vector<Point_3D> &points_gnp, const vector<Point_3D> &vertices, const TrFaceL &triangle)
{
    //Triangle points
    Point_3D P1 = Get_point(triangle.v1, points_gnp, vertices);
    Point_3D P2 = Get_point(triangle.v2, points_gnp, vertices);
    Point_3D P3 = Get_point(triangle.v3, points_gnp, vertices);
    
    //Calculate points P, Q and R = PxQ
    Point_3D P = P1 - P3;
    Point_3D Q = P2 - P3;
    Point_3D R = P.cross(Q);
    
    //Calculate squared norms
    double Pn = P.dot(P);
    double Qn = Q.dot(Q);
    double Rn = R.dot(R);
    
    //Calculate circumcircle
    Point_3D C = (Q*Pn - P*Qn).cross(R)/(2*Rn);
    C = C + P3;
    
    return C;
}
//There are some negative indices to refer to the vertices vectors while all positive indices refer to the point_list vector
//so this function deals with that and returns the proper point
Point_3D Triangulation::Get_point(const long int &index, const vector<Point_3D> &points_gnp, const vector<Point_3D> &vertices)
{
    //If the index is equal or greater to zero, then the index remains the same
    if (index >= 0)
        return points_gnp[index];
    //If the index is negative, add one and change the sign to get the correct index in vertices
    else
        return vertices[-1-index];
}
//This function adds a triangle to the vector edges as the three edges that make up the triangle
int Triangulation::Add_triangle_as_edges(const TrFaceL &triangle, vector<EdgeL> &edges)
{
    //Add edges to last vector of edges
    //Edge 1,2
    edges.push_back(EdgeL(triangle.v1, triangle.v2));
    //Edge 1,3
    edges.push_back(EdgeL(triangle.v2, triangle.v3));
    //Edge 2,3
    edges.push_back(EdgeL(triangle.v3, triangle.v1));
    
    return 1;
}
//This function finds the edges that will make the new triangles in the triangulation
//then the new triangles are added to the triangulation with the edges found and a given point
int Triangulation::Add_new_triangles(const long int &vertex, vector<TrFaceL> &triangles, vector<EdgeL> &bad_triangles_edges)
{
    //Vector that counts the repetitions of each edge
    vector<int> repetitions_count(bad_triangles_edges.size(), 0);
    
    //Scan each edge in bad triangles and compare it with the rest
    for (int j = 0; j < (int)bad_triangles_edges.size()-1; j++) {
        
        //Scan all edges after edge j
        for (int k=j+1; k < (int)bad_triangles_edges.size(); k++) {
            
            //Check if the two edges are the same
            if (bad_triangles_edges[j] == bad_triangles_edges[k]) {
                
                //If there are the same edge, increase the count for edges j and k
                repetitions_count[j] = repetitions_count[j] + 1;
                repetitions_count[k] = repetitions_count[k] + 1;
                
            }
        }
    }
    
    //Scan the vector of repetitions_count and add a new triangle if an edge has 0 repetitions
    for (int i = 0; i < (int)repetitions_count.size(); i++) {
        
        //Check if the count is zero
        if (!repetitions_count[i]) {
            
            //Add a new triangle using the two points of the edge and the current point
            triangles.push_back(TrFaceL(bad_triangles_edges[i].v1, bad_triangles_edges[i].v2, vertex));
        }
    }
    
    return 1;
}
//This function removes the edges that have vertices of the super triangle and
//generates the vector of edges that make up the triangulation
int Triangulation::Final_triangulation_edges(const vector<TrFaceL> &triangles, GNP &gnp_i)
{
    //Vector to store the triangulation edges using the consecutive numbering for vertices
    vector<EdgeL> edges;
    
    //Scan all triangles looking for valid edges
    for (int i = 0; i < (int)triangles.size(); i++) {
        
        //Check if there are any valid edges and add them to the vector of edges
        if (triangles[i].v1 >= 0 && triangles[i].v2 >= 0) {
            
            //Add the edge v1-v2
            edges.push_back(EdgeL(triangles[i].v1, triangles[i].v2));
        }
        if (triangles[i].v2 >= 0 && triangles[i].v3 >= 0) {
            
            //Add the edge v2-v3
            edges.push_back(EdgeL(triangles[i].v2, triangles[i].v3));
        }
        if (triangles[i].v3 >= 0 && triangles[i].v1 >= 0) {
            
            //Add the edge v3-v1
            edges.push_back(EdgeL(triangles[i].v3, triangles[i].v1));
        }
    }
    
    //Vector that counts the repetitions of each edge
    vector<int> repetitions_count(edges.size(), 0);
    
    //Scan each edge in bad triangles and compare it with the rest
    for (int i = 0; i < (int)edges.size()-1; i++) {
        for (int j = i+1; j < (int)edges.size(); j++) {
            
            //If the edge is repeated, then increase the count
            if (edges[i] == edges[j]) {
                repetitions_count[i]++;
                repetitions_count[j]++;
            }
        }
    }
    
    //Add only those edges that have count 0
    for (int i = 0; i < (int)edges.size(); i++) {
        
        //Check the count of edge i
        if (!repetitions_count[i]) {
            
            //Count is zero, so add it to the final triangulation vector in the GNP variable
            gnp_i.triangulation.push_back(edges[i]);
        }
    }
    
    return 1;
}
int Triangulation::Generate_3d_trangulation(const long int &last_cnt_node, const vector<Point_3D> &point_list, const vector<vector<long int> > &structure, const vector<Point_3D> &point_list_gnp, const vector<long int> &gnp_junctions, GCH &hybrid)
{
    
    return 1;
}
