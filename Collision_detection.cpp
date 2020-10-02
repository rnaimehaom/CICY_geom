//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GS networks
//OBJECTIVE:    Detect interpenetrations (collisions) of GS (polyhedra)
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "Collision_detection.h"

//This function determines whethere there is interpenetration of two paralellepipeds
int Collision_detection::GJK(const GNP &gnp1, const GNP &gnp_new, vector<Point_3D> &simplex, bool &p_flag, bool &t_flag)
{
    //Get a point in the Minkowski sum
    Point_3D D = gnp1.vertices[0] - gnp_new.vertices[0];
    
    //Initialize the simplex with the support point of direction D
    Point_3D S = Support_AB(D, gnp1.vertices, gnp_new.vertices);
    simplex.push_back(S);
    
    //Initialize the search direction with the vector from S to the origin
    D = S*(-1);
    
    bool terminate = false;
    
    while (!terminate) {
        
        //Get a new support point
        Point_3D A = Support_AB(D, gnp1.vertices, gnp_new.vertices);
        
        //Check if A is the origin
        if (A.length2() < Zero) {
            
            //Set the flag for touching to true as the polyhedrons are touching at a vertex
            t_flag = true;
            
            //Set the penetration flag to false, as touching is not penetration
            p_flag = false;
            
            //Terminate the function
            return 1;
        }
        //hout<<"A="<<A.str()<<endl;
        
        //Check if there is no intersection
        if (A.dot(D) < Zero) {
            
            //This conditions means there is no penetration, thus set the penetration flag to false
            p_flag = false;
            
            //Terminate the function
            return 1;
        }
        
        //If this part of the code is reached, then update the simplex and check if it
        //contains the origin
        if (Is_origin_in_simplex(simplex, A, D, terminate)) {
            
            //Add the last vertex to the simplex
            simplex.push_back(A);
            if (D.length2() < Zero) {
                //Set the flag for touching to true as the polyhedrons are touching
                t_flag = true;
            }
            
            //If the simplex contains the origin, set the penetration flag to true
            p_flag = true;
            
            //Terminate the function
            return 1;
        }
        //If the simplex was not found to contain the origin, continue the loop with the
        //updated direction D as with the next support point A the simplex might
        //enclose the origin
        
        //If D is a zero vector, then the origin is located at a face, edge or point
        //Thus check for D being a zero vector
        //hout<<"D="<<D.str()<<endl;
        if (D.length2() < Zero) {
            //Set the flag for touching to true as the polyhedrons are touching
            t_flag = true;
            
            //Set the penetration flag to false, as touching is not penetration
            p_flag = false;
            
            //Terminate the function
            return 1;
        }
        //hout<<"END IT"<<endl<<endl;
    }
    
    return 1;
}
//Function to obtain the support point in a given direction in the Minkowski sum
Point_3D Collision_detection::Support_AB(const Point_3D &D, const Point_3D verticesA[], const Point_3D verticesB[])
{
    return (Support_map(D, verticesA) - Support_map(D*(-1), verticesB) );
}
//Function for the suport point in a given direction
Point_3D Collision_detection::Support_map(const Point_3D &D, const Point_3D vertices[])
{
    //Initialize distance with first index
    double dist = D.dot(vertices[0]);
    
    //Initialize index of maximum dot product with 0 (first index)
    int idx = 0;
    
    //Find the maximum distance in the rest of the vertices
    //The array of vertices has a fixed size of 8, so that value is used here
    //It will need to change if other types of polyhedrons are used
    for (int i = 1; i < 8; i++) {
        
        //Calculate the new distance
        double new_dist =  D.dot(vertices[i]);
        
        //Check if update is needed
        if (new_dist > dist) {
            //Update variables
            dist = new_dist;
            idx = i;
        }
    }
    
    //Return the vertex with the maximum distance, which is the one with index idx
    return vertices[idx];
}
//Check if the simplex contains the origin or update it and the search direction
bool Collision_detection::Is_origin_in_simplex(vector<Point_3D> &simplex, Point_3D &A, Point_3D &D, bool &terminate) {
    
    if (simplex.size() == 1) {
        
        //Calculate the vector from A to B, where B is the point in the simplex
        Point_3D AB = simplex[0] - A;
        
        //Calculate the vector from A to the origin
        Point_3D AO = A*(-1);
        
        //hout<<"simplex0={B="<<simplex[0].str()<<"}"<<endl;
        //Check if there is intersection
        if (AB.dot(AO) > Zero) {
            
            //Add A to the simplex in the first position
            simplex.push_back(simplex[0]); //Add another B
            simplex[0] = A; //Turn the first element to A
            //hout<<"simplex={A="<<simplex[0].str()<<", B="<<simplex[1].str()<<"}"<<endl;
            
            //Calculate new direction
            D = (AB.cross(AO)).cross(AB);
        }
        else {
            //Origin is outside Minkowski sum, so update terminate flag
            //There is no intersection
            terminate =  true;
        }
        
        //Return false as either way we cannot say that the origin is inside the simplex
        return false;
    }
    else if (simplex.size() == 2) {
        //Go to the case where the simplex has two elements before (potentially)
        //adding the new vertex A
        return Update_simplex_case2(simplex, A, D, terminate);
    }
    else if (simplex.size() == 3) {
        //Go to the case where the simplex has three elements before (potentially)
        //adding the new vertex A
        return Update_simplex_case3(simplex, A, D, terminate);
    }
    else {
        //hout<<"There was a problem in the code, simplex has "<<simplex.size()<<" elments"<<endl;
        //Also return true to terminate the function GJK
        return true;
    }
}
//This is the case when the simplex already has two points and the support point A is added
bool Collision_detection::Update_simplex_case2(vector<Point_3D> &simplex, Point_3D &A, Point_3D &D, bool &terminate)
{
    
    //cout<<"case 2"<<endl;
    
    //Note that simplex = {B, C}
    Point_3D B = simplex[0];
    Point_3D C = simplex[1];
    //hout<<"simplex0={B="<<simplex[0].str()<<", C="<<simplex[1].str()<<"}"<<endl;
    
    //Precompute some vectors
    Point_3D AB = B-A;
    Point_3D AC = C-A;
    Point_3D AO = A*(-1);
    //ABC = ABxAC
    Point_3D ABC = AB.cross(AC);
    
    //Check where is the origin
    if ( (ABC.cross(AC)).dot(AO) > Zero) {
        
        if (AC.dot(AO) > Zero) {
            //Replace B by A
            simplex[0] = A;
            //Update direction
            D = (AC.cross(AO)).cross(AC);
        }
        else {
            //go to common if-statement
            return Common_if_case2(simplex, A, AB, AO, D, terminate);
        }
    }
    else {
        if ( (AB.cross(ABC)).dot(AO) > Zero) {
            //go to common if-statement
            return Common_if_case2(simplex, A, AB, AO, D, terminate);
        }
        else {
            //Compute D as only the sign will chenge below
            D = ABC;
            
            if (ABC.dot(AO) > Zero) {
                //Update the simplex adding A as the first element
                //Add another C, here simplex = {B, C, C}
                simplex.push_back(simplex.back());
                //Swap B to the middle, here simplex = {B, B, C}
                simplex[1] = simplex[0];
                //Swap A to the front, here simplex = {A, B, C}
                simplex[0] = A;
                //hout<<"simplex={A="<<simplex[0].str()<<", B="<<simplex[1].str()<<", C="<<simplex[2].str()<<"}"<<endl;
            }
            else {
                //Update the simplex adding A as the first element and swapping B and C
                //Add another B, here simplex = {B, C, B}
                simplex.push_back(simplex[0]);
                //Swap the first B with A, here simplex = {A, C, B}
                simplex[0] = A;
                //hout<<"simplex={A="<<simplex[0].str()<<", C="<<simplex[1].str()<<", B="<<simplex[2].str()<<"}"<<endl;
                
                //Update direction
                D = D*(-1);
            }
        }
    }
    
    //Return false as either way we cannot say that the origin is inside the simplex
    return false;
}
//Common if in case 2
bool Collision_detection::Common_if_case2(vector<Point_3D> &simplex, Point_3D &A, Point_3D &AB, Point_3D &AO, Point_3D &D, bool &terminate)
{
    if (AB.dot(AO) > Zero) {
        
        //Update the simplex adding A as the first element, and swapping C by B
        //Swap C by B, here simplex = {B, B}
        simplex[1] = simplex[0];
        //Swap the first B by A, here simplex = {A, B}
        simplex[0] = A;
        
        //Update direction
        D = (AB.cross(AO)).cross(AB);
    }
    else {
        //Origin is outside Minkowski sum, so update terminate flag
        terminate = true;
        //hout<<"Origin is outside Minkowski sum, so update terminate flag (case 2 common)"<<endl;
    }
    
    //Return false as either way we cannot say that the origin is inside the simple
    return false;
}
//This is the case when the simplex already has three points and the support point A is added
bool Collision_detection::Update_simplex_case3(vector<Point_3D> &simplex, Point_3D &A, Point_3D &D, bool &terminate) {
    
    //hout<<"case 3"<<endl;
    
    //Check in which half plane the origin is, and take the triangle in that plane and
    //call case 2
    //Note that if the origin is in the three half-planes that define the volume of the
    //tetrahedron formed by the four points (three in the simplex and A), then the
    //origin is inside and the loop also terminates
    
    //Calculate some important vectors, considering that simplex = {B, C, E}
    Point_3D AB = simplex[0] - A;
    Point_3D AC = simplex[1] - A;
    Point_3D AE = simplex[2] - A;
    Point_3D AO = A*(-1);
    
    //hout<<"simplex0={B="<<simplex[0].str()<<", C="<<simplex[1].str()<<", E="<<simplex[2].str()<<"}"<<endl;
    
    //Check if the origin is outside the tetrahedron in the direction of face AEB
    Point_3D AEB = AE.cross(AB);
    
    //Check for the origin being at the bundary of the Minkowski sum,
    //which means both polyhedra are touching
    if (AEB.length2() < Zero) {
        //When the perpendicular vector is (close to) zero, then there is a touch
        //So set the direction vector equal to zero to avoid numerical error
        D.set(0,0,0);
        
        //The origin is part of the Minkowski sum, however, there is no penetration
        //so return false
        return false;
    }
    
    //hout<<"AEB.dot(AO)="<<AEBdot(AO)<<" AEB="<<AEB.str()<<endl;
    if (AEB.dot(AO) > Zero) {
        //Remove C from the simplex and search in the direction AEB
        //hout<<"Remove C from the simplex and search in the direction AEB"<<endl;
        
        //Here simplex = {B, E, E}
        simplex[1] = simplex[2];
        //Here simplex = {B, E}
        simplex.pop_back();
        //hout<<"simplex={B="<<simplex[0].str()<<", E="<<simplex[1].str()<<"}"<<endl;
        
        //Call case 2
        return Update_simplex_case2(simplex, A, AEB, terminate);
    }
    else {
        //The origin might be inside the tetrahedron, or in other of the faces
        
        //Check if the origin is outside the tetrahedron in the direction of face ABC
        Point_3D ABC = AB.cross(AC);
        
        //Check for the origin being at the bundary of the Minkowski sum,
        //which means both polyhedra are touching
        if (ABC.length2() < Zero) {
            //When the perpendicular vector is (close to) zero, then there is a touch
            //So set the direction vector equal to zero to avoid numerical error
            D.set(0,0,0);
            
            //The origin is part of the Minkowski sum, however, there is no penetration
            //so return false
            return false;
        }
        
        if (ABC.dot(AO) > Zero) {
            //Remove E from the simplex and search in the direction ABC
            //hout<<"Remove E from the simplex and search in the direction ABC"<<endl;
            
            //Swap B and C
            //Here simplex = {B, C, B}
            simplex[2] = simplex[0];
            //Here simplex = {C, C, B}
            simplex[0] = simplex[1];
            //Here simplex = {C, B, B}
            simplex[1] = simplex[2];
            
            //Here simplex = {C, B}
            simplex.pop_back();
            
            //Call case 2
            return Update_simplex_case2(simplex, A, ABC, terminate);
        }
        else {
            //The origin might be inside the tetrahedron, or in other of the faces
            
            //Check if the origin is outside the tetrahedron in the direction of face ACE
            Point_3D ACE = AC.cross(AE);
            
            //Check for the origin being at the bundary of the Minkowski sum,
            //which means both polyhedra are touching
            if (ACE.length2() < Zero) {
                //When the perpendicular vector is (close to) zero, then there is a touch
                //So set the direction vector equal to zero to avoid numerical error
                D.set(0,0,0);
                
                //The origin is part of the Minkowski sum, however, there is no penetration
                //so return false
                return false;
            }
            
            if (ACE.dot(AO) > Zero) {
                //Remove B from the simplex and search in the direction ACE
                //hout<<"Remove B from the simplex and search in the direction ACE"<<endl;
                
                //Here simplex = {C, C, E}
                simplex[0] = simplex [1];
                //Here simplex = {C, E, E}
                simplex[1] = simplex [2];
                //hout<<"simplex={C="<<simplex[0].str()<<", E="<<simplex[1].str()<<"}"<<endl;
                
                //Here simplex = {C, E}
                simplex.pop_back();
                
                //Call case 2
                return Update_simplex_case2(simplex, A, ACE, terminate);
            }
            else {
                //The origin is inside the three faces with common vertex A
                //We already know that the origin is on the halfplane of the remaining
                //triangle BCE where A is located
                //Thus, the origin is in the region contained by the four faces of the
                //tetrahedron, i.e., the origin is inside the tetrahedron
                //Update terminate flag (actually this does not matter but is done to keep
                //consistency)
                terminate = true;
                
                //The origin is found inside the simplex, thus return true
                return true;
            }
        }
    }
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//Extended Polytope Algorithm (EPA)
int Collision_detection::EPA(const Point_3D verticesA[], const Point_3D verticesB[], vector<Point_3D> &simplex, Point_3D &normal, double &PD) {
    
    //Generate the faces of the original simplex (which is a tetrahedron)
    vector<TrFace> Faces(4);
    Faces[0].set(0,1,2);
    Faces[1].set(0,1,3);
    Faces[2].set(0,2,3);
    Faces[3].set(1,2,3);
    
    //Generate a vector of Normals to the faces (in the direction from the origin to the faces)
    //And a vector of distances from the origin to the faces
    vector<Point_3D> Normals(4);
    vector<double> distances(4);
    
    //cout<<"normal_and_distance_to_origin"<<endl;
    for (int i = 0; i < (int)Faces.size(); i++) {
        
        Normal_and_distance_to_origin(simplex, Faces[i], Normals[i], distances[i]);
        //cout<<"dist="<<distances[i]<<endl;
    }
    
    //Enter the loop to expand the polytope
    while (true) {
        
        //Find the face of the simplex closest to the origin
        //Actually I only need the normal and distance
        Point_3D N;
        double PD_;
        //hout<<"find_closest_face"<<endl;
        Find_closest_face(Faces, Normals, distances, N, PD_);
        
        //Get new support point in direction opposite to normal
        Point_3D S = Support_AB(N, verticesA, verticesB);
        //hout<<"S="<<S.str()<<endl;
        
        double dist = S.dot(N);
        
        //hout<<"PD="<<PD_<<" Normal="<<N.str()<<" dist="<<dist<<endl;
        
        //Get the distance of the plane to the new support point
        //This is easily done by turning the normal into a unit vector
        if (dist-PD_ < Zero) {
            
            //If the condition is true, then the distance of the plane to the new support point
            //is the same as the minimum distance of the origin to the plane, so
            //the closest face has been found and terminate the function
            normal = N;
            PD = PD_;
            return 1;
        }
        
        //Add the new support point to the simplex
        simplex.push_back(S);
        /*hout<<"simplex.size="<<simplex.size()<<endl;
        for (int i = 0; i < (int)simplex.size(); i++) {
            hout<<"simplex["<<i<<"]="<<simplex[i].str()<<endl;
        }*/
        int id_S = (int)simplex.size() - 1;
        
        //Reconstruct the simplex
        //hout<<"reconstruct_simplex"<<endl;
        Reconstruct_simplex(simplex, id_S, Faces, Normals, distances);
        
        //hout<<"END while"<<endl<<endl;
    }
    
    return 1;
}
//Calculate the normal unit vector of a given face and its distance to the origin
void Collision_detection::Normal_and_distance_to_origin(const vector<Point_3D> &simplex, const TrFace &f, Point_3D &normal, double &distance) {
    
    //Get vector normal to face
    Point_3D N = (simplex[f.v2] - simplex[f.v1]).cross(simplex[f.v3] - simplex[f.v1]);
    //cout<<"r1="<<(simplex[f.v2] - simplex[f.v1]).str()<<" r2="<<(simplex[f.v3] - simplex[f.v1]).str()<<endl;
    //cout<<"N cross="<<N.str();
    //Make unit
    N.make_unit();
    //cout<<"N unit="<<N.str()<<endl;
    
    //Update the normal going from the origin towards the face
    //The vector v1 goes from the origin towards the face f,
    //then, the vector -v1 goes from the face towards the origin
    //If the dot product is positive, then both vectors "go towards the same direction"
    //So if the dot product of N and -v1 is positive, then N needs to be reversed
    //cout<<"dir check="<<dot(N, simplex[f.v1]*(-1.0))<<endl;
    //Also, save the value of the dot procduct as it is used to calculate the distance
    //from the origin to the face f
    double signed_dist = N.dot(simplex[f.v1]*(-1.0));
    normal = (signed_dist > Zero)? N*(-1.0) : N;
    
    //Calculate distance to origin
    //distance = abs(dot(normal, simplex[f.v1]));
    distance = abs(signed_dist);
    
    //cout<<"Face="<<f.str()<<" N="<<normal.str()<<" dist="<<distance<<endl;
}
//This function finds the face that is closest to the origin
//The output of this function is actually the normal vector and distance to origin of the closest face
void Collision_detection::Find_closest_face(const vector<TrFace> &Faces, const vector<Point_3D> &Normals, const vector<double> &distances, Point_3D &normal, double &PD) {
    
    //Initialize the minimum distance (i.e., the penetration depth PD) with the first face
    PD = distances[0];
    
    //cout<<"Faces.size="<<Faces.size()<<" Normals.size="<<Normals.size()<<" distances.size="<<distances.size()<<endl;
    
    //Set the normal equal to the first element in the vector of normals to keep
    //consistency as the distance was initialized with the first element in the
    //vector of distances
    normal = Normals[0];
    
    int idx = 0;
    
    for (int i = 1 ; i < (int)distances.size(); i++) {
        
        if (distances[i] < PD) {
            //Update minimum distance and its index
            PD = distances[i];
            normal = Normals[i];
            idx = i;
        }
    }
    //cout<<"Closest Face="<<Faces[idx].str()<<endl;
}
//Reconstruct the polyhedron
void Collision_detection::Reconstruct_simplex(const vector<Point_3D> &simplex, const int &S, vector<TrFace> &Faces, vector<Point_3D> &Normals, vector<double> &distances) {
    
    //Vector to store the indices of the faces to remove
    vector<int> idxs;
    
    //Find the faces that are seen by the support point S, these face will be removed
    for (int i = 0; i < (int)Normals.size(); i++) {
        
        //The face is seen by the support point, S, if the vector that goes from the face to S
        //is in the same direction as the normal of the face
        //cout<<"Faces["<<i<<"]="<<Faces[i].str()<<" dot(S-v1,N)="<<dot(simplex[S]-simplex[Faces[i].v1], Normals[i])<<endl;
        if (Normals[i].dot(simplex[S]-simplex[Faces[i].v1]) > Zero ) {
            
            //Add the index of faces to remove
            idxs.push_back(i);
        }
    }
    //cout<<"idxs.size="<<idxs.size()<<endl;
    
    //Vector of edges to store non-shared edges, initialied with the edges of the first face
    vector<Edge> edges(3);
    edges[0].set(Faces[idxs[0]].v1, Faces[idxs[0]].v2);
    edges[1].set(Faces[idxs[0]].v2, Faces[idxs[0]].v3);
    edges[2].set(Faces[idxs[0]].v3, Faces[idxs[0]].v1);
    
    //Check if more than one face was found
    for (int i = 1; i < (int)idxs.size(); i++) {
        
        //Find if any of the edges of the current face is repeated
        
        //Edge 1
        Edge e1 = Edge(Faces[idxs[i]].v1, Faces[idxs[i]].v2);
        int id = Is_edge_in_edges(e1, edges);
        if (id == - 1) {
            //Edge not found, so add it
            edges.push_back(e1);
        }
        else {
            //Edge found, so delete it
            edges.erase(edges.begin()+id);
        }
        
        //Edge 2
        Edge e2 = Edge(Faces[idxs[i]].v2, Faces[idxs[i]].v3);
        id = Is_edge_in_edges(e2, edges);
        if (id == - 1) {
            //Edge not found, so add it
            edges.push_back(e2);
        }
        else {
            //Edge found, so delete it
            edges.erase(edges.begin()+id);
        }
        
        //Edge 3
        Edge e3 = Edge(Faces[idxs[i]].v3, Faces[idxs[i]].v1);
        id = Is_edge_in_edges(e3, edges);
        if (id == - 1) {
            //Edge not found, so add it
            edges.push_back(e3);
        }
        else {
            //Edge found, so delete it
            edges.erase(edges.begin()+id);
        }
    }
    
    //Remove the faces from the polytope
    //Loop goes in rever order in order to avoid deleting the wrong faces
    //as deleting an element will change the index of all elements after the one deleted
    for (int i = (int)idxs.size()-1; i >= 0; i--) {
        //cout<<"Deleted Face="<<Faces[idxs[i]].str()<<endl;
        Faces.erase(Faces.begin()+idxs[i]);
        //Also erase the same elements from the vectors of normals and distances
        Normals.erase(Normals.begin()+idxs[i]);
        distances.erase(distances.begin()+idxs[i]);
    }
    
    //Make new faces using the support point and all the edges that were not shared by
    //the faces that were deleted
    for (int i = 0; i < (int)edges.size(); i++) {
        
        //Make a new face
        Faces.push_back(TrFace(S, edges[i].v1, edges[i].v2));
        
        //Also calculate the distance of the face to the origin and its normal vector
        Point_3D N;
        double dist;
        Normal_and_distance_to_origin(simplex, Faces.back(), N, dist);
        
        //Add the distance to the origin and normal to the corresponding vectors
        Normals.push_back(N);
        distances.push_back(dist);
        //cout<<"New Face="<<Faces.back().str()<<" new_dist="<<distances.back()<<endl;
    }
    
}
//Check if the given edge is in the vector of edges
int Collision_detection::Is_edge_in_edges(const Edge &e, const vector<Edge> &edges) {
    
    for (int i = 0; i < (int)edges.size(); i++) {
        if (edges[i] == e) {
            //A common edge was found, so return the current index
            return i;
        }
    }
    
    //If no repeated index was found return -1 (an invalid index)
    return -1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//This function calculates the distance between two GNPs and finds the closest points on each one
int Collision_detection::Distance_and_direction_estimation(const GNP &gnpA, const GNP &gnpB, Point_3D &N, double &dist)
{
    //Get a point V fron inside the Minkowski sum (e.g., gnpA.center - gnpB.center is
    //inside the Mikonski sum)
    //Then get the vector VO, which goes from V towards the origin.
    //Thus VO = - V = gnpB.center - gnpA.center
    Point_3D VO = gnpB.center - gnpA.center;
    
    //Variables to store the four points of the Minkowski sum that are closest to the origin
    int idxsA[2], idxsB[2];
    
    //Find the two furthest points in A (their indices) in the direction of VO
    Support_map_two_points(VO, gnpA.vertices, idxsA);
    
    //Find the two furthest points in B (their indices) in the direction of -VO
    Support_map_two_points(VO*(-1), gnpB.vertices, idxsB);
    
    //Generate the four points in the Minkowski sum that are closest to the index
    //Closest point to origin in Minkowski sum is in index 0
    Point_3D tetrahedron[4];
    tetrahedron[0] = gnpA.vertices[idxsA[0]] - gnpB.vertices[idxsB[0]];
    tetrahedron[1] = gnpA.vertices[idxsA[0]] - gnpB.vertices[idxsB[1]];
    tetrahedron[2] = gnpA.vertices[idxsA[1]] - gnpB.vertices[idxsB[0]];
    tetrahedron[3] = gnpA.vertices[idxsA[1]] - gnpB.vertices[idxsB[1]];
    
    //From the simplex int the Minkowski sum closest to the origin
    vector<int> simplex;
    if (!Find_simplex_closest_to_origin(tetrahedron, simplex)) {
        hout<<"Error in Distance_estimation when calling Find_simplex_closest_to_origin."<<endl;
        return 0;
    }
    
    //Determine the distance from simplex to origin and direction vector
    if (!Distance_and_direction_from_simplex_to_origin(tetrahedron, simplex, N, dist)) {
        hout<<"Error in Distance_estimation when calling Distance_and_direction_from_simplex_to_origin."<<endl;
        return 0;
    }
    
    return 1;
}
//Function for finding two suport points (their indices) in a given direction
void Collision_detection::Support_map_two_points(const Point_3D &D, const Point_3D vertices[], int idxs[])
{
    //Initialize the first distance with the vertex in the first index
    double dist1 = D.dot(vertices[0]);
    
    //Initialize the first index of maximum dot product with 0 (first index)
    idxs[0] = 0;
    
    //Initialize the second distance with the vertex in the second index
    double dist2 = D.dot(vertices[1]);
    
    //Initialize the second index of maximum dot product depending on the value of idxs[0]
    if (dist1 > dist2) {
        idxs[1] = 1;
    }
    else {
        //Swap distances and indices
        double tmp = dist1;
        dist1 = dist2;
        dist2 = tmp;
        idxs[0] = 1;
        idxs[1] = 0;
    }
    
    //Find the maximum distance in the rest of the vertices
    //The array of vertices has a fixed size of 8, so that value is used here
    //It will need to change if other types of polyhedrons are used
    for (int i = 2; i < 8; i++) {
        
        //Calculate the new distance
        double new_dist =  D.dot(vertices[i]);
        
        //Check if update is needed
        //First check if new_dist is a new maximum
        if (new_dist > dist1) {
            
            //What was stored as dist1, becomes dist2 (i.e., the second maximum)
            dist2 = dist1;
            idxs[1] = idxs[0];
            
            //new_dist becomes the maximum distance
            dist1 = new_dist;
            idxs[0] = i;
        }
        //Second, if new_dist is not a new maximum, check if it is largest that the second maximum
        else if (new_dist > dist2) {
            
            //Only update the second maximum
            dist2 = new_dist;
            idxs[1] = i;
        }
    }
}
//This function finds the simplex (point, edge or face) from a tetrahedron that is closest to the origin
//The simplex is stored in a vector as the indices of the tetrahedron (which is an array of points)
//In this function, the vertices of the tetrahedon are labeled as follows:
//tetrahedron[0] = A
//tetrahedron[1] = B
//tetrahedron[2] = C
//tetrahedron[3] = D
int Collision_detection::Find_simplex_closest_to_origin(const Point_3D tetrahedron[4], vector<int> &simplex)
{
    //Rename tetrahedron vertices to facilitate the writing and reading of the code
    Point_3D A = tetrahedron[0];
    Point_3D B = tetrahedron[1];
    Point_3D C = tetrahedron[2];
    Point_3D D = tetrahedron[3];
    
    //Get a vector from tetrahedron to origin
    Point_3D AO = A*(-1);
    
    //The resulting simplex will for sure have the point A because by construction it is the closest
    //or one of the closest points to the origin
    simplex.push_back(0);
    
    //Determine if the origin is closest to any of the faces
    //Face BCD is ignored because by construction, it cannot be the clossest face to the origin
    //as it does not contain A. By construction A is the clossest point to the origin.
    //So a face that does not contain A cannot be clossest to the origin
    
    //Get a normal to ABC going outside the tetrahedron
    Point_3D ABC = Normal_to_face_ouside_tetrahedron(A, B, C, D);
    
    //If the four points are in the same plane, then ABC is orthogonal to AD and their dot
    //product is zero, in such case, all points are in the simplex,
    //However, only three points are necessary to calculate the distance of the Minkowski sum to the origin
    if (abs(ABC.dot(D-A)) < Zero) {
        
        //Set the three first point in the tetrahedron array as the simplex
        //The first is already added, so just add the next two
        simplex.push_back(1);
        simplex.push_back(2);
        
        //Terminate the function
        return 1;
    }
    
    //Get a normal to ACD going outside the tetrahedron
    Point_3D ACD = Normal_to_face_ouside_tetrahedron(A, C, D, B);
    
    //Get a normal to ABD going outside the tetrahedron
    Point_3D ABD = Normal_to_face_ouside_tetrahedron(A, B, D, C);
    
    //Variable to count the positive dot products
    int n = 0;
    
    //Variable to identify the positive dot products
    vector<int> idxs;
    
    //Carlculate the dot products of the normals to the faces with the vector AO, and check the signs
    //Increase the counter when the dot product is positive
    //Add an identifier when the dot product is positive
    if (AO.dot(ABC) > Zero) {
        n++;
        idxs.push_back(0);
    }
    if (AO.dot(ACD) > Zero) {
        n++;
        idxs.push_back(1);
    }
    if (AO.dot(ABD) > Zero) {
        n++;
        idxs.push_back(2);
    }
    
    //If the three dot products were positive, the the origin is closest to vertex A
    if (n == 3) {
        //There is nothing to do as A was already added to the simplex
        //Just terminate the function
        return 1;
    }
    //If two dot products were positive, then the simplex is an edge
    else if (n == 2) {
        
        //Determine the second vertex of the edge
        if (idxs[0] == 0) {
            if (idxs[1] == 1) {
                
                //This is the case when the origin is above planes for faces ABC and ACD
                //Those faces share the edge AC, thus this is the simplex
                //Since A is already in the simplex, add C (2)
                simplex.push_back(2);
            }
            else if (idxs[1] == 2) {
                
                //This is the case when the origin is above planes for faces ABC and ABD
                //Those faces share the edge AB, thus this is the simplex
                //Since A is already in the simplex, add B (1)
                simplex.push_back(1);
            }
            else {
                //This should no happen, so return error
                hout<<"Error in Find_simplex_closest_to_origin. Vector idxs has an invalid combination:"<<endl;
                for (size_t i = 0; i < idxs.size(); i++) {
                    hout<<"idxs["<<i<<"]="<<idxs[i]<<endl;
                }
                hout<<endl;
                return 0;
            }
            
            return 1;
        }
        else if (idxs[0] == 1 && idxs[2] == 2) {
            
            //This is the case when the origin is above planes for faces ACD and ABD
            //Those faces share the edge AD, thus this is the simplex
            //Since A is already in the simplex, add D (3) and terminate the function
            simplex.push_back(3);
            return 1;
        }
        else {
            //This should no happen, so return error
            hout<<"Error in Find_simplex_closest_to_origin. Vector idxs has an invalid combination:"<<endl;
            for (size_t i = 0; i < idxs.size(); i++) {
                hout<<"idxs["<<i<<"]="<<idxs[i]<<endl;
            }
            hout<<endl;
            return 0;
        }
    }
    //If only one dot product was positive, then the simplex is a face
    else if(n == 1) {
        
        //Determine the face from idxs
        if (idxs[0] == 0) {
            //Origin is above face ABC
            //Add B and C
            simplex.push_back(1);
            simplex.push_back(2);
        }
        else if (idxs[0] == 1) {
            //Origin is above face ACD
            //Add C and D
            simplex.push_back(2);
            simplex.push_back(3);
        }
        else if (idxs[0] == 2) {
            //Origin is above face ABD
            //Add B and D
            simplex.push_back(1);
            simplex.push_back(3);
        }
        else {
            //This should no happen, so return error
            hout<<"Error in Find_simplex_closest_to_origin. The first element in idxs has an invalid value:"<<idxs[0]<<". Valid values are 0, 1, or 2 only."<<endl;
            return 0;
        }
        
        return 1;
    }
    else {
    
        //This should not happen
        hout<<"Error in Find_simplex_closest_to_origin. There were n="<<n<<" positive dot products. Valid values of n can only be 1, 2, or 3."<<endl;
        return 0;
    }
}
//
Point_3D Collision_detection::Normal_to_face_ouside_tetrahedron(const Point_3D &A, const Point_3D &B, const Point_3D &C, const Point_3D &not_in_face)
{
    //Get a normal to face ABC
    Point_3D ABC = (C - A).cross(B - A);
    
    //Make sure ABC goes outside the tetrahedron, e.g., in the oposite direction of AD,
    //where D is not_in_face
    Point_3D AD = not_in_face - A;
    
    //Calculate dot product
    double dot_tmp = ABC.dot(AD);
    if (dot_tmp > Zero) {
        //Reverse the direction is ABC and AD go in the same direction, i.e.,
        //ABC goes inside the tetrahedron
        ABC = ABC*(-1);
    }
    
    return ABC;
}
int Collision_detection::Distance_and_direction_from_simplex_to_origin(const Point_3D tetrahedron[4], const vector<int> &simplex, Point_3D &N, double &dist)
{
    //Vector AO is used a few times
    Point_3D AO = tetrahedron[simplex[0]]*(-1);
    
    //Calculate the distance from simplex to origin depending on the type of simplex (vertex, edge, or face)
    if (simplex.size() == 3) {
        //Simplex is a face
        
        //Get the normal vector
        N = (tetrahedron[simplex[1]] - tetrahedron[simplex[0]]).cross(tetrahedron[simplex[2]] - tetrahedron[simplex[0]]);
        //Make the vector unitary
        N.make_unit();
        
        //Make sure it goes towards the origin, e.g., in the same direction as AO
        double signed_dist = N.dot(AO*(-1));
        N = (signed_dist > Zero)? N*(-1.0) : N;
        
        //Calculate distance to origin
        dist = abs(signed_dist);
    }
    else if (simplex.size() == 2) {
        //Simplex is an edge
        
        //Calculate the projection of AO on AB, where simplez = {A,B}
        Point_3D AB = tetrahedron[simplex[1]] - tetrahedron[simplex[0]];
        double dot_tmp = AO.dot(AB);
        //This is actually the projection squared, this is the quantity needed to obtain dist
        double proj2 = dot_tmp*dot_tmp/AB.length2();
        
        //The distance of the edge to the origin is calculated using Pythagoras
        dist = sqrt(AO.length2() - proj2);
        
        //The direction vector is the triple cross product
        N = (AB.cross(AO)).cross(AB);
        
    }
    else if (simplex.size() == 1) {
        //Simplex is a vertex
        
        //The direction vector goes from the vertex to origin, i.e., along AO
        N = AO.unit();
        
        //The distance is the length of the vector from origin to point
        dist = tetrahedron[simplex[0]].length();
    }
    else {
        //This should not happen
        hout<<"Error in Distance_and_direction_from_simplex_to_origin. Simplex has an invalid number of elements: "<<simplex.size()<<" elements. Simplex can only have 1, 2, o 3 elements."<<endl;
        return 0;
    }
    return 1;
}
